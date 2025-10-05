/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "MeshManager.hpp"
#include "MeshBody.hpp"
#include "MeshLevel.hpp"
#include "mesh/LogLevelsInfo.hpp"

#include "generators/CellBlockManagerABC.hpp"
#include "generators/InternalMeshGenerator.hpp"
#include "generators/ExternalMeshGeneratorBase.hpp"
#include "generators/ParticleMeshGenerator.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/PartitionerManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "common/TimingMacros.hpp"

#include <unordered_set>

namespace geos
{

using namespace dataRepository;

MeshManager::MeshManager( string const & name,
                          Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );

  addLogLevel< logInfo::ImportFields >();
}

MeshManager::~MeshManager()
{}

Group * MeshManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< MeshGeneratorBase > mesh =
    MeshGeneratorBase::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &this->registerGroup< MeshGeneratorBase >( childName, std::move( mesh ) );
}


void MeshManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from MeshGeneratorBase here
  for( auto & catalogIter: MeshGeneratorBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


void MeshManager::generateMeshes( DomainPartition & domain )
{
  // First, check if there's any work to do at all.
  // Count the number of mesh generators registered.
  int numMeshGenerators = 0;
  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase const & ) { ++numMeshGenerators; } );

  // If there are no mesh generators, this is a meshless run.
  // This function has nothing to do, so return immediately.
  if( numMeshGenerators == 0 )
  {
    GEOS_LOG_RANK_0( "No mesh generators found in MeshManager. Assuming meshless simulation." );
    return;
  }

  // Temporary partitioner manager for unit tests (will be deleted at end of function)
  std::unique_ptr< PartitionerManager > tempPartitionerManager;

  // If no partitioner exists, create a default one based on mesh type
  // This might happen in tests where there's no ProblemManager parent
  if( !domain.hasPartitioner()  )
  {
    // Check if the domain has a parent (ProblemManager)
    PartitionerManager * partitionerManager = nullptr;

    if( domain.hasParent() )
    {
      Group & problemManager = domain.getParent();
      if( problemManager.hasGroup( ProblemManager::groupKeysStruct().partitionerManager.key() ) )
      {
        partitionerManager = &problemManager.getGroup< PartitionerManager >( ProblemManager::groupKeysStruct().partitionerManager.key() );
      }
    }

    // If no PartitionerManager found, create a temporary one (for unit tests)
    if( partitionerManager == nullptr )
    {
      GEOS_LOG_RANK_0( "No PartitionerManager available (likely a unit test). Creating temporary PartitionerManager." );
      tempPartitionerManager = std::make_unique< PartitionerManager >( "tempPartitionerManager", &domain );
      partitionerManager = tempPartitionerManager.get();
    }

    bool hasInternalMesh = false;
    bool hasExternalMesh = false;
    bool hasParticleMesh = false;

    // Check what type of mesh generators are being used
    forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase const & meshGen )
    {
      if( dynamic_cast< InternalMeshGenerator const * >( &meshGen ) != nullptr )
      {
        hasInternalMesh = true;
      }
      else if( dynamic_cast< ExternalMeshGeneratorBase const * >( &meshGen ) != nullptr )
      {
        hasExternalMesh = true;
      }
      else if( dynamic_cast< ParticleMeshGenerator const * >( &meshGen ) != nullptr )
      {
        hasParticleMesh = true;
      }
    } );

    // Create default partitioner based on mesh type
    if( hasParticleMesh )
    {
      GEOS_LOG_RANK_0( "No partitioner specified. Creating default ParticleCartesianPartitioner for particle mesh." );
      partitionerManager->createChild( "ParticleCartesian", "defaultPartitioner" );
    }
    else if( hasInternalMesh && hasExternalMesh )
    {
      GEOS_ERROR( "Both internal and external meshes detected, but no partitioner specified." );
    }
    else if( hasInternalMesh )
    {
      GEOS_LOG_RANK_0( "No partitioner specified. Creating default CartesianPartitioner for internal mesh." );
      partitionerManager->createChild( "Cartesian", "defaultPartitioner" );
    }
    else if( hasExternalMesh )
    {
      GEOS_LOG_RANK_0( "No partitioner specified. Creating default ParMetisPartitioner for external mesh." );
      partitionerManager->createChild( "ParMetis", "defaultPartitioner" );
    }
    else
    {
      GEOS_ERROR( "No partitioner defined." );
    }
  }

  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    MeshBody & meshBody = domain.getMeshBodies().registerGroup< MeshBody >( meshGen.getName() );
    meshBody.createMeshLevel( 0 );

    // Get the partitioner (either from ProblemManager or temporary)
    PartitionerBase * partitioner = nullptr;
    if( domain.hasPartitioner() )
    {
      partitioner = &domain.getPartitioner();
    }
    else if( tempPartitionerManager != nullptr )
    {
      partitioner = &tempPartitionerManager->getPartitioner();
    }

    if( partitioner != nullptr )
    {
      meshGen.generateMesh( meshBody, *partitioner );
    }
    else
    {
      GEOS_ERROR( "No partitioner available for mesh generation" );
    }

    if( !meshBody.hasParticles() )
    {
      CellBlockManagerABC const & cellBlockManager = meshBody.getCellBlockManager();

      meshBody.setGlobalLengthScale( std::max( cellBlockManager.getGlobalLength(), cellBlockManager.getGlobalOffset() ) );
    }
  } );
}
void MeshManager::generateMeshLevels( DomainPartition & domain )
{
  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    string const & meshName = meshGen.getName();
    domain.getMeshBodies().registerGroup< MeshBody >( meshName ).createMeshLevel( MeshBody::groupStructKeys::baseDiscretizationString() );
  } );
}

/**
 * @brief Collect a set of material field names registered in a subregion.
 * @param subRegion the target subregion
 * @return a set of wrapper names
 */
std::unordered_set< string > getMaterialWrapperNames( ElementSubRegionBase const & subRegion )
{
  using namespace constitutive;
  std::unordered_set< string > materialWrapperNames;
  subRegion.getConstitutiveModels().forSubGroups< ConstitutiveBase >( [&]( ConstitutiveBase const & material )
  {
    material.forWrappers( [&]( WrapperBase const & wrapper )
    {
      if( wrapper.sizedFromParent() )
      {
        materialWrapperNames.insert( ConstitutiveBase::makeFieldName( material.getName(), wrapper.getName() ) );
      }
    } );
  } );
  return materialWrapperNames;
}

void MeshManager::importFields( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  forSubGroups< MeshGeneratorBase >( [&domain]( MeshGeneratorBase const & generator )
  {
    if( !domain.hasMeshBody( generator.getName() ) )
    {
      return;
    }
    else if( domain.getMeshBody( generator.getName() ).hasParticles() ) // field import is not currently compatible with particle mesh
                                                                        // bodies
    {
      return;
    }

    GEOS_LOG_RANK_0( GEOS_FMT( "{}: importing field data from mesh dataset", generator.getName() ) );
    MeshBody & meshBody = domain.getMeshBody( generator.getName() );
    meshBody.forMeshLevels( [&]( MeshLevel & meshLevel )
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "  mesh level = {}", meshLevel.getName() ) );
      FieldIdentifiers fieldsToBeSync;
      meshLevel.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >(
        [&]( localIndex,
             localIndex,
             ElementRegionBase const & region,
             CellElementSubRegion & subRegion )
      {
        GEOS_LOG_RANK_0( GEOS_FMT( "  volumic fields on {}/{}", region.getName(), subRegion.getName() ) );
        importFields( generator, region.getName(), subRegion, MeshGeneratorBase::Block::VOLUMIC, generator.getVolumicFieldsMapping(), fieldsToBeSync );
      } );
      meshLevel.getElemManager().forElementSubRegionsComplete< FaceElementSubRegion >(
        [&]( localIndex,
             localIndex,
             ElementRegionBase const & region,
             FaceElementSubRegion & subRegion )
      {
        GEOS_LOG_RANK_0( GEOS_FMT( "  surfaic fields on {}/{}", region.getName(), subRegion.getName() ) );
        importFields( generator, region.getName(), subRegion, MeshGeneratorBase::Block::SURFACIC, generator.getSurfacicFieldsMapping(), fieldsToBeSync );
      } );
      CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, meshLevel, domain.getNeighbors(), false ); // TODO Validate this.
    } );
  } );

  forSubGroups< MeshGeneratorBase >( []( MeshGeneratorBase & generator )
  {
    generator.freeResources();
  } );
}

void MeshManager::importFields( MeshGeneratorBase const & generator,
                                string const & regionName,
                                ElementSubRegionBase & subRegion,
                                MeshGeneratorBase::Block const block,
                                std::map< string, string > const & fieldsMapping,
                                FieldIdentifiers & fieldsToBeSync )
{
  std::unordered_set< string > const materialWrapperNames = getMaterialWrapperNames( subRegion );
  // Writing properties
  for( auto const & pair : fieldsMapping )
  {
    string const & meshFieldName = pair.first;
    string const & geosFieldName = pair.second;
    // Find destination
    if( !subRegion.hasWrapper( geosFieldName ) )
    {
      // Skip - the user may have not enabled a particular physics model/solver on this destination region.
      GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::ImportFields,
                                      GEOS_FMT( "    Skipping import of {} -> {} (field not found)", meshFieldName, geosFieldName ),
                                      generator );

      continue;
    }

    // Now that we know that the subRegion has this wrapper,
    // we can add the geosFieldName to the list of fields to synchronize
    fieldsToBeSync.addElementFields( { geosFieldName }, { regionName } );
    WrapperBase & wrapper = subRegion.getWrapperBase( geosFieldName );
    GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::ImportFields,
                                    GEOS_FMT( "    {} -> {}", meshFieldName, geosFieldName ),
                                    generator );

    bool const isMaterialField = materialWrapperNames.count( geosFieldName ) > 0 && wrapper.numArrayDims() > 1;
    generator.importFieldOnArray( block, subRegion.getName(), meshFieldName, isMaterialField, wrapper );
  }
}

} /* namespace geos */
