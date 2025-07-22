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

#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "generators/CellBlockManagerABC.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
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
  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    MeshBody & meshBody = domain.getMeshBodies().registerGroup< MeshBody >( meshGen.getName() );
    meshBody.createMeshLevel( 0 );
    SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );

    meshGen.generateMesh( meshBody, partition );

    if( !meshBody.hasParticles() )
    {
      CellBlockManagerABC const & cellBlockManager = meshBody.getCellBlockManager();

      meshBody.setGlobalLengthScale( cellBlockManager.getGlobalLength() );
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
    printf("here0\n");
    if( !domain.hasMeshBody( generator.getName() ) )
    {
    printf("here1\n");
      return;
    }
    else if( domain.getMeshBody( generator.getName() ).hasParticles() ) // field import is not currently compatible with particle mesh
                                                                        // bodies
    {
    printf("here2\n");
      return;
    }

    printf("here3\n");
    //GEOS_LOG_RANK_0( GEOS_FMT( "{}: importing field data from mesh dataset", generator.getName() ) );
    printf("here4\n");
    MeshBody & meshBody = domain.getMeshBody( generator.getName() );
    printf("here5\n");
    meshBody.forMeshLevels( [&]( MeshLevel & meshLevel )
    {
    printf("here6\n");

    printf("meshLevel.getName()=%s\n", meshLevel.getName().c_str());
      //GEOS_LOG_RANK_0( GEOS_FMT( "  mesh level = {}", meshLevel.getName() ) );

    printf("here7\n");
      FieldIdentifiers fieldsToBeSync;
    printf("here7bis\n");
      meshLevel.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >(
        [&]( localIndex,
             localIndex,
             ElementRegionBase const & region,
             CellElementSubRegion & subRegion )
      {
    printf("here8\n");
    printf("region.getName()=%s\n", region.getName().c_str());
    printf("subRegion.getName()=%s\n",subRegion.getName().c_str());
    printf("still here8\n");
    //GEOS_LOG_RANK_0( GEOS_FMT( "  volumic fields on {}/{}", region.getName(), subRegion.getName() ) );
    printf("here9\n");
        importFields( generator, region.getName(), subRegion, MeshGeneratorBase::Block::VOLUMIC, generator.getVolumicFieldsMapping(), fieldsToBeSync );
        printf("afterimportfield\n");
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
      printf("herecomm\n");
      CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, meshLevel, domain.getNeighbors(), false ); // TODO Validate this.
      printf("hereaftercomm\n");
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
  printf("beforeloopimportfield\n");
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
  printf("afterloopimportfield\n");

    // Now that we know that the subRegion has this wrapper,
    // we can add the geosFieldName to the list of fields to synchronize
    printf("ligne1\n");
    fieldsToBeSync.addElementFields( { geosFieldName }, { regionName } );
    printf("ligne2\n");
    WrapperBase & wrapper = subRegion.getWrapperBase( geosFieldName );
    printf("ligne3\n");
    GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::ImportFields,
                                    GEOS_FMT( "    {} -> {}", meshFieldName, geosFieldName ),
                                    generator );

    bool const isMaterialField = materialWrapperNames.count( geosFieldName ) > 0 && wrapper.numArrayDims() > 1;
    printf("ligne4\n");
    generator.importFieldOnArray( block, subRegion.getName(), meshFieldName, isMaterialField, wrapper );
    printf("ligne5\n");
  }
}

} /* namespace geos */
