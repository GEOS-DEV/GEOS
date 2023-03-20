/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "MeshManager.hpp"
#include "MeshBody.hpp"
#include "MeshLevel.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "generators/MeshGeneratorBase.hpp"
#include "generators/CellBlockManager.hpp"
#include "common/TimingMacros.hpp"

#include <unordered_set>

namespace geosx
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
  GEOSX_LOG_RANK_0( "Adding Mesh: " << childKey << ", " << childName );
  std::unique_ptr< MeshGeneratorBase > solver = MeshGeneratorBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< MeshGeneratorBase >( childName, std::move( solver ) );
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
    MeshLevel & meshLevel = meshBody.getBaseDiscretization();
    CellBlockManager & cellBlockManager = meshBody.registerGroup< CellBlockManager >( keys::cellManager );

    meshGen.generateMesh( cellBlockManager );
    // std::tuple< std::unique_ptr< CellBlockManagerABC >, PartitionDescriptor const &, real64 > tup =  meshGen.generateMesh( );
    // std::unique_ptr< CellBlockManagerABC > cellBlockManager( std::move( std::get<0>( tup ) ) );
    // string name = keys::cellManager;
    // CellBlockManagerABC & cbm = meshBody.registerGroup< CellBlockManagerABC >( name,  std::move( cellBlockManager ) );
    // PartitionDescriptor const & partitionDescriptor = std::get<1>( tup );
    // real64 globalLength = std::get<2>( tup );
    PartitionDescriptor const & partitionDescriptor = cellBlockManager.getPartitionDescriptor();
    real64 globalLength = cellBlockManager.getGlobalLength();
    meshBody.setGlobalLengthScale( globalLength );
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    //PartitionDescriptor & partitionDescriptor = cellBlockManager.getPartitionDescriptor();
    if ( partitionDescriptor.hasMetisNeighborList() )
    {
      domain.getMetisNeighborList() = partitionDescriptor.getMetisNeighborList();
    }
    else
    {
      SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );
      partition = partitionDescriptor.getSpatialPartition();
    }
    meshBody.setGlobalLengthScale( globalLength );
  } );
}


void MeshManager::generateMeshLevels( DomainPartition & domain )
{
  this->forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & meshGen )
  {
    if( dynamicCast< InternalWellGenerator * >( &meshGen ) )
    {
      return;
    }

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
  GEOSX_MARK_FUNCTION;
  forSubGroups< MeshGeneratorBase >( [&]( MeshGeneratorBase & generator )
  {
    if( !domain.hasMeshBody( generator.getName() ) )
      return;

    FieldIdentifiers fieldsToBeSync;
    std::map< string, string > fieldNamesMapping = generator.getFieldsMapping();

    MeshLevel & meshLevel = domain.getMeshBody( generator.getName() ).getBaseDiscretization();
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementSubRegionsComplete< CellElementSubRegion >( [&]( localIndex,
                                                                           localIndex,
                                                                           ElementRegionBase const & region,
                                                                           CellElementSubRegion & subRegion )
    {
      std::unordered_set< string > const materialWrapperNames = getMaterialWrapperNames( subRegion );
      // Writing properties
      for( auto const & pair : fieldNamesMapping )
      {
        string const & meshFieldName = pair.first;
        string const & geosxFieldName = pair.second;
        // Find destination
        if( !subRegion.hasWrapper( geosxFieldName ) )
        {
          // Skip - the user may have not enabled a particular physics model/solver on this dstRegion.
          GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Skipping import of {} -> {} on {}/{} (field not found)",
                                                meshFieldName, geosxFieldName, region.getName(), subRegion.getName() ) );

          continue;
        }

        // Now that we know that the subRegion has this wrapper, we can add the geosxFieldName to the list of fields to
        // synchronize
        fieldsToBeSync.addElementFields( {geosxFieldName}, {region.getName()} );
        WrapperBase & wrapper = subRegion.getWrapperBase( geosxFieldName );
        GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Importing field {} -> {} on {}/{}",
                                              meshFieldName, geosxFieldName,
                                              region.getName(), subRegion.getName() ) );

        bool const isMaterialField = materialWrapperNames.count( geosxFieldName ) > 0 && wrapper.numArrayDims() > 1;
        generator.importFieldsOnArray( region.getName(), meshFieldName, isMaterialField, wrapper );
      }
    } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, meshLevel, domain.getNeighbors(), false );
    generator.freeResources();
  } );
}

} /* namespace geosx */
