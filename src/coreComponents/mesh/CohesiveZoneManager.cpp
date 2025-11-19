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

#include <vector>

#include "CohesiveZoneManager.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "constitutive/ConstitutiveManager.hpp"

#include "mesh/MeshManager.hpp"
#include "schema/schemaUtilities.hpp"

namespace geos
{

using namespace dataRepository;

class SpatialPartition;

// *********************************************************************************************************************
/**
 * @return
 */

CohesiveZoneManager::CohesiveZoneManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
  this->registerGroup< Group >( CohesiveZoneManager::groupKeyStruct::cohesiveZoneRegionsGroup() );
}

CohesiveZoneManager::~CohesiveZoneManager()
{
  // TODO Auto-generated destructor stub
}

localIndex CohesiveZoneManager::numRegions() const
{
  localIndex numRegions = 0;
  this->forCohesiveZoneRegions< CohesiveZoneRegionBase >( [&]( CohesiveZoneRegionBase const & )
  {
    numRegions += 1;
  } );
  return numRegions;
}

void CohesiveZoneManager::resize( integer_array const & numNodes,
                                  string_array const & regionNames )
{
  localIndex const n_regions = LvArray::integerConversion< localIndex >( regionNames.size() );
  for( localIndex reg = 0; reg < n_regions; ++reg )
  {
    this->getRegion( reg ).resize( numNodes[reg] );
  }
}

Group * CohesiveZoneManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  Group & cohesiveZoneRegions = this->getGroup( CohesiveZoneManager::groupKeyStruct::cohesiveZoneRegionsGroup() );
  return &cohesiveZoneRegions.registerGroup( childName,
                                             CatalogInterface::factory( childKey, getDataContext(),
                                                                        childName, &cohesiveZoneRegions ) );
}

void CohesiveZoneManager::expandObjectCatalogs()
{
  ObjectManagerBase::CatalogInterface::CatalogType const & catalog = ObjectManagerBase::getCatalog();
  for( ObjectManagerBase::CatalogInterface::CatalogType::const_iterator iter = catalog.begin();
       iter!=catalog.end();
       ++iter )
  {
    string const key = iter->first;
    if( key.find( "CohesiveZoneRegion" ) != string::npos )
    {
      this->createChild( key, key );
    }
  }
}

void CohesiveZoneManager::setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                               xmlWrapper::xmlNode schemaParent,
                                               integer documentationType )
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child( "xsd:choice" );
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child( "xsd:choice" );
    targetChoiceNode.append_attribute( "minOccurs" ) = "0";
    targetChoiceNode.append_attribute( "maxOccurs" ) = "unbounded";
  }

  std::set< string > names;
  this->forCohesiveZoneRegions( [&]( CohesiveZoneRegionBase & cohesiveZoneRegion )
  {
    names.insert( cohesiveZoneRegion.getName() );
  } );

  for( string const & name: names )
  {
    schemaUtilities::SchemaConstruction( getRegion( name ), schemaRoot, targetChoiceNode, documentationType );
  }
}

void CohesiveZoneManager::syncCohesiveZones()
{
  // this->forCohesiveZoneRegions< CohesiveZoneRegionBase >( [&]( CohesiveZoneRegionBase const & czRegion )
  // {
  //   arrayView1d< globalIndex > const czGlobalID = czRegion.getGloblaID();
  //   arrayView2d< real64 > const czReferencePosition = czRegion.getReferencePosition();
  //   arrayView3d< real64 > const czReferenceSurfaceNormal = czRegion.getReferenceSurfaceNormal();

  //   // Sync initial grid positions with rank 0
  //   int const numRanks = MpiWrapper::commSize();
  //   int rank = MpiWrapper::commRank( MPI_COMM_GEOS );

  //   array1d< MPI_Request > mpiRequestIndices( numRanks );
  //   array1d< MPI_Request > mpiRequestPositions( numRanks );
  //   array1d< MPI_Request > mpiRequestSurfaceNormals( numRanks );

  //   array1d< MPI_Status > mpiStatusIndices( numRanks );
  //   array1d< MPI_Status > mpiStatusPositions( numRanks );
  //   array1d< MPI_Status > mpiStatusSurfaceNormals( numRanks );

  //   array1d< int > gridNodeIndices;
  //   if( rank != 0 )
  //   {
  //     // Send list of indices to overwrite nodal positions
  //     mpiRequestIndices[rank] = MPI_REQUEST_NULL;
  //     MpiWrapper::iSend( localCohesiveGridNodes.data(),
  //                        localCohesiveGridNodes.size(),
  //                        0,
  //                        0,
  //                        MPI_COMM_GEOS,
  //                        &mpiRequestIndices[rank] );

  //     // Send nodal positions
  //     mpiRequestPositions[rank] = MPI_REQUEST_NULL;
  //     MpiWrapper::iSend( m_referenceCohesiveGridNodePositions.data(),
  //                       m_referenceCohesiveGridNodePositions.size(),
  //                       0,
  //                       1,
  //                       MPI_COMM_GEOS,
  //                       &mpiRequestPositions[rank] );

  //     // Send nodal positions
  //     mpiRequestSurfaceNormals[rank] = MPI_REQUEST_NULL;
  //     MpiWrapper::iSend( m_referenceCohesiveGridNodePartitioningSurfaceNormals.data(),
  //                       m_referenceCohesiveGridNodePartitioningSurfaceNormals.size(),
  //                       0,
  //                       1,
  //                       MPI_COMM_GEOS,
  //                       &mpiRequestSurfaceNormals[rank] );
  //   }
  //   else
  //   {
  //     for( int r = 1; r < numRanks; r++ )
  //     {
  //       MpiWrapper::recv( gridNodeIndices,
  //                         r,
  //                         0,
  //                         MPI_COMM_GEOS,
  //                         &mpiStatusIndices[r] );

  //       array2d< real64 > gridNodePositions( numCohesiveNodes, 3 );
  //       mpiRequestPositions[r] = MPI_REQUEST_NULL;
  //       MpiWrapper::iRecv( gridNodePositions.data(),
  //                         gridNodePositions.size(),
  //                         r,
  //                         1,
  //                         MPI_COMM_GEOS,
  //                         &mpiRequestPositions[r] );

  //       MpiWrapper::wait( &mpiRequestPositions[r], &mpiStatusPositions[r] );

  //       // Send nodal surface normals
  //       array2d< real64 > gridNodeSurfaceNormals( numCohesiveNodes, 3 );
  //       mpiRequestSurfaceNormals[r] = MPI_REQUEST_NULL;
  //       MpiWrapper::iRecv( gridNodeSurfaceNormals.data(),
  //                         gridNodeSurfaceNormals.size(),
  //                         r,
  //                         1,
  //                         MPI_COMM_GEOS,
  //                         &mpiRequestSurfaceNormals[r] );

  //       MpiWrapper::wait( &mpiRequestSurfaceNormals[r], &mpiStatusSurfaceNormals[r] );

  //       // Combine grid positions
  //       for( int g = 0; g < gridNodeIndices.size(); ++g )
  //       {
  //         for( int k = 0; k < 3; ++k )
  //         {
  //           referenceCohesiveGridNodePositions[gridNodeIndices[g]][k] = gridNodePositions[gridNodeIndices[g]][k];
  //           // TODO: Check if we should be taking the max here for the partitioning surface normal?
  //           referenceCohesiveGridNodePartitioningSurfaceNormals[gridNodeIndices[g]][k] = gridNodeSurfaceNormals[gridNodeIndices[g]][k];
  //         }
  //       }
  //     }
  //   }
  // } );


}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CohesiveZoneManager, string const &, Group * const )
}
