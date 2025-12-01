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
#include "physicsSolvers/solidMechanics/SolidMechanicsMPM.hpp"
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

void CohesiveZoneManager::initializeReferenceConfiguration( NodeManager & nodeManager,
                                                            ParticleManager & particleManager )
{
  GEOS_UNUSED_VAR( particleManager );

  // Grid fields
  // int const numNodes = nodeManager.size();
  auto & nodeGlobalToLocalMap = nodeManager.globalToLocalMap();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const gridExplicitSurfaceNormal = nodeManager.getReference< array2d< real64 > >( SolidMechanicsMPM::viewKeyStruct::gridExplicitSurfaceNormalString() );

  this->forCohesiveZoneRegions< CohesiveZoneRegionBase >( [&]( CohesiveZoneRegionBase & czRegion )
  {
    SortedArrayView< globalIndex const > const globalID = czRegion.getGlobalID();
    arrayView2d< real64 > referencePosition = czRegion.getReferencePosition();
    arrayView2d< real64 > referencePartitioningSurfaceNormal = czRegion.getReferencePartitioningSurfaceNormal();
    
    array1d< localIndex > indexInGlobalCZArray; // Perhaps this should be initialized to size given we know the count of local cohesive grid nodes from beginning of initialization
    array1d< localIndex > localNodalIndex;

    // Following in current form must be done in serial
    // Unsure how best to parallelize it without race conditions or loss of efficiency
    // This is also only be performde once at CZ initialization
    int const numCohesiveNodes = czRegion.size();
    for( int g  = 0; g < numCohesiveNodes; ++g )
    {
      // Check if global index of cohesive zone is on this partition
      // if( nodeGlobalToLocalMap.contains( globalID[g] ) ) // C++ 20
      if( nodeGlobalToLocalMap.find( globalID[g] ) != nodeGlobalToLocalMap.end() )
      {
        indexInGlobalCZArray.emplace_back( g );
        localNodalIndex.emplace_back( nodeGlobalToLocalMap.at( globalID[g] ) );
      }
    }

    // Load nodal reference fields into reference arrays to perform global sync
    arrayView1d< localIndex > indexInGlobalCZArrayView = indexInGlobalCZArray;
    arrayView1d< localIndex > localNodalIndexView = localNodalIndex;
    forAll< serialPolicy >( indexInGlobalCZArray.size(), [=] GEOS_HOST ( localIndex const g ) 
    {
      localIndex czIndex = indexInGlobalCZArrayView[g];
      localIndex nodeIndex = localNodalIndexView[g];
      LvArray::tensorOps::copy< 3 >( referencePosition[czIndex], gridPosition[nodeIndex] );
      LvArray::tensorOps::copy< 3 >( referencePartitioningSurfaceNormal[czIndex], gridExplicitSurfaceNormal[nodeIndex] );
    } );

    // Sync initial reference nodal values with rank 0
    int const numRanks = MpiWrapper::commSize();
    int rank = MpiWrapper::commRank( MPI_COMM_GEOS );

    array1d< MPI_Request > mpiRequestIndices( numRanks );
    array1d< MPI_Request > mpiRequestPositions( numRanks );
    array1d< MPI_Request > mpiRequestPartitioningSurfaceNormals( numRanks );

    array1d< MPI_Status > mpiStatusIndices( numRanks );
    array1d< MPI_Status > mpiStatusPositions( numRanks );
    array1d< MPI_Status > mpiStatusPartitioningSurfaceNormals( numRanks );

    array1d< int > nodalIndices;
    if( rank != 0 )
    {
      // Send list of indices to overwrite nodal positions
      mpiRequestIndices[rank] = MPI_REQUEST_NULL;
      MpiWrapper::iSend( indexInGlobalCZArray.data(),
                         indexInGlobalCZArray.size(),
                         0,
                         0,
                         MPI_COMM_GEOS,
                         &mpiRequestIndices[rank] );

      // Send nodal positions
      mpiRequestPositions[rank] = MPI_REQUEST_NULL;
      MpiWrapper::iSend( referencePosition.data(),
                         referencePosition.size(),
                         0,
                         1,
                         MPI_COMM_GEOS,
                         &mpiRequestPositions[rank] );

      // Send nodal positions
      mpiRequestPartitioningSurfaceNormals[rank] = MPI_REQUEST_NULL;
      MpiWrapper::iSend( referencePartitioningSurfaceNormal.data(),
                         referencePartitioningSurfaceNormal.size(),
                         0,
                         1,
                         MPI_COMM_GEOS,
                         &mpiRequestPartitioningSurfaceNormals[rank] );
    }
    else
    {
      for( int r = 1; r < numRanks; r++ )
      {
        MpiWrapper::recv( nodalIndices,
                          r,
                          0,
                          MPI_COMM_GEOS,
                          &mpiStatusIndices[r] );

        array2d< real64 > nodalReferencePositions( numCohesiveNodes, 3 );
        mpiRequestPositions[r] = MPI_REQUEST_NULL;
        MpiWrapper::iRecv( nodalReferencePositions.data(),
                           nodalReferencePositions.size(),
                           r,
                           1,
                           MPI_COMM_GEOS,
                           &mpiRequestPositions[r] );

        MpiWrapper::wait( &mpiRequestPositions[r], &mpiStatusPositions[r] );

        // Send nodal surface normals
        array2d< real64 > nodalReferencePartitioningSurfaceNormals( numCohesiveNodes, 3 );
        mpiRequestPartitioningSurfaceNormals[r] = MPI_REQUEST_NULL;
        MpiWrapper::iRecv( nodalReferencePartitioningSurfaceNormals.data(),
                           nodalReferencePartitioningSurfaceNormals.size(),
                           r,
                           1,
                           MPI_COMM_GEOS,
                           &mpiRequestPartitioningSurfaceNormals[r] );

        MpiWrapper::wait( &mpiRequestPartitioningSurfaceNormals[r], &mpiStatusPartitioningSurfaceNormals[r] );

        // Combine grid positions
        for( int g = 0; g < nodalIndices.size(); ++g )
        {
          localIndex n = nodalIndices[g];
          LvArray::tensorOps::copy< 3 >( referencePosition[n], nodalReferencePositions[n] );
          // TODO: Check if we should be taking the max here for the partitioning surface normal?
          LvArray::tensorOps::copy< 3 >( referencePartitioningSurfaceNormal[n], nodalReferencePartitioningSurfaceNormals[n] );
        }
      }
    }

    // Wait for everything to complete
    MpiWrapper::barrier();

    // Scatter complete reference values back to other ranks
    MpiWrapper::bcast( referencePosition.data(),
                       referencePosition.size(),
                       0,
                       MPI_COMM_GEOS );

    MpiWrapper::bcast( referencePartitioningSurfaceNormal.data(),
                       referencePartitioningSurfaceNormal.size(),
                       0,
                       MPI_COMM_GEOS );
  } );

  // // Re-zero grid mass
  // // This may or may not be inside the above lambda for cz regions
  // forAll< serialPolicy >( numNodes, [=] GEOS_HOST ( localIndex const g )
  // {
  //   globalIndex const mappedNode = localToGlobalMap[ g ];

  //   for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
  //   {
  //     gridMass[g][fieldIndex] = 0.0;
  //   }
  // } );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CohesiveZoneManager, string const &, Group * const )
}
