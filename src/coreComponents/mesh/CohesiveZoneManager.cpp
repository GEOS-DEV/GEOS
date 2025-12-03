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


void CohesiveZoneManager::initializeReferenceConfiguration( int const numDims,
                                                            real64 const smallMass,
                                                            int const numVelocityFields,
                                                            stdVector< array2d< localIndex > > const & m_mappedNodes,
                                                            stdVector< array2d< real64 > > const & m_shapeFunctionValues, 
                                                            NodeManager & nodeManager,
                                                            ParticleManager & particleManager )
{
  GEOS_UNUSED_VAR( particleManager );

  // Grid fields
  arrayView1d< globalIndex > localToGlobalMap = nodeManager.localToGlobalMap();
  auto & nodeGlobalToLocalMap = nodeManager.globalToLocalMap();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const gridPosition = nodeManager.referencePosition();
  arrayView2d< real64 const > const gridExplicitSurfaceNormal = nodeManager.getReference< array2d< real64 > >( SolidMechanicsMPM::viewKeyStruct::gridExplicitSurfaceNormalString() );

  this->forCohesiveZoneRegions< CohesiveZoneRegion >( [&]( CohesiveZoneRegion & czRegion )
  {
    SortedArrayView< globalIndex const > globalID = czRegion.getGlobalID();
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
      // if( nodeGlobalToLocalMap.contains( globalID[g] ) ) // Only available with C++ 20
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

    // Initialize state variables
    arrayView1d< real64 > czDamage = czRegion.getDamage();
    arrayView1d< real64 > czMaxNormalDisplacement = czRegion.getMaxNormalDisplacement();
    arrayView1d< real64 > czMaxTangentialDisplacement = czRegion.getMaxTangentialDisplacement();

    forAll< serialPolicy >( numCohesiveNodes, [=]( localIndex const & g ){
      czDamage[g] = 0.0;
      czMaxNormalDisplacement[g] = 0.0;
      czMaxTangentialDisplacement[g] = 0.0;
    } );
  } );

  // Temporarily hardcoded for now, but ideally should look up the node being mapped to and using the velocity field figure out which czRegion to map to
  CohesiveZoneRegionBase & czRegion = this->getRegion( "cz1" );
  int const numCohesiveNodes = czRegion.size();
  SortedArrayView< globalIndex const > globalID = czRegion.getGlobalID();
  arrayView2d< real64 const > const czReferencePosition = czRegion.getReferencePosition();

  // Flag cohesive particles and store reference mapping data
  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< int const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView3d< real64 const > const particleDeformationGradient = subRegion.getField< fields::mpm::particleDeformationGradient >();

    arrayView1d< int > particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();
    arrayView2d< globalIndex > particleReferenceMappedNodes = subRegion.getField< fields::mpm::particleReferenceMappedNodes >();
    arrayView2d< real64 > particleReferenceShapeFunctionValues = subRegion.getField< fields::mpm::particleReferenceShapeFunctionValues >();

    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();

    arrayView2d< real64 > const particleCohesiveReferenceSurfaceNormal = subRegion.getField< fields::mpm::particleCohesiveReferenceSurfaceNormal >();
    arrayView2d< real64 > const particleCohesiveReferenceSurfacePosition = subRegion.getField< fields::mpm::particleCohesiveReferenceSurfacePosition >();
    arrayView3d< real64 > const particleCohesiveReferenceDeformationGradient = subRegion.getField< fields::mpm::particleCohesiveReferenceDeformationGradient >();

    // Get views to mapping arrays
    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
    arrayView2d< real64 const > const shapeFunctionValues = m_shapeFunctionValues[subRegionIndex];

    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(), [=] GEOS_HOST_DEVICE ( localIndex const pp )
    {
      localIndex const p = activeParticleIndices[pp];

      // Store reference deformation gradient of particles to update surface normals and positions
      LvArray::tensorOps::copy< 3, 3 >( particleCohesiveReferenceDeformationGradient[p], particleDeformationGradient[p] );
      LvArray::tensorOps::copy< 3 >( particleCohesiveReferenceSurfaceNormal[p], particleSurfaceNormal[p] );
      LvArray::tensorOps::copy< 3 >( particleCohesiveReferenceSurfacePosition[p], particleSurfacePosition[p] );

      for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
      {
        globalIndex const mappedNode = localToGlobalMap[ mappedNodes[pp][g] ];

        if( globalID.contains( mappedNode ) )
        {
          if( static_cast< SolidMechanicsMPM::SurfaceFlag >( particleSurfaceFlag[p] ) == SolidMechanicsMPM::SurfaceFlag::Cohesive )
          {
            particleCohesiveZoneFlag[p] = 1;
          }

          // Store particle shape function data for future cohesive computations
          for( int gg = 0; gg < numberOfVerticesPerParticle * 8; gg++ )
          {
            particleReferenceMappedNodes[p][gg] = localToGlobalMap[ mappedNodes[pp][gg] ];
            particleReferenceShapeFunctionValues[p][gg] = shapeFunctionValues[pp][gg];
          }

          // If at least one corner is involved no need to look further
          break;
        }
      }
    } );
    ++subRegionIndex;
  } );

  // Now we need to compute the grid area at cohesive initialization
  // We map the surface position of each particle ( vector from particle center to interface surface ), this is also the particle surface
  // normal direction
  array2d< real64 > tempGridMassLocal( numCohesiveNodes, numVelocityFields );
  array1d< real64 > tempGridVolumeLocal( numCohesiveNodes );
  array3d< real64 > tempGridParticleSurfaceNormalLocal( numCohesiveNodes, numVelocityFields, 3 );
  array3d< real64 > tempGridSurfacePositionLocal( numCohesiveNodes, numVelocityFields, 3 );

  // Initialize temporary grid fields to zero
  for( int g = 0; g < numCohesiveNodes; ++g )
  {
    tempGridVolumeLocal[g] = 0.0;
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      tempGridMassLocal[g][fieldIndex] = 0.0;
      for( int i = 0; i < numDims; ++i )
      {
        tempGridParticleSurfaceNormalLocal[g][fieldIndex][i] = 0.0;
        tempGridSurfacePositionLocal[g][fieldIndex][i] = 0.0;
      }
    }
  }

  // Map relevant quantities to cz grid for cz node initialization
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
    arrayView2d< real64 const > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();

    arrayView1d< int const > const particleCohesiveZoneFlag = subRegion.getField< fields::mpm::particleCohesiveZoneFlag >();
    arrayView2d< globalIndex const > const particleReferenceMappedNodes = subRegion.getField< fields::mpm::particleReferenceMappedNodes >();
    arrayView2d< real64 const > const particleReferenceShapeFunctionValues = subRegion.getField< fields::mpm::particleReferenceShapeFunctionValues >();
    arrayView2d< int const > const particleCohesiveFieldMapping = subRegion.getField< fields::mpm::particleCohesiveFieldMapping >();

    int const numberOfVerticesPerParticle = subRegion.numberOfVerticesPerParticle();

    auto globalToLocalMap = nodeManager.globalToLocalMap();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();
    forAll< serialPolicy >( activeParticleIndices.size(),
                            [=,
                             &tempGridMassLocal,
                             &tempGridVolumeLocal,
                             &tempGridParticleSurfaceNormalLocal,
                             &tempGridSurfacePositionLocal] GEOS_HOST ( localIndex const pp )
      {
        localIndex const p = activeParticleIndices[pp];

        // Map particle displacement to grid
        for( int g = 0; g < 8 * numberOfVerticesPerParticle; ++g )
        {
          globalIndex const mappedNode = particleReferenceMappedNodes[p][g];

          if( globalID.contains( mappedNode ) )
          {
            real64 shapeFunctionValue = particleReferenceShapeFunctionValues[p][g];

            // CC: TODO must be a better way to find index in temp arrays
            localIndex nodeIndex = 0;
            for( int n = 0; n < numCohesiveNodes; ++n )
            {
              if( globalID[n] == mappedNode )
              {
                nodeIndex = n;
                break;
              }
            }

            localIndex const fieldIndex = particleCohesiveFieldMapping[p][g];

            if( particleCohesiveZoneFlag[p] == 1 )
            {
              tempGridMassLocal[nodeIndex][fieldIndex] += particleMass[p] * shapeFunctionValue;

              for( int i  = 0; i < numDims; ++i )
              {
                tempGridParticleSurfaceNormalLocal[nodeIndex][fieldIndex][i] += particleMass[p] * particleSurfaceNormal[p][i] * shapeFunctionValue;
                tempGridSurfacePositionLocal[nodeIndex][fieldIndex][i] += particleMass[p] *
                                                                          ( particlePosition[p][i] - czReferencePosition[nodeIndex][i] + particleSurfacePosition[p][i] ) *
                                                                          shapeFunctionValue;
              }

            }
            tempGridVolumeLocal[nodeIndex] += particleVolume[p] * shapeFunctionValue;
          }
        }
      } );
  } );

  // Sync cz nodal fields
  array2d< real64 > tempGridMassGlobal( numCohesiveNodes, numVelocityFields );
  array1d< real64 > tempGridVolumeGlobal( numCohesiveNodes );
  array3d< real64 > tempGridParticleSurfaceNormalGlobal( numCohesiveNodes, numVelocityFields, 3 );
  array3d< real64 > tempGridSurfacePositionGlobal( numCohesiveNodes, numVelocityFields, 3 );

  MpiWrapper::allReduce( tempGridMassLocal,
                         tempGridMassGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridVolumeLocal,
                         tempGridVolumeGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridParticleSurfaceNormalLocal,
                         tempGridParticleSurfaceNormalGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  MpiWrapper::allReduce( tempGridSurfacePositionLocal,
                         tempGridSurfacePositionGlobal,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  // Normalize cz nodal fields
  forAll< serialPolicy >( numCohesiveNodes, [=, &tempGridParticleSurfaceNormalGlobal, &tempGridSurfacePositionGlobal] GEOS_HOST ( localIndex const g )
  {
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      if( tempGridMassGlobal[g][fieldIndex] > smallMass )
      {
        for( int i = 0; i < numDims; ++i ) // Does this need to be using numDims? 
        {
          tempGridParticleSurfaceNormalGlobal[g][fieldIndex][i] /= tempGridMassGlobal[g][fieldIndex];
          tempGridSurfacePositionGlobal[g][fieldIndex][i] /= tempGridMassGlobal[g][fieldIndex];
        }
      }
    }
  } );


}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CohesiveZoneManager, string const &, Group * const )
}
