/**
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

#include "VTKSuperCellPartitioning.hpp"

// GEOS mesh includes
#include "mesh/generators/VTKMeshGeneratorTools.hpp"


// GEOS ParMETIS/Scotch (conditional)
#ifdef GEOS_USE_PARMETIS
#include "mesh/generators/ParMETISInterface.hpp"
#endif

// GEOS common includes
#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"
#include "common/format/StringUtilities.hpp"

// LvArray
#include "LvArray/src/ArrayOfArrays.hpp"

// VTK includes
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkCellData.h>
#include <vtkPartitionedDataSet.h>
#include <vtkExtractCells.h>

#include <unordered_set>


namespace geos
{

namespace vtk
{

// =============================================================================================
// SECTION 1: SUPER-CELL IDENTIFICATION (rank 0 only)
// =============================================================================================
pmet_idx_t computeSuperCellWeight( vtkIdType numCells, integer fractureWeight )
{
  if( numCells > 1 )
  {
    return numCells + fractureWeight;  // Fracture super-cell
  }
  else
  {
    return numCells;  // Regular cell (weight = 1)
  }
}

SuperCellInfo tagCellsWithSuperCellIds(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  stdMap< string, ArrayOfArrays< vtkIdType, int64_t > > const & fractureNeighbors,
  integer fractureWeight )
{
  GEOS_LOG_RANK_0( "TAGGING 3D CELLS WITH SUPER-CELL IDs" );

  vtkIdType const numLocalCells = cells3D->GetNumberOfCells();

  vtkIdTypeArray * globalIds =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetGlobalIds() );

  GEOS_ERROR_IF( !globalIds, "3D mesh missing global IDs" );

// -----------------------------------------------------------------------
// Step 1: Build 3D cell connectivity graph via fractures
// -----------------------------------------------------------------------
  std::map< vtkIdType, std::set< vtkIdType > > fractureGraph;
  vtkIdType totalFracturePairs = 0;

  for( auto const & [fractureName, neighbors] : fractureNeighbors )
  {
    vtkIdType numFracCells = neighbors.size();

    for( vtkIdType i = 0; i < numFracCells; ++i )
    {
      auto neighborList = neighbors[i];

      if( neighborList.size() >= 2 )
      {
        vtkIdType gidA = neighborList[0];
        vtkIdType gidB = neighborList[1];

        fractureGraph[gidA].insert( gidB );
        fractureGraph[gidB].insert( gidA );

        totalFracturePairs++;
      }
    }

    GEOS_LOG_RANK_0( GEOS_FMT( "Fracture '{}': processed {} fracture elements",
                               fractureName, numFracCells ) );
  }

  GEOS_LOG_RANK_0( GEOS_FMT( "Built fracture graph with {} 3D cells having fracture connections, {} fracture pairs",
                             fractureGraph.size(), totalFracturePairs ) );

  // -----------------------------------------------------------------------
  // Step 2: Find the maximum global ID from ALL cells (not just fracture cells)
  // -----------------------------------------------------------------------
  vtkIdType maxGlobalId = 0;

  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType gid = globalIds->GetValue( i );
    maxGlobalId = std::max( maxGlobalId, gid );
  }

  // -----------------------------------------------------------------------
  // Step 3: Find connected components using DFS
  // -----------------------------------------------------------------------
  std::map< vtkIdType, vtkIdType > cellToSuperCell;
  std::set< vtkIdType > visited;

  // Start super-cell IDs AFTER the maximum global ID to avoid collisions
  vtkIdType nextSuperCellId = maxGlobalId + 1;

  std::function< void(vtkIdType, vtkIdType, std::vector< vtkIdType > &) > dfs =
    [&]( vtkIdType cell, vtkIdType superCellId, std::vector< vtkIdType > & component )
  {
    if( visited.count( cell ) )
      return;
    // Mark this cell as visited
    visited.insert( cell );
    // Assign this cell to the current super-cell
    cellToSuperCell[cell] = superCellId;
    component.push_back( cell );

    // Recursively visit all neighbors
    if( fractureGraph.count( cell ) )
    {
      for( vtkIdType neighbor : fractureGraph.at( cell ) )
      {
        dfs( neighbor, superCellId, component );
      }
    }
  };

  std::map< vtkIdType, std::vector< vtkIdType > > superCellComponents;

  for( auto const & [cell, neighbors] : fractureGraph )
  {
    if( !visited.count( cell ) )
    {
      std::vector< vtkIdType > component;
      dfs( cell, nextSuperCellId, component );

      superCellComponents[nextSuperCellId] = component;
      nextSuperCellId++;
    }
  }

  // -----------------------------------------------------------------------
  // Step 4: Build SuperCellInfo
  // -----------------------------------------------------------------------
  SuperCellInfo info;

  vtkIdType numCellsInSuperCells = 0;
  vtkIdType numSuperCellsCreated = 0;
  vtkIdType largestSuperCellSize = 0;

  for( auto const & [scId, members] : superCellComponents )
  {
    info.superCellToOriginalCells[scId] = members;
    info.vertexWeights[scId] = computeSuperCellWeight( members.size(), fractureWeight );

    if( members.size() > 1 )
    {
      info.atomicSuperCells.insert( scId );
    }

    numCellsInSuperCells += members.size();
    numSuperCellsCreated++;
    largestSuperCellSize = std::max( largestSuperCellSize,
                                     static_cast< vtkIdType >( members.size() ) );
  }

  // -----------------------------------------------------------------------
  // Step 5: Tag cells with SuperCellIds
  // -----------------------------------------------------------------------
  vtkNew< vtkIdTypeArray > superCellIdArray;
  superCellIdArray->SetName( "SuperCellId" );
  superCellIdArray->SetNumberOfTuples( numLocalCells );

  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType globalId = globalIds->GetValue( i );

    if( cellToSuperCell.count( globalId ) )
    {
      // Cell is part of a fracture-connected super-cell
      superCellIdArray->SetValue( i, cellToSuperCell.at( globalId ) );
    }
    else
    {
      // Regular cell - use its own global ID as super-cell ID
      // This is safe because super-cell IDs start at maxGlobalId + 1
      superCellIdArray->SetValue( i, globalId );
    }
  }

  cells3D->GetCellData()->AddArray( superCellIdArray );

  // -----------------------------------------------------------------------
  // Step 6: Report statistics
  // -----------------------------------------------------------------------
  vtkIdType numRegularCells = numLocalCells - numCellsInSuperCells;
  vtkIdType totalSuperCells = numSuperCellsCreated + numRegularCells;
  vtkIdType cellReduction = numLocalCells - totalSuperCells;

  GEOS_LOG_RANK_0( "SUPER-CELL TAGGING SUMMARY" );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Total 3D cells:                   {:8}", numLocalCells ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Cells in super-cells:             {:8}", numCellsInSuperCells ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Regular cells (no fractures):     {:8}", numRegularCells ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Number of super-cells created:    {:8}", numSuperCellsCreated ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Total super-cells (incl regular): {:8}", totalSuperCells ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Cell reduction:                   {:8}", cellReduction ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Largest super-cell size:          {:8} cells", largestSuperCellSize ) );

  return info;
}


// =============================================================================================
// SECTION 2: SUPER-CELL RECONSTRUCTION (after redistribution)
// =============================================================================================
SuperCellInfo reconstructSuperCellInfo( vtkSmartPointer< vtkUnstructuredGrid > mesh, integer fractureWeight )
{
  SuperCellInfo info;

  vtkIdTypeArray * superCellIdArray =
    vtkIdTypeArray::SafeDownCast( mesh->GetCellData()->GetArray( "SuperCellId" ) );

  if( !superCellIdArray )
  {
    return info;  // No super-cells, return empty info
  }

  vtkIdTypeArray * globalIds =
    vtkIdTypeArray::SafeDownCast( mesh->GetCellData()->GetGlobalIds() );

  GEOS_ERROR_IF( !globalIds, "Mesh missing global IDs" );

  // Build map: super-cell ID -> vector of cell global IDs
  std::map< vtkIdType, std::vector< vtkIdType > > localSuperCells;

  for( vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i )
  {
    vtkIdType scId = superCellIdArray->GetValue( i );
    vtkIdType globalId = globalIds->GetValue( i );
    localSuperCells[scId].push_back( globalId );
  }

// Populate SuperCellInfo
  for( auto const & [scId, cells] : localSuperCells )
  {
    info.superCellToOriginalCells[scId] = cells;
    info.vertexWeights[scId] = computeSuperCellWeight( cells.size(), fractureWeight );

    if( cells.size() > 1 )
    {
      info.atomicSuperCells.insert( scId );
    }
  }

  return info;
}


// =============================================================================================
// SECTION 3: INITIAL REDISTRIBUTION (preserving super-cells)
// =============================================================================================
vtkSmartPointer< vtkDataSet >
redistributeBySuperCellBlocks( vtkSmartPointer< vtkUnstructuredGrid > cells3D,
                               MPI_Comm comm )
{
  GEOS_MARK_FUNCTION;

  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );

  vtkSmartPointer< vtkPartitionedDataSet > partitionedMesh =
    vtkSmartPointer< vtkPartitionedDataSet >::New();

  if( rank == 0 )
  {
    vtkIdTypeArray * superCellIdArray =
      vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetArray( "SuperCellId" ) );

    GEOS_ERROR_IF( !superCellIdArray, "SuperCellId array not found" );

    vtkIdType numCells = cells3D->GetNumberOfCells();

    // Build super-cell metadata
    std::map< vtkIdType, std::vector< vtkIdType > > superCellToLocalCells;

    for( vtkIdType i = 0; i < numCells; ++i )
    {
      vtkIdType scId = superCellIdArray->GetValue( i );
      superCellToLocalCells[scId].push_back( i );
    }
    vtkIdType numSuperCells = superCellToLocalCells.size();

    // -----------------------------------------------------------------------
    // Step 1: Pre-compute ALL cell centroids
    // -----------------------------------------------------------------------
    std::vector< std::array< double, 3 > > allCellCentroids( numCells );
    vtkPoints * meshPoints = cells3D->GetPoints();

    for( vtkIdType cellIdx = 0; cellIdx < numCells; ++cellIdx )
    {
      vtkIdType npts;
      const vtkIdType * pts;
      cells3D->GetCellPoints( cellIdx, npts, pts );

      double centroid[3] = {0.0, 0.0, 0.0};
      for( vtkIdType i = 0; i < npts; ++i )
      {
        double pt[3];
        meshPoints->GetPoint( pts[i], pt );
        centroid[0] += pt[0];
        centroid[1] += pt[1];
        centroid[2] += pt[2];
      }

      allCellCentroids[cellIdx][0] = centroid[0] / npts;
      allCellCentroids[cellIdx][1] = centroid[1] / npts;
      allCellCentroids[cellIdx][2] = centroid[2] / npts;
    }

    // -----------------------------------------------------------------------
    // Step 2: Compute super-cell centroids from pre-computed values
    // -----------------------------------------------------------------------
    struct SuperCellPartitionInfo
    {
      vtkIdType scId;
      std::array< double, 3 > centroid;
      std::vector< vtkIdType > cellIndices;
    };

    std::vector< SuperCellPartitionInfo > superCells;
    superCells.reserve( numSuperCells );

    for( auto const & [scId, cellIndices] : superCellToLocalCells )
    {
      std::array< double, 3 > centroid = {0.0, 0.0, 0.0};

      // Fast: just sum pre-computed cell centroids
      for( vtkIdType cellIdx : cellIndices )
      {
        centroid[0] += allCellCentroids[cellIdx][0];
        centroid[1] += allCellCentroids[cellIdx][1];
        centroid[2] += allCellCentroids[cellIdx][2];
      }

      // Average across cells in super-cell
      centroid[0] /= cellIndices.size();
      centroid[1] /= cellIndices.size();
      centroid[2] /= cellIndices.size();

      superCells.push_back( SuperCellPartitionInfo{ scId, centroid, cellIndices } );
    }

    // Free cell centroids (no longer needed)
    allCellCentroids.clear();
    allCellCentroids.shrink_to_fit();

    // -----------------------------------------------------------------------
    // Step 3: Find bounding box for Morton code normalization
    // -----------------------------------------------------------------------
    double minCoord[3] = {std::numeric_limits< double >::max(),
                          std::numeric_limits< double >::max(),
                          std::numeric_limits< double >::max()};
    double maxCoord[3] = {std::numeric_limits< double >::lowest(),
                          std::numeric_limits< double >::lowest(),
                          std::numeric_limits< double >::lowest()};

    for( auto const & sc : superCells )
    {
      for( int d = 0; d < 3; ++d )
      {
        minCoord[d] = std::min( minCoord[d], sc.centroid[d] );
        maxCoord[d] = std::max( maxCoord[d], sc.centroid[d] );
      }
    }

    GEOS_LOG_RANK_0( GEOS_FMT( "Bounding box: X=[{:.3f}, {:.3f}], Y=[{:.3f}, {:.3f}], Z=[{:.3f}, {:.3f}]",
                               minCoord[0], maxCoord[0], minCoord[1], maxCoord[1], minCoord[2], maxCoord[2] ));

    // -----------------------------------------------------------------------
    // Step 4: Sort by Morton algorith to ensure spatial locality
    // -----------------------------------------------------------------------
    auto computeMorton = []( std::array< double, 3 > const & centroid,
                             double bounds_min[3],
                             double bounds_max[3] ) -> uint64_t
    {
      // Normalize coordinates to [0, 1]
      auto normalize = [&]( double val, int dim ) -> uint32_t
      {
        double range = bounds_max[dim] - bounds_min[dim];
        if( range < 1e-10 )
          return 0;
        double norm = (val - bounds_min[dim]) / range;
        norm = std::max( 0.0, std::min( 1.0, norm ) );
        return static_cast< uint32_t >( norm * ((1u << 21) - 1) ); // 21 bits per dimension
      };

      uint32_t x = normalize( centroid[0], 0 );
      uint32_t y = normalize( centroid[1], 1 );
      uint32_t z = normalize( centroid[2], 2 );

      // Interleave bits (Morton encoding)
      uint64_t code = 0;
      for( int i = 0; i < 21; ++i )
      {
        code |= ((x & (1u << i)) ? (1ull << (3*i)) : 0);
        code |= ((y & (1u << i)) ? (1ull << (3*i + 1)) : 0);
        code |= ((z & (1u << i)) ? (1ull << (3*i + 2)) : 0);
      }
      return code;
    };


    // Sort by Morton
    std::sort( superCells.begin(), superCells.end(),
               [&]( SuperCellPartitionInfo const & a, SuperCellPartitionInfo const & b )
    {
      return computeMorton( a.centroid, minCoord, maxCoord ) <
      computeMorton( b.centroid, minCoord, maxCoord );
    } );

    GEOS_LOG_RANK_0( "Super-cells sorted by spatial locality (Z-order curve)" );

    // -----------------------------------------------------------------------
    // Step 5: Assign sorted super-cells to ranks in contiguous blocks
    // -----------------------------------------------------------------------
    array1d< int64_t > cellPartitions( numCells );
    std::vector< vtkIdType > cellsPerRank( numRanks, 0 );

    vtkIdType superCellsPerRank = (numSuperCells + numRanks - 1) / numRanks;

    for( vtkIdType scIdx = 0; scIdx < numSuperCells; ++scIdx )
    {
      int targetRank = std::min( static_cast< int >(scIdx / superCellsPerRank), numRanks - 1 );

      // All cells in this super-cell go to the same rank
      for( vtkIdType cellIdx : superCells[scIdx].cellIndices )
      {
        cellPartitions[cellIdx] = targetRank;
        cellsPerRank[targetRank]++;
      }
    }

    // -----------------------------------------------------------------------
    // Step 6: Build partitions
    // -----------------------------------------------------------------------
    partitionedMesh->SetNumberOfPartitions( numRanks );

    for( int r = 0; r < numRanks; ++r )
    {
      vtkNew< vtkIdList > cellsForRank;

      for( vtkIdType i = 0; i < numCells; ++i )
      {
        if( cellPartitions[i] == r )
        {
          cellsForRank->InsertNextId( i );
        }
      }

      if( cellsForRank->GetNumberOfIds() == 0 )
      {
        vtkNew< vtkUnstructuredGrid > emptyPart;
        partitionedMesh->SetPartition( r, emptyPart );
        continue;
      }

      vtkNew< vtkExtractCells > extractor;
      extractor->SetInputData( cells3D );
      extractor->SetCellList( cellsForRank );
      extractor->Update();

      vtkSmartPointer< vtkUnstructuredGrid > partition =
        vtkUnstructuredGrid::SafeDownCast( extractor->GetOutput() );

      partitionedMesh->SetPartition( r, partition );
    }

  }
  else
  {
    // Other ranks: create empty partitioned dataset
    partitionedMesh->SetNumberOfPartitions( numRanks );
    for( int r = 0; r < numRanks; ++r )
    {
      vtkNew< vtkUnstructuredGrid > emptyPart;
      partitionedMesh->SetPartition( r, emptyPart );
    }
  }


  // -----------------------------------------------------------------------
  // Step 7: Redistribute using VTK
  // -----------------------------------------------------------------------
  vtkSmartPointer< vtkDataSet > result = vtk::redistribute( *partitionedMesh, comm );

  partitionedMesh = nullptr;
  if( rank == 0 )
  {
    cells3D = nullptr;
  }

  vtkIdType localCells = result->GetNumberOfCells();

  // Verify SuperCellId array survived redistribution
  if( localCells > 0 )
  {
    vtkIdTypeArray * resultSuperCellIdArray =
      vtkIdTypeArray::SafeDownCast( result->GetCellData()->GetArray( "SuperCellId" ) );

    GEOS_ERROR_IF( !resultSuperCellIdArray,
                   GEOS_FMT( "Rank {}: SuperCellId array lost during redistribution", rank ) );
  }

  return result;
}


// =============================================================================================
// SECTION 4: SUPER-CELL GRAPH BUILDING
// =============================================================================================
std::pair< ArrayOfArrays< pmet_idx_t, pmet_idx_t >, array1d< pmet_idx_t > >
buildSuperCellGraph(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const & baseGraph,
  arrayView1d< pmet_idx_t const > const & baseElemDist,
  SuperCellInfo const & info,
  MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );

  vtkIdTypeArray * superCellIdArray =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetArray( "SuperCellId" ) );
  GEOS_ERROR_IF( superCellIdArray == nullptr, "SuperCellId array not found" );

  vtkIdTypeArray * cellGlobalIds =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetGlobalIds() );
  GEOS_ERROR_IF( cellGlobalIds == nullptr, "Cell global IDs not found" );

  vtkIdType numLocalCells = cells3D->GetNumberOfCells();

  std::map< vtkIdType, std::vector< vtkIdType > > superCellToLocalCells;
  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType superCellId = superCellIdArray->GetValue( i );
    superCellToLocalCells[superCellId].push_back( i );
  }

  localIndex numLocalSuperCells = superCellToLocalCells.size();

  // Build super-cell distribution
  array1d< pmet_idx_t > superElemDist( numRanks + 1 );
  pmet_idx_t localSuperCellCount = numLocalSuperCells;
  MpiWrapper::allgather( &localSuperCellCount, 1, superElemDist.data(), 1, comm );

  pmet_idx_t temp = superElemDist[0];
  superElemDist[0] = 0;
  for( int r = 1; r <= numRanks; ++r )
  {
    pmet_idx_t next = superElemDist[r];
    superElemDist[r] = superElemDist[r-1] + temp;
    temp = next;
  }

  pmet_idx_t myGlobalStart = superElemDist[rank];
  pmet_idx_t myGlobalEnd = superElemDist[rank + 1];

  std::vector< vtkIdType > orderedSuperCellIds;
  orderedSuperCellIds.reserve( numLocalSuperCells );
  for( auto const & [superCellId, localCells] : superCellToLocalCells )
  {
    orderedSuperCellIds.push_back( superCellId );
  }
  std::sort( orderedSuperCellIds.begin(), orderedSuperCellIds.end() );

  std::map< vtkIdType, pmet_idx_t > superCellIdToGlobalIdx;
  for( pmet_idx_t i = 0; i < numLocalSuperCells; ++i )
  {
    superCellIdToGlobalIdx[orderedSuperCellIds[i]] = myGlobalStart + i;
  }

  // -----------------------------------------------------------------------
  // Step 1: Build mappings for exchange
  // -----------------------------------------------------------------------
  // Local maps
  std::unordered_map< pmet_idx_t, vtkIdType > myGlobalCellToSuperCell;
  myGlobalCellToSuperCell.reserve( numLocalCells );

  std::unordered_map< pmet_idx_t, vtkIdType > myParmetisToVtk;
  myParmetisToVtk.reserve( numLocalCells );

  pmet_idx_t myParmetisStart = baseElemDist[rank];
  pmet_idx_t myParmetisEnd = baseElemDist[rank + 1];

  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType vtkGlobalId = cellGlobalIds->GetValue( i );
    vtkIdType superCellId = superCellIdArray->GetValue( i );

    myGlobalCellToSuperCell[vtkGlobalId] = superCellId;
    myParmetisToVtk[myParmetisStart + i] = vtkGlobalId;
  }

  // -----------------------------------------------------------------------
  // Step 2: Identify unique ghost ParMETIS indices
  // -----------------------------------------------------------------------
  std::set< pmet_idx_t > ghostParmetisSet;

  for( localIndex i = 0; i < baseGraph.size(); ++i )
  {
    auto neighbors = baseGraph[i];
    for( localIndex j = 0; j < neighbors.size(); ++j )
    {
      pmet_idx_t nbrIdx = neighbors[j];
      if( nbrIdx < myParmetisStart || nbrIdx >= myParmetisEnd )
      {
        ghostParmetisSet.insert( nbrIdx );
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 3: Exchange ghost mappings
  // -----------------------------------------------------------------------
  // Prepare local data to send
  std::vector< pmet_idx_t > localParmetisIndices;
  std::vector< vtkIdType > localVtkIds;
  std::vector< vtkIdType > localSuperCellIds;

  localParmetisIndices.reserve( numLocalCells );
  localVtkIds.reserve( numLocalCells );
  localSuperCellIds.reserve( numLocalCells );

  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    localParmetisIndices.push_back( myParmetisStart + i );
    localVtkIds.push_back( cellGlobalIds->GetValue( i ) );
    localSuperCellIds.push_back( superCellIdArray->GetValue( i ) );
  }

  // Gather sizes
  int localCount = static_cast< int >( localParmetisIndices.size() );
  array1d< int > allCounts( numRanks );
  MpiWrapper::allgather( &localCount, 1, allCounts.data(), 1, comm );

  array1d< int > displs( numRanks + 1 );
  displs[0] = 0;
  for( int r = 0; r < numRanks; ++r )
  {
    displs[r + 1] = displs[r] + allCounts[r];
  }

  int totalMappings = displs[numRanks];

  // Allocate and gather
  array1d< pmet_idx_t > allParmetisIndices( totalMappings );
  array1d< vtkIdType > allVtkIds( totalMappings );
  array1d< vtkIdType > allSuperCellIds( totalMappings );

  MpiWrapper::allgatherv( localParmetisIndices.data(), localCount,
                          allParmetisIndices.data(), allCounts.data(), displs.data(), comm );
  MpiWrapper::allgatherv( localVtkIds.data(), localCount,
                          allVtkIds.data(), allCounts.data(), displs.data(), comm );
  MpiWrapper::allgatherv( localSuperCellIds.data(), localCount,
                          allSuperCellIds.data(), allCounts.data(), displs.data(), comm );

  // Build lookup tables from gathered data (only for ghosts + local)
  std::unordered_map< pmet_idx_t, vtkIdType > parmetisToVtk;
  std::unordered_map< vtkIdType, vtkIdType > vtkToSuperCell;

  parmetisToVtk.reserve( ghostParmetisSet.size() + numLocalCells );
  vtkToSuperCell.reserve( ghostParmetisSet.size() + numLocalCells );

  for( int i = 0; i < totalMappings; ++i )
  {
    pmet_idx_t parmetisIdx = allParmetisIndices[i];

    // Only store if it's local or a ghost we need
    if( (parmetisIdx >= myParmetisStart && parmetisIdx < myParmetisEnd) ||
        ghostParmetisSet.count( parmetisIdx ) > 0 )
    {
      vtkIdType vtkId = allVtkIds[i];
      vtkIdType superCellId = allSuperCellIds[i];

      parmetisToVtk[parmetisIdx] = vtkId;
      vtkToSuperCell[vtkId] = superCellId;
    }
  }

  // Cleanup
  allParmetisIndices.clear();
  allVtkIds.clear();
  allSuperCellIds.clear();

  // -----------------------------------------------------------------------
  // Step 4: Exchange super-cell global indices
  // -----------------------------------------------------------------------
  std::vector< vtkIdType > sendSCIds;
  std::vector< pmet_idx_t > sendSCGlobalIndices;

  for( auto const & [scId, gIdx] : superCellIdToGlobalIdx )
  {
    sendSCIds.push_back( scId );
    sendSCGlobalIndices.push_back( gIdx );
  }

  int localSCCount = static_cast< int >( sendSCIds.size() );
  array1d< int > allSCCounts( numRanks );
  MpiWrapper::allgather( &localSCCount, 1, allSCCounts.data(), 1, comm );

  array1d< int > scDispls( numRanks + 1 );
  scDispls[0] = 0;
  for( int r = 0; r < numRanks; ++r )
  {
    scDispls[r + 1] = scDispls[r] + allSCCounts[r];
  }

  int totalSCMappings = scDispls[numRanks];

  array1d< vtkIdType > allSCIds( totalSCMappings );
  array1d< pmet_idx_t > allSCGlobalIndices( totalSCMappings );

  MpiWrapper::allgatherv( sendSCIds.data(), localSCCount,
                          allSCIds.data(), allSCCounts.data(), scDispls.data(), comm );
  MpiWrapper::allgatherv( sendSCGlobalIndices.data(), localSCCount,
                          allSCGlobalIndices.data(), allSCCounts.data(), scDispls.data(), comm );

  // Update local map with all super-cell mappings
  for( int i = 0; i < totalSCMappings; ++i )
  {
    superCellIdToGlobalIdx[allSCIds[i]] = allSCGlobalIndices[i];
  }

  // -----------------------------------------------------------------------
  // Step 5: Build super-cell graph edges
  // -----------------------------------------------------------------------
  array1d< pmet_idx_t > superVertexWeights( numLocalSuperCells );
  std::vector< std::set< pmet_idx_t > > neighborSets( numLocalSuperCells );

  for( localIndex localSuperIdx = 0; localSuperIdx < numLocalSuperCells; ++localSuperIdx )
  {
    vtkIdType superCellId = orderedSuperCellIds[localSuperIdx];
    auto const & localCells = superCellToLocalCells.at( superCellId );

    auto itWeight = info.vertexWeights.find( superCellId );
    superVertexWeights[localSuperIdx] = (itWeight != info.vertexWeights.end())
                                        ? itWeight->second : localCells.size();

    for( vtkIdType cellLocalIdx : localCells )
    {
      auto neighbors = baseGraph[cellLocalIdx];
      for( localIndex j = 0; j < neighbors.size(); ++j )
      {
        pmet_idx_t nbrParmetis = neighbors[j];

        auto itVtk = parmetisToVtk.find( nbrParmetis );
        if( itVtk == parmetisToVtk.end() )
          continue;

        vtkIdType nbrVtkId = itVtk->second;

        auto itSC = vtkToSuperCell.find( nbrVtkId );
        if( itSC == vtkToSuperCell.end() )
          continue;

        vtkIdType nbrSuperCellId = itSC->second;
        if( nbrSuperCellId == superCellId )
          continue;  // Skip self-loops

        auto itGlobal = superCellIdToGlobalIdx.find( nbrSuperCellId );
        if( itGlobal != superCellIdToGlobalIdx.end() )
        {
          neighborSets[localSuperIdx].insert( itGlobal->second );
        }
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 6: Symmetrize
  // -----------------------------------------------------------------------
  // Collect ALL edges to send (to any rank)
  std::vector< pmet_idx_t > myEdgesSrc;
  std::vector< pmet_idx_t > myEdgesDst;

  for( localIndex i = 0; i < numLocalSuperCells; ++i )
  {
    pmet_idx_t globalI = myGlobalStart + i;
    for( pmet_idx_t nbrIdx : neighborSets[i] )
    {
      if( nbrIdx < myGlobalStart || nbrIdx >= myGlobalEnd )
      {
        myEdgesSrc.push_back( nbrIdx );  // destination of reverse edge
        myEdgesDst.push_back( globalI ); // source of reverse edge
      }
    }
  }

  int myEdgeCount = static_cast< int >( myEdgesSrc.size() );

  // Gather counts from all ranks
  array1d< int > allEdgeCounts( numRanks );
  MpiWrapper::allgather( &myEdgeCount, 1, allEdgeCounts.data(), 1, comm );

  // Compute displacements for symmetrization
  array1d< int > edgeDispls( numRanks + 1 );
  edgeDispls[0] = 0;
  for( int r = 0; r < numRanks; ++r )
  {
    edgeDispls[r + 1] = edgeDispls[r] + allEdgeCounts[r];
  }

  int totalEdges = edgeDispls[numRanks];

  // Gather all edges
  array1d< pmet_idx_t > allEdgeSrc( totalEdges > 0 ? totalEdges : 1 );
  array1d< pmet_idx_t > allEdgeDst( totalEdges > 0 ? totalEdges : 1 );

  if( myEdgeCount > 0 )
  {
    MpiWrapper::allgatherv( myEdgesSrc.data(), myEdgeCount,
                            allEdgeSrc.data(), allEdgeCounts.data(), edgeDispls.data(), comm );
    MpiWrapper::allgatherv( myEdgesDst.data(), myEdgeCount,
                            allEdgeDst.data(), allEdgeCounts.data(), edgeDispls.data(), comm );
  }
  else
  {
    // Handle empty send
    pmet_idx_t dummy = 0;
    MpiWrapper::allgatherv( &dummy, 0,
                            allEdgeSrc.data(), allEdgeCounts.data(), edgeDispls.data(), comm );
    MpiWrapper::allgatherv( &dummy, 0,
                            allEdgeDst.data(), allEdgeCounts.data(), edgeDispls.data(), comm );
  }

  // Add reverse edges for those I own
  int reverseEdgesAdded = 0;
  for( int i = 0; i < totalEdges; ++i )
  {
    pmet_idx_t dst = allEdgeSrc[i];
    pmet_idx_t src = allEdgeDst[i];

    if( dst >= myGlobalStart && dst < myGlobalEnd )
    {
      localIndex localIdx = dst - myGlobalStart;
      if( neighborSets[localIdx].insert( src ).second )
      {
        reverseEdgesAdded++;
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 7: Build final ArrayOfArrays
  // -----------------------------------------------------------------------
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > superGraph;
  array1d< pmet_idx_t > offsets( numLocalSuperCells + 1 );
  offsets[0] = 0;

  for( localIndex i = 0; i < numLocalSuperCells; ++i )
  {
    offsets[i + 1] = offsets[i] + neighborSets[i].size();
  }

  superGraph.resizeFromOffsets( numLocalSuperCells, offsets.data() );

  for( localIndex i = 0; i < numLocalSuperCells; ++i )
  {
    superGraph.appendToArray( i, neighborSets[i].begin(), neighborSets[i].end() );
  }

  // -----------------------------------------------------------------------
  // Step 8: Statistics
  // -----------------------------------------------------------------------
  localIndex totalSuperEdges = 0;
  for( localIndex i = 0; i < superGraph.size(); ++i )
  {
    totalSuperEdges += superGraph.sizeOfArray( i );
  }

  localIndex totalBaseEdges = 0;
  for( localIndex i = 0; i < baseGraph.size(); ++i )
  {
    totalBaseEdges += baseGraph.sizeOfArray( i );
  }

  vtkIdType globalSuperCells = MpiWrapper::sum(
    static_cast< vtkIdType >(numLocalSuperCells), comm );
  vtkIdType globalSuperEdges = MpiWrapper::sum(
    static_cast< vtkIdType >(totalSuperEdges), comm );
  vtkIdType globalBaseCells = MpiWrapper::sum(
    static_cast< vtkIdType >(numLocalCells), comm );
  vtkIdType globalBaseEdges = MpiWrapper::sum(
    static_cast< vtkIdType >(totalBaseEdges), comm );

  GEOS_LOG_RANK_0( "SUPER-CELL GRAPH:" );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Super-graph nodes:     {:>10L}", globalSuperCells ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Super-graph edges:     {:>10L}", globalSuperEdges ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Base graph nodes:      {:>10L}", globalBaseCells ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Base graph edges:      {:>10L}", globalBaseEdges ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "  Node reduction:        {:>10L} cells ({:.1f}%)",
                           globalBaseCells - globalSuperCells,
                           100.0 * (globalBaseCells - globalSuperCells) / globalBaseCells ) );
  return std::make_pair( std::move( superGraph ), std::move( superVertexWeights ) );
}

void validateSuperCellGraph(
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const & superGraph,
  arrayView1d< pmet_idx_t const > const & superElemDist,
  arrayView1d< pmet_idx_t const > const & vertexWeights,
  MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );

  GEOS_LOG_RANK_0( "Running graph integrity checks..." );

  pmet_idx_t myStart = superElemDist[rank];
  pmet_idx_t myEnd = superElemDist[rank + 1];

  // -----------------------------------------------------------------------
  // Check 1: No self-loops
  // -----------------------------------------------------------------------
  {
    for( localIndex i = 0; i < superGraph.size(); ++i )
    {
      pmet_idx_t myGlobalId = myStart + i;

      for( localIndex j = 0; j < superGraph.sizeOfArray( i ); ++j )
      {
        GEOS_ERROR_IF( superGraph[i][j] == myGlobalId,
                       GEOS_FMT( "Rank {}: Super-cell {} (global {}) has self-loop",
                                 rank, i, myGlobalId ) );
      }
    }
  }

  // -----------------------------------------------------------------------
  // Check 2: All neighbors in valid global range
  // -----------------------------------------------------------------------
  {
    pmet_idx_t globalMin = 0;
    pmet_idx_t globalMax = superElemDist[numRanks] - 1;

    for( localIndex i = 0; i < superGraph.size(); ++i )
    {
      for( localIndex j = 0; j < superGraph.sizeOfArray( i ); ++j )
      {
        pmet_idx_t neighbor = superGraph[i][j];

        GEOS_ERROR_IF( neighbor< globalMin || neighbor > globalMax,
                       GEOS_FMT( "Rank {}: Super-cell {} has out-of-range neighbor {} (valid: [{}, {}])",
                                 rank, i, neighbor, globalMin, globalMax ) );
      }
    }
  }

  // -----------------------------------------------------------------------
  // Check 3: No duplicate edges
  // -----------------------------------------------------------------------
  {
    for( localIndex i = 0; i < superGraph.size(); ++i )
    {
      std::set< pmet_idx_t > uniqueNeighbors;

      for( localIndex j = 0; j < superGraph.sizeOfArray( i ); ++j )
      {
        pmet_idx_t neighbor = superGraph[i][j];

        GEOS_ERROR_IF( !uniqueNeighbors.insert( neighbor ).second,
                       GEOS_FMT( "Rank {}: Super-cell {} has duplicate neighbor {}",
                                 rank, i, neighbor ) );
      }
    }
  }

  // -----------------------------------------------------------------------
  // Check 4: Isolated vertices
  // -----------------------------------------------------------------------
  {
    localIndex isolated = 0;

    for( localIndex i = 0; i < superGraph.size(); ++i )
    {
      if( superGraph.sizeOfArray( i ) == 0 )
      {
        isolated++;
        if( isolated <= 5 )
        {
          pmet_idx_t globalId = myStart + i;
          GEOS_LOG_RANK( GEOS_FMT( "WARNING: Super-cell {} (global {}) has no neighbors (isolated)",
                                   i, globalId ) );
        }
      }
    }
    GEOS_WARNING_IF( isolated > 0, GEOS_FMT( "Found {} isolated vertices ", isolated ) );
  }

  // -----------------------------------------------------------------------
  // Check 5: Verify weights
  // -----------------------------------------------------------------------
  {
    if( !vertexWeights.empty() )
    {
      GEOS_ERROR_IF( vertexWeights.size() != superGraph.size(),
                     GEOS_FMT( "Rank {}: Vertex weights size ({}) != graph size ({})",
                               rank, vertexWeights.size(), superGraph.size() ) );

      for( localIndex i = 0; i < vertexWeights.size(); ++i )
      {
        GEOS_ERROR_IF( vertexWeights[i] <= 0,
                       GEOS_FMT( "Rank {}: Super-cell {} has invalid weight {}",
                                 rank, i, vertexWeights[i] ) );
      }
    }
  }

  // -----------------------------------------------------------------------
  // Check 6: Graph symmetry (local edges only)
  // -----------------------------------------------------------------------
  GEOS_LOG_RANK_0( "Checking local graph symmetry..." );
  {
    std::unordered_set< uint64_t > localOutgoingNeighbors;

    // Build edge set
    for( localIndex i = 0; i < superGraph.size(); ++i )
    {
      pmet_idx_t globalI = myStart + i;
      auto neighbors = superGraph[i];

      for( localIndex j = 0; j < neighbors.size(); ++j )
      {
        pmet_idx_t globalJ = neighbors[j];

        if( globalJ >= myStart && globalJ < myEnd )
        {
          uint64_t edgeKey = (static_cast< uint64_t >( globalI ) << 32) |
                             static_cast< uint64_t >( globalJ );
          localOutgoingNeighbors.insert( edgeKey );
        }
      }
    }

    // Verify symmetry
    for( localIndex i = 0; i < superGraph.size(); ++i )
    {
      pmet_idx_t globalI = myStart + i;
      auto neighbors = superGraph[i];

      for( localIndex j = 0; j < neighbors.size(); ++j )
      {
        pmet_idx_t globalJ = neighbors[j];

        if( globalJ >= myStart && globalJ < myEnd )
        {
          uint64_t reverseKey = (static_cast< uint64_t >( globalJ ) << 32) |
                                static_cast< uint64_t >( globalI );

          GEOS_ERROR_IF( localOutgoingNeighbors.find( reverseKey ) == localOutgoingNeighbors.end(),
                         GEOS_FMT( "Rank {}: Asymmetric edge ({} -> {}) - reverse edge missing",
                                   rank, globalI, globalJ ) );
        }
      }
    }
  }

}


// =============================================================================================
// SECTION 5: PARTITION UNPACKING (after ParMETIS)
// =============================================================================================
array1d< int64_t >
unpackSuperCellPartitioning(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  array1d< int64_t > const & superPartitioning,
  stdMap< vtkIdType, localIndex > const & superCellIdToLocalIdx,
  MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );

  GEOS_LOG_RANK_0( " UNPACKING SUPER-CELL PARTITIONING" );

  // -----------------------------------------------------------------------
  // Step 1: Get super-cell ID array
  // -----------------------------------------------------------------------
  vtkIdTypeArray * superCellIdArray =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetArray( "SuperCellId" ) );

  GEOS_ERROR_IF( superCellIdArray == nullptr,
                 GEOS_FMT( "Rank {}: SuperCellId array not found", rank ) );

  GEOS_ERROR_IF( static_cast< size_t >(superPartitioning.size()) != superCellIdToLocalIdx.size(),
                 GEOS_FMT( "Rank {}: Super-cell partitioning size ({}) doesn't match number of local super-cells ({})",
                           rank, superPartitioning.size(), superCellIdToLocalIdx.size() ) );

  // -----------------------------------------------------------------------
  // Step 2: Assign cells to ranks based on super-cell partitioning
  // -----------------------------------------------------------------------
  vtkIdType const numCells = cells3D->GetNumberOfCells();
  array1d< int64_t > cellPartitioning( numCells );

  for( vtkIdType i = 0; i < numCells; ++i )
  {
    vtkIdType superCellId = superCellIdArray->GetValue( i );
    auto it = superCellIdToLocalIdx.find( superCellId );

    GEOS_ERROR_IF( it == superCellIdToLocalIdx.end(),
                   GEOS_FMT( "Rank {}: Cell {} has unknown super-cell ID {}",
                             rank, i, superCellId ) );

    cellPartitioning[i] = superPartitioning[it->second];
  }

  // -----------------------------------------------------------------------
  // Step 3: Validate that super-cells weren't split
  // -----------------------------------------------------------------------
  stdMap< vtkIdType, std::set< int64_t > > superCellToRanks;

  for( vtkIdType i = 0; i < numCells; ++i )  // ← Reuse numCells (no redeclaration)
  {
    vtkIdType const scId = superCellIdArray->GetValue( i );
    superCellToRanks.get_inserted( scId ).insert( cellPartitioning[i] );
  }

  vtkIdType numSplitSuperCells = 0;
  for( auto const & [superCellId, ranks] : superCellToRanks )
  {
    if( ranks.size() > 1 )
    {
      ++numSplitSuperCells;
    }
  }

  vtkIdType totalSplitSuperCells = MpiWrapper::sum( numSplitSuperCells, comm );

  GEOS_ERROR_IF( totalSplitSuperCells > 0,
                 GEOS_FMT( "Partitioning failed: {:L} super-cells were split across ranks!",
                           totalSplitSuperCells ) );

  GEOS_LOG_RANK_0( "Super-cell partitioning validated successfully" );

  return cellPartitioning;
}

} // namespace vtk
} // namespace geos
