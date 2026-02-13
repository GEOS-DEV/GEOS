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



namespace geos
{

namespace vtk
{

// =============================================================================================
// SECTION 1: SUPER-CELL IDENTIFICATION (rank 0 only)
// =============================================================================================
SuperCellInfo tagCellsWithSuperCellIds(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  stdMap< string, ArrayOfArrays< vtkIdType, int64_t > > const & fractureNeighbors,
  MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );

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

    GEOS_LOG_RANK( "Fracture '" << fractureName
                                << "': processed " << numFracCells << " fracture elements" );
  }

  GEOS_LOG_RANK( "Rank " << rank << ": Built fracture graph with "
                         << fractureGraph.size() << " 3D cells having fracture connections, "
                         << totalFracturePairs << " fracture pairs" );

  // -----------------------------------------------------------------------
  // Step 2: Find the maximum global ID to avoid collisions
  // -----------------------------------------------------------------------

  vtkIdType maxGlobalId = 0;
  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType gid = globalIds->GetValue( i );
    maxGlobalId = std::max( maxGlobalId, gid );
  }

  GEOS_LOG_RANK( "Max global ID on rank " << rank << ": " << maxGlobalId );

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

    visited.insert( cell );
    cellToSuperCell[cell] = superCellId;
    component.push_back( cell );

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

  GEOS_LOG_RANK( "Found " << superCellComponents.size() << " connected components via DFS" );

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
    info.vertexWeights[scId] = members.size();

    if( members.size() > 1 )
    {
      info.atomicSuperCells.insert( scId );
    }

    numCellsInSuperCells += members.size();
    numSuperCellsCreated++;
    largestSuperCellSize = std::max( largestSuperCellSize,
                                     static_cast< vtkIdType >(members.size()) );
  }

  // -----------------------------------------------------------------------
  // Step 5: Tag cells with SuperCellIds
  // -----------------------------------------------------------------------
  vtkNew< vtkIdTypeArray > superCellIdArray;
  superCellIdArray->SetName( "SuperCellId" );
  superCellIdArray->SetNumberOfTuples( numLocalCells );

  vtkIdType numRegularCells = 0;

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
      numRegularCells++;
    }
  }

  cells3D->GetCellData()->AddArray( superCellIdArray );

  GEOS_LOG_RANK( "Tagged " << numLocalCells << " cells with SuperCellIds" );

  // -----------------------------------------------------------------------
  // Step 6: Report statistics
  // -----------------------------------------------------------------------

  vtkIdType globalNumCellsInSuperCells = numCellsInSuperCells;
  vtkIdType globalNumSuperCells = numSuperCellsCreated;
  vtkIdType globalLargestSize = largestSuperCellSize;
  vtkIdType globalRegularCells = numRegularCells;
  vtkIdType globalTotalCells = numLocalCells;

  vtkIdType globalTotalSuperCells = globalNumSuperCells + globalRegularCells;
  vtkIdType cellReduction = globalTotalCells - globalTotalSuperCells;

  GEOS_LOG_RANK_0( "SUPER-CELL TAGGING SUMMARY" );
  GEOS_LOG_RANK_0( "  Total 3D cells:                   "
                   << std::setw( 8 ) << globalTotalCells << std::setw( 8 )  );
  GEOS_LOG_RANK_0( "  Cells in super-cells:             "
                   << std::setw( 8 ) << globalNumCellsInSuperCells << std::setw( 8 )  );
  GEOS_LOG_RANK_0( "  Regular cells (no fractures):     "
                   << std::setw( 8 ) << globalRegularCells << std::setw( 8 ));
  GEOS_LOG_RANK_0( "  Number of super-cells created:    "
                   << std::setw( 8 ) << globalNumSuperCells << std::setw( 8 ) );
  GEOS_LOG_RANK_0( "  Total super-cells (incl regular): "
                   << std::setw( 8 ) << globalTotalSuperCells << std::setw( 8 ) );
  GEOS_LOG_RANK_0( "  Cell reduction:                   "
                   << std::setw( 8 ) << cellReduction << std::setw( 8 ) );
  GEOS_LOG_RANK_0( "  Largest super-cell size:          "
                   << std::setw( 8 ) << globalLargestSize << " cells" );

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

  // Build map: super-cell ID → vector of cell global IDs
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
    //info.vertexWeights[scId] = cells.size();

    // Apply same weighting strategy as initial tagging
    if( cells.size() > 1 )
    {
      // Fracture super-cell: weight = 10x
      info.vertexWeights[scId] = cells.size() + fractureWeight;
      info.atomicSuperCells.insert( scId );
    }
    else
    {
      // Regular cell: weight = 1
      info.vertexWeights[scId] = cells.size();    // = 1
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

  // Only rank 0 has the mesh and does the partitioning
  if( rank == 0 )
  {
    // Get SuperCellId array
    vtkIdTypeArray * superCellIdArray =
      vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetArray( "SuperCellId" ) );

    if( !superCellIdArray )
    {
      GEOS_LOG_RANK_0( "No SuperCellId array found - using simple redistribution" );
      // Fall back: just split evenly by cell index
      vtkIdType numCells = cells3D->GetNumberOfCells();

      partitionedMesh->SetNumberOfPartitions( numRanks );

      for( int r = 0; r < numRanks; ++r )
      {
        vtkNew< vtkIdList > cellsForRank;
        for( vtkIdType i = 0; i < numCells; ++i )
        {
          if( i % numRanks == r )
          {
            cellsForRank->InsertNextId( i );
          }
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
      // Build super-cell metadata
      std::map< vtkIdType, std::vector< vtkIdType > > superCellToLocalCells;
      vtkIdType numCells = cells3D->GetNumberOfCells();

      for( vtkIdType i = 0; i < numCells; ++i )
      {
        vtkIdType scId = superCellIdArray->GetValue( i );
        superCellToLocalCells[scId].push_back( i );
      }

      vtkIdType numSuperCells = superCellToLocalCells.size();
      GEOS_LOG_RANK_0( GEOS_FMT( "Redistributing {} super-cells ({} total cells) across {} ranks",
                                 numSuperCells, numCells, numRanks ));

      // Compute centroid for each super-cell
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
        // Compute centroid of this super-cell
        std::array< double, 3 > centroid = {0.0, 0.0, 0.0};
        vtkIdType totalPoints = 0;

        for( vtkIdType cellIdx : cellIndices )
        {
          vtkCell * cell = cells3D->GetCell( cellIdx );
          vtkPoints * points = cell->GetPoints();
          vtkIdType nPts = points->GetNumberOfPoints();

          for( vtkIdType p = 0; p < nPts; ++p )
          {
            double pt[3];
            points->GetPoint( p, pt );
            centroid[0] += pt[0];
            centroid[1] += pt[1];
            centroid[2] += pt[2];
            totalPoints++;
          }
        }

        if( totalPoints > 0 )
        {
          centroid[0] /= totalPoints;
          centroid[1] /= totalPoints;
          centroid[2] /= totalPoints;
        }

        superCells.push_back( SuperCellPartitionInfo{ scId, centroid, cellIndices } );
      }

      // Sort super-cells spatially using Z-order (Morton) curve
      // This provides good spatial locality in 3D
      auto computeMortonCode = []( std::array< double, 3 > const & centroid,
                                   double minCoord[3],
                                   double maxCoord[3] ) -> uint64_t
      {
        // Normalize coordinates to [0, 1]
        auto normalize = [&]( double val, int dim ) -> uint32_t
        {
          double range = maxCoord[dim] - minCoord[dim];
          if( range < 1e-10 )
            return 0;
          double norm = (val - minCoord[dim]) / range;
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

      // Find bounding box
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

      // Sort by Morton code for spatial locality
      std::sort( superCells.begin(), superCells.end(),
                 [&]( SuperCellPartitionInfo const & a, SuperCellPartitionInfo const & b )
      {
        return computeMortonCode( a.centroid, minCoord, maxCoord ) <
        computeMortonCode( b.centroid, minCoord, maxCoord );
      } );

      GEOS_LOG_RANK_0( "Super-cells sorted by spatial locality (Z-order curve)" );

      // Assign sorted super-cells to ranks in contiguous blocks
      array1d< int64_t > cellPartitions( numCells );
      std::vector< int > cellsPerRank( numRanks, 0 );

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


      // Build partitioned dataset - extract cells for each rank
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

        vtkNew< vtkExtractCells > extractor;
        extractor->SetInputData( cells3D );
        extractor->SetCellList( cellsForRank );
        extractor->Update();

        vtkSmartPointer< vtkUnstructuredGrid > partition =
          vtkUnstructuredGrid::SafeDownCast( extractor->GetOutput() );
        partitionedMesh->SetPartition( r, partition );
      }
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


  vtkSmartPointer< vtkDataSet > result = vtk::redistribute( *partitionedMesh, comm );

  partitionedMesh = nullptr;
  if( rank == 0 )
    cells3D = nullptr;

  vtkIdType localCells = result->GetNumberOfCells();

  // Verify SuperCellId array survived redistribution
  if( localCells > 0 )
  {
    vtkIdTypeArray * resultSuperCellIdArray =
      vtkIdTypeArray::SafeDownCast( result->GetCellData()->GetArray( "SuperCellId" ) );

    if( !resultSuperCellIdArray )
    {
      GEOS_ERROR( "SuperCellId array was lost during redistribution!" );
    }
  }

  return result;
}



// =============================================================================================
// SECTION 4: SUPER-CELL GRAPH BUILDING (for ParMETIS)
// =============================================================================================
std::pair< ArrayOfArrays< pmet_idx_t, pmet_idx_t >, array1d< pmet_idx_t > >
buildSuperCellGraph(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const & baseGraph,
  arrayView1d< pmet_idx_t const > const & baseElemDist,
  SuperCellInfo const & info,
  globalIndex const GEOS_UNUSED_PARAM( localStart ),
  MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );

  std::unordered_map< pmet_idx_t, vtkIdType > globalCellIdToSuperCellId;
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > superGraph;
  array1d< pmet_idx_t > superVertexWeights;

  // -----------------------------------------------------------------------
  // Step 1: Get super-cell ID array and validate
  // -----------------------------------------------------------------------

  vtkIdTypeArray * superCellIdArray =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetArray( "SuperCellId" ) );

  GEOS_ERROR_IF( superCellIdArray == nullptr,
                 "Rank " << rank << ": SuperCellId array not found" );

  vtkIdTypeArray * cellGlobalIds =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetGlobalIds() );

  GEOS_ERROR_IF( cellGlobalIds == nullptr,
                 "Rank " << rank << ": Cell global IDs not found" );

  vtkIdType numLocalCells = cells3D->GetNumberOfCells();

  GEOS_ERROR_IF( baseGraph.size() != numLocalCells,
                 "Rank " << rank << ": Base graph size (" << baseGraph.size()
                         << ") != mesh cell count (" << numLocalCells << ")" );

  // -----------------------------------------------------------------------
  // Step 2: Map local cells to super-cells
  // -----------------------------------------------------------------------
  std::map< vtkIdType, std::vector< vtkIdType > > superCellToLocalCells;

  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType superCellId = superCellIdArray->GetValue( i );
    superCellToLocalCells[superCellId].push_back( i );
  }

  localIndex numLocalSuperCells = superCellToLocalCells.size();

  // -----------------------------------------------------------------------
  // Step 3: Assign GLOBAL super-cell indices
  // -----------------------------------------------------------------------
  array1d< pmet_idx_t > superElemDist( numRanks + 1 );
  pmet_idx_t localSuperCellCount = numLocalSuperCells;

  MPI_Allgather( &localSuperCellCount, 1, MPI_LONG_LONG_INT,
                 superElemDist.data(), 1, MPI_LONG_LONG_INT, comm );

  pmet_idx_t temp = superElemDist[0];
  superElemDist[0] = 0;
  for( int r = 1; r <= numRanks; ++r )
  {
    pmet_idx_t next = superElemDist[r];
    superElemDist[r] = superElemDist[r-1] + temp;
    temp = next;
  }

  pmet_idx_t myGlobalStart = superElemDist[rank];

  // Build ordered list FIRST, then assign indices
  std::vector< vtkIdType > orderedSuperCellIds;
  orderedSuperCellIds.reserve( numLocalSuperCells );

  for( auto const & [superCellId, localCells] : superCellToLocalCells )
  {
    orderedSuperCellIds.push_back( superCellId );
  }

  // Sort to ensure deterministic ordering
  std::sort( orderedSuperCellIds.begin(), orderedSuperCellIds.end() );

  // Now assign global indices in sorted order
  std::map< vtkIdType, pmet_idx_t > superCellIdToGlobalIdx;
  for( pmet_idx_t i = 0; i < numLocalSuperCells; ++i )
  {
    vtkIdType superCellId = orderedSuperCellIds[i];
    pmet_idx_t globalIdx = myGlobalStart + i;

    superCellIdToGlobalIdx[superCellId] = globalIdx;
  }

  GEOS_ERROR_IF( superCellIdToGlobalIdx.size() != static_cast< size_t >( numLocalSuperCells ),
                 "Rank " << rank << ": Mapping size mismatch: "
                         << superCellIdToGlobalIdx.size() << " != " << numLocalSuperCells );


  // -----------------------------------------------------------------------
  // Step 4: Build GLOBAL map (globalCellId -> SuperCellId) via MPI
  // -----------------------------------------------------------------------
  {
    std::vector< pmet_idx_t > sendGlobalIds;
    std::vector< vtkIdType > sendSuperCellIds;

    sendGlobalIds.reserve( numLocalCells );
    sendSuperCellIds.reserve( numLocalCells );

    for( vtkIdType i = 0; i < numLocalCells; ++i )
    {
      sendGlobalIds.push_back( cellGlobalIds->GetValue( i ) );
      sendSuperCellIds.push_back( superCellIdArray->GetValue( i ) );
    }

    // Gather counts from all ranks
    int localCount = sendGlobalIds.size();
    std::vector< int > allCounts( numRanks );
    MPI_Allgather( &localCount, 1, MPI_INT, allCounts.data(), 1, MPI_INT, comm );

    std::vector< int > displs( numRanks + 1, 0 );
    for( int r = 0; r < numRanks; ++r )
    {
      displs[r+1] = displs[r] + allCounts[r];
    }

    int totalMappings = displs[numRanks];
    {
      std::vector< pmet_idx_t > allGlobalIds( totalMappings );
      std::vector< vtkIdType > allSuperCellIds( totalMappings );

      // All-gather the mappings
      MPI_Allgatherv( sendGlobalIds.data(), localCount, MPI_LONG_LONG_INT,
                      allGlobalIds.data(), allCounts.data(), displs.data(),
                      MPI_LONG_LONG_INT, comm );

      MPI_Allgatherv( sendSuperCellIds.data(), localCount, MPI_LONG_LONG_INT,
                      allSuperCellIds.data(), allCounts.data(), displs.data(),
                      MPI_LONG_LONG_INT, comm );

      // Build GLOBAL map (includes cells from ALL ranks)
      globalCellIdToSuperCellId.reserve( totalMappings );
      for( int i = 0; i < totalMappings; ++i )
      {
        globalCellIdToSuperCellId[allGlobalIds[i]] = allSuperCellIds[i];
      }
    }
  } // Step 4 scope ends here

  // -----------------------------------------------------------------------
  // Step 5: Exchange super-cell global indices
  // -----------------------------------------------------------------------
  std::vector< vtkIdType > sendSCIds;
  std::vector< pmet_idx_t > sendSCGlobalIndices;

  for( auto const & [scId, gIdx] : superCellIdToGlobalIdx )
  {
    sendSCIds.push_back( scId );
    sendSCGlobalIndices.push_back( gIdx );
  }

  int localSCCount = sendSCIds.size();
  std::vector< int > allSCCounts( numRanks );
  MPI_Allgather( &localSCCount, 1, MPI_INT, allSCCounts.data(), 1, MPI_INT, comm );

  std::vector< int > scDispls( numRanks + 1, 0 );
  for( int r = 0; r < numRanks; ++r )
  {
    scDispls[r+1] = scDispls[r] + allSCCounts[r];
  }

  int totalSCMappings = scDispls[numRanks];
  std::vector< vtkIdType > allSCIds( totalSCMappings );
  std::vector< pmet_idx_t > allSCGlobalIndices( totalSCMappings );

  MPI_Allgatherv( sendSCIds.data(), localSCCount, MPI_LONG_LONG_INT,
                  allSCIds.data(), allSCCounts.data(), scDispls.data(),
                  MPI_LONG_LONG_INT, comm );

  MPI_Allgatherv( sendSCGlobalIndices.data(), localSCCount, MPI_LONG_LONG_INT,
                  allSCGlobalIndices.data(), allSCCounts.data(), scDispls.data(),
                  MPI_LONG_LONG_INT, comm );

  // Add ALL super-cell mappings to our local map
  for( int i = 0; i < totalSCMappings; ++i )
  {
    superCellIdToGlobalIdx[allSCIds[i]] = allSCGlobalIndices[i];
  }

  // Verify coverage
  pmet_idx_t totalGlobalSuperCells = superElemDist[numRanks];
  std::set< pmet_idx_t > assignedIndices;

  for( int i = 0; i < totalSCMappings; ++i )
  {
    assignedIndices.insert( allSCGlobalIndices[i] );
  }

  GEOS_LOG_RANK_0( "Global super-cell index coverage: "
                   << assignedIndices.size() << " / " << totalGlobalSuperCells );

// -----------------------------------------------------------------------
// Step 6: Build super-cell graph edges
// -----------------------------------------------------------------------

  GEOS_LOG_RANK_0( "Building super-cell graph edges..." );

// -----------------------------------------------------------------------
// Step 6.1: Build Index Translator AND Step 6.4: Build edges
// -----------------------------------------------------------------------

  { // START BIG SCOPE

    std::unordered_map< pmet_idx_t, vtkIdType > parmetisToVtkId;
    std::unordered_map< vtkIdType, pmet_idx_t > vtkToParmetisId;

    // Build translator maps
    {
      pmet_idx_t myParmetisStart = baseElemDist[rank];

      for( vtkIdType i = 0; i < numLocalCells; ++i )
      {
        pmet_idx_t parmetisIdx = myParmetisStart + i;
        vtkIdType vtkGlobalId = cellGlobalIds->GetValue( i );

        parmetisToVtkId[parmetisIdx] = vtkGlobalId;
        vtkToParmetisId[vtkGlobalId] = parmetisIdx;
      }

      // Exchange mappings via AllGather
      std::vector< pmet_idx_t > sendParmetisIndices;
      std::vector< vtkIdType > sendVtkIds;

      sendParmetisIndices.reserve( numLocalCells );
      sendVtkIds.reserve( numLocalCells );

      for( vtkIdType i = 0; i < numLocalCells; ++i )
      {
        sendParmetisIndices.push_back( myParmetisStart + i );
        sendVtkIds.push_back( cellGlobalIds->GetValue( i ) );
      }

      int translatorLocalCount = numLocalCells;
      std::vector< int > translatorCounts( numRanks );
      MPI_Allgather( &translatorLocalCount, 1, MPI_INT,
                     translatorCounts.data(), 1, MPI_INT, comm );

      std::vector< int > translatorDispls( numRanks + 1, 0 );
      for( int r = 0; r < numRanks; ++r )
      {
        translatorDispls[r+1] = translatorDispls[r] + translatorCounts[r];
      }

      int translatorTotalMappings = translatorDispls[numRanks];

      // Nested scope for temporary vectors
      {
        std::vector< pmet_idx_t > allParmetisIndices( translatorTotalMappings );
        std::vector< vtkIdType > allVtkIds( translatorTotalMappings );

        MPI_Allgatherv( sendParmetisIndices.data(), translatorLocalCount, MPI_LONG_LONG_INT,
                        allParmetisIndices.data(), translatorCounts.data(),
                        translatorDispls.data(), MPI_LONG_LONG_INT, comm );

        MPI_Allgatherv( sendVtkIds.data(), translatorLocalCount, MPI_LONG_LONG_INT,
                        allVtkIds.data(), translatorCounts.data(),
                        translatorDispls.data(), MPI_LONG_LONG_INT, comm );

        parmetisToVtkId.reserve( translatorTotalMappings );
        vtkToParmetisId.reserve( translatorTotalMappings );

        for( int i = 0; i < translatorTotalMappings; ++i )
        {
          parmetisToVtkId[allParmetisIndices[i]] = allVtkIds[i];
          vtkToParmetisId[allVtkIds[i]] = allParmetisIndices[i];
        }
      }

      pmet_idx_t totalCells = baseElemDist[numRanks];
      GEOS_ERROR_IF( static_cast< pmet_idx_t >( parmetisToVtkId.size() ) != totalCells,
                     "Rank " << rank << ": Translator size mismatch: "
                             << parmetisToVtkId.size() << " != " << totalCells );
    }

    // -----------------------------------------------------------------------
    // Step 6.2: Verify base graph index range
    // -----------------------------------------------------------------------
    {
      pmet_idx_t totalCells = baseElemDist[numRanks];
      pmet_idx_t localMinNeighbor = std::numeric_limits< pmet_idx_t >::max();
      pmet_idx_t localMaxNeighbor = 0;

      for( localIndex i = 0; i < baseGraph.size(); ++i )
      {
        auto neighbors = baseGraph[i];
        for( localIndex j = 0; j < neighbors.size(); ++j )
        {
          pmet_idx_t neighborIdx = neighbors[j];
          if( neighbors.size() > 0 )
          {
            localMinNeighbor = std::min( localMinNeighbor, neighborIdx );
            localMaxNeighbor = std::max( localMaxNeighbor, neighborIdx );
          }
        }
      }

      pmet_idx_t globalMinNeighbor = MpiWrapper::min( localMinNeighbor, comm );
      pmet_idx_t globalMaxNeighbor = MpiWrapper::max( localMaxNeighbor, comm );

      GEOS_ERROR_IF( globalMaxNeighbor >= totalCells || globalMinNeighbor < 0,
                     "Base graph contains invalid ParMETIS indices! "
                     << "Range [" << globalMinNeighbor << ", " << globalMaxNeighbor
                     << "] exceeds valid ParMETIS range [0, " << (totalCells - 1) << "]" );

    }

    // -----------------------------------------------------------------------
    // Step 6.3: Initialize super-cell data structures
    // -----------------------------------------------------------------------
    superVertexWeights.resize( numLocalSuperCells );
    std::vector< std::set< pmet_idx_t > > neighborSets( numLocalSuperCells );

    // Statistics
    localIndex ghostNeighborCount = 0;
    localIndex localNeighborCount = 0;
    localIndex selfLoopCount = 0;
    localIndex translationFailures = 0;
    localIndex mapLookupFailures = 0;

    pmet_idx_t myStart = superElemDist[rank];
    pmet_idx_t myEnd = superElemDist[rank + 1];

    // -----------------------------------------------------------------------
    // Step 6.4: Build super-cell edges using correct index translation
    // -----------------------------------------------------------------------
    for( localIndex localSuperIdx = 0; localSuperIdx < numLocalSuperCells; ++localSuperIdx )
    {
      vtkIdType superCellId = orderedSuperCellIds[localSuperIdx];
      auto const & localCells = superCellToLocalCells.at( superCellId );

      // Set vertex weight
      auto itWeight = info.vertexWeights.find( superCellId );
      superVertexWeights[localSuperIdx] = (itWeight != info.vertexWeights.end())
                                      ? itWeight->second
                                      : localCells.size();

      // Collect neighbors from all 3D cells in this super-cell
      for( vtkIdType cellLocalIdx : localCells )
      {
        auto neighbors = baseGraph[cellLocalIdx];

        for( localIndex j = 0; j < neighbors.size(); ++j )
        {
          pmet_idx_t neighborParmetisIdx = neighbors[j];

          // Translate ParMETIS index -> VTK global ID
          auto itTranslate = parmetisToVtkId.find( neighborParmetisIdx );

          if( itTranslate == parmetisToVtkId.end() )
          {
            translationFailures++;
            if( translationFailures <= 3 )
            {
              GEOS_ERROR( "Rank " << rank << ": ParMETIS index " << neighborParmetisIdx
                                  << " not in translation map!" );
            }
            continue;
          }

          vtkIdType neighborVtkGlobalId = itTranslate->second;

          // Look up neighbor's SuperCellId using VTK global ID
          auto itSuperCell = globalCellIdToSuperCellId.find( neighborVtkGlobalId );

          if( itSuperCell == globalCellIdToSuperCellId.end() )
          {
            mapLookupFailures++;
            if( mapLookupFailures <= 3 )
            {
              GEOS_LOG( "Rank " << rank << ": VTK global ID " << neighborVtkGlobalId
                                << " (from ParMETIS idx " << neighborParmetisIdx
                                << ") not in globalCellIdToSuperCellId map!" );
            }
            continue;
          }

          vtkIdType neighborSuperCellId = itSuperCell->second;

          // Skip self-loops
          if( neighborSuperCellId == superCellId )
          {
            selfLoopCount++;
            continue;
          }

          // Convert neighbor's SuperCellId -> global super-cell index
          auto itGlobalIdx = superCellIdToGlobalIdx.find( neighborSuperCellId );

          if( itGlobalIdx == superCellIdToGlobalIdx.end() )
          {
            mapLookupFailures++;
            if( mapLookupFailures <= 3 )
            {
              GEOS_LOG( "Rank " << rank << ": Neighbor SuperCellId " << neighborSuperCellId
                                << " not in superCellIdToGlobalIdx map!" );
            }
            continue;
          }

          pmet_idx_t neighborGlobalSuperIdx = itGlobalIdx->second;

          // Add edge
          neighborSets[localSuperIdx].insert( neighborGlobalSuperIdx );

          // Track local vs ghost neighbors
          if( neighborGlobalSuperIdx >= myStart && neighborGlobalSuperIdx < myEnd )
          {
            localNeighborCount++;
          }
          else
          {
            ghostNeighborCount++;
          }
        }
      }
    }

    // -----------------------------------------------------------------------
    // Step 6.5: Report statistics
    // -----------------------------------------------------------------------

    localIndex totalUniqueNeighbors = 0;
    for( auto const & nset : neighborSets )
    {
      totalUniqueNeighbors += nset.size();
    }

    vtkIdType globalTranslationFailures = MpiWrapper::sum(
      static_cast< vtkIdType >( translationFailures ), comm );
    vtkIdType globalMapFailures = MpiWrapper::sum(
      static_cast< vtkIdType >( mapLookupFailures ), comm );

    GEOS_ERROR_IF( globalTranslationFailures > 0,
                   "Graph building failed: " << globalTranslationFailures
                                             << " ParMETIS indices could not be translated to VTK IDs!" );

    GEOS_ERROR_IF( globalMapFailures > 0,
                   "Graph building failed: " << globalMapFailures
                                             << " lookups failed in super-cell mapping!" );



    // -----------------------------------------------------------------------
    // Step 6.6: Symmetrize graph BEFORE building ArrayOfArrays
    // -----------------------------------------------------------------------

    { // Scope for symmetrization vectors
      // Pass 1: Collect all LOCAL edges (i -> j where we own i)
      std::vector< pmet_idx_t > localSrc, localDst;

      for( localIndex i = 0; i < numLocalSuperCells; ++i )
      {
        pmet_idx_t globalI = myStart + i;

        for( pmet_idx_t neighborIdx : neighborSets[i] )
        {
          localSrc.push_back( globalI );
          localDst.push_back( neighborIdx );
        }
      }

      // Gather edge counts
      int localEdgeCount = static_cast< int >( localSrc.size() );
      std::vector< int > edgeCountsSymm( numRanks );
      MPI_Allgather( &localEdgeCount, 1, MPI_INT, edgeCountsSymm.data(), 1, MPI_INT, comm );

      std::vector< int > displsSymm( numRanks + 1, 0 );
      for( int r = 0; r < numRanks; ++r )
      {
        displsSymm[r+1] = displsSymm[r] + edgeCountsSymm[r];
      }

      int totalEdgesSymm = displsSymm[numRanks];

      // Nested scope for the large all-gathered vectors
      {
        std::vector< pmet_idx_t > allSrcNodes( totalEdgesSymm );
        std::vector< pmet_idx_t > allDstNodes( totalEdgesSymm );

        // Gather ALL edges from ALL ranks
        MPI_Allgatherv( localSrc.data(), localEdgeCount, MPI_LONG_LONG_INT,
                        allSrcNodes.data(), edgeCountsSymm.data(), displsSymm.data(),
                        MPI_LONG_LONG_INT, comm );
        MPI_Allgatherv( localDst.data(), localEdgeCount, MPI_LONG_LONG_INT,
                        allDstNodes.data(), edgeCountsSymm.data(), displsSymm.data(),
                        MPI_LONG_LONG_INT, comm );

        // Pass 2: Add REVERSE edges to neighborSets
        // For each edge (A -> B), ensure edge (B -> A) exists on rank owning B

        int reverseEdgesAdded = 0;

        for( int e = 0; e < totalEdgesSymm; ++e )
        {
          pmet_idx_t src = allSrcNodes[e]; // Source of original edge
          pmet_idx_t dst = allDstNodes[e]; // Destination of original edge

          // Check if we OWN the DESTINATION (where reverse edge should be added)
          if( dst >= myStart && dst < myEnd )
          {
            localIndex localDstIdx = dst - myStart;

            // Add REVERSE edge: dst → src
            if( neighborSets[localDstIdx].insert( src ).second )
            {
              reverseEdgesAdded++;
            }
          }
        }

        if( reverseEdgesAdded > 0 )
        {
          GEOS_LOG_RANK( "Added " << reverseEdgesAdded << " reverse edges for symmetry" );
        }
      }
    }


    // -----------------------------------------------------------------------
    // Step 7: NOW build ArrayOfArrays from SYMMETRIZED neighborSets
    // -----------------------------------------------------------------------

    // Recompute total edges after symmetrization
    localIndex totalSymEdgesLocal = 0;
    for( localIndex i = 0; i < numLocalSuperCells; ++i )
    {
      totalSymEdgesLocal += neighborSets[i].size();
    }

    // Build offsets
    array1d< pmet_idx_t > offsets( numLocalSuperCells + 1 );
    offsets[0] = 0;
    for( localIndex i = 0; i < numLocalSuperCells; ++i )
    {
      offsets[i + 1] = offsets[i] + neighborSets[i].size();
    }

    // Build superGraph
    superGraph.resizeFromOffsets( numLocalSuperCells, offsets.data() );

    // Populate from symmetrized neighborSets
    for( localIndex i = 0; i < numLocalSuperCells; ++i )
    {
      superGraph.appendToArray( i, neighborSets[i].begin(), neighborSets[i].end() );
    }

    // Verify strict allocation
    GEOS_ERROR_IF( superGraph.valueCapacity() != offsets[numLocalSuperCells],
                   "Graph allocation mismatch after symmetrization: capacity="
                   << superGraph.valueCapacity() << ", expected=" << offsets[numLocalSuperCells] );

  } // END BIG SCOPE - parmetisToVtkId, vtkToParmetisId, neighborSets freed here

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

  vtkIdType globalSuperCells = MpiWrapper::sum( static_cast< vtkIdType >(numLocalSuperCells), comm );
  vtkIdType globalSuperEdges = MpiWrapper::sum( static_cast< vtkIdType >(totalSuperEdges), comm );
  vtkIdType globalBaseCells = MpiWrapper::sum( static_cast< vtkIdType >(numLocalCells), comm );
  vtkIdType globalBaseEdges = MpiWrapper::sum( static_cast< vtkIdType >(totalBaseEdges), comm );


  GEOS_LOG_RANK_0( "SUPER-CELL GRAPH:" );
  GEOS_LOG_RANK_0( "  Super-graph nodes:     " << std::setw( 10 ) << globalSuperCells );
  GEOS_LOG_RANK_0( "  Super-graph edges:     " << std::setw( 10 ) << globalSuperEdges );
  GEOS_LOG_RANK_0( "  Base graph nodes:      " << std::setw( 10 ) << globalBaseCells );
  GEOS_LOG_RANK_0( "  Base graph edges:      " << std::setw( 10 ) << globalBaseEdges );
  GEOS_LOG_RANK_0( "  Node reduction:        " << std::setw( 10 ) << (globalBaseCells - globalSuperCells) << " cells" );

  // -----------------------------------------------------------------------
  // Step 9: Validation
  // -----------------------------------------------------------------------
  pmet_idx_t minNeighbor = std::numeric_limits< pmet_idx_t >::max();
  pmet_idx_t maxNeighbor = 0;
  localIndex outOfRangeCount = 0;

  for( localIndex i = 0; i < superGraph.size(); ++i )
  {
    auto neighbors = superGraph[i];
    for( localIndex j = 0; j < neighbors.size(); ++j )
    {
      pmet_idx_t neighborIdx = neighbors[j];

      if( neighbors.size() > 0 )
      {
        minNeighbor = std::min( minNeighbor, neighborIdx );
        maxNeighbor = std::max( maxNeighbor, neighborIdx );
      }

      if( neighborIdx < 0 || neighborIdx >= globalSuperCells )
      {
        outOfRangeCount++;
        if( outOfRangeCount <= 3 )
        {
          GEOS_LOG_RANK( " ERROR: Super-cell " << i
                                               << " has out-of-range neighbor " << neighborIdx
                                               << " (valid: [0, " << globalSuperCells << "))" );
        }
      }
    }
  }
  vtkIdType globalOutOfRange = MpiWrapper::sum( static_cast< vtkIdType >(outOfRangeCount), comm );
  GEOS_ERROR_IF( globalOutOfRange > 0,
                 "Super-cell graph has " << globalOutOfRange << " out-of-range neighbor indices!" );

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

  localIndex numErrors = 0;
  localIndex numWarnings = 0;

  // -----------------------------------------------------------------------
  // Check 1: No self-loops
  // -----------------------------------------------------------------------
  localIndex selfLoops = 0;
  for( localIndex i = 0; i < superGraph.size(); ++i )
  {
    pmet_idx_t myGlobalId = superElemDist[rank] + i;

    for( localIndex j = 0; j < superGraph.sizeOfArray( i ); ++j )
    {
      if( superGraph[i][j] == myGlobalId )
      {
        selfLoops++;
        if( selfLoops <= 3 )
        {
          GEOS_LOG( "Rank " << rank << ": ERROR: Super-cell " << i
                            << " has self-loop (neighbor = " << myGlobalId << ")" );
        }
      }
    }
  }
  if( selfLoops > 0 )
  {
    GEOS_LOG( "Rank " << rank << ": Found " << selfLoops << " self-loops!" );
    numErrors += selfLoops;
  }

  // -----------------------------------------------------------------------
  // Check 2: All neighbors in valid global range
  // -----------------------------------------------------------------------
  localIndex outOfRange = 0;
  pmet_idx_t globalMin = 0;
  pmet_idx_t globalMax = superElemDist[numRanks] - 1;

  for( localIndex i = 0; i < superGraph.size(); ++i )
  {
    for( localIndex j = 0; j < superGraph.sizeOfArray( i ); ++j )
    {
      pmet_idx_t neighbor = superGraph[i][j];

      if( neighbor < globalMin || neighbor > globalMax )
      {
        outOfRange++;
        if( outOfRange <= 3 )
        {
          GEOS_LOG( "Rank " << rank << ": ERROR: Super-cell " << i
                            << " has out-of-range neighbor " << neighbor
                            << " (valid: [" << globalMin << ", " << globalMax << "])" );
        }
      }
    }
  }
  if( outOfRange > 0 )
  {
    GEOS_LOG( "Rank " << rank << ": Found " << outOfRange << " out-of-range neighbors!" );
    numErrors += outOfRange;
  }

  // -----------------------------------------------------------------------
  // Check 3: Duplicate edges
  // -----------------------------------------------------------------------
  localIndex duplicates = 0;
  for( localIndex i = 0; i < superGraph.size(); ++i )
  {
    std::set< pmet_idx_t > uniqueNeighbors;

    for( localIndex j = 0; j < superGraph.sizeOfArray( i ); ++j )
    {
      pmet_idx_t neighbor = superGraph[i][j];

      if( !uniqueNeighbors.insert( neighbor ).second )
      {
        duplicates++;
        if( duplicates <= 3 )
        {
          GEOS_LOG( "Rank " << rank << ": WARNING: Super-cell " << i
                            << " has duplicate neighbor " << neighbor );
        }
      }
    }
  }
  if( duplicates > 0 )
  {
    GEOS_LOG( "Rank " << rank << ": Found " << duplicates << " duplicate edges!" );
    numWarnings += duplicates;
  }

  // -----------------------------------------------------------------------
  // Check 4: Isolated vertices (no neighbors)
  // -----------------------------------------------------------------------
  localIndex isolated = 0;
  std::vector< localIndex > isolatedIndices;

  for( localIndex i = 0; i < superGraph.size(); ++i )
  {
    if( superGraph.sizeOfArray( i ) == 0 )
    {
      isolated++;
      isolatedIndices.push_back( i );

      if( isolated <= 3 )
      {
        GEOS_LOG( "Rank " << rank << ": WARNING: Super-cell " << i
                          << " (global " << (superElemDist[rank] + i)
                          << ") has no neighbors (isolated)" );
      }
    }
  }

  if( isolated > 0 )
  {
    GEOS_LOG( "Rank " << rank << ": Found " << isolated << " isolated vertices" );
    numWarnings += isolated;
  }
  // -----------------------------------------------------------------------
  // Check 5: Verify weights match graph size
  // -----------------------------------------------------------------------
  if( !vertexWeights.empty() )
  {
    if( vertexWeights.size() != superGraph.size() )
    {
      GEOS_LOG( "Rank " << rank << ": ERROR: Vertex weights size ("
                        << vertexWeights.size() << ") != graph size ("
                        << superGraph.size() << ")" );
      numErrors++;
    }

    // Check for zero or negative weights
    localIndex badWeights = 0;
    for( localIndex i = 0; i < vertexWeights.size(); ++i )
    {
      if( vertexWeights[i] <= 0 )
      {
        badWeights++;
        if( badWeights <= 3 )
        {
          GEOS_LOG( "Rank " << rank << ": ERROR: Super-cell " << i
                            << " has invalid weight " << vertexWeights[i] );
        }
      }
    }
    if( badWeights > 0 )
    {
      GEOS_LOG( "Rank " << rank << ": Found " << badWeights << " invalid weights!" );
      numErrors += badWeights;
    }
  }


  // -----------------------------------------------------------------------
  // Check 6: Graph symmetry (REQUIRED by ParMETIS/METIS)
  // -----------------------------------------------------------------------
  GEOS_LOG_RANK_0( "Checking graph symmetry..." );

  // Build global edge list (i → j)
  std::set< std::pair< pmet_idx_t, pmet_idx_t > > localEdges;

  pmet_idx_t myStart = superElemDist[rank];

  for( localIndex i = 0; i < superGraph.size(); ++i )
  {
    pmet_idx_t globalI = myStart + i;
    auto neighbors = superGraph[i];

    for( localIndex j = 0; j < neighbors.size(); ++j )
    {
      pmet_idx_t globalJ = neighbors[j];
      localEdges.insert( {globalI, globalJ} );
    }
  }

  // Gather all edges on rank 0
  int localEdgeCount = static_cast< int >( localEdges.size() );
  std::vector< int > edgeCounts( numRanks );
  MPI_Gather( &localEdgeCount, 1, MPI_INT,
              edgeCounts.data(), 1, MPI_INT, 0, comm );

  int asymmetricCount = 0;

  if( rank == 0 )
  {
    std::vector< int > displsEdge( numRanks + 1, 0 );
    for( int r = 0; r < numRanks; ++r )
    {
      displsEdge[r+1] = displsEdge[r] + edgeCounts[r];
    }

    int totalEdges = displsEdge[numRanks];
    std::vector< pmet_idx_t > allEdgesI( totalEdges );
    std::vector< pmet_idx_t > allEdgesJ( totalEdges );

    // Copy local edges
    int offset = 0;
    for( auto const & [i, j] : localEdges )
    {
      allEdgesI[offset] = i;
      allEdgesJ[offset] = j;
      offset++;
    }

    // Gather from other ranks
    MPI_Gatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                 allEdgesI.data(), edgeCounts.data(), displsEdge.data(),
                 MPI_LONG_LONG_INT, 0, comm );
    MPI_Gatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                 allEdgesJ.data(), edgeCounts.data(), displsEdge.data(),
                 MPI_LONG_LONG_INT, 0, comm );

    // Check symmetry
    std::set< std::pair< pmet_idx_t, pmet_idx_t > > allEdges;
    for( int e = 0; e < totalEdges; ++e )
    {
      allEdges.insert( {allEdgesI[e], allEdgesJ[e]} );
    }

    for( auto const & [i, j] : allEdges )
    {
      // Check if reverse edge exists
      if( allEdges.find( {j, i} ) == allEdges.end() )
      {
        asymmetricCount++;
        if( asymmetricCount <= 5 )
        {
          GEOS_LOG_RANK_0( "  ASYMMETRIC: Edge (" << i << " → " << j
                                                  << ") exists but reverse (" << j << " → " << i << ") doesn't" );
        }
      }
    }

    if( asymmetricCount > 0 )
    {
      GEOS_LOG_RANK_0( " Graph has " << asymmetricCount
                                     << " asymmetric edges (will be symmetrized before ParMETIS)" );
      numWarnings += asymmetricCount;
    }
    else
    {
      GEOS_LOG_RANK_0( "Graph is symmetric (" << allEdges.size() << " edges checked)" );
    }
  }
  else
  {
    // Other ranks send their edges
    std::vector< pmet_idx_t > sendI, sendJ;
    for( auto const & [i, j] : localEdges )
    {
      sendI.push_back( i );
      sendJ.push_back( j );
    }

    MPI_Gatherv( sendI.data(), localEdgeCount, MPI_LONG_LONG_INT,
                 nullptr, nullptr, nullptr, MPI_LONG_LONG_INT, 0, comm );
    MPI_Gatherv( sendJ.data(), localEdgeCount, MPI_LONG_LONG_INT,
                 nullptr, nullptr, nullptr, MPI_LONG_LONG_INT, 0, comm );
  }

  MPI_Barrier( comm );

  // -----------------------------------------------------------------------
  // Global summary
  // -----------------------------------------------------------------------
  long long globalErrors = 0;
  long long globalWarnings = 0;
  long long localErr = numErrors;
  long long localWarn = numWarnings;

  MPI_Allreduce( &localErr, &globalErrors, 1, MPI_LONG_LONG_INT, MPI_SUM, comm );
  MPI_Allreduce( &localWarn, &globalWarnings, 1, MPI_LONG_LONG_INT, MPI_SUM, comm );

  GEOS_THROW_IF( globalErrors > 0,
                 GEOS_FMT( "Graph has integrity errors - cannot proceed with partitioning (total errors: {})",
                           globalErrors ), InputError );
  if( globalWarnings > 0 )
  {
    GEOS_LOG_RANK_0( "\n GRAPH INTEGRITY CHECK PASSED WITH WARNINGS" );
    GEOS_LOG_RANK_0( "   Total warnings: " << globalWarnings );
    if( asymmetricCount > 0 )
    {
      GEOS_LOG_RANK_0( "   (Graph will be symmetrized automatically)" );
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
                 "Rank " << rank << ": SuperCellId array not found" );

  GEOS_ERROR_IF( static_cast< size_t >(superPartitioning.size()) != superCellIdToLocalIdx.size(),
                 "Rank " << rank << ": Super-cell partitioning size ("
                         << superPartitioning.size() << ") doesn't match number of local super-cells ("
                         << superCellIdToLocalIdx.size() << ")" );

  // -----------------------------------------------------------------------
  // Step 2: Create partitioning array for original cells
  // -----------------------------------------------------------------------

  array1d< int64_t > cellPartitioning( cells3D->GetNumberOfCells() );

  // For each original cell, look up its super-cell and assign the same partition
  for( vtkIdType i = 0; i < cells3D->GetNumberOfCells(); ++i )
  {
    vtkIdType superCellId = superCellIdArray->GetValue( i );

    auto it = superCellIdToLocalIdx.find( superCellId );
    GEOS_ERROR_IF( it == superCellIdToLocalIdx.end(),
                   "Rank " << rank << ": Cell " << i
                           << " has super-cell ID " << superCellId
                           << " which is not in local super-cell map" );

    localIndex superIdx = it->second;

    GEOS_ERROR_IF( superIdx >= superPartitioning.size(),
                   "Rank " << rank << ": Super-cell index " << superIdx
                           << " out of range [0, " << superPartitioning.size() << ")" );

    int64_t targetRank = superPartitioning[superIdx];
    cellPartitioning[i] = targetRank;
  }

  //GEOS_LOG_RANK( "Assigned partitions to " << cells3D->GetNumberOfCells() << " cells" );

  // -----------------------------------------------------------------------
  // Step 3: Verify - All cells in same super-cell go to same rank
  // -----------------------------------------------------------------------

  stdMap< vtkIdType, std::set< int64_t > > superCellToRanks;

  vtkIdType const numCells2 = cells3D->GetNumberOfCells();
  for( vtkIdType i = 0; i < numCells2; ++i )
  {
    vtkIdType const scId = superCellIdArray->GetValue( i );
    int64_t const targetRank = cellPartitioning[i];

    superCellToRanks.get_inserted( scId ).insert( targetRank );
  }

  // Check for splits
  vtkIdType numSplitSuperCells = 0;

  for( auto const & [superCellId, ranks] : superCellToRanks )
  {
    if( ranks.size() > 1 )
    {
      GEOS_ERROR( "Rank " << rank << ": Super-cell " << superCellId
                          << " was split across " << ranks.size() << " ranks: {"
                          << stringutilities::join( std::vector< int64_t >( ranks.begin(), ranks.end()), ", " )
                          << "}" );
      ++numSplitSuperCells;
    }
  }

  GEOS_ERROR_IF( numSplitSuperCells > 0,
                 "Rank " << rank << ": " << numSplitSuperCells
                         << " super-cells were incorrectly split!" );


  // -----------------------------------------------------------------------
  // Step 4: Global statistics
  // -----------------------------------------------------------------------

  vtkIdType totalSplitSuperCells = MpiWrapper::sum( numSplitSuperCells, comm );

  if( totalSplitSuperCells > 0 )
  {
    GEOS_ERROR( totalSplitSuperCells << " super-cells were split globally!" );
  }

  // Count cells going to each rank
  std::map< int64_t, vtkIdType > cellsPerRank;
  for( vtkIdType i = 0; i < cellPartitioning.size(); ++i )
  {
    cellsPerRank[cellPartitioning[i]]++;
  }

  return cellPartitioning;
}


} // namespace vtk

} // namespace geos
