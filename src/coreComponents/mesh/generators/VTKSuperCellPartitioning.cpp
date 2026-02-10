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

/**
 * # Super-Cell Partitioning for Fracture Meshes
 *
 * ## Problem Statement
 *
 * When a fracture element connects two 3D cells:
 *
 *     Cell A
 *     ══════  Fracture surface (2D)
 *     Cell B
 *
 * Both cells MUST reside on the same MPI rank to:
 * - Enable proper fracture-to-3D neighbor assignment
 * - Avoid ghost layer complications
 * - Maintain topological consistency for contact mechanics
 *
 * Standard graph partitioning (ParMETIS/PTScotch) treats cells independently,
 * potentially assigning Cell A and Cell B to different ranks.
 *
 * ## Solution: Super-Cell Graph Coarsening
 *
 * We treat fracture-connected cell groups as atomic "super-cells" that cannot be split:
 *
 * 1. **Tag cells** (rank 0): Identify connected components via fractures, assign SuperCellId
 * 2. **Distribute**: Initial mesh distribution preserves SuperCellId array
 * 3. **Coarsen graph**: Build graph where nodes = super-cells (not cells)
 * 4. **Partition**: ParMETIS partitions super-cells (using vertex weights for load balancing)
 * 5. **Unpack**: Map super-cell assignments back to individual cells
 * 6. **Verify**: Ensure no super-cells were split
 *
 * After partitioning, all cells with the same SuperCellId are on the same rank.
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
  //GEOS_LOG_RANK( "Processing " << numLocalCells << " 3D cells" );

  vtkIdTypeArray * globalIds =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetGlobalIds() );

  GEOS_ERROR_IF( !globalIds,"3D mesh missing global IDs" );

  // -----------------------------------------------------------------------
  // Step 1: Build 3D cell connectivity graph via fractures
  // -----------------------------------------------------------------------
  
  // Map: 3D cell global ID → set of neighboring 3D cell global IDs (via fractures)
  std::map< vtkIdType, std::set< vtkIdType > > fractureGraph;

  vtkIdType totalFracturePairs = 0;
  
  // Build graph from fracture neighbor information
  for( auto const & [fractureName, neighbors] : fractureNeighbors )
  {
    vtkIdType numFracCells = neighbors.size();
    
    for( vtkIdType i = 0; i < numFracCells; ++i )
    {
      auto neighborList = neighbors[i];
      
      // Each fracture element connects (typically) 2 3D cells
      if( neighborList.size() >= 2 )
      {
        vtkIdType gidA = neighborList[0];
        vtkIdType gidB = neighborList[1];
        
        // Add bidirectional edge between the 3D cells
        fractureGraph[gidA].insert( gidB );
        fractureGraph[gidB].insert( gidA );
        
        totalFracturePairs++;
      }
      // Handle boundary fractures (only 1 neighbor) - these cells stay regular
    }
    
    GEOS_LOG_RANK( "Fracture '" << fractureName
                           << "': processed " << numFracCells << " fracture elements" );
  }

  GEOS_LOG_RANK( "Rank " << rank << ": Built fracture graph with "
                         << fractureGraph.size() << " 3D cells having fracture connections, "
                         << totalFracturePairs << " fracture pairs" );

  // -----------------------------------------------------------------------
  // Step 2: Find connected components using DFS (3D cells only)
  // -----------------------------------------------------------------------

  std::map< vtkIdType, vtkIdType > cellToSuperCell;
  std::set< vtkIdType > visited;
  vtkIdType nextSuperCellId = 0;

  // DFS to find connected component
  std::function< void(vtkIdType, vtkIdType, std::vector< vtkIdType > &) > dfs =
    [&]( vtkIdType cell, vtkIdType superCellId, std::vector< vtkIdType > & component )
  {
if( visited.count( cell ) )
      return;

    visited.insert( cell );
    cellToSuperCell[cell] = superCellId;
    component.push_back( cell );

    // Visit all fracture-connected 3D neighbors
    if( fractureGraph.count( cell ) )
    {
      for( vtkIdType neighbor : fractureGraph.at( cell ) )
      {
        dfs( neighbor, superCellId, component );
      }
    }
  };

  // Find all connected components
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
  // Step 3: Build SuperCellInfo (3D cells only)
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
  // Step 4: Tag 3D cells with SuperCellIds
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
      // Cell is part of a super-cell (has fracture neighbors)
      superCellIdArray->SetValue( i, cellToSuperCell.at( globalId ) );
    }
    else
    {
      // Regular cell (no fractures) - becomes its own super-cell
      superCellIdArray->SetValue( i, globalId );
      numRegularCells++;
    }
  }

  cells3D->GetCellData()->AddArray( superCellIdArray );

  GEOS_LOG_RANK( "Tagged " << numLocalCells 
                         << " cells with SuperCellIds" );

// -----------------------------------------------------------------------
  // Step 5: Report statistics
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
                   << std::setw( 8 ) << globalRegularCells << std::setw( 8));
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
SuperCellInfo reconstructSuperCellInfo( vtkSmartPointer< vtkUnstructuredGrid > mesh )
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
    info.vertexWeights[scId] = cells.size();

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

  GEOS_LOG_RANK_0( "INITIAL REDISTRIBUTION (Super-Cell Aware)" );

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

      // Assign each super-cell to a rank (round-robin)
      array1d< int64_t > cellPartitions( numCells );
      int targetRank = 0;
      
      for( auto const & [scId, cellIndices] : superCellToLocalCells )
      {
        // All cells in this super-cell go to the same rank
        for( vtkIdType cellIdx : cellIndices )
        {
          cellPartitions[cellIdx] = targetRank;
        }
        targetRank = (targetRank + 1) % numRanks;
      }

      // Count cells per rank for logging
      std::vector< int > cellsPerRank( numRanks, 0 );
      for( vtkIdType i = 0; i < numCells; ++i )
      {
        cellsPerRank[cellPartitions[i]]++;
      }

      GEOS_LOG_RANK_0( "Cell distribution plan:" );
      for( int r = 0; r < numRanks; ++r )
      {
        GEOS_LOG_RANK_0( GEOS_FMT( "  Rank {} will receive {} cells", r, cellsPerRank[r] ));
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
        
        GEOS_LOG_RANK_0( GEOS_FMT( "  Partition {}: {} cells extracted", r, partition->GetNumberOfCells() ));
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

  // All ranks participate in redistribution
  vtkSmartPointer< vtkDataSet > result = vtk::redistribute( *partitionedMesh, comm );
  
  vtkIdType localCells = result->GetNumberOfCells();
  //GEOS_LOG_RANK( GEOS_FMT( "Received {} cells after redistribution", localCells ));
  
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

  //GEOS_LOG_RANK( "Building super-cell graph from " << numLocalCells << " local cells" );

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
  //GEOS_LOG_RANK( "Local super-cells: " << numLocalSuperCells << " (from " << numLocalCells << " cells)" );

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
  //pmet_idx_t myGlobalEnd = superElemDist[rank + 1];
  //GEOS_LOG_RANK( "Super-cell global range: [" << myGlobalStart << ", " << myGlobalEnd << ")" );

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
  // Step 4: Build GLOBAL map (globalCellId → SuperCellId) via MPI
  // -----------------------------------------------------------------------
  // Collect local mappings
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
  std::map< pmet_idx_t, vtkIdType > globalCellIdToSuperCellId;
  for( int i = 0; i < totalMappings; ++i )
  {
    globalCellIdToSuperCellId[allGlobalIds[i]] = allSuperCellIds[i];
  }

  // -----------------------------------------------------------------------
  // Step 4.5: Build ParMETIS global index -> VTK global ID mapping
  // -----------------------------------------------------------------------
  // Build local mapping: position in allGlobalIds -> VTK global ID
  std::vector< vtkIdType > globalParmetisToVtk( totalMappings );
  for( int i = 0; i < totalMappings; ++i )
  {
    globalParmetisToVtk[i] = allGlobalIds[i];
  }


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

  //GEOS_LOG_RANK( "Total super-cell ID→index mappings: " << superCellIdToGlobalIdx.size() << " (all ranks)" );

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
// Step 6.1: Build Index Translator
// -----------------------------------------------------------------------

std::map< pmet_idx_t, vtkIdType > parmetisToVtkId;
std::map< vtkIdType, pmet_idx_t > vtkToParmetisId;

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
  std::vector< pmet_idx_t > allParmetisIndices( translatorTotalMappings );
  std::vector< vtkIdType > allVtkIds( translatorTotalMappings );
  
  MPI_Allgatherv( sendParmetisIndices.data(), translatorLocalCount, MPI_LONG_LONG_INT,
                  allParmetisIndices.data(), translatorCounts.data(), 
                  translatorDispls.data(), MPI_LONG_LONG_INT, comm );
  
  MPI_Allgatherv( sendVtkIds.data(), translatorLocalCount, MPI_LONG_LONG_INT,
                  allVtkIds.data(), translatorCounts.data(), 
                  translatorDispls.data(), MPI_LONG_LONG_INT, comm );
  
  for( int i = 0; i < translatorTotalMappings; ++i )
  {
    parmetisToVtkId[allParmetisIndices[i]] = allVtkIds[i];
    vtkToParmetisId[allVtkIds[i]] = allParmetisIndices[i];
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

array1d< pmet_idx_t > superVertexWeights( numLocalSuperCells );
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

  // Set vertex weight (number of 3D cells in this super-cell)
  auto itWeight = info.vertexWeights.find( superCellId );
  superVertexWeights[localSuperIdx] = (itWeight != info.vertexWeights.end())
                                      ? itWeight->second
                                      : localCells.size();

  // Collect neighbors from all 3D cells in this super-cell
  for( vtkIdType cellLocalIdx : localCells )
  {
    auto neighbors = baseGraph[cellLocalIdx];  // Returns ParMETIS indices

    for( localIndex j = 0; j < neighbors.size(); ++j )
    {
      pmet_idx_t neighborParmetisIdx = neighbors[j];

      // CRITICAL: Translate ParMETIS index -> VTK global ID
      auto itTranslate = parmetisToVtkId.find( neighborParmetisIdx );
      
      if( itTranslate == parmetisToVtkId.end() )
      {
        translationFailures++;
        if( translationFailures <= 3 )
        {
          GEOS_ERROR( "Rank " << rank << ": ParMETIS index " << neighborParmetisIdx  << " not in translation map!" );
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

      // Skip self-loops (internal edges within same super-cell)
      if( neighborSuperCellId == superCellId )
      {
        selfLoopCount++;
        continue;
      }

      // Convert neighbor's SuperCellId → global super-cell index
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

      // Add edge to neighbor set (automatically handles duplicates via std::set)
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

//GEOS_LOG_RANK( "Super-cell graph edge building complete:" );
//GEOS_LOG_RANK( "  Local neighbor edges:       " << localNeighborCount );
//GEOS_LOG_RANK( "  Ghost neighbor edges:       " << ghostNeighborCount );
//GEOS_LOG_RANK( "  Self-loops skipped:         " << selfLoopCount );
//GEOS_LOG_RANK( "  Translation failures:       " << translationFailures );
//GEOS_LOG_RANK( "  Map lookup failures:        " << mapLookupFailures );

localIndex totalUniqueNeighbors = 0;
for( auto const & nset : neighborSets )
{
  totalUniqueNeighbors += nset.size();
}

// Error if we had failures
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
  std::vector< pmet_idx_t > allSrcNodes( totalEdgesSymm );
  std::vector< pmet_idx_t > allDstNodes( totalEdgesSymm );

  // Gather ALL edges from ALL ranks
  MPI_Allgatherv( localSrc.data(), localEdgeCount, MPI_LONG_LONG_INT,
                  allSrcNodes.data(), edgeCountsSymm.data(), displsSymm.data(),
                  MPI_LONG_LONG_INT, comm );
  MPI_Allgatherv( localDst.data(), localEdgeCount, MPI_LONG_LONG_INT,
                  allDstNodes.data(), edgeCountsSymm.data(), displsSymm.data(),
                  MPI_LONG_LONG_INT, comm );

  //GEOS_LOG_RANK( "Collected " << totalEdgesSymm << " total edges from all ranks" );

  // Pass 2: Add REVERSE edges to neighborSets
  // For each edge (A -> B), ensure edge (B -> A) exists on rank owning B

  int reverseEdgesAdded = 0;

  for( int e = 0; e < totalEdgesSymm; ++e )
  {
    pmet_idx_t src = allSrcNodes[e];  // Source of original edge
    pmet_idx_t dst = allDstNodes[e];  // Destination of original edge
    
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

  if (reverseEdgesAdded > 0 )
  {
    GEOS_LOG_RANK( "Added " << reverseEdgesAdded << " reverse edges for symmetry" );
  }
  

  // -----------------------------------------------------------------------
  // Verify symmetry
  // -----------------------------------------------------------------------

  // Rebuild edge list after symmetrization
  std::vector< pmet_idx_t > verifyLocalSrc, verifyLocalDst;
  for( localIndex i = 0; i < numLocalSuperCells; ++i )
  {
    pmet_idx_t globalI = myStart + i;
    for( pmet_idx_t neighborIdx : neighborSets[i] )
    {
      verifyLocalSrc.push_back( globalI );
      verifyLocalDst.push_back( neighborIdx );
    }
  }

  int verifyLocalCount = static_cast< int >( verifyLocalSrc.size() );
  std::vector< int > verifyEdgeCounts( numRanks );
  MPI_Allgather( &verifyLocalCount, 1, MPI_INT, verifyEdgeCounts.data(), 1, MPI_INT, comm );

  std::vector< int > verifyDispls( numRanks + 1, 0 );
  for( int r = 0; r < numRanks; ++r )
  {
    verifyDispls[r+1] = verifyDispls[r] + verifyEdgeCounts[r];
  }

  int totalVerifyEdges = verifyDispls[numRanks];
  std::vector< pmet_idx_t > allVerifySrc( totalVerifyEdges );
  std::vector< pmet_idx_t > allVerifyDst( totalVerifyEdges );

  MPI_Allgatherv( verifyLocalSrc.data(), verifyLocalCount, MPI_LONG_LONG_INT,
                  allVerifySrc.data(), verifyEdgeCounts.data(), verifyDispls.data(),
                  MPI_LONG_LONG_INT, comm );
  MPI_Allgatherv( verifyLocalDst.data(), verifyLocalCount, MPI_LONG_LONG_INT,
                  allVerifyDst.data(), verifyEdgeCounts.data(), verifyDispls.data(),
                  MPI_LONG_LONG_INT, comm );

  // Build global edge set
  std::set< std::pair< pmet_idx_t, pmet_idx_t > > edgeSet;
  for( int e = 0; e < totalVerifyEdges; ++e )
  {
    edgeSet.insert( {allVerifySrc[e], allVerifyDst[e]} );
  }

  // Check symmetry
  int asymmetricLocal = 0;
  for( auto const & [src, dst] : edgeSet )
  {
    if( edgeSet.find( {dst, src} ) == edgeSet.end() )
    {
      asymmetricLocal++;
      if( asymmetricLocal <= 5 )
      {
        GEOS_LOG_RANK_0( "Missing reverse: (" << src << " → " << dst 
                         << ") exists but (" << dst << " → " << src << ") doesn't" );
      }
    }
  }

  int asymmetricGlobal = MpiWrapper::sum( asymmetricLocal, comm );

  if( asymmetricGlobal > 0 )
  {
    GEOS_ERROR( "Graph still has " << asymmetricGlobal << " asymmetric edges after symmetrization!" );
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

  // Create ArrayOfArrays with STRICT allocation
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > superGraph;
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
    GEOS_LOG_RANK_0( "  Super-graph nodes:     " << std::setw( 10 ) << globalSuperCells  );
    GEOS_LOG_RANK_0( "  Super-graph edges:     " << std::setw( 10 ) << globalSuperEdges );
    GEOS_LOG_RANK_0( "  Base graph nodes:      " << std::setw( 10 ) << globalBaseCells  );
    GEOS_LOG_RANK_0( "  Base graph edges:      " << std::setw( 10 ) << globalBaseEdges  );
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

  if( minNeighbor != std::numeric_limits< pmet_idx_t >::max() )
  {
    GEOS_LOG_RANK( "Neighbor index range: ["
                   << minNeighbor << ", " << maxNeighbor
                   << "], out-of-range: " << outOfRangeCount );
  }

  GEOS_LOG_RANK_0( "Expected global range: [0, " << globalSuperCells << ")" );
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
  
  if( globalErrors > 0 )
  {
    GEOS_LOG_RANK_0( "\nGRAPH INTEGRITY CHECK FAILED!" );
    GEOS_LOG_RANK_0( "   Total errors: " << globalErrors );
    GEOS_THROW( "Graph has integrity errors - cannot proceed with partitioning", 
                std::runtime_error );
  }
  else if( globalWarnings > 0 )
  {
    GEOS_LOG_RANK_0( "\n GRAPH INTEGRITY CHECK PASSED WITH WARNINGS" );
    GEOS_LOG_RANK_0( "   Total warnings: " << globalWarnings );
    if( asymmetricCount > 0 )
    {
      GEOS_LOG_RANK_0( "   (Graph will be symmetrized automatically)" );
    }
  }
  else
  {
    GEOS_LOG_RANK_0( "\n GRAPH INTEGRITY CHECK PASSED" );
    GEOS_LOG_RANK_0( "   No errors or warnings found" );
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
    GEOS_ERROR(  totalSplitSuperCells << " super-cells were split globally!" );
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