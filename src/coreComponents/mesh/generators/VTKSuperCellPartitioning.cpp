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
    return numCells;  // Regular cell
  }
}

SuperCellInfo tagCellsWithSuperCellIds(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  stdMap< string, ArrayOfArrays< localIndex, int64_t > > const & fractureNeighbors,

  integer fractureWeight )
{
  vtkIdType const numLocalCells = cells3D->GetNumberOfCells();

  vtkIdTypeArray * globalIds =
    vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetGlobalIds() );

  GEOS_ERROR_IF( !globalIds, "3D mesh missing global IDs" );

  // -----------------------------------------------------------------------
  // Step 1: Build 3D cell connectivity graph via fractures
  // -----------------------------------------------------------------------
  stdMap< vtkIdType, std::set< vtkIdType > > fractureGraph;
  vtkIdType totalFractureElements = 0;

  for( auto const & [fractureName, neighbors] : fractureNeighbors )
  {
    totalFractureElements += static_cast< vtkIdType >( neighbors.size() );

    for( vtkIdType i = 0; i < neighbors.size(); ++i )
    {
      auto neighborList = neighbors[i];

      GEOS_ERROR_IF( neighborList.size() != 2,
                     GEOS_FMT( "Fracture '{}' element {} has {} 3D neighbors (expected exactly 2). ",
                               fractureName, i, neighborList.size() ) );

      vtkIdType gidA = neighborList[0];
      vtkIdType gidB = neighborList[1];

      fractureGraph.get_inserted( gidA ).insert( gidB );
      fractureGraph.get_inserted( gidB ).insert( gidA );
    }
  }

  GEOS_LOG_RANK_0( GEOS_FMT( "Built fracture graph: {} 3D cells connected via {} fracture elements",
                             fractureGraph.size(), totalFractureElements ) );

  // -----------------------------------------------------------------------
  // Step 2: Find the maximum global ID from ALL cells (not just fracture cells)
  // -----------------------------------------------------------------------
  vtkIdType const maxGlobalId = numLocalCells > 0
                                ? globalIds->GetValueRange( 0 )[1]
                                : 0;

  // -----------------------------------------------------------------------
  // Step 3: Find connected components using DFS
  // -----------------------------------------------------------------------
  stdMap< vtkIdType, vtkIdType > cellToSuperCell;
  std::set< vtkIdType > visited;

  // Start super-cell IDs AFTER the maximum global ID to avoid collisions
  vtkIdType nextSuperCellId = maxGlobalId + 1;

  std::function< void(vtkIdType, vtkIdType, stdVector< vtkIdType > &) > dfs =
    [&]( vtkIdType cell, vtkIdType superCellId, stdVector< vtkIdType > & component )
  {
    if( visited.count( cell ) )
      return;
    visited.insert( cell );
    cellToSuperCell.try_emplace( cell, superCellId );
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

  stdMap< vtkIdType, stdVector< vtkIdType > > superCellComponents;

  for( auto const & [cell, neighbors] : fractureGraph )
  {
    if( !visited.count( cell ) )
    {
      stdVector< vtkIdType > component;
      dfs( cell, nextSuperCellId, component );

      superCellComponents.try_emplace( nextSuperCellId, component );
      nextSuperCellId++;
    }
  }

  // -----------------------------------------------------------------------
  // Step 4: Build SuperCellInfo
  // -----------------------------------------------------------------------
  SuperCellInfo info;

  vtkIdType numSuperCellsCreated = 0;
  vtkIdType largestSuperCellSize = 0;

  for( auto const & [scId, members] : superCellComponents )
  {
    info.superCellToOriginalCells.try_emplace( scId, members );
    info.vertexWeights.try_emplace( scId, computeSuperCellWeight( members.size(), fractureWeight ) );

    if( members.size() > 1 )
    {
      // Mark as atomic super-cell (cannot be split during partitioning)
      info.atomicSuperCells.insert( scId );
    }

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
  GEOS_LOG_RANK_0( GEOS_FMT( "Super-cells: {} created (largest: {} cells)",
                             numSuperCellsCreated, largestSuperCellSize ) );

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
  stdMap< vtkIdType, stdVector< vtkIdType > > localSuperCells;

  for( vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i )
  {
    vtkIdType scId = superCellIdArray->GetValue( i );
    vtkIdType globalId = globalIds->GetValue( i );
    localSuperCells.get_inserted( scId ).push_back( globalId );
  }

  for( auto const & [scId, cells] : localSuperCells )
  {
    info.superCellToOriginalCells.try_emplace( scId, cells );
    info.vertexWeights.try_emplace( scId, computeSuperCellWeight( cells.size(), fractureWeight ) );

    if( cells.size() > 1 )
    {
      info.atomicSuperCells.insert( scId );
    }
  }

  return info;
}


// =============================================================================================
// SECTION 3: INITIAL REDISTRIBUTION
// =============================================================================================
vtkSmartPointer< vtkDataSet >
redistributeBySuperCellBlocks( vtkSmartPointer< vtkUnstructuredGrid > cells3D,
                               MPI_Comm comm,
                               ScatterMethod scatterMethod,
                               arrayView1d< integer const > cartesianPartitions )
{
  GEOS_MARK_FUNCTION;

  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );

  if( scatterMethod == ScatterMethod::kdtree )
  {
    GEOS_LOG_RANK_0( "scatterMethod=kdtree is not supported with fractures (cannot preserve "
                     "super-cell atomicity). Automatically falling back to rcb." );
    scatterMethod = ScatterMethod::rcb;
  }

  // Per-cell rank assignment is computed on rank 0 (empty elsewhere) then consumed by
  // scatterByRankAssignment. Super-cell atomicity is enforced by assigning one rank
  // per super-cell and propagating it to every cell of that super-cell.
  stdVector< integer > cellRanks;

  if( rank == 0 )
  {
    vtkIdTypeArray * superCellIdArray =
      vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetArray( "SuperCellId" ) );
    GEOS_ERROR_IF( !superCellIdArray, "SuperCellId array not found" );

    vtkIdType const numCells = cells3D->GetNumberOfCells();

    // Group cells by super-cell ID. The iteration order of the map gives a stable
    // 0..numSuperCells-1 indexing used by the virtual "atom mesh" below.
    stdMap< vtkIdType, stdVector< vtkIdType > > superCellToCells;
    for( vtkIdType i = 0; i < numCells; ++i )
    {
      superCellToCells.get_inserted( superCellIdArray->GetValue( i ) ).push_back( i );
    }
    vtkIdType const numSuperCells = superCellToCells.size();

    GEOS_LOG_RANK_0( GEOS_FMT( "Initial distribution: {} super-cells across {} ranks (method={})",
                               numSuperCells, numRanks,
                               EnumStrings< ScatterMethod >::toString( scatterMethod ) ) );

    // Build a virtual mesh with one VTK_VERTEX cell per super-cell, placed at the
    // super-cell's centroid (average over all points of all member cells). This lets
    // computeCellRanks() see super-cells as atomic and partition them as such.
    vtkNew< vtkPoints > atomPoints;
    atomPoints->SetNumberOfPoints( numSuperCells );

    vtkNew< vtkUnstructuredGrid > atomMesh;
    atomMesh->SetPoints( atomPoints );
    atomMesh->Allocate( numSuperCells );

    stdVector< integer > cellToAtom( numCells );
    vtkPoints * meshPoints = cells3D->GetPoints();

    vtkIdType atomIdx = 0;
    for( auto const & kv : superCellToCells )
    {
      stdVector< vtkIdType > const & cellIndices = kv.second;

      real64 cx = 0.0, cy = 0.0, cz = 0.0;
      vtkIdType totalPoints = 0;
      for( vtkIdType cellIdx : cellIndices )
      {
        vtkIdType npts;
        vtkIdType const * pts;
        cells3D->GetCellPoints( cellIdx, npts, pts );
        for( vtkIdType i = 0; i < npts; ++i )
        {
          real64 p[3];
          meshPoints->GetPoint( pts[i], p );
          cx += p[0];
          cy += p[1];
          cz += p[2];
        }
        totalPoints += npts;
        cellToAtom[cellIdx] = static_cast< integer >( atomIdx );
      }

      real64 const inv = (totalPoints > 0) ? 1.0 / totalPoints : 0.0;
      atomPoints->SetPoint( atomIdx, cx * inv, cy * inv, cz * inv );

      vtkIdType const pid = atomIdx;
      atomMesh->InsertNextCell( VTK_VERTEX, 1, &pid );

      ++atomIdx;
    }

    // One rank per super-cell, then expand to per-cell so atomicity is preserved.
    stdVector< integer > const atomRanks =
      computeCellRanks( scatterMethod, *atomMesh, cartesianPartitions, numRanks );

    cellRanks.resize( numCells );
    for( vtkIdType c = 0; c < numCells; ++c )
    {
      cellRanks[c] = atomRanks[ cellToAtom[c] ];
    }
  }

  // All ranks ship cells
  vtkSmartPointer< vtkUnstructuredGrid > result =
    scatterByRankAssignment( cells3D.Get(), std::move( cellRanks ), comm );

  if( rank == 0 )
  {
    cells3D = nullptr;
  }

  // SuperCellId must survive redistribution
  if( result->GetNumberOfCells() > 0 )
  {
    GEOS_ERROR_IF( !vtkIdTypeArray::SafeDownCast( result->GetCellData()->GetArray( "SuperCellId" ) ),
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

  stdMap< vtkIdType, stdVector< vtkIdType > > superCellToLocalCells;
  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType superCellId = superCellIdArray->GetValue( i );
    superCellToLocalCells.get_inserted( superCellId ).push_back( i );
  }

  localIndex numLocalSuperCells = superCellToLocalCells.size();

  // Build super-cell global index distribution across ranks
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

  pmet_idx_t localGlobalStart = superElemDist[rank];
  pmet_idx_t localGlobalEnd = superElemDist[rank + 1];

  stdVector< vtkIdType > orderedSuperCellIds;
  orderedSuperCellIds.reserve( numLocalSuperCells );
  for( auto const & [superCellId, localCells] : superCellToLocalCells )
  {
    orderedSuperCellIds.push_back( superCellId );
  }
  std::sort( orderedSuperCellIds.begin(), orderedSuperCellIds.end() );

  stdMap< vtkIdType, pmet_idx_t > superCellIdToGlobalIdx;
  for( pmet_idx_t i = 0; i < numLocalSuperCells; ++i )
  {
    superCellIdToGlobalIdx.try_emplace( orderedSuperCellIds[i], localGlobalStart + i );
  }

  // -----------------------------------------------------------------------
  // Step 1: Build mappings for exchange
  // -----------------------------------------------------------------------
  // Local maps
  stdUnorderedMap< pmet_idx_t, vtkIdType > localGlobalCellToSuperCell;
  localGlobalCellToSuperCell.reserve( numLocalCells );

  stdUnorderedMap< pmet_idx_t, vtkIdType > localParmetisToVtk;
  localParmetisToVtk.reserve( numLocalCells );

  pmet_idx_t localParmetisStart = baseElemDist[rank];
  pmet_idx_t localParmetisEnd = baseElemDist[rank + 1];

  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    vtkIdType vtkGlobalId = cellGlobalIds->GetValue( i );
    vtkIdType superCellId = superCellIdArray->GetValue( i );

    localGlobalCellToSuperCell.try_emplace( vtkGlobalId, superCellId );
    localParmetisToVtk.try_emplace( localParmetisStart + i, vtkGlobalId );
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
      if( nbrIdx < localParmetisStart || nbrIdx >= localParmetisEnd )
      {
        ghostParmetisSet.insert( nbrIdx );
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 3: Exchange ghost mappings
  // -----------------------------------------------------------------------
  // Prepare local data to send
  array1d< pmet_idx_t > localParmetisIndices;
  array1d< vtkIdType > localVtkIds;
  array1d< vtkIdType > localSuperCellIds;

  localParmetisIndices.reserve( numLocalCells );
  localVtkIds.reserve( numLocalCells );
  localSuperCellIds.reserve( numLocalCells );

  for( vtkIdType i = 0; i < numLocalCells; ++i )
  {
    localParmetisIndices.emplace_back( localParmetisStart + i );
    localVtkIds.emplace_back( cellGlobalIds->GetValue( i ) );
    localSuperCellIds.emplace_back( superCellIdArray->GetValue( i ) );
  }

  // Gather all mappings
  array1d< pmet_idx_t > allParmetisIndices;
  array1d< vtkIdType > allVtkIds;
  array1d< vtkIdType > allSuperCellIds;

  MpiWrapper::allGatherv( localParmetisIndices.toViewConst(), allParmetisIndices, comm );
  MpiWrapper::allGatherv( localVtkIds.toViewConst(), allVtkIds, comm );
  MpiWrapper::allGatherv( localSuperCellIds.toViewConst(), allSuperCellIds, comm );

  // Build lookup tables from gathered data (only for ghosts + local)
  stdUnorderedMap< pmet_idx_t, vtkIdType > parmetisToVtk;
  stdUnorderedMap< vtkIdType, vtkIdType > vtkToSuperCell;

  parmetisToVtk.reserve( ghostParmetisSet.size() + numLocalCells );
  vtkToSuperCell.reserve( ghostParmetisSet.size() + numLocalCells );

  for( localIndex i = 0; i < allParmetisIndices.size(); ++i )
  {
    pmet_idx_t const parmetisIdx = allParmetisIndices[i];

    // Only store if it's local or a ghost we need
    if( (parmetisIdx >= localParmetisStart && parmetisIdx < localParmetisEnd) ||
        ghostParmetisSet.count( parmetisIdx ) > 0 )
    {
      vtkIdType const vtkId = allVtkIds[i];
      vtkIdType const superCellId = allSuperCellIds[i];

      parmetisToVtk.try_emplace( parmetisIdx, vtkId );
      vtkToSuperCell.try_emplace( vtkId, superCellId );
    }
  }

  // -----------------------------------------------------------------------
  // Step 4: Exchange super-cell global indices
  // -----------------------------------------------------------------------
  array1d< vtkIdType > sendSCIds;
  array1d< pmet_idx_t > sendSCGlobalIndices;

  sendSCIds.reserve( superCellIdToGlobalIdx.size() );
  sendSCGlobalIndices.reserve( superCellIdToGlobalIdx.size() );

  for( auto const & [scId, gIdx] : superCellIdToGlobalIdx )
  {
    sendSCIds.emplace_back( scId );
    sendSCGlobalIndices.emplace_back( gIdx );
  }

  array1d< vtkIdType > allSCIds;
  array1d< pmet_idx_t > allSCGlobalIndices;

  MpiWrapper::allGatherv( sendSCIds.toViewConst(), allSCIds, comm );
  MpiWrapper::allGatherv( sendSCGlobalIndices.toViewConst(), allSCGlobalIndices, comm );

  // Update local map with all super-cell mappings
  for( localIndex i = 0; i < allSCIds.size(); ++i )
  {
    superCellIdToGlobalIdx.insert_or_assign( allSCIds[i], allSCGlobalIndices[i] );
  }

  // -----------------------------------------------------------------------
  // Step 5: Build super-cell graph edges
  // -----------------------------------------------------------------------
  array1d< pmet_idx_t > superVertexWeights( numLocalSuperCells );
  stdVector< std::set< pmet_idx_t > > neighborSets( numLocalSuperCells );

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
  stdVector< pmet_idx_t > localEdgesSrc;
  stdVector< pmet_idx_t > localEdgesDst;

  for( localIndex i = 0; i < numLocalSuperCells; ++i )
  {
    pmet_idx_t globalI = localGlobalStart + i;
    for( pmet_idx_t nbrIdx : neighborSets[i] )
    {
      if( nbrIdx < localGlobalStart || nbrIdx >= localGlobalEnd )
      {
        localEdgesSrc.push_back( nbrIdx );  // destination of reverse edge
        localEdgesDst.push_back( globalI ); // source of reverse edge
      }
    }
  }

  int localEdgeCount = static_cast< int >( localEdgesSrc.size() );

  // Gather counts from all ranks
  array1d< int > allEdgeCounts( numRanks );
  MpiWrapper::allgather( &localEdgeCount, 1, allEdgeCounts.data(), 1, comm );

  // Compute displacements for symmetrization
  array1d< int > edgeDispls( numRanks + 1 );
  edgeDispls[0] = 0;
  std::partial_sum( allEdgeCounts.begin(), allEdgeCounts.end(),
                    edgeDispls.begin() + 1 );

  int totalEdges = edgeDispls[numRanks];

  // Gather all edges
  array1d< pmet_idx_t > allEdgeSrc( totalEdges > 0 ? totalEdges : 1 );
  array1d< pmet_idx_t > allEdgeDst( totalEdges > 0 ? totalEdges : 1 );

  if( localEdgeCount > 0 )
  {
    MpiWrapper::allgatherv( localEdgesSrc.data(), localEdgeCount,
                            allEdgeSrc.data(), allEdgeCounts.data(), edgeDispls.data(), comm );
    MpiWrapper::allgatherv( localEdgesDst.data(), localEdgeCount,
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

  pmet_idx_t localStart = superElemDist[rank];
  pmet_idx_t localEnd = superElemDist[rank + 1];

  // -----------------------------------------------------------------------
  // Check 1: No self-loops
  // -----------------------------------------------------------------------
  {
    for( localIndex i = 0; i < superGraph.size(); ++i )
    {
      pmet_idx_t localGlobalId = localStart + i;

      for( localIndex j = 0; j < superGraph.sizeOfArray( i ); ++j )
      {
        GEOS_ERROR_IF( superGraph[i][j] == localGlobalId,
                       GEOS_FMT( "Rank {}: Super-cell {} (global {}) has self-loop",
                                 rank, i, localGlobalId ) );
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
          GEOS_LOG_RANK( GEOS_FMT( "WARNING: Super-cell {} (global {}) has no neighbors (isolated)",
                                   i, localStart + i ) );
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
  // Check 6: Graph symmetry for local edges
  // -----------------------------------------------------------------------
  {
    // Pack a pair of (<32-bit) ids into a single uint64_t hash key.
    auto const edgeKey = []( pmet_idx_t a, pmet_idx_t b ) -> uint64_t
    {
      return (static_cast< uint64_t >( a ) << 32) | static_cast< uint64_t >( b );
    };

    // Iterate every edge whose endpoints are both local to this rank.
    auto const forEachLocalEdge = [&]( auto && fn )
    {
      for( localIndex i = 0; i < superGraph.size(); ++i )
      {
        pmet_idx_t const globalI = localStart + i;
        auto neighbors = superGraph[i];
        for( localIndex j = 0; j < neighbors.size(); ++j )
        {
          pmet_idx_t const globalJ = neighbors[j];
          if( globalJ >= localStart && globalJ < localEnd )
          {
            fn( globalI, globalJ );
          }
        }
      }
    };

    std::unordered_set< uint64_t > localOutgoingNeighbors;

    forEachLocalEdge( [&]( pmet_idx_t globalI, pmet_idx_t globalJ )
    {
      localOutgoingNeighbors.insert( edgeKey( globalI, globalJ ) );
    } );

    forEachLocalEdge( [&]( pmet_idx_t globalI, pmet_idx_t globalJ )
    {
      GEOS_ERROR_IF( localOutgoingNeighbors.find( edgeKey( globalJ, globalI ) ) == localOutgoingNeighbors.end(),
                     GEOS_FMT( "Rank {}: Asymmetric edge ({} -> {}) - reverse edge missing",
                               rank, globalI, globalJ ) );
    } );
  }

}


// =============================================================================================
// SECTION 5: PARTITION EXPANSION (after ParMETIS)
// =============================================================================================
array1d< int64_t >
expandSuperCellPartitioningToCells(
  vtkSmartPointer< vtkUnstructuredGrid > cells3D,
  array1d< int64_t > const & superPartitioning,
  stdMap< vtkIdType, localIndex > const & superCellIdToLocalIdx,
  MPI_Comm comm )
{
  int const rank = MpiWrapper::commRank( comm );

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
  vtkIdType numSplitSuperCells = 0;

  for( vtkIdType i = 0; i < numCells; ++i )
  {
    vtkIdType const scId = superCellIdArray->GetValue( i );
    auto & ranks = superCellToRanks.get_inserted( scId );
    // A super-cell becomes "split" the first time we add a second distinct rank.
    numSplitSuperCells += ( ranks.insert( cellPartitioning[i] ).second && ranks.size() == 2 );
  }

  vtkIdType totalSplitSuperCells = MpiWrapper::sum( numSplitSuperCells, comm );

  GEOS_ERROR_IF( totalSplitSuperCells > 0,
                 GEOS_FMT( "Partitioning failed: {:L} super-cells were split across ranks!",
                           totalSplitSuperCells ) );

  return cellPartitioning;
}

} // namespace vtk
} // namespace geos
