/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Coarsening.cpp
 */

#include "Coarsening.hpp"

#include "linearAlgebra/multiscale/mesh/MeshData.hpp"
#include "linearAlgebra/multiscale/mesh/MeshLevel.hpp"
#include "linearAlgebra/multiscale/mesh/MeshUtils.hpp"
#include "PartitionerBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{
namespace multiscale
{
namespace coarsening
{

namespace
{

void fillBasicCellData( MeshObjectManager const & fineCellManager,
                        MeshObjectManager & coarseCellManager )
{
  arrayView1d< localIndex const > const coarseCellIndex = fineCellManager.getField< fields::multiscale::CoarseCellLocalIndex >().toViewConst();
  arrayView1d< integer const > const fineGhostRank = fineCellManager.ghostRank();
  arrayView1d< integer const > const fineIsExternal = fineCellManager.isExternal();
  arrayView1d< integer const > const fineDomainBoundary = fineCellManager.getDomainBoundaryIndicator();

  arrayView1d< integer > const coarseGhostRank = coarseCellManager.ghostRank();
  arrayView1d< integer > const coarseIsExternal = coarseCellManager.isExternal();
  arrayView1d< integer > const coarseDomainBoundary = coarseCellManager.getDomainBoundaryIndicator();

  forAll< parallelHostPolicy >( fineCellManager.size(), [=]( localIndex const icf )
  {
    localIndex const icc = coarseCellIndex[icf];
    RAJA::atomicMax( parallelHostAtomic{}, &coarseGhostRank[icc], fineGhostRank[icf] );
    RAJA::atomicMax( parallelHostAtomic{}, &coarseIsExternal[icc], fineIsExternal[icf] );
    RAJA::atomicMax( parallelHostAtomic{}, &coarseDomainBoundary[icc], fineDomainBoundary[icf] );
  } );
}

void buildCellLocalToGlobalMaps( std::set< globalIndex > const & ghostGlobalIndices,
                                 globalIndex const rankOffset,
                                 MeshObjectManager & coarseCellManager )
{
  arrayView1d< globalIndex > const coarseLocalToGlobal = coarseCellManager.localToGlobalMap();
  {
    localIndex icc = 0;
    for(; icc < coarseCellManager.numOwnedObjects(); ++icc )
    {
      coarseLocalToGlobal[icc] = rankOffset + icc;
    }
    for( globalIndex const coarseGlobalIndex : ghostGlobalIndices )
    {
      coarseLocalToGlobal[icc++] = coarseGlobalIndex;
    }
  }
  coarseCellManager.constructGlobalToLocalMap();
  coarseCellManager.setMaxGlobalIndex();
}

void fillBasicNodeData( MeshObjectManager const & fineNodeManager,
                        MeshObjectManager & coarseNodeManager )
{
  arrayView1d< localIndex const > const fineNodeIndex = coarseNodeManager.getField< fields::multiscale::FineNodeLocalIndex >().toViewConst();

  meshUtils::fillArrayBySrcIndex< parallelHostPolicy >( fineNodeIndex,
                                                        coarseNodeManager.ghostRank(),
                                                        fineNodeManager.ghostRank() );
  meshUtils::fillArrayBySrcIndex< parallelHostPolicy >( fineNodeIndex,
                                                        coarseNodeManager.isExternal(),
                                                        fineNodeManager.isExternal() );
  meshUtils::fillArrayBySrcIndex< parallelHostPolicy >( fineNodeIndex,
                                                        coarseNodeManager.getDomainBoundaryIndicator(),
                                                        fineNodeManager.getDomainBoundaryIndicator() );
}

void buildNodeToCellMap( MeshObjectManager const & fineCellManager,
                         MeshObjectManager const & fineNodeManager,
                         MeshObjectManager & coarseNodeManager )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< localIndex const > const coarseCellLocalIndex = fineCellManager.getField< fields::multiscale::CoarseCellLocalIndex >();
  arrayView1d< localIndex const > const fineNodeLocalIndex = coarseNodeManager.getField< fields::multiscale::FineNodeLocalIndex >();
  MeshObjectManager::MapViewConst const fineNodeToCell = fineNodeManager.toDualRelation().toViewConst();

  // First pass: count length of each sub-array in order to do exact allocation
  array1d< localIndex > rowCounts( coarseNodeManager.size() );
  forAll< parallelHostPolicy >( coarseNodeManager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const inc )
  {
    localIndex count = 0;
    meshUtils::forUniqueNeighborValues< 256 >( fineNodeLocalIndex[inc],
                                               fineNodeToCell,
                                               coarseCellLocalIndex,
                                               [&]( localIndex const )
    {
      ++count;
    } );
    rowCounts[inc] = count;
  } );

  // Resize the map
  MeshObjectManager::MapType & coarseNodeToCell = coarseNodeManager.toDualRelation();
  coarseNodeToCell.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Second pass: fill the map
  forAll< parallelHostPolicy >( coarseNodeManager.size(), [=, coarseNodeToCell = coarseNodeToCell.toView()]( localIndex const inc )
  {
    meshUtils::forUniqueNeighborValues< 256 >( fineNodeLocalIndex[inc],
                                               fineNodeToCell,
                                               coarseCellLocalIndex,
                                               [&]( localIndex const icc )
    {
      coarseNodeToCell.insertIntoSet( inc, icc );
    } );
  } );
}

void buildCellToNodeMap( MeshObjectManager const & coarseNodeManager,
                         MeshObjectManager & coarseCellManager )
{
  GEOS_MARK_FUNCTION;

  MeshObjectManager::MapViewConst const nodeToCell = coarseNodeManager.toDualRelation().toViewConst();
  MeshObjectManager::MapType & cellToNode = coarseCellManager.toDualRelation();

  // First pass: count the row lengths in transpose map
  array1d< localIndex > rowCounts( cellToNode.size() );
  forAll< parallelHostPolicy >( nodeToCell.size(), [=, rowCounts = rowCounts.toView()]( localIndex const inc )
  {
    for( localIndex const icc : nodeToCell[inc] )
    {
      RAJA::atomicInc< parallelHostAtomic >( &rowCounts[icc] );
    }
  } );

  // Create and resize the temporary map
  ArrayOfArrays< localIndex > cellToNodeTemp;
  cellToNodeTemp.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Second pass: fill the map
  forAll< parallelHostPolicy >( nodeToCell.size(), [=, cellToNode = cellToNodeTemp.toView()]( localIndex const inc )
  {
    for( localIndex const icc : nodeToCell[inc] )
    {
      cellToNode.emplaceBackAtomic< parallelHostAtomic >( icc, inc );
    }
  } );

  // Move the temp map and sort the entries
  cellToNode.assimilate< parallelHostPolicy >( std::move( cellToNodeTemp ),
                                               LvArray::sortedArrayManipulation::UNSORTED_NO_DUPLICATES );
}

void buildCoarseCells( multiscale::MeshLevel & fineMesh,
                       multiscale::MeshLevel & coarseMesh,
                       LinearSolverParameters::Multiscale::Coarsening const & params )
{
  GEOS_MARK_FUNCTION;

  MeshObjectManager & fineCellManager = fineMesh.cellManager();
  MeshObjectManager & coarseCellManager = coarseMesh.cellManager();

  // Allocate arrays to hold partitioning info (local and global)
  arrayView1d< localIndex > const coarseLocalIndex =
    fineCellManager.registerField< fields::multiscale::CoarseCellLocalIndex >( coarseMesh.name() ).referenceAsView();
  arrayView1d< globalIndex > const coarseGlobalIndex =
    fineCellManager.registerField< fields::multiscale::CoarseCellGlobalIndex >( coarseMesh.name() ).referenceAsView();

  // Generate the partitioning locally
  std::unique_ptr< PartitionerBase > partitioner = PartitionerBase::create( params );
  localIndex const numLocalCoarseCells = partitioner->generate( fineMesh, coarseLocalIndex );

  // Compute global number of partitions
  globalIndex const rankOffset = MpiWrapper::prefixSum< globalIndex >( numLocalCoarseCells );

  // Fill in partition global index for locally owned cells
  forAll< parallelHostPolicy >( fineCellManager.numOwnedObjects(),
                                [coarseLocalIndex = coarseLocalIndex.toViewConst(),
                                 coarseGlobalIndex, rankOffset]( localIndex const icf )
  {
    coarseGlobalIndex[icf] = coarseLocalIndex[icf] + rankOffset;
  } );

  // Synchronize partition global index across ranks
  string_array fieldNames;
  fieldNames.emplace_back( fields::multiscale::CoarseCellGlobalIndex::key() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames, fineCellManager, fineMesh.domain()->getNeighbors(), false );

  // Scan ghosted cells and collect all new partition global indices
  std::set< globalIndex > ghostGlobalIndices;
  for( localIndex icf = fineCellManager.numOwnedObjects(); icf < fineCellManager.size(); ++icf )
  {
    ghostGlobalIndices.insert( coarseGlobalIndex[icf] );
  }
  localIndex const numPresentCoarseCells = numLocalCoarseCells + LvArray::integerConversion< localIndex >( ghostGlobalIndices.size() );

  // Resize and start populating coarse cell manager
  coarseCellManager.resize( numPresentCoarseCells );
  coarseCellManager.setNumOwnedObjects( numLocalCoarseCells );

  // Compute cartesian indices for the coarse cells
  partitioner->setCoarseData( coarseMesh );

  // Populate coarse local-global maps
  buildCellLocalToGlobalMaps( ghostGlobalIndices, rankOffset, coarseCellManager );

  // finish filling the partition array for ghosted fine cells
  forRange< parallelHostPolicy >( fineCellManager.numOwnedObjects(), fineCellManager.size(),
                                  [coarseGlobalIndex = coarseGlobalIndex.toViewConst(),
                                   coarseLocalIndex, &coarseCellManager]( localIndex const icf )
  {
    coarseLocalIndex[icf] = coarseCellManager.globalToLocalMap( coarseGlobalIndex[icf] );
  } );

  fillBasicCellData( fineCellManager, coarseCellManager );

  // Populate neighbor data and sets
  std::vector< int > neighborRanks = fineMesh.domain()->getNeighborRanks();
  for( int const rank : neighborRanks )
  {
    coarseCellManager.addNeighbor( rank );
  }
  meshUtils::copySets( fineCellManager,
                       fields::multiscale::CoarseCellLocalIndex::key(),
                       coarseCellManager );
  meshUtils::copyNeighborData( fineCellManager,
                               fields::multiscale::CoarseCellLocalIndex::key(),
                               neighborRanks,
                               coarseCellManager,
                               meshUtils::mapIndexArrayUnique< localIndex > );
}

array1d< localIndex >
findCoarseNodes( DomainPartition & domain,
                 MeshObjectManager & fineNodeManager,
                 MeshObjectManager const & fineCellManager,
                 arrayView1d< string const > const & boundaryNodeSets )
{
  // Build and sync an adjacency map of fine nodes to coarse subdomains (including global boundaries)
  ArrayOfSets< globalIndex > nodeToSubdomainLocal =
    meshUtils::buildFineObjectToSubdomainMap( fineNodeManager,
                                              fineCellManager.getField< fields::multiscale::CoarseCellGlobalIndex >().toViewConst(),
                                              boundaryNodeSets );

  ArrayOfArrays< globalIndex > nodeToSubdomainArray;
  nodeToSubdomainArray.assimilate( std::move( nodeToSubdomainLocal ) );

  {
    meshUtils::ScopedDataRegistrar reg( fineNodeManager, "nodeToCoarseSubdomain", nodeToSubdomainArray );
    reg.sync( domain );
  }

  ArrayOfSets< globalIndex > nodeToSubdomain;
  nodeToSubdomain.assimilate< parallelHostPolicy >( std::move( nodeToSubdomainArray ),
                                                    LvArray::sortedArrayManipulation::SORTED_UNIQUE );

  array1d< localIndex > candidates =
    meshUtils::findCoarseNodesByDualPartition( fineNodeManager.toDualRelation().toViewConst(),
                                               fineCellManager.toDualRelation().toViewConst(),
                                               nodeToSubdomain.toViewConst(), 3 );

  // In some rare cases, ranks may still disagree on which nodes are identified as coarse.
  // In this case we let the owner of the node decide and exchange decisions.
  // This is unfortunate, since purely local analysis works in most cases.

  array1d< integer > isCoarseNode( fineNodeManager.size() );
  forAll< parallelHostPolicy >( candidates.size(), [candidates = candidates.toViewConst(),
                                                    isCoarseNode = isCoarseNode.toView()]( localIndex const i )
  {
    isCoarseNode[candidates[i]] = 1;
  } );

  {
    meshUtils::ScopedDataRegistrar reg( fineNodeManager, "isCoarseNode", isCoarseNode );
    reg.sync( domain );
  }

  array1d< localIndex > coarseNodes;
  coarseNodes.reserve( candidates.size() );
  for( localIndex i = 0; i < isCoarseNode.size(); ++i )
  {
    if( isCoarseNode[i] == 1 )
    {
      coarseNodes.emplace_back( i );
    }
  }
  return coarseNodes;
}

void buildCoarseNodes( multiscale::MeshLevel & fineMesh,
                       multiscale::MeshLevel & coarseMesh,
                       arrayView1d< string const > const & boundaryNodeSets )
{
  GEOS_MARK_FUNCTION;

  MeshObjectManager & fineNodeManager = fineMesh.nodeManager();
  MeshObjectManager & fineCellManager = fineMesh.cellManager();
  MeshObjectManager & coarseNodeManager = coarseMesh.nodeManager();

  // Find all locally present coarse nodes
  array1d< localIndex > const coarseNodes = findCoarseNodes( *fineMesh.domain(),
                                                             fineNodeManager,
                                                             fineCellManager,
                                                             boundaryNodeSets );

  // Reorder them to have all local nodes precede ghosted (stable partition to preserve order)
  arrayView1d< integer const > const nodeGhostRank = fineNodeManager.ghostRank();
  auto const localEnd = std::stable_partition( coarseNodes.begin(), coarseNodes.end(),
                                               [=]( localIndex const inf ){ return nodeGhostRank[inf] < 0; } );

  localIndex const numLocalCoarseNodes =
    LvArray::integerConversion< localIndex >( std::distance( coarseNodes.begin(), localEnd ) );
  globalIndex const firstLocalNodeIndex = MpiWrapper::prefixSum< globalIndex >( numLocalCoarseNodes );

  coarseNodeManager.resize( coarseNodes.size() );
  coarseNodeManager.setNumOwnedObjects( numLocalCoarseNodes );

  // Finally, build coarse node maps
  arrayView1d< localIndex > const coarseNodeLocalIndex =
    fineNodeManager.registerField< fields::multiscale::CoarseNodeLocalIndex >( coarseMesh.name() ).reference().toView();
  arrayView1d< globalIndex > const coarseNodeGlobalIndex =
    fineNodeManager.registerField< fields::multiscale::CoarseNodeGlobalIndex >( coarseMesh.name() ).reference().toView();
  arrayView1d< localIndex > const fineNodeLocalIndex =
    coarseNodeManager.registerField< fields::multiscale::FineNodeLocalIndex >( coarseMesh.name() ).reference().toView();
  arrayView1d< globalIndex > const coarseNodeLocalToGlobal =
    coarseNodeManager.localToGlobalMap();

  // Fill the local part
  forAll< parallelHostPolicy >( numLocalCoarseNodes, [=, coarseNodes = coarseNodes.toViewConst()]( localIndex const i )
  {
    localIndex const inf = coarseNodes[i];
    coarseNodeLocalIndex[inf] = i;
    coarseNodeGlobalIndex[inf] = firstLocalNodeIndex + i;
    coarseNodeLocalToGlobal[i] = firstLocalNodeIndex + i;
    fineNodeLocalIndex[i] = inf;
  } );

  // Sync across ranks
  string_array fieldNames;
  fieldNames.emplace_back( fields::multiscale::CoarseNodeGlobalIndex::key() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames, fineNodeManager, fineMesh.domain()->getNeighbors(), false );

  // Fill the ghosted part
  forAll< parallelHostPolicy >( coarseNodes.size() - numLocalCoarseNodes,
                                [=, coarseNodes = coarseNodes.toViewConst()]( localIndex const k )
  {
    localIndex const i = k + numLocalCoarseNodes;
    localIndex const inf = coarseNodes[i];
    coarseNodeLocalIndex[inf] = i;
    coarseNodeLocalToGlobal[i] = coarseNodeGlobalIndex[inf];
    fineNodeLocalIndex[i] = inf;
  } );

  coarseNodeManager.constructGlobalToLocalMap();
  coarseNodeManager.setMaxGlobalIndex();
  fillBasicNodeData( fineNodeManager, coarseNodeManager );

  // Populate neighbor data and sets
  std::vector< int > neighborRanks = fineMesh.domain()->getNeighborRanks();
  for( int const rank : neighborRanks )
  {
    coarseNodeManager.addNeighbor( rank );
  }
  meshUtils::copySets( fineNodeManager,
                       fields::multiscale::CoarseNodeLocalIndex::key(),
                       coarseNodeManager );
  meshUtils::copyNeighborData( fineNodeManager,
                               fields::multiscale::CoarseNodeLocalIndex::key(),
                               neighborRanks,
                               coarseNodeManager,
                               meshUtils::mapIndexArray< localIndex > );
}

} // namespace

void buildCoarseMesh( MeshLevel & fineMesh,
                      MeshLevel & coarseMesh,
                      LinearSolverParameters::Multiscale::Coarsening const & params,
                      array1d< string > const & boundaryNodeSets )
{
  buildCoarseCells( fineMesh, coarseMesh, params );
  buildCoarseNodes( fineMesh, coarseMesh, boundaryNodeSets );
  buildNodeToCellMap( fineMesh.cellManager(), fineMesh.nodeManager(), coarseMesh.nodeManager() );
  buildCellToNodeMap( coarseMesh.nodeManager(), coarseMesh.cellManager() );
}

} // namespace coarsening
} // namespace multiscale
} // namespace geos
