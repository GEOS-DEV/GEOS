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
 * @file SemistructuredPartitioner.cpp
 */

#include "SemistructuredPartitioner.hpp"

#include "linearAlgebra/multiscale/mesh/MeshUtils.hpp"
#include "linearAlgebra/multiscale/mesh/MeshData.hpp"
#include "common/MpiWrapper.hpp"

#ifdef GEOSX_USE_METIS
#include "linearAlgebra/multiscale/mesh/coarsening/MetisInterface.hpp"
#endif

#ifdef GEOSX_USE_SCOTCH
#include "linearAlgebra/multiscale/mesh/coarsening/ScotchInterface.hpp"
#endif

namespace geos
{

namespace multiscale
{

template< typename FUNC >
CRSMatrix< int64_t, int64_t, int64_t >
buildLayerGraph( multiscale::MeshLevel const & mesh,
                 arrayView1d< localIndex const > const & layerCells,
                 std::unordered_map< integer, localIndex > const & structIndexToLayerCell,
                 integer const layerIndex,
                 integer const minCommonNodes,
                 FUNC && weightFunc )
{
  localIndex const numCells = layerCells.size();
  arrayView1d< integer const > const cellGhostRank = mesh.cellManager().ghostRank().toViewConst();
  MeshObjectManager::MapViewConst const cellToNode = mesh.cellManager().toDualRelation().toViewConst();
  MeshObjectManager::MapViewConst const nodeToCell = mesh.nodeManager().toDualRelation().toViewConst();

  auto const structIndex = mesh.cellManager().getField< fields::StructuredIndex >().toViewConst();

  // Count exact length of each row (METIS does not accept holes in CRS graph)
  array1d< int64_t > rowCounts( numCells );
  forAll< parallelHostPolicy >( numCells, [=, rowCounts = rowCounts.toView()] ( localIndex const i )
  {
    localIndex const cellIdx = layerCells[i];
    GEOS_ASSERT_EQ( structIndex[cellIdx][1], layerIndex );
    localIndex numUniqueNeighborCells = 0;
    meshUtils::forUniqueNeighbors< 512 >( cellIdx, cellToNode, nodeToCell,
                                          [&]( localIndex const nbrIdx,
                                               localIndex const numCommonNodes )
    {
      if( numCommonNodes >= minCommonNodes && nbrIdx != cellIdx && cellGhostRank[nbrIdx] < 0 && structIndex[nbrIdx][1] == layerIndex )
      {
        ++numUniqueNeighborCells;
      }
    } );
    rowCounts[i] = numUniqueNeighborCells;
  } );

  CRSMatrix< int64_t, int64_t, int64_t > graph;
  graph.resizeFromRowCapacities< parallelHostPolicy >( numCells, numCells, rowCounts.data() );

  // Fill the graph
  forAll< parallelHostPolicy >( numCells, [=, rowCounts = rowCounts.toView(), graph = graph.toView()] ( localIndex const i )
  {
    localIndex const cellIdx = layerCells[i];
    meshUtils::forUniqueNeighbors< 512 >( cellIdx, cellToNode, nodeToCell, [&]( localIndex const nbrIdx, localIndex const numCommonNodes )
    {
      if( numCommonNodes >= minCommonNodes && nbrIdx != cellIdx && cellGhostRank[nbrIdx] < 0 && structIndex[nbrIdx][1] == layerIndex )
      {
        graph.insertNonZero( LvArray::integerConversion< int64_t >( i ),
                             LvArray::integerConversion< int64_t >( structIndexToLayerCell.at( structIndex[nbrIdx][0] ) ),
                             weightFunc( cellIdx, nbrIdx ) );
      }
    } );
  } );

  return graph;
}

localIndex SemistructuredPartitioner::generate( MeshLevel const & mesh,
                                                arrayView1d< localIndex > const & partition )
{
  GEOS_MARK_FUNCTION;

  MeshObjectManager const & cellManager = mesh.cellManager();
  localIndex const numCells = cellManager.numOwnedObjects();

  auto const structIndex = cellManager.getField< fields::StructuredIndex >().toViewConst();
  GEOS_ERROR_IF_LT_MSG( structIndex.size( 1 ), 2, "Not enough structured indices" );

  RAJA::ReduceMin< parallelHostReduce, integer > loZIndex( std::numeric_limits< integer >::max() );
  RAJA::ReduceMax< parallelHostReduce, integer > hiZIndex( std::numeric_limits< integer >::min() );

  forAll< parallelHostPolicy >( numCells, [=]( localIndex const i )
  {
    GEOS_ASSERT_GE( structIndex[i][0], 0 );
    GEOS_ASSERT_GE( structIndex[i][1], 0 );
    loZIndex.min( structIndex[i][1] );
    hiZIndex.max( structIndex[i][1] );
  } );

  // Special treatment for ranks that don't have a piece of the mesh
  if( numCells == 0 )
  {
    loZIndex = 0;
    hiZIndex = -1;
  }

  integer const numCellsZ = hiZIndex - loZIndex + 1;
  integer const numCellsA = numCells / numCellsZ;
  GEOS_ASSERT_EQ( numCells % numCellsZ, 0 );

  // When semi-coarsening is enabled, only coarsen in z-direction,
  // until every rank has reduced its mesh to 1 cell in z-direction.
  if( m_params.structured.semicoarsening && MpiWrapper::max( numCellsZ ) > 1 )
  {
    m_params.ratio[0] = 1.0;
  }

  // Compute coarse grid sizes
  m_numPart[0] = static_cast< integer >( std::ceil( numCellsA / m_params.ratio[0] ) ); // rounded up
  m_numPart[1] = static_cast< integer >( std::ceil( numCellsZ / m_params.ratio[1] ) ); // rounded up

  if( m_numPart[1] == 0 )
  {
    return 0; // exit to avoid metis call on empty graph
  }

  // Compute (can/should we parallelize this?):
  // - a list of cells in one specific layer
  // - a map of area index to index in the above list
  integer const targetLayer = loZIndex;
  arrayView1d< integer const > const ghostRank = cellManager.ghostRank();

  array1d< localIndex > layerCells;
  layerCells.reserve( numCellsA );

  std::unordered_map< integer, localIndex > indexToLayer;
  indexToLayer.reserve( numCellsA );

  forAll< serialPolicy >( numCells, [&]( localIndex const i )
  {
    if( structIndex[i][1] == targetLayer && ghostRank[i] < 0 )
    {
      indexToLayer.emplace( structIndex[i][0], layerCells.size() );
      layerCells.emplace_back( i );
    }
  } );

  // Special handling for another trivial case (to avoid METIS problems)
  array1d< int64_t > areaPart( numCellsA );
  if( m_numPart[0] == numCellsA )
  {
    areaPart.setValues< parallelHostPolicy >( 1 );
    RAJA::exclusive_scan_inplace< parallelHostPolicy >( areaPart, RAJA::operators::plus< localIndex >{} );
  }
  else
  {
    // Compute a graph for area partitioning
    CRSMatrix< int64_t, int64_t, int64_t > const graph = [&]()
    {
      if( m_params.graph.preserveRegions && mesh.cellManager().hasWrapper( fields::multiscale::OrigElementRegion::key() ) )
      {
        arrayView1d< localIndex const > const region = mesh.cellManager().getField< fields::multiscale::OrigElementRegion >();
        auto weight = [region]( auto const i, auto const j ){ return 1.0 + 1000.0 * ( region[i] == region[j] ); };
        return buildLayerGraph( mesh, layerCells, indexToLayer, targetLayer, m_params.graph.minCommonNodes, weight );
      }
      else
      {
        auto weight = []( localIndex, localIndex ){ return 1.0; };
        return buildLayerGraph( mesh, layerCells, indexToLayer, targetLayer, m_params.graph.minCommonNodes, weight );
      }
    }();

    switch( m_params.graph.method )
    {
      using Method = LinearSolverParameters::Multiscale::Coarsening::Graph::Method;
      case Method::metis:
      {
#ifdef GEOSX_USE_METIS
        metis::partition( graph.toViewConst(), m_params.graph.metis, m_numPart[0], areaPart.toView() );
#else
        GEOS_THROW( "Attempted to use METIS partition, but GEOSX is not built with METIS support.", InputError );
#endif
        break;
      }
      case Method::scotch:
      {
#ifdef GEOSX_USE_METIS
        scotch::partition( graph.toViewConst(), m_params.graph.scotch, m_numPart[0], areaPart.toView() );
#else
        GEOS_THROW( "Attempted to use Scotch partition, but GEOSX is not built with Scotch support.", InputError );
#endif
        break;
      }
      default:
      {
        GEOS_ERROR( "Unrecognized graph partitioner type" );
      }
    }
  }

  integer const zRatio = ( numCellsZ + m_numPart[1] - 1 ) / m_numPart[1]; // rounded up

  // Compute cartesian coarse cell indices
  forAll< parallelHostPolicy >( numCells, [&loZIndex, &indexToLayer, numPartZ = m_numPart[1],
                                           structIndex, zRatio, partition,
                                           layerCells = layerCells.toViewConst(),
                                           areaPart = areaPart.toViewConst()]( localIndex const i )
  {
    localIndex const partIdxA = LvArray::integerConversion< localIndex >( areaPart[indexToLayer.at( structIndex[i][0] )] );
    localIndex const partIdxZ = ( structIndex[i][1] - loZIndex ) / zRatio;
    partition[i] = partIdxA * numPartZ + partIdxZ;
  } );

  return m_numPart[0] * m_numPart[1]; // TODO
}

void SemistructuredPartitioner::setCoarseData( MeshLevel & coarseMesh ) const
{
  array2d< integer > & structIndex =
    coarseMesh.cellManager().registerField< fields::StructuredIndex >( {} ).reference();
  structIndex.resizeDimension< 1 >( 2 );

  forAll< parallelHostPolicy >( coarseMesh.cellManager().numOwnedObjects(),
                                [structIndex = structIndex.toView(), numPartZ = m_numPart[1]]( localIndex const i )
  {
    structIndex[i][0] = LvArray::integerConversion< integer >( i / numPartZ );
    structIndex[i][1] = LvArray::integerConversion< integer >( i % numPartZ );
  } );
}


} // namespace multiscale
} // geosx
