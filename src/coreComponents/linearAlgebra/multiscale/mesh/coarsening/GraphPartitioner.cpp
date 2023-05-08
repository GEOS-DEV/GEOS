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
 * @file GraphPartitioner.cpp
 */

#include "GraphPartitioner.hpp"

#include "linearAlgebra/multiscale/mesh/MeshUtils.hpp"
#include "linearAlgebra/multiscale/mesh/MeshData.hpp"

#ifdef GEOSX_USE_METIS
#include "linearAlgebra/multiscale/mesh/coarsening/MetisInterface.hpp"
#endif

#ifdef GEOSX_USE_SCOTCH
#include "linearAlgebra/multiscale/mesh/coarsening/ScotchInterface.hpp"
#endif

#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"

namespace geos
{
namespace multiscale
{

template< typename FUNC >
CRSMatrix< int64_t, int64_t, int64_t >
buildCellGraph( multiscale::MeshLevel const & mesh,
                integer const minCommonNodes,
                FUNC && weightFunc )
{
  localIndex const numCells = mesh.cellManager().numOwnedObjects();
  arrayView1d< integer const > const cellGhostRank = mesh.cellManager().ghostRank().toViewConst();
  MeshObjectManager::MapViewConst const cellToNode = mesh.cellManager().toDualRelation().toViewConst();
  MeshObjectManager::MapViewConst const nodeToCell = mesh.nodeManager().toDualRelation().toViewConst();

  // Count exact length of each row (METIS does not accept holes in CRS graph)
  array1d< int64_t > rowCounts( numCells );
  forAll< parallelHostPolicy >( numCells, [=, rowCounts = rowCounts.toView()] ( localIndex const cellIdx )
  {
    localIndex numUniqueNeighborCells = 0;
    meshUtils::forUniqueNeighbors< 512 >( cellIdx, cellToNode, nodeToCell,
                                          [&]( localIndex const nbrIdx,
                                               localIndex const numCommonNodes )
    {
      if( numCommonNodes >= minCommonNodes && nbrIdx != cellIdx && cellGhostRank[nbrIdx] < 0 )
      {
        ++numUniqueNeighborCells;
      }
    } );
    rowCounts[cellIdx] = numUniqueNeighborCells;
  } );

  CRSMatrix< int64_t, int64_t, int64_t > graph;
  graph.resizeFromRowCapacities< parallelHostPolicy >( numCells, numCells, rowCounts.data() );

  // Fill the graph
  forAll< parallelHostPolicy >( numCells, [=, rowCounts = rowCounts.toView(), graph = graph.toView()] ( localIndex const cellIdx )
  {
    meshUtils::forUniqueNeighbors< 512 >( cellIdx, cellToNode, nodeToCell, [&]( localIndex const nbrIdx, localIndex const numCommonNodes )
    {
      if( numCommonNodes >= minCommonNodes && nbrIdx != cellIdx && cellGhostRank[nbrIdx] < 0 )
      {
        graph.insertNonZero( LvArray::integerConversion< int64_t >( cellIdx ),
                             LvArray::integerConversion< int64_t >( nbrIdx ),
                             weightFunc( cellIdx, nbrIdx ) );
      }
    } );
  } );

  return graph;
}

static CRSMatrix< int64_t, int64_t, int64_t >
buildCellGraphFromMatrix( CRSMatrixView< real64 const, localIndex const > const & mat,
                          integer const multiplier )
{
  array1d< int64_t > rowSizes( mat.numRows() );
  array1d< real64 > diag( mat.numRows() );
  forAll< parallelHostPolicy >( mat.numRows(), [rowSizes = rowSizes.toView(),
                                                diag = diag.toView(), mat]( localIndex const i )
  {
    rowSizes[i] = mat.numNonZeros( i ) - 1;
    auto const columns = mat.getColumns( i );
    auto const values = mat.getEntries( i );
    for( localIndex k = 0; k < columns.size(); ++k )
    {
      if( columns[k] == i )
      {
        diag[i] = LvArray::math::abs( values[k] );
      }
    }
    GEOS_ASSERT_GT( diag[i], 0.0 );
  } );

  // Now build the symmetric weighted graph
  CRSMatrix< int64_t, int64_t, int64_t > graph;
  graph.resizeFromRowCapacities< parallelHostPolicy >( mat.numRows(), mat.numColumns(), rowSizes.data() );

  forAll< parallelHostPolicy >( mat.numRows(), [graph = graph.toView(), mat, diag, multiplier]( localIndex const i )
  {
    auto const icols = mat.getColumns( i );
    auto const ivals = mat.getEntries( i );
    for( localIndex k = 0; k < icols.size(); ++k )
    {
      localIndex const j = icols[k];
      real64 const aij = !isZero( ivals[k] ) ? ivals[k] : -diag[i] / ( icols.size() - 1 );
      if( j != i )
      {
        auto const jcols = mat.getColumns( j );
        auto const jvals = mat.getEntries( j );

        auto const d = LvArray::sortedArrayManipulation::find( jcols.dataIfContiguous(), mat.numNonZeros( j ), i );
        localIndex const m = LvArray::integerConversion< localIndex >( d );
        GEOS_ASSERT_GT( mat.numNonZeros( j ), m );

        real64 const aji = !isZero( jvals[m] ) ? jvals[m] : -diag[j] / ( jcols.size() - 1 );

        // w = gamma * 0.5 * |a_ij + a_ji| / sqrt(a_ii * a_jj)
        real64 const value = LvArray::math::abs( aij + aji ) / 2;
        real64 const normalizer = LvArray::math::sqrt( diag[i] * diag[j] );
        real64 const weight = multiplier * ( value / normalizer );

        graph.insertNonZero( i, j, static_cast< int64_t >( weight ) );
      }
    }
  } );

  return graph;
}

localIndex GraphPartitioner::generate( multiscale::MeshLevel const & mesh,
                                       arrayView1d< localIndex > const & partition )
{
  GEOS_MARK_FUNCTION;

  real64 ratio = 1.0;
  for( real64 const r : m_params.ratio )
  {
    ratio *= r;
  }

  localIndex const numCells = mesh.cellManager().numOwnedObjects();
  localIndex const numPart = static_cast< localIndex >( std::ceil( numCells / ratio ) ); // rounded up

  // Special handling for a trivial case
  if( numPart <= 1 )
  {
    partition.zero();
    return numPart;
  }

  // Special handling for another trivial case (to avoid METIS problems)
  if( numPart == numCells )
  {
    partition.setValues< parallelHostPolicy >( 1 );
    RAJA::exclusive_scan_inplace< parallelHostPolicy >( partition, RAJA::operators::plus< localIndex >{} );
    return numPart;
  }

  CRSMatrix< int64_t, int64_t, int64_t > const graph = [&]()
  {
    auto const & params = m_params.graph;
    if( params.matrixWeights > 0 && params.localMatrix != nullptr )
    {
      GEOS_ERROR_IF_NE_MSG( params.localMatrix->numRows(), numCells, "Invalid local matrix provided" );
      return buildCellGraphFromMatrix( params.localMatrix->toViewConst(), params.matrixWeights );
    }
    else if( m_params.graph.preserveRegions && mesh.cellManager().hasWrapper( fields::multiscale::OrigElementRegion::key() ) )
    {
      arrayView1d< localIndex const > const region = mesh.cellManager().getField< fields::multiscale::OrigElementRegion >();
      auto weight = [region]( auto const i, auto const j ){ return 1.0 + 1000.0 * ( region[i] == region[j] ); };
      return buildCellGraph( mesh, m_params.graph.minCommonNodes, weight );
    }
    else
    {
      auto weight = []( localIndex, localIndex ){ return 1.0; };
      return buildCellGraph( mesh, m_params.graph.minCommonNodes, weight );
    }
  }();

  array1d< int64_t > part( numCells );
  switch( m_params.graph.method )
  {
    using Method = LinearSolverParameters::Multiscale::Coarsening::Graph::Method;
    case Method::metis:
    {
#ifdef GEOSX_USE_METIS
      metis::partition( graph.toViewConst(), m_params.graph.metis, numPart, part.toView() );
#else
      GEOS_THROW( "Attempted to use METIS partition, but GEOSX is not built with METIS support.", InputError );
#endif
      break;
    }
    case Method::scotch:
    {
#ifdef GEOSX_USE_SCOTCH
      scotch::partition( graph.toViewConst(), m_params.graph.scotch, numPart, part.toView() );
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

  partition.move( LvArray::MemorySpace::host, true );
  std::copy( part.begin(), part.end(), partition.begin() );
  return numPart;
}

} // namespace multiscale
} // namespace geos
