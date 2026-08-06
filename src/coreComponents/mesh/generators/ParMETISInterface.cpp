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

/**
 * @file ParMETISInterface.cpp
 */

#include "ParMETISInterface.hpp"

#include "common/GEOS_RAJA_Interface.hpp"

#include <parmetis.h>

#include <numeric>

#define GEOS_PARMETIS_CHECK( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_NE_MSG( ierr, METIS_OK, "Error in call to:\n" #call ); \
  } while( false )

namespace geos
{
namespace parmetis
{

static_assert( std::is_same< idx_t, pmet_idx_t >::value, "Non-matching index types. ParMETIS must be built with 64-bit indices." );

ArrayOfArrays< idx_t, idx_t >
meshToDual( ArrayOfArraysView< idx_t const, idx_t > const & elemToNodes,
            arrayView1d< idx_t const > const & elemDist,
            MPI_Comm comm,
            int const minCommonNodes )
{
  idx_t const numElems = elemToNodes.size();

  // `parmetis` awaits the arrays to be allocated as two continuous arrays: one for values, the other for offsets.
  // Our `ArrayOfArrays` allows to reserve some extra space for further element insertion,
  // but this is not compatible with what `parmetis` requires.
  GEOS_ASSERT_EQ_MSG( std::accumulate( elemToNodes.getSizes(), elemToNodes.getSizes() + numElems, 0 ),
                      elemToNodes.valueCapacity(),
                      "Internal error. The element to nodes mapping must be strictly allocated for compatibility with a third party library." );

  idx_t numflag = 0;
  idx_t ncommonnodes = minCommonNodes;
  idx_t * xadj;
  idx_t * adjncy;

  // ParMETIS does not handle ranks with zero local elements (it crashes on NULL array pointers).
  // We provide dummy non-NULL arrays for empty ranks so the collective call can proceed.
  idx_t dummyOffset = 0;
  idx_t dummyValue = 0;
  idx_t * eptr = numElems > 0 ? const_cast< idx_t * >( elemToNodes.getOffsets() ) : &dummyOffset;
  idx_t * eind = numElems > 0 ? const_cast< idx_t * >( elemToNodes.getValues() ) : &dummyValue;

  // Technical UB if ParMETIS writes into these arrays; in practice we discard them right after
  GEOS_PARMETIS_CHECK( ParMETIS_V3_Mesh2Dual( const_cast< idx_t * >( elemDist.data() ),
                                              eptr,
                                              eind,
                                              &numflag, &ncommonnodes, &xadj, &adjncy, &comm ) );

  ArrayOfArrays< idx_t, idx_t > graph;

  if( numElems > 0 )
  {
    graph.resizeFromOffsets( numElems, xadj );

    // There is no way to direct-copy values into ArrayOfArrays without UB (casting away const)
    forAll< parallelHostPolicy >( numElems, [xadj, adjncy, graph = graph.toView()]( localIndex const k )
    {
      graph.appendToArray( k, adjncy + xadj[k], adjncy + xadj[k+1] );
    } );

    METIS_Free( xadj );
    METIS_Free( adjncy );
  }

  return graph;
}

array1d< idx_t >
partition( ArrayOfArraysView< idx_t const, idx_t > const & graph,
           arrayView1d< idx_t const > const & vertDist,
           idx_t const numParts,
           MPI_Comm comm,
           int const numRefinements )
{
  GEOS_ERROR_IF( numParts <= 0, "Number of partitions must be strictly positive" );

  array1d< idx_t > part( graph.size() ); // all 0 by default
  if( numParts == 1 )
  {
    return part;
  }

  // Compute tpwgts parameters (target partition weights)
  array1d< real_t > tpwgts( numParts );
  tpwgts.setValues< serialPolicy >( 1.0f / static_cast< real_t >( numParts ) );

  // Set other ParMETIS parameters
  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t npart = numParts;
  idx_t options[4] = { 1, 0, 2022, PARMETIS_PSR_UNCOUPLED };
  idx_t edgecut = 0;
  real_t ubvec = 1.05;

  // ParMETIS does not handle ranks with zero local vertices (it crashes on NULL adjncy pointer).
  // We provide dummy non-NULL pointers for empty ranks so the collective call can proceed.
  idx_t dummyAdj = 0;
  idx_t dummyPart = 0;
  idx_t * xadj = const_cast< idx_t * >( graph.getOffsets() );
  idx_t * adjncy = graph.size() > 0 && graph.getValues() != nullptr
                   ? const_cast< idx_t * >( graph.getValues() )
                   : &dummyAdj;
  idx_t * partData = graph.size() > 0 ? part.data() : &dummyPart;

  // Technical UB if ParMETIS writes into these arrays; in practice we discard them right after
  GEOS_PARMETIS_CHECK( ParMETIS_V3_PartKway( const_cast< idx_t * >( vertDist.data() ),
                                             xadj,
                                             adjncy,
                                             nullptr, nullptr, &wgtflag,
                                             &numflag, &ncon, &npart, tpwgts.data(),
                                             &ubvec, options, &edgecut, partData, &comm ) );

  for( int iter = 0; iter < numRefinements; ++iter )
  {
    GEOS_PARMETIS_CHECK( ParMETIS_V3_RefineKway( const_cast< idx_t * >( vertDist.data() ),
                                                 xadj,
                                                 adjncy,
                                                 nullptr, nullptr, &wgtflag,
                                                 &numflag, &ncon, &npart, tpwgts.data(),
                                                 &ubvec, options, &edgecut, partData, &comm ) );
  }

  return part;
}

array1d< idx_t >
partitionWeighted( ArrayOfArraysView< idx_t const, idx_t > const & graph,
                   arrayView1d< idx_t const > const & vertexWeights,
                   arrayView1d< idx_t const > const & vertDist,
                   idx_t const numParts,
                   MPI_Comm comm,
                   int const numRefinements )
{
  GEOS_ERROR_IF( numParts <= 0, "Number of partitions must be strictly positive" );

  array1d< idx_t > part( graph.size() );
  if( numParts == 1 )
  {
    return part;
  }

  array1d< real_t > tpwgts( numParts );
  tpwgts.setValues< serialPolicy >( 1.0f / static_cast< real_t >( numParts ) );

  idx_t wgtflag = 2;  // vertex weights only
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t npart = numParts;

  // Options: [use_defaults, log_level, seed, coupling]
  idx_t options[4] = { 1, 0, 2022, PARMETIS_PSR_UNCOUPLED };

  idx_t edgecut = 0;
  real_t ubvec = 1.05;

  GEOS_PARMETIS_CHECK( ParMETIS_V3_PartKway(
                         const_cast< idx_t * >( vertDist.data() ),
                         const_cast< idx_t * >( graph.getOffsets() ),
                         const_cast< idx_t * >( graph.getValues() ),
                         const_cast< idx_t * >( vertexWeights.data() ),
                         nullptr, // edge weights
                         &wgtflag,
                         &numflag, &ncon, &npart, tpwgts.data(),
                         &ubvec, options, &edgecut, part.data(), &comm ) );

  for( int iter = 0; iter < numRefinements; ++iter )
  {
    GEOS_PARMETIS_CHECK( ParMETIS_V3_RefineKway(
                           const_cast< idx_t * >( vertDist.data() ),
                           const_cast< idx_t * >( graph.getOffsets() ),
                           const_cast< idx_t * >( graph.getValues() ),
                           const_cast< idx_t * >( vertexWeights.data() ),
                           nullptr,
                           &wgtflag,
                           &numflag, &ncon, &npart, tpwgts.data(),
                           &ubvec, options, &edgecut, part.data(), &comm ) );
  }

  return part;
}

array1d< idx_t >
partitionWeighted( ArrayOfArraysView< idx_t const, idx_t > const & graph,
                   arrayView1d< idx_t const > const & edgeWeights,
                   arrayView2d< idx_t const > const & vertexWeights,
                   arrayView1d< idx_t const > const & vertDist,
                   idx_t const numParts,
                   arrayView1d< real64 const > const & imbalance,
                   MPI_Comm comm,
                   int const numRefinements,
                   idx_t const seed )
{
  GEOS_ERROR_IF( numParts <= 0, "Number of partitions must be strictly positive" );

  idx_t const numVertices = graph.size();
  idx_t const numConstraints = vertexWeights.size( 1 );
  idx_t const numAdjacencies = numVertices > 0 ? graph.getOffsets()[numVertices] : 0;
  GEOS_ERROR_IF_NE_MSG( vertexWeights.size( 0 ), numVertices,
                        "ParMETIS vertex-weight row count must match the local graph size" );
  GEOS_ERROR_IF( numConstraints <= 0, "ParMETIS requires at least one balance constraint" );
  GEOS_ERROR_IF_NE_MSG( imbalance.size(), numConstraints,
                        "ParMETIS requires one imbalance tolerance per constraint" );
  GEOS_ERROR_IF_NE_MSG( edgeWeights.size(), numAdjacencies,
                        "ParMETIS edge weights must match the CSR adjacency array" );
  for( idx_t i = 0; i < numAdjacencies; ++i )
  {
    GEOS_ERROR_IF( edgeWeights[i] <= 0, "ParMETIS edge weights must be positive" );
  }
  for( idx_t i = 0; i < numVertices; ++i )
  {
    for( idx_t k = 0; k < numConstraints; ++k )
    {
      GEOS_ERROR_IF( vertexWeights( i, k ) < 0, "ParMETIS vertex weights must be nonnegative" );
    }
  }

  array1d< idx_t > part( numVertices );
  if( numParts == 1 )
  {
    return part;
  }

  array1d< real_t > targetWeights( numParts * numConstraints );
  targetWeights.setValues< serialPolicy >( 1.0f / static_cast< real_t >( numParts ) );

  array1d< real_t > imbalanceVector( numConstraints );
  for( idx_t k = 0; k < numConstraints; ++k )
  {
    GEOS_ERROR_IF( imbalance[k] < 0.0,
                   GEOS_FMT( "ParMETIS imbalance tolerance {} is negative", imbalance[k] ) );
    imbalanceVector[k] = static_cast< real_t >( 1.0 + imbalance[k] );
  }

  idx_t wgtflag = 3;
  idx_t numflag = 0;
  idx_t ncon = numConstraints;
  idx_t npart = numParts;
  idx_t options[4] = { 1, 0, seed, PARMETIS_PSR_UNCOUPLED };
  idx_t edgecut = 0;

  idx_t dummyOffset = 0;
  idx_t dummyAdjacency = 0;
  idx_t dummyWeight = 0;
  idx_t dummyPart = 0;
  idx_t * xadj = numVertices > 0 ? const_cast< idx_t * >( graph.getOffsets() ) : &dummyOffset;
  idx_t * adjncy = numAdjacencies > 0 ? const_cast< idx_t * >( graph.getValues() ) : &dummyAdjacency;
  idx_t * adjacencyWeights = numAdjacencies > 0 ? const_cast< idx_t * >( edgeWeights.data() ) : &dummyWeight;
  idx_t * vertexWeightData = numVertices > 0 ? const_cast< idx_t * >( vertexWeights.data() ) : &dummyWeight;
  idx_t * partData = numVertices > 0 ? part.data() : &dummyPart;

  GEOS_PARMETIS_CHECK( ParMETIS_V3_PartKway(
                         const_cast< idx_t * >( vertDist.data() ),
                         xadj,
                         adjncy,
                         vertexWeightData,
                         adjacencyWeights,
                         &wgtflag,
                         &numflag, &ncon, &npart, targetWeights.data(),
                         imbalanceVector.data(), options, &edgecut, partData, &comm ) );

  for( int iter = 0; iter < numRefinements; ++iter )
  {
    GEOS_PARMETIS_CHECK( ParMETIS_V3_RefineKway(
                           const_cast< idx_t * >( vertDist.data() ),
                           xadj,
                           adjncy,
                           vertexWeightData,
                           adjacencyWeights,
                           &wgtflag,
                           &numflag, &ncon, &npart, targetWeights.data(),
                           imbalanceVector.data(), options, &edgecut, partData, &comm ) );
  }

  return part;
}

} // namespace parmetis
} // namespace geos
