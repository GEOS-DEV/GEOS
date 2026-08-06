/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testMETISInterface.cpp
 */

#include "mesh/generators/METISInterface.hpp"

#include <gtest/gtest.h>

namespace geos
{
namespace
{

TEST( METISInterface, WeightedMultiConstraintPartitionIsBalancedAndDeterministic )
{
  constexpr pmet_idx_t numVertices = 8;
  constexpr pmet_idx_t numConstraints = 3;

  // A path with two additional cross links exercises nonuniform edge weights
  // while retaining an unambiguous, balanced two-way partition.
  pmet_idx_t const capacities[numVertices] = { 1, 3, 2, 3, 3, 2, 3, 1 };
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > graph;
  graph.resizeFromCapacities< serialPolicy >( numVertices, capacities );

  auto addUndirectedEdge = [&graph]( pmet_idx_t const first, pmet_idx_t const second )
  {
    graph.emplaceBack( first, second );
    graph.emplaceBack( second, first );
  };
  for( pmet_idx_t vertex = 0; vertex + 1 < numVertices; ++vertex )
  {
    addUndirectedEdge( vertex, vertex + 1 );
  }
  addUndirectedEdge( 1, 3 );
  addUndirectedEdge( 4, 6 );

  array1d< pmet_idx_t > edgeWeights( graph.toViewConst().getOffsets()[numVertices] );
  edgeWeights.setValues< serialPolicy >( 1 );
  for( pmet_idx_t edge = 0; edge < edgeWeights.size(); ++edge )
  {
    // Use nonuniform positive weights to verify that the weighted API path is active.
    edgeWeights[edge] += edge % 3;
  }

  array2d< pmet_idx_t > vertexWeights( numVertices, numConstraints );
  for( pmet_idx_t vertex = 0; vertex < numVertices; ++vertex )
  {
    vertexWeights( vertex, 0 ) = 1;
    vertexWeights( vertex, 1 ) = 1 + vertex % 2;
    vertexWeights( vertex, 2 ) = 2 - vertex % 2;
  }

  array1d< real64 > imbalance( numConstraints );
  imbalance.setValues< serialPolicy >( 0.05 );

  array1d< pmet_idx_t > const first = metis::partitionWeighted( graph.toViewConst(),
                                                                edgeWeights.toViewConst(),
                                                                vertexWeights.toViewConst(),
                                                                2,
                                                                imbalance.toViewConst(),
                                                                739 );
  array1d< pmet_idx_t > const second = metis::partitionWeighted( graph.toViewConst(),
                                                                 edgeWeights.toViewConst(),
                                                                 vertexWeights.toViewConst(),
                                                                 2,
                                                                 imbalance.toViewConst(),
                                                                 739 );

  pmet_idx_t counts[2] = { 0, 0 };
  pmet_idx_t loads[2][numConstraints] = {};
  for( pmet_idx_t vertex = 0; vertex < numVertices; ++vertex )
  {
    ASSERT_GE( first[vertex], 0 );
    ASSERT_LT( first[vertex], 2 );
    EXPECT_EQ( first[vertex], second[vertex] );
    ++counts[first[vertex]];
    for( pmet_idx_t constraint = 0; constraint < numConstraints; ++constraint )
    {
      loads[first[vertex]][constraint] += vertexWeights( vertex, constraint );
    }
  }

  EXPECT_GT( counts[0], 0 );
  EXPECT_GT( counts[1], 0 );
  for( pmet_idx_t constraint = 0; constraint < numConstraints; ++constraint )
  {
    EXPECT_EQ( loads[0][constraint], loads[1][constraint] );
  }
}

} // namespace
} // namespace geos
