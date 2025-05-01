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
 * @file testColoring.cpp
 */

#include "../graphs/GraphTools.hpp"
#include "../graphs/GraphColoringBase.hpp"
#include "../graphs/RLFGraphColoring.hpp"

#include <gtest/gtest.h>


using namespace geos;
using namespace graph;
TEST( GraphColoringTest, CountPositiveDistinctColors )
{
  std::vector< int > colors = {1, -1, 3, 2, 1, 4, 5, 3, -1, 0};
  EXPECT_EQ( GraphColoringBase::getNumberOfColors( colors ), 6 );
}


TEST( GraphColoringTest, CartesianDecomposition3D6 )
{
  idx_t const nx = 3, ny = 4, nz = 3;
  auto [xadj, adjncy] = generateGraphCartPartitionning3D6( nx, ny, nz );
  geos::graph::RLFGraphColoring graphColoring;
  std::vector< int > colors = graphColoring.colorGraph( xadj, adjncy );

  EXPECT_TRUE( graphColoring.isColoringValid( xadj, adjncy, colors ));
  EXPECT_EQ( graphColoring.getNumberOfColors( colors ), 2 );
}


TEST( GraphColoringTest, CartesianDecomposition3D26 )
{
  idx_t const nx = 3, ny = 4, nz = 3;
  auto [xadj, adjncy] = generateGraphCartPartitionning3D26( nx, ny, nz );
  geos::graph::RLFGraphColoring graphColoring;
  std::vector< int > colors = graphColoring.colorGraph( xadj, adjncy );
  EXPECT_TRUE( graphColoring.isColoringValid( xadj, adjncy, colors ));
  EXPECT_EQ( graphColoring.getNumberOfColors( colors ), 8 );
}


TEST( GraphColoringTest, RandomGraphs )
{
  size_t const iterations = 10;
  for( size_t i = 0; i < iterations; ++i )
  {
    size_t num_nodes = rand() % 100 + 5;     // between 5 and 104
    size_t num_edges = rand() % (num_nodes * 6 + 1) + num_nodes;     // between num_nodes and num_nodes * 6
    auto [xadj, adjncy] = generateGraphRandom( num_nodes, num_edges );
    geos::graph::RLFGraphColoring graphColoring;
    std::vector< int > colors = graphColoring.colorGraph( xadj, adjncy );
    EXPECT_TRUE( graphColoring.isColoringValid( xadj, adjncy, colors ));
  }
}


TEST( GraphColoringTest, InvalidColoring )
{
  // Create a simple graph with 4 nodes and 3 edges
  std::vector< idx_t > xadj = {0, 1, 2, 3, 3};
  std::vector< idx_t > adjncy = {1, 0, 2, 1};

  // Intentionally create an invalid coloring where two adjacent nodes have the same color
  std::vector< int > colors = {0, 0, 1, 1};

  // Check if the coloring is valid (should fail)
  EXPECT_FALSE( GraphColoringBase::isColoringValid( xadj, adjncy, colors ));
}
