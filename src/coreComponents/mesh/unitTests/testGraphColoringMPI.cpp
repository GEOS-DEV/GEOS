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
#include "../graphs/GraphToolsMPI.hpp"

#include "../graphs/GraphColoringBase.hpp"
#include "../graphs/ZoltanGraphColoring.hpp"
#include "../graphs/RLFGraphColoringMPI.hpp"

#include <gtest/gtest.h>
#include <iostream>
#include <mpi.h>

using namespace geos;
using namespace graph;

class GraphColoringTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  }

  int rank;
};

void runColoringTest( GraphColoringBase & graphColoring, const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy, int expectedNumberOfColors )
{
  auto [localXadj, localAdjncy] = scatterGraphData( xadj, adjncy, MPI_COMM_WORLD );
  int color = graphColoring.colorGraph( localAdjncy );

  EXPECT_TRUE( graphColoring.isColoringValid( localAdjncy, color, MPI_COMM_WORLD ));
  if( expectedNumberOfColors>0 )
  {
    EXPECT_EQ( graphColoring.getNumberOfColors( color, MPI_COMM_WORLD ), expectedNumberOfColors );
  }
}

TEST_F( GraphColoringTest, CartesianDecomposition3D6 )
{
  ZoltanGraphColoring zoltanColoring;
  RLFGraphColoringMPI rlfColoringMPI;

  std::vector< camp::idx_t > xadj;
  std::vector< camp::idx_t > adjncy;

  if( rank == 0 )
  {
    idx_t const nx = 4, ny = 2, nz = 1;
    std::tie( xadj, adjncy ) = generateGraphCartPartitionning3D6( nx, ny, nz );
  }

  runColoringTest( zoltanColoring, xadj, adjncy, 2 );
  runColoringTest( rlfColoringMPI, xadj, adjncy, 2 );
}

TEST_F( GraphColoringTest, CartesianDecomposition3D26 )
{
  ZoltanGraphColoring zoltanColoring;
  RLFGraphColoringMPI rlfColoringMPI;

  std::vector< camp::idx_t > xadj;
  std::vector< camp::idx_t > adjncy;

  if( rank == 0 )
  {
    idx_t const nx = 2, ny = 2, nz = 2;
    std::tie( xadj, adjncy ) = generateGraphCartPartitionning3D26( nx, ny, nz );
  }

  runColoringTest( zoltanColoring, xadj, adjncy, 8 );
  runColoringTest( rlfColoringMPI, xadj, adjncy, 8 );
}

TEST_F( GraphColoringTest, RandomGraphs )
{
  ZoltanGraphColoring zoltanColoring;
  RLFGraphColoringMPI rlfColoringMPI;

  size_t const iterations = 10;
  for( size_t i = 0; i < iterations; ++i )
  {
    std::vector< camp::idx_t > xadj;
    std::vector< camp::idx_t > adjncy;

    if( rank == 0 )
    {
      size_t num_nodes = 8;
      size_t num_edges = rand() % (num_nodes * 3 + 1) + num_nodes;
      std::tie( xadj, adjncy ) = generateGraphRandom( num_nodes, num_edges );
    }

    runColoringTest( zoltanColoring, xadj, adjncy, 0 );
    runColoringTest( rlfColoringMPI, xadj, adjncy, 0 );
  }
}

TEST_F( GraphColoringTest, CountPositiveDistinctColors )
{
  std::vector< int > colors = {1, -1, 3, 2, 1, 4, 5, 3};
  EXPECT_EQ( GraphColoringBase::getNumberOfColors( colors, MPI_COMM_WORLD ), 6 );
}

int main( int argc, char * *argv )
{
  // Initialize MPI
  MPI_Init( &argc, &argv );

  // Initialize Google Test
  ::testing::InitGoogleTest( &argc, argv );

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Suppress output from non-root ranks
  if( rank != 0 )
  {
    ::testing::TestEventListeners & listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release( listeners.default_result_printer());
  }

  // Run all tests
  int result = RUN_ALL_TESTS();

  // Finalize MPI
  MPI_Finalize();

  return result;
}
