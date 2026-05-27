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
 * @file testGraphColoringMPI.cpp
 */

#include "../graphs/GraphTools.hpp"
#include "../graphs/GraphToolsMPI.hpp"

#ifdef GEOS_USE_TRILINOS
#include "../graphs/ZoltanGraphColoring.hpp"
#endif
#include "../graphs/RLFGraphColoringMPI.hpp"

#include "common/MpiWrapper.hpp"

#include <gtest/gtest.h>
#include <iostream>



using namespace geos;
using namespace graph;

class GraphColoringTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  }

  int rank;
};


void runColoringTest( GraphColoringBase & graphColoring, const stdVector< size_t > & xadj, const stdVector< size_t > & adjncy, int expectedNumberOfColors )
{
  auto [localXadj, localAdjncy] = scatterGraphData( xadj, adjncy, MPI_COMM_GEOS );
  int color = graphColoring.colorGraph( localAdjncy );

  EXPECT_TRUE( graphColoring.isColoringValid( localAdjncy, color, MPI_COMM_GEOS ));
  if( expectedNumberOfColors>0 )
  {
    EXPECT_EQ( graphColoring.getNumberOfColors( color, MPI_COMM_GEOS ), expectedNumberOfColors );
  }
}

TEST_F( GraphColoringTest, CartesianDecomposition3D6 )
{
#ifdef GEOS_USE_TRILINOS
  ZoltanGraphColoring zoltanColoring;
#endif
  RLFGraphColoringMPI rlfColoringMPI;

  stdVector< size_t > xadj;
  stdVector< size_t > adjncy;

  if( rank == 0 )
  {
    size_t const nx = 2, ny = 2, nz = 1;
    std::tie( xadj, adjncy ) = generateGraphCartPartitioning3D6( nx, ny, nz );
  }

#ifdef GEOS_USE_TRILINOS
  runColoringTest( zoltanColoring, xadj, adjncy, 2 );
#endif
  runColoringTest( rlfColoringMPI, xadj, adjncy, 2 );
}


TEST_F( GraphColoringTest, CartesianDecomposition3D26 )
{
#ifdef GEOS_USE_TRILINOS
  ZoltanGraphColoring zoltanColoring;
#endif
  RLFGraphColoringMPI rlfColoringMPI;

  stdVector< size_t > xadj;
  stdVector< size_t > adjncy;

  if( rank == 0 )
  {
    size_t const nx = 2, ny = 2, nz = 1;
    std::tie( xadj, adjncy ) = generateGraphCartPartitioning3D26( nx, ny, nz );
  }

#ifdef GEOS_USE_TRILINOS
  runColoringTest( zoltanColoring, xadj, adjncy, 4 );
#endif
  runColoringTest( rlfColoringMPI, xadj, adjncy, 4 );
}


TEST_F( GraphColoringTest, RandomGraphs )
{
#ifdef GEOS_USE_TRILINOS
  ZoltanGraphColoring zoltanColoring;
#endif
  RLFGraphColoringMPI rlfColoringMPI;

  size_t const iterations = 10;
  for( size_t i = 0; i < iterations; ++i )
  {
    stdVector< size_t > xadj;
    stdVector< size_t > adjncy;

    if( rank == 0 )
    {
      size_t num_nodes = 4;
      size_t num_edges = rand() % (num_nodes * 3 + 1) + num_nodes;
      std::tie( xadj, adjncy ) = generateGraphRandom( num_nodes, num_edges );
    }

#ifdef GEOS_USE_TRILINOS
    runColoringTest( zoltanColoring, xadj, adjncy, 0 );
#endif
    runColoringTest( rlfColoringMPI, xadj, adjncy, 0 );
  }
}


TEST_F( GraphColoringTest, CountPositiveDistinctColors )
{
  stdVector< int > colors = {1, -1, 3, 2, 1, 4, 5, 3};
  EXPECT_EQ( GraphColoringBase::getNumberOfColors( colors, MPI_COMM_GEOS ), 6 );
}


int main( int argc, char * *argv )
{
  // Initialize MPI
  MpiWrapper::init( &argc, &argv );
  MPI_COMM_GEOS = MpiWrapper::commDup( MPI_COMM_WORLD );

  // Initialize Google Test
  ::testing::InitGoogleTest( &argc, argv );

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  // Suppress output from non-root ranks
  if( rank != 0 )
  {
    ::testing::TestEventListeners & listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release( listeners.default_result_printer());
  }

  // Run all tests
  int result = RUN_ALL_TESTS();

  // Finalize MPI
  MpiWrapper::finalize();

  return result;
}
