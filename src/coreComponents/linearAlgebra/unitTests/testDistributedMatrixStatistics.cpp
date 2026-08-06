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

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "linearAlgebra/utilities/DistributedMatrixStatistics.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

using namespace geos;

TEST( DistributedMatrixStatistics, TwoRankPattern )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) != 2 )
  {
    GTEST_SKIP() << "This pattern exercises a two-rank halo exchange";
  }

  CRSMatrix< real64, globalIndex > matrix( 2, 4, 3 );
  real64 const values[3] = { 1.0, 1.0, 1.0 };
  if( rank == 0 )
  {
    globalIndex const row0Columns[2] = { 0, 2 };
    globalIndex const row1Columns[3] = { 0, 1, 3 };
    matrix.insertNonZeros( 0, row0Columns, values, 2 );
    matrix.insertNonZeros( 1, row1Columns, values, 3 );
  }
  else
  {
    globalIndex const row0Columns[2] = { 0, 2 };
    globalIndex const row1Columns[3] = { 1, 2, 3 };
    matrix.insertNonZeros( 0, row0Columns, values, 2 );
    matrix.insertNonZeros( 1, row1Columns, values, 3 );
  }

  DistributedMatrixStatistics const statistics =
    computeDistributedMatrixStatistics( matrix.toViewConst(), 2 * rank, MPI_COMM_GEOS );

  ASSERT_EQ( statistics.ranks.size(), 2 );
  for( int metricRank = 0; metricRank < 2; ++metricRank )
  {
    MatrixRankStatistics const & rankStatistics = statistics.ranks[metricRank];
    EXPECT_EQ( rankStatistics.firstRow, 2 * metricRank );
    EXPECT_EQ( rankStatistics.numRows, 2 );
    EXPECT_EQ( rankStatistics.numNonzeros, 5 );
    EXPECT_EQ( rankStatistics.numOffRankNonzeros, 2 );
    EXPECT_EQ( rankStatistics.numHaloReceiveDofs, 2 );
    EXPECT_EQ( rankStatistics.numHaloSendDofs, 2 );
    EXPECT_EQ( rankStatistics.numReceiveNeighbors, 1 );
    EXPECT_EQ( rankStatistics.numSendNeighbors, 1 );
    EXPECT_EQ( rankStatistics.maxRowNonzeros, 3 );
    EXPECT_EQ( rankStatistics.numEmptyRows, 0 );
  }
}

#ifdef GEOS_USE_HYPRE
TEST( DistributedMatrixStatistics, StreamedMatrixMarketOutput )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) != 2 )
  {
    GTEST_SKIP() << "This output test exercises two distributed row blocks";
  }

  CRSMatrix< real64, globalIndex > localMatrix( 2, 4, 3 );
  real64 const values[3] = { 1.0, 2.0, 3.0 };
  if( rank == 0 )
  {
    globalIndex const row0Columns[2] = { 0, 2 };
    globalIndex const row1Columns[3] = { 0, 1, 3 };
    localMatrix.insertNonZeros( 0, row0Columns, values, 2 );
    localMatrix.insertNonZeros( 1, row1Columns, values, 3 );
  }
  else
  {
    globalIndex const row0Columns[2] = { 0, 2 };
    globalIndex const row1Columns[3] = { 1, 2, 3 };
    localMatrix.insertNonZeros( 0, row0Columns, values, 2 );
    localMatrix.insertNonZeros( 1, row1Columns, values, 3 );
  }

  string const filename = "testDistributedMatrixStatistics_output.mtx";
  if( rank == 0 )
  {
    std::remove( filename.c_str() );
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );

  HypreMatrix matrix;
  matrix.create( localMatrix.toViewConst(), 2, MPI_COMM_GEOS );
  matrix.write( filename, LAIOutputFormat::MATRIX_MARKET );

  if( rank == 0 )
  {
    std::ifstream input( filename );
    ASSERT_TRUE( input );

    string line;
    ASSERT_TRUE( static_cast< bool >( std::getline( input, line ) ) );
    EXPECT_EQ( line, "%%MatrixMarket matrix coordinate real general" );
    ASSERT_TRUE( static_cast< bool >( std::getline( input, line ) ) );

    std::istringstream dimensions( line );
    globalIndex rows = 0;
    globalIndex columns = 0;
    globalIndex nonzeros = 0;
    dimensions >> rows >> columns >> nonzeros;
    EXPECT_EQ( rows, 4 );
    EXPECT_EQ( columns, 4 );
    EXPECT_EQ( nonzeros, 10 );

    std::map< std::pair< globalIndex, globalIndex >, real64 > entries;
    globalIndex row = 0;
    globalIndex column = 0;
    real64 value = 0.0;
    while( input >> row >> column >> value )
    {
      auto const [it, inserted] = entries.emplace( std::make_pair( row - 1, column - 1 ), value );
      GEOS_UNUSED_VAR( it );
      EXPECT_TRUE( inserted );
    }
    EXPECT_EQ( entries.size(), nonzeros );

    std::map< std::pair< globalIndex, globalIndex >, real64 > const expected =
    {
      { { 0, 0 }, 1.0 }, { { 0, 2 }, 2.0 },
      { { 1, 0 }, 1.0 }, { { 1, 1 }, 2.0 }, { { 1, 3 }, 3.0 },
      { { 2, 0 }, 1.0 }, { { 2, 2 }, 2.0 },
      { { 3, 1 }, 1.0 }, { { 3, 2 }, 2.0 }, { { 3, 3 }, 3.0 }
    };
    EXPECT_EQ( entries, expected );
    std::remove( filename.c_str() );
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );
}
#endif

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
