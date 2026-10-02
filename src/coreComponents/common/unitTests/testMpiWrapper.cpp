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

#include "common/MpiWrapper.hpp"
#include "common/MpiChunkedCommunication.hpp"
#include "common/initializeEnvironment.hpp"

#include <gtest/gtest.h>

#include <set>

using namespace geos;


class MpiTestScope
{
public:

  MpiTestScope( int argc, char * argv[] )
  {
    ::testing::InitGoogleTest( &argc, argv );
    geos::setupEnvironment( argc, argv );
  }

  ~MpiTestScope()
  {
    geos::cleanupEnvironment();
  }

};

template< typename FIRST, typename SECOND >
struct PairTestCase
{
  string_view testName;
  stdVector< stdVector< MpiWrapper::PairType< FIRST, SECOND > > > pairs; // Pair set for each rank
  MpiWrapper::PairType< FIRST, SECOND > expectedMin; // Expected result for Min reduction
  MpiWrapper::PairType< FIRST, SECOND > expectedMax; // Expected result for Max reduction
};

template< typename FIRST, typename SECOND >
string reductionTestFailureMsg( PairTestCase< FIRST, SECOND > const & testCase,
                                MpiWrapper::PairType< FIRST, SECOND > pair,
                                int rankId,
                                string_view opName )
{
  return GEOS_FMT( "Test case '{}', rank {}:\n  Error in {} reduction, incorrect first value in pair ({}, {}).",
                   testCase.testName, rankId, opName, pair.first, pair.second );
}

/**
 * @name Testing the reduction of single pair per ranks
 */
namespace pairTestCases
{

template< typename FIRST, typename SECOND >
void runAllTestCases( int const rankId )
{
  stdVector< PairTestCase< FIRST, SECOND > > const testCases = {
    {
      "Basic case with distinct values",
      { { { 9, 23 }, { 18, 886 }, { 19, 26 }, { 8, 549 } } }, // pairs
      { 8, 549 }, // expectedMin
      { 19, 26 } // expectedMax
    },
    {
      "Identical second values",
      { { { 23, 12 }, { 886, 12 }, { -26, 12 }, { 549, 12 } } }, // pairs
      { -26, 12 }, // expectedMin
      { 886, 12 } // expectedMax
    },
    {
      "Identical first values",
      { { { 12, 23 }, { 12, 886 }, { 12, -26 }, { 12, 549 } } }, // pairs
      { 12, -26 }, // expectedMin
      { 12, 886 } // expectedMax
    },
  };
  for( auto const & testCase : testCases )
  {
    MpiWrapper::PairType< FIRST, SECOND > min = MpiWrapper::min( testCase.pairs[0][rankId] );
    EXPECT_EQ( min.first, testCase.expectedMin.first ) << reductionTestFailureMsg( testCase, min, rankId, "Min" );
    EXPECT_EQ( min.second, testCase.expectedMin.second ) << reductionTestFailureMsg( testCase, min, rankId, "Min" );
    MpiWrapper::PairType< FIRST, SECOND > max = MpiWrapper::max( testCase.pairs[0][rankId] );
    EXPECT_EQ( max.first, testCase.expectedMax.first ) << reductionTestFailureMsg( testCase, max, rankId, "Max" );
    EXPECT_EQ( max.second, testCase.expectedMax.second ) << reductionTestFailureMsg( testCase, max, rankId, "Max" );
  }
}

TEST( MpiWrapperTesting, MpiPairReductionTest )
{
  int const rankId = MpiWrapper::commRank();
  int const nbRanks = MpiWrapper::commSize();
  if( nbRanks==4 )
  { // here, we only want to test the MIN/MAXLOC behaviour
    // basic MPI types
    runAllTestCases< int32_t, int32_t >( rankId );
    runAllTestCases< double, int32_t >( rankId );

    // custom MPI types
    runAllTestCases< int64_t, int64_t >( rankId );
    runAllTestCases< double, int64_t >( rankId );
    runAllTestCases< double, double >( rankId );
  }
}

} /* namespace pairTestCases */

/**
 * @name Testing the reduction of sets over ranks
 */
namespace pairSetsTestCases
{

template< typename FIRST, typename SECOND >
void runTestCase( PairTestCase< FIRST, SECOND > const & testCase, int rankId, int nbRanks )
{
  // Create local pairs for 1 or n ranks
  stdVector< MpiWrapper::PairType< FIRST, SECOND > > localPairs;
  if( nbRanks == 1 )
  {
    // Special case when we test on only one rank : all rank data are combined
    for( auto const & rankPairs : testCase.pairs )
    {
      localPairs.insert( localPairs.end(), rankPairs.begin(), rankPairs.end());
    }
  }
  else
  {
    // MPI testing : each ranks get its own data
    ASSERT_EQ( testCase.pairs.size(), nbRanks )
      << "Test case " << testCase.testName << ": The number of ranks does not match the test case data.";
    localPairs = testCase.pairs[rankId];
  }

  MpiWrapper::PairType< FIRST, SECOND > min = MpiWrapper::min< FIRST, SECOND >( localPairs );
  EXPECT_EQ( min.first, testCase.expectedMin.first ) << reductionTestFailureMsg( testCase, min, rankId, "Min" );
  EXPECT_EQ( min.second, testCase.expectedMin.second ) << reductionTestFailureMsg( testCase, min, rankId, "Min" );

  MpiWrapper::PairType< FIRST, SECOND > max = MpiWrapper::max< FIRST, SECOND >( localPairs );
  EXPECT_EQ( max.first, testCase.expectedMax.first ) << reductionTestFailureMsg( testCase, max, rankId, "Max" );
  EXPECT_EQ( max.second, testCase.expectedMax.second ) << reductionTestFailureMsg( testCase, max, rankId, "Max" );
}

template< typename FIRST, typename SECOND >
void runAllTestCases( int rankId, int nbRanks )
{
  // the set of cases is designed for 4 ranks
  stdVector< PairTestCase< FIRST, SECOND > > const testCases = {
    {
      "Basic case with distinct values",
      { { { 55, 549 }, { 34, 886 } }, { { 68, 23 }, { 10, 6 } }, { { 59, -781 }, { 64, 823 } }, { { 15, 543 }, { 22, 363 } } }, // pairs
      { 10, 6 }, // expectedMin
      { 68, 23 } // expectedMax
    },
    {
      "Distinct values but the first two ranks has no data",
      { { }, { }, { { 59, -781 }, { 64, 823 } }, { { 15, 543 }, { 22, 363 } } }, // pairs
      { 15, 543 }, // expectedMin
      { 64, 823 } // expectedMax
    },
    {
      "Distinct values, but with different data count and the last two ranks has no data",
      { { { 55, 549 }, { 34, 886 }, { 42, 23 }, { 10, 6 }, { 59, -781 }, { 64, 823 } }, { { 68, 23 }, { 22, 363 } }, { }, { } }, // pairs
      { 10, 6 }, // expectedMin
      { 68, 23 } // expectedMax
    },
    {
      "No rank have any data",
      { { }, { }, { }, { } }, // pairs
      { std::numeric_limits< FIRST >::max(), std::numeric_limits< SECOND >::max() }, // expectedMin
      { std::numeric_limits< FIRST >::lowest(), std::numeric_limits< SECOND >::lowest() }, // expectedMax
    },
    {
      "All pairs have the same value but different indices",
      { { { 12, 549 }, { 12, 886 } }, { { 12, 23 }, { 12, -6 } }, { { 12, 781 }, { 12, 823 } }, { { 12, 543 }, { 12, 363 } } }, // pairs
      { 12, -6 }, // expectedMin
      { 12, 886 } // expectedMax
    },
    {
      "All pairs have the same value but different indices, variation",
      { { { 12, -6 }, { 12, 886 } }, { { 12, 23 }, { 12, 549 } }, { { 12, 781 }, { 12, 823 } }, { { 12, 543 }, { 12, 363 } } }, // pairs
      { 12, -6 }, // expectedMin
      { 12, 886 } // expectedMax
    },
    {
      "All pairs have the same value except one",
      { { { 12, 549 }, { 12, 886 } }, { { 3, 23 }, { 12, 6 } }, { { 12, 781 }, { 12, 823 } }, { { 12, 543 }, { 12, 363 } } }, // pairs
      { 3, 23 }, // expectedMin
      { 12, 886 } // expectedMax
    },
    {
      "A unique value for one rank",
      { { { 50, 50 }, { 50, 50 } }, { { 75, -23 }, { 50, 50 } }, { { 50, 50 }, { 50, 50 } }, { { 50, 50 }, { 50, 50 } } }, // pairs
      { 50, 50 }, // expectedMin
      { 75, -23 } // expectedMax
    },
    {
      "All pairs are identical",
      { { { 123, 456 }, { 123, 456 } }, { { 123, 456 }, { 123, 456 } }, { { 123, 456 }, { 123, 456 } }, { { 123, 456 }, { 123, 456 } } }, // pairs
      { 123, 456 }, // expectedMin
      { 123, 456 } // expectedMax
    },
  };
  for( auto const & testCase : testCases )
  {
    runTestCase< FIRST, SECOND >( testCase, rankId, nbRanks );
  }
}

TEST( MpiWrapperTesting, MpiPairSetReductionTest )
{
  int const rankId = MpiWrapper::commRank();
  int const nbRanks = MpiWrapper::commSize();

  // basic MPI types
  runAllTestCases< int32_t, int32_t >( rankId, nbRanks );
  runAllTestCases< double, int32_t >( rankId, nbRanks );

  // custom MPI types
  runAllTestCases< int64_t, int64_t >( rankId, nbRanks );
  runAllTestCases< double, int64_t >( rankId, nbRanks );
  runAllTestCases< double, double >( rankId, nbRanks );
}

} /* namespace pairSetsTestCases */

TEST( MpiWrapperTesting, GatherCopiesTypedValuesAndHonorsDisplacements )
{
  int const rank = MpiWrapper::commRank();
  int const size = MpiWrapper::commSize();
  int const value = 17 + rank;
  stdVector< int > gathered( size, -1 );
  EXPECT_EQ( MpiWrapper::allgather( &value, 1, gathered.data(), 1, MPI_COMM_GEOS ), MPI_SUCCESS );
  for( int r = 0; r < size; ++r )
    EXPECT_EQ( gathered[r], 17 + r );

  stdVector< int > counts( size, 1 );
  stdVector< int > offsets( size );
  stdVector< int > sparse( 2 * size + 1, -1 );
  for( int r = 0; r < size; ++r )
    offsets[r] = 2 * r + 1;
  EXPECT_EQ( MpiWrapper::allgatherv( &value, 1, sparse.data(), counts.data(), offsets.data(), MPI_COMM_GEOS ), MPI_SUCCESS );
  for( int r = 0; r < size; ++r )
  {
    EXPECT_EQ( sparse[2 * r], -1 );
    EXPECT_EQ( sparse[2 * r + 1], 17 + r );
  }
  EXPECT_EQ( sparse.back(), -1 );

  EXPECT_EQ( MpiWrapper::gather( value, gathered, 0, MPI_COMM_GEOS ), MPI_SUCCESS );
  if( rank == 0 )
  {
    for( int r = 0; r < size; ++r )
      EXPECT_EQ( gathered[r], 17 + r );
  }
}

TEST( MpiWrapperTesting, PrefixSumUsesTheSuppliedCommunicator )
{
  int const rank = MpiWrapper::commRank();
  int const size = MpiWrapper::commSize();
  MPI_Comm reverse = MpiWrapper::commSplit( MPI_COMM_GEOS, 0, size - 1 - rank );
  EXPECT_EQ( MpiWrapper::prefixSum< globalIndex >( 1, reverse ), MpiWrapper::commRank( reverse ) );
  MpiWrapper::commFree( reverse );
}

TEST( MpiWrapperTesting, SparseExchangeDeliversOnlyAddressedBuffers )
{
  int const rank = MpiWrapper::commRank();
  int const size = MpiWrapper::commSize();
  // Each rank sends to itself, to the next rank and to the rank three ahead.
  // Buffer contents identify the sender, the receiver and a length that varies.
  auto payload = []( int from, int to )
  {
    stdVector< char > bytes( 3 + ( from + 2 * to ) % 5 );
    bytes[0] = static_cast< char >( from );
    bytes[1] = static_cast< char >( to );
    for( std::size_t i = 2; i < bytes.size(); ++i )
      bytes[i] = static_cast< char >( from + to + i );
    return bytes;
  };
  stdMap< int, stdVector< char > > outgoing;
  for( int const shift : { 0, 1, 3 } )
  {
    int const to = ( rank + shift ) % size;
    outgoing.get_inserted( to ) = payload( rank, to );
  }
  auto const incoming = mpi::sparseExchange( outgoing, MPI_COMM_GEOS, 7401, 2 );
  std::set< int > expected;
  for( int const shift : { 0, 1, 3 } )
    expected.insert( ( rank - shift % size + size ) % size );
  ASSERT_EQ( incoming.size(), expected.size() );
  for( int const from : expected )
  {
    auto const found = incoming.find( from );
    ASSERT_NE( found, incoming.end() ) << "from rank " << from;
    EXPECT_EQ( found->second, payload( from, rank ) ) << "from rank " << from;
  }
  // A rank without peers still takes part.
  stdMap< int, stdVector< char > > none;
  EXPECT_TRUE( mpi::sparseExchange( none, MPI_COMM_GEOS ).empty() );
}

int main( int argc, char * argv[] )
{
  MpiTestScope testScope{ argc, argv };
  return RUN_ALL_TESTS();
}
