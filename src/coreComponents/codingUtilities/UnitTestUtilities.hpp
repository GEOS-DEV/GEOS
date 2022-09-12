/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file UnitTestUtilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP_
#define GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP_

#include "common/DataTypes.hpp"

#include "common/MpiWrapper.hpp"

#include <gtest/gtest.h>

#define SKIP_TEST_IF( COND, REASON )  \
  do                                  \
  {                                   \
    if( COND )                        \
    {                                 \
      GTEST_SKIP_( ": " REASON );     \
    }                                 \
  } while( false )

#define SKIP_TEST_IN_SERIAL( REASON ) \
  SKIP_TEST_IF( geosx::MpiWrapper::commSize() == 1, REASON )

#define SKIP_TEST_IN_PARALLEL( REASON ) \
  SKIP_TEST_IF( geosx::MpiWrapper::commSize() != 1, REASON )

namespace geosx
{

namespace testing
{

/**
 * @brief Returns the expected value depending on the MPI context.
 * @tparam T Type of the expected value.
 * @param[in] expectedSerial Expected value for serial case.
 * @param[in] expectedParallel Expected values for parallel MPI cases, for the current MPI rank.
 *        The length of the list should match the size of the MPI communicator @p comm.
 *        The @p i^th element of the list will be chosen for MPI rank @p i.
 * @param[in] comm The MPI_Comm communicator that the function will act on.
 * @return The expected value.
 *
 * @note This function is meant to be used to run the same test in serial or parallel environments.
 */
template< typename T >
T expected( T expectedSerial,
            std::initializer_list< T > expectedParallel,
            MPI_Comm const & comm = MPI_COMM_GEOSX )
{
  int const mpiSize = MpiWrapper::commSize( comm );
  if( mpiSize == 1 )
  {
    return expectedSerial;
  }
  else
  {
    GEOSX_ASSERT( expectedParallel.size() == std::size_t( mpiSize ) );
    std::vector< T > tmp( expectedParallel );
    return tmp[MpiWrapper::commRank( comm )];
  }
}


constexpr real64 DEFAULT_ABS_TOL = 1E-12;
constexpr real64 DEFAULT_REL_TOL = std::numeric_limits< real64 >::epsilon();

::testing::AssertionResult checkRelativeErrorFormat( const char *, const char *, const char *, const char *,
                                                     real64 const v1, real64 const v2, real64 const relTol, real64 const absTol )
{
  real64 const delta = std::abs( v1 - v2 );
  real64 const value = std::max( std::abs( v1 ), std::abs( v2 ) );
  if( delta > absTol && delta > relTol * value )
  {
    return ::testing::AssertionFailure() << std::scientific << std::setprecision( 5 )
                                         << " relative error: " << delta / value
                                         << " (" << v1 << " vs " << v2 << "),"
                                         << " exceeds " << relTol <<". "
                                         << " absolute error: " << delta << " exeeds "
                                         << absTol <<std::endl;
  }
  return ::testing::AssertionSuccess();
}

::testing::AssertionResult checkRelativeErrorFormat( char const * a, char const * b, char const * c,
                                                     real64 const v1, real64 const v2, real64 const relTol )
{ return checkRelativeErrorFormat( a, b, c, "DEFAULT_ABS_TOL", v1, v2, relTol, DEFAULT_ABS_TOL ); }

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol )
{ EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol ); }

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol, real64 const absTol )
{ EXPECT_PRED_FORMAT4( checkRelativeErrorFormat, v1, v2, relTol, absTol ); }

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol, string const & name )
{
  SCOPED_TRACE( name );
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol, real64 const absTol, string const & name )
{
  SCOPED_TRACE( name );
  EXPECT_PRED_FORMAT4( checkRelativeErrorFormat, v1, v2, relTol, absTol );
}

template< typename COL_INDEX, typename VALUE >
void compareMatrixRow( VALUE const relTol, VALUE const absTol,
                       localIndex const length1, COL_INDEX const * const indices1, VALUE const * const values1,
                       localIndex const length2, COL_INDEX const * const indices2, VALUE const * const values2 )
{
  EXPECT_EQ( length1, length2 );

  for( localIndex j1 = 0, j2 = 0; j1 < length1 && j2 < length2; ++j1, ++j2 )
  {
    while( j1 < length1 && j2 < length2 && indices1[j1] != indices1[j2] )
    {
      while( j1 < length1 && indices1[j1] < indices2[j2] )
      {
        ADD_FAILURE() << "column " << indices1[j1] << ") in matrix 1 does not have a match";
      }
      while( j2 < length2 && indices2[j2] < indices1[j1] )
      {
        ADD_FAILURE() << "column " << indices2[j2] << ") in matrix 2 does not have a match";
      }
    }
    if( j1 < length1 && j2 < length2 )
    {
      SCOPED_TRACE( "Column " + std::to_string( indices1[j1] ));

      checkRelativeError( values1[j1], values2[j1], relTol, absTol );
    }
  }
}

template< typename T, typename COL_INDEX >
void compareLocalMatrices( CRSMatrixView< T const, COL_INDEX const > const & matrix1,
                           CRSMatrixView< T const, COL_INDEX const > const & matrix2,
                           real64 const relTol = DEFAULT_REL_TOL,
                           real64 const absTol = DEFAULT_ABS_TOL,
                           globalIndex const rowOffset = 0 )
{
  ASSERT_EQ( matrix1.numRows(), matrix2.numRows() );
  ASSERT_EQ( matrix1.numColumns(), matrix2.numColumns() );

  matrix1.move( LvArray::MemorySpace::host, false );
  matrix2.move( LvArray::MemorySpace::host, false );

  // check the accuracy across local rows
  for( localIndex i = 0; i < matrix1.numRows(); ++i )
  {
    SCOPED_TRACE( GEOSX_FMT( "Row {}", i + rowOffset ) );
    compareMatrixRow( relTol, absTol,
                      matrix1.numNonZeros( i ),
                      matrix1.getColumns( i ).dataIfContiguous(),
                      matrix1.getEntries( i ).dataIfContiguous(),
                      matrix2.numNonZeros( i ),
                      matrix2.getColumns( i ).dataIfContiguous(),
                      matrix2.getEntries( i ).dataIfContiguous() );
  }
}

template< typename MATRIX >
void compareMatrices( MATRIX const & matrix1,
                      MATRIX const & matrix2,
                      real64 const relTol = DEFAULT_REL_TOL,
                      real64 const absTol = DEFAULT_ABS_TOL )
{
  ASSERT_EQ( matrix1.numGlobalRows(), matrix2.numGlobalRows() );
  ASSERT_EQ( matrix1.numGlobalCols(), matrix2.numGlobalCols() );

  ASSERT_EQ( matrix1.numLocalRows(), matrix2.numLocalRows() );
  ASSERT_EQ( matrix1.numLocalCols(), matrix2.numLocalCols() );

  CRSMatrix< real64, globalIndex > const mat1 = matrix1.extract();
  CRSMatrix< real64, globalIndex > const mat2 = matrix2.extract();

  compareLocalMatrices( mat1.toViewConst(), mat2.toViewConst(), relTol, absTol, matrix1.ilower() );
}

} // namespace testing

} // namespace geosx

#endif //GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP_
