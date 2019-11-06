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

/*
 * @file UnitTestUtilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP
#define GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP

#include "gtest/gtest.h"

#include "common/DataTypes.hpp"

namespace geosx
{

namespace testing
{

constexpr real64 DEFAULT_ABS_TOL = 1E-13;

::testing::AssertionResult checkRelativeErrorFormat( const char *, const char *, const char *, const char *,
                                                     real64 const v1, real64 const v2, real64 const relTol, real64 const absTol )
{
  real64 const delta = std::abs( v1 - v2 );
  real64 const value = std::max( std::abs(v1), std::abs(v2) );
  if (delta > absTol && delta > relTol * value)
  {
    return ::testing::AssertionFailure() << std::scientific << std::setprecision(5)
                                         << " relative error: " << delta / value
                                         << " (" << v1 << " vs " << v2 << "),"
                                         << " exceeds " << relTol << std::endl;
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
  SCOPED_TRACE(name);
  EXPECT_PRED_FORMAT4( checkRelativeErrorFormat, v1, v2, relTol, absTol );
}

void compareMatrixRow( globalIndex rowNumber, real64 relTol,
                       localIndex numRowEntries1, globalIndex * indices1, real64 * values1,
                       localIndex numRowEntries2, globalIndex * indices2, real64 * values2 )
{
  SCOPED_TRACE( "Row " + std::to_string( rowNumber ));

  EXPECT_EQ( numRowEntries1, numRowEntries2 );

  for( localIndex j1 = 0, j2 = 0 ; j1 < numRowEntries1 && j2 < numRowEntries2 ; ++j1, ++j2 )
  {
    while( j1 < numRowEntries1 && j2 < numRowEntries2 && indices1[j1] != indices1[j2] )
    {
      while( j1 < numRowEntries1 && indices1[j1] < indices2[j2] )
      {
        ADD_FAILURE() << "column " << indices1[j1] << ") in matrix 1 does not have a match";
      }
      while( j2 < numRowEntries2 && indices2[j2] < indices1[j1] )
      {
        ADD_FAILURE() << "column " << indices2[j2] << ") in matrix 2 does not have a match";
      }
    }
    if( j1 < numRowEntries1 && j2 < numRowEntries2 )
    {
      SCOPED_TRACE( "Column " + std::to_string( indices1[j1] ));

      checkRelativeError( values1[j1], values2[j1], relTol );
    }
  }
}

template< typename MATRIX >
void compareMatrices( MATRIX const & matrix1,
                      MATRIX const & matrix2,
                      real64 relTol )
{
  ASSERT_EQ( matrix1.localRows(), matrix2.localRows());

  array1d< globalIndex > indices1, indices2;
  array1d< real64 > values1, values2;

  // check the accuracy across local rows
  for( globalIndex i = matrix1.ilower() ; i < matrix1.iupper() ; ++i )
  {
    matrix1.getRowCopy( i, indices1, values1 );
    matrix2.getRowCopy( i, indices2, values2 );

    compareMatrixRow( i, relTol,
                      indices1.size(), indices1.data(), values1.data(),
                      indices2.size(), indices2.data(), values2.data() );
  }
}

}

}

#endif //GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP
