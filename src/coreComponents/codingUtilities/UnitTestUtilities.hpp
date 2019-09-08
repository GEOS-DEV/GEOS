/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * @file UnitTestUtilities.hpp
 */

#ifndef GEOSX_UNITTESTUTILITIES_HPP
#define GEOSX_UNITTESTUTILITIES_HPP

#include "gtest/gtest.h"

#include "common/DataTypes.hpp"

namespace geosx
{

namespace testing
{

template< typename T >
::testing::AssertionResult checkRelativeErrorFormat( const char *, const char *, const char *,
                                                     T v1, T v2, T relTol )
{
  T const delta = std::abs( v1 - v2 );
  T const value = std::max( std::abs( v1 ), std::abs( v2 ));

  if( delta < 1E-13 )
  {
    return ::testing::AssertionSuccess();
  }

  if( delta > relTol * value )
  {
    return ::testing::AssertionFailure() << std::scientific << std::setprecision( 5 )
                                         << " relative error: " << delta / value
                                         << " (" << v1 << " vs " << v2 << "),"
                                         << " exceeds " << relTol << std::endl;
  }
  return ::testing::AssertionSuccess();
}

template< typename T >
void checkRelativeError( T v1, T v2, T relTol )
{
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

template< typename T >
void checkRelativeError( T v1, T v2, T relTol, string const & name )
{
  SCOPED_TRACE( name );
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
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

#endif //GEOSX_UNITTESTUTILITIES_HPP
