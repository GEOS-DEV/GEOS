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

/**
 * @file DofManagerHelpers.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_DENSELAHELPERS_HPP
#define GEOSX_LINEARALGEBRA_DENSELAHELPERS_HPP

#include <random>

// BLAS and LAPACK function declaration
#ifdef USE_LAPACK
#include "BlasLapackFunctions.h"
#endif


namespace geosx
{

// unnamed namespace to avoid needless external linkage
namespace
{

// Random device and random number generator seed integer array used
// to populate a vector/matrix with random coefficients
static std::random_device rd;
static std::mt19937 gen( rd());
static std::uniform_int_distribution< int > dis( 0, 4095 );
static std::uniform_int_distribution< int > disOdd( 0, 2047 );
static int ISEED[] = {dis( gen ), dis( gen ), dis( gen ), disOdd( gen )*2 + 1};

// --------------------- HELPERS --------------------

template< int USD >
static real64 determinant3by3( arraySlice2d< real64 const, USD > const & A )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == 3 && A.size( 1 ) == 3,
                    "Matrix size must be 3x3" );

  real64 const det =
    A( 0, 0 ) * ( A( 1, 1 ) * A( 2, 2 ) - A( 1, 2 ) * A( 2, 1 ) ) +
    A( 0, 1 ) * ( A( 1, 2 ) * A( 2, 0 ) - A( 1, 0 ) * A( 2, 2 ) ) +
    A( 0, 2 ) * ( A( 1, 0 ) * A( 2, 1 ) - A( 1, 1 ) * A( 2, 0 ) );

  return det;
}

template< int USD >
static real64 determinant4by4( arraySlice2d< real64 const, USD > const & A )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == 4 && A.size( 1 ) == 4,
                    "Matrix size must be 4x4" );

  real64 const det =
    A( 0, 0 ) * ( A( 1, 1 ) * ( A( 2, 2 ) * A( 3, 3 ) - A( 3, 2 ) * A( 2, 3 ) ) -
                  A( 1, 2 ) * ( A( 2, 1 ) * A( 3, 3 ) - A( 3, 1 ) * A( 2, 3 ) ) +
                  A( 1, 3 ) * ( A( 2, 1 ) * A( 3, 2 ) - A( 3, 1 ) * A( 2, 2 ) )
                  ) -
    A( 0, 1 ) * ( A( 1, 0 ) * ( A( 2, 2 ) * A( 3, 3 ) - A( 3, 2 ) * A( 2, 3 ) ) -
                  A( 1, 2 ) * ( A( 2, 0 ) * A( 3, 3 ) - A( 3, 0 ) * A( 2, 3 ) ) +
                  A( 1, 3 ) * ( A( 2, 0 ) * A( 3, 2 ) - A( 3, 0 ) * A( 2, 2 ) )
                  ) +
    A( 0, 2 ) * ( A( 1, 0 ) * ( A( 2, 1 ) * A( 3, 3 ) - A( 3, 1 ) * A( 2, 3 ) ) -
                  A( 1, 1 ) * ( A( 2, 0 ) * A( 3, 3 ) - A( 3, 0 ) * A( 2, 3 ) ) +
                  A( 1, 3 ) * ( A( 2, 0 ) * A( 3, 1 ) - A( 3, 0 ) * A( 2, 1 ) )
                  ) -
    A( 0, 3 ) * ( A( 1, 0 ) * ( A( 2, 1 ) * A( 3, 2 ) - A( 3, 1 ) * A( 2, 2 ) ) -
                  A( 1, 1 ) * ( A( 2, 0 ) * A( 3, 2 ) - A( 3, 0 ) * A( 2, 2 ) ) +
                  A( 1, 2 ) * ( A( 2, 0 ) * A( 3, 1 ) - A( 3, 0 ) * A( 2, 1 ) )
                  );

  return det;
}

template< int USD >
static void matrixInverse2by2( arraySlice2d< real64 const, USD > const & A,
                               arraySlice2d< real64, USD > const & Ainv,
                               real64 const & oneOverDetA )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == 2 && A.size( 1 ) == 2,
                    "Matrix size must be 2x2" );

  Ainv( 0, 0 ) =  A( 1, 1 ) * oneOverDetA;
  Ainv( 0, 1 ) = -A( 0, 1 ) * oneOverDetA;
  Ainv( 1, 0 ) = -A( 1, 0 ) * oneOverDetA;
  Ainv( 1, 1 ) =  A( 0, 0 ) * oneOverDetA;
}

template< int USD >
static void matrixInverse3by3( arraySlice2d< real64 const, USD > const & A,
                               arraySlice2d< real64, USD > const & Ainv,
                               real64 const & oneOverDetA )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == 3 && A.size( 1 ) == 3,
                    "Matrix size must be 3x3" );

  Ainv( 0, 0 ) = ( A( 1, 1 ) * A( 2, 2 ) -
                   A( 1, 2 ) * A( 2, 1 ) ) * oneOverDetA;
  Ainv( 0, 1 ) = ( A( 0, 2 ) * A( 2, 1 ) -
                   A( 0, 1 ) * A( 2, 2 ) ) * oneOverDetA;
  Ainv( 0, 2 ) = ( A( 0, 1 ) * A( 1, 2 ) -
                   A( 0, 2 ) * A( 1, 1 ) ) * oneOverDetA;
  Ainv( 1, 0 ) = ( A( 1, 2 ) * A( 2, 0 ) -
                   A( 1, 0 ) * A( 2, 2 ) ) * oneOverDetA;
  Ainv( 1, 1 ) = ( A( 0, 0 ) * A( 2, 2 ) -
                   A( 0, 2 ) * A( 2, 0 ) ) * oneOverDetA;
  Ainv( 1, 2 ) = ( A( 0, 2 ) * A( 1, 0 ) -
                   A( 0, 0 ) * A( 1, 2 ) ) * oneOverDetA;
  Ainv( 2, 0 ) = ( A( 1, 0 ) * A( 2, 1 ) -
                   A( 1, 1 ) * A( 2, 0 ) ) * oneOverDetA;
  Ainv( 2, 1 ) = ( A( 0, 1 ) * A( 2, 0 ) -
                   A( 0, 0 ) * A( 2, 1 ) ) * oneOverDetA;
  Ainv( 2, 2 ) = ( A( 0, 0 ) * A( 1, 1 ) -
                   A( 0, 1 ) * A( 1, 0 ) ) * oneOverDetA;
}


template< typename MatType, typename VecType >
static real64 determinantFromLUFactors( MatType const & LUFactor,
                                        VecType const & ipiv,
                                        localIndex const offset )
{
  real64 det = 1.0;
  for( int i = 0; i < LUFactor.size( 0 ); ++i )
  {
    if( ipiv( i ) != i + offset )
    {
      det *= -LUFactor( i, i );
    }
    else
    {
      det *= LUFactor( i, i );
    }
  }
  return det;
}

#ifndef USE_LAPACK

// TODO: see if there is any value in keeping these functions, not sure

static void matrixFactorize( arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & LUFactor,
                             arraySlice1d< localIndex > const & ipiv )
{
  for( localIndex i = 0; i < LUFactor.size( 0 ); ++i )
  {

    // if needed, use pivoting
    localIndex piv = i;
    real64 a = std::fabs( LUFactor( i, i ) );
    for( localIndex j = i+1; j < LUFactor.size( 0 ); j++ )
    {
      real64 const b = std::fabs( LUFactor( j, i ) );
      if( b > a )
      {
        a = b;
        piv = j;
      }
    }
    ipiv( i ) = piv;

    // check if we have to swap rows i and piv
    if( piv != i )
    {
      for( int j = 0; j < LUFactor.size( 1 ); j++ )
      {
        real64 const tmp   = LUFactor( i, j );
        LUFactor( i, j )   = LUFactor( piv, j );
        LUFactor( piv, j ) = tmp;
      }
    }

    GEOSX_ASSERT_MSG( std::fabs( LUFactor( i, i ) )  > std::numeric_limits< real64 >::epsilon(),
                      "Found a zero pivot!" );

    // compute the entries of L and U
    real64 const A_ii_inv = 1.0 / LUFactor( i, i );
    for( localIndex j = i+1; j < LUFactor.size( 0 ); j++ )
    {
      LUFactor( j, i ) *= A_ii_inv;
    }
    for( localIndex k = i+1; k < LUFactor.size( 1 ); k++ )
    {
      real64 const A_ik = LUFactor( i, k );
      for( localIndex j = i+1; j < LUFactor.size( 0 ); j++ )
      {
        LUFactor( j, k ) -= A_ik * LUFactor( j, i );
      }
    }
  }
}


static void matrixTFactorize( arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & LUFactor,
                              arraySlice1d< localIndex > const & ipiv )
{
  for( localIndex i = 0; i < LUFactor.size( 1 ); ++i )
  {

    // if needed, use pivoting
    localIndex piv = i;
    real64 a = std::fabs( LUFactor( i, i ) );
    for( localIndex j = i+1; j < LUFactor.size( 1 ); j++ )
    {
      real64 const b = std::fabs( LUFactor( i, j ) );
      if( b > a )
      {
        a = b;
        piv = j;
      }
    }
    ipiv( i ) = piv;

    // check if we have to swap rows i and piv
    if( piv != i )
    {
      for( int j = 0; j < LUFactor.size( 0 ); j++ )
      {
        real64 const tmp   = LUFactor( j, i );
        LUFactor( j, i )   = LUFactor( j, piv );
        LUFactor( j, piv ) = tmp;
      }
    }

    GEOSX_ASSERT_MSG( std::fabs( LUFactor( i, i ) )  > std::numeric_limits< real64 >::epsilon(),
                      "Found a zero pivot!" );

    // compute the entries of L and U
    real64 const A_ii_inv = 1.0 / LUFactor( i, i );
    for( localIndex j = i+1; j < LUFactor.size( 1 ); j++ )
    {
      LUFactor( i, j ) *= A_ii_inv;
    }
    for( localIndex k = i+1; k < LUFactor.size( 0 ); k++ )
    {
      real64 const A_ki = LUFactor( k, i );
      for( localIndex j = i+1; j < LUFactor.size( 1 ); j++ )
      {
        LUFactor( k, j ) -= A_ki * LUFactor( i, j );
      }
    }
  }
}

static void matrixInverseFromLUFactors( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & LUFactor,
                                        arraySlice1d< localIndex const > const & ipiv,
                                        arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & Ainv )
{
  // TODO: add comments

  for( localIndex k = 0; k < Ainv.size( 0 ); k++ )
  {
    Ainv( k, k ) = 1.0 / LUFactor( k, k );
    real64 const minus_Ainv_kk = -Ainv( k, k );

    for( localIndex i = 0; i < k; i++ )
    {
      Ainv( i, k ) = LUFactor( i, k ) * minus_Ainv_kk;
    }
    for( localIndex j = k-1; j >= 0; j-- )
    {
      Ainv( j, k ) /= LUFactor( j, j );
      real64 const Ainv_jk = Ainv( j, k );
      for( localIndex i = 0; i < j; i++ )
      {
        Ainv( i, k ) -= LUFactor( i, j ) * Ainv_jk;
      }
    }
  }

  localIndex const l = Ainv.size( 0 ) - 1;
  for( localIndex j = 0; j < l; j++ )
  {
    real64 const minus_L_lj = -LUFactor( l, j );
    for( localIndex i = 0; i <= j; i++ )
    {
      Ainv( i, j ) += Ainv( i, l ) * minus_L_lj;
    }
    for( localIndex i = j+1; i < Ainv.size( 0 ); i++ )
    {
      Ainv( i, j ) = Ainv( i, l ) * minus_L_lj;
    }
  }

  for( localIndex k = Ainv.size( 0 ) - 2; k >= 0; k-- )
  {
    for( localIndex j = 0; j < k; j++ )
    {
      real64 const L_kj = LUFactor( k, j );
      for( localIndex i = 0; i < Ainv.size( 0 ); i++ )
      {
        Ainv( i, j ) -= Ainv( i, k ) * L_kj;
      }
    }
  }

  for( localIndex k = Ainv.size( 0 ) - 1; k >= 0; k-- )
  {
    localIndex const piv_k = ipiv( k );
    if( k != piv_k )
    {
      for( localIndex i = 0; i < Ainv.size( 0 ); i++ )
      {
        real64 const tmp = Ainv( i, k );
        Ainv( i, k ) = Ainv( i, piv_k );
        Ainv( i, piv_k ) = tmp;
      }
    }
  }
}


static void matrixTInverseFromLUFactors( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & LUFactor,
                                         arraySlice1d< localIndex const > const & ipiv,
                                         arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & Ainv )
{
  // TODO: add comments

  for( localIndex k = 0; k < Ainv.size( 0 ); k++ )
  {
    Ainv( k, k ) = 1.0 / LUFactor( k, k );
    real64 const minus_Ainv_kk = -Ainv( k, k );

    for( localIndex i = 0; i < k; i++ )
    {
      Ainv( k, i ) = LUFactor( k, i ) * minus_Ainv_kk;
    }
    for( localIndex j = k-1; j >= 0; j-- )
    {
      Ainv( k, j ) /= LUFactor( j, j );
      real64 const Ainv_kj = Ainv( k, j );
      for( localIndex i = 0; i < j; i++ )
      {
        Ainv( k, i ) -= LUFactor( j, i ) * Ainv_kj;
      }
    }
  }

  localIndex const l = Ainv.size( 0 ) - 1;
  for( localIndex j = 0; j < l; j++ )
  {
    real64 const minus_L_jl = -LUFactor( j, l );
    for( localIndex i = 0; i <= j; i++ )
    {
      Ainv( j, i ) += Ainv( l, i ) * minus_L_jl;
    }
    for( localIndex i = j+1; i < Ainv.size( 0 ); i++ )
    {
      Ainv( j, i ) = Ainv( l, i ) * minus_L_jl;
    }
  }

  for( localIndex k = Ainv.size( 0 ) - 2; k >= 0; k-- )
  {
    for( localIndex j = 0; j < k; j++ )
    {
      real64 const L_jk = LUFactor( j, k );
      for( localIndex i = 0; i < Ainv.size( 0 ); i++ )
      {
        Ainv( j, i ) -= Ainv( k, i ) * L_jk;
      }
    }
  }

  for( localIndex k = Ainv.size( 0 ) - 1; k >= 0; k-- )
  {
    localIndex const piv_k = ipiv( k );
    if( k != piv_k )
    {
      for( localIndex i = 0; i < Ainv.size( 0 ); i++ )
      {
        real64 const tmp = Ainv( k, i );
        Ainv( k, i ) = Ainv( piv_k, i );
        Ainv( piv_k, i ) = tmp;
      }
    }
  }
}

#endif

// --------------------- LAPACK WRAPPERS --------------------

// TODO: decide if we want to remove these LAPACK operations, they seem less important
// than those left in DenseLA.cpp

#ifdef USE_LAPACK

static void vectorNorm1Lapack( arraySlice1d< real64 const > const & X, real64 & norm )
{
  int const INCX = 1;
  int const N = integer_conversion< int >( X.size() );
  norm = GEOSX_dasum( &N, X.dataIfContiguous(), &INCX );
}

static void vectorNorm2Lapack( arraySlice1d< real64 const > const & X, real64 & norm )
{
  int const INCX = 1;
  int const N = integer_conversion< int >( X.size() );
  norm = GEOSX_dnrm2( &N, X.dataIfContiguous(), &INCX );
}

static void vectorNormInfLapack( arraySlice1d< real64 const > const & X, real64 & norm )
{
  int const INCX = 1;
  int const N = integer_conversion< int >( X.size() );
  int ind = GEOSX_idamax( &N, X.dataIfContiguous(), &INCX );
  ind -= 1; // Fortran convention, subtract 1
  norm = std::abs( X( ind ) );
}

template< int USD >
static void matrixNormLapack( arraySlice2d< real64 const, USD > const & A,
                              char const c,
                              bool const transpose,
                              real64 & norm )
{
  int const M = transpose
              ? integer_conversion< int >( A.size( 0 ) )
              : integer_conversion< int >( A.size( 1 ) );
  int const N = transpose
              ? integer_conversion< int >( A.size( 1 ) )
              : integer_conversion< int >( A.size( 0 ) );

  array1d< double > temp;
  double * work = nullptr;
  if( c == 'I' )
  {
    temp.resize( N );
    work = temp.data();
  }

  norm = GEOSX_dlange( &c, &N, &M, A.dataIfContiguous(), &N, work );
}

static void vectorVectorAddLapack( arraySlice1d< real64 const > const & X,
                                   arraySlice1d< real64 > const & Y,
                                   real64 const & alpha )
{
  int const INCX = 1;
  int const INCY = 1;
  int const N = integer_conversion< int >( X.size() );
  GEOSX_daxpy( &N, &alpha, X.dataIfContiguous(), &INCX, Y.dataIfContiguous(), &INCY );
}

template< int USD >
static void matrixMatrixAddLapack( arraySlice2d< real64 const, USD > const & A,
                                   arraySlice2d< real64, USD > const & B,
                                   real64 const & alpha )
{
  int const INCX = 1;
  int const INCY = 1;
  int const N = integer_conversion< int >( A.size() );
  GEOSX_daxpy( &N, &alpha, A.dataIfContiguous(), &INCX, B.dataIfContiguous(), &INCY );
}

template< int USD >
static void matrixScaleLapack( real64 const & alpha,
                               arraySlice2d< real64, USD > const & A )
{
  int const INCX = 1;
  int const N = integer_conversion< int >( A.size() );
  GEOSX_dscal( &N, &alpha, A.dataIfContiguous(), &INCX );
}

static void vectorScaleLapack( real64 const & alpha,
                               arraySlice1d< real64 > const & X )
{
  int const INCX = 1;
  int const N = integer_conversion< int >( X.size() );
  GEOSX_dscal( &N, &alpha, X.dataIfContiguous(), &INCX );
}

static void vectorDotLapack( arraySlice1d< real64 const > const & X,
                             arraySlice1d< real64 const > const & Y,
                             real64 & dotProduct )
{
  int const INCX = 1;
  int const INCY = 1;
  int const N = integer_conversion< int >( X.size() );
  dotProduct = GEOSX_ddot( &N, X.dataIfContiguous(), &INCX, Y.dataIfContiguous(), &INCY );
}

static void vectorCopyLapack( array1d< real64 > const & X,
                              array1d< real64 > & Y )
{
  int const INCX = 1;
  int const INCY = 1;
  int const N = integer_conversion< int >( X.size() );
  GEOSX_dcopy( &N, X.data(), &INCX, Y.data(), &INCY );
}

template< int USD >
static void matrixCopyLapack( arraySlice2d< real64 const, USD > const & A,
                              arraySlice2d< real64, USD > const & B )
{
  int const INCX = 1;
  int const INCY = 1;
  int const N = integer_conversion< int >( A.size() );
  GEOSX_dcopy( &N, A.dataIfContiguous(), &INCX, B.dataIfContiguous(), &INCY );
}

#endif

} // namespace

} // namespace geosx

#endif // GEOSX_LINEARALGEBRA_DENSELAHELPERS_HPP
