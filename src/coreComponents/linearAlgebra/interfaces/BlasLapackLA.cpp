/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BlasLapack.cpp
 */

// Include the corresponding header file.
#include "BlasLapackLA.hpp"

#include <random>

// BLAS and LAPACK function declaration
#include "BlasLapackFunctions.h"

// Put everything under the geosx namespace.
namespace geosx
{

// Random device and random number generator seed integer array used
// to populate a vector/matrix with random coefficients
static std::random_device rd;
static std::mt19937 gen( rd());
static std::uniform_int_distribution< int > dis( 0, 4095 );
static std::uniform_int_distribution< int > disOdd( 0, 2047 );
static int ISEED[] = {dis( gen ), dis( gen ), dis( gen ), disOdd( gen )*2 + 1};

real64 BlasLapackLA::vectorNorm1( arraySlice1d< real64 const > const & x )
{
  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( x.size() );
  return GEOSX_dasum( &N, x.dataIfContiguous(), &INCX );
}

real64 BlasLapackLA::vectorNorm2( arraySlice1d< real64 const > const & x )
{
  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( x.size() );
  return GEOSX_dnrm2( &N, x.dataIfContiguous(), &INCX );
}

real64 BlasLapackLA::vectorNormInf( arraySlice1d< real64 const > const & x )
{
  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( x.size() );
  int ind = GEOSX_idamax( &N, x.dataIfContiguous(), &INCX );
  ind -= 1; // Fortran convention, subtract 1
  return std::abs( x( ind ) );
}

namespace detail
{

template< int USD >
static real64 determinant( arraySlice2d< real64 const, USD > const & a )
{
  // --- check that matrix is square
  GEOSX_ASSERT_MSG( a.size( 0 ) == a.size( 1 ) && a.size( 0 ) > 0,
                    "Matrix must be square with order greater than zero" );

  real64 det;
  switch( a.size( 0 ) )
  {
    case 1:
    {
      det = a( 0, 0 );
      break;
    }
    case 2:
    {
      det = a( 0, 0 ) * a( 1, 1 ) - a( 0, 1 ) * a( 1, 0 );
      break;
    }
    case 3:
    {
      det =
        a( 0, 0 ) * ( a( 1, 1 ) * a( 2, 2 ) - a( 1, 2 ) * a( 2, 1 ) ) +
        a( 0, 1 ) * ( a( 1, 2 ) * a( 2, 0 ) - a( 1, 0 ) * a( 2, 2 ) ) +
        a( 0, 2 ) * ( a( 1, 0 ) * a( 2, 1 ) - a( 1, 1 ) * a( 2, 0 ) );
      break;
    }
    case 4:
    {
      det =
        a( 0, 0 ) * ( a( 1, 1 ) * ( a( 2, 2 ) * a( 3, 3 ) - a( 3, 2 ) * a( 2, 3 ) ) -
                      a( 1, 2 ) * ( a( 2, 1 ) * a( 3, 3 ) - a( 3, 1 ) * a( 2, 3 ) ) +
                      a( 1, 3 ) * ( a( 2, 1 ) * a( 3, 2 ) - a( 3, 1 ) * a( 2, 2 ) )
                      ) -
        a( 0, 1 ) * ( a( 1, 0 ) * ( a( 2, 2 ) * a( 3, 3 ) - a( 3, 2 ) * a( 2, 3 ) ) -
                      a( 1, 2 ) * ( a( 2, 0 ) * a( 3, 3 ) - a( 3, 0 ) * a( 2, 3 ) ) +
                      a( 1, 3 ) * ( a( 2, 0 ) * a( 3, 2 ) - a( 3, 0 ) * a( 2, 2 ) )
                      ) +
        a( 0, 2 ) * ( a( 1, 0 ) * ( a( 2, 1 ) * a( 3, 3 ) - a( 3, 1 ) * a( 2, 3 ) ) -
                      a( 1, 1 ) * ( a( 2, 0 ) * a( 3, 3 ) - a( 3, 0 ) * a( 2, 3 ) ) +
                      a( 1, 3 ) * ( a( 2, 0 ) * a( 3, 1 ) - a( 3, 0 ) * a( 2, 1 ) )
                      ) -
        a( 0, 3 ) * ( a( 1, 0 ) * ( a( 2, 1 ) * a( 3, 2 ) - a( 3, 1 ) * a( 2, 2 ) ) -
                      a( 1, 1 ) * ( a( 2, 0 ) * a( 3, 2 ) - a( 3, 0 ) * a( 2, 2 ) ) +
                      a( 1, 2 ) * ( a( 2, 0 ) * a( 3, 1 ) - a( 3, 0 ) * a( 2, 1 ) )
                      );
      break;
    }
    default:
    {
      // Compute the determinant via LU factorization
      int const NN = LvArray::integerConversion< int >( a.size( 0 ) );
      int INFO;
      array1d< int > IPIV( NN );

      array2d< double > LUFactor( a.size( 0 ), a.size( 1 ) );

      for( localIndex i = 0; i < a.size( 0 ); ++i )
      {
        for( localIndex j = 0; j < a.size( 1 ); ++j )
        {
          LUFactor( i, j ) = a( i, j );
        }
      }

      // We compute the LU factors for the transpose matrix, i.e. choosing the
      // LAPACK_COL_MAJOR ordering, to avoid transposition/copy requires for
      // LAPACK_ROW_MAJOR ordering.
      GEOSX_dgetrf( &NN, &NN, LUFactor.data(), &NN, IPIV.data(), &INFO );

      GEOSX_ASSERT_MSG( INFO == 0, "LAPACK dgetrf error code: " << INFO );

      det = 1.0;
      for( int i = 0; i < NN; ++i )
      {
        if( IPIV[i] != i + 1 ) //IPIV is based on Fortran convention (counting from 1)
        {
          det *= -LUFactor( i, i );
        }
        else
        {
          det *= LUFactor( i, i );
        }
      }

      break;
    }
  }

  return det;
}

template< int USD >
static real64 matrixNorm( arraySlice2d< real64 const, USD > const & a,
                          char const norm )
{
  int const M = LvArray::integerConversion< int >( a.size( 0 ) );
  int const N = LvArray::integerConversion< int >( a.size( 1 ) );

  array1d< double > temp;
  double * WORK = nullptr;
  if( norm == 'I' )
  {
    temp.resize( N );
    WORK = temp.data();
  }

  return GEOSX_dlange( &norm, &N, &M, a.dataIfContiguous(), &N, WORK );
}

template< int USD >
void matrixMatrixAdd( arraySlice2d< real64 const, USD > const & a,
                      arraySlice2d< real64, USD > const & b,
                      real64 const alpha )
{

  GEOSX_ASSERT_MSG( a.size( 0 ) == b.size( 0 ) &&
                    a.size( 1 ) == b.size( 1 ),
                    "Matrix dimensions not compatible for sum" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = LvArray::integerConversion< int >( a.size() );
  GEOSX_daxpy( &N, &alpha, a.dataIfContiguous(), &INCX, b.dataIfContiguous(), &INCY );
}

template< int USD >
void matrixScale( real64 const alpha,
                  arraySlice2d< real64, USD > const & a )
{
  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( a.size() );
  GEOSX_dscal( &N, &alpha, a.dataIfContiguous(), &INCX );
}

template< int USD >
void matrixRand( arraySlice2d< real64, USD > const & a,
                 BlasLapackLA::RandomNumberDistribution const & idist )
{
  int const IDIST = static_cast< int >(idist);
  int const NN = LvArray::integerConversion< int >( a.size() );
  GEOSX_ASSERT_MSG( NN > 0, "The matrix cannot be empty" );
  GEOSX_dlarnv( &IDIST, ISEED, &NN, a.dataIfContiguous() );
}

template< int USD >
void matrixCopy( arraySlice2d< real64 const, USD > const & a,
                 arraySlice2d< real64, USD > const & b )
{
  GEOSX_ASSERT_MSG( a.size( 0 ) == b.size( 0 ) &&
                    a.size( 1 ) == b.size( 1 ),
                    "Matrix dimensions not compatible for copying" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = LvArray::integerConversion< int >( a.size() );
  GEOSX_dcopy( &N, a.dataIfContiguous(), &INCX, b.dataIfContiguous(), &INCY );
}

template< int USD >
void matrixInverse( arraySlice2d< real64 const, USD > const & a,
                    arraySlice2d< real64, USD > const & ainv,
                    real64 & detA )
{
  // --- Check that source matrix is square
  int const NN = LvArray::integerConversion< int >( a.size( 0 ));
  GEOSX_ASSERT_MSG( NN > 0 &&
                    NN == a.size( 1 ),
                    "Matrix must be square" );

  // --- Check that inverse matrix has appropriate dimension
  GEOSX_ASSERT_MSG( ainv.size( 0 ) == NN &&
                    ainv.size( 1 ) == NN,
                    "Inverse matrix has wrong dimensions" );

  // --- Check if matrix is singular by computing the determinant
  //     note: if order greater than 3 we compute the determinant by
  //           first constructing the LU factors, later reused for calculating
  //           the inverse.
  array1d< int > IPIV;
  array1d< double > INV_WORK;

  if( NN <= 3 )
  {
    detA = determinant( a.toSliceConst());
  }
  else
  {
    // Copy A in Ainv
    matrixCopy( a, ainv );

    // Declare workspace for permutations and scratch array
    IPIV.resize( NN );
    INV_WORK.resize( NN );

    // Compute determinant (not done calling directly the function determinant
    // (avoid computing twice LUFactors, currently stored in Ainv, needed for
    // computing the inverse). Again we compute the LU factors for the
    // transpose matrix, i.e. choosing the LAPACK_COL_MAJOR ordering, to
    // avoid transposition/copy requires for LAPACK_ROW_MAJOR ordering.
    int INFO;
    GEOSX_dgetrf( &NN, &NN, ainv.dataIfContiguous(), &NN, IPIV.data(), &INFO );

    GEOSX_ASSERT_MSG( INFO == 0, "LAPACK dgetrf error code: " << INFO );

    detA = 1.0;
    for( int i = 0; i < NN; ++i )
    {
      if( IPIV[i] != i + 1 ) //IPIV is based on Fortran convention (counting from 1)
      {
        detA *= -ainv( i, i );
      }
      else
      {
        detA *= ainv( i, i );
      }
    }
  }

  // Check if matrix is singular
  GEOSX_ASSERT_MSG( std::abs( detA ) > 0, "Matrix is singular" );

  real64 oneOverDetA = 1. / detA;

  // --- Compute inverse
  switch( NN )
  {
    case 1:
    {
      ainv( 0, 0 ) = oneOverDetA;
      break;
    }

    case 2:
    {
      ainv( 0, 0 ) =  a( 1, 1 ) * oneOverDetA;
      ainv( 0, 1 ) = -a( 0, 1 ) * oneOverDetA;
      ainv( 1, 0 ) = -a( 1, 0 ) * oneOverDetA;
      ainv( 1, 1 ) =  a( 0, 0 ) * oneOverDetA;
      break;
    }

    case 3:
    {
      ainv( 0, 0 ) = ( a( 1, 1 ) * a( 2, 2 ) -
                       a( 1, 2 ) * a( 2, 1 ) ) * oneOverDetA;
      ainv( 0, 1 ) = ( a( 0, 2 ) * a( 2, 1 ) -
                       a( 0, 1 ) * a( 2, 2 ) ) * oneOverDetA;
      ainv( 0, 2 ) = ( a( 0, 1 ) * a( 1, 2 ) -
                       a( 0, 2 ) * a( 1, 1 ) ) * oneOverDetA;
      ainv( 1, 0 ) = ( a( 1, 2 ) * a( 2, 0 ) -
                       a( 1, 0 ) * a( 2, 2 ) ) * oneOverDetA;
      ainv( 1, 1 ) = ( a( 0, 0 ) * a( 2, 2 ) -
                       a( 0, 2 ) * a( 2, 0 ) ) * oneOverDetA;
      ainv( 1, 2 ) = ( a( 0, 2 ) * a( 1, 0 ) -
                       a( 0, 0 ) * a( 1, 2 ) ) * oneOverDetA;
      ainv( 2, 0 ) = ( a( 1, 0 ) * a( 2, 1 ) -
                       a( 1, 1 ) * a( 2, 0 ) ) * oneOverDetA;
      ainv( 2, 1 ) = ( a( 0, 1 ) * a( 2, 0 ) -
                       a( 0, 0 ) * a( 2, 1 ) ) * oneOverDetA;
      ainv( 2, 2 ) = ( a( 0, 0 ) * a( 1, 1 ) -
                       a( 0, 1 ) * a( 1, 0 ) ) * oneOverDetA;
      break;
    }
    default:
    {
      // Invert (LAPACK function DGETRI). The LU factors computed for the
      // transpose matrix stored in Ainv are used.
      int INFO;
      GEOSX_dgetri( &NN, ainv.dataIfContiguous(), &NN, IPIV.data(), INV_WORK.data(), &NN, &INFO );

      GEOSX_ASSERT_MSG( INFO == 0, "LAPACK dgetri error code: " << INFO );

      break;
    }
  }
}

template< int USD >
void matrixInverse( arraySlice2d< real64 const, USD > const & a,
                    arraySlice2d< real64, USD > const & ainv )
{
  real64 detA;
  matrixInverse( a, ainv, detA );
}

} // namespace detail

real64 BlasLapackLA::determinant( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a )
{
  return detail::determinant( a );
}

real64 BlasLapackLA::determinant( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a )
{
  return detail::determinant( a );
}

real64 BlasLapackLA::matrixNormInf( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a )
{
  // For row-major, computed as one-norm of the transpose matrix
  return detail::matrixNorm( a, '1' );
}

real64 BlasLapackLA::matrixNormInf( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a )
{
  return detail::matrixNorm( a, 'I' );
}

real64 BlasLapackLA::matrixNorm1( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a )
{
  // For row-major, computed as infinity-norm of the transpose matrix
  return detail::matrixNorm( a, 'I' );
}

real64 BlasLapackLA::matrixNorm1( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a )
{
  return detail::matrixNorm( a, '1' );
}

real64 BlasLapackLA::matrixNormFrobenius( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a )
{
  // For row-major, computed using the transpose matrix
  return detail::matrixNorm( a, 'F' );
}

real64 BlasLapackLA::matrixNormFrobenius( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a )
{
  // For row-major, computed using the transpose matrix
  return detail::matrixNorm( a, 'F' );
}

void BlasLapackLA::vectorVectorAdd( arraySlice1d< real64 const > const & x,
                                    arraySlice1d< real64 > const & y,
                                    real64 const alpha )
{
  GEOSX_ASSERT_MSG( x.size() == y.size(),
                    "Vector dimensions not compatible for sum" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = LvArray::integerConversion< int >( x.size() );
  GEOSX_daxpy( &N, &alpha, x.dataIfContiguous(), &INCX, y.dataIfContiguous(), &INCY );
}

void BlasLapackLA::matrixMatrixAdd( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                    arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & b,
                                    real64 const alpha )
{
  detail::matrixMatrixAdd( a, b, alpha );
}

void BlasLapackLA::matrixMatrixAdd( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a,
                                    arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & b,
                                    real64 const alpha )
{
  detail::matrixMatrixAdd( a, b, alpha );
}

void BlasLapackLA::vectorScale( real64 const alpha,
                                arraySlice1d< real64 > const & x )
{
  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( x.size() );
  GEOSX_dscal( &N, &alpha, x.dataIfContiguous(), &INCX );
}

void BlasLapackLA::matrixScale( real64 const alpha, arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & a )
{
  detail::matrixScale( alpha, a );
}

void BlasLapackLA::matrixScale( real64 const alpha, arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & a )
{
  detail::matrixScale( alpha, a );
}

real64 BlasLapackLA::vectorDot( arraySlice1d< real64 const > const & x,
                                arraySlice1d< real64 const > const & y )
{
  GEOSX_ASSERT_MSG( x.size() == y.size(), "Vector dimensions not compatible for dot product" );
  int const INCX = 1;
  int const INCY = 1;
  int const N = LvArray::integerConversion< int >( x.size() );
  return GEOSX_ddot( &N, x.dataIfContiguous(), &INCX, y.dataIfContiguous(), &INCY );

}

void BlasLapackLA::matrixVectorMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                         arraySlice1d< real64 const > const & x,
                                         arraySlice1d< real64 > const & y,
                                         real64 const alpha,
                                         real64 const beta )
{
  GEOSX_ASSERT_MSG( a.size( 1 ) == x.size() && a.size( 0 ) == y.size(),
                    "Matrix, source vector and destination vector not compatible" );

  int const M = LvArray::integerConversion< int >( a.size( 0 ) );
  int const N = 1;
  int const K = LvArray::integerConversion< int >( a.size( 1 ) );

  // A*X = Y is computed as X^T * A^T = Y^T, i.e. accessing the transpose
  // matrix using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, x.dataIfContiguous(), &N, a.dataIfContiguous(), &K, &beta, y.dataIfContiguous(), &N );
}

void BlasLapackLA::matrixTVectorMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                          arraySlice1d< real64 const > const & x,
                                          arraySlice1d< real64 > const & y,
                                          real64 const alpha,
                                          real64 const beta )
{
  GEOSX_ASSERT_MSG( a.size( 0 ) == x.size() && a.size( 1 ) == y.size(),
                    "Matrix, source vector and destination vector not compatible" );

  int const M = LvArray::integerConversion< int >( a.size( 1 ) );
  int const N = 1;
  int const K = LvArray::integerConversion< int >( a.size( 0 ) );

  // A^T*X = Y is computed as X^T * A = Y^T, i.e. accessing the transpose
  // matrix using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'T';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, x.dataIfContiguous(), &N, a.dataIfContiguous(), &M, &beta, y.dataIfContiguous(), &N );
}

void BlasLapackLA::matrixMatrixMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                         arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & b,
                                         arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & c,
                                         real64 const alpha,
                                         real64 const beta )
{

  GEOSX_ASSERT_MSG( c.size( 0 ) == a.size( 0 ) &&
                    c.size( 1 ) == b.size( 1 ) &&
                    a.size( 1 ) == b.size( 0 ),
                    "Matrix dimensions not compatible for product" );

  int const M = LvArray::integerConversion< int >( a.size( 0 ) );
  int const N = LvArray::integerConversion< int >( b.size( 1 ) );
  int const K = LvArray::integerConversion< int >( a.size( 1 ) );

  // A*B = C is computed as B^T * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, b.dataIfContiguous(), &N, a.dataIfContiguous(), &K, &beta, c.dataIfContiguous(), &N );
}

void BlasLapackLA::matrixTMatrixMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                          arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & b,
                                          arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & c,
                                          real64 const alpha,
                                          real64 const beta )
{

  GEOSX_ASSERT_MSG( c.size( 0 ) == a.size( 1 ) &&
                    c.size( 1 ) == b.size( 1 ) &&
                    a.size( 0 ) == b.size( 0 ),
                    "Matrix dimensions not compatible for product" );

  int const M = LvArray::integerConversion< int >( a.size( 1 ) );
  int const N = LvArray::integerConversion< int >( b.size( 1 ) );
  int const K = LvArray::integerConversion< int >( a.size( 0 ) );

  // A^T*B = C is computed as B^T * A = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'N';
  char const TRANS2 = 'T';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, b.dataIfContiguous(), &N, a.dataIfContiguous(), &M, &beta, c.dataIfContiguous(), &N );
}

void BlasLapackLA::matrixMatrixTMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                          arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & b,
                                          arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & c,
                                          real64 const alpha,
                                          real64 const beta )
{

  GEOSX_ASSERT_MSG( c.size( 0 ) == a.size( 0 ) &&
                    c.size( 1 ) == b.size( 0 ) &&
                    a.size( 1 ) == b.size( 1 ),
                    "Matrix dimensions not compatible for product" );

  int const M = LvArray::integerConversion< int >( a.size( 0 ) );
  int const N = LvArray::integerConversion< int >( b.size( 0 ) );
  int const K = LvArray::integerConversion< int >( a.size( 1 ) );

  // A*B^T = C is computed as B * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'T';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, b.dataIfContiguous(), &K, a.dataIfContiguous(), &K, &beta, c.dataIfContiguous(), &N );
}

void BlasLapackLA::matrixTMatrixTMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                           arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & b,
                                           arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & c,
                                           real64 const alpha,
                                           real64 const beta )
{

  GEOSX_ASSERT_MSG( c.size( 0 ) == a.size( 1 ) &&
                    c.size( 1 ) == b.size( 0 ) &&
                    a.size( 0 ) == b.size( 1 ),
                    "Matrix dimensions not compatible for product" );

  int const M = LvArray::integerConversion< int >( a.size( 1 ) );
  int const N = LvArray::integerConversion< int >( b.size( 0 ) );
  int const K = LvArray::integerConversion< int >( a.size( 0 ) );

  // A^T*B^T = C is computed as B * A = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'T';
  char const TRANS2 = 'T';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, b.dataIfContiguous(), &K, a.dataIfContiguous(), &M, &beta, c.dataIfContiguous(), &N );

  return;
}

void BlasLapackLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                  arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & ainv,
                                  real64 & detA )
{
  detail::matrixInverse( a, ainv, detA );
}

void BlasLapackLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a,
                                  arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & ainv,
                                  real64 & detA )
{
  detail::matrixInverse( a, ainv, detA );
}

void BlasLapackLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                                  arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & ainv )
{
  detail::matrixInverse( a, ainv );
}

void BlasLapackLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a,
                                  arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & ainv )
{
  detail::matrixInverse( a, ainv );
}

void BlasLapackLA::vectorCopy( arraySlice1d< real64 const > const & x,
                               arraySlice1d< real64 > const & y )
{
  GEOSX_ASSERT_MSG( x.size() == y.size(),
                    "Vector dimensions not compatible for copying" );

  int const INCX = 1;
  int const INCY = 1;
  int const N = LvArray::integerConversion< int >( x.size() );
  GEOSX_dcopy( &N, x.dataIfContiguous(), &INCX, y.dataIfContiguous(), &INCY );
}

void BlasLapackLA::matrixCopy( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                               arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & b )
{
  return detail::matrixCopy( a, b );
}

void BlasLapackLA::matrixCopy( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a,
                               arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & b )
{
  return detail::matrixCopy( a, b );
}

void BlasLapackLA::setRandomNumberGeneratorSeed( arraySlice1d< int const > const & seed )
{
  // Error checking
  GEOSX_ASSERT_MSG( seed.size() >= 4, "Seed array must have size at least four" );

  GEOSX_ASSERT_MSG( 0 <= seed( 0 ) && seed( 0 ) < 4096 &&
                    0 <= seed( 1 ) && seed( 1 ) < 4096 &&
                    0 <= seed( 2 ) && seed( 2 ) < 4096 &&
                    0 <= seed( 3 ) && seed( 3 ) < 4096,
                    "Seed array integer entries must be in interval [0,4096)" );

  GEOSX_ASSERT_MSG( seed( 3 ) % 2 > 0, "Seed array 4th element must be odd" );

  for( int i = 0; i < 4; ++i )
  {
    ISEED[i] = seed[i];
  }
}

void BlasLapackLA::getRandomNumberGeneratorSeed( arraySlice1d< int > const & seed )
{
  // Error checking
  GEOSX_ASSERT_MSG( seed.size() >= 4, "Seed array must have size at least four" );
  for( int i = 0; i < 4; ++i )
  {
    seed[i] = ISEED[i];
  }
}

void BlasLapackLA::vectorRand( arraySlice1d< real64 > const & x,
                               RandomNumberDistribution const & idist )
{

  int IDIST = static_cast< int >(idist);
  int const N = LvArray::integerConversion< int >( x.size() );
  GEOSX_ASSERT_MSG( N > 0, "The vector cannot be empty" );
  GEOSX_dlarnv( &IDIST, ISEED, &N, x.dataIfContiguous());
}

void BlasLapackLA::matrixRand( arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & a,
                               RandomNumberDistribution const & idist )
{
  detail::matrixRand( a, idist );
}

void BlasLapackLA::matrixRand( arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & a,
                               RandomNumberDistribution const & idist )
{
  detail::matrixRand( a, idist );
}

void BlasLapackLA::matrixSVD( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & a,
                              arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & u,
                              arraySlice1d< real64 > const & s,
                              arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & vt )
{
  int const minDim = (a.size( 0 ) < a.size( 1 ))
                   ? LvArray::integerConversion< int >( a.size( 0 ) )
                   : LvArray::integerConversion< int >( a.size( 1 ) );

  GEOSX_ASSERT_MSG( a.size( 0 ) == u.size( 0 ) && minDim == u.size( 1 ),
                    "The matrices A and U have an incompatible size" );

  GEOSX_ASSERT_MSG( minDim == vt.size( 0 ) && a.size( 1 ) == vt.size( 1 ),
                    "The matrices A and V have an incompatible size" );

  GEOSX_ASSERT_MSG( s.size() == minDim,
                    "The matrix A and vector S have an incompatible size" );

  // make a copy of A, since dgesvd destroys contents
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > ACOPY( a.size( 0 ), a.size( 1 ) );
  for( int i = 0; i < a.size( 0 ); ++i )
  {
    for( int j = 0; j < a.size( 1 ); ++j )
    {
      ACOPY( i, j ) = a( i, j );
    }
  }

  // define the arguments of dgesvd
  int const M     = LvArray::integerConversion< int >( a.size( 0 ) );
  int const N     = LvArray::integerConversion< int >( a.size( 1 ) );
  int const LDA   = M;
  int const LDU   = M;
  int const LDVT  = minDim;
  int LWORK = 0;
  int INFO  = 0;
  double WKOPT = 0.0;

  // 1) query and allocate the optimal workspace
  LWORK = -1;
  GEOSX_dgesvd( "S", "S",
                &M, &N, ACOPY.data(), &LDA,
                s.dataIfContiguous(), u.dataIfContiguous(), &LDU, vt.dataIfContiguous(), &LDVT,
                &WKOPT, &LWORK, &INFO );

  LWORK = static_cast< int >( WKOPT );
  array1d< real64 > WORK( LWORK );

  // 2) compute svd
  GEOSX_dgesvd( "S", "S",
                &M, &N, ACOPY.data(), &LDA,
                s.dataIfContiguous(), u.dataIfContiguous(), &LDU, vt.dataIfContiguous(), &LDVT,
                WORK.data(), &LWORK, &INFO );

  GEOSX_ASSERT_MSG( INFO == 0, "The algorithm computing SVD failed to converge." );
}

void BlasLapackLA::matrixSVD( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & a,
                              arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & u,
                              arraySlice1d< real64 > const & s,
                              arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & vt )
{
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > AT( a.size( 0 ), a.size( 1 ) );
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > UT( u.size( 0 ), u.size( 1 ) );
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > V( vt.size( 0 ), vt.size( 1 ) );

  // convert A to a column major format
  for( int i = 0; i < a.size( 0 ); ++i )
  {
    for( int j = 0; j < a.size( 1 ); ++j )
    {
      AT( i, j ) = a( i, j );
    }
  }

  matrixSVD( AT.toSliceConst(), UT.toSlice(), s, V.toSlice() );

  // convert U and VT back to row-major format
  for( int i = 0; i < u.size( 0 ); ++i )
  {
    for( int j = 0; j < u.size( 1 ); ++j )
    {
      u( i, j ) = UT( i, j );
    }
  }
  for( int i = 0; i < vt.size( 0 ); ++i )
  {
    for( int j = 0; j < vt.size( 1 ); ++j )
    {
      vt( i, j ) = V( i, j );
    }
  }
}

} // end geosx namespace
