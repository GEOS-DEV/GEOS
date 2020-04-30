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
 * @file DenseLA.cpp
 */

// Include the corresponding header file.
#include "DenseLA.hpp"

//#define USE_LAPACK

#include "DenseLAHelpers.hpp"
#include "BlasLapackFunctions.h"

// Put everything under the geosx namespace.
namespace geosx
{

real64 DenseLA::vectorNorm1( arraySlice1d< real64 const > const & X )
{
#ifdef USE_LAPACK

  real64 norm = 0;
  vectorNorm1Lapack( X, norm );
  return norm;

#else

  real64 norm = 0;
  for( localIndex i = 0; i < X.size(); ++i )
  {
    norm += std::fabs( X( i ) );
  }
  return norm;

#endif

}

real64 DenseLA::vectorNorm2( arraySlice1d< real64 const > const & X )
{
#ifdef USE_LAPACK

  real64 norm = 0;
  vectorNorm2Lapack( X, norm );
  return norm;

#else

  real64 norm = 0;
  for( localIndex i = 0; i < X.size(); ++i )
  {
    norm += X( i ) * X( i );
  }
  return sqrt( norm );

#endif
}

real64 DenseLA::vectorNormInf( arraySlice1d< real64 const > const & X )
{
#ifdef USE_LAPACK

  real64 norm = 0;
  vectorNormInfLapack( X, norm );
  return norm;

#else

  real64 norm = 0;
  real64 absEntry = 0;
  for( localIndex i = 0; i < X.size(); ++i )
  {
    absEntry = std::fabs( X( i ) );
    if( norm < absEntry )
    {
      norm = absEntry;
    }
  }
  return norm;

#endif
}


real64 DenseLA::determinant( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A )
{
  // --- check that matrix is square
  GEOSX_ASSERT_MSG( A.size( 0 ) == A.size( 1 ) &&
                    A.size( 0 ) > 0,
                    "Matrix must be square with order greater than zero" );

  real64 det;
  switch( A.size( 0 ) )
  {
    case 1:
    {
      det = A( 0, 0 );
      break;
    }
    case 2:
    {
      det = A( 0, 0 ) * A( 1, 1 ) - A( 0, 1 ) * A( 1, 0 );
      break;
    }
    case 3:
    {
      det = determinant3by3( A );
      break;
    }
    case 4:
    {
      det = determinant4by4( A );
      break;
    }
    default:
    {
      // Compute the determinant via LU factorization of A^t

#ifdef USE_LAPACK

      // TODO: put this in this helpers

      int const NN = integer_conversion< int >( A.size( 1 ) );
      int info;

      array1d< int > ipiv( NN );
      array2d< double, MatrixLayout::ROW_MAJOR_PERM > LUFactor( A.size( 0 ), A.size( 1 ) );
      matrixCopy( A, LUFactor );

      // Lapack will compute the LU factors for the transpose matrix, i.e. choosing the
      // LAPACK_COL_MAJOR ordering, to avoid transposition/copy required for LAPACK_ROW_MAJOR ordering.
      GEOSX_dgetrf( &NN, &NN, LUFactor.data(), &NN, ipiv.data(), &info );

      GEOSX_ASSERT_MSG( info == 0, "LAPACK dgetrf error code: " << info );

      // we need this offset because ipiv is based on Fortran convention (counting from 1)
      localIndex const offset = 1;

#else

      // 1. A is const, so copy its entries into a new matrix, LUFactor
      array1d< localIndex > ipiv( A.size( 1 ) );
      array2d< real64, MatrixLayout::ROW_MAJOR_PERM > LUFactor( A.size( 0 ), A.size( 1 ) );
      matrixCopy( A, LUFactor );

      // 2. Compute the LU decomposition of the matrix A^t, stored in LUFactor
      // TODO: see if this can be made faster
      matrixTFactorize( LUFactor, ipiv );

      // we don't need any offset when we do not use lapack
      localIndex const offset = 0;

#endif

      // 3. The LU factors are ready, we can compute the determinant
      det = determinantFromLUFactors( LUFactor, ipiv, offset );
      break;
    }
  }
  return det;
}

real64 DenseLA::determinant( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A )
{
  // --- check that matrix is square
  GEOSX_ASSERT_MSG( A.size( 0 ) == A.size( 1 ) &&
                    A.size( 0 ) > 0,
                    "Matrix must be square with order greater than zero" );

  real64 det;
  switch( A.size( 0 ) )
  {
    case 1:
    {
      det = A( 0, 0 );
      break;
    }
    case 2:
    {
      det = A( 0, 0 ) * A( 1, 1 ) - A( 0, 1 ) * A( 1, 0 );
      break;
    }
    case 3:
    {
      det = determinant3by3( A );
      break;
    }
    case 4:
    {
      det = determinant4by4( A );
      break;
    }
    default:
    {
      // Compute the determinant via LU factorization of A

#ifdef USE_LAPACK

      // TODO: put this is in the helpers

      int const NN = integer_conversion< int >( A.size( 0 ) );
      int info;

      array1d< int > ipiv( NN );
      array2d< double, MatrixLayout::COL_MAJOR_PERM > LUFactor( A.size( 0 ), A.size( 1 ) );
      matrixCopy( A, LUFactor );

      GEOSX_dgetrf( &NN, &NN, LUFactor.data(), &NN, ipiv.data(), &info );

      GEOSX_ASSERT_MSG( info == 0, "LAPACK dgetrf error code: " << info );

      // we need this offset because ipiv is based on Fortran convention (counting from 1)
      localIndex const offset = 1;

#else

      // 1. A is const, so copy its entries into a new matrix, LUFactor
      array1d< localIndex > ipiv( A.size( 0 ) );
      array2d< real64, MatrixLayout::COL_MAJOR_PERM > LUFactor( A.size( 0 ), A.size( 1 ) );
      matrixCopy( A, LUFactor );

      // 2. Compute the LU factors of the matrix A, stored in LUFactor
      // TODO: see if this can be made faster
      matrixFactorize( LUFactor, ipiv );

      // we don't need any offset when we do not use lapack
      localIndex const offset = 0;

#endif

      // 3. The LU factors are ready, we can compute the determinant
      det = determinantFromLUFactors( LUFactor, ipiv, offset );
      break;
    }
  }
  return det;
}

real64 DenseLA::matrixNormInf( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A )
{
#ifdef USE_LAPACK

  // For row-major, computed as one-norm of the transpose matrix
  real64 norm = 0;
  matrixNormLapack( A, '1', true, norm );
  return norm;

#else

  real64 norm = 0;
  real64 rowSum = 0;

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    rowSum = 0;
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      rowSum += std::fabs( A( i, j ) );
    }

    if( norm < rowSum )
    {
      norm = rowSum;
    }
  }
  return norm;

#endif
}

real64 DenseLA::matrixNormInf( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A )
{
#ifdef USE_LAPACK

  real64 norm = 0;
  matrixNormLapack( A, 'I', false, norm );
  return norm;

#else

  array1d< real64 > rowSum;
  rowSum.resize( A.size( 0 ) );
  rowSum = 0;

  for( localIndex j = 0; j < A.size( 1 ); ++j )
  {
    for( localIndex i = 0; i < A.size( 0 ); ++i )
    {
      rowSum( i ) += std::fabs( A( i, j ) );
    }
  }
  return *std::max_element( rowSum.begin(), rowSum.end());

#endif
}

real64 DenseLA::matrixNorm1( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A )
{
#ifdef USE_LAPACK

  // For row-major, computed as infinity-norm of the transpose matrix
  real64 norm = 0;
  matrixNormLapack( A, 'I', true, norm );
  return norm;

#else

  array1d< real64 > columnSum;
  columnSum.resize( A.size( 1 ) );
  columnSum = 0;

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      columnSum( j ) += std::fabs( A( i, j ) );
    }
  }
  return *std::max_element( columnSum.begin(), columnSum.end());

#endif
}

real64 DenseLA::matrixNorm1( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A )
{
#ifdef USE_LAPACK

  real64 norm = 0;
  matrixNormLapack( A, '1', false, norm );
  return norm;

#else

  real64 norm = 0;
  real64 columnSum = 0;

  for( localIndex j = 0; j < A.size( 1 ); ++j )
  {
    columnSum = 0;
    for( localIndex i = 0; i < A.size( 0 ); ++i )
    {
      columnSum += std::fabs( A( i, j ) );
    }

    if( norm < columnSum )
    {
      norm = columnSum;
    }
  }
  return norm;

#endif
}


real64 DenseLA::matrixNormFrobenius( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A )
{
#ifdef USE_LAPACK

  // For row-major, computed using the transpose matrix
  real64 norm = 0;
  matrixNormLapack( A, 'F', true, norm );
  return norm;

#else

  real64 norm = 0;

  for( localIndex j = 0; j < A.size( 1 ); ++j )
  {
    for( localIndex i = 0; i < A.size( 0 ); ++i )
    {
      norm += A( i, j ) * A( i, j );
    }
  }
  return sqrt( norm );

#endif
}

real64 DenseLA::matrixNormFrobenius( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A )
{
#ifdef USE_LAPACK

  real64 norm = 0;
  matrixNormLapack( A, 'F', false, norm );
  return norm;

#else

  real64 norm = 0;

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      norm += A( i, j ) * A( i, j );
    }
  }
  return sqrt( norm );

#endif
}

void DenseLA::vectorVectorAdd( arraySlice1d< real64 const > const & X,
                               arraySlice1d< real64 > const & Y,
                               real64 const alpha )
{
  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                    "Vector dimensions not compatible for sum" );

#ifdef USE_LAPACK

  vectorVectorAddLapack( X, Y, alpha );

#else

  for( localIndex i = 0; i < X.size(); ++i )
  {
    Y( i ) += alpha * X( i );
  }

#endif
}

void DenseLA::matrixMatrixAdd( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                               arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & B,
                               real64 const alpha )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == B.size( 0 ) &&
                    A.size( 1 ) == B.size( 1 ),
                    "Matrix dimensions not compatible for sum" );

#ifdef USE_LAPACK

  matrixMatrixAddLapack( A, B, alpha );

#else

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      B( i, j ) += alpha * A( i, j );
    }
  }

#endif
}

void DenseLA::matrixMatrixAdd( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                               arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & B,
                               real64 const alpha )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == B.size( 0 ) &&
                    A.size( 1 ) == B.size( 1 ),
                    "Matrix dimensions not compatible for sum" );

#ifdef USE_LAPACK

  matrixMatrixAddLapack( A, B, alpha );

#else

  for( localIndex j = 0; j < A.size( 1 ); ++j )
  {
    for( localIndex i = 0; i < A.size( 0 ); ++i )
    {
      B( i, j ) += alpha * A( i, j );
    }
  }

#endif
}

void DenseLA::vectorScale( real64 const & alpha,
                           arraySlice1d< real64 > const & X )
{
#ifdef USE_LAPACK

  vectorScaleLapack( alpha, X );

#else

  for( localIndex i = 0; i < X.size(); ++i )
  {
    X( i ) *= alpha;
  }

#endif
}

void DenseLA::matrixScale( real64 const & alpha,
                           arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & A )
{
#ifdef USE_LAPACK

  matrixScaleLapack( alpha, A );

#else

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      A( i, j ) *= alpha;
    }
  }

#endif
}

void DenseLA::matrixScale( real64 const & alpha,
                           arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & A )
{
#ifdef USE_LAPACK

  matrixScaleLapack( alpha, A );

#else

  for( localIndex j = 0; j < A.size( 1 ); ++j )
  {
    for( localIndex i = 0; i < A.size( 0 ); ++i )
    {
      A( i, j ) *= alpha;
    }
  }

#endif
}

real64 DenseLA::vectorDot( arraySlice1d< real64 const > const & X,
                           arraySlice1d< real64 const > const & Y )
{
  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                    "Vector dimensions not compatible for dot product" );

#ifdef USE_LAPACK

  real64 dotProduct = 0;
  vectorDotLapack( X, Y, dotProduct );
  return dotProduct;

#else

  real64 dotProduct = 0;

  for( localIndex i = 0; i < X.size(); ++i )
  {
    dotProduct += X( i ) * Y( i );
  }
  return dotProduct;

#endif
}

void DenseLA::matrixVectorMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                    arraySlice1d< real64 const > const & X,
                                    arraySlice1d< real64 > const & Y,
                                    real64 const alpha,
                                    real64 const beta )
{
  GEOSX_ASSERT_MSG( A.size( 1 ) == X.size() &&
                    A.size( 0 ) == Y.size(),
                    "Matrix, source vector and destination vector not compatible" );

#ifdef USE_LAPACK

  // TODO: move to the helper file

  int const M = integer_conversion< int >( A.size( 0 ) );
  int const N = 1;
  int const K = integer_conversion< int >( A.size( 1 ) );

  // A*X = Y is computed as X^T * A^T = Y^T, i.e. accessing the transpose
  // matrix using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, X.dataIfContiguous(), &N, A.dataIfContiguous(), &K, &beta, Y.dataIfContiguous(), &N );

#else

  vectorScale( beta, Y );

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    real64 sum_A_ij_X_j = 0;
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      sum_A_ij_X_j += A( i, j ) * X( j );
    }
    Y( i ) += alpha * sum_A_ij_X_j;
  }

#endif
}

void DenseLA::matrixTVectorMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                     arraySlice1d< real64 const > const & X,
                                     arraySlice1d< real64 > const & Y,
                                     real64 const alpha,
                                     real64 const beta )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == X.size() &&
                    A.size( 1 ) == Y.size(),
                    "Matrix, source vector and destination vector not compatible" );

#ifdef USE_LAPACK

  // TODO: move to the helper file

  int const M = integer_conversion< int >( A.size( 1 ) );
  int const N = 1;
  int const K = integer_conversion< int >( A.size( 0 ) );

  // A^T*X = Y is computed as X^T * A = Y^T, i.e. accessing the transpose
  // matrix using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'T';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, X.dataIfContiguous(), &N, A.dataIfContiguous(), &M, &beta, Y.dataIfContiguous(), &N );

#else

  vectorScale( beta, Y );

  for( localIndex j = 0; j < A.size( 0 ); ++j )
  {
    real64 const alpha_X_j = alpha * X( j );
    for( localIndex i = 0; i < A.size( 1 ); ++i )
    {
      Y( i ) += alpha_X_j * A( j, i );
    }
  }

#endif
}

void DenseLA::matrixMatrixMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                                    arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                                    real64 const alpha,
                                    real64 const beta )
{
  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                    C.size( 1 ) == B.size( 1 ) &&
                    A.size( 1 ) == B.size( 0 ),
                    "Matrix dimensions not compatible for product" );

#ifdef USE_LAPACK

  // TODO: move to the helper file

  int const M = integer_conversion< int >( A.size( 0 ) );
  int const N = integer_conversion< int >( B.size( 1 ) );
  int const K = integer_conversion< int >( A.size( 1 ) );

  // A*B = C is computed as B^T * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, B.dataIfContiguous(), &N, A.dataIfContiguous(), &K, &beta, C.dataIfContiguous(), &N );


#else

  // TODO: double-check loop order, see if the loops need to be interchanged
  // TODO: see if using the raw pointer instead of the accessors is more efficient

  // compute beta * C
  matrixScale( beta, C );

  // add alpha * A * B
  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex k = 0; k < A.size( 1 ); ++k )
    {
      real64 const alpha_A_ik = alpha * A( i, k );
      for( localIndex j = 0; j < B.size( 1 ); ++j )
      {
        C( i, j ) +=  alpha_A_ik * B( k, j );
      }
    }
  }

#endif
}

void DenseLA::matrixTMatrixMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                     arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                                     arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                                     real64 const alpha,
                                     real64 const beta )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 1 ) &&
                    C.size( 1 ) == B.size( 1 ) &&
                    A.size( 0 ) == B.size( 0 ),
                    "Matrix dimensions not compatible for product" );

#ifdef USE_LAPACK

  // TODO: move to the helper file

  int const M = integer_conversion< int >( A.size( 1 ) );
  int const N = integer_conversion< int >( B.size( 1 ) );
  int const K = integer_conversion< int >( A.size( 0 ) );

  // A^T*B = C is computed as B^T * A = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'N';
  char const TRANS2 = 'T';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, B.dataIfContiguous(), &N, A.dataIfContiguous(), &M, &beta, C.dataIfContiguous(), &N );

#else

  // TODO: double-check loop order, see if the loops need to be interchanged
  // TODO: see if using the raw pointer instead of the accessors is more efficient

  // compute beta * C
  matrixScale( beta, C );

  // add alpha * A * B
  for( localIndex k = 0; k < A.size( 0 ); ++k )
  {
    for( localIndex i= 0; i < A.size( 1 ); ++i )
    {
      real64 const alpha_A_ki = alpha * A( k, i );
      for( localIndex j = 0; j < B.size( 1 ); ++j )
      {
        C( i, j ) += alpha_A_ki * B( k, j );
      }
    }
  }

#endif
}

void DenseLA::matrixMatrixTMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                     arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                                     arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                                     real64 const alpha,
                                     real64 const beta )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                    C.size( 1 ) == B.size( 0 ) &&
                    A.size( 1 ) == B.size( 1 ),
                    "Matrix dimensions not compatible for product" );

#ifdef USE_LAPACK

  // TODO: move to the helper file

  int const M = integer_conversion< int >( A.size( 0 ) );
  int const N = integer_conversion< int >( B.size( 0 ) );
  int const K = integer_conversion< int >( A.size( 1 ) );

  // A*B^T = C is computed as B * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'T';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, B.dataIfContiguous(), &K, A.dataIfContiguous(), &K, &beta, C.dataIfContiguous(), &N );

#else

  // TODO: double-check loop order, see if the loops need to be interchanged
  // TODO: see if using the raw pointer instead of the accessors is more efficient

  // compute beta * C
  matrixScale( beta, C );

  // add alpha * A * B
  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < B.size( 0 ); ++j )
    {
      real64 sum_a_ik_b_jk = 0;
      for( localIndex k = 0; k < A.size( 1 ); ++k )
      {
        sum_a_ik_b_jk += A( i, k ) * B( j, k );
      }
      C( i, j ) += alpha * sum_a_ik_b_jk;
    }
  }

#endif
}

void DenseLA::matrixTMatrixTMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                      arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                                      arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                                      real64 const alpha,
                                      real64 const beta )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 1 ) &&
                    C.size( 1 ) == B.size( 0 ) &&
                    A.size( 0 ) == B.size( 1 ),
                    "Matrix dimensions not compatible for product" );

#ifdef USE_LAPACK

  int const M = integer_conversion< int >( A.size( 1 ) );
  int const N = integer_conversion< int >( B.size( 0 ) );
  int const K = integer_conversion< int >( A.size( 0 ) );

  // A^T*B^T = C is computed as B * A = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'T';
  char const TRANS2 = 'T';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, B.dataIfContiguous(), &K, A.dataIfContiguous(), &M, &beta, C.dataIfContiguous(), &N );

#else

  // TODO: what could be a good loop ordering for this one??

  // compute beta * C
  matrixScale( beta, C );

  // add alpha * A * B
  for( localIndex i = 0; i < A.size( 1 ); ++i )
  {
    for( localIndex j = 0; j < B.size( 0 ); ++j )
    {
      real64 sum_a_ki_b_jk = 0;
      for( localIndex k = 0; k < A.size( 0 ); ++k )
      {
        sum_a_ki_b_jk += A( k, i ) * B( j, k );
      }
      C( i, j ) += alpha * sum_a_ki_b_jk;
    }
  }

#endif
}

void DenseLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                             arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & Ainv,
                             real64 & detA )
{
  // --- Check that source matrix is square
  GEOSX_ASSERT_MSG( A.size( 0 ) > 0 &&
                    A.size( 0 )  == A.size( 1 ),
                    "Matrix must be square" );

  // --- Check that inverse matrix has appropriate dimension
  GEOSX_ASSERT_MSG( Ainv.size( 0 ) == A.size( 0 ) &&
                    Ainv.size( 1 ) == A.size( 0 ),
                    "Inverse matrix has wrong dimensions" );

#ifdef USE_LAPACK

  int const NN = integer_conversion< int >( A.size( 0 ) );
  int info;

  array1d< int > ipiv;
  array1d< double > invWork;

#else

  array1d< localIndex > ipiv;
  array2d< real64, MatrixLayout::ROW_MAJOR_PERM > LUFactor;

#endif

  // --- Check if matrix is singular by computing the determinant
  //     note: if order greater than 3 we compute the determinant by
  //           first constructing the LU factors, later reused for calculating
  //           the inverse.

  if( A.size( 0 ) <= 3 )
  {
    detA = determinant( A.toSliceConst());
  }
  else
  {
#ifdef USE_LAPACK

    // 1. Copy A in Ainv
    ipiv.resize( A.size( 0 ) );
    matrixCopy( A, Ainv );

    // 2.  Compute determinant (not done calling directly the function determinant
    // (avoid computing twice LUFactors, currently stored in Ainv, needed for
    // computing the inverse). We compute the LU factors for the
    // transpose matrix, i.e. choosing the LAPACK_COL_MAJOR ordering, to
    // avoid transposition/copy requires for LAPACK_ROW_MAJOR ordering.
    GEOSX_dgetrf( &NN, &NN, Ainv.dataIfContiguous(), &NN, ipiv.data(), &info );

    GEOSX_ASSERT_MSG( info == 0, "LAPACK dgetrf error code: " << info );

    // 3. The LU factors are ready, we can compute the determinant
    localIndex const offset = 1;
    detA = determinantFromLUFactors( Ainv, ipiv, offset );

#else

    // 1. Copy A into LUFactor
    ipiv.resize( A.size( 0 ) );
    LUFactor.resize( A.size( 0 ), A.size( 1 ) );
    matrixCopy( A, LUFactor );

    // 2. Compute the LU factors of the matrix A^t
    matrixTFactorize( LUFactor, ipiv );

    // 3. The LU factors are ready, we can compute the determinant
    localIndex const offset = 0;
    detA = determinantFromLUFactors( LUFactor, ipiv, offset );

#endif

  }

  // Check if matrix is singular
  GEOSX_ASSERT_MSG( std::fabs( detA ) >
                    std::numeric_limits< real64 >::epsilon() *
                    matrixNormFrobenius( A ),
                    "Matrix is singular" );
  real64 const oneOverDetA = 1. / detA;

  // --- Compute inverse
  switch( A.size( 0 ) )
  {
    case 1:
    {
      Ainv( 0, 0 ) = oneOverDetA;
      break;
    }

    case 2:
    {
      matrixInverse2by2( A, Ainv, oneOverDetA );
      break;
    }

    case 3:
    {
      matrixInverse3by3( A, Ainv, oneOverDetA );
      break;
    }
    default:
    {

#ifdef USE_LAPACK

      invWork.resize( A.size( 0 ) );

      // Invert (LAPACK function DGETRI). The LU factors computed for the
      // transpose matrix stored in Ainv are used.
      GEOSX_dgetri( &NN, Ainv.dataIfContiguous(), &NN, ipiv.data(), invWork.data(), &NN, &info );

      GEOSX_ASSERT_MSG( info == 0, "LAPACK dgetri error code: " << info );

#else

      matrixTInverseFromLUFactors( LUFactor, ipiv, Ainv );

#endif

      break;
    }
  }
}

void DenseLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                             arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & Ainv,
                             real64 & detA )
{
  // --- Check that source matrix is square
  GEOSX_ASSERT_MSG( A.size( 0 ) > 0 &&
                    A.size( 0 )  == A.size( 1 ),
                    "Matrix must be square" );

  // --- Check that inverse matrix has appropriate dimension
  GEOSX_ASSERT_MSG( Ainv.size( 0 ) == A.size( 0 ) &&
                    Ainv.size( 1 ) == A.size( 0 ),
                    "Inverse matrix has wrong dimensions" );

#ifdef USE_LAPACK

  int const NN = integer_conversion< int >( A.size( 0 ) );
  int info;

  array1d< int > ipiv;
  array1d< double > invWork;

#else

  array1d< localIndex > ipiv;
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > LUFactor;

#endif

  // --- Check if matrix is singular by computing the determinant
  //     note: if order greater than 3 we compute the determinant by
  //           first constructing the LU factors, later reused for calculating
  //           the inverse.

  if( A.size( 0 ) <= 3 )
  {
    detA = determinant( A.toSliceConst());
  }
  else
  {

#ifdef USE_LAPACK

    // 1. Copy A in Ainv
    ipiv.resize( A.size( 0 ) );
    matrixCopy( A, Ainv );

    // 2.  Compute determinant (not done calling directly the function determinant
    // (avoid computing twice LUFactors, currently stored in Ainv, needed for
    // computing the inverse).
    GEOSX_dgetrf( &NN, &NN, Ainv.dataIfContiguous(), &NN, ipiv.data(), &info );

    GEOSX_ASSERT_MSG( info == 0, "LAPACK dgetrf error code: " << info );

    // 3. The LU factors are ready, we can compute the determinant
    localIndex const offset = 1;
    detA = determinantFromLUFactors( Ainv, ipiv, offset );

#else

    // 1. Copy A into LUFactor
    ipiv.resize( A.size( 0 ) );
    LUFactor.resize( A.size( 0 ), A.size( 1 ) );
    matrixCopy( A, LUFactor );

    // 2. Compute the LU factors of the matrix A, stored in LUFactor
    matrixFactorize( LUFactor, ipiv );

    // 3. The LU factors are ready, we can compute the determinant
    localIndex const offset = 0;
    detA = determinantFromLUFactors( LUFactor, ipiv, offset );

#endif
  }

  // Check if matrix is singular
  GEOSX_ASSERT_MSG( std::fabs( detA ) >
                    std::numeric_limits< real64 >::epsilon() *
                    matrixNormFrobenius( A ),
                    "Matrix is singular" );
  real64 const oneOverDetA = 1. / detA;

  // --- Compute inverse
  switch( A.size( 0 ) )
  {
    case 1:
    {
      Ainv( 0, 0 ) = oneOverDetA;
      break;
    }

    case 2:
    {
      matrixInverse2by2( A, Ainv, oneOverDetA );
      break;
    }

    case 3:
    {
      matrixInverse3by3( A, Ainv, oneOverDetA );
      break;
    }
    default:
    {

#ifdef USE_LAPACK

      invWork.resize( A.size( 0 ) );

      // Invert (LAPACK function DGETRI). The LU factors stored in Ainv are used.
      GEOSX_dgetri( &NN, Ainv.dataIfContiguous(), &NN, ipiv.data(), invWork.data(), &NN, &info );

      GEOSX_ASSERT_MSG( info == 0, "LAPACK dgetri error code: " << info );

#else

      matrixInverseFromLUFactors( LUFactor, ipiv, Ainv );

#endif

      break;
    }
  }
}

void DenseLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                             arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & Ainv )
{
  real64 detA = 0;
  matrixInverse( A, Ainv, detA );
}

void DenseLA::matrixInverse( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                             arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & Ainv )
{
  real64 detA = 0;
  matrixInverse( A, Ainv, detA );
}

void DenseLA::vectorCopy( array1d< real64 > const & X,
                          array1d< real64 > & Y )
{
  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                    "Vector dimensions not compatible for copying" );

#ifdef USE_LAPACK

  vectorCopyLapack( X, Y );

#else

  for( localIndex i = 0; i < X.size(); ++i )
  {
    Y( i ) = X( i );
  }

#endif
}

void DenseLA::matrixCopy( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                          arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & B )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == B.size( 0 ) &&
                    A.size( 1 ) == B.size( 1 ),
                    "Matrix dimensions not compatible for copying" );

#ifdef USE_LAPACK

  matrixCopyLapack( A, B );

#else

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      B( i, j ) = A( i, j );
    }
  }

#endif
}

void DenseLA::matrixCopy( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                          arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & B )
{
  GEOSX_ASSERT_MSG( A.size( 0 ) == B.size( 0 ) &&
                    A.size( 1 ) == B.size( 1 ),
                    "Matrix dimensions not compatible for copying" );

#ifdef USE_LAPACK

  matrixCopyLapack( A, B );

#else

  for( localIndex j = 0; j < A.size( 1 ); ++j )
  {
    for( localIndex i = 0; i < A.size( 0 ); ++i )
    {
      B( i, j ) = A( i, j );
    }
  }

#endif
}

// below, this is work in progress

void DenseLA::setRandomNumberGeneratorSeed( arraySlice1d< int const > const & seed )
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

void DenseLA::getRandomNumberGeneratorSeed( arraySlice1d< int > const & seed )
{
  // Error checking
  GEOSX_ASSERT_MSG( seed.size() >= 4, "Seed array must have size at least four" );
  for( int i = 0; i < 4; ++i )
  {
    seed[i] = ISEED[i];
  }
}

void DenseLA::vectorRand( arraySlice1d< real64 > const & X,
                          RandomNumberDistribution const & idist )
{

  int IDIST = static_cast< int >(idist);
  int const N = integer_conversion< int >( X.size() );
  GEOSX_ASSERT_MSG( N > 0, "The vector cannot be empty" );
  GEOSX_dlarnv( &IDIST, ISEED, &N, X.dataIfContiguous());
}

void DenseLA::matrixRand( arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & A,
                          RandomNumberDistribution const & idist )
{
  int const IDIST = static_cast< int >(idist);
  int const NN = integer_conversion< int >( A.size() );
  GEOSX_ASSERT_MSG( NN > 0, "The matrix cannot be empty" );
  GEOSX_dlarnv( &IDIST, ISEED, &NN, A.dataIfContiguous() );
}

void DenseLA::matrixRand( arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & A,
                          RandomNumberDistribution const & idist )
{
  int const IDIST = static_cast< int >(idist);
  int const NN = integer_conversion< int >( A.size() );
  GEOSX_ASSERT_MSG( NN > 0, "The matrix cannot be empty" );
  GEOSX_dlarnv( &IDIST, ISEED, &NN, A.dataIfContiguous() );
}

void DenseLA::matrixSVD( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                         arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & U,
                         arraySlice1d< real64 > const & S,
                         arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & VT )
{
  int const minDim = (A.size( 0 ) < A.size( 1 ))
                   ? integer_conversion< int >( A.size( 0 ) )
                   : integer_conversion< int >( A.size( 1 ) );

  GEOSX_ASSERT_MSG( A.size( 0 ) == U.size( 0 ) && minDim == U.size( 1 ),
                    "The matrices A and U have an incompatible size" );

  GEOSX_ASSERT_MSG( minDim == VT.size( 0 ) && A.size( 1 ) == VT.size( 1 ),
                    "The matrices A and V have an incompatible size" );

  GEOSX_ASSERT_MSG( S.size() == minDim,
                    "The matrix A and vector S have an incompatible size" );

  // make a copy of A, since dgesvd destroys contents
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > Acopy( A.size( 0 ), A.size( 1 ) );
  matrixCopy( A, Acopy );

  // define the arguments of dgesvd
  int const M     = integer_conversion< int >( A.size( 0 ) );
  int const N     = integer_conversion< int >( A.size( 1 ) );
  int const LDA   = M;
  int const LDU   = M;
  int const LDVT  = minDim;
  int lwork = 0;
  int info  = 0;
  double wkopt = 0.0;

  // 1) query and allocate the optimal workspace
  lwork = -1;
  GEOSX_dgesvd( "S", "S",
                &M, &N, Acopy.data(), &LDA,
                S.dataIfContiguous(), U.dataIfContiguous(), &LDU, VT.dataIfContiguous(), &LDVT,
                &wkopt, &lwork, &info );

  lwork = static_cast< int >( wkopt );
  array1d< real64 > work( lwork );

  // 2) compute svd
  GEOSX_dgesvd( "S", "S",
                &M, &N, Acopy.data(), &LDA,
                S.dataIfContiguous(), U.dataIfContiguous(), &LDU, VT.dataIfContiguous(), &LDVT,
                work.data(), &lwork, &info );

  GEOSX_ASSERT_MSG( info == 0, "The algorithm computing SVD failed to converge." );
}

void DenseLA::matrixSVD( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                         arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & U,
                         arraySlice1d< real64 > const & S,
                         arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & VT )
{
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > AT( A.size( 0 ), A.size( 1 ) );
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > UT( U.size( 0 ), U.size( 1 ) );
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > V( VT.size( 0 ), VT.size( 1 ) );

  // convert A to a column major format
  for( int i = 0; i < A.size( 0 ); ++i )
  {
    for( int j = 0; j < A.size( 1 ); ++j )
    {
      AT( i, j ) = A( i, j );
    }
  }

  matrixSVD( AT.toSliceConst(), UT.toSlice(), S, V.toSlice() );

  // convert U and VT back to row-major format
  for( int i = 0; i < U.size( 0 ); ++i )
  {
    for( int j = 0; j < U.size( 1 ); ++j )
    {
      U( i, j ) = UT( i, j );
    }
  }
  for( int i = 0; i < VT.size( 0 ); ++i )
  {
    for( int j = 0; j < VT.size( 1 ); ++j )
    {
      VT( i, j ) = V( i, j );
    }
  }
}

} // end geosx namespace
