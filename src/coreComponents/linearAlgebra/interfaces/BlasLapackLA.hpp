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
 * @file BlasLapackLA.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKLA_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKLA_HPP_

#include "common/DataTypes.hpp"
#include "common/Logger.hpp"
#include "linearAlgebra/common.hpp"

// BLAS and LAPACK function declaration
#include "BlasLapackFunctions.h"

namespace geosx
{

namespace BlasLapackLA
{

/**
 * \enum  RandomNumberDistribution
 * \brief This enum class specifies the type of distribution for
 *        generating random real numbers.
 */
enum class RandomNumberDistribution : int
{
  UNIFORM_01 = 1,     /**< uniform distribution (0,1); */
  UNIFORM_m1p1 = 2,   /**< uniform distribution (-1,1); */
  NORMAL_01 = 3       /**< normal distribution (0,1); */
};

/**
 * @brief Returns the 1-norm of the vector.
 *
 * @param [in] X GEOSX array1d.
 */
real64 vectorNorm1( arraySlice1d< real64 const > const & X );

/**
 * @brief Returns the two norm of the vector.
 *
 * @param [in] X GEOSX array1d.
 */

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static real64 vectorNorm2( arraySlice1d< real64 const > const & X )
{
#if !defined(__CUDA_ARCH__)

  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( X.size() );
  return GEOSX_dnrm2( &N, X.dataIfContiguous(), &INCX );

#else

  // TODO: template on dimensions, move to LvArray

  real64 norm = 0;
  for( localIndex i = 0; i < X.size(); ++i )
  {
    norm = norm + X( i ) * X( i );
  }
  return sqrt( norm );

#endif
}


/**
 * @brief Returns the infinity-norm of the vector.
 *
 * @param [in] X GEOSX array1d.
 */
real64 vectorNormInf( arraySlice1d< real64 const > const & X );

/**
 * @brief Returns the determinant of a square matrix.
 *
 * @param [in] A GEOSX array2d.
 *
 * @note
 * This function is hardcoded for square matrices up to order four.
 * For dimensions larger than four, the determinant is computed using
 * LAPACK's function DGETRF. To avoid matrix transposition/copy
 * due to the row major ordering used in GEOSX for array2d, the determinant
 * is computed for the transpose matrix, i.e. assuming column major
 * ordering, for best performance.
 */
real64 determinant( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A );

/**
 * @copydoc determinant( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const & )
 */
real64 determinant( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A );

/**
 * @brief Returns the infinity norm of the matrix.
 *
 * @param [in] A GEOSX array2d.
 *
 * @note
 * Row major ordering is used for GEOSX array2d. Since LAPACK native
 * routines are using a column major ordering (Fortran), the infinity
 * norm is computed as the one norm of the transpose matrix, i.e. assuming
 * column major ordering, for best performance.
 */
real64 matrixNormInf( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A );

/**
 * @copydoc matrixNormInf( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const & )
 */
real64 matrixNormInf( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A );

/**
 * @brief Returns the one norm of the matrix.
 *
 * @param [in] A GEOSX array2d.
 *
 * @note
 * Row major ordering is used for GEOSX array2d. Since LAPACK native
 * routines are using a column major ordering (Fortran), the one norm
 * is computed as the infinity norm of the transpose matrix, i.e. assuming
 * column major ordering, for best performance.
 */
real64 matrixNorm1( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A );

/**
 * @copydoc matrixNorm1( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const & )
 */
real64 matrixNorm1( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A );

/**
 * @brief Returns the Frobenius norm of the matrix.
 *
 * @param [in] A GEOSX array2d.
 *
 * @note
 * Row major ordering is used for GEOSX array2d. Since LAPACK native
 * routines are using a column major ordering (Fortran), the one norm
 * is computed for the transpose matrix, i.e. assuming column major
 * ordering, for best performance.
 */
real64 matrixNormFrobenius( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A );

/**
 * @copydoc matrixNormFrobenius( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const & )
 */
real64 matrixNormFrobenius( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A );

/**
 * @brief Vector-Vector sum;
 * \p Y = \p alpha*\p X + \p Y.
 *
 * Computes (\p alpha*\p X + \p Y) and overwrites the result on
 * \p Y, with optional scaling.
 *
 * @param [in]     X     GEOSX array1d.
 * @param [in,out] Y     GEOSX array1d.
 * @param [in]     alpha Optional scalar to multiply with \p X.
 *
 * @warning
 * Assumes that \p X and \p Y have the same size.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static void vectorVectorAdd( arraySlice1d< real64 const > const & X,
                             arraySlice1d< real64 > const & Y,
                             real64 const alpha = 1. )
{
  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                    "Vector dimensions not compatible for sum" );

#if !defined(__CUDA_ARCH__)

  int const INCX = 1;
  int const INCY = 1;
  int const N = LvArray::integerConversion< int >( X.size() );
  GEOSX_daxpy( &N, &alpha, X.dataIfContiguous(), &INCX, Y.dataIfContiguous(), &INCY );

#else

  // TODO: template on dimensions, move to LvArray

  for( localIndex i = 0; i < X.size(); ++i )
  {
    Y( i ) = Y( i ) + alpha * X( i );
  }

#endif
}


/**
 * @brief Matrix-Matrix sum;
 * \p B = \p alpha*\p A + \p B.
 *
 * Computes (\p alpha*\p A + \p B) and overwrites the result on p B, with
 * optional scaling.
 *
 * @param [in]     A     GEOSX array2d.
 * @param [in,out] B     GEOSX array2d.
 * @param [in]     alpha Optional scalar to multiply with \p A.
 *
 * @warning
 * Assumes that \p A and \p B have the same size.
 */
void matrixMatrixAdd( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                      arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & B,
                      real64 const alpha = 1. );

/**
 * @copydoc matrixMatrixAdd( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const & A,
                             arraySlice2d<real64, MatrixLayout::ROW_MAJOR> const & B,
                             real64 const alpha = 1. )
 */
void matrixMatrixAdd( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                      arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & B,
                      real64 const alpha = 1. );

/**
 * @brief In-place scalar-vector product;
 * \p X = \p alpha*\p X.
 *
 * @param [in]     alpha Scalar to multiply with \p X.
 * @param [in,out] X     GEOSX array1d.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static void vectorScale( real64 const & alpha,
                         arraySlice1d< real64 > const & X )
{
#if !defined(__CUDA_ARCH__)

  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( X.size() );
  GEOSX_dscal( &N, &alpha, X.dataIfContiguous(), &INCX );

#else

  // TODO: template on dimensions, move to LvArray

  for( localIndex i = 0; i < X.size(); ++i )
  {
    X( i ) = X( i ) * alpha;
  }

#endif
}


/**
 * @brief In-place scalar-matrix product;
 * \p A = \p alpha*\p A.
 *
 * @param [in]     alpha Scalar to multiply with <tt>A</tt>.
 * @param [in,out] A     GEOSX array2d.
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static void matrixScale( real64 const & alpha,
                         arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & A )
{
#if !defined(__CUDA_ARCH__)

  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( A.size() );
  GEOSX_dscal( &N, &alpha, A.dataIfContiguous(), &INCX );

#else

  // TODO: template on dimensions, move to LvArray

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < A.size( 1 ); ++j )
    {
      A( i, j ) = A( i, j ) * alpha;
    }
  }

#endif
}


/**
 * @copydoc static void matrixScale( real64 const alpha, arraySlice2d<real64, MatrixLayout::ROW_MAJOR> const & )
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static void matrixScale( real64 const & alpha,
                         arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & A )
{
#if !defined(__CUDA_ARCH__)

  int const INCX = 1;
  int const N = LvArray::integerConversion< int >( A.size() );
  GEOSX_dscal( &N, &alpha, A.dataIfContiguous(), &INCX );

#else

  // TODO: template on dimensions, move to LvArray

  for( localIndex j = 0; j < A.size( 1 ); ++j )
  {
    for( localIndex i = 0; i < A.size( 0 ); ++i )
    {
      A( i, j ) = A( i, j ) * alpha;
    }
  }

#endif
}
  

/**
 * @brief Returns the dot product of two vectors.
 *
 * @param [in] X GEOSX array1d.
 * @param [in] Y GEOSX array1d.
 *
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static real64 vectorDot( arraySlice1d< real64 const > const & X,
                         arraySlice1d< real64 const > const & Y )
{
  GEOSX_ASSERT_MSG( X.size() == Y.size(),
                    "Vector dimensions not compatible for dot product" );

#if !defined(__CUDA_ARCH__)

  int const INCX = 1;
  int const INCY = 1;
  int const N = LvArray::integerConversion< int >( X.size() );
  return GEOSX_ddot( &N, X.dataIfContiguous(), &INCX, Y.dataIfContiguous(), &INCY );

#else

  // TODO: template on dimensions, move to LvArray

  real64 dotProduct = 0;

  for( localIndex i = 0; i < X.size(); ++i )
  {
    dotProduct = dotProduct + X( i ) * Y( i );
  }
  return dotProduct;

#endif
}

/**
 * @brief Matrix-Vector product;
 * \p Y = \p alpha*\p A*\p X + \p beta*\p Y.
 *
 * Computes matrix-vector product with optional scaling and accumulation.
 *
 * @param [in]     A     GEOSX array2d.
 * @param [in]     X     GEOSX array1d.
 * @param [in,out] Y     GEOSX array1d.
 * @param [in]     alpha Optional scalar to multiply with \p A*\p X.
 * @param [in]     beta  Optional parameter to control the accumulation.
 *
 * @warning
 * Assumes that \p X and \p Y have compatible sizes with \p A.
 */
void matrixVectorMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                           arraySlice1d< real64 const > const & X,
                           arraySlice1d< real64 > const & Y,
                           real64 const alpha = 1.0,
                           real64 const beta = 0.0 );

/**
 * @brief transpose(Matrix)-Vector product;
 * \p Y = \p alpha*\p transpose(A)*\p X + \p beta*\p Y.
 *
 * Computes transpose(matrix)-vector product with optional scaling and
 * accumulation.
 *
 * @param [in]     A     GEOSX array2d.
 * @param [in]     X     GEOSX array1d.
 * @param [in,out] Y     GEOSX array1d.
 * @param [in]     alpha Optional scalar to multiply with
 *                       \p transpose(A)*\p X.
 * @param [in]     beta  Optional parameter to control the accumulation.
 *
 * @warning
 * Assumes that \p X and \p Y have compatible sizes with \p transpose(A).
 */
void matrixTVectorMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                            arraySlice1d< real64 const > const & X,
                            arraySlice1d< real64 > const & Y,
                            real64 const alpha = 1.0,
                            real64 const beta = 0.0 );

/**
 * @brief Matrix-Matrix product;
 * \p C = \p alpha*\p A*\p B + \p beta*\p C.
 *
 * Computes matrix-matrix product with optional scaling and accumulation.
 *
 * @param [in]     A     GEOSX array2d.
 * @param [in]     B     GEOSX array2d.
 * @param [in,out] C     GEOSX array1d.
 * @param [in]     alpha Optional scalar to multiply with \p A*\p B.
 * @param [in]     beta  Optional parameter to control the accumulation.
 *
 * @warning
 * Assumes that \p A and \p B have compatible sizes and that \p C already
 * has the right size.
 *
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static void matrixMatrixMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                  arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                                  arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                                  real64 const alpha = 1.0,
                                  real64 const beta = 0.0 )
{
  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                    C.size( 1 ) == B.size( 1 ) &&
                    A.size( 1 ) == B.size( 0 ),
                    "Matrix dimensions not compatible for product" );

#if !defined(__CUDA_ARCH__)

  // TODO: move to the helper file

  int const M = LvArray::integerConversion< int >( A.size( 0 ) );
  int const N = LvArray::integerConversion< int >( B.size( 1 ) );
  int const K = LvArray::integerConversion< int >( A.size( 1 ) );

  // A*B = C is computed as B^T * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout
  char const TRANS1 = 'N';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, B.dataIfContiguous(), &N, A.dataIfContiguous(), &K, &beta, C.dataIfContiguous(), &N );


#else

  // TODO: template on dimensions, move to LvArray

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < B.size( 1 ); ++j )
    {
      real64 sum_A_ik_B_kj = 0;
      for( localIndex k = 0; k < A.size( 1 ); ++k )
      {
        sum_A_ik_B_kj = sum_A_ik_B_kj + A( i, k ) * B( k, j );
      }
      C( i, j ) = beta * C( i, j ) + alpha * sum_A_ik_B_kj;
    }
  }

#endif
}


/**
 * @brief transpose(Matrix)-Matrix product;
 * \p C = \p alpha*\p transpose(A)*\p B + \p beta*\p C.
 *
 * Computes transpose(matrix)-matrix product with optional scaling
 * and accumulation.
 *
 * @param [in]     A     GEOSX array2d.
 * @param [in]     B     GEOSX array2d.
 * @param [in,out] C     GEOSX array1d.
 * @param [in]     alpha Optional scalar to multiply with
 *                       \p transpose(A)*\p B.
 * @param [in]     beta  Optional parameter to control the accumulation.
 *
 * @warning
 * Assumes that \p transpose(A) and \p B have compatible sizes and that
 * \p C already has the right size.
 *
 */
void matrixTMatrixMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                            arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                            arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                            real64 const alpha = 1.0,
                            real64 const beta = 0.0 );

/**
 * @brief Matrix-transpose(Matrix) product;
 * \p C = \p alpha*\p A*\p transpose(B) + \p beta*\p C.
 *
 * Computes matrix-transpose(matrix) product with optional scaling
 * and accumulation.
 *
 * @param [in]     A     GEOSX array2d.
 * @param [in]     B     GEOSX array2d.
 * @param [in,out] C     GEOSX array1d.
 * @param [in]     alpha Optional scalar to multiply with
 *                       \p A*transpose(B).
 * @param [in]     beta  Optional parameter to control the accumulation.
 *
 * @warning
 * Assumes that \p A and \p transpose(B) have compatible sizes and that
 * \p C already has the right size.
 *
 */

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static void matrixMatrixTMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                                   arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                                   arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                                   real64 const alpha = 1.0,
                                   real64 const beta = 0.0 )
{

  GEOSX_ASSERT_MSG( C.size( 0 ) == A.size( 0 ) &&
                    C.size( 1 ) == B.size( 0 ) &&
                    A.size( 1 ) == B.size( 1 ),
                    "Matrix dimensions not compatible for product" );

#if !defined(__CUDA_ARCH__)

  // TODO: move to the helper file

  int const M = LvArray::integerConversion< int >( A.size( 0 ) );
  int const N = LvArray::integerConversion< int >( B.size( 0 ) );
  int const K = LvArray::integerConversion< int >( A.size( 1 ) );

  // A*B^T = C is computed as B * A^T = C^T, i.e. accessing the transpose
  // matrices using a column-major layout

  char const TRANS1 = 'T';
  char const TRANS2 = 'N';

  GEOSX_dgemm( &TRANS1, &TRANS2, &N, &M, &K, &alpha, B.dataIfContiguous(), &K, A.dataIfContiguous(), &K, &beta, C.dataIfContiguous(), &N );

#else

  // TODO: template on dimensions, move to LvArray

  for( localIndex i = 0; i < A.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < B.size( 0 ); ++j )
    {
      real64 sum_a_ik_b_jk = 0;
      for( localIndex k = 0; k < A.size( 1 ); ++k )
      {
        sum_a_ik_b_jk = sum_a_ik_b_jk + A( i, k ) * B( j, k );
      }
      C( i, j ) = beta * C( i, j ) + alpha * sum_a_ik_b_jk;
    }
  }

#endif
}

/**
 * @brief transpose(Matrix)-transpose(Matrix) product;
 * \p C = \p alpha*\p transpose(A)*\p transpose(B) + \p beta*\p C.
 *
 * Computes transpose(matrix)-transpose(matrix) product with optional
 * scaling and accumulation.
 *
 * @param [in]     A     GEOSX array2d.
 * @param [in]     B     GEOSX array2d.
 * @param [in,out] C     GEOSX array1d.
 * @param [in]     alpha Optional scalar to multiply with
 *                       \p tranpose(A)*transpose(B).
 * @param [in]     beta  Optional parameter to control the accumulation.
 *
 * @warning
 * Assumes that \p tranpose(A) and \p transpose(B) have compatible sizes
 * and that \p C already has the right size.
 *
 */
void matrixTMatrixTMultiply( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                             arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & B,
                             arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & C,
                             real64 const alpha = 1.0,
                             real64 const beta = 0.0 );

/**
 * @brief Computes the inverse matrix;
 * \p Ainv = \p inverse(A).
 *
 * Computes the inverse of the square matrix \p A and stores it in \p Ainv,
 * and also returns the determinant of \p A.
 *
 * @param [in]  A    GEOSX array2d.
 * @param [out] Ainv GEOSX array2d.
 * @param [out] detA Determinant of matrix \p A.
 *
 * @warning
 * Assumes \p Ainv already has the same size as \p A.
 *
 * @note This function is hardcoded for matrices up to order three.
 * For dimensions larger than three, the function calls LAPACK
 * functions DGETRF and DGETRI. Because of the row major
 * ordering used by GEOSX array2d, the inverse of \p A is practically
 * computed as the transpose matrix of the transpose matrix inverse using
 * lapack operation based on column-major ordering. This removes the need
 * for any copy/transposition that would be required operating with the
 * row-major layout.
 */
void matrixInverse( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                    arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & Ainv,
                    real64 & detA );

/**
 * @copydoc matrixInverse( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const &,
                           arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> &,
                           real64 & )
 */
void matrixInverse( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                    arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & Ainv,
                    real64 & detA );

/**
 * @brief Computes the inverse matrix;
 * \p Ainv = \p inverse(A).
 *
 * Computes the inverse of the square matrix \p A and stores it in \p Ainv.
 *
 * @param [in]  A    GEOSX array2d.
 * @param [out] Ainv GEOSX array2d.
 *
 * @warning
 * Assumes \p Ainv already has the same size as \p A.
 *
 * @note This function is hardcoded for matrices up to order three.
 * For dimensions larger than three, the function calls LAPACK
 * functions DGETRF and DGETRI. Because of the row major
 * ordering used by GEOSX array2d, the inverse of \p A is practically
 * computed as the transpose matrix of the transpose matrix inverse using
 * lapack operation based on column-major ordering. This removes the need
 * for any copy/transposition that would be required operating with the
 * row-major layout.
 */
void matrixInverse( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                    arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & Ainv );

/**
 * @copydoc matrixInverse( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const & A,
                           arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> & Ainv )
 */
void matrixInverse( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                    arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & Ainv );

/**
 * @brief Vector copy;
 * \p Y = \p X.
 *
 * @param [in]  X GEOSX array1d.
 * @param [out] Y GEOSX array1d.
 *
 * @warning
 * Assumes that \p X and \p Y have the same size.
 *
 */
void vectorCopy( array1d< real64 > const & X,
                 array1d< real64 > & Y );

/**
 * @brief Matrix copy;
 * \p B = \p A.
 *
 * @param [in]  A GEOSX array2d.
 * @param [out] B GEOSX array2d.
 *
 * @warning
 * Assumes that \p A and \p B have the same size.
 *
 */
void matrixCopy( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                 arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & B );

/**
 * @copydoc matrixCopy( arraySlice2d<real64 const, MatrixLayout::ROW_MAJOR> const & A,
                        arraySlice2d<real64, MatrixLayout::ROW_MAJOR> const & B )
 */
void matrixCopy( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                 arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & B );

/**
 * @brief Set the random number generator seed.
 *
 * @param [in] seed GEOSX array1d, dimension (4). The elements must be
 *             between 0 and 4095, and seed(4) must be odd.
 */
void setRandomNumberGeneratorSeed( arraySlice1d< int const > const & seed );

/**
 * @brief Get the random number generator seed.
 *
 * @param [out] seed GEOSX array1d, dimension (4).
 */
void getRandomNumberGeneratorSeed( arraySlice1d< int > const & seed );

/**
 * @brief Sets vector entries to random real numbers.
 *
 * Sets vector entries to random real numbers from a uniform or normal
 * distribution without specifying the seed of the random number generator.
 *
 * @param [out]    X     GEOSX array1d.
 * @param [in]     idist Optional RandomNumberDistribution enum value
 *                       specifying the distribution of the random numbers.
 */
void vectorRand( arraySlice1d< real64 > const & X,
                 RandomNumberDistribution const & idist = RandomNumberDistribution::UNIFORM_01 );

/**
 * @brief Sets matrix entries to random real numbers.
 *
 * Sets matrix entries to random real numbers from a uniform or normal
 * distribution without specifying the seed of the random number generator.
 *
 * @param [out]    X     GEOSX array1d.
 * @param [in]     idist Optional RandomNumberDistribution enum value
 *                       specifying the distribution of the random numbers.
 */
void matrixRand( arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & A,
                 RandomNumberDistribution const & idist = RandomNumberDistribution::UNIFORM_01 );

/**
 * @copydoc matrixRand( arraySlice2d<real64, MatrixLayout::ROW_MAJOR> const & A,
                        RandomNumberDistribution const & idist = RandomNumberDistribution::UNIFORM_01 )
 */
void matrixRand( arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & A,
                 RandomNumberDistribution const & idist = RandomNumberDistribution::UNIFORM_01 );


/**
 * @brief Computes the singular value decomposition of A
 *
 * If size(A) = (M,N), this function expects:
 * size(U) = (M,min(M,N)),
 * size(VT) = (min(M,N),N), and
 * size(S) = min(M,N)
 * On exit, S contains the singular values of A
 * and U contains an orthonormal basis of range(A)
 *
 * @param [in]    A GEOSX array2d.
 * @param [out]   U GEOSX array2d.
 * @param [out]   S GEOSX array1d.
 * @param [out]   VT GEOSX array2d.
 */
void matrixSVD( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & A,
                arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & U,
                arraySlice1d< real64 > const & S,
                arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & VT );

/**
 * @copydoc matrixSVD( arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const &,
                       arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const &,
                       arraySlice1d< real64 > const &,
                       arraySlice2d< real64, MatrixLayout::ROW_MAJOR > const & )
 */
void matrixSVD( arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & A,
                arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & U,
                arraySlice1d< real64 > const & S,
                arraySlice2d< real64, MatrixLayout::COL_MAJOR > const & VT );


} // namespace BlasLapackLA

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_BLASLAPACKLA_HPP_*/
