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


namespace geosx
{

/**
 * \class BlasLapackLA
 * \brief This class contains a collection of BLAS and LAPACK linear
 *        algebra operations for GEOSX array1d and array2d
 */
class BlasLapackLA
{

public:

  /**
    * \enum  RandomNumberDistribution
    * \brief This enum class specifies the type of distribution for
    *        generating random real numbers.
    */
  enum class RandomNumberDistribution : int
  {
    UNIFORM_01 = 1,   /**< uniform distribution (0,1); */
    UNIFORM_m1p1 = 2, /**< uniform distribution (-1,1); */
    NORMAL_01 = 3     /**< normal distribution (0,1); */
  };

  /**
   * @brief Returns the 1-norm of the vector.
   *
   * @param [in] X GEOSX array1d.
   */
  static real64 vectorNorm1( array1d<real64> const & X ) ;

  /**
   * @brief Returns the two norm of the vector.
   *
   * @param [in] X GEOSX array1d.
   */
  static real64 vectorNorm2( array1d<real64> const & X ) ;

  /**
   * @brief Returns the infinity-norm of the vector.
   *
   * @param [in] X GEOSX array1d.
   */
  static real64 vectorNormInf( array1d<real64> const & X ) ;

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
  static real64 determinant( array2d<real64> const & A ) ;

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
  static real64 matrixNormInf(array2d<real64> const & A) ;

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
  static real64 matrixNorm1(array2d<real64> const & A) ;

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
  static real64 matrixNormFrobenius(array2d<real64> const & A) ;

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
  static void vectorVectorAdd( array1d<real64> const & X,
                               array1d<real64> & Y,
                               real64 const alpha = 1. );

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
  static void matrixMatrixAdd( array2d<real64> const & A,
                               array2d<real64> & B,
                               real64 const alpha = 1. ) ;

  /**
   * @brief In-place scalar-vector product;
   * \p X = \p alpha*\p X.
   *
   * @param [in]     alpha Scalar to multiply with \p X.
   * @param [in,out] X     GEOSX array1d.
   */
  static void vectorScale( real64 alpha,
                           array1d<real64> & X ) ;

  /**
   * @brief In-place scalar-matrix product;
   * \p A = \p alpha*\p A.
   *
   * @param [in]     alpha Scalar to multiply with <tt>A</tt>.
   * @param [in,out] A     GEOSX array2d.
   */
  static void matrixScale( real64 alpha,
                           array2d<real64> & A ) ;

  /**
   * @brief Returns the dot product of two vectors.
   *
   * @param [in] X GEOSX array1d.
   * @param [in] Y GEOSX array1d.
   *
   */
  static real64 vectorDot( array1d<real64> const & X,
                           array1d<real64> const & Y) ;

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
  static void matrixVectorMultiply(array2d<real64> const & A,
                                   array1d<real64> const & X,
                                   array1d<real64>  & Y,
                                   real64 const alpha=1.,
                                   real64 const beta=0.) ;

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
  static void matrixTVectorMultiply(array2d<real64> const & A,
                                    array1d<real64> const & X,
                                    array1d<real64>  & Y,
                                    real64 const alpha=1.,
                                    real64 const beta=0.) ;

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
  static void matrixMatrixMultiply( array2d<real64> const & A,
                                    array2d<real64> const & B,
                                    array2d<real64> & C,
                                    real64 const alpha=1.,
                                    real64 const beta=0.) ;

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
  static void matrixTMatrixMultiply( array2d<real64> const & A,
                                     array2d<real64> const & B,
                                     array2d<real64> & C,
                                     real64 const alpha=1.,
                                     real64 const beta=0.);

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
  static void matrixMatrixTMultiply( array2d<real64> const & A,
                                     array2d<real64> const & B,
                                     array2d<real64> & C,
                                     real64 const alpha=1.,
                                     real64 const beta=0.) ;

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
  static void matrixTMatrixTMultiply( array2d<real64> const & A,
                                      array2d<real64> const & B,
                                      array2d<real64> & C,
                                      real64 const alpha=1.,
                                      real64 const beta=0.) ;

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
  static void matrixInverse( array2d<real64> const & A,
                             array2d<real64> & Ainv ) ;

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
  static void matrixInverse( array2d<real64> const & A,
                             array2d<real64> & Ainv,
                             real64 & detA) ;

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
  static void vectorCopy( array1d<real64> const & X,
                          array1d<real64> & Y ) ;

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
  static void matrixCopy( array2d<real64> const & A,
                          array2d<real64> & B ) ;

  /**
   * @brief Set the random number generator seed.
   *
   * @param [in] seed GEOSX array1d, dimension (4). The elements must be
   *             between 0 and 4095, and seed(4) must be odd.
   */
  static void setRandomNumberGeneratorSeed( array1d<int> const & seed);

  /**
   * @brief Get the random number generator seed.
   *
   * @param [out] seed GEOSX array1d, dimension (4).
   */
  static void getRandomNumberGeneratorSeed( array1d<int> & seed);

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
  static void vectorRand( array1d<real64> & X,
                          RandomNumberDistribution const & idist = RandomNumberDistribution::UNIFORM_01) ;

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
  static void matrixRand( array2d<real64> & A,
                          RandomNumberDistribution const & idist = RandomNumberDistribution::UNIFORM_01) ;


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
  static void matrixSVD( array2d<real64> const & A,
	                 array2d<real64> & U,
	   	         array1d<real64> & S,
		         array2d<real64> & VT ) ;

  
};

}

#endif /*GEOSX_LINEARALGEBRA_BLASLAPACKLA_HPP_*/
