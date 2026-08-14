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

/**
 * @file denseLASolvers.hpp
 */
#ifndef GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_
#define GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_

#include "common/DataTypes.hpp"
#include "denseLinearAlgebra/common/layouts.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/logger/Logger.hpp"

#include <complex>

namespace geos
{

namespace denseLinearAlgebra
{

namespace details
{

constexpr real64 singularMatrixTolerance = 1e2*LvArray::NumericLimits< real64 >::epsilon;

/**
 * @brief Solves a 2x2 linear system A * x = b.
 *
 * This function solves a linear system of the form A * x = b, where A is a 2x2 matrix,
 * b is a 2x1 vector, and x is the solution vector. The function checks the sizes
 * of the inputs to ensure they conform to the expected dimensions. It also checks that
 * the determinant of matrix A is not near zero to avoid solving a singular system.
 *
 * @tparam MATRIX_TYPE The type of the matrix A. Must support indexing with `A[i][j]`.
 * @tparam RHS_TYPE The type of the right-hand side vector b. Must support indexing with `b[i]`.
 * @tparam SOL_TYPE The type of the solution vector x. Must support indexing with `x[i]`.
 *
 * @param[in] A The 2x2 matrix representing the system of equations. Must have size 2x2.
 * @param[in] b The 2-element vector representing the right-hand side of the equation.
 * @param[out] x The 2-element vector that will store the solution to the system.
 * @return bool that specifies whether the solve succeeded (1) or not (0).
 */
template< typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveTwoByTwoSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  LvArray::tensorOps::internal::checkSizes< 2, 2 >( A );
  LvArray::tensorOps::internal::checkSizes< 2 >( b );
  LvArray::tensorOps::internal::checkSizes< 2 >( x );

  real64 const detA = LvArray::tensorOps::determinant< 2 >( A );

  if( LvArray::math::abs( detA ) < singularMatrixTolerance )
    return false;

  real64 const invA = 1.0 / detA;

  x[0] = ( A[1][1] * b[0] - A[0][1] * b[1] ) * invA;
  x[1] = ( A[0][0] * b[1] - A[1][0] * b[0] ) * invA;

  return true;
}

/**
 * @brief Solves the linear system `A * X = B`, where A is a 2×2 matrix and X and B are 2×M matrices.
 *
 * @details
 * This function solves a set of linear systems of the form `A * X = B`, where:
 * - `A` is a 2×2 coefficient matrix,
 * - `B` is a 2×M matrix of right-hand sides,
 * - `X` is the solution matrix, also 2×M.
 *
 * On input, the matrix `X` contains the right-hand side values (i.e., `B`). On output,
 * it is overwritten with the solution vectors (i.e., `X`). The function checks for
 * singularity using a tolerance and returns `false` if the system cannot be solved.
 *
 * @tparam M The number of columns in the rhs matrix.
 * @tparam MATRIX_TYPE The type of the matrix A. Must support indexing with `A[i][j]`.
 * @tparam SOL_TYPE The type of the rhs/solution matrix x. Must support indexing with `X[i][j]`.
 *
 * @param[in] A The 2x2 matrix representing the system of equations. Must have size 2x2.
 * @param[in,out] X The 2xM matrix that will store the right-hand-side and solution to the system.
 * @return bool that specifies whether the solve succeeded (1) or not (0).
 */
template< std::ptrdiff_t M,
          typename MATRIX_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveTwoByTwoSystem( MATRIX_TYPE const & A, SOL_TYPE && X )
{
  LvArray::tensorOps::internal::checkSizes< 2, 2 >( A );
  LvArray::tensorOps::internal::checkSizes< 2, M >( X );

  real64 const detA = LvArray::tensorOps::determinant< 2 >( A );

  if( LvArray::math::abs( detA ) < singularMatrixTolerance )
    return false;

  real64 const invA = 1.0 / detA;

  for( std::ptrdiff_t j = 0; j < M; j++ )
  {
    real64 const b0 = X[0][j];
    real64 const b1 = X[1][j];
    X[0][j] = ( A[1][1] * b0 - A[0][1] * b1 ) * invA;
    X[1][j] = ( A[0][0] * b1 - A[1][0] * b0 ) * invA;
  }
  return true;
}

/**
 * @brief Solves a 3x3 linear system A * x = b.
 *
 * This function solves a linear system of the form A * x = b, where A is a 3x3 matrix,
 * b is a 3x1 vector, and x is the solution vector. The function checks the sizes
 * of the inputs to ensure they conform to the expected dimensions. It also checks that
 * the determinant of matrix A is not near zero to avoid solving a singular system.
 *
 * @tparam MATRIX_TYPE The type of the matrix A. Must support indexing with `A[i][j]`.
 * @tparam RHS_TYPE The type of the right-hand side vector b. Must support indexing with `b[i]`.
 * @tparam SOL_TYPE The type of the solution vector x. Must support indexing with `x[i]`.
 *
 * @param[in] A The 3x3 matrix representing the system of equations. Must have size 3x3.
 * @param[in] b The 3-element vector representing the right-hand side of the equation.
 * @param[out] x The 3-element vector that will store the solution to the system.
 * @return bool that specifies whether the solve succeeded (1) or not (0).
 */
template< typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveThreeByThreeSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  LvArray::tensorOps::internal::checkSizes< 3, 3 >( A );
  LvArray::tensorOps::internal::checkSizes< 3 >( b );
  LvArray::tensorOps::internal::checkSizes< 3 >( x );

  real64 const detA = LvArray::tensorOps::determinant< 3 >( A );

  if( LvArray::math::abs( detA ) < singularMatrixTolerance )
    return false;

  real64 const invA = 1.0 / detA;

  real64 const detX0 = b[0] * ( A[1][1] * A[2][2] - A[2][1] * A[1][2] ) -
                       b[1] * ( A[0][1] * A[2][2] - A[0][2] * A[2][1] ) +
                       b[2] * ( A[0][1] * A[1][2] - A[0][2] * A[1][1] );

  real64 const detX1 = A[0][0] * ( b[1] * A[2][2] - b[2] * A[1][2] ) -
                       A[1][0] * ( b[0] * A[2][2] - b[2] * A[0][2] ) +
                       A[2][0] * ( b[0] * A[1][2] - b[1] * A[0][2] );

  real64 const detX2 = A[0][0] * ( A[1][1] * b[2] - A[2][1] * b[1] ) -
                       A[1][0] * ( A[0][1] * b[2] - A[2][1] * b[0] ) +
                       A[2][0] * ( A[0][1] * b[1] - A[1][1] * b[0] );

  x[0] = detX0 * invA;
  x[1] = detX1 * invA;
  x[2] = detX2 * invA;

  return true;
}

/**
 * @brief Solves the linear system `A * X = B`, where A is a 3x3 matrix and X and B are 3×M matrices.
 *
 * @details
 * This function solves a set of linear systems of the form `A * X = B`, where:
 * - `A` is a 3x3 coefficient matrix,
 * - `B` is a 3×M matrix of right-hand sides,
 * - `X` is the solution matrix, also 3×M.
 *
 * On input, the matrix `X` contains the right-hand side values (i.e., `B`). On output,
 * it is overwritten with the solution vectors (i.e., `X`). The function checks for
 * singularity using a tolerance and returns `false` if the system cannot be solved.
 *
 * @tparam M The number of columns in the rhs matrix.
 * @tparam MATRIX_TYPE The type of the matrix A. Must support indexing with `A[i][j]`.
 * @tparam SOL_TYPE The type of the rhs/solution matrix x. Must support indexing with `X[i][j]`.
 *
 * @param[in] A The 3x3 matrix representing the system of equations. Must have size 3x3.
 * @param[in,out] X The 3xM matrix that will store the right-hand-side and solution to the system.
 * @return bool that specifies whether the solve succeeded (1) or not (0).
 */
template< std::ptrdiff_t M,
          typename MATRIX_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveThreeByThreeSystem( MATRIX_TYPE const & A, SOL_TYPE && X )
{
  LvArray::tensorOps::internal::checkSizes< 3, 3 >( A );
  LvArray::tensorOps::internal::checkSizes< 3, M >( X );

  real64 const detA = LvArray::tensorOps::determinant< 3 >( A );

  if( LvArray::math::abs( detA ) < singularMatrixTolerance )
    return false;

  real64 const invA = 1.0 / detA;

  for( std::ptrdiff_t j = 0; j < M; j++ )
  {
    real64 const b0 = X[0][j];
    real64 const b1 = X[1][j];
    real64 const b2 = X[2][j];

    real64 const detX0 = b0 * ( A[1][1] * A[2][2] - A[2][1] * A[1][2] ) -
                         b1 * ( A[0][1] * A[2][2] - A[0][2] * A[2][1] ) +
                         b2 * ( A[0][1] * A[1][2] - A[0][2] * A[1][1] );

    real64 const detX1 = A[0][0] * ( b1 * A[2][2] - b2 * A[1][2] ) -
                         A[1][0] * ( b0 * A[2][2] - b2 * A[0][2] ) +
                         A[2][0] * ( b0 * A[1][2] - b1 * A[0][2] );

    real64 const detX2 = A[0][0] * ( A[1][1] * b2 - A[2][1] * b1 ) -
                         A[1][0] * ( A[0][1] * b2 - A[2][1] * b0 ) +
                         A[2][0] * ( A[0][1] * b1 - A[1][1] * b0 );

    X[0][j] = detX0 * invA;
    X[1][j] = detX1 * invA;
    X[2][j] = detX2 * invA;
  }
  return true;
}

/**
 * @brief Solves a linear system where the matrix is upper triangular using back substitution.
 *
 * This function solves the linear system `Ax = b`, where `A` is an upper triangular matrix, using
 * back substitution. The solution `x` is computed and stored in the provided output vector.
 *
 * @tparam N The size of the square matrix `A`.
 * @tparam MATRIX_TYPE The type of the matrix `A`.
 * @tparam RHS_TYPE The type of the right-hand side vector `b`.
 * @tparam SOL_TYPE The type of the solution vector `x`.
 * @param[in] A The upper triangular matrix representing the coefficients of the system.
 * @param[in] b The right-hand side vector. It is used to compute the solution.
 * @param[out] x The solution vector. The result of solving the system `Ax = b` using back substitution.
 */
template< std::ptrdiff_t N,
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
void solveUpperTriangularSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  for( std::ptrdiff_t i = N - 1; i >= 0; --i )
  {
    real64 sum = b[i];
    for( std::ptrdiff_t j = i + 1; j < N; ++j )
    {
      sum -= A[i][j] * x[j];
    }
    x[i] = sum / A[i][i];
  }
}

/**
 * @brief Solves a linear system with an upper triangular coefficient matrix.
 *
 * @details
 * This function performs back-substitution to solve the linear system `A * X = B`,
 * where `A` is an upper triangular matrix of size NxN, and `X` is a
 * solution matrix of size `NxM`. On input, `X` contains the right-hand side vectors (B),
 * and on output, it is overwritten with the solution vectors (X).
 *
 * The function assumes that all diagonal elements of `A` are non-zero.
 *
 * @tparam N Number of rows and columns in the square matrix `A`.
 * @tparam M Number of columns in the right-hand side and solution matrix `X`.
 * @tparam MATRIX_TYPE Type of the matrix `A`. Must support indexing via `A[i][j]`.
 * @tparam SOL_TYPE Type of the matrix `X`. Must support indexing via `X[i][j]`.
 *
 * @param[in] A Upper triangular matrix of size NxN.
 * @param[in,out] X Matrix of size NxM. On input, contains the right-hand sides;
 *                  on output, contains the solutions.
 */
template< std::ptrdiff_t N,
          std::ptrdiff_t M,
          typename MATRIX_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
void solveUpperTriangularSystem( MATRIX_TYPE const & A, SOL_TYPE && X )
{
  for( std::ptrdiff_t i = N - 1; i >= 0; --i )
  {
    real64 const invAii = 1.0 / A[i][i];
    for( std::ptrdiff_t k = 0; k < M; ++k )
    {
      real64 sum = X[i][k];
      for( std::ptrdiff_t j = i + 1; j < N; ++j )
      {
        sum -= A[i][j] * X[j][k];
      }
      X[i][k] = sum * invAii;
    }
  }
}

/**
 * @brief Solves a linear system using Gaussian elimination.
 *
 * This function performs Gaussian elimination on the given matrix `A` and right-hand side vector `b`.
 * It transforms the matrix `A` into an upper triangular matrix and then solves for the solution `x`
 * using back substitution.
 *
 * @tparam N The size of the square matrix `A`.
 * @tparam MATRIX_TYPE The type of the matrix `A`.
 * @tparam RHS_TYPE The type of the right-hand side vector `b`.
 * @tparam SOL_TYPE The type of the solution vector `x`.
 * @param[in,out] A The matrix to be transformed into an upper triangular matrix. Modified in place.
 * @param[in,out] b The right-hand side vector. Modified in place to reflect the transformed system.
 * @param[out] x The solution vector. The result of solving the system `Ax = b`.
 * @return bool that specifies whether the solve succeeded (1) or not (0).
 */
template< std::ptrdiff_t N,
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveGaussianElimination( MATRIX_TYPE & A, RHS_TYPE & b, SOL_TYPE && x )
{
  static_assert( N > 0, "N must be greater than 0." );
  LvArray::tensorOps::internal::checkSizes< N, N >( A );
  LvArray::tensorOps::internal::checkSizes< N >( b );
  LvArray::tensorOps::internal::checkSizes< N >( x );

  // Step 1: Transform  into an upper triangular matrix

  // 1.a. Find the pivot
  for( std::ptrdiff_t i = 0; i < N; ++i )
  {
    std::ptrdiff_t max_row = i;
    for( std::ptrdiff_t k = i + 1; k < N; ++k )
    {
      if( LvArray::math::abs( A[k][i] ) > LvArray::math::abs( A[max_row][i] ))
      {
        max_row = k;
      }
    }

    // 1.b. Swap rows
    for( std::ptrdiff_t k = i; k < N; ++k )
    {
      // std::swap( A[i][k], A[max_row][k] );
      real64 const temp = A[max_row][k];
      A[max_row][k] = A[i][k];
      A[i][k] = temp;
    }
    // std::swap( b[i], b[max_row] ); cannot be done on device
    real64 const temp = b[i];
    b[i] =  b[max_row];
    b[max_row] = temp;


    if( LvArray::math::abs( A[i][i] ) < singularMatrixTolerance )
      return false;

    // 1.c Eliminate entries below the pivot
    for( std::ptrdiff_t k = i + 1; k < N; ++k )
    {
      real64 const scaling = A[k][i] / A[i][i];
      for( std::ptrdiff_t j = i; j < N; ++j )
      {
        A[k][j] -= scaling * A[i][j];
      }
      b[k] -= scaling * b[i];
    }
  }

  // Step 2: Backward substitution
  solveUpperTriangularSystem< N >( A, b, std::forward< SOL_TYPE >( x ) );

  return true;
}

/**
 * @brief Solves a linear system using Gaussian elimination.
 *
 * @details
 * This function solves the system `A * X = B`, where:
 * - `A` is an `NxN` matrix,
 * - `X` is an `NxM` matrix containing the right-hand sides on input,
 *   and the solutions on output.
 *
 * The function performs Gaussian elimination with partial pivoting to transform
 * the matrix into upper triangular form, followed by back-substitution to compute
 * the solution. It returns `false` if the matrix is singular (i.e., pivot is below
 * a tolerance threshold).
 *
 * @tparam N Number of rows and columns in the square matrix `A`.
 * @tparam M Number of columns in the right-hand side and solution matrix `X`.
 * @tparam MATRIX_TYPE Type of the matrix `A`. Must support indexing via `A[i][j]`.
 * @tparam SOL_TYPE Type of the matrix `X`. Must support indexing via `X[i][j]`.
 *
 * @param[in,out] A The `NxN` coefficient matrix. Modified in-place to upper triangular form.
 * @param[in,out] X The `NxM` matrix. On input, contains the right-hand sides; on output, the solutions.
 * @return `true` if the system was successfully solved; `false` if the matrix is singular.
 */
template< std::ptrdiff_t N,
          std::ptrdiff_t M,
          typename MATRIX_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
bool solveGaussianElimination( MATRIX_TYPE & A, SOL_TYPE && X )
{
  static_assert( N > 0, "N must be greater than 0." );
  static_assert( M > 0, "N must be greater than 0." );
  LvArray::tensorOps::internal::checkSizes< N, N >( A );
  LvArray::tensorOps::internal::checkSizes< N, M >( X );

  // Step 1: Transform  into an upper triangular matrix

  // 1.a. Find the pivot
  for( std::ptrdiff_t i = 0; i < N; ++i )
  {
    std::ptrdiff_t max_row = i;
    for( std::ptrdiff_t k = i + 1; k < N; ++k )
    {
      if( LvArray::math::abs( A[k][i] ) > LvArray::math::abs( A[max_row][i] ))
      {
        max_row = k;
      }
    }

    // 1.b. Swap rows
    for( std::ptrdiff_t k = i; k < N; ++k )
    {
      // std::swap( A[i][k], A[max_row][k] );
      real64 const temp = A[max_row][k];
      A[max_row][k] = A[i][k];
      A[i][k] = temp;
    }

    if( LvArray::math::abs( A[i][i] ) < singularMatrixTolerance )
      return false;

    // std::swap( X[i], X[max_row] ); cannot be done on device
    for( std::ptrdiff_t j = 0; j < M; j++ )
    {
      real64 const temp = X[i][j];
      X[i][j] = X[max_row][j];
      X[max_row][j] = temp;
    }

    // 1.c Eliminate entries below the pivot
    for( std::ptrdiff_t k = i + 1; k < N; ++k )
    {
      real64 const scaling = A[k][i] / A[i][i];
      for( std::ptrdiff_t j = i; j < N; ++j )
      {
        A[k][j] -= scaling * A[i][j];
      }
      for( std::ptrdiff_t j = 0; j < M; j++ )
      {
        X[k][j] -= scaling * X[i][j];
      }
    }
  }

  // Step 2: Backward substitution
  solveUpperTriangularSystem< N, M >( A, std::forward< SOL_TYPE >( X ) );

  return true;
}

} // details namespace

/**
 * @brief Solves a linear system using the most appropriate method based on the size of the system.
 *
 * This function determines the appropriate method for solving a linear system `Ax = b` based on
 * the size of the matrix `A`. For 2x2 and 3x3 systems, specialized solvers are used. For larger systems,
 * Gaussian elimination is employed. The matrix and the rhs are modified by the function.
 *
 * @tparam N The size of the square matrix `A`.
 * @tparam MATRIX_TYPE The type of the matrix `A`.
 * @tparam RHS_TYPE The type of the right-hand side vector `b`.
 * @tparam SOL_TYPE The type of the solution vector `x`.
 * @tparam MODIFY_MATRIX boolean flag indicating whether the input matrix `A` and vector `b` should be modified.
 *                       If `1`, the matrix `A` and vector `b` are modified in place. If `0`, copies of
 *                       `A` and `b` are made, and the original data is left unchanged.
 * @param[in] A The matrix representing the coefficients of the system.
 * @param[in] b The right-hand side vector.
 * @param[out] x The solution vector. The result of solving the system `Ax = b`.
 * @return bool that specifies whether the solve succeeded (1) or not (0).
 */
template< std::ptrdiff_t N,
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE,
          bool MODIFY_MATRIX = true >
GEOS_HOST_DEVICE
inline
bool solve( MATRIX_TYPE & A, RHS_TYPE & b, SOL_TYPE && x )
{
  static_assert( N > 0, "N must be greater than 0." );
  static_assert( N < 11, "N cannot be larger than 10" );
  LvArray::tensorOps::internal::checkSizes< N, N >( A );
  LvArray::tensorOps::internal::checkSizes< N >( b );
  LvArray::tensorOps::internal::checkSizes< N >( x );

  if constexpr ( N == 2 )
  {
    return details::solveTwoByTwoSystem( A, b, std::forward< SOL_TYPE >( x ) );
  }
  else if constexpr ( N == 3 )
  {
    return details::solveThreeByThreeSystem( A, b, std::forward< SOL_TYPE >( x ) );
  }
  else
  {
    if constexpr ( MODIFY_MATRIX )
    {
      return details::solveGaussianElimination< N >( A, b, std::forward< SOL_TYPE >( x ) );
    }
    else
    {
      real64 A_copy[N][N]{};
      real64 b_copy[N]{};

      for( std::ptrdiff_t i=0; i < N; ++i )
      {
        b_copy[i] = b[i];
        for( std::ptrdiff_t j=0; j < N; ++j )
        {
          A_copy[i][j] = A[i][j];
        }
      }
      return details::solveGaussianElimination< N >( A_copy, b_copy, std::forward< SOL_TYPE >( x ) );
    }
  }
}

/**
 * @brief Solves a multiple column linear system using the most appropriate method based on the size of the system.
 *
 * This function determines the appropriate method for solving a linear system `AX = B` based on
 * the size of the matrix `A`. For 2x2 and 3x3 systems, specialized solvers are used. For larger systems,
 * Gaussian elimination is employed. The matrix and the rhs are modified by the function.
 *
 * @tparam N The size of the square matrix `A`. Also the number of rows in B and X
 * @tparam M The number of columns of the `B` and `X`
 * @tparam MATRIX_TYPE The type of the matrix `A`.
 * @tparam RHS_TYPE The type of the right-hand side matrix `B`.
 * @tparam SOL_TYPE The type of the solution matrix `X`.
 * @tparam MODIFY_MATRIX boolean flag indicating whether the input matrix `A` and matrix `B` should be modified.
 *                       If `1`, the matrix `A` and matrix `B` are modified in place. If `0`, copies of
 *                       `A` and `B` are made, and the original data is left unchanged.
 *
 * @param[in] A The matrix representing the coefficients of the system (NxN).
 * @param[in] B The right-hand side matrix (NxM).
 * @param[out] X The solution vector. The result of solving the system `AX = B` (NxM).
 * @return bool that specifies whether the solve succeeded (1) or not (0).
 */
template< std::ptrdiff_t N,
          std::ptrdiff_t M,
          typename MATRIX_TYPE,
          typename SOL_TYPE,
          bool MODIFY_MATRIX = true >
GEOS_HOST_DEVICE
inline
bool solve( MATRIX_TYPE & A, SOL_TYPE && X )
{
  static_assert( N > 0, "N must be greater than 0." );
  static_assert( N < 11, "N cannot be larger than 10" );
  static_assert( M > 0, "M must be greater than 0." );
  static_assert( M < 12, "M cannot be larger than 11" );
  LvArray::tensorOps::internal::checkSizes< N, N >( A );
  LvArray::tensorOps::internal::checkSizes< N, M >( X );

  if constexpr ( N == 2 )
  {
    return details::solveTwoByTwoSystem< M >( A, std::forward< SOL_TYPE >( X ) );
  }
  else if constexpr ( N == 3 )
  {
    return details::solveThreeByThreeSystem< M >( A, std::forward< SOL_TYPE >( X ) );
  }
  else if constexpr ( MODIFY_MATRIX )
  {
    return details::solveGaussianElimination< N, M >( A, std::forward< SOL_TYPE >( X ) );
  }
  else
  {
    real64 A_copy[N][N];

    for( std::ptrdiff_t i=0; i < N; ++i )
    {
      for( std::ptrdiff_t j=0; j < N; ++j )
      {
        A_copy[i][j] = A[i][j];
      }
    }
    return details::solveGaussianElimination< N, M >( A_copy, std::forward< SOL_TYPE >( X ) );
  }
}

/**
 * @brief Solves a linear system `A * X = B` using Gaussian elimination.
 *
 * @details
 * This function solves a system of linear equations where:
 * - `A` is an `NxN` coefficient matrix,
 * - `B` is an `NxM` right-hand side matrix,
 * - `X` is the solution matrix of size `NxM`.
 *
 * The function first copies the contents of `B` into `X`, then calls the
 * Gaussian elimination solver. The matrix `A` may be modified in-place depending
 * on the `MODIFY_MATRIX` flag.
 *
 * @tparam N Number of rows and columns in the square matrix `A`.
 * @tparam M Number of columns in the right-hand side and solution matrices.
 * @tparam MATRIX_TYPE Type of the matrix `A`. Must support indexing via `A[i][j]`.
 * @tparam RHS_TYPE Type of the matrix `B`. Must support indexing via `B[i][j]`.
 * @tparam SOL_TYPE Type of the matrix `X`. Must support indexing via `X[i][j]`.
 * @tparam MODIFY_MATRIX If true, allows in-place modification of matrix `A`.
 *
 * @param[in,out] A The `NxN` coefficient matrix. May be modified if `MODIFY_MATRIX` is true.
 * @param[in] B The `NxM` right-hand side matrix.
 * @param[out] X The `NxM` solution matrix. On input, overwritten with `B`; on output, contains the solution.
 * @return `true` if the system was successfully solved; `false` if the matrix is singular.
 */
template< std::ptrdiff_t N,
          std::ptrdiff_t M,
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE,
          bool MODIFY_MATRIX = true >
GEOS_HOST_DEVICE
inline
bool solve( MATRIX_TYPE & A, RHS_TYPE & B, SOL_TYPE && X )
{
  static_assert( N > 0, "N must be greater than 0." );
  static_assert( N < 11, "N cannot be larger than 10" );
  static_assert( M > 0, "M must be greater than 0." );
  static_assert( M < 12, "M cannot be larger than 11" );
  LvArray::tensorOps::internal::checkSizes< N, N >( A );
  LvArray::tensorOps::internal::checkSizes< N, M >( B );
  LvArray::tensorOps::internal::checkSizes< N, M >( X );

  // Copy the RHS onto the solution
  for( std::ptrdiff_t i=0; i < N; ++i )
  {
    for( std::ptrdiff_t j=0; j < M; ++j )
    {
      X[i][j] = B[i][j];
    }
  }
  return solve< N, M, MATRIX_TYPE, SOL_TYPE, MODIFY_MATRIX >( A, std::forward< SOL_TYPE >( X ) );
}

} // denseLinearAlgebra

} // geos


#endif /*GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_*/
