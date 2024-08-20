/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
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

#include <complex>

namespace geos
{

namespace denseLinearAlgebra
{

namespace internal
{  
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
 */
template< typename MATRIX_TYPE, 
          typename RHS_TYPE, 
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
void solveTwoByTwoSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x)
{
    LvArray::tensorOps::internal::checkSizes< 2, 2 >( A );
    LvArray::tensorOps::internal::checkSizes< 2 >( b );
    LvArray::tensorOps::internal::checkSizes< 2 >( x );

    real64 const detA = A[0][0] * A[1][1] - A[0][1] * A[1][0];

    GEOS_ERROR_IF_LT_MSG( LvArray::math::abs(detA), LvArray::NumericLimits< real64 >::epsilon;, "Singular system." );

    x[0] = (A[1][1] * b[0] - A[0][1] * b[1] ) / detA;
    x[1] = (A[0][0] * b[1] - A[1][0] * b[0] ) / detA;
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
 */
template< typename MATRIX_TYPE, 
          typename RHS_TYPE, 
          typename SOL_TYPE >
GEOS_HOST_DEVICE
inline
void solveThreeByThreeSystem( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x)
{
    LvArray::tensorOps::internal::checkSizes< 3, 3 >( A );
    LvArray::tensorOps::internal::checkSizes< 3 >( b );
    LvArray::tensorOps::internal::checkSizes< 3 >( x );

    real64 const detA = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                        A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                        A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    GEOS_ERROR_IF_LT_MSG( LvArray::math::abs(detA), LvArray::NumericLimits< real64 >::epsilon;, "Singular system." );

    real64 const detX0 = b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                         b[1] * (A[0][1] * A[2][2] - A[0][2] * A[2][1]) +
                         b[2] * (A[0][1] * A[1][2] - A[0][2] * A[1][1]);

    real64 const detX1 = A[0][0] * (b[1] * A[2][2] - b[2] * A[2][1]) -
                         A[0][1] * (b[0] * A[2][2] - b[2] * A[2][0]) +
                         A[0][2] * (b[0] * A[1][2] - b[1] * A[1][0]);

    real64 const detX2 = A[0][0] * (A[1][1] * b[2] - A[1][2] * b[1]) -
                         A[0][1] * (A[1][0] * b[2] - A[1][2] * b[0]) +
                         A[0][2] * (A[1][0] * b[1] - A[1][1] * b[0]);

    x[0] = detX0 / detA;
    x[1] = detX1 / detA;
    x[2] = detX2 / detA;
}


template< std::ptrdiff_t N, 
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
void solveUpperTriangularSystem( MATRIX_TYPE const & A, RHS_TYPE & b, SOL_TYPE && x )
{
  for (std::ptrdiff_t i = N - 1; i >= 0; --i)
  {
    real64 sum = b[i];
    for (std::ptrdiff_t j = i + 1; j < N; ++j)
    {
      sum -= A[i][j] * x[j];
    }
    x[i] = sum / A[i][i];
  }
}

template< std::ptrdiff_t N, 
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE 
inline
void solveGaussianElimination( MATRIX_TYPE & A, RHS_TYPE & b, SOL_TYPE && x )
{
  static_assert( N > 0, "N must be greater than 0." );
  internal::checkSizes< N, N >( matrix );
  internal::checkSizes< N >( b ); 
  internal::checkSizes< N >( x );
  

  // Step 1: Transform  into an upper triangular matrix 
  
  // 1.a. Find the pivot
  for (std::ptrdiff_t i = 0; i < N; ++i)
  {
    std::ptrdiff_t max_row = i;
    for (std::ptrdiff_t k = i + 1; k < N; ++k)
    {
      if (std::abs(A[k][i]) > std::abs(A[max_row][i]))
      {
        max_row = k;
      }
    }

    // 1.b. Swap rows
    for (std::ptrdiff_t k = i; k < N; ++k)
    {
      std::swap(A[i][k], A[max_row][k]);
    }
    std::swap(b[i], b[max_row]);

    GEOS_ERROR_IF_LT_MSG( LvArray::math::abs(A[i][i]), LvArray::NumericLimits< real64 >::epsilon, "Singular matrix." );

    // 1.c Eliminate entries below the pivot
    for (std::ptrdiff_t k = i + 1; k < N; ++k)
    {
      real64 const scaling = A[k][i] / A[i][i];
      for (std::ptrdiff_t j = i; j < N; ++j)
      {
        A[k][j] -= scaling * A[i][j];
      }
      b[k] -= scaling * b[i];
    }
  }

  // Step 2: Backward substitution
  solveUpperTriangularSystem<N>( A, b, std::forward<N>(x) )
}

/**
 * Const version of the function
 */
template< std::ptrdiff_t N, 
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE 
inline
void solveGaussianElimination( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  real64[N][N] A_copy{};
  real64[N] b_copy{};

  for(std::ptrdiff_t i=0; i < N; ++j)
  {
    b_copy[i] = b[i];
    for( std::ptrdiff_t j=0; j < N; ++j )
    {
      A_copy[i][j] = A[i][j];
    };
  };
  
  solveGaussianElimination( A_copy, b_copy, std::forward<SOL_TYPE>(x) );
}

}; // internal namespace


/**
 * 
 */
template< std::ptrdiff_t N, 
          typename MATRIX_TYPE,
          typename RHS_TYPE,
          typename SOL_TYPE >
GEOS_HOST_DEVICE 
inline
void solve( MATRIX_TYPE const & A, RHS_TYPE const & b, SOL_TYPE && x )
{
  static_assert( N > 0, "N must be greater than 0." );
  internal::checkSizes< N, N >( A );
  internal::checkSizes< N >( b ); 
  internal::checkSizes< N >( x );
  
  if constexpr ( N == 2 )
  {
    internal::solveTwoByTwoSystem( A, b, std::forward<SOL_TYPE>(x) );
  }
  else if constexpr( N == 3 )
  {
    internal::solveThreeByThreeSystem( A, b, std::forward<SOL_TYPE>(x) );
  }
  else
  {
    internal::solveGaussianElimination< N >( A, b, std::forward<SOL_TYPE>(x) );
  }
}

};

};


#endif /*GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_*/
