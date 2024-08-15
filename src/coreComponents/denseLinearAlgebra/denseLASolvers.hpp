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

    real64 const A_inv[2][2] = { { A[1][1] / detA, -A[0][1] / detA },
                                 { -A[1][0] / detA, A[0][0] / detA } };

    x[0] = A_inv[0][0] * b[0] + A_inv[0][1] * b[1];
    x[1] = A_inv[1][0] * b[0] + A_inv[1][1] * b[1];
}

};

};


#endif /*GEOS_DENSELINEARALGEBRA_DENSELASOLVERS_HPP_*/
