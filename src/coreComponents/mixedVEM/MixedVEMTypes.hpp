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
 * @file MixedVEMTypes.hpp
 *
 * Sizes and symmetric-tensor algebra shared by the lowest-order mixed VEM for
 * linear elasticity (Hellinger-Reissner), following Dassi, Lovadina & Visinoni,
 * CMAME 364 (2020) 112910.
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMTYPES_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMTYPES_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "denseLinearAlgebra/common/layouts.hpp"

namespace geos
{

namespace mixedVEM
{

/// dim T_h(f) = 6: two constant tangential, one in-plane rotational, three linear normal modes.
constexpr integer NUM_FACE_DOF = 6;

/// dim RM(E) = 6: three translations and three infinitesimal rotations.
constexpr integer NUM_RM_DOF = 6;

/// dim [P_0(E)]_s^{3x3} = 6: the constant symmetric stresses reproduced by Pi_E.
constexpr integer NUM_SYM_COMP = 6;

/// Off-diagonal scaling that makes the symmetric tensor basis orthonormal for ":".
constexpr real64 INV_SQRT_2 = 0.7071067811865475244;

/// Row-major dense block, the layout every element operator is written into.
using MatrixSlice = arraySlice2d< real64, MatrixLayout::ROW_MAJOR >;

/// Read-only row-major dense block.
using MatrixSliceConst = arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR >;

/**
 * @brief Coefficients of a tensor on the orthonormal symmetric basis, v_a = pi_a : A.
 * @param[in] A the second-order tensor (only its symmetric part contributes)
 * @param[out] v the six coefficients (xx, yy, zz, yz, xz, xy) with shears scaled by sqrt(2)
 *
 * The basis is pi_0 = e1 x e1, pi_1 = e2 x e2, pi_2 = e3 x e3,
 * pi_3 = (e2 x e3 + e3 x e2)/sqrt(2), pi_4 = (e1 x e3 + e3 x e1)/sqrt(2),
 * pi_5 = (e1 x e2 + e2 x e1)/sqrt(2), so that pi_a : pi_b = delta_ab.
 */
GEOS_HOST_DEVICE
inline void projectOnSymBasis( real64 const (&A)[3][3],
                               real64 (& v)[NUM_SYM_COMP] )
{
  v[0] = A[0][0];
  v[1] = A[1][1];
  v[2] = A[2][2];
  v[3] = INV_SQRT_2 * ( A[1][2] + A[2][1] );
  v[4] = INV_SQRT_2 * ( A[0][2] + A[2][0] );
  v[5] = INV_SQRT_2 * ( A[0][1] + A[1][0] );
}

/**
 * @brief Rebuild the symmetric tensor A = sum_a v_a pi_a from its coefficients.
 * @param[in] v the six coefficients
 * @param[out] A the symmetric tensor
 */
GEOS_HOST_DEVICE
inline void expandSymBasis( real64 const (&v)[NUM_SYM_COMP],
                            real64 (& A)[3][3] )
{
  real64 const s3 = INV_SQRT_2 * v[3];
  real64 const s4 = INV_SQRT_2 * v[4];
  real64 const s5 = INV_SQRT_2 * v[5];

  A[0][0] = v[0]; A[0][1] = s5;   A[0][2] = s4;
  A[1][0] = s5;   A[1][1] = v[1]; A[1][2] = s3;
  A[2][0] = s4;   A[2][1] = s3;   A[2][2] = v[2];
}

/**
 * @brief Traction (sum_a v_a pi_a) n of a symmetric tensor given by its coefficients.
 * @param[in] v the six coefficients
 * @param[in] n the direction
 * @param[out] t the resulting traction vector
 */
GEOS_HOST_DEVICE
inline void symBasisTraction( real64 const (&v)[NUM_SYM_COMP],
                              real64 const (&n)[3],
                              real64 (& t)[3] )
{
  real64 const s3 = INV_SQRT_2 * v[3];
  real64 const s4 = INV_SQRT_2 * v[4];
  real64 const s5 = INV_SQRT_2 * v[5];

  t[0] = v[0] * n[0] + s5 * n[1] + s4 * n[2];
  t[1] = s5 * n[0] + v[1] * n[1] + s3 * n[2];
  t[2] = s4 * n[0] + s3 * n[1] + v[2] * n[2];
}

/**
 * @brief Isotropic compliance D = C^{-1} on the orthonormal symmetric basis.
 * @param[in] lambda the first Lame parameter
 * @param[in] mu the shear modulus
 * @param[out] compliance the 6x6 matrix of a(sigma,tau) = (D sigma, tau)
 *
 * D sigma = sigma/(2 mu) - lambda tr(sigma) I / (2 mu (3 lambda + 2 mu)).
 */
GEOS_HOST_DEVICE
inline void makeIsotropicCompliance( real64 const lambda,
                                     real64 const mu,
                                     real64 (& compliance)[NUM_SYM_COMP][NUM_SYM_COMP] )
{
  real64 const den = 2.0 * mu * ( 3.0 * lambda + 2.0 * mu );

  real64 const a = ( lambda + mu ) / ( mu * ( 3.0 * lambda + 2.0 * mu ) );
  real64 const b = -lambda / den;
  real64 const c = 0.5 / mu;

  for( integer i = 0; i < NUM_SYM_COMP; ++i )
  {
    for( integer j = 0; j < NUM_SYM_COMP; ++j )
    {
      compliance[i][j] = 0.0;
    }
  }

  // normal block: the shear entries carry no lambda because the basis is orthonormal
  compliance[0][0] = a; compliance[0][1] = b; compliance[0][2] = b;
  compliance[1][0] = b; compliance[1][1] = a; compliance[1][2] = b;
  compliance[2][0] = b; compliance[2][1] = b; compliance[2][2] = a;

  compliance[3][3] = c;
  compliance[4][4] = c;
  compliance[5][5] = c;
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMTYPES_HPP_
