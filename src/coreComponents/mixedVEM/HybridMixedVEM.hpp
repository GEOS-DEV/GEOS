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
 * @file HybridMixedVEM.hpp
 *
 * Static condensation of the element saddle point block onto the interface multiplier.
 *
 * Breaking the stress space elementwise and restoring traction continuity with a
 * multiplier lambda on the mesh skeleton turns the global problem into
 *
 *   [ K^  B^T  -C^T ] [ sigma^ ]   [  g ]
 *   [ B^   0     0  ] [   u    ] = [ -f ]
 *   [ C    0     0  ] [ lambda ]   [  0 ].
 *
 * The element unknowns are local to E, so they are eliminated through
 *
 *   M_E = [ K_E  B_E^T ],   S_E = ( M_E^{-1} )_{sigma sigma},
 *         [ B_E   0    ]
 *
 * leaving the interface problem H lambda = h with H = sum_E C_E S_E C_E^T. Because
 * K_E is symmetric positive definite and B_E has full row rank,
 *
 *   S_E = K_E^{-1} - W G^{-1} W^T,   W = K_E^{-1} B_E^T,   G = B_E K_E^{-1} B_E^T,
 *
 * is symmetric positive semidefinite with the six dimensional kernel range(W), so H is
 * symmetric positive semidefinite and becomes definite once the Dirichlet multipliers
 * are eliminated. That is the whole point of the option: a conjugate gradient solve on
 * 6 N_f^int unknowns replaces an indefinite solve on 6 N_f + 6 N_c.
 *
 * With the degree of freedom convention of this package the traction of E on f is
 * s_{E,f} sum_j sigma_j phi_j, and lambda_j = int_f u . phi_j df is the matching moment
 * of the displacement trace, so the continuity operator is simply C_{E,f} = s_{E,f} I_6.
 */

#ifndef GEOS_MIXEDVEM_HYBRIDMIXEDVEM_HPP_
#define GEOS_MIXEDVEM_HYBRIDMIXEDVEM_HPP_

#include "mixedVEM/MixedVEMElementOperators.hpp"

namespace geos
{

namespace mixedVEM
{

/**
 * @brief In place Cholesky factorization A = L L^T of a symmetric positive definite matrix.
 * @param[in,out] matrix the matrix, overwritten by the lower factor L
 * @param[in] n the order of the matrix
 * @return false if a non positive pivot is met, in which case @p matrix is left partial
 */
GEOS_HOST_DEVICE
inline bool choleskyFactorize( MatrixSlice const & matrix,
                               integer const n )
{
  for( integer i = 0; i < n; ++i )
  {
    for( integer j = 0; j <= i; ++j )
    {
      real64 sum = matrix( i, j );
      for( integer k = 0; k < j; ++k )
      {
        sum -= matrix( i, k ) * matrix( j, k );
      }

      if( i == j )
      {
        if( !( sum > 0.0 ) )
        {
          return false;
        }
        matrix( i, i ) = LvArray::math::sqrt( sum );
      }
      else
      {
        matrix( i, j ) = sum / matrix( j, j );
      }
    }

    for( integer j = i + 1; j < n; ++j )
    {
      matrix( i, j ) = 0.0;
    }
  }

  return true;
}

/**
 * @brief Solve L L^T x = b in place given the lower Cholesky factor.
 * @param[in] factor the lower factor L
 * @param[in] n the order of the matrix
 * @param[in,out] x the right hand side, overwritten by the solution
 */
GEOS_HOST_DEVICE
inline void choleskySolve( MatrixSliceConst const & factor,
                           integer const n,
                           real64 * const x )
{
  for( integer i = 0; i < n; ++i )
  {
    real64 sum = x[i];
    for( integer k = 0; k < i; ++k )
    {
      sum -= factor( i, k ) * x[k];
    }
    x[i] = sum / factor( i, i );
  }

  for( integer i = n - 1; i >= 0; --i )
  {
    real64 sum = x[i];
    for( integer k = i + 1; k < n; ++k )
    {
      sum -= factor( k, i ) * x[k];
    }
    x[i] = sum / factor( i, i );
  }
}

/**
 * @brief Invert a small symmetric positive definite matrix in place.
 * @param[in,out] matrix the matrix, overwritten by its inverse
 * @return false if the matrix is not positive definite
 */
GEOS_HOST_DEVICE
inline bool invertRigidMotionGram( real64 (& matrix)[NUM_RM_DOF][NUM_RM_DOF] )
{
  real64 factor[NUM_RM_DOF][NUM_RM_DOF] = { { 0.0 } };

  for( integer i = 0; i < NUM_RM_DOF; ++i )
  {
    for( integer j = 0; j <= i; ++j )
    {
      real64 sum = matrix[i][j];
      for( integer k = 0; k < j; ++k )
      {
        sum -= factor[i][k] * factor[j][k];
      }

      if( i == j )
      {
        if( !( sum > 0.0 ) )
        {
          return false;
        }
        factor[i][i] = LvArray::math::sqrt( sum );
      }
      else
      {
        factor[i][j] = sum / factor[j][j];
      }
    }
  }

  // inverse of the lower factor, then A^{-1} = L^{-T} L^{-1}
  real64 inverseFactor[NUM_RM_DOF][NUM_RM_DOF] = { { 0.0 } };
  for( integer i = 0; i < NUM_RM_DOF; ++i )
  {
    inverseFactor[i][i] = 1.0 / factor[i][i];
    for( integer j = 0; j < i; ++j )
    {
      real64 sum = 0.0;
      for( integer k = j; k < i; ++k )
      {
        sum += factor[i][k] * inverseFactor[k][j];
      }
      inverseFactor[i][j] = -sum / factor[i][i];
    }
  }

  for( integer i = 0; i < NUM_RM_DOF; ++i )
  {
    for( integer j = 0; j <= i; ++j )
    {
      real64 sum = 0.0;
      for( integer k = i; k < NUM_RM_DOF; ++k )
      {
        sum += inverseFactor[k][i] * inverseFactor[k][j];
      }
      matrix[i][j] = sum;
      matrix[j][i] = sum;
    }
  }

  return true;
}

/**
 * @brief Eliminate the element unknowns and form the local Schur complement S_E.
 * @param[in] stiffness the matrix K_E
 * @param[in] divergence the matrix B_E
 * @param[in] numFaces the number of faces of the element
 * @param[in,out] factorization a (6 numFaces) x (6 numFaces) scratch block
 * @param[out] couplingTranspose W^T = B_E K_E^{-1}, of size 6 x (6 numFaces)
 * @param[out] inverseDivGram the inverse of the divergence Gram matrix G = B_E K_E^{-1} B_E^T
 * @param[out] schur the (6 numFaces) x (6 numFaces) matrix S_E
 * @return false if K_E or G is not positive definite
 *
 * S_E = K_E^{-1} - W G^{-1} W^T is the stress-stress block of M_E^{-1}. It is symmetric
 * positive semidefinite and annihilates range(W), the six dimensional space conjugate to
 * the rigid body motions. Every triangular solve runs along a contiguous row, using the
 * symmetry of K_E^{-1} to write columns as rows.
 */
GEOS_HOST_DEVICE
inline bool computeLocalCondensation( MatrixSliceConst const & stiffness,
                                      MatrixSliceConst const & divergence,
                                      integer const numFaces,
                                      MatrixSlice const & factorization,
                                      MatrixSlice const & couplingTranspose,
                                      real64 (& inverseDivGram)[NUM_RM_DOF][NUM_RM_DOF],
                                      MatrixSlice const & schur )
{
  integer const numStressDof = NUM_FACE_DOF * numFaces;

  for( integer i = 0; i < numStressDof; ++i )
  {
    for( integer j = 0; j < numStressDof; ++j )
    {
      factorization( i, j ) = stiffness( i, j );
    }
  }

  if( !choleskyFactorize( factorization, numStressDof ) )
  {
    return false;
  }

  MatrixSliceConst const factor = factorization.toSliceConst();

  // K_E^{-1} is symmetric, so solving for column j and storing it in row j is the same
  for( integer j = 0; j < numStressDof; ++j )
  {
    for( integer i = 0; i < numStressDof; ++i )
    {
      schur( j, i ) = ( i == j ) ? 1.0 : 0.0;
    }
    choleskySolve( factor, numStressDof, &schur( j, 0 ) );
  }

  // W^T = B_E K_E^{-1}, one row per rigid body motion
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    for( integer i = 0; i < numStressDof; ++i )
    {
      couplingTranspose( k, i ) = divergence( k, i );
    }
    choleskySolve( factor, numStressDof, &couplingTranspose( k, 0 ) );
  }

  // G = B_E W, inverted in place
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    for( integer l = 0; l < NUM_RM_DOF; ++l )
    {
      real64 sum = 0.0;
      for( integer i = 0; i < numStressDof; ++i )
      {
        sum += divergence( k, i ) * couplingTranspose( l, i );
      }
      inverseDivGram[k][l] = sum;
    }
  }

  if( !invertRigidMotionGram( inverseDivGram ) )
  {
    return false;
  }

  // S_E = K_E^{-1} - W G^{-1} W^T, written into the lower triangle and mirrored
  for( integer i = 0; i < numStressDof; ++i )
  {
    real64 row[NUM_RM_DOF];
    for( integer k = 0; k < NUM_RM_DOF; ++k )
    {
      real64 sum = 0.0;
      for( integer l = 0; l < NUM_RM_DOF; ++l )
      {
        sum += inverseDivGram[k][l] * couplingTranspose( l, i );
      }
      row[k] = sum;
    }

    for( integer j = 0; j <= i; ++j )
    {
      real64 correction = 0.0;
      for( integer k = 0; k < NUM_RM_DOF; ++k )
      {
        correction += row[k] * couplingTranspose( k, j );
      }

      real64 const value = 0.5 * ( schur( i, j ) + schur( j, i ) ) - correction;
      schur( i, j ) = value;
      schur( j, i ) = value;
    }
  }

  return true;
}

/**
 * @brief Apply the continuity operator, S_E -> C_E S_E C_E^T.
 * @param[in] faceGeom the geometry of the faces of the element
 * @param[in] numFaces the number of faces of the element
 * @param[in,out] schur the matrix S_E, overwritten by the interface contribution
 *
 * C_{E,f} = s_{E,f} I_6, so the operator is the symmetric diagonal scaling by the
 * outward signs and the result can be scattered row by row without repacking.
 */
GEOS_HOST_DEVICE
inline void applyContinuityOperator( FaceGeometry const * const faceGeom,
                                     integer const numFaces,
                                     MatrixSlice const & schur )
{
  integer const numStressDof = NUM_FACE_DOF * numFaces;

  for( integer i = 0; i < numStressDof; ++i )
  {
    real64 const rowSign = faceGeom[i / NUM_FACE_DOF].outwardSign;

    for( integer j = 0; j < numStressDof; ++j )
    {
      schur( i, j ) *= rowSign * faceGeom[j / NUM_FACE_DOF].outwardSign;
    }
  }
}

/**
 * @brief Recover the element stress and displacement from the interface multiplier.
 * @param[in] schur the matrix S_E, before the continuity scaling
 * @param[in] couplingTranspose the block W^T = B_E K_E^{-1}
 * @param[in] inverseDivGram the inverse of the divergence Gram matrix G
 * @param[in] faceGeom the geometry of the faces of the element
 * @param[in] numFaces the number of faces of the element
 * @param[in] multiplier the 6 numFaces values of lambda on the faces of the element
 * @param[in] stressRhs the element stress right hand side g_E
 * @param[in] bodyForce the six components of the element load f_E
 * @param[out] stress the 6 numFaces stress degrees of freedom sigma_E
 * @param[out] displacement the six rigid body motion coefficients u_E
 *
 * [ sigma_E ; u_E ] = M_E^{-1} [ g_E + C_E^T lambda ; -f_E ], which is elementwise
 * independent, so the recovery is a perfectly parallel pass over the mesh.
 */
GEOS_HOST_DEVICE
inline void recoverElementSolution( MatrixSliceConst const & schur,
                                    MatrixSliceConst const & couplingTranspose,
                                    real64 const (&inverseDivGram)[NUM_RM_DOF][NUM_RM_DOF],
                                    FaceGeometry const * const faceGeom,
                                    integer const numFaces,
                                    real64 const * const multiplier,
                                    real64 const * const stressRhs,
                                    real64 const * const bodyForce,
                                    real64 * const stress,
                                    real64 * const displacement )
{
  integer const numStressDof = NUM_FACE_DOF * numFaces;

  // u_E = G^{-1} ( W^T ( g_E + C_E^T lambda ) + f_E )
  real64 reduced[NUM_RM_DOF];
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    real64 sum = bodyForce[k];
    for( integer i = 0; i < numStressDof; ++i )
    {
      real64 const rhs = stressRhs[i] + faceGeom[i / NUM_FACE_DOF].outwardSign * multiplier[i];
      sum += couplingTranspose( k, i ) * rhs;
    }
    reduced[k] = sum;
  }

  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    real64 sum = 0.0;
    for( integer l = 0; l < NUM_RM_DOF; ++l )
    {
      sum += inverseDivGram[k][l] * reduced[l];
    }
    displacement[k] = sum;
  }

  // sigma_E = S_E ( g_E + C_E^T lambda ) - W G^{-1} f_E
  real64 scaledForce[NUM_RM_DOF];
  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    real64 sum = 0.0;
    for( integer l = 0; l < NUM_RM_DOF; ++l )
    {
      sum += inverseDivGram[k][l] * bodyForce[l];
    }
    scaledForce[k] = sum;
  }

  for( integer i = 0; i < numStressDof; ++i )
  {
    real64 sum = 0.0;
    for( integer j = 0; j < numStressDof; ++j )
    {
      sum += schur( i, j ) * ( stressRhs[j] + faceGeom[j / NUM_FACE_DOF].outwardSign * multiplier[j] );
    }

    for( integer k = 0; k < NUM_RM_DOF; ++k )
    {
      sum -= couplingTranspose( k, i ) * scaledForce[k];
    }

    stress[i] = sum;
  }
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_HYBRIDMIXEDVEM_HPP_
