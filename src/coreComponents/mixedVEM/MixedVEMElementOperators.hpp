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
 * @file MixedVEMElementOperators.hpp
 *
 * Element operators of the lowest-order mixed VEM for linear elasticity.
 *
 * With {phi_j}_{j=1}^{6 nf} the local virtual stress basis and {R_i}_{i=1}^{6} the
 * rigid body motions, the element contributes the saddle-point block
 *
 *   [ K_E  B_E^T ] [ sigma_E ]   [  g_E ]
 *   [ B_E   0    ] [   u_E   ] = [ -f_E ],
 *
 *   (B_E)_ij = (div phi_j, R_i)_E,
 *   (K_E)_ij = a_E(Pi_E phi_j, Pi_E phi_i) + s_E((I-Pi_E) phi_j, (I-Pi_E) phi_i).
 *
 * phi_j is virtual, but its face traction lies in T_h(f) and its divergence lies in
 * RM(E), and both operators are built from those two representations alone.
 *
 * All integrals are polynomial moments of the mesh entities and are evaluated in
 * closed form, so B_E and K_E are exact and no quadrature loop appears below. The
 * outputs are contiguous row-major blocks already expressed in the face-intrinsic
 * degree of freedom convention, so assembly is a plain index scatter.
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMELEMENTOPERATORS_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMELEMENTOPERATORS_HPP_

#include "mixedVEM/MixedVEMFaceBasis.hpp"

namespace geos
{

namespace mixedVEM
{

/**
 * @brief Build the divergence operator B_E.
 * @param[in] faceGeom the geometry of the @p numFaces faces of the element
 * @param[in] numFaces the number of faces
 * @param[in] elemCenter the point x_E
 * @param[out] divergence the 6 x (6 numFaces) matrix B_E
 *
 * (B_E)_ij = (div phi_j, R_i)_E = int_{dE} (phi_j n_E) . R_i df by parts, so with
 * R_i = e_i the rows carry the zeroth face moments and with R_{3+i} = e_i ^ (x - x_E)
 * they carry the rotational moments int_f (x - x_E) ^ phi_j.
 */
GEOS_HOST_DEVICE
inline void computeDivergenceOperator( FaceGeometry const * const faceGeom,
                                       integer const numFaces,
                                       real64 const (&elemCenter)[3],
                                       MatrixSlice const & divergence )
{
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    FaceGeometry const & geom = faceGeom[lf];

    real64 mean[NUM_FACE_DOF][3];
    real64 N[NUM_FACE_DOF][3][3];

    computeFaceBasisMeans( geom, mean );
    computeFaceBasisMoments( geom, elemCenter, N );

    for( integer j = 0; j < NUM_FACE_DOF; ++j )
    {
      integer const col = NUM_FACE_DOF * lf + j;

      real64 b[3];
      faceBasisRotationalMoment( N[j], b );

      // the outward sign turns the face-intrinsic traction into phi_j n_E
      for( integer i = 0; i < 3; ++i )
      {
        divergence( i, col ) = geom.outwardSign * mean[j][i];
        divergence( 3 + i, col ) = geom.outwardSign * b[i];
      }
    }
  }
}

/**
 * @brief Build the divergence reconstruction D_E, the RM(E) coefficients of div phi_j.
 * @param[in] divergence the matrix B_E
 * @param[in] numFaces the number of faces
 * @param[in] moments the element moments
 * @param[out] divReconstruction the 6 x (6 numFaces) matrix D_E
 *
 * div phi_j = alpha_j + omega_j ^ (x - x_E) lies in RM(E), so B_E = W D_E with W the
 * Gram matrix W_ik = (R_i, R_k)_E of the rigid body motions,
 *
 *   W = [ |E| I   C  ],   C = -[m1]_x,   A_E = tr(M) I - M,
 *       [  C^T   A_E ]
 *
 * with m1 and M the first and second element moments. Eliminating alpha leaves the
 * 3x3 Schur complement S = A_E - (|m1|^2 I - m1 @ m1)/|E|, which reduces to the
 * system (8) of Proposition 3.1 when x_E is the exact barycenter and m1 vanishes.
 * Keeping C makes the reconstruction exact for any choice of x_E.
 */
GEOS_HOST_DEVICE
inline void computeDivergenceReconstruction( MatrixSliceConst const & divergence,
                                             integer const numFaces,
                                             ElementMoments const & moments,
                                             MatrixSlice const & divReconstruction )
{
  real64 inertia[3][3];
  computeInertia( moments, inertia );

  real64 const (&m1)[3] = moments.firstMoment;
  real64 const invVolume = 1.0 / moments.volume;
  real64 const m1Squared = LvArray::tensorOps::AiBi< 3 >( m1, m1 );

  real64 schur[3][3];
  for( integer i = 0; i < 3; ++i )
  {
    for( integer k = 0; k < 3; ++k )
    {
      schur[i][k] = inertia[i][k] - invVolume * ( ( ( i == k ) ? m1Squared : 0.0 ) - m1[i] * m1[k] );
    }
  }

  real64 inverseSchur[3][3];
  LvArray::tensorOps::invert< 3 >( inverseSchur, schur );

  integer const numStressDof = NUM_FACE_DOF * numFaces;

  for( integer j = 0; j < numStressDof; ++j )
  {
    real64 b1[3], b2[3];
    for( integer i = 0; i < 3; ++i )
    {
      b1[i] = divergence( i, j );
      b2[i] = divergence( 3 + i, j );
    }

    // omega solves S omega = b2 - m1 ^ b1 / |E|
    real64 rhs[3];
    LvArray::tensorOps::crossProduct( rhs, m1, b1 );
    for( integer i = 0; i < 3; ++i )
    {
      rhs[i] = b2[i] - invVolume * rhs[i];
    }

    real64 omega[3] = { 0.0, 0.0, 0.0 };
    for( integer i = 0; i < 3; ++i )
    {
      for( integer k = 0; k < 3; ++k )
      {
        omega[i] += inverseSchur[i][k] * rhs[k];
      }
    }

    // alpha = ( b1 + m1 ^ omega ) / |E|
    real64 alpha[3];
    LvArray::tensorOps::crossProduct( alpha, m1, omega );
    for( integer i = 0; i < 3; ++i )
    {
      divReconstruction( i, j ) = invVolume * ( b1[i] + alpha[i] );
      divReconstruction( 3 + i, j ) = omega[i];
    }
  }
}

/**
 * @brief Moments of the rigid body motions against the projection basis.
 * @param[in] moments the element moments
 * @param[out] rmMoments the 6x6 matrix G with G_ai = int_E p_a . R_i dE
 *
 * p_a is the linear field with eps(p_a) = pi_a, and because pi_a is symmetric
 * p_a(x) = pi_a (x - x_E). Hence G_ai = pi_a : int_E R_i @ (x - x_E) dE, which the
 * first and second element moments give exactly.
 */
GEOS_HOST_DEVICE
inline void computeRigidMotionMoments( ElementMoments const & moments,
                                       real64 (& rmMoments)[NUM_SYM_COMP][NUM_RM_DOF] )
{
  real64 const (&M)[3][3] = moments.secondMoment;
  real64 const (&m1)[3] = moments.firstMoment;

  real64 Z[NUM_RM_DOF][3][3];

  // translations: int_E e_k @ (x - x_E) dE, which vanishes when x_E is the barycenter
  for( integer k = 0; k < 3; ++k )
  {
    for( integer p = 0; p < 3; ++p )
    {
      for( integer q = 0; q < 3; ++q )
      {
        Z[k][p][q] = ( p == k ) ? m1[q] : 0.0;
      }
    }
  }

  // rotations: int_E (e_k ^ (x - x_E)) @ (x - x_E) dE = [e_k]_x M
  for( integer q = 0; q < 3; ++q )
  {
    Z[3][0][q] = 0.0;      Z[3][1][q] = -M[2][q]; Z[3][2][q] = M[1][q];
    Z[4][0][q] = M[2][q];  Z[4][1][q] = 0.0;      Z[4][2][q] = -M[0][q];
    Z[5][0][q] = -M[1][q]; Z[5][1][q] = M[0][q];  Z[5][2][q] = 0.0;
  }

  for( integer i = 0; i < NUM_RM_DOF; ++i )
  {
    real64 column[NUM_SYM_COMP];
    projectOnSymBasis( Z[i], column );

    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      rmMoments[a][i] = column[a];
    }
  }
}

/**
 * @brief Build the constant stress projection P_E of Pi_E.
 * @param[in] faceGeom the geometry of the faces of the element
 * @param[in] numFaces the number of faces
 * @param[in] elemCenter the point x_E
 * @param[in] moments the element moments
 * @param[in] divReconstruction the matrix D_E
 * @param[out] projection the 6 x (6 numFaces) matrix P_E
 *
 * Integrating (13) by parts gives
 *   |E| (Pi_E phi_j)_a = -int_E div phi_j . p_a dE + int_f (phi_j n_E) . p_a df,
 * the first term being -G D_E because div phi_j expands on RM(E), and the second
 * being the symmetric part of the face moment N_j because p_a(x) = pi_a (x - x_E).
 */
GEOS_HOST_DEVICE
inline void computeProjectionOperator( FaceGeometry const * const faceGeom,
                                       integer const numFaces,
                                       real64 const (&elemCenter)[3],
                                       ElementMoments const & moments,
                                       MatrixSliceConst const & divReconstruction,
                                       MatrixSlice const & projection )
{
  real64 rmMoments[NUM_SYM_COMP][NUM_RM_DOF];
  computeRigidMotionMoments( moments, rmMoments );

  real64 const invVolume = 1.0 / moments.volume;
  integer const numStressDof = NUM_FACE_DOF * numFaces;

  // volume term: -G D_E
  for( integer a = 0; a < NUM_SYM_COMP; ++a )
  {
    for( integer j = 0; j < numStressDof; ++j )
    {
      real64 value = 0.0;
      for( integer i = 0; i < NUM_RM_DOF; ++i )
      {
        value += rmMoments[a][i] * divReconstruction( i, j );
      }
      projection( a, j ) = -value;
    }
  }

  // face term: block diagonal, one 6x6 per face
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    FaceGeometry const & geom = faceGeom[lf];

    real64 N[NUM_FACE_DOF][3][3];
    computeFaceBasisMoments( geom, elemCenter, N );

    for( integer j = 0; j < NUM_FACE_DOF; ++j )
    {
      real64 column[NUM_SYM_COMP];
      projectOnSymBasis( N[j], column );

      integer const col = NUM_FACE_DOF * lf + j;
      for( integer a = 0; a < NUM_SYM_COMP; ++a )
      {
        projection( a, col ) += geom.outwardSign * column[a];
      }
    }
  }

  // M_E = |E| I because the projection basis is orthonormal for ":"
  for( integer a = 0; a < NUM_SYM_COMP; ++a )
  {
    for( integer j = 0; j < numStressDof; ++j )
    {
      projection( a, j ) *= invVolume;
    }
  }
}

/**
 * @brief Traction map Lambda_n with (Lambda_n)_{pa} = (pi_a n)_p.
 * @param[in] normal the direction n
 * @param[out] tractionMap the 3x6 matrix
 */
GEOS_HOST_DEVICE
inline void computeTractionMap( real64 const (&normal)[3],
                                real64 (& tractionMap)[3][NUM_SYM_COMP] )
{
  for( integer a = 0; a < NUM_SYM_COMP; ++a )
  {
    real64 coefficients[NUM_SYM_COMP] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    coefficients[a] = 1.0;

    real64 traction[3];
    symBasisTraction( coefficients, normal, traction );

    for( integer p = 0; p < 3; ++p )
    {
      tractionMap[p][a] = traction[p];
    }
  }
}

/**
 * @brief Build the element stiffness K_E, consistency plus stabilization.
 * @param[in] faceGeom the geometry of the faces of the element
 * @param[in] numFaces the number of faces
 * @param[in] volume the element volume |E|
 * @param[in] diameter the element diameter h_E
 * @param[in] compliance the 6x6 matrix of D = C^{-1}
 * @param[in] projection the matrix P_E
 * @param[in,out] workspace a 6 x (6 numFaces) scratch block
 * @param[out] stiffness the (6 numFaces) x (6 numFaces) matrix K_E
 *
 * a_E^h(sigma,tau) = a_E(Pi_E sigma, Pi_E tau) + s_E((I-Pi_E) sigma, (I-Pi_E) tau)
 * with the boundary stabilization s_E(sigma,tau) = kappa_E h_E int_{dE} (sigma n).(tau n).
 * Expanding the residual (I-Pi_E) phi_j n on a face into the projected traction minus
 * phi_j collects every quadratic term into a single 6x6 weight,
 *   W = |E| D + kappa_E h_E sum_f |f| Lambda_{n_f}^T Lambda_{n_f},
 * so the dense part of K_E is the one product P_E^T W P_E; the remaining cross and
 * Gram terms only touch the six columns of each face. The outward sign appears twice
 * in the stabilization and cancels.
 */
GEOS_HOST_DEVICE
inline void computeStiffness( FaceGeometry const * const faceGeom,
                              integer const numFaces,
                              real64 const volume,
                              real64 const diameter,
                              real64 const (&compliance)[NUM_SYM_COMP][NUM_SYM_COMP],
                              MatrixSliceConst const & projection,
                              MatrixSlice const & workspace,
                              MatrixSlice const & stiffness )
{
  integer const numStressDof = NUM_FACE_DOF * numFaces;

  // kappa_E is the deviatoric compliance 1/(2 mu). The paper takes half the trace of D,
  // equation (15); measured on a tetrahedral column this choice is slightly the better of
  // the two, and it carries no lambda, so the stabilization does not drift as the material
  // approaches incompressibility
  real64 const stabScale = compliance[3][3] * diameter;

  real64 weight[NUM_SYM_COMP][NUM_SYM_COMP];
  for( integer a = 0; a < NUM_SYM_COMP; ++a )
  {
    for( integer b = 0; b < NUM_SYM_COMP; ++b )
    {
      weight[a][b] = volume * compliance[a][b];
    }
  }

  for( integer lf = 0; lf < numFaces; ++lf )
  {
    FaceGeometry const & geom = faceGeom[lf];

    real64 tractionMap[3][NUM_SYM_COMP];
    computeTractionMap( geom.normal, tractionMap );

    real64 const w = stabScale * geom.area;
    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      for( integer b = 0; b < NUM_SYM_COMP; ++b )
      {
        real64 value = 0.0;
        for( integer p = 0; p < 3; ++p )
        {
          value += tractionMap[p][a] * tractionMap[p][b];
        }
        weight[a][b] += w * value;
      }
    }
  }

  // workspace = W P_E, then K_E = P_E^T workspace with a contiguous inner loop
  for( integer a = 0; a < NUM_SYM_COMP; ++a )
  {
    for( integer j = 0; j < numStressDof; ++j )
    {
      real64 value = 0.0;
      for( integer b = 0; b < NUM_SYM_COMP; ++b )
      {
        value += weight[a][b] * projection( b, j );
      }
      workspace( a, j ) = value;
    }
  }

  for( integer i = 0; i < numStressDof; ++i )
  {
    real64 column[NUM_SYM_COMP];
    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      column[a] = projection( a, i );
    }

    for( integer j = 0; j < numStressDof; ++j )
    {
      stiffness( i, j ) = 0.0;
    }

    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      real64 const c = column[a];
      for( integer j = 0; j < numStressDof; ++j )
      {
        stiffness( i, j ) += c * workspace( a, j );
      }
    }
  }

  // cross and Gram terms of the stabilization, confined to the columns of each face
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    FaceGeometry const & geom = faceGeom[lf];

    real64 tractionMap[3][NUM_SYM_COMP];
    computeTractionMap( geom.normal, tractionMap );

    real64 mean[NUM_FACE_DOF][3];
    computeFaceBasisMeans( geom, mean );

    real64 gram[NUM_FACE_DOF][NUM_FACE_DOF];
    computeFaceBasisGram( geom, gram );

    // cross weight Lambda_n^T int_f phi_k df
    real64 cross[NUM_SYM_COMP][NUM_FACE_DOF];
    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      for( integer k = 0; k < NUM_FACE_DOF; ++k )
      {
        real64 value = 0.0;
        for( integer p = 0; p < 3; ++p )
        {
          value += tractionMap[p][a] * mean[k][p];
        }
        cross[a][k] = stabScale * value;
      }
    }

    integer const offset = NUM_FACE_DOF * lf;

    for( integer i = 0; i < numStressDof; ++i )
    {
      for( integer k = 0; k < NUM_FACE_DOF; ++k )
      {
        real64 value = 0.0;
        for( integer a = 0; a < NUM_SYM_COMP; ++a )
        {
          value += projection( a, i ) * cross[a][k];
        }
        stiffness( i, offset + k ) -= value;
        stiffness( offset + k, i ) -= value;
      }
    }

    for( integer k = 0; k < NUM_FACE_DOF; ++k )
    {
      for( integer l = 0; l < NUM_FACE_DOF; ++l )
      {
        stiffness( offset + k, offset + l ) += stabScale * gram[k][l];
      }
    }
  }
}

/**
 * @brief Build every element operator of the mixed VEM in one pass.
 * @param[in] faceGeom the geometry of the faces of the element, already oriented
 * @param[in] numFaces the number of faces
 * @param[in] elemCenter the point x_E
 * @param[in] diameter the element diameter h_E
 * @param[in] moments the element moments
 * @param[in] compliance the 6x6 matrix of D = C^{-1}
 * @param[out] divergence the 6 x (6 numFaces) matrix B_E
 * @param[out] divReconstruction the 6 x (6 numFaces) matrix D_E
 * @param[out] projection the 6 x (6 numFaces) matrix P_E
 * @param[in,out] workspace a 6 x (6 numFaces) scratch block, reusable across elements
 * @param[out] stiffness the (6 numFaces) x (6 numFaces) matrix K_E
 */
GEOS_HOST_DEVICE
inline void computeElementOperators( FaceGeometry const * const faceGeom,
                                     integer const numFaces,
                                     real64 const (&elemCenter)[3],
                                     real64 const diameter,
                                     ElementMoments const & moments,
                                     real64 const (&compliance)[NUM_SYM_COMP][NUM_SYM_COMP],
                                     MatrixSlice const & divergence,
                                     MatrixSlice const & divReconstruction,
                                     MatrixSlice const & projection,
                                     MatrixSlice const & workspace,
                                     MatrixSlice const & stiffness )
{
  computeDivergenceOperator( faceGeom, numFaces, elemCenter, divergence );

  computeDivergenceReconstruction( divergence.toSliceConst(), numFaces, moments, divReconstruction );

  computeProjectionOperator( faceGeom, numFaces, elemCenter, moments,
                             divReconstruction.toSliceConst(), projection );

  computeStiffness( faceGeom, numFaces, moments.volume, diameter, compliance,
                    projection.toSliceConst(), workspace, stiffness );
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMELEMENTOPERATORS_HPP_
