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
 * Element operators of Dassi, Lovadina, Visinoni (2020). With tau_j the basis of Sigma_h(E), dual to
 * the face moments (5)-(6), and r_i the basis of RM(E), equation (2),
 *
 *   (B_E)_ij = (div tau_j, r_i)_E,   (M_E)_ij = a_E^h(tau_j, tau_i), equation (14).
 *
 * Every integral is a closed-form moment, so B_E and M_E are exact.
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
 * (B_E)_ij = int_{dE} (tau_j n) . r_i df, equation (9): zeroth face moments for r_i = e_i,
 * rotational moments int_f (x - x_E) ^ psi_k for r_{3+i} = e_i ^ (x - x_E).
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

      // s_{E,f} turns the face unknown into the outward traction tau_j n
      for( integer i = 0; i < 3; ++i )
      {
        divergence( i, col ) = geom.outwardSign * mean[j][i];
        divergence( 3 + i, col ) = geom.outwardSign * b[i];
      }
    }
  }
}

/**
 * @brief Coefficients (alpha_E, omega_E) of div tau_j = alpha_E + omega_E ^ (x - x_E), Proposition 3.1.
 * @param[in] divergence the matrix B_E
 * @param[in] numFaces the number of faces
 * @param[in] moments the element moments
 * @param[out] divReconstruction the 6 x (6 n_f^E) matrix, rows (alpha_E, omega_E)
 *
 * B_E = W_RM (alpha_E, omega_E) with W_RM = (r_i, r_k)_E. The first moment m1 of x_E is kept,
 * so this reduces to (7)-(8) when x_E is the barycenter and holds for any x_E otherwise.
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
  real64 const m1Squared = LvArray::tensorOps::l2NormSquared< 3 >( m1 );

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
 * @brief G_ai = int_E p_a . r_i dE, with p_a(x) = pi_a (x - x_E) so that eps(p_a) = pi_a.
 * @param[in] moments the element moments
 * @param[out] rmMoments the 6x6 matrix G
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
 * @brief Matrix P_E of Pi_E, equation (13).
 * @param[in] faceGeom the geometry of the faces of the element
 * @param[in] numFaces the number of faces
 * @param[in] elemCenter the point x_E
 * @param[in] moments the element moments
 * @param[in] divReconstruction the coefficients (alpha_E, omega_E)
 * @param[out] projection the 6 x (6 n_f^E) matrix P_E
 *
 * |E| (Pi_E tau_j)_a = -int_E div tau_j . p_a dE + int_{dE} (tau_j n) . p_a df.
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

  // volume term: -G (alpha_E, omega_E)
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

  // int_E pi_a : pi_b dE = |E| delta_ab
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
 * @brief Element matrix M_E of a_E^h, equation (14).
 * @param[in] faceGeom the geometry of the faces of the element
 * @param[in] numFaces the number of faces
 * @param[in] volume the element volume |E|
 * @param[in] diameter the element diameter h_E
 * @param[in] stabilizationLength h_E, or |E| / |dE| with consistency weights
 * @param[in] compliance the 6x6 matrix of D = C^{-1}
 * @param[in] projection the matrix P_E
 * @param[in,out] workspace a 6 x (6 n_f^E) scratch block
 * @param[out] complianceMatrix the (6 n_f^E) x (6 n_f^E) matrix M_E
 *
 * s_E is (15) with kappa_E = 1/(2 mu). Expanding (I - Pi_E) tau n gathers the quadratic terms in
 *   W = |E| D + kappa_E h sum_f |f| Lambda_{n_f}^T Lambda_{n_f},
 * so M_E = P_E^T W P_E plus cross and G_f terms local to each face.
 *
 * Option |E| / |dE| adds, with C = |E| kappa_E P_E^T P_E, T_E the unknowns of a constant stress and
 * R_c = I - T_E P_E^c (P_E^c the constant columns),
 *   sum_f sigma_{f,m}^T C_{f,mm} sigma_{f,m} + (R_c sigma)^T blkdiag_f(C_{f,cc}) (R_c sigma),
 * which vanishes on constant stresses.
 */
GEOS_HOST_DEVICE
inline void computeComplianceMatrix( FaceGeometry const * const faceGeom,
                              integer const numFaces,
                              real64 const volume,
                              real64 const diameter,
                              StabilizationLength const stabilizationLength,
                              real64 const (&compliance)[NUM_SYM_COMP][NUM_SYM_COMP],
                              MatrixSliceConst const & projection,
                              MatrixSlice const & workspace,
                              MatrixSlice const & complianceMatrix )
{
  integer const numStressDof = NUM_FACE_DOF * numFaces;

  // kappa_E = 1/(2 mu); h is h_E, equation (15), or |E| / |dE|
  real64 surfaceArea = 0.0;
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    surfaceArea += faceGeom[lf].area;
  }
  real64 const length = ( stabilizationLength != StabilizationLength::elementDiameter )
                        ? volume / surfaceArea
                        : diameter;
  real64 const stabScale = compliance[3][3] * length;

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

  // workspace = W P_E, then M_E = P_E^T workspace with a contiguous inner loop
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

  // P_E^T W P_E is symmetric, so only its lower triangle is formed
  for( integer i = 0; i < numStressDof; ++i )
  {
    real64 column[NUM_SYM_COMP];
    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      column[a] = projection( a, i );
    }

    for( integer j = 0; j <= i; ++j )
    {
      complianceMatrix( i, j ) = 0.0;
    }

    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      real64 const c = column[a];
      for( integer j = 0; j <= i; ++j )
      {
        complianceMatrix( i, j ) += c * workspace( a, j );
      }
    }
  }

  for( integer i = 0; i < numStressDof; ++i )
  {
    for( integer j = 0; j < i; ++j )
    {
      complianceMatrix( j, i ) = complianceMatrix( i, j );
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

    // cross weight Lambda_n^T int_f psi_k df
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
      // the projection column is read once instead of once per face mode
      real64 column[NUM_SYM_COMP];
      for( integer a = 0; a < NUM_SYM_COMP; ++a )
      {
        column[a] = projection( a, i );
      }

      for( integer k = 0; k < NUM_FACE_DOF; ++k )
      {
        real64 value = 0.0;
        for( integer a = 0; a < NUM_SYM_COMP; ++a )
        {
          value += column[a] * cross[a][k];
        }
        complianceMatrix( i, offset + k ) -= value;
        complianceMatrix( offset + k, i ) -= value;
      }
    }

    for( integer k = 0; k < NUM_FACE_DOF; ++k )
    {
      for( integer l = 0; l < NUM_FACE_DOF; ++l )
      {
        complianceMatrix( offset + k, offset + l ) += stabScale * gram[k][l];
      }
    }
  }

  if( stabilizationLength != StabilizationLength::hydraulicRadius )
  {
    return;
  }

  // C = |E| kappa_E P_E^T P_E, only its face blocks enter
  real64 const consistencyScale = volume * compliance[3][3];

  // C_{f,mm} on sigma, and C_{f,cc}, the identity part of R_c^T C_cc R_c
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    integer const offset = NUM_FACE_DOF * lf;
    for( integer k = 0; k < NUM_FACE_DOF; ++k )
    {
      for( integer l = 0; l < NUM_FACE_DOF; ++l )
      {
        if( ( k < 3 ) != ( l < 3 ) )
        {
          continue;
        }
        real64 value = 0.0;
        for( integer a = 0; a < NUM_SYM_COMP; ++a )
        {
          value += projection( a, offset + k ) * projection( a, offset + l );
        }
        complianceMatrix( offset + k, offset + l ) += consistencyScale * value;
      }
    }
  }

  // workspace = T_E^T C_cc on the constant columns, K = T_E^T C_cc T_E
  real64 K[NUM_SYM_COMP][NUM_SYM_COMP] = {};
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    FaceGeometry const & geom = faceGeom[lf];
    integer const offset = NUM_FACE_DOF * lf;

    real64 tractionMap[3][NUM_SYM_COMP];
    computeTractionMap( geom.normal, tractionMap );

    for( integer l = 0; l < 3; ++l )
    {
      for( integer a = 0; a < NUM_SYM_COMP; ++a )
      {
        real64 value = 0.0;
        for( integer k = 0; k < 3; ++k )
        {
          real64 c = 0.0;
          for( integer b = 0; b < NUM_SYM_COMP; ++b )
          {
            c += projection( b, offset + k ) * projection( b, offset + l );
          }
          value += geom.area * tractionMap[k][a] * consistencyScale * c;
        }
        workspace( a, offset + l ) = value;
      }
      for( integer a = 0; a < NUM_SYM_COMP; ++a )
      {
        for( integer b = 0; b < NUM_SYM_COMP; ++b )
        {
          K[a][b] += workspace( a, offset + l ) * geom.area * tractionMap[l][b];
        }
      }
    }
  }

  // the rest of R_c^T C_cc R_c, -C_cc T_E P^c - (P^c)^T T_E^T C_cc + (P^c)^T K P^c, on constant tractions
  for( integer j = 0; j < numStressDof; ++j )
  {
    if( j % NUM_FACE_DOF >= 3 )
    {
      continue;
    }

    real64 kp[NUM_SYM_COMP];
    for( integer a = 0; a < NUM_SYM_COMP; ++a )
    {
      kp[a] = -workspace( a, j );
      for( integer b = 0; b < NUM_SYM_COMP; ++b )
      {
        kp[a] += K[a][b] * projection( b, j );
      }
    }

    for( integer i = 0; i < numStressDof; ++i )
    {
      if( i % NUM_FACE_DOF >= 3 )
      {
        continue;
      }
      real64 value = 0.0;
      for( integer a = 0; a < NUM_SYM_COMP; ++a )
      {
        value += projection( a, i ) * kp[a] - workspace( a, i ) * projection( a, j );
      }
      complianceMatrix( i, j ) += value;
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
 * @param[out] divReconstruction the coefficients (alpha_E, omega_E) of Proposition 3.1
 * @param[out] projection the 6 x (6 n_f^E) matrix P_E
 * @param[in,out] workspace a 6 x (6 n_f^E) scratch block
 * @param[out] complianceMatrix the (6 numFaces) x (6 numFaces) matrix M_E
 */
GEOS_HOST_DEVICE
inline void computeElementOperators( FaceGeometry const * const faceGeom,
                                     integer const numFaces,
                                     real64 const (&elemCenter)[3],
                                     real64 const diameter,
                                     StabilizationLength const stabilizationLength,
                                     ElementMoments const & moments,
                                     real64 const (&compliance)[NUM_SYM_COMP][NUM_SYM_COMP],
                                     MatrixSlice const & divergence,
                                     MatrixSlice const & divReconstruction,
                                     MatrixSlice const & projection,
                                     MatrixSlice const & workspace,
                                     MatrixSlice const & complianceMatrix )
{
  computeDivergenceOperator( faceGeom, numFaces, elemCenter, divergence );

  computeDivergenceReconstruction( divergence.toSliceConst(), numFaces, moments, divReconstruction );

  computeProjectionOperator( faceGeom, numFaces, elemCenter, moments,
                             divReconstruction.toSliceConst(), projection );

  computeComplianceMatrix( faceGeom, numFaces, moments.volume, diameter, stabilizationLength, compliance,
                    projection.toSliceConst(), workspace, complianceMatrix );
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMELEMENTOPERATORS_HPP_
