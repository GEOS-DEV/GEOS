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
 * @file MixedVEMFaceBasis.hpp
 *
 * The six traction modes spanning T_h(f) and the face moments they generate.
 *
 * T_h(f) is the trace of the rigid body motions on the face, so with r = x - x_f and
 * L = sqrt(m20+m02) the polar radius of gyration, the basis is written directly in the
 * global frame as
 *
 *   phi_k = e_k/|f|,  phi_{3+k} = (e_k ^ r)/(|f| L),   k = 0, 1, 2.
 *
 * The three rotational modes span the same space as the local {n ^ r, xt n, yt n} used
 * elsewhere in the literature, since e_k ^ r resolves into an in plane rotation about n
 * and a normal component varying linearly across the face.
 *
 * Writing all six in the global frame removes the arbitrary choice of t1, turns a
 * prescribed traction or displacement into its plain Cartesian components, and makes a
 * global rigid motion produce the same rotational multiplier values on every face. That
 * last property is what lets an unknown-based algebraic multigrid reproduce the whole
 * near null space from its own per-unknown constants, on a mesh of any orientation.
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMFACEBASIS_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMFACEBASIS_HPP_

#include "mixedVEM/MixedVEMGeometry.hpp"

namespace geos
{

namespace mixedVEM
{

/**
 * @brief Cross product of a Cartesian axis with a vector, (e_k ^ v).
 * @param[in] k the axis index
 * @param[in] v the vector
 * @param[out] out the result
 */
GEOS_HOST_DEVICE
inline void axisCross( integer const k,
                       real64 const (&v)[3],
                       real64 (& out)[3] )
{
  integer const a = ( k + 1 ) % 3;
  integer const b = ( k + 2 ) % 3;

  out[k] = 0.0;
  out[a] = -v[b];
  out[b] = v[a];
}

/**
 * @brief Second moment tensor M_pq = int_f r_p r_q df of the face, in the global frame.
 * @param[in] geom the face geometry
 * @param[out] M the tensor
 */
GEOS_HOST_DEVICE
inline void computeFaceSecondMoment( FaceGeometry const & geom,
                                     real64 (& M)[3][3] )
{
  for( integer p = 0; p < 3; ++p )
  {
    for( integer q = 0; q < 3; ++q )
    {
      M[p][q] = geom.m20 * geom.t1[p] * geom.t1[q]
                + geom.m11 * ( geom.t1[p] * geom.t2[q] + geom.t2[p] * geom.t1[q] )
                + geom.m02 * geom.t2[p] * geom.t2[q];
    }
  }
}

/**
 * @brief Evaluate the six traction modes of T_h(f) at a point of the face.
 * @param[in] geom the face geometry
 * @param[in] x the evaluation point, assumed to lie in the plane of the face
 * @param[out] phi the six vector values, phi[j][i]
 */
GEOS_HOST_DEVICE
inline void evaluateFaceBasis( FaceGeometry const & geom,
                               real64 const (&x)[3],
                               real64 (& phi)[NUM_FACE_DOF][3] )
{
  real64 r[3];
  LvArray::tensorOps::copy< 3 >( r, x );
  LvArray::tensorOps::subtract< 3 >( r, geom.center );

  real64 const invArea = 1.0 / geom.area;
  real64 const sL = invArea / LvArray::math::sqrt( geom.m20 + geom.m02 );

  for( integer i = 0; i < 3; ++i )
  {
    phi[0][i] = 0.0;
    phi[1][i] = 0.0;
    phi[2][i] = 0.0;
    phi[i][i] = invArea;
  }

  for( integer k = 0; k < 3; ++k )
  {
    real64 w[3];
    axisCross( k, r, w );

    for( integer i = 0; i < 3; ++i )
    {
      phi[3 + k][i] = sL * w[i];
    }
  }
}

/**
 * @brief Zeroth moments int_f phi_j df of the traction modes.
 * @param[in] geom the face geometry
 * @param[out] mean the six vectors, mean[j][i]
 *
 * The three linear modes integrate to zero because x_f is the exact centroid,
 * so only the constant modes survive: int_f phi_j = (e_x, e_y, e_z, 0, 0, 0).
 */
GEOS_HOST_DEVICE
inline void computeFaceBasisMeans( FaceGeometry const & geom,
                                   real64 (& mean)[NUM_FACE_DOF][3] )
{
  GEOS_UNUSED_VAR( geom );

  for( integer i = 0; i < 3; ++i )
  {
    mean[0][i] = 0.0;
    mean[1][i] = 0.0;
    mean[2][i] = 0.0;
    mean[i][i] = 1.0;
    mean[3][i] = 0.0;
    mean[4][i] = 0.0;
    mean[5][i] = 0.0;
  }
}

/**
 * @brief First moments N_j = int_f phi_j @ (x - x_E) df of the traction modes.
 * @param[in] geom the face geometry
 * @param[in] elemCenter the point x_E
 * @param[out] N the six tensors, N[j][p][q]
 *
 * Writing x - x_E = df + xt t1 + yt t2 with df = x_f - x_E and using
 * int_f xt = int_f yt = 0, every entry reduces to the second moments of the face.
 * N_j carries both element operators: its skew part gives the rotational rows of
 * B_E, its symmetric part gives the face contribution to the projection.
 */
GEOS_HOST_DEVICE
inline void computeFaceBasisMoments( FaceGeometry const & geom,
                                     real64 const (&elemCenter)[3],
                                     real64 (& N)[NUM_FACE_DOF][3][3] )
{
  real64 df[3];
  LvArray::tensorOps::copy< 3 >( df, geom.center );
  LvArray::tensorOps::subtract< 3 >( df, elemCenter );

  real64 const sL = 1.0 / LvArray::math::sqrt( geom.m20 + geom.m02 );

  real64 M[3][3];
  computeFaceSecondMoment( geom, M );

  for( integer p = 0; p < 3; ++p )
  {
    for( integer q = 0; q < 3; ++q )
    {
      // constant modes: only the offset df survives
      N[0][p][q] = ( p == 0 ) ? df[q] : 0.0;
      N[1][p][q] = ( p == 1 ) ? df[q] : 0.0;
      N[2][p][q] = ( p == 2 ) ? df[q] : 0.0;
    }
  }

  // rotational modes: df drops out with int_f r, leaving sL (e_k ^ M_q) on each column
  for( integer q = 0; q < 3; ++q )
  {
    real64 const column[3] = { M[0][q], M[1][q], M[2][q] };

    for( integer k = 0; k < 3; ++k )
    {
      real64 w[3];
      axisCross( k, column, w );

      for( integer p = 0; p < 3; ++p )
      {
        N[3 + k][p][q] = sL * w[p];
      }
    }
  }
}

/**
 * @brief Rotational moment int_f (x - x_E) ^ phi_j df obtained from N_j.
 * @param[in] N the first moment tensor of one mode
 * @param[out] b the resulting vector
 *
 * int_f (x - x_E) ^ phi_j = -axial(N_j), the skew part of N_j read as a vector.
 */
GEOS_HOST_DEVICE
inline void faceBasisRotationalMoment( real64 const (&N)[3][3],
                                       real64 (& b)[3] )
{
  b[0] = N[2][1] - N[1][2];
  b[1] = N[0][2] - N[2][0];
  b[2] = N[1][0] - N[0][1];
}

/**
 * @brief Gram matrix int_f phi_i . phi_j df of the traction modes.
 * @param[in] geom the face geometry
 * @param[out] gram the 6x6 matrix
 *
 * The constant block is |f|^{-1} I and the two blocks do not mix, because x_f is the
 * centroid and int_f r vanishes. The rotational block follows from
 * (e_i ^ r) . (e_j ^ r) = delta_ij |r|^2 - r_i r_j averaged over the face, the second
 * moments carried by the geometry already being normalized by |f|.
 */
GEOS_HOST_DEVICE
inline void computeFaceBasisGram( FaceGeometry const & geom,
                                  real64 (& gram)[NUM_FACE_DOF][NUM_FACE_DOF] )
{
  real64 const invArea = 1.0 / geom.area;
  real64 const trace = geom.m20 + geom.m02;
  real64 const sL2 = invArea / trace;

  real64 M[3][3];
  computeFaceSecondMoment( geom, M );

  for( integer i = 0; i < NUM_FACE_DOF; ++i )
  {
    for( integer j = 0; j < NUM_FACE_DOF; ++j )
    {
      gram[i][j] = 0.0;
    }
  }

  for( integer i = 0; i < 3; ++i )
  {
    gram[i][i] = invArea;

    for( integer j = 0; j < 3; ++j )
    {
      gram[3 + i][3 + j] = sL2 * ( ( ( i == j ) ? trace : 0.0 ) - M[i][j] );
    }
  }
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMFACEBASIS_HPP_
