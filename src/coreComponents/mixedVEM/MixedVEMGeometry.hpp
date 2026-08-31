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
 * @file MixedVEMGeometry.hpp
 *
 * Face frames, face moments and element moments of the lowest-order mixed VEM.
 *
 * Every integral needed by the element operators is a polynomial moment of the
 * fan decomposition, so it is evaluated in closed form: no quadrature loop
 * survives into the operator kernels.
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMGEOMETRY_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMGEOMETRY_HPP_

#include "mixedVEM/MixedVEMTypes.hpp"
#include "common/logger/Logger.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace mixedVEM
{

/**
 * @struct FaceGeometry
 * @brief The face-intrinsic data defining T_h(f) plus this element's orientation of it.
 *
 * The frame {t1, t2, n} and the moments are built from the face alone, so the two
 * elements sharing f produce identical data and their stress degrees of freedom are
 * directly identified at assembly. The only element-dependent entry is @p outwardSign.
 */
struct FaceGeometry
{
  /// Area centroid x_f, the origin of the local coordinates (xt, yt).
  real64 center[3];

  /// Reference unit normal n of the face, shared by both neighbours.
  real64 normal[3];

  /// First in-plane tangent, from the first non-degenerate edge of the node loop.
  real64 t1[3];

  /// Second in-plane tangent, t2 = n x t1.
  real64 t2[3];

  /// Area |f|.
  real64 area;

  /// m20 = |f|^{-1} int_f xt^2 df.
  real64 m20;

  /// m02 = |f|^{-1} int_f yt^2 df.
  real64 m02;

  /// m11 = |f|^{-1} int_f xt yt df.
  real64 m11;

  /// +-1: orientation of the node loop relative to n, right-handed with {t1, t2}.
  real64 loopOrientation;

  /// s_{E,f} = sign( n . (x_f - x_E) ): orientation of n relative to the outward normal.
  real64 outwardSign;
};

/**
 * @struct ElementMoments
 * @brief Volume and centred moments of an element, relative to the point x_E.
 */
struct ElementMoments
{
  /// Volume |E|.
  real64 volume;

  /// int_E (x - x_E) dE, zero when x_E is the exact barycenter.
  real64 firstMoment[3];

  /// int_E (x - x_E) x (x - x_E) dE.
  real64 secondMoment[3][3];
};

/**
 * @brief Zero the accumulators of an ElementMoments.
 * @param[out] moments the moments to reset
 */
GEOS_HOST_DEVICE
inline void initialize( ElementMoments & moments )
{
  moments.volume = 0.0;
  LvArray::tensorOps::fill< 3 >( moments.firstMoment, 0.0 );
  LvArray::tensorOps::fill< 3, 3 >( moments.secondMoment, 0.0 );
}

/**
 * @brief Build the face frame, area, centroid and second moments of a polygonal face.
 * @tparam COORDS callable with signature (integer node, integer dim) -> real64
 * @param[in] X the coordinates of the face nodes, in loop order
 * @param[in] numNodes the number of face nodes
 * @param[in] referenceNormal the face-intrinsic unit normal n
 * @param[out] geom the resulting face geometry, with @p outwardSign left untouched
 *
 * The polygon is fanned from the node average and every triangle contributes its
 * signed area n . [(x_k - xh) ^ (x_{k+1} - xh)] / 2, so re-entrant polygons cancel
 * correctly. The moments are then taken about the exact centroid, which is what makes
 * int_f xt df = int_f yt df = 0 and collapses the face integrals to closed form.
 */
template< typename COORDS >
GEOS_HOST_DEVICE
void computeFaceGeometry( COORDS const & X,
                          integer const numNodes,
                          real64 const (&referenceNormal)[3],
                          FaceGeometry & geom )
{
  LvArray::tensorOps::copy< 3 >( geom.normal, referenceNormal );
  LvArray::tensorOps::normalize< 3 >( geom.normal );

  // The two elements sharing f must agree on the frame, otherwise their stress degrees of
  // freedom mean different things. A mesh face normal is only oriented with respect to one
  // of the two, and in parallel each rank may orient it its own way, so both the sign and
  // the tangent are fixed here from the direction alone: the sign makes the dominant
  // component positive, and t1 comes from the global axis least aligned with n. Both are
  // pure functions of n, hence identical on every rank and independent of the node loop.
  integer dominant = 0;
  integer weakest = 0;
  for( integer i = 1; i < 3; ++i )
  {
    if( LvArray::math::abs( geom.normal[i] ) > LvArray::math::abs( geom.normal[dominant] ) )
    {
      dominant = i;
    }
    if( LvArray::math::abs( geom.normal[i] ) < LvArray::math::abs( geom.normal[weakest] ) )
    {
      weakest = i;
    }
  }

  if( geom.normal[dominant] < 0.0 )
  {
    LvArray::tensorOps::scale< 3 >( geom.normal, -1.0 );
  }

  real64 axis[3] = { 0.0, 0.0, 0.0 };
  axis[weakest] = 1.0;

  LvArray::tensorOps::copy< 3 >( geom.t1, axis );
  LvArray::tensorOps::scaledAdd< 3 >( geom.t1, geom.normal, -LvArray::tensorOps::AiBi< 3 >( axis, geom.normal ) );
  LvArray::tensorOps::normalize< 3 >( geom.t1 );

  // fan apex: the node average, always inside the plane of the face
  real64 apex[3] = { 0.0, 0.0, 0.0 };
  for( integer k = 0; k < numNodes; ++k )
  {
    for( integer i = 0; i < 3; ++i )
    {
      apex[i] += X( k, i );
    }
  }
  LvArray::tensorOps::scale< 3 >( apex, 1.0 / numNodes );

  // signed area and area-weighted centroid of the triangle fan
  real64 signedArea = 0.0;
  real64 centroid[3] = { 0.0, 0.0, 0.0 };

  for( integer k = 0; k < numNodes; ++k )
  {
    integer const kp1 = ( k + 1 ) % numNodes;

    real64 e0[3], e1[3], cross[3];
    for( integer i = 0; i < 3; ++i )
    {
      e0[i] = X( k, i ) - apex[i];
      e1[i] = X( kp1, i ) - apex[i];
    }
    LvArray::tensorOps::crossProduct( cross, e0, e1 );

    real64 const areaK = 0.5 * LvArray::tensorOps::AiBi< 3 >( cross, geom.normal );

    signedArea += areaK;
    for( integer i = 0; i < 3; ++i )
    {
      centroid[i] += areaK * ( apex[i] + X( k, i ) + X( kp1, i ) ) / 3.0;
    }
  }

  geom.loopOrientation = ( signedArea < 0.0 ) ? -1.0 : 1.0;
  geom.area = LvArray::math::abs( signedArea );

  LvArray::tensorOps::copy< 3 >( geom.center, centroid );
  LvArray::tensorOps::scale< 3 >( geom.center, 1.0 / signedArea );

  LvArray::tensorOps::crossProduct( geom.t2, geom.normal, geom.t1 );
  LvArray::tensorOps::normalize< 3 >( geom.t2 );

  // re-orthogonalize t1 = t2 ^ n to remove floating point drift
  LvArray::tensorOps::crossProduct( geom.t1, geom.t2, geom.normal );
  LvArray::tensorOps::normalize< 3 >( geom.t1 );

  // second moments about x_f, fanned from x_f so that the linear terms drop out
  real64 m20 = 0.0, m02 = 0.0, m11 = 0.0;

  for( integer k = 0; k < numNodes; ++k )
  {
    integer const kp1 = ( k + 1 ) % numNodes;

    real64 r0[3], r1[3];
    for( integer i = 0; i < 3; ++i )
    {
      r0[i] = X( k, i ) - geom.center[i];
      r1[i] = X( kp1, i ) - geom.center[i];
    }

    real64 const x1 = LvArray::tensorOps::AiBi< 3 >( r0, geom.t1 );
    real64 const y1 = LvArray::tensorOps::AiBi< 3 >( r0, geom.t2 );
    real64 const x2 = LvArray::tensorOps::AiBi< 3 >( r1, geom.t1 );
    real64 const y2 = LvArray::tensorOps::AiBi< 3 >( r1, geom.t2 );

    // T_h(f) is defined on a planar face, so any out of plane node invalidates the moments
    GEOS_ASSERT_MSG( LvArray::math::abs( LvArray::tensorOps::AiBi< 3 >( r0, geom.normal ) )
                     <= 1e-8 * LvArray::math::sqrt( geom.area ),
                     "MixedVEM: face nodes are not coplanar." );

    // signed twice-area of the triangle (x_f, x_k, x_{k+1}) in the local frame
    real64 const jac = x1 * y2 - x2 * y1;

    m20 += jac * ( x1 * x1 + x1 * x2 + x2 * x2 ) / 12.0;
    m02 += jac * ( y1 * y1 + y1 * y2 + y2 * y2 ) / 12.0;
    m11 += jac * ( x1 * y1 / 12.0 + x1 * y2 / 24.0 + x2 * y1 / 24.0 + x2 * y2 / 12.0 );
  }

  // normalize by the signed area so that m20, m02 stay positive for either loop orientation
  real64 const invArea = geom.loopOrientation / geom.area;

  geom.m20 = m20 * invArea;
  geom.m02 = m02 * invArea;
  geom.m11 = m11 * invArea;

  geom.outwardSign = 1.0;
}

/**
 * @brief Orient a face geometry with respect to an element center.
 * @param[in] elemCenter the point x_E
 * @param[in,out] geom the face geometry whose @p outwardSign is set
 *
 * s_{E,f} = sign( n . (x_f - x_E) ), valid for elements star-shaped about x_E (A1).
 */
GEOS_HOST_DEVICE
inline void orientFace( real64 const (&elemCenter)[3],
                        FaceGeometry & geom )
{
  real64 df[3];
  LvArray::tensorOps::copy< 3 >( df, geom.center );
  LvArray::tensorOps::subtract< 3 >( df, elemCenter );

  geom.outwardSign = ( LvArray::tensorOps::AiBi< 3 >( df, geom.normal ) < 0.0 ) ? -1.0 : 1.0;
}

/**
 * @brief Add the tetrahedral fan of one face to the element moments.
 * @tparam COORDS callable with signature (integer node, integer dim) -> real64
 * @param[in] X the coordinates of the face nodes, in loop order
 * @param[in] numNodes the number of face nodes
 * @param[in] geom the geometry of that face
 * @param[in] elemCenter the point x_E
 * @param[in,out] moments the accumulated element moments
 *
 * Each tetrahedron (x_E, x_f, x_k, x_{k+1}) has one vertex at r = 0, so with
 * a = x_f - x_E, b = x_k - x_E, c = x_{k+1} - x_E,
 *   int_T r dT   = V (a + b + c)/4,
 *   int_T r@r dT = (V/20) [ a@a + b@b + c@c + (a+b+c)@(a+b+c) ],
 * which are exact, not quadrature approximations.
 */
template< typename COORDS >
GEOS_HOST_DEVICE
void accumulateElementMoments( COORDS const & X,
                               integer const numNodes,
                               FaceGeometry const & geom,
                               real64 const (&elemCenter)[3],
                               ElementMoments & moments )
{
  real64 a[3];
  LvArray::tensorOps::copy< 3 >( a, geom.center );
  LvArray::tensorOps::subtract< 3 >( a, elemCenter );

  // the loop orientation and the outward sign together make every fan volume positive
  real64 const orientation = geom.loopOrientation * geom.outwardSign;

  for( integer k = 0; k < numNodes; ++k )
  {
    integer const kp1 = ( k + 1 ) % numNodes;

    real64 b[3], c[3], cross[3];
    for( integer i = 0; i < 3; ++i )
    {
      b[i] = X( k, i ) - elemCenter[i];
      c[i] = X( kp1, i ) - elemCenter[i];
    }

    LvArray::tensorOps::crossProduct( cross, b, c );
    real64 const volume = orientation * LvArray::tensorOps::AiBi< 3 >( a, cross ) / 6.0;

    real64 s[3];
    for( integer i = 0; i < 3; ++i )
    {
      s[i] = a[i] + b[i] + c[i];
    }

    moments.volume += volume;

    for( integer i = 0; i < 3; ++i )
    {
      moments.firstMoment[i] += 0.25 * volume * s[i];
    }

    real64 const w = volume / 20.0;
    for( integer i = 0; i < 3; ++i )
    {
      for( integer j = 0; j < 3; ++j )
      {
        moments.secondMoment[i][j] += w * ( a[i] * a[j] + b[i] * b[j] + c[i] * c[j] + s[i] * s[j] );
      }
    }
  }
}

/**
 * @brief Inertia matrix A_E of the rotational part of the divergence reconstruction.
 * @param[in] moments the element moments
 * @param[out] inertia the symmetric positive definite matrix A_E
 *
 * A_E omega = int_E (x - x_E) ^ [ omega ^ (x - x_E) ] dE = [ tr(M) I - M ] omega,
 * with M the second moment; this is the 3x3 system (8) of Proposition 3.1.
 */
GEOS_HOST_DEVICE
inline void computeInertia( ElementMoments const & moments,
                            real64 (& inertia)[3][3] )
{
  real64 const trace = moments.secondMoment[0][0] + moments.secondMoment[1][1] + moments.secondMoment[2][2];

  for( integer i = 0; i < 3; ++i )
  {
    for( integer j = 0; j < 3; ++j )
    {
      inertia[i][j] = ( ( i == j ) ? trace : 0.0 ) - moments.secondMoment[i][j];
    }
  }
}

/**
 * @brief Diameter h_E of an element, used to scale the stabilization.
 * @tparam COORDS callable with signature (integer node, integer dim) -> real64
 * @param[in] X the coordinates of the element nodes
 * @param[in] numNodes the number of element nodes
 * @return the largest distance between two nodes
 */
template< typename COORDS >
GEOS_HOST_DEVICE
real64 computeElementDiameter( COORDS const & X,
                               integer const numNodes )
{
  real64 maxDistanceSquared = 0.0;

  for( integer k = 0; k < numNodes; ++k )
  {
    for( integer l = k + 1; l < numNodes; ++l )
    {
      real64 distanceSquared = 0.0;
      for( integer i = 0; i < 3; ++i )
      {
        real64 const d = X( k, i ) - X( l, i );
        distanceSquared += d * d;
      }
      maxDistanceSquared = LvArray::math::max( maxDistanceSquared, distanceSquared );
    }
  }

  return LvArray::math::sqrt( maxDistanceSquared );
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMGEOMETRY_HPP_
