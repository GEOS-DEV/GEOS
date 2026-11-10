/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */

/**
 * @file MPMGeometryMeasures.hpp
 * @brief Physical interface-area and nodal-volume measures for native 2D MPM.
 *
 * These host/device helpers integrate the bilinear (Q1) nodal hat over a
 * reconstructed interface in the meridional or Cartesian grid plane. They
 * convert a geometric interface line into physical surface area using either
 * a constant out-of-plane thickness or the cylindrical full-ring measure.
 * This geometric integration is independent of the mechanical particle
 * interpolation scheme (single-point B-spline, CPDI, or CPTI).
 */
#ifndef GEOS_SOLIDMECHANICS_MPMGEOMETRYMEASURES_HPP
#define GEOS_SOLIDMECHANICS_MPMGEOMETRYMEASURES_HPP

#include "common/DataTypes.hpp"
#include <cmath>
#include <cfloat>

namespace geos
{
namespace mpmMeasures
{

/**
 * @brief Integrates one nodal Q1 hat over an interface line, returning physical area.
 *
 * Work in coordinates relative to the reference grid node. Its compact
 * rectangular Q1 support is |x| < hx, |y| < hy, with
 *
 *   N_I(x,y) = (1 - |x|/hx) (1 - |y|/hy).
 *
 * The interface is the infinite line n . (x - position) = 0. Clipping it
 * against the Q1 support gives a finite segment Gamma_I. The result is
 *
 *   Cartesian 2D: integral_Gamma_I thickness * N_I ds,
 *   axisymmetric: integral_Gamma_I 2*pi*r * N_I ds,
 *
 * where r = rI + x and rI = measureScale is the absolute radius of the
 * reference node. In the Cartesian case measureScale is instead the physical
 * out-of-plane thickness. Thus the result is a physical crack/contact area,
 * not the unweighted interface length; the two differ even in native 2D.
 *
 * The line is parameterized by its closest point to the node plus distance
 * along a unit tangent. First clip the parameter interval to the nodal hat
 * and, in axisymmetry, to r >= 0. Then split the surviving segment wherever
 * x or y changes sign, because the Q1 absolute-value factors change slope.
 * On each segment N_I is quadratic in arclength; multiplication by radius
 * makes the axisymmetric integrand cubic. Two-point Gauss-Legendre integration
 * is therefore exact for the polynomial on each such segment.
 *
 * @tparam AXISYMMETRIC Use a full 2*pi*r ring measure if true, otherwise
 *                      a constant plane-strain thickness.
 * @tparam NORMAL Array-like type supplying the two in-plane normal components.
 * @tparam POSITION Array-like type supplying an in-plane point on the line,
 *                  relative to the reference node.
 * @param hx Positive radial/x nodal support half-width (grid spacing).
 * @param hy Positive axial/y nodal support half-width (grid spacing).
 * @param measureScale Plane-strain thickness or absolute reference-node radius.
 * @param normal Line normal; its magnitude need not be unity, but must be nonzero.
 * @param position In-plane position on the interface relative to the node.
 * @return Physical nodal interface area, or zero for an empty/invalid support.
 *
 * @note The plane-strain path requires a strictly positive thickness. The
 *       axisymmetric path admits rI = 0; it uses no artificial radius floor.
 *       Both paths assume positive grid spacings, and the radius passed to an
 *       axisymmetric call is the physical radius of the reference node.
 */
template< bool AXISYMMETRIC, class NORMAL, class POSITION >
GEOS_HOST_DEVICE real64 interfaceLineIntegral( real64 const hx,
                                              real64 const hy,
                                              real64 const measureScale,
                                              NORMAL const & normal,
                                              POSITION const & position )
{
  // Construct a unit line normal. Degenerate normals or nonpositive spacing
  // cannot define a meaningful nodal line-area integral.
  real64 const norm = ::sqrt( normal[0] * normal[0] + normal[1] * normal[1] );
  if( !(norm > 0 && hx > 0 && hy > 0) ) return 0;

  // A Cartesian calculation needs physical thickness. Axisymmetry instead
  // interprets measureScale as node radius, for which zero is permitted.
  if constexpr( !AXISYMMETRIC )
  {
    if( !(measureScale > 0) ) return 0;
  }

  // Signed distance and the closest point of the infinite line to this node.
  // A unit tangent makes the line parameter t equal to physical arclength.
  real64 const n[2] = { normal[0] / norm, normal[1] / norm };
  real64 const distance = n[0] * position[0] + n[1] * position[1];
  real64 const center[2] = { distance * n[0], distance * n[1] };
  real64 const tangent[2] = { -n[1], n[0] }, h[2] = { hx, hy };
  real64 lo = -DBL_MAX, hi = DBL_MAX;

  // Intersect center+t*tangent with the two strips [-h[d], h[d]]. A line
  // parallel to a support face is either outside that strip or imposes no
  // additional bound on t. Touching a face alone has zero interface measure.
  for( int d = 0; d < 2; ++d )
  {
    if( ::fabs( tangent[d] ) < DBL_EPSILON )
    {
      if( ::fabs( center[d] ) >= h[d] ) return 0;
      continue;
    }
    real64 a = (-h[d] - center[d]) / tangent[d], b = (h[d] - center[d]) / tangent[d];
    if( a > b )
    {
      real64 const tmp = a;
      a = b;
      b = tmp;
    }
    if( a > lo ) lo = a;
    if( b < hi ) hi = b;
  }

  // Only the nonnegative-radius half-plane represents physical material.
  // Clipping in arclength coordinates retains support at the symmetry axis:
  // a node with rI = 0 can still integrate a finite interval at r > 0.
  // There is no division by radius or regularization parameter here.
  if constexpr( AXISYMMETRIC )
  {
    if( ::fabs( tangent[0] ) < DBL_EPSILON )
    {
      if( !(measureScale + center[0] > 0) ) return 0;
    }
    else
    {
      real64 const axis = -(measureScale + center[0]) / tangent[0];
      if( tangent[0] > 0 && axis > lo ) lo = axis;
      if( tangent[0] < 0 && axis < hi ) hi = axis;
    }
  }

  // The clipped line has no nodal area when its parameter interval is empty.
  if( !(hi > lo) ) return 0;

  // Split at x = 0 and y = 0, the knots of the tensor-product Q1 hat.
  // At most two interior knots are possible, so four break points suffice.
  real64 breaks[4] = { lo, hi, 0, 0 };
  int count = 2;
  for( int d = 0; d < 2; ++d ) if( ::fabs( tangent[d] ) >= DBL_EPSILON )
  {
    real64 const t = -center[d] / tangent[d];
    if( t > lo && t < hi ) breaks[count++] = t;
  }

  // Sort the endpoints and interior knots into successive intervals on which
  // both absolute-value signs, and hence the polynomial factors, are fixed.
  for( int a = 1; a < count; ++a ) for( int b = a; b > 0 && breaks[b] < breaks[b-1]; --b )
  {
    real64 const tmp = breaks[b];
    breaks[b] = breaks[b-1];
    breaks[b-1] = tmp;
  }

  // Two symmetric Gauss points integrate up to cubic order exactly. The
  // product of the two Q1 hat factors is quadratic along a line; the
  // axisymmetric r = rI+x factor raises its degree by one, not beyond cubic.
  real64 integral = 0;
  for( int a = 1; a < count; ++a )
  {
    real64 const mid = 0.5 * (breaks[a] + breaks[a-1]), half = 0.5 * (breaks[a] - breaks[a-1]);
    for( int sign = -1; sign <= 1; sign += 2 )
    {
      real64 const t = mid + sign * half / ::sqrt( 3.0 );
      real64 weight = 1;
      for( int d = 0; d < 2; ++d ) weight *= 1 - ::fabs( center[d] + t * tangent[d] ) / h[d];

      // Multiply the Q1 value by physical radius for the axisymmetric path.
      // The positive-radius restriction was imposed by the preceding clip.
      if constexpr( AXISYMMETRIC ) weight *= measureScale + center[0] + t * tangent[0];
      integral += half * weight;
    }
  }

  // Convert the Q1-weighted meridional line integral into physical area:
  // either a constant thickness or the full 2*pi circumferential factor.
  return integral * (AXISYMMETRIC ? 6.283185307179586476925286766559 : measureScale);
}

/**
 * @brief Returns Q1-weighted physical interface area for a Cartesian 2D grid.
 *
 * The line geometry is two-dimensional, but its physical surface area is
 * thickness times the Q1-weighted line length. This wrapper selects the
 * constant-thickness branch of interfaceLineIntegral().
 *
 * @param hx Positive x-direction nodal support half-width.
 * @param hy Positive y-direction nodal support half-width.
 * @param thickness Positive out-of-plane physical thickness.
 * @param normal Meridional/Cartesian line normal.
 * @param position Point on the interface relative to the reference node.
 */
template< class NORMAL, class POSITION >
GEOS_HOST_DEVICE real64 lineArea( real64 const hx,
                                 real64 const hy,
                                 real64 const thickness,
                                 NORMAL const & normal,
                                 POSITION const & position )
{
  return interfaceLineIntegral< false >( hx, hy, thickness, normal, position );
}

/**
 * @brief Returns the full-ring Q1-weighted area of an axisymmetric interface.
 *
 * Uses the radial position r = rI+x along the clipped line to evaluate
 * integral_Gamma 2*pi*r*N_I ds. Unlike the plane-strain wrapper, its third
 * argument is the *node radius*, not a thickness or regularization length.
 * The integration itself clips the line at the physical symmetry axis.
 *
 * @param hx Positive radial grid spacing.
 * @param hy Positive axial grid spacing.
 * @param rI Absolute reference-node radius; zero is allowed.
 * @param normal Meridional line normal.
 * @param position Point on the interface relative to the reference node.
 */
template< class NORMAL, class POSITION >
GEOS_HOST_DEVICE real64 axisymmetricLineArea( real64 const hx,
                                             real64 const hy,
                                             real64 const rI,
                                             NORMAL const & normal,
                                             POSITION const & position )
{
  return interfaceLineIntegral< true >( hx, hy, rI, normal, position );
}

/**
 * @brief Integrates the Q1 nodal hat over its full, nonnegative-radius ring volume.
 *
 * The axial Q1 factor integrates to hy over [-hy, hy]. The remaining radial
 * integral is obtained by passing a horizontal interface through the node to
 * axisymmetricLineArea(), whose normal is (0,1). Thus
 *
 *   V_I = integral_{r>=0} 2*pi*r*N_I(r,y) dr dy
 *       = hy * integral_{y=0,r>=0} 2*pi*r*N_I(r,0) dr.
 *
 * At the axis (rI = 0) the result is 2*pi*hx^2*hy/6. For a node whose
 * support lies wholly outside the axis (rI >= hx), it is 2*pi*rI*hx*hy.
 * Intermediate node radii are handled by the same exact clipped integration.
 *
 * @param hx Positive radial grid spacing.
 * @param hy Positive axial grid spacing.
 * @param rI Physical reference-node radius, assumed nonnegative.
 * @return Full-ring physical volume associated with the nodal Q1 support.
 */
GEOS_HOST_DEVICE inline real64 axisymmetricNodalVolume( real64 const hx,
                                                        real64 const hy,
                                                        real64 const rI )
{
  // At y = 0 the axial hat equals unity; multiplication by hy supplies
  // its complete axial integral without a second quadrature dimension.
  real64 const normal[2] = { 0, 1 }, position[2] = { 0, 0 };
  return hy * axisymmetricLineArea( hx, hy, rI, normal, position );
}

} // namespace mpmMeasures
} // namespace geos

#endif
