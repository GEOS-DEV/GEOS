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


/**
 * @brief Integrate a 3D Q1 nodal hat over its intersection with a plane.
 *
 * The plane n.(x-position)=0 cuts the nodal support box [-hx,hx] x
 * [-hy,hy] x [-hz,hz]. Intersect its 12 edges using the segment parameter
 * t in [0,1], deduplicate shared-vertex hits, and sort the resulting
 * convex polygon about its centroid (which lies inside the cut).
 *
 * Q1 is the product of three piecewise-linear hats. Clip the polygon at
 * x=0, y=0, and z=0, giving regions with a single cubic polynomial.
 * Fan-triangulate each region and use a six-point degree-four triangle
 * rule with the physical Jacobian |(b-a) x (c-a)|. The rule's weights
 * sum to 1/2, so the result is physical area, not reference area.
 * The quadrature is exact on nondegenerate planar cuts up to roundoff.
 *
 * The plane normal need not be unit length. No allocation, sorting library,
 * or global geometry search is needed, so this routine can run on device.
 * This is an idealized plane-area operator, independent of how the
 * particle surface was detected or which mechanical shape functions map it.
 *
 * @return Integral of N_I dA, or zero if the plane misses the support.
 */
template< class EXTENTS, class NORMAL, class POSITION >
GEOS_HOST_DEVICE real64 nodalPlaneArea3D( EXTENTS const & h,
                                         NORMAL const & normal,
                                         POSITION const & position )
{
  for( int d = 0; d < 3; ++d )
  {
    if( !(h[d] > 0.0 && h[d] < DBL_MAX) ) return 0.0;
  }

  real64 const normalLength = ::sqrt( normal[0] * normal[0] + normal[1] * normal[1] +
                                      normal[2] * normal[2] );
  if( !(normalLength > 0.0 && normalLength < DBL_MAX) ) return 0.0;

  real64 const corners[8][3] = {
    {-h[0],-h[1],-h[2]}, {h[0],-h[1],-h[2]},
    { h[0], h[1],-h[2]}, {-h[0],h[1],-h[2]},
    {-h[0],-h[1], h[2]}, {h[0],-h[1], h[2]},
    { h[0], h[1], h[2]}, {-h[0],h[1],h[2]}
  };
  int const edges[12][2] = {
    {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6},
    {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}
  };

  real64 const maxH = ::fmax(h[0], ::fmax(h[1], h[2]));
  real64 const coordinateTol = 64.0 * DBL_EPSILON * maxH;
  real64 const signedTol = coordinateTol * normalLength;
  real64 const distanceTol2 = coordinateTol * coordinateTol;
  real64 signedDistance[8];
  for( int k = 0; k < 8; ++k )
  {
    signedDistance[k] = 0.0;
    for( int d = 0; d < 3; ++d )
      signedDistance[k] += normal[d] * (corners[k][d] - position[d]);
  }

  // Every accepted point belongs to the finite support box. A plane through
  // a box vertex/edge is seen by several edges but must contribute once.
  real64 vertices[12][3] = {};
  int count = 0;
  for( int edge = 0; edge < 12; ++edge )
  {
    int const a = edges[edge][0], b = edges[edge][1];
    real64 const sa = signedDistance[a], sb = signedDistance[b];
    real64 parameters[2] = {};
    int nParameters = 0;
    if( ::fabs(sa) <= signedTol ) parameters[nParameters++] = 0.0;
    if( ::fabs(sb) <= signedTol ) parameters[nParameters++] = 1.0;
    if( (sa < -signedTol && sb > signedTol) || (sa > signedTol && sb < -signedTol) )
      parameters[nParameters++] = sa / (sa-sb);

    for( int q = 0; q < nParameters; ++q )
    {
      real64 const t = parameters[q];
      real64 candidate[3];
      for( int d = 0; d < 3; ++d )
        candidate[d] = corners[a][d] + t * (corners[b][d] - corners[a][d]);

      bool duplicate = false;
      for( int k = 0; k < count; ++k )
      {
        real64 squared = 0.0;
        for( int d = 0; d < 3; ++d )
          squared += (candidate[d] - vertices[k][d]) * (candidate[d] - vertices[k][d]);
        if( squared <= distanceTol2 ) { duplicate = true; break; }
      }
      if( !duplicate && count < 12 )
      {
        for( int d = 0; d < 3; ++d ) vertices[count][d] = candidate[d];
        ++count;
      }
    }
  }
  if( count < 3 ) return 0.0;

  // The arithmetic mean of a convex polygon's vertices is an interior
  // point. Use it only to sort polar angles; unlike the source surface
  // position, it cannot lie outside a valid support intersection.
  real64 center[3] = {};
  for( int k = 0; k < count; ++k )
    for( int d = 0; d < 3; ++d ) center[d] += vertices[k][d] / count;
  real64 tangent[3];
  for( int d = 0; d < 3; ++d ) tangent[d] = vertices[0][d] - center[d];
  real64 bitangent[3] = {
    normal[1]*tangent[2] - normal[2]*tangent[1],
    normal[2]*tangent[0] - normal[0]*tangent[2],
    normal[0]*tangent[1] - normal[1]*tangent[0]
  };
  real64 records[12][4];
  for( int k = 0; k < count; ++k )
  {
    real64 u = 0.0, v = 0.0;
    for( int d = 0; d < 3; ++d )
    {
      real64 const delta = vertices[k][d] - center[d];
      u += delta * tangent[d];
      v += delta * bitangent[d];
    }
    records[k][0] = ::atan2(v, u);
    for( int d = 0; d < 3; ++d ) records[k][d+1] = vertices[k][d];
  }
  // Direct real-valued comparisons avoid truncating angular differences to
  // int (the old qsort comparator) and work in host/device builds.
  for( int k = 1; k < count; ++k )
  {
    real64 current[4];
    for( int d = 0; d < 4; ++d ) current[d] = records[k][d];
    int j = k;
    while( j > 0 && records[j-1][0] > current[0] )
    {
      for( int d = 0; d < 4; ++d ) records[j][d] = records[j-1][d];
      --j;
    }
    for( int d = 0; d < 4; ++d ) records[j][d] = current[d];
  }

  // Positive degree-four rule on the reference triangle. The weights
  // already sum to 1/2; no additional reference-area factor is needed.
  constexpr real64 barycentric[6][2] = {
    {0.445948490915965, 0.445948490915965},
    {0.445948490915965, 0.108103018168070},
    {0.108103018168070, 0.445948490915965},
    {0.091576213509771, 0.091576213509771},
    {0.091576213509771, 0.816847572980459},
    {0.816847572980459, 0.091576213509771}
  };
  constexpr real64 quadratureWeight[6] = {
    0.111690794839005, 0.111690794839005, 0.111690794839005,
    0.054975871827661, 0.054975871827661, 0.054975871827661
  };

  // A planar cut can coincide exactly with one of the Q1 knot planes.
  // In that case both half-space clips would contain the SAME 2D patch;
  // assign it to one side only to avoid double counting its area.
  bool knotPlane[3] = {};
  for( int d = 0; d < 3; ++d )
  {
    knotPlane[d] = true;
    for( int k = 0; k < count; ++k )
      if( ::fabs(records[k][d+1]) > coordinateTol ) knotPlane[d] = false;
  }

  real64 area = 0.0;
  for( int octant = 0; octant < 8; ++octant )
  {
    bool repeatedKnotPlane = false;
    for( int d = 0; d < 3; ++d )
      repeatedKnotPlane |= knotPlane[d] && (octant & (1 << d));
    if( repeatedKnotPlane ) continue;

    // Clip the original planar polygon at each Q1 knot plane. Within a
    // surviving octant N_I is a cubic polynomial in planar coordinates.
    real64 polygon[16][3] = {};
    int size = count;
    for( int k = 0; k < count; ++k )
      for( int d = 0; d < 3; ++d ) polygon[k][d] = records[k][d+1];
    for( int d = 0; d < 3 && size >= 3; ++d )
    {
      real64 clipped[16][3] = {};
      int newSize = 0;
      real64 const sign = (octant & (1 << d)) ? 1.0 : -1.0;
      for( int k = 0; k < size; ++k )
      {
        int const next = (k+1)%size;
        real64 const a = sign*polygon[k][d], b = sign*polygon[next][d];
        if( a >= 0.0 && newSize < 16 )
        {
          for( int j = 0; j < 3; ++j ) clipped[newSize][j] = polygon[k][j];
          ++newSize;
        }
        if( (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) )
        {
          real64 const t = a / (a-b);
          if( newSize < 16 )
          {
            for( int j = 0; j < 3; ++j )
              clipped[newSize][j] = polygon[k][j] + t*(polygon[next][j]-polygon[k][j]);
            ++newSize;
          }
        }
      }
      size = newSize;
      for( int k = 0; k < size; ++k )
        for( int j = 0; j < 3; ++j ) polygon[k][j] = clipped[k][j];
    }

    for( int k = 1; k+1 < size; ++k )
    {
      real64 const * a = polygon[0];
      real64 const * b = polygon[k];
      real64 const * c = polygon[k+1];
      real64 const u[3] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
      real64 const v[3] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
      real64 const cross[3] = {
        u[1]*v[2]-u[2]*v[1],
        u[2]*v[0]-u[0]*v[2],
        u[0]*v[1]-u[1]*v[0]
      };
      real64 const jacobian = ::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
      if( !(jacobian > 0.0) ) continue;

      for( int q = 0; q < 6; ++q )
      {
        real64 const eta = barycentric[q][0], xi = barycentric[q][1];
        real64 const w = 1.0-eta-xi;
        real64 shape = 1.0;
        for( int d = 0; d < 3; ++d )
        {
          real64 const x = eta*a[d] + xi*b[d] + w*c[d];
          shape *= 1.0-::fabs(x)/h[d];
        }
        area += jacobian * quadratureWeight[q] * shape;
      }
    }
  }
  return ::fmax(0.0, area);
}

} // namespace mpmMeasures
} // namespace geos

#endif
