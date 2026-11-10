/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */

/**
 * @file MPMAxisymmetricSurfaces.hpp
 * @brief Ring-weighted Q1 geometry and surface-reconstruction kernels for axisymmetric MPM.
 *
 * These host/device helpers use the meridional radius and axial coordinate;
 * geometric integration includes the full 2*pi*r measure. Mechanical particle
 * interpolation and its chosen CPDI/CPTI/B-spline stencil remain independent.
 */
#ifndef GEOS_SOLIDMECHANICS_MPMAXISYMMETRICSURFACES_HPP
#define GEOS_SOLIDMECHANICS_MPMAXISYMMETRICSURFACES_HPP
#include "common/DataTypes.hpp"
#include "mesh/ParticleGeometry.hpp"
#include <cmath>
#include <cfloat>
namespace geos
{
namespace axisymmetricSurface
{
// Full azimuthal angle. All geometric support integrals include this factor;
// normalized particle-to-node weights cancel it between numerator and denominator.
constexpr real64 twoPi = 6.283185307179586476925286766559;

/**
 * @brief Project an affine particle domain onto the axisymmetric geometric Q1 grid.
 *
 * Unlike mechanical CPDI/CPTI or B-spline interpolation, surface color and
 * surface-normal geometry always use a piecewise-bilinear Q1 basis on the
 * meridional (r,y) grid. Integrating the same Q1 basis over either a triangle
 * or a parallelogram avoids introducing a geometry-dependent surface stencil.
 *
 * The returned nodal contributions have the normalized ring-volume weight
 *
 *   S_I = integral_domain r N_I dA / (r_p A_p),
 *
 * and the axisymmetric geometric-gradient weights
 *
 *   B_{r,I} = integral_domain (r dN_I/dr + N_I) dA / (r_p A_p),
 *   B_{y,I} = integral_domain r dN_I/dy dA / (r_p A_p).
 *
 * The radial term includes the cylindrical N/r contribution. The common 2*pi
 * factor cancels in these ratios. The caller converts the dimensionless
 * weights to actual material volume and geometric-normal contributions.
 *
 * @param simplex Whether the affine domain is triangular (otherwise a quad).
 * @param center Particle center in meridional coordinates; center[0] is radius.
 * @param geometry Generators used to construct the current affine domain.
 * @param origin Local grid origin, including the cells accessible to this rank.
 * @param h Local grid spacing in the meridional directions.
 * @param cells Number of local cells in each direction.
 * @param emit Callback receiving (node_i, node_j, S, B_r, B_y).
 * @return False for nonpositive normalization or an out-of-grid vertex,
 *         before any callback is invoked; true after the cell integration.
 *
 * @note The domain and local cell map must use the same r-y coordinates.
 *       The denominator uses center radius times meridional domain measure;
 *       its validity is checked, rather than repairing degenerate geometry.
 */
template< class POINT, class MATRIX, class ORIGIN, class SPACING, class CELLS, class EMIT >
GEOS_HOST_DEVICE bool projectGeometry( bool const simplex,
                                       POINT const & center,
                                       MATRIX const & geometry,
                                       ORIGIN const & origin,
                                       SPACING const & h,
                                       CELLS const & cells,
                                       EMIT const & emit )
{
  using namespace particleGeometry;

  // Construct the native affine simplex/quad without changing its geometry.
  // For a centered affine domain, r_p*A_p is the meridional first radial moment.
  auto const domain = affine::Domain::fromGenerators( simplex ? Shape::Triangle : Shape::Parallelogram, geometry );
  real64 const denominator = center[0] * domain.measure(); // 2*pi cancels in normalized weights
  if( !(denominator > 0) )
    return false;

  // Reject any domain that extends beyond local-grid coverage before emitting
  // contributions. A small index-space tolerance accepts roundoff at faces.
  real64 vertices[8][3];
  domain.vertices( center, vertices );
  int const count = domain.verticesCount();
  int first[2], last[2];
  for( int d = 0; d < 2; ++d )
  {
    real64 lo = DBL_MAX, hi = -DBL_MAX;
    for( int v = 0; v < count; ++v )
    {
      lo = ::fmin( lo, vertices[v][d] );
      hi = ::fmax( hi, vertices[v][d] );
    }
    real64 const a = (lo - origin[d]) / h[d], b = (hi - origin[d]) / h[d];
    real64 const tolerance = 64 * DBL_EPSILON * (cells[d] + 1);
    if( a < -tolerance || b > cells[d] + tolerance )
      return false;

    first[d] = static_cast< int >( ::floor( ::fmax( 0.0, a ) ) );
    last[d] = static_cast< int >( ::floor( ::fmin( static_cast< real64 >( cells[d] ), b ) ) );
    if( last[d] >= cells[d] )
      last[d] = cells[d]-1;
  }

  // Three-point Gauss-Legendre quadrature in each unit-square coordinate.
  // The triangle's Duffy map supplies the (1-t) Jacobian separately below.
  constexpr real64 q[3] = { 0.11270166537925831148, 0.5, 0.88729833462074168852 };
  constexpr real64 w[3] = { 5.0 / 18, 4.0 / 9, 5.0 / 18 };
  for( int i = first[0]; i <= last[0]; ++i )
    for( int j = first[1]; j <= last[1]; ++j )
    {
      // Intersect this particle polygon with the current grid cell. Working
      // cell by cell is necessary because Q1 gradients change across faces.
      real64 polygon[10][2] = {};
      int size = count;
      for( int v = 0; v < count; ++v )
        for( int d = 0; d < 2; ++d )
          polygon[v][d] = vertices[v][d];

      real64 const cell[2] = { origin[0] + i * h[0], origin[1] + j * h[1] };

      // Four successive half-plane clips: low/high radius, low/high axial.
      // Sign-adjusted signed distances keep vertices on each cell face and
      // interpolate only when an edge strictly crosses that face.
      for( int face = 0; face < 4 && size > 0; ++face )
      {
        int const d = face / 2;
        real64 const sign = face % 2 == 0 ? -1 : 1;
        real64 const boundary = cell[d] + (face % 2) * h[d];
        real64 clipped[10][2] = {};
        int nextSize = 0;
        for( int v = 0; v < size; ++v )
        {
          int const next = (v + 1) % size;
          real64 const a = sign * (polygon[v][d] - boundary), b = sign * (polygon[next][d] - boundary);
          if( a <= 0 )
          {
            clipped[nextSize][0] = polygon[v][0];
            clipped[nextSize++][1] = polygon[v][1];
          }
          if( (a < 0 && b > 0) || (a > 0 && b < 0) )
          {
            real64 const t = a / (a - b);
            for( int k = 0; k < 2; ++k )
              clipped[nextSize][k] = polygon[v][k] + t * (polygon[next][k] - polygon[v][k]);

            ++nextSize;
          }
        }
        size = nextSize;
        for( int v = 0; v < size; ++v )
          for( int k = 0; k < 2; ++k )
            polygon[v][k] = clipped[v][k];
      }

      // Triangulate the clipped convex polygon as a fan about vertex zero.
      // Accumulate the four Q1-node weights before calling emit, so each
      // node receives at most one contribution from this particle/cell pair.
      real64 S[4] = {}, B[4][2] = {};
      for( int v = 1; v + 1 < size; ++v )
      {
        real64 const ax = polygon[v][0] - polygon[0][0], ay = polygon[v][1] - polygon[0][1];
        real64 const bx = polygon[v + 1][0] - polygon[0][0], by = polygon[v + 1][1] - polygon[0][1];
        real64 const jac = ::fabs( ax * by - ay * bx ) / denominator;
        for( int u = 0; u < 3; ++u )
          for( int z = 0; z < 3; ++z )
          {
            // Duffy coordinates span the reference triangle: (t,(1-t)q_z).
            // Multiply by the signed-free affine area Jacobian and the
            // radial coordinate to account for the revolved ring measure.
            real64 const t = q[u], s = (1 - t) * q[z];
            real64 const r = polygon[0][0] + t * ax + s * bx, y = polygon[0][1] + t * ay + s * by;
            real64 const xhat = (r - cell[0]) / h[0], yhat = (y - cell[1]) / h[1], weight = jac * w[u] * w[z] * (1 - t);
            for( int a = 0; a < 2; ++a )
              for( int b = 0; b < 2; ++b )
              {
                int const node = 2 * a + b;
                real64 const nx = a ? xhat : 1 - xhat, ny = b ? yhat : 1 - yhat;
                S[node] += weight * r * nx * ny;
                // r*dN/dr + N is the radial cylindrical-divergence term;
                // the axial term is r*dN/dy. Both use the same ring measure.
                B[node][0] += weight * (r * (a ? 1 : -1) / h[0] * ny + nx * ny);
                B[node][1] += weight * r * nx * (b ? 1 : -1) / h[1];
              }
          }
      }

      // Skip a node with no positive mapped ring volume; a valid geometric
      // projection can nevertheless have empty intersections with some cells.
      for( int a = 0; a < 2; ++a )
        for( int b = 0; b < 2; ++b )
          if( S[2 * a + b] > 0 )
            emit( i + a, j + b, S[2 * a + b], B[2 * a + b][0], B[2 * a + b][1] );
    }
  return true;
}

/**
 * @brief Gather a nodal interface plane into a particle-relative offset.
 *
 * A nodal surface position is a displacement from its node, whereas the
 * particle needs a displacement relative to its own center. Project the
 * node-to-particle displacement plus that nodal offset onto the normal.
 * Averaging the nodal offset vectors directly would bias a plane because
 * the ring-weighted mapping centroid is not the meridional area centroid.
 *
 * @param weight Particle-to-node interpolation weight.
 * @param node Absolute meridional node position.
 * @param particle Absolute meridional particle position.
 * @param normal Nodal surface normal; its magnitude need not be unity.
 * @param offset Nodal plane offset vector relative to the node.
 * @param result Accumulator for the particle-relative offset; not cleared.
 */
template< class NODE, class POINT, class NORMAL, class POSITION, class OUTPUT >
GEOS_HOST_DEVICE void gatherOffset( real64 const weight,
                                    NODE const & node,
                                    POINT const & particle,
                                    NORMAL const & normal,
                                    POSITION const & offset,
                                    OUTPUT && result )
{
  real64 const norm2 = normal[0] * normal[0] + normal[1] * normal[1];
  if( !(norm2 > 0) )
    return;

  // (n dot (node + offset - particle))/|n|^2 is the normal-coordinate
  // displacement. Multiplying back by n removes sensitivity to |n|.
  real64 const distance = ((node[0] - particle[0] + offset[0]) * normal[0] +
                           (node[1] - particle[1] + offset[1]) * normal[1]) / norm2;
  for( int d = 0; d < 2; ++d )
    result[d] += weight * distance * normal[d];
}

/**
 * @brief Describe one geometric Q1 nodal support in local meridional coordinates.
 *
 * The unclipped hat support occupies [-h_r,h_r] x [-h_y,h_y] relative
 * to a node at (radius,y). Bounds lo/hi may be shortened at the symmetry
 * axis and at physical symmetry faces. Every physical measure incorporates
 * the circumference 2*pi*(radius+x) once; no inactive third-grid direction
 * or user-specified plane-strain thickness enters these integrals.
 */
struct Support
{
  real64 h[2], radius, lo[2], hi[2];

  /**
   * @brief Evaluate the physical ring-measure density 2*pi*r*N_Q1.
   *
   * @param x Local radial displacement from the support node.
   * @param y Local axial displacement from the support node.
   * @return Weighted circumference at this point (before integration).
   * @note Callers integrate only within the clipped support lo/hi.
   */
  GEOS_HOST_DEVICE real64 weight( real64 const x, real64 const y ) const
  {
    return (1 - ::fabs( x ) / h[0] ) * (1 - ::fabs( y ) / h[1] ) * twoPi * (radius + x);
  }

  /**
   * @brief Integrate the full Q1-weighted physical ring volume of this support.
   *
   * Since the Q1 hat separates in r and y, the integral factors into a
   * radial first-moment integral and an axial hat integral. Splitting each
   * direction at its zero-offset knot makes both pieces polynomial and
   * exactly integrable; the final 2*pi converts meridional measure to volume.
   *
   * @return Nonnegative physical support volume, or zero for empty bounds.
   */
  GEOS_HOST_DEVICE real64 volume() const
  {
    if( !(hi[0] > lo[0] && hi[1] > lo[1]) )
      return 0;

    real64 radial = 0, axial = 0;
    for( int sign = -1; sign <= 1; sign += 2 )
    {
      // On either side of the knot, |x| = sign*x. The antiderivative
      // contains terms through x^3 because ring radius is radius+x.
      real64 const a = sign < 0 ? lo[0] : ::fmax( lo[0], 0.0 );
      real64 const b = sign < 0 ? ::fmin( hi[0], 0.0 ) : hi[0];
      if( b > a )
        radial += radius * (b - a) + (1 - sign * radius / h[0] ) * (b * b - a * a) / 2
                  - sign * (b * b * b - a * a * a) / (3 * h[0]);

      // The axial hat integral is linear on each side of y=0.
      real64 const c = sign < 0 ? lo[1] : ::fmax( lo[1], 0.0 );
      real64 const d = sign < 0 ? ::fmin( hi[1], 0.0 ) : hi[1];
      if( d > c )
        axial += d - c - sign * (d * d - c * c) / (2 * h[1]);
    }
    return ::fmax( 0.0, twoPi * radial * axial );
  }
};

/**
 * @brief Construct a Q1 surface-integration support clipped to physical bounds.
 *
 * The axis always clips local radial offsets at x=-radius, ensuring r>=0.
 * Additional clipping at lower/upper domain faces occurs only for faces
 * marked as symmetry boundaries; arbitrary non-symmetry boundaries do not
 * modify this geometric support here.
 *
 * @param hx Radial Q1 support half-width (grid spacing).
 * @param hy Axial Q1 support half-width (grid spacing).
 * @param radius Radius of the grid node.
 * @param y Axial coordinate of the grid node.
 * @param lower Lower physical bounds in r and y.
 * @param upper Upper physical bounds in r and y.
 * @param symmetry Four lower/upper symmetry-face flags, ordered by direction.
 */
template< class LOWER, class UPPER, class SYMMETRY >
GEOS_HOST_DEVICE Support support( real64 const hx,
                                  real64 const hy,
                                  real64 const radius,
                                  real64 const y,
                                  LOWER const & lower,
                                  UPPER const & upper,
                                  SYMMETRY const & symmetry )
{
  Support s{ { hx, hy }, radius, { -hx, -hy }, { hx, hy } };
  s.lo[0] = ::fmax( s.lo[0], -radius );
  real64 const node[2] = { radius, y };
  for( int d = 0; d < 2; ++d )
  {
    if( symmetry[2 * d] )
      s.lo[d] = ::fmax( s.lo[d], lower[d] - node[d] );

    if( symmetry[2 * d + 1] )
      s.hi[d] = ::fmin( s.hi[d], upper[d] - node[d] );
  }
  return s;
}

/**
 * @brief Integrate the ring-weighted Q1 support on one side of a local plane.
 *
 * Returns the physical material-side support volume
 *
 *   V(c) = integral_support 2*pi*(radius+x)*N_Q1(x,y)
 *                            * H(c - normal dot (x,y)) dx dy.
 *
 * Each of the four rectangles split at the Q1 hat knots is clipped by
 * normal dot (x,y) <= offset. A triangle fan partitions each clipped
 * polygon. Three-point Gauss integration in each Duffy coordinate handles
 * the polynomial integrand without subtracting nearly equal full volumes.
 *
 * @param s Clipped geometric support in node-relative coordinates.
 * @param normal Meridional plane normal; need not be unit length.
 * @param offset Signed plane constant using the same scale as normal.
 * @return Nonnegative physical volume on the included side of the plane.
 */
template< class NORMAL >
GEOS_HOST_DEVICE real64 halfSpaceVolume( Support const & s, NORMAL const & normal, real64 const offset )
{
  if( !(s.hi[0] > s.lo[0] && s.hi[1] > s.lo[1]) )
    return 0;

  constexpr real64 q[3] = { 0.11270166537925831148, 0.5, 0.88729833462074168852 };
  constexpr real64 w[3] = { 5.0 / 18, 4.0 / 9, 5.0 / 18 };
  real64 volume = 0;

  // Partition around x=0 and y=0: within each rectangle the absolute
  // values in the tensor-product Q1 hat have a fixed algebraic sign.
  for( int i = 0; i < 2; ++i )
    for( int j = 0; j < 2; ++j )
    {
      real64 const x0 = i == 0 ? s.lo[0] : ::fmax( 0.0, s.lo[0] );
      real64 const x1 = i == 0 ? ::fmin( 0.0, s.hi[0] ) : s.hi[0];
      real64 const y0 = j == 0 ? s.lo[1] : ::fmax( 0.0, s.lo[1] );
      real64 const y1 = j == 0 ? ::fmin( 0.0, s.hi[1] ) : s.hi[1];
      if( !(x1 > x0 && y1 > y0) )
        continue;

      // Sutherland-Hodgman clipping of the Q1 quadrant against the plane.
      // Keep points on the included side, then add an intersection whenever
      // successive corner values straddle the plane strictly.
      real64 const corners[4][2] = { { x0, y0 }, { x1, y0 }, { x1, y1 }, { x0, y1 } };
      real64 polygon[6][2] = {};
      int count = 0;
      for( int k = 0; k < 4; ++k )
      {
        int const next = (k + 1) % 4;
        real64 const a = normal[0] * corners[k][0] + normal[1] * corners[k][1] - offset;
        real64 const b = normal[0] * corners[next][0] + normal[1] * corners[next][1] - offset;
        if( a <= 0 )
        {
          polygon[count][0] = corners[k][0];
          polygon[count++][1] = corners[k][1];
        }
        if( (a < 0 && b > 0) || (a > 0 && b < 0) )
        {
          real64 const t = a / (a - b);
          polygon[count][0] = corners[k][0] + t * (corners[next][0] - corners[k][0]);
          polygon[count++][1] = corners[k][1] + t * (corners[next][1] - corners[k][1]);
        }
      }

      // Triangulate the clipped polygon about vertex zero. A Duffy map
      // provides a positive Jacobian factor (1-t) for each 3x3 quadrature
      // point, where Support::weight supplies the complete 2*pi*r factor.
      for( int k = 1; k + 1 < count; ++k )
      {
        real64 const ax = polygon[k][0] - polygon[0][0], ay = polygon[k][1] - polygon[0][1];
        real64 const bx = polygon[k + 1][0] - polygon[0][0], by = polygon[k + 1][1] - polygon[0][1];
        real64 const jac = ::fabs( ax * by - ay * bx );
        for( int u = 0; u < 3; ++u )
          for( int v = 0; v < 3; ++v )
          {
            real64 const t = q[u], z = (1 - t) * q[v];
            volume += jac * w[u] * w[v] * (1 - t) * s.weight( polygon[0][0] + t * ax + z * bx,
                                                             polygon[0][1] + t * ay + z * by );
          }
      }
    }
  return ::fmax( 0.0, volume );
}

/**
 * @brief Invert the clipped Q1 support-volume function for a plane offset.
 *
 * For a prescribed meridional normal, halfSpaceVolume is nondecreasing in
 * the plane constant. Bisect a symmetric interval for a fixed 48 iterations
 * to recover the interface offset whose enclosed ring volume approaches the
 * supplied nodal material volume.
 *
 * @param s Geometric support of the node, including axis/symmetry clipping.
 * @param normal Plane normal, in the same coordinates used for the support.
 * @param volume Target physical material-side support volume.
 * @return Approximate signed plane constant normal dot (x,y).
 *
 * @note The implemented initial bound uses the full Q1 half-widths h[]
 *       rather than re-extremizing the clipped support. This routine does
 *       not independently verify that the requested volume is bracketed.
 */
template< class NORMAL >
GEOS_HOST_DEVICE real64 offsetFromVolume( Support const & s, NORMAL const & normal, real64 const volume )
{
  real64 bound = ::fabs( normal[0] ) * s.h[0] + ::fabs( normal[1] ) * s.h[1];
  real64 a = -bound, b = bound;
  for( int iter = 0; iter < 48; ++iter )
  {
    real64 const mid = (a + b) / 2;
    if( halfSpaceVolume( s, normal, mid ) < volume )
      a = mid;
    else
      b = mid;
  }
  return (a + b) / 2;
}

/**
 * @brief Integrate the ring-weighted surface area of a meridional interface.
 *
 * The plane normal and its node-relative position specify a line in the
 * r-y cross-section. Revolving that line gives a surface with weighted area
 *
 *   A_I = integral_line 2*pi*r*N_Q1 ds.
 *
 * First intersect the line with the clipped rectangular support, then split
 * its parameter interval at x=0 and y=0 (the Q1 hat-function knots). On
 * each segment the integrand is cubic or lower, so two-point Gauss-Legendre
 * quadrature integrates it exactly for the given straight-line geometry.
 *
 * @param s Q1 support, including axis/symmetry truncation.
 * @param normal Unnormalized meridional interface normal.
 * @param position Node-relative point on the interface plane.
 * @return Nonnegative Q1-weighted physical surface area.
 */
template< class NORMAL, class POSITION >
GEOS_HOST_DEVICE real64 lineArea( Support const & s, NORMAL const & normal, POSITION const & position )
{
  real64 const norm = ::sqrt( normal[0] * normal[0] + normal[1] * normal[1] );
  if( !(norm > 0 && s.hi[0] > s.lo[0] && s.hi[1] > s.lo[1]) )
    return 0;

  // Normalize n and parameterize the line as center + t*tangent. Choosing
  // center = (n dot position)*n makes it the closest point to the node.
  real64 const n[2] = { normal[0] / norm, normal[1] / norm };
  real64 const distance = n[0] * position[0] + n[1] * position[1];
  real64 const center[2] = { distance * n[0], distance * n[1] }, tangent[2] = { -n[1], n[0] };
  real64 lo = -DBL_MAX, hi = DBL_MAX;
  for( int d = 0; d < 2; ++d )
  {
    // Intersect the line parameter with each rectangular coordinate bound.
    // If it is parallel to a face, reject it only when outside that slab.
    if( ::fabs( tangent[d] ) < DBL_EPSILON )
    {
      if( center[d] < s.lo[d] || center[d] > s.hi[d] )
        return 0;
    }
    else
    {
      real64 a = (s.lo[d] - center[d]) / tangent[d], b = (s.hi[d] - center[d]) / tangent[d];
      if( a > b )
      {
        real64 const tmp = a;
        a = b;
        b = tmp;
      }
      lo = ::fmax( lo, a );
      hi = ::fmin( hi, b );
    }
  }
  if( !(hi > lo) )
    return 0;

  // Split wherever the line crosses a Q1 knot. Sorting the at-most-four
  // cuts permits a common two-point quadrature on every smooth segment.
  real64 cuts[4] = { lo, hi, 0, 0 };
  int count = 2;
  for( int d = 0; d < 2; ++d )
    if( ::fabs( tangent[d] ) >= DBL_EPSILON )
    {
      real64 const t = -center[d] / tangent[d];
      if( t > lo && t < hi )
        cuts[count++] = t;
    }

  for( int i = 1; i < count; ++i )
    for( int j = i; j > 0 && cuts[j] < cuts[j-1]; --j )
    {
      real64 const tmp = cuts[j];
      cuts[j] = cuts[j-1];
      cuts[j-1] = tmp;
    }

  // Tangent has unit length, so dt is physical meridional arc length.
  // Support::weight already includes the complete ring circumference.
  real64 area = 0;
  for( int i = 1; i < count; ++i )
  {
    real64 const mid = (cuts[i] + cuts[i-1]) / 2, half = (cuts[i] - cuts[i-1]) / 2;
    for( int sign = -1; sign <= 1; sign += 2 )
    {
      real64 const t = mid + sign * half / ::sqrt( 3.0 );
      area += half * s.weight( center[0] + t * tangent[0], center[1] + t * tangent[1] );
    }
  }
  return ::fmax( 0.0, area );
}

/**
 * @brief Fit a meridional interface normal from two weighted particle classes.
 *
 * visit(callback) iterates dimensionless local samples
 * (x/hx, y/hy, label=-1/+1, area/(hx*hy)); callers can include reflected
 * axis samples without changing their scalar phase labels. Fit a signed
 * tanh classifier f=tanh((phi_r*x+phi_y*y+phi_0)/2) to the two classes by
 * weighted nonlinear least squares with a small slope penalty.
 *
 * A Gauss-Newton-like 3x3 normal system is solved by partially pivoted
 * Gauss-Jordan elimination. Backtracking accepts only a non-increasing
 * objective. The fitted slope is converted from normalized coordinates to
 * a unit normal in physical r-y coordinates; the intercept is not returned
 * because the interface offset is recovered separately from nodal volume.
 *
 * @param visit Repeatable callback supplier for weighted signed samples.
 * @param hx Radial coordinate scale (grid spacing).
 * @param hy Axial coordinate scale (grid spacing).
 * @param maxIter Maximum nonlinear fitting iterations.
 * @param tolerance Convergence tolerance for normal and intercept increments.
 * @param normal In/out meridional normal; supplies initial direction and is
 *               overwritten only when a nondegenerate fit is available.
 * @return False for one-sided data, a singular system or a degenerate slope;
 *         true for a finite normalized fitted slope.
 */
template< class VISIT >
GEOS_HOST_DEVICE bool logisticNormal( VISIT const & visit,
                                      real64 const hx,
                                      real64 const hy,
                                      int const maxIter,
                                      real64 const tolerance,
                                      real64 ( &normal )[3] )
{
  // Initialize the dimensionless slope so a physical normal has consistent
  // radial and axial scaling even for unequal meridional cell dimensions.
  // A weak L2 slope penalty regularizes the weighted classification fit.
  real64 const length = ::sqrt( hx * hy ), penalty = 1.e-7;
  real64 phi[3] = { normal[0] * hx / length, normal[1] * hy / length, 0 };
  real64 previous[2] = { normal[0], normal[1] };
  real64 classWeight[2] = {};
  visit( [&]( real64, real64, real64 label, real64 weight )
  {
    classWeight[label > 0 ? 1 : 0] += weight;
  } );
  if( !(classWeight[0] > 0 && classWeight[1] > 0) )
    return false;

  // The residual is label - tanh(linear_predictor/2). Compute this complete
  // weighted objective again for every line-search trial; sample storage is
  // owned by the supplied visit callback rather than duplicated locally.
  auto objective = [&]( real64 const ( &p )[3] )
  {
    real64 value = penalty * (p[0] * p[0] + p[1] * p[1]);
    visit( [&]( real64 x, real64 y, real64 label, real64 weight )
    {
      real64 const residual = label - ::tanh( (p[0] * x + p[1] * y + p[2]) / 2 );
      value += weight * residual * residual;
    } );
    return value;
  };

  for( int iter = 0; iter < maxIter; ++iter )
  {
    // Assemble the three-variable weighted Gauss-Newton system. Deriv is
    // d(tanh(eta/2))/d(eta), so each normal-matrix term is w*deriv^2*v*v^T;
    // the right-hand side is w*deriv*(label-f)*v minus slope regularization.
    real64 system[3][4] = { { penalty, 0, 0, -penalty * phi[0] },
                            { 0, penalty, 0, -penalty * phi[1] }, { 0, 0, 0, 0 } };
    visit( [&]( real64 x, real64 y, real64 label, real64 weight )
    {
      real64 const v[3] = { x, y, 1 }, f = ::tanh( (phi[0] * x + phi[1] * y + phi[2]) / 2 ), deriv = (1 - f * f) / 2;
      for( int i = 0; i < 3; ++i )
      {
        system[i][3] += weight * deriv * (label - f) * v[i];
        for( int j = 0; j < 3; ++j )
          system[i][j] += weight * deriv * deriv * v[i] * v[j];
      }
    } );

    // Solve the 3x3 augmented system with partial row pivoting and
    // Gauss-Jordan elimination. Scale the singularity threshold by the
    // largest diagonal entry rather than using a dimensional fixed cutoff.
    real64 matrixScale = 0;
    for( int i = 0; i < 3; ++i )
      matrixScale = ::fmax( matrixScale, ::fabs( system[i][i] ) );

    for( int k = 0; k < 3; ++k )
    {
      int pivot = k;
      for( int i = k + 1; i < 3; ++i )
        if( ::fabs( system[i][k] ) > ::fabs( system[pivot][k] ) )
          pivot = i;

      if( ::fabs( system[pivot][k] ) <= 64 * DBL_EPSILON * matrixScale )
        return false;

      for( int j = k; j < 4; ++j )
      {
        real64 const tmp = system[k][j];
        system[k][j] = system[pivot][j];
        system[pivot][j] = tmp;
      }
      real64 const diagonal = system[k][k];
      for( int j = k; j < 4; ++j )
        system[k][j] /= diagonal;

      for( int i = 0; i < 3; ++i )
        if( i != k )
        {
          real64 const factor = system[i][k];
          for( int j = k; j < 4; ++j )
            system[i][j] -= factor * system[k][j];
        }
    }

    // Accept the Gauss-Newton increment only when the full nonlinear
    // objective does not increase. If sixteen halvings fail, retain the
    // previous iterate and return its final direction below.
    real64 const oldObjective = objective( phi );
    real64 step = 1, trial[3] = {};
    for( int backtrack = 0; backtrack < 16; ++backtrack )
    {
      for( int i = 0; i < 3; ++i )
        trial[i] = phi[i] + step * system[i][3];

      if( objective( trial ) <= oldObjective )
        break;

      step /= 2;
    }
    if( objective( trial ) > oldObjective )
      break;

    for( int i = 0; i < 3; ++i )
      phi[i] = trial[i];

    // Convert the dimensionless fit slope back to physical coordinates.
    // Stop when both the normal direction and intercept update stabilize.
    real64 const n[2] = { phi[0] / hx, phi[1] / hy }, norm = ::sqrt( n[0] * n[0] + n[1] * n[1] );
    if( !(norm > 0 && norm < DBL_MAX) )
      return false;

    real64 const change = ::hypot( n[0] / norm - previous[0], n[1] / norm - previous[1] );
    previous[0] = n[0] / norm;
    previous[1] = n[1] / norm;
    if( iter > 2 && change < tolerance && step * ::fabs( system[2][3] ) < tolerance )
      break;
  }

  // Only the unit orientation is exported; the caller subsequently inverts
  // material volume on this normal to set a consistent plane position.
  real64 const norm = ::hypot( phi[0] / hx, phi[1] / hy );
  if( !(norm > 0 && norm < DBL_MAX) )
    return false;

  normal[0] = phi[0] / (hx * norm);
  normal[1] = phi[1] / (hy * norm);
  normal[2] = 0;
  return true;
}

} // namespace axisymmetricSurface
} // namespace geos
#endif
