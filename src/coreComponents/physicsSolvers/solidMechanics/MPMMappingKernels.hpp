/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */

/**
 * @file MPMMappingKernels.hpp
 * @brief Fixed-dimensional particle-to-grid interpolation kernels for Q1,
 *        cubic B-spline, CPDI, and CPTI particle geometries.
 *
 * These helpers produce node indices, interpolation weights, spatial
 * gradients, and (for axisymmetry) hoop terms. They do not transfer mass,
 * momentum, or stresses. The solver performs those operations only after
 * it has accepted a valid mapping.
 */
#ifndef GEOS_SOLIDMECHANICS_MPMMAPPINGKERNELS_HPP
#define GEOS_SOLIDMECHANICS_MPMMAPPINGKERNELS_HPP

#include "mesh/ParticleGeometry.hpp"
#include <type_traits>

namespace geos
{
namespace mpmMapping
{

/**
 * @brief Invoke a fixed-dimensional launch using a compile-time dimension tag.
 *
 * The solver selects the 2D or 3D specialization once on the host; kernels
 * dispatched by the callback can then size their stencils and unroll loops
 * without a run-time dimension branch for each particle.
 *
 * @tparam LAUNCH Callable accepting std::integral_constant<int, DIM>.
 * @param dimension Spatial dimension (2 selects 2D; all other values select 3D).
 * @param launch Callable forwarded to the selected specialization.
 * @return The selected call's result, preserving its declared return category.
 */
template< class LAUNCH >
decltype(auto) dispatch( int const dimension, LAUNCH && launch )
{
  if( dimension == 2 ) return launch( std::integral_constant< int, 2 >{} );
  return launch( std::integral_constant< int, 3 >{} );
}

/**
 * @brief Locate a normalized coordinate in an inclusive range of grid cells.
 *
 * Here u = (x - origin) / h, first and last are legal integer cell indices,
 * and fraction is the local position within the selected cell. Coordinates
 * within a floating-point tolerance of the outer support endpoints are
 * snapped to the nearest endpoint; coordinates genuinely outside fail.
 *
 * @return False for an empty cell range, a nonfinite/out-of-range coordinate,
 *         or an unsupported position; true with cell/fraction populated.
 */
GEOS_HOST_DEVICE inline bool cellCoordinate( real64 const u, int const first,
                                             int const last, int & cell, real64 & fraction )
{
  // Reject nonfinite or grossly escaped coordinates before converting to an
  // integer cell index. The broad check is also safe for NaN, because the
  // conjunctive interval comparison then evaluates false.
  if( last < first || !( u >= first - 1.0 && u <= last + 2.0 ) ) return false;

  // The admissible normalized support is [first, last+1], with a small
  // magnitude-dependent tolerance for coordinates on a numerical boundary.
  real64 const tolerance = 64 * DBL_EPSILON * ( ::fabs( u ) > 1 ? ::fabs( u ) : 1 );
  if( u < first - tolerance || u > last + 1.0 + tolerance ) return false;

  // Floor normally chooses the containing cell. At either extreme, retain
  // the endpoint shape function by snapping its cell fraction to 0 or 1.
  cell = static_cast< int >( ::floor( u ) );
  fraction = u - cell;
  if( cell < first ) { cell = first; fraction = 0; }
  if( cell > last ) { cell = last; fraction = 1; }
  return true;
}

/**
 * @brief Physical-face locations used for one-sided cubic interpolation.
 *
 * Entries are cell indices in the local patch coordinate system, not MPI
 * partition boundaries. An index of -1 disables closure on that face.
 * Buffer and ghost nodes remain allocated; closure changes interpolation
 * weights instead of deleting or overwriting grid degrees of freedom.
 */
struct CubicBoundary
{
  int lower[3] = { -1, -1, -1 };
  int upper[3] = { -1, -1, -1 };
};

/**
 * @brief Replace cubic B-spline exterior controls by a physical-face basis.
 *
 * A boundary control outside the physical domain is eliminated by linear
 * continuation from the first two physical controls: u[-1] = 2*u[0] - u[1],
 * with the corresponding construction at the upper face. The same linear
 * combination is applied to shape functions and their physical derivatives,
 * before taking tensor products. Grid ghost values are not overwritten.
 *
 * For the first cell, s in [0,1], the nonzero lower-face weights are
 * (1-s+s^3/6, s-s^3/3, s^3/6). They reproduce a linear field, sum to one,
 * retain the boundary control value at s=0, and meet the interior cubic basis
 * with C2 continuity. They do not constrain tangential velocity components.
 * If both faces bound the same cell, the routine uses a linear two-node basis.
 *
 * @param lower Apply the lower-face closure.
 * @param upper Apply the upper-face closure.
 * @param s Physical-cell fraction in [0,1].
 * @param invH Inverse physical grid spacing used to scale derivatives.
 * @param weight Four stencil weights ordered by increasing grid-node index.
 * @param derivative Corresponding derivatives with respect to position.
 */
GEOS_HOST_DEVICE inline void closeCubicBoundary( bool const lower, bool const upper,
                                                 real64 const s, real64 const invH,
                                                 real64 ( &weight )[4], real64 ( &derivative )[4] )
{
  // Evaluate the combined basis polynomials directly instead of subtracting
  // nearly equal cardinal weights. Near r=0 that subtraction could lose a
  // physically small O(r/h) positive weight through roundoff.
  if( lower && upper )
  {
    // The physical cell spans both faces: only its two endpoint controls
    // remain, with constant derivatives and zero exterior support.
    weight[0] = 0.0;
    weight[1] = 1 - s;
    weight[2] = s;
    weight[3] = 0.0;
    derivative[0] = 0.0;
    derivative[1] = -invH;
    derivative[2] = invH;
    derivative[3] = 0.0;
  }
  else if( lower )
  {
    // At the lower face, slot 0 is the eliminated exterior control. The
    // remaining three polynomials and their derivatives form the closed basis.
    real64 const s2 = s * s, s3 = s2 * s;
    weight[0] = 0.0;
    weight[1] = 1 - s + s3 / 6;
    weight[2] = s - s3 / 3;
    weight[3] = s3 / 6;
    derivative[0] = 0.0;
    derivative[1] = ( -1 + s2 / 2 ) * invH;
    derivative[2] = ( 1 - s2 ) * invH;
    derivative[3] = ( s2 / 2 ) * invH;
  }
  else if( upper )
  {
    // Reflect the same closure about the upper face. t measures distance
    // inward from the upper boundary, and slot 3 is the exterior control.
    real64 const t = 1 - s, t2 = t * t, t3 = t2 * t;
    weight[0] = t3 / 6;
    weight[1] = t - t3 / 3;
    weight[2] = 1 - t + t3 / 6;
    weight[3] = 0.0;
    derivative[0] = ( -t2 / 2 ) * invH;
    derivative[1] = ( -1 + t2 ) * invH;
    derivative[2] = ( 1 - t2 / 2 ) * invH;
    derivative[3] = 0.0;
  }
}

/**
 * @brief Map one point to Q1 or cubic B-spline nodal interpolation data.
 *
 * Each spatial coordinate is located within a legal cell, and its 1D basis
 * weights/derivatives are formed before assembling the tensor product.
 * Q1 uses 2^DIM nodes (4 in 2D, 8 in 3D); cubic uses 4^DIM nodes (16 or 64).
 * The physical grid map supplies each node's local index. Gradients have
 * three components even for DIM=2, where the inactive z component is zero.
 *
 * For interior stencils, the final entry receives the arithmetic residual
 * needed to enforce partition of unity and zero summed gradients. Boundary-
 * closed cubic stencils instead retain small near-axis values by assigning
 * those residuals to their largest weight.
 *
 * @tparam DIM Compile-time active grid dimension (2 or 3).
 * @tparam CUBIC True for cubic B-splines, false for Q1.
 * @param boundary Optional physical-face indices used only for cubic closure.
 * @return False when spacing or support location is invalid; true on success.
 *         The caller must use output arrays only after a true return.
 */
template< int DIM, bool CUBIC, class POINT, class MAP, class NODES, class VALUES, class GRADIENTS >
GEOS_HOST_DEVICE bool point( POINT const & position, MAP const & map,
                            real64 const (&origin)[3], real64 const (&spacing)[3],
                            localIndex const (&cells)[3], NODES && nodes,
                            VALUES && values, GRADIENTS && gradients,
                            CubicBoundary const & boundary = CubicBoundary{} )
{
  // Compile-time stencil sizes avoid variable-length particle temporaries.
  // base[d] is the first grid-node index touched by the 1D support in d.
  constexpr int WIDTH = CUBIC ? 4 : 2;
  constexpr int COUNT = CUBIC ? ( 1 << ( 2 * DIM ) ) : ( 1 << DIM );
  int base[DIM];
  bool closed = false;
  real64 weight[DIM][WIDTH], derivative[DIM][WIDTH];
  // Reject a point unless every active direction has enough allocated
  // support nodes; cubic uses a one-cell wider footprint on either side.
  for( int d = 0; d < DIM; ++d )
  {
    if( !( spacing[d] > 0 && spacing[d] <= DBL_MAX ) ) return false;
    real64 s;
    int cell;
    if( !cellCoordinate( ( position[d] - origin[d] ) / spacing[d], CUBIC ? 1 : 0,
                         cells[d] - ( CUBIC ? 2 : 1 ), cell, s ) ) return false;
    base[d] = cell - ( CUBIC ? 1 : 0 );
    real64 const invH = 1 / spacing[d];
    if constexpr( CUBIC )
    {
      // Cardinal cubic weights and spatial derivatives in this direction.
      // The four entries correspond to node indices cell-1 through cell+2.
      real64 const t = 1 - s, s2 = s * s, s3 = s2 * s;
      weight[d][0] = t*t*t / 6;
      weight[d][1] = ( 3*s3 - 6*s2 + 4 ) / 6;
      weight[d][2] = ( -3*s3 + 3*s2 + 3*s + 1 ) / 6;
      weight[d][3] = s3 / 6;
      derivative[d][0] = -0.5*t*t*invH;
      derivative[d][1] = ( 1.5*s2 - 2*s )*invH;
      derivative[d][2] = ( -1.5*s2 + s + 0.5 )*invH;
      derivative[d][3] = 0.5*s2*invH;
      // Closure is driven by a physical domain face, not by the local MPI
      // partition edge. It modifies these 1D bases before tensor products.
      bool const lower = cell == boundary.lower[d];
      bool const upper = cell + 1 == boundary.upper[d];
      if( lower || upper )
      {
        closed = true;
        // Subtract the physical cell origin before scaling. At the axis this
        // preserves tiny positive radii lost by (r-ghostOrigin)/h - cell.
        s = ( position[d] - ( origin[d] + cell * spacing[d] ) ) * invH;
        s = s < 0.0 ? 0.0 : ( s > 1.0 ? 1.0 : s );
        closeCubicBoundary( lower, upper, s, invH, weight[d], derivative[d] );
      }
    }
    else
    {
      // Linear Q1 shape functions have two endpoint controls. Face closure
      // applies only to the cubic formulation.
      weight[d][0] = 1 - s;
      weight[d][1] = s;
      derivative[d][0] = -invH;
      derivative[d][1] = invH;
    }
  }
  // Enumerate tensor-product tuples in flattened stencil order (last
  // dimension varies fastest). For native 2D, index[2] remains zero so the
  // three-index grid map still addresses the single inactive plane.
  real64 sum = 0, gradientSum[DIM] = {};
  int largest = 0;
  for( int n = 0; n < COUNT; ++n )
  {
    int digit[DIM], index[3] = {}, remaining = n;
    for( int d = DIM - 1; d >= 0; --d )
    {
      digit[d] = remaining % WIDTH;
      remaining /= WIDTH;
      index[d] = base[d] + digit[d];
    }
    nodes[n] = map[index[0]][index[1]][index[2]];

    // Tensor products provide N_I. For an unclosed stencil, use the final
    // coefficient as a residual to cancel the accumulated roundoff.
    real64 value = 1;
    for( int d = 0; d < DIM; ++d ) value *= weight[d][digit[d]];
    values[n] = !closed && n == COUNT - 1 ? 1 - sum : value;
    if( values[n] > values[largest] ) largest = n;
    sum += values[n];
    if constexpr( DIM == 2 ) gradients[n][2] = 0;

    // dN_I/dx_d is the derivative factor in d times the ordinary weights
    // in every other direction. Correct the last interior entry so these
    // gradients sum to zero despite roundoff in the tensor products.
    for( int d = 0; d < DIM; ++d )
    {
      real64 gradient = derivative[d][digit[d]];
      for( int k = 0; k < DIM; ++k ) if( k != d ) gradient *= weight[k][digit[k]];
      gradients[n][d] = !closed && n == COUNT - 1 ? -gradientSum[d] : gradient;
      gradientSum[d] += gradients[n][d];
    }
  }
  if( closed )
  {
    // Closed faces have exact zero weights in their exterior slots, which
    // must remain zero. Apply roundoff corrections to the largest positive
    // contribution instead of disturbing small near-axis support values.
    values[largest] += 1 - sum;
    for( int d = 0; d < DIM; ++d ) gradients[largest][d] -= gradientSum[d];
  }
  return true;
}

/**
 * @brief Construct affine CPDI/CPTI mapping by Q1 evaluation at domain vertices.
 *
 * The generator matrix defines one temporary parallelogram/parallelepiped
 * (CPDI) or triangle/tetrahedron (CPTI). For each material-domain vertex,
 * evaluate the background-grid Q1 shape functions and combine them with the
 * material basis's geometrically determined gradient. Box gradients use
 * averages of the corner basis; simplex gradients are constant barycentric
 * gradients. The geometry used for positions and derivatives is identical.
 *
 * Each of the CORNERS vertices contributes SUPPORT=2^DIM entries to a flat
 * output stencil, even when two vertices address the same grid node. The
 * caller may later coalesce identical nodes. Physical volume/thickness does
 * not enter this geometrical mapping calculation.
 *
 * AXISYMMETRIC selects meridional ring-moment weights and a hoop-gradient
 * channel, while still using the two-dimensional Q1 background grid.
 *
 * @tparam DIM Active spatial dimension (2 or 3).
 * @tparam SIMPLEX Triangle/tetrahedron for true; affine box for false.
 * @tparam AXISYMMETRIC Use radial ring weighting for DIM=2.
 * @return False for a degenerate generator matrix, an invalid axisymmetric
 *         center radius, or an unmappable domain vertex. Some earlier
 *         vertices may already have written output on failure; discard it.
 */
template< int DIM, bool SIMPLEX, bool AXISYMMETRIC=false, class POINT, class MATRIX, class MAP,
          class NODES, class VALUES, class GRADIENTS >
GEOS_HOST_DEVICE bool domain( POINT const & center, MATRIX const & generators,
                             MAP const & map, real64 const (&origin)[3],
                             real64 const (&spacing)[3], localIndex const (&cells)[3],
                             NODES && nodes, VALUES && values, GRADIENTS && gradients )
{
  // Axisymmetric mapping lives entirely in the (r,y) plane, and normalized
  // ring moments below divide by the strictly positive particle-center radius.
  static_assert( !AXISYMMETRIC || DIM == 2, "Axisymmetric mapping is meridional" );
  if constexpr( AXISYMMETRIC )
    if( !( center[0] > 0 && center[0] <= DBL_MAX ) ) return false;

  // Each material-domain vertex samples a complete Q1 background-grid stencil.
  constexpr int CORNERS = SIMPLEX ? DIM + 1 : ( 1 << DIM );
  constexpr int SUPPORT = 1 << DIM;
  constexpr auto SHAPE = DIM == 2
    ? ( SIMPLEX ? particleGeometry::Shape::Triangle : particleGeometry::Shape::Parallelogram )
    : ( SIMPLEX ? particleGeometry::Shape::Tetrahedron : particleGeometry::Shape::Parallelepiped );
  // Reuse this one affine domain for its vertex offsets and reference-basis
  // gradients. The reciprocal arrays also supply the box generator inverses.
  auto const geometry = particleGeometry::affine::Domain::fromGenerators( SHAPE, generators );
  real64 reciprocal[DIM + 1][DIM];
  if constexpr( DIM == 2 )
  {
    if( !particleGeometry::affine::triangleGradients( generators, reciprocal ) ) return false;
  }
  else
  {
    if( !particleGeometry::affine::tetrahedronGradients( generators, reciprocal ) ) return false;
  }
  // Domain vertices carry material-basis data. A failure at any vertex is a
  // failure of the entire domain mapping, even if earlier outputs were filled.
  for( int a = 0; a < CORNERS; ++a )
  {
    real64 offset[3], vertex[DIM];
    geometry.vertexOffset( a, offset );
    for( int d = 0; d < DIM; ++d ) vertex[d] = center[d] + offset[d];
    // Evaluate nodal Q1 weights at this physical domain vertex. The point
    // gradients dw are not used for CPDI/CPTI: domain gradients come from the
    // derivatives of the material-domain shape functions themselves.
    localIndex support[SUPPORT];
    real64 w[SUPPORT], dw[SUPPORT][3];
    if( !point< DIM, false >( vertex, map, origin, spacing, cells, support, w, dw ) ) return false;

    real64 basisGradient[DIM] = {};
    if constexpr( SIMPLEX )
    {
      // A linear simplex's barycentric gradient is constant over its volume,
      // so the vertex's material gradient is one row of the reciprocal map.
      for( int d = 0; d < DIM; ++d ) basisGradient[d] = reciprocal[a][d];
    }
    else
    {
      // Average the affine box basis derivative at its vertices. The signs
      // follow the geometry's cyclic XY ordering, then its two Z planes;
      // reciprocal[k+1] expresses generator-direction derivatives in space.
      int const signs[8][3] = { { -1, -1, -1 }, { 1, -1, -1 }, { 1, 1, -1 }, { -1, 1, -1 },
                                { -1, -1, 1 }, { 1, -1, 1 }, { 1, 1, 1 }, { -1, 1, 1 } };
      for( int d = 0; d < DIM; ++d )
        for( int k = 0; k < DIM; ++k )
          basisGradient[d] += signs[a][k] * reciprocal[k + 1][d] / ( 1 << ( DIM - 1 ) );
    }
    // Cartesian integration assigns each vertex the same reference measure.
    // For a ring, the radial first moment changes the normalized corner
    // weight by offset_r/(12*r_center). The 2*pi factor cancels in this ratio.
    real64 cornerWeight = 1.0 / CORNERS;
    if constexpr( AXISYMMETRIC )
    {
      cornerWeight += offset[0] / ( 12 * center[0] );
      if constexpr( !SIMPLEX )
      {
        // For a meridional box, r-weighted spatial derivatives contain
        // cross terms when either generator has a radial component. Constant
        // simplex barycentric gradients need no such correction.
        int const signs[4][2] = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } };
        int const s = signs[a][0], t = signs[a][1];
        for( int d = 0; d < 2; ++d )
          basisGradient[d] = s * reciprocal[1][d] * ( 0.5 + t * generators[1][0] / ( 12 * center[0] ) )
                             + t * reciprocal[2][d] * ( 0.5 + s * generators[0][0] / ( 12 * center[0] ) );
      }
    }
    // Combine material vertex a with each of its Q1 grid nodes. Store all
    // vertex-node entries, including repeated node IDs, in a fixed ordering.
    // Per-vertex residuals enforce the intended sums without a separate pass.
    real64 sum = 0, gradientSum[DIM] = {}, hoopSum = 0;
    for( int b = 0; b < SUPPORT; ++b )
    {
      int const n = a * SUPPORT + b;
      nodes[n] = support[b];
      values[n] = b == SUPPORT - 1 ? cornerWeight - sum : w[b] * cornerWeight;
      sum += values[n];
      if constexpr( DIM == 2 ) gradients[n][2] = 0;
      if constexpr( AXISYMMETRIC )
      {
        // The third gradient slot represents the axisymmetric hoop term,
        // not an inactive z derivative. Each material corner contributes
        // 1/(CORNERS*r_center), distributed by its Q1 weights.
        real64 const hoopWeight = 1.0 / ( CORNERS * center[0] );
        gradients[n][2] = b == SUPPORT - 1 ? hoopWeight - hoopSum : w[b] * hoopWeight;
        hoopSum += gradients[n][2];
      }
      // Meridional/Cartesian gradients sum to the material-corner gradient;
      // the final Q1 node absorbs roundoff just as it does for the weights.
      for( int d = 0; d < DIM; ++d )
      {
        gradients[n][d] = b == SUPPORT - 1 ? basisGradient[d] - gradientSum[d] : w[b] * basisGradient[d];
        gradientSum[d] += gradients[n][d];
      }
    }
  }
  return true;
}

} // namespace mpmMapping
} // namespace geos
#endif
