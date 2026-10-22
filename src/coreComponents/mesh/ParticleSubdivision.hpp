/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */
#ifndef GEOS_MESH_PARTICLESUBDIVISION_HPP
#define GEOS_MESH_PARTICLESUBDIVISION_HPP

#include "mesh/ParticleGeometry.hpp"

namespace geos
{
namespace particleGeometry
{
namespace subdivision
{

/** A child expressed in its parent's centered generator coordinates.
 * Generators are rows: G_child = A G_parent, dx_child = G_parent^T q.
 * Apply the SAME transform to G0 and G, even after an FOnly reset. The
 * fraction partitions physical mass/volume independently of geometry measure.
 * This is temporary topology-update data, never a persistent particle field.
 */
struct Child
{
  real64 matrix[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
  real64 center[3] = {};
  real64 fraction = 1.0;

  GEOS_HOST_DEVICE affine::Domain geometry( affine::Domain const & parent ) const
  {
    affine::Domain result;
    result.shape = parent.shape;
    for( int a = 0; a < 3; ++a ) for( int i = 0; i < 3; ++i )
    {
      result.generators[a][i] = 0.0;
      for( int b = 0; b < 3; ++b )
        result.generators[a][i] += matrix[a][b] * parent.generators[b][i];
    }
    return result;
  }

  GEOS_HOST_DEVICE void offset( affine::Domain const & parent, real64 (&dx)[3] ) const
  {
    for( int i = 0; i < 3; ++i )
    {
      dx[i] = 0.0;
      for( int a = 0; a < 3; ++a ) dx[i] += center[a] * parent.generators[a][i];
    }
  }
};

/** Compose grandchild-in-child with child-in-parent; avoid reconstructing F. */
GEOS_HOST_DEVICE inline Child compose( Child const & parent, Child const & local )
{
  Child result;
  result.fraction = parent.fraction * local.fraction;
  for( int a = 0; a < 3; ++a )
  {
    result.center[a] = parent.center[a];
    for( int b = 0; b < 3; ++b ) result.center[a] += local.center[b] * parent.matrix[b][a];
    for( int b = 0; b < 3; ++b )
    {
      result.matrix[a][b] = 0.0;
      for( int k = 0; k < 3; ++k ) result.matrix[a][b] += local.matrix[a][k] * parent.matrix[k][b];
    }
  }
  return result;
}

/** Bounded leaf set; only selected parents need to retain one on the host. */
struct Plan
{
  int count = 1;
  Child children[8];
  bool unresolved = false;
};

/** Retain the established box child order (positive side first, last axis fastest). */
GEOS_HOST_DEVICE inline Plan boxPlan( int const mask, int const dimensions )
{
  Plan plan;
  int directions = 0;
  for( int a = 0; a < dimensions; ++a ) if( mask & (1 << a) ) ++directions;
  plan.count = 1 << directions;
  for( int c = 0; c < plan.count; ++c )
  {
    int bit = directions;
    plan.children[c].fraction = 1.0 / plan.count;
    for( int a = 0; a < dimensions; ++a ) if( mask & (1 << a) )
    {
      --bit;
      plan.children[c].matrix[a][a] = 0.5;
      plan.children[c].center[a] = c & (1 << bit) ? -0.25 : 0.25;
    }
  }
  return plan;
}

/** Diameter in the active spatial dimensions; plane strain excludes thickness. */
GEOS_HOST_DEVICE inline real64 diameter( affine::Domain const & domain, int const dimensions )
{
  real64 largestSquared = 0.0;
  for( int a = 0; a < domain.verticesCount(); ++a )
    for( int b = a + 1; b < domain.verticesCount(); ++b )
    {
      real64 va[3], vb[3];
      domain.vertexOffset( a, va );
      domain.vertexOffset( b, vb );
      real64 lengthSquared = 0.0;
      for( int i = 0; i < dimensions; ++i ) lengthSquared += (va[i]-vb[i]) * (va[i]-vb[i]);
      if( lengthSquared > largestSquared ) largestSquared = lengthSquared;
    }
  return ::sqrt( largestSquared );
}

/** Bisect a triangle/tet edge in canonical coordinates. Replacing one endpoint
 * at a time by its midpoint preserves orientation and gives exact half measure.
 * For wedges, the same base-triangle split is applied to both end faces.
 */
GEOS_HOST_DEVICE inline Child simplexChild( int const dimensions, int const edgeA,
                                           int const edgeB, int const side )
{
  real64 vertices[4][3] = {};
  for( int v = 0; v <= dimensions; ++v ) for( int a = 0; a < dimensions; ++a )
    vertices[v][a] = (v == a+1 ? 1.0 : 0.0) - 1.0 / (dimensions+1);
  int const replaced = side == 0 ? edgeB : edgeA;
  for( int a = 0; a < dimensions; ++a )
    vertices[replaced][a] = 0.5 * (vertices[edgeA][a] + vertices[edgeB][a]);

  Child child;
  child.fraction = 0.5;
  for( int a = 0; a < dimensions; ++a )
  {
    for( int v = 0; v <= dimensions; ++v ) child.center[a] += vertices[v][a] / (dimensions+1);
    for( int b = 0; b < dimensions; ++b ) child.matrix[a][b] = vertices[a+1][b] - vertices[0][b];
  }
  return child;
}

/** Shape-preserving binary split; ties use canonical edge order. Longest edge
 * is evaluated on CURRENT geometry, not G0. A wedge may split its thickness
 * only in 3D; its representation remains an affine triangular prism.
 */
GEOS_HOST_DEVICE inline bool bisect( affine::Domain const & domain, int const dimensions,
                                    Child (&children)[2] )
{
  bool const box = domain.shape == Shape::Parallelepiped || domain.shape == Shape::Parallelogram;
  bool const tet = domain.shape == Shape::Tetrahedron;
  bool const wedge = domain.shape == Shape::TriangleWedge;
  if( !box && !tet && !wedge && domain.shape != Shape::Triangle ) return false;

  int direction = -1;
  int edgeA = 0, edgeB = 1;
  real64 longestSquared = -1.0;
  if( !box )
  {
    int const vertices = tet ? 4 : 3;
    for( int a = 0; a < vertices; ++a ) for( int b = a+1; b < vertices; ++b )
    {
      real64 va[3], vb[3];
      domain.vertexOffset( a, va );
      domain.vertexOffset( b, vb );
      real64 lengthSquared = 0.0;
      for( int i = 0; i < dimensions; ++i ) lengthSquared += (va[i]-vb[i]) * (va[i]-vb[i]);
      if( lengthSquared > longestSquared )
      {
        longestSquared = lengthSquared;
        edgeA = a; edgeB = b;
      }
    }
  }
  for( int a = 0; a < dimensions; ++a ) if( box || (wedge && a == 2) )
  {
    real64 lengthSquared = 0.0;
    for( int i = 0; i < dimensions; ++i ) lengthSquared += domain.generators[a][i] * domain.generators[a][i];
    if( lengthSquared > longestSquared ) { longestSquared = lengthSquared; direction = a; }
  }
  if( !(longestSquared > 0.0) ) return false;
  for( int side = 0; side < 2; ++side )
  {
    children[side] = Child{};
    if( direction >= 0 )
    {
      children[side].matrix[direction][direction] = 0.5;
      children[side].center[direction] = side == 0 ? 0.25 : -0.25;
      children[side].fraction = 0.5;
    }
    else children[side] = simplexChild( tet ? 3 : 2, edgeA, edgeB, side );
  }
  return true;
}

/** Refine the largest violating leaf, with at most 2^d leaves per pass.
 * A tet's diameter need not decrease on the first bisection. Always inspect
 * the resulting leaves and report unresolved size at the work limit.
 */
GEOS_HOST_DEVICE inline Plan sizePlan( affine::Domain const & parent, int const dimensions,
                                      real64 const maximumDiameter )
{
  Plan plan;
  if( !(maximumDiameter > 0.0) ) return plan;
  int const capacity = 1 << dimensions;
  while( true )
  {
    int selected = -1;
    real64 largest = maximumDiameter * (1.0 + 1.0e-12);
    for( int c = 0; c < plan.count; ++c )
    {
      real64 const size = diameter( plan.children[c].geometry( parent ), dimensions );
      if( size > largest ) { largest = size; selected = c; }
    }
    if( selected < 0 ) return plan;
    if( plan.count == capacity ) { plan.unresolved = true; return plan; }
    Child local[2];
    Child const leaf = plan.children[selected];
    if( !bisect( leaf.geometry( parent ), dimensions, local ) )
    { plan.unresolved = true; return plan; }
    plan.children[selected] = compose( leaf, local[0] );
    plan.children[plan.count++] = compose( leaf, local[1] );
  }
}

} // namespace subdivision
} // namespace particleGeometry
} // namespace geos
#endif
