/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */
#ifndef GEOS_MESH_PARTICLEREFERENCEGEOMETRY_HPP_
#define GEOS_MESH_PARTICLEREFERENCEGEOMETRY_HPP_

#include "common/DataTypes.hpp"
#include <cmath>

namespace geos
{
namespace particleGeometry
{

/** Geometry is independent of ParticleType (the interpolation algorithm).
 * Boxes use half edges; wedges/tetrahedra use full edges. The existing solver
 * R-vectors ALWAYS retain their old equivalent-box meaning. These generators
 * are separate reference data, not additional mechanical degrees of freedom.
 */
enum class Shape : integer { Parallelepiped = 0, TriangleWedge = 1, Tetrahedron = 2 };
constexpr char const * shapeKey() { return "particleReferenceGeometry"; }
constexpr char const * vectorsKey() { return "particleReferenceGeometryVectors"; }

GEOS_HOST_DEVICE inline int vertexCount( integer const shape )
{ return shape == 0 ? 8 : ( shape == 1 ? 6 : ( shape == 2 ? 4 : 0 ) ); }

// All matrix arguments store vectors in ROWS, like particleRVectors.
template< class MATRIX >
GEOS_HOST_DEVICE inline real64 determinant( MATRIX const & g )
{
  return g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])
       - g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])
       + g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]);
}

template< class MATRIX >
GEOS_HOST_DEVICE inline real64 signedMeasure( integer const shape, MATRIX const & g )
{
  return determinant( g ) * ( shape == 0 ? 8.0 : ( shape == 1 ? 0.5 : 1.0/6.0 ) );
}

template< class MATRIX, class FTYPE >
GEOS_HOST_DEVICE inline void deform( MATRIX const & reference, FTYPE const & F, real64 (&current)[3][3] )
{
  for( int a=0; a<3; ++a )
    for( int i=0; i<3; ++i )
    {
      current[a][i] = 0.0;
      for( int j=0; j<3; ++j ) current[a][i] += F[i][j]*reference[a][j];
    }
}

/** Return vertices in the geometry convention, independent of the backend:
 * standard hex signs; bottom triangle (origin, +edge1, +edge2), then its top;
 * tetrahedron (origin, +edge1, +edge2, +edge3). Connectivity permutations below
 * supply outward-facing VTK/Silo cells. Neighboring particles do not share nodes.
 */
template< class POINT, class MATRIX >
GEOS_HOST_DEVICE inline void vertices( integer const shape, POINT const & center,
                                       MATRIX const & g, real64 (&v)[8][3] )
{
  if( shape == 0 )
  {
    int const signs[8][3] = { {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
                             {-1,-1,1}, {1,-1,1}, {1,1,1}, {-1,1,1} };
    for( int k=0; k<8; ++k )
      for( int i=0; i<3; ++i )
        v[k][i] = center[i]+signs[k][0]*g[0][i]+signs[k][1]*g[1][i]+signs[k][2]*g[2][i];
  }
  else if( shape == 1 )
  {
    for( int i=0; i<3; ++i )
    {
      v[0][i] = center[i]-(g[0][i]+g[1][i])/3.0-g[2][i]/2.0;
      v[1][i] = v[0][i]+g[0][i]; v[2][i] = v[0][i]+g[1][i];
      for( int k=0; k<3; ++k ) v[k+3][i] = v[k][i]+g[2][i];
    }
  }
  else if( shape == 2 )
  {
    for( int i=0; i<3; ++i )
    {
      v[0][i] = center[i]-(g[0][i]+g[1][i]+g[2][i])/4.0;
      for( int k=0; k<3; ++k ) v[k+1][i] = v[0][i]+g[k][i];
    }
  }
}

/** Change the reference WITHOUT moving the current vertices.
 * Rnew = Rold * Fold^T * Fnew^{-T}. Only call for numerical/reference resets,
 * never for genuine deformation or a physical TransformParticles rotation.
 * Failure leaves reference unchanged. A singular reference reset must not
 * silently destroy native geometry. Scaling the inverse avoids overflow for
 * otherwise well-conditioned matrices with very small/large entries.
 */
template< class OLD_F, class NEW_F, class MATRIX >
GEOS_HOST_DEVICE inline bool rebase( OLD_F const & oldF, NEW_F const & newF, MATRIX && reference )
{
  real64 scale = 0.0;
  for( int i=0; i<3; ++i ) for( int j=0; j<3; ++j )
  {
    real64 const a = ::fabs( newF[i][j] );
    if( !std::isfinite( a ) ) return false;
    if( a > scale ) scale = a;
  }
  if( !( scale > 0.0 ) ) return false;
  real64 f[3][3];
  for( int i=0; i<3; ++i ) for( int j=0; j<3; ++j ) f[i][j] = newF[i][j]/scale;
  real64 const det = determinant( f );
  if( !std::isfinite( det ) || ::fabs( det ) < 1.0e-15 ) return false;
  real64 inv[3][3];
  for( int i=0; i<3; ++i ) for( int j=0; j<3; ++j )
    inv[i][j] = ( f[(j+1)%3][(i+1)%3]*f[(j+2)%3][(i+2)%3]
                 -f[(j+1)%3][(i+2)%3]*f[(j+2)%3][(i+1)%3] ) / det / scale;
  real64 current[3][3], next[3][3];
  deform( reference, oldF, current );
  for( int a=0; a<3; ++a ) for( int i=0; i<3; ++i )
  {
    next[a][i] = 0.0;
    for( int j=0; j<3; ++j ) next[a][i] += inv[i][j]*current[a][j];
    if( !std::isfinite( next[a][i] ) ) return false;
  }
  for( int a=0; a<3; ++a ) for( int i=0; i<3; ++i ) reference[a][i] = next[a][i];
  return true;
}

/** Output connectivity only. VTK linear wedges (through 9.6.2) require
 * (0,1,2) to point OUTWARD, away from the other triangle. Our positively
 * oriented generators point toward the other triangle, so reverse both faces.
 * Keep this backend convention out of the reusable reference geometry itself.
 */
GEOS_HOST_DEVICE inline int vtkVertex( integer const shape, int const i )
{
  int const wedge[6] = {0,2,1,3,5,4};
  return shape == 1 ? wedge[i] : i;
}

/** Canonical vertex at each Silo index. Compose the VTK wedge orientation
 * above with VisIt's TranslateSiloWedgeToVTKWedge permutation. Tetrahedra
 * also need an odd permutation because Silo and VTK orientations differ.
 */
GEOS_HOST_DEVICE inline int siloVertex( integer const shape, int const i )
{
  int const prism[6] = {5,2,0,3,4,1};
  int const tet[4] = {1,0,2,3};
  return shape == 1 ? prism[i] : ( shape == 2 ? tet[i] : i );
}

} // namespace particleGeometry
} // namespace geos
#endif
