/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */
#ifndef GEOS_MESH_PARTICLEGEOMETRY_HPP
#define GEOS_MESH_PARTICLEGEOMETRY_HPP

#include "common/DataTypes.hpp"
#include <cmath>
#include <cfloat>

namespace geos
{
namespace particleGeometry
{

// Geometry and interpolation are independent. The file schema is version 2. The 2D shapes are geometry-kernel support only;
// they do not enable a native-2D mesh or axisymmetric mechanics. CPTI uses 3D tetrahedra.
enum class Shape : integer
{
  Parallelepiped = 0, TriangleWedge = 1, Tetrahedron = 2,
  Parallelogram = 3, Triangle = 4
};

namespace affine
{

GEOS_HOST_DEVICE inline int dimension( Shape const shape )
{ return shape == Shape::Parallelogram || shape == Shape::Triangle ? 2 : 3; }

GEOS_HOST_DEVICE inline int vertexCount( Shape const shape )
{
  switch( shape )
  {
    case Shape::Parallelepiped: return 8;
    case Shape::TriangleWedge: return 6;
    case Shape::Tetrahedron: return 4;
    case Shape::Parallelogram: return 4;
    case Shape::Triangle: return 3;
  }
  return 0;
}

template< class MATRIX >
GEOS_HOST_DEVICE inline real64 determinant( MATRIX const & g )
{
  return g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])
       - g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])
       + g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0]);
}

// Generators are rows in storage: G_current = G_reference F^T.
// All generators in this interface are FULL edges, including box generators.
template< class MATRIX, class FTYPE, class OUTPUT >
GEOS_HOST_DEVICE inline void deform( MATRIX const & reference, FTYPE const & F, OUTPUT && current )
{
  for( int a=0; a<3; ++a ) for( int i=0; i<3; ++i )
    current[a][i] = reference[a][0]*F[i][0] + reference[a][1]*F[i][1] + reference[a][2]*F[i][2];
}

/** Stack-local affine geometry; never an additional per-particle field.
 * x_vertex = x_center + sum_a xi_vertex[a] generators[a]. The reference
 * vertices below are centered on the geometric centroid. Particle mass and
 * physical volume remain independent solver state, not inferred from a scaled
 * integration domain. No dynamic allocation or virtual dispatch is required.
 */
struct Domain
{
  Shape shape = Shape::Parallelepiped;
  real64 generators[3][3] = {};

  GEOS_HOST_DEVICE Domain() = default;

  template< class MATRIX >
  GEOS_HOST_DEVICE static Domain fromGenerators( Shape const shape, MATRIX const & g )
  {
    Domain result;
    result.shape = shape;
    for( int a=0; a<3; ++a ) for( int i=0; i<3; ++i ) result.generators[a][i] = g[a][i];
    return result;
  }

  GEOS_HOST_DEVICE int verticesCount() const { return vertexCount( shape ); }

  GEOS_HOST_DEVICE void vertexOffset( int const k, real64 (&offset)[3] ) const
  {
    real64 xi[3] = {};
    if( shape == Shape::Parallelepiped || shape == Shape::Parallelogram )
    {
      int const signs[8][3] = { {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
                               {-1,-1,1}, {1,-1,1}, {1,1,1}, {-1,1,1} };
      for( int a=0; a<dimension( shape ); ++a ) xi[a] = 0.5*signs[k][a];
    }
    else if( shape == Shape::Tetrahedron )
    {
      for( int a=0; a<3; ++a ) xi[a] = ( k == a+1 ? 1.0 : 0.0 ) - 0.25;
    }
    else
    {
      int const triangleVertex = k % 3;
      xi[0] = ( triangleVertex == 1 ? 1.0 : 0.0 ) - 1.0/3.0;
      xi[1] = ( triangleVertex == 2 ? 1.0 : 0.0 ) - 1.0/3.0;
      xi[2] = shape == Shape::TriangleWedge ? ( k < 3 ? -0.5 : 0.5 ) : 0.0;
    }
    for( int i=0; i<3; ++i )
      offset[i] = xi[0]*generators[0][i] + xi[1]*generators[1][i] + xi[2]*generators[2][i];
  }

  template< class POINT >
  GEOS_HOST_DEVICE void vertices( POINT const & center, real64 (&v)[8][3] ) const
  {
    for( int k=0; k<verticesCount(); ++k )
    {
      real64 offset[3];
      vertexOffset( k, offset );
      for( int i=0; i<3; ++i ) v[k][i] = center[i] + offset[i];
    }
  }

  GEOS_HOST_DEVICE real64 signedMeasure() const
  {
    if( dimension( shape ) == 2 )
      return ( generators[0][0]*generators[1][1] - generators[0][1]*generators[1][0] )
             * ( shape == Shape::Triangle ? 0.5 : 1.0 );
    real64 const factor = shape == Shape::TriangleWedge ? 0.5 : ( shape == Shape::Tetrahedron ? 1.0/6.0 : 1.0 );
    return factor*determinant( generators );
  }

  GEOS_HOST_DEVICE real64 measure() const { return ::fabs( signedMeasure() ); }

  GEOS_HOST_DEVICE real64 diameter() const
  {
    real64 largest = 0.0;
    // For centered boxes, opposite corners give the diameter. This is also
    // the cheap path used by the current CPDI scaling/splitting criteria.
    if( shape == Shape::Parallelepiped || shape == Shape::Parallelogram )
    {
      // Only one vertex from each opposite pair is needed.
      int const count = dimension( shape ) == 2 ? 2 : 4;
      for( int k=0; k<count; ++k )
      {
        real64 x[3]; vertexOffset( k, x );
        real64 const d2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        if( d2 > largest ) largest = d2;
      }
      return 2.0*::sqrt( largest );
    }
    for( int a=0; a<verticesCount(); ++a ) for( int b=a+1; b<verticesCount(); ++b )
    {
      real64 x[3], y[3]; vertexOffset( a, x ); vertexOffset( b, y );
      real64 d2 = 0.0;
      for( int i=0; i<3; ++i ) d2 += (x[i]-y[i])*(x[i]-y[i]);
      if( d2 > largest ) largest = d2;
    }
    return ::sqrt( largest );
  }

  /** Parameter t for the first boundary hit x=c+t*direction.
   * Map the ray into generator coordinates, then intersect the centered box
   * bounds or simplex barycentric half spaces. Direction need not be unit.
   * The wedge combines triangle half spaces with its two thickness planes.
   */
  template<class VECTOR>
  GEOS_HOST_DEVICE real64 rayDistance(VECTOR const & direction) const
  {
    real64 q[3]={};
    if(dimension(shape)==2)
    {
      real64 det=generators[0][0]*generators[1][1]-generators[0][1]*generators[1][0];
      if(!(::fabs(det)>0)) return 0;
      q[0]=(direction[0]*generators[1][1]-direction[1]*generators[1][0])/det;
      q[1]=(generators[0][0]*direction[1]-generators[0][1]*direction[0])/det;
    }
    else
    {
      real64 det=determinant(generators);
      if(!(::fabs(det)>0)) return 0;
      for(int a=0;a<3;++a) for(int i=0;i<3;++i)
        q[a]+=direction[i]*(generators[(a+1)%3][(i+1)%3]*generators[(a+2)%3][(i+2)%3]
                           -generators[(a+1)%3][(i+2)%3]*generators[(a+2)%3][(i+1)%3])/det;
    }
    real64 t=1e300;
    if(shape==Shape::Parallelepiped || shape==Shape::Parallelogram)
    {
      for(int a=0;a<dimension(shape);++a) if(::fabs(q[a])>0 && 0.5/::fabs(q[a])<t) t=0.5/::fabs(q[a]);
    }
    else
    {
      int n=shape==Shape::Tetrahedron ? 3 : 2;
      real64 c=1.0/(n+1), sum=0;
      for(int a=0;a<n;++a) { sum+=q[a]; if(q[a]<0 && -c/q[a]<t) t=-c/q[a]; }
      if(sum>0 && c/sum<t) t=c/sum;
      if(shape==Shape::TriangleWedge && ::fabs(q[2])>0 && 0.5/::fabs(q[2])<t) t=0.5/::fabs(q[2]);
    }
    return t<1e300 ? t : 0;
  }

  // Return one axis-bisected box child and its offset from the parent center.
  // mask selects generator directions; positiveMask selects each child's side.
  // The caller transfers mass/history and preserves the existing child order.
  GEOS_HOST_DEVICE bool boxChild( int const mask, int const positiveMask,
                                 Domain & child, real64 (&offset)[3] ) const
  {
    if( shape != Shape::Parallelepiped && shape != Shape::Parallelogram ) return false;
    child = *this;
    for( int i=0; i<3; ++i ) offset[i] = 0.0;
    for( int a=0; a<dimension( shape ); ++a ) if( mask & (1 << a) )
      for( int i=0; i<3; ++i )
      {
        child.generators[a][i] *= 0.5;
        offset[i] += ( positiveMask & (1 << a) ? 0.5 : -0.5 ) * child.generators[a][i];
      }
    return true;
  }
};

/** Constant barycentric gradients for a positively oriented tetrahedron.
 * G rows are edges from vertex 0. The other three vertex gradients are the
 * reciprocal edge vectors; vertex 0 closes their sum to zero. Normalize first
 * to avoid cubing the dimensional length scale in the determinant test.
 * Integration volume is geometric: physical particle volume is not an input.
 */
template< class MATRIX >
GEOS_HOST_DEVICE inline bool tetrahedronGradients( MATRIX const & g, real64 (&gradient)[4][3] )
{
  real64 scale = 0.0;
  for( int a=0; a<3; ++a ) for( int i=0; i<3; ++i )
  {
    real64 const value = ::fabs( g[a][i] );
    if( !( value <= DBL_MAX ) ) return false;
    if( value > scale ) scale = value;
  }
  if( !( scale > 0.0 ) ) return false;
  real64 normalized[3][3];
  for( int a=0; a<3; ++a ) for( int i=0; i<3; ++i ) normalized[a][i] = g[a][i]/scale;
  real64 const det = determinant( normalized );
  if( !( det > 0.0 ) ) return false;
  for( int i=0; i<3; ++i )
  {
    gradient[0][i] = 0.0;
    for( int a=0; a<3; ++a )
    {
      gradient[a+1][i] = ( ( normalized[(a+1)%3][(i+1)%3]*normalized[(a+2)%3][(i+2)%3]
                             - normalized[(a+1)%3][(i+2)%3]*normalized[(a+2)%3][(i+1)%3] ) / det ) / scale;
      if( !( ::fabs( gradient[a+1][i] ) <= DBL_MAX ) ) return false;
      gradient[0][i] -= gradient[a+1][i];
    }
    if( !( ::fabs( gradient[0][i] ) <= DBL_MAX ) ) return false;
  }
  return true;
}

enum class InputError { None, Shape, Volume, Generators, Measure, Wedge };

inline char const * inputErrorMessage( InputError const e )
{
  switch( e )
  {
    case InputError::None: return "valid particle geometry";
    case InputError::Shape: return "GeometryShape must be integer 0 (box), 1 (wedge), or 2 (tet)";
    case InputError::Volume: return "ParticleVolume must be finite and positive";
    case InputError::Generators: return "geometry generators must be finite with positive orientation";
    case InputError::Measure: return "initial geometry measure must agree with ParticleVolume";
    case InputError::Wedge: return "wedge must be an in-plane triangle extruded in z";
  }
  return "invalid particle geometry";
}

// One file contract: explicit physical volume and full generators. No inferred
// box, overloaded volume column, or second geometry definition is accepted.
template< class MATRIX >
inline InputError validateFileGeometry( real64 shapeValue, real64 volume, MATRIX const & g )
{
  if( !std::isfinite(shapeValue) || shapeValue<0 || shapeValue>2 || shapeValue > ::floor(shapeValue) ) return InputError::Shape;
  if( !std::isfinite(volume) || volume<=0 ) return InputError::Volume;
  Shape shape=static_cast<Shape>(static_cast<integer>(shapeValue));
  for(int a=0;a<3;++a) for(int i=0;i<3;++i) if(!std::isfinite(g[a][i])) return InputError::Generators;
  real64 v=Domain::fromGenerators(shape,g).signedMeasure();
  if(!std::isfinite(v) || v<=0) return InputError::Generators;
  if(::fabs(v-volume)>1e-9*volume) return InputError::Measure;
  if(shape==Shape::TriangleWedge)
  {
    real64 extent=0;
    for(int a=0;a<3;++a) for(int i=0;i<3;++i) if(::fabs(g[a][i])>extent) extent=::fabs(g[a][i]);
    if(::fabs(g[0][2])>1e-12*extent || ::fabs(g[1][2])>1e-12*extent ||
       ::fabs(g[2][0])>1e-12*extent || ::fabs(g[2][1])>1e-12*extent) return InputError::Wedge;
  }
  return InputError::None;
}

// Physical convection is independent of constitutive F resets. Use the same
// physical increment as the constitutive update, with an alias-safe temporary.
template< class MATRIX, class INCREMENT >
GEOS_HOST_DEVICE inline void advance( MATRIX && current, INCREMENT const & increment )
{
  real64 next[3][3]; deform(current,increment,next);
  for(int a=0;a<3;++a) for(int i=0;i<3;++i) current[a][i]=next[a][i];
}

template< class REFERENCE, class CURRENT, class FTYPE >
GEOS_HOST_DEVICE inline void reset( REFERENCE const & reference, CURRENT && current,
                                    FTYPE const & newF, bool const resetGeometry )
{
  if(resetGeometry) deform(reference,newF,current);
}

} // namespace affine

constexpr char const * shapeKey() { return "particleGeometryShape"; }
constexpr char const * vectorsKey() { return "particleGeometryReference"; }
constexpr char const * currentKey() { return "particleGeometryCurrent"; }
inline int vertexCount(integer s) { return s>=0 && s<=4 ? affine::vertexCount(static_cast<Shape>(s)) : 0; }
// Backend connectivity permutations are not geometric state.
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
