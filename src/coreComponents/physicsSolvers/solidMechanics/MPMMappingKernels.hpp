/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */
#ifndef GEOS_SOLIDMECHANICS_MPMMAPPINGKERNELS_HPP
#define GEOS_SOLIDMECHANICS_MPMMAPPINGKERNELS_HPP

#include "mesh/ParticleGeometry.hpp"
#include <type_traits>

namespace geos
{
namespace mpmMapping
{

/** The timestep is shared. Only the host chooses which fixed-dimension launch
 * to instantiate; particle loops never test the simulation's dimension.
 */
template< class LAUNCH >
decltype(auto) dispatch( int const dimension, LAUNCH && launch )
{
  if( dimension == 2 ) return launch( std::integral_constant<int,2>{} );
  return launch( std::integral_constant<int,3>{} );
}

/** Check finite coordinates before conversion to an integer. Roundoff at an
 * outer support boundary is clamped, but genuinely escaped particles fail.
 */
GEOS_HOST_DEVICE inline bool cellCoordinate( real64 const u, int const first,
                                             int const last, int & cell, real64 & fraction )
{
  if( last < first || !(u >= first-1.0 && u <= last+2.0) ) return false;
  real64 const tolerance = 64*DBL_EPSILON*(::fabs(u)>1 ? ::fabs(u) : 1);
  if( u < first-tolerance || u > last+1.0+tolerance ) return false;
  cell = static_cast<int>(::floor(u));
  fraction = u-cell;
  if( cell < first ) { cell=first; fraction=0; }
  if( cell > last ) { cell=last; fraction=1; }
  return true;
}

/** Tensor-product Q1 or cubic B-spline support: 4/8 or 16/64 entries.
 * Local arrays and loops use DIM. The shared framework output gradient keeps
 * three components, with an explicitly zero inactive z component in 2D.
 */
template< int DIM, bool CUBIC, class POINT, class MAP, class NODES, class VALUES, class GRADIENTS >
GEOS_HOST_DEVICE bool point( POINT const & position, MAP const & map,
                            real64 const (&origin)[3], real64 const (&spacing)[3],
                            localIndex const (&cells)[3], NODES && nodes,
                            VALUES && values, GRADIENTS && gradients )
{
  constexpr int WIDTH = CUBIC ? 4 : 2;
  constexpr int COUNT = CUBIC ? (1 << (2*DIM)) : (1 << DIM);
  int base[DIM];
  real64 weight[DIM][WIDTH], derivative[DIM][WIDTH];
  for( int d=0; d<DIM; ++d )
  {
    if( !(spacing[d]>0 && spacing[d]<=DBL_MAX) ) return false;
    real64 s;
    int cell;
    if( !cellCoordinate((position[d]-origin[d])/spacing[d], CUBIC ? 1 : 0,
                         cells[d]-(CUBIC ? 2 : 1), cell, s) ) return false;
    base[d]=cell-(CUBIC ? 1 : 0);
    real64 const invH=1/spacing[d];
    if constexpr( CUBIC )
    {
      real64 const t=1-s, s2=s*s, s3=s2*s;
      weight[d][0]=t*t*t/6; weight[d][1]=(3*s3-6*s2+4)/6;
      weight[d][2]=(-3*s3+3*s2+3*s+1)/6; weight[d][3]=s3/6;
      derivative[d][0]=-0.5*t*t*invH; derivative[d][1]=(1.5*s2-2*s)*invH;
      derivative[d][2]=(-1.5*s2+s+0.5)*invH; derivative[d][3]=0.5*s2*invH;
    }
    else
    {
      weight[d][0]=1-s; weight[d][1]=s;
      derivative[d][0]=-invH; derivative[d][1]=invH;
    }
  }
  real64 sum=0, gradientSum[DIM]={};
  for( int n=0; n<COUNT; ++n )
  {
    int digit[DIM], index[3]={}, remaining=n;
    for( int d=DIM-1; d>=0; --d )
    {
      digit[d]=remaining%WIDTH; remaining/=WIDTH; index[d]=base[d]+digit[d];
    }
    nodes[n]=map[index[0]][index[1]][index[2]];
    real64 value=1;
    for( int d=0; d<DIM; ++d ) value*=weight[d][digit[d]];
    values[n]=n==COUNT-1 ? 1-sum : value;
    sum+=values[n];
    if constexpr( DIM==2 ) gradients[n][2]=0;
    for( int d=0; d<DIM; ++d )
    {
      real64 gradient=derivative[d][digit[d]];
      for( int k=0; k<DIM; ++k ) if(k!=d) gradient*=weight[k][digit[k]];
      gradients[n][d]=n==COUNT-1 ? -gradientSum[d] : gradient;
      gradientSum[d]+=gradients[n][d];
    }
  }
  return true;
}

/** Affine CPDI/CPTI uses the same temporary domain for vertex positions and
 * geometric derivatives. Physical volume and thickness do not enter mapping.
 * CPDI: averaged Q1 corner gradients of the material box basis.
 * CPTI: constant simplex barycentric gradients (triangle or tetrahedron).
 */
template< int DIM, bool SIMPLEX, bool AXISYMMETRIC=false, class POINT, class MATRIX, class MAP,
          class NODES, class VALUES, class GRADIENTS >
GEOS_HOST_DEVICE bool domain( POINT const & center, MATRIX const & generators,
                             MAP const & map, real64 const (&origin)[3],
                             real64 const (&spacing)[3], localIndex const (&cells)[3],
                             NODES && nodes, VALUES && values, GRADIENTS && gradients )
{
  static_assert(!AXISYMMETRIC || DIM==2, "Axisymmetric mapping is meridional");
  if constexpr(AXISYMMETRIC) if(!(center[0]>0 && center[0]<=DBL_MAX)) return false;
  constexpr int CORNERS = SIMPLEX ? DIM+1 : (1<<DIM);
  constexpr int SUPPORT = 1<<DIM;
  constexpr auto SHAPE = DIM==2
    ? (SIMPLEX ? particleGeometry::Shape::Triangle : particleGeometry::Shape::Parallelogram)
    : (SIMPLEX ? particleGeometry::Shape::Tetrahedron : particleGeometry::Shape::Parallelepiped);
  auto const geometry=particleGeometry::affine::Domain::fromGenerators(SHAPE,generators);
  real64 reciprocal[DIM+1][DIM];
  if constexpr( DIM==2 )
  {
    if(!particleGeometry::affine::triangleGradients(generators,reciprocal)) return false;
  }
  else
  {
    if(!particleGeometry::affine::tetrahedronGradients(generators,reciprocal)) return false;
  }
  for( int a=0; a<CORNERS; ++a )
  {
    real64 offset[3], vertex[DIM]; geometry.vertexOffset(a,offset);
    for(int d=0; d<DIM; ++d) vertex[d]=center[d]+offset[d];
    localIndex support[SUPPORT]; real64 w[SUPPORT], dw[SUPPORT][3];
    if(!point<DIM,false>(vertex,map,origin,spacing,cells,support,w,dw)) return false;
    real64 basisGradient[DIM]={};
    if constexpr( SIMPLEX )
    {
      for(int d=0; d<DIM; ++d) basisGradient[d]=reciprocal[a][d];
    }
    else
    {
      // Canonical box vertex order is cyclic in XY, then the two Z planes.
      int const signs[8][3]={{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
                             {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
      for(int d=0; d<DIM; ++d) for(int k=0; k<DIM; ++k)
        basisGradient[d]+=signs[a][k]*reciprocal[k+1][d]/(1<<(DIM-1));
    }
    // Integrate the affine material basis against r dA. On a triangle the
    // barycentric gradients are constant; on a box their radial moments matter.
    real64 cornerWeight=1.0/CORNERS;
    if constexpr(AXISYMMETRIC)
    {
      cornerWeight += offset[0]/(12*center[0]);
      if constexpr(!SIMPLEX)
      {
        int const signs[4][2]={{-1,-1},{1,-1},{1,1},{-1,1}};
        int const s=signs[a][0], t=signs[a][1];
        for(int d=0;d<2;++d)
          basisGradient[d]=s*reciprocal[1][d]*(0.5+t*generators[1][0]/(12*center[0]))
                          +t*reciprocal[2][d]*(0.5+s*generators[0][0]/(12*center[0]));
      }
    }
    real64 sum=0, gradientSum[DIM]={}, hoopSum=0;
    for(int b=0; b<SUPPORT; ++b)
    {
      int const n=a*SUPPORT+b;
      nodes[n]=support[b];
      values[n]=b==SUPPORT-1 ? cornerWeight-sum : w[b]*cornerWeight;
      sum+=values[n];
      if constexpr(DIM==2) gradients[n][2]=0;
      if constexpr(AXISYMMETRIC)
      {
        real64 const hoopWeight=1.0/(CORNERS*center[0]);
        gradients[n][2]=b==SUPPORT-1 ? hoopWeight-hoopSum : w[b]*hoopWeight;
        hoopSum+=gradients[n][2];
      }
      for(int d=0; d<DIM; ++d)
      {
        gradients[n][d]=b==SUPPORT-1 ? basisGradient[d]-gradientSum[d] : w[b]*basisGradient[d];
        gradientSum[d]+=gradients[n][d];
      }
    }
  }
  return true;
}

} // namespace mpmMapping
} // namespace geos
#endif
