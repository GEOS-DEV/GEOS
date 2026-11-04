/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
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
constexpr real64 twoPi = 6.283185307179586476925286766559;

/** Integrate the geometric Q1 basis over a quad/triangle, splitting at grid
 * cells. This is deliberately independent of the mechanical CPDI/CPTI
 * interpolant: a conforming triangle tiling must fill exactly the same Q1
 * nodal support as quads. emit(i,j,S,B_r,B_y) receives normalized ring volume
 * and boundary-normal weights B=(<dN/dr>+<N/r>,<dN/dy>).
 * Returns false before emitting anything if geometry leaves the local grid.
 */
template<class POINT,class MATRIX,class ORIGIN,class SPACING,class CELLS,class EMIT>
GEOS_HOST_DEVICE bool projectGeometry( bool const simplex, POINT const & center, MATRIX const & geometry,
                                       ORIGIN const & origin, SPACING const & h, CELLS const & cells,
                                       EMIT const & emit )
{
  using namespace particleGeometry;
  auto const domain=affine::Domain::fromGenerators(simplex?Shape::Triangle:Shape::Parallelogram,geometry);
  real64 const denominator=center[0]*domain.measure(); // 2*pi cancels in normalized weights
  if( !(denominator>0) ) return false;
  real64 vertices[8][3]; domain.vertices(center,vertices);
  int const count=domain.verticesCount();
  int first[2],last[2];
  for(int d=0;d<2;++d)
  {
    real64 lo=DBL_MAX,hi=-DBL_MAX;
    for(int v=0;v<count;++v){lo=::fmin(lo,vertices[v][d]);hi=::fmax(hi,vertices[v][d]);}
    real64 const a=(lo-origin[d])/h[d],b=(hi-origin[d])/h[d];
    real64 const tolerance=64*DBL_EPSILON*(cells[d]+1);
    if(a < -tolerance || b > cells[d]+tolerance) return false;
    first[d]=static_cast<int>(::floor(::fmax(0.0,a)));
    last[d]=static_cast<int>(::floor(::fmin(static_cast<real64>(cells[d]),b)));
    if(last[d]>=cells[d])last[d]=cells[d]-1;
  }
  constexpr real64 q[3]={0.11270166537925831148,0.5,0.88729833462074168852};
  constexpr real64 w[3]={5.0/18,4.0/9,5.0/18};
  for(int i=first[0];i<=last[0];++i)for(int j=first[1];j<=last[1];++j)
  {
    real64 polygon[10][2]={};int size=count;
    for(int v=0;v<count;++v)for(int d=0;d<2;++d)polygon[v][d]=vertices[v][d];
    real64 const cell[2]={origin[0]+i*h[0],origin[1]+j*h[1]};
    for(int face=0;face<4 && size>0;++face)
    {
      int const d=face/2;real64 const sign=face%2==0?-1:1;
      real64 const boundary=cell[d]+(face%2)*h[d];
      real64 clipped[10][2]={};int nextSize=0;
      for(int v=0;v<size;++v)
      {
        int const next=(v+1)%size;
        real64 const a=sign*(polygon[v][d]-boundary),b=sign*(polygon[next][d]-boundary);
        if(a<=0){clipped[nextSize][0]=polygon[v][0];clipped[nextSize++][1]=polygon[v][1];}
        if((a<0 && b>0)||(a>0 && b<0))
        {real64 const t=a/(a-b);for(int k=0;k<2;++k)clipped[nextSize][k]=polygon[v][k]+t*(polygon[next][k]-polygon[v][k]);++nextSize;}
      }
      size=nextSize;for(int v=0;v<size;++v)for(int k=0;k<2;++k)polygon[v][k]=clipped[v][k];
    }
    real64 S[4]={},B[4][2]={};
    for(int v=1;v+1<size;++v)
    {
      real64 const ax=polygon[v][0]-polygon[0][0],ay=polygon[v][1]-polygon[0][1];
      real64 const bx=polygon[v+1][0]-polygon[0][0],by=polygon[v+1][1]-polygon[0][1];
      real64 const jac=::fabs(ax*by-ay*bx)/denominator;
      for(int u=0;u<3;++u)for(int z=0;z<3;++z)
      {
        real64 const t=q[u],s=(1-t)*q[z];
        real64 const r=polygon[0][0]+t*ax+s*bx,y=polygon[0][1]+t*ay+s*by;
        real64 const xhat=(r-cell[0])/h[0],yhat=(y-cell[1])/h[1],weight=jac*w[u]*w[z]*(1-t);
        for(int a=0;a<2;++a)for(int b=0;b<2;++b)
        {
          int const node=2*a+b;real64 const nx=a?xhat:1-xhat,ny=b?yhat:1-yhat;
          S[node]+=weight*r*nx*ny;
          B[node][0]+=weight*(r*(a?1:-1)/h[0]*ny+nx*ny);
          B[node][1]+=weight*r*nx*(b?1:-1)/h[1];
        }
      }
    }
    for(int a=0;a<2;++a)for(int b=0;b<2;++b)
      if(S[2*a+b]>0)emit(i+a,j+b,S[2*a+b],B[2*a+b][0],B[2*a+b][1]);
  }
  return true;
}

/// Gather a nodal plane as a particle-relative normal offset. The cylindrical
/// mapping centroid differs from the meridional area centroid, so averaging
/// node-relative offsets directly biases even a perfectly planar surface.
template<class NODE,class POINT,class NORMAL,class POSITION,class OUTPUT>
GEOS_HOST_DEVICE void gatherOffset(real64 const weight,NODE const & node,POINT const & particle,
                                   NORMAL const & normal,POSITION const & offset,OUTPUT && result)
{
  real64 const norm2=normal[0]*normal[0]+normal[1]*normal[1];
  if(!(norm2>0))return;
  real64 const distance=((node[0]-particle[0]+offset[0])*normal[0]+
                         (node[1]-particle[1]+offset[1])*normal[1])/norm2;
  for(int d=0;d<2;++d)result[d]+=weight*distance*normal[d];
}

/** Geometric Q1 support, clipped at the axis and physical symmetry planes.
 * Offsets are relative to the node. The inactive geometric coordinate has
 * no role; every volume/area includes the complete circumference once.
 */
struct Support
{
  real64 h[2], radius, lo[2], hi[2];

  GEOS_HOST_DEVICE real64 weight( real64 const x, real64 const y ) const
  { return (1-::fabs(x)/h[0])*(1-::fabs(y)/h[1])*twoPi*(radius+x); }

  GEOS_HOST_DEVICE real64 volume() const
  {
    if( !(hi[0]>lo[0] && hi[1]>lo[1]) ) return 0;
    real64 radial=0, axial=0;
    for( int sign=-1; sign<=1; sign+=2 )
    {
      real64 const a=sign<0 ? lo[0] : ::fmax(lo[0],0.0);
      real64 const b=sign<0 ? ::fmin(hi[0],0.0) : hi[0];
      if( b>a ) radial += radius*(b-a)+(1-sign*radius/h[0])*(b*b-a*a)/2
                           -sign*(b*b*b-a*a*a)/(3*h[0]);
      real64 const c=sign<0 ? lo[1] : ::fmax(lo[1],0.0);
      real64 const d=sign<0 ? ::fmin(hi[1],0.0) : hi[1];
      if( d>c ) axial += d-c-sign*(d*d-c*c)/(2*h[1]);
    }
    return ::fmax(0.0,twoPi*radial*axial);
  }
};

template< class LOWER, class UPPER, class SYMMETRY >
GEOS_HOST_DEVICE Support support( real64 const hx, real64 const hy,
                                  real64 const radius, real64 const y,
                                  LOWER const & lower, UPPER const & upper,
                                  SYMMETRY const & symmetry )
{
  Support s{{hx,hy},radius,{-hx,-hy},{hx,hy}};
  s.lo[0]=::fmax(s.lo[0],-radius);
  real64 const node[2]={radius,y};
  for( int d=0; d<2; ++d )
  {
    if( symmetry[2*d] ) s.lo[d]=::fmax(s.lo[d],lower[d]-node[d]);
    if( symmetry[2*d+1] ) s.hi[d]=::fmin(s.hi[d],upper[d]-node[d]);
  }
  return s;
}

/** Integral of 2*pi*r*N over a clipped material half-space n.x <= offset.
 * Split the support at both hat-function knots, clip each rectangle by the
 * plane, then integrate its cubic polynomial exactly on a triangle fan.
 * A positive 3x3 Duffy-Gauss rule avoids subtracting nearly equal volumes.
 */
template< class NORMAL >
GEOS_HOST_DEVICE real64 halfSpaceVolume( Support const & s, NORMAL const & normal, real64 const offset )
{
  if( !(s.hi[0]>s.lo[0] && s.hi[1]>s.lo[1]) ) return 0;
  constexpr real64 q[3]={0.11270166537925831148,0.5,0.88729833462074168852};
  constexpr real64 w[3]={5.0/18,4.0/9,5.0/18};
  real64 volume=0;
  for( int i=0; i<2; ++i ) for( int j=0; j<2; ++j )
  {
    real64 const x0=i==0 ? s.lo[0] : ::fmax(0.0,s.lo[0]);
    real64 const x1=i==0 ? ::fmin(0.0,s.hi[0]) : s.hi[0];
    real64 const y0=j==0 ? s.lo[1] : ::fmax(0.0,s.lo[1]);
    real64 const y1=j==0 ? ::fmin(0.0,s.hi[1]) : s.hi[1];
    if( !(x1>x0 && y1>y0) ) continue;
    real64 const corners[4][2]={{x0,y0},{x1,y0},{x1,y1},{x0,y1}};
    real64 polygon[6][2]={}; int count=0;
    for( int k=0; k<4; ++k )
    {
      int const next=(k+1)%4;
      real64 const a=normal[0]*corners[k][0]+normal[1]*corners[k][1]-offset;
      real64 const b=normal[0]*corners[next][0]+normal[1]*corners[next][1]-offset;
      if( a<=0 ) { polygon[count][0]=corners[k][0]; polygon[count++][1]=corners[k][1]; }
      if( (a<0 && b>0) || (a>0 && b<0) )
      {
        real64 const t=a/(a-b);
        polygon[count][0]=corners[k][0]+t*(corners[next][0]-corners[k][0]);
        polygon[count++][1]=corners[k][1]+t*(corners[next][1]-corners[k][1]);
      }
    }
    for( int k=1; k+1<count; ++k )
    {
      real64 const ax=polygon[k][0]-polygon[0][0], ay=polygon[k][1]-polygon[0][1];
      real64 const bx=polygon[k+1][0]-polygon[0][0], by=polygon[k+1][1]-polygon[0][1];
      real64 const jac=::fabs(ax*by-ay*bx);
      for( int u=0; u<3; ++u ) for( int v=0; v<3; ++v )
      {
        real64 const t=q[u], z=(1-t)*q[v];
        volume+=jac*w[u]*w[v]*(1-t)*s.weight(polygon[0][0]+t*ax+z*bx,
                                                          polygon[0][1]+t*ay+z*by);
      }
    }
  }
  return ::fmax(0.0,volume);
}

template< class NORMAL >
GEOS_HOST_DEVICE real64 offsetFromVolume( Support const & s, NORMAL const & normal, real64 const volume )
{
  real64 bound=::fabs(normal[0])*s.h[0]+::fabs(normal[1])*s.h[1];
  real64 a=-bound, b=bound;
  for( int iter=0; iter<48; ++iter )
  {
    real64 const mid=(a+b)/2;
    if( halfSpaceVolume(s,normal,mid)<volume ) a=mid; else b=mid;
  }
  return (a+b)/2;
}

/// Exact Q1-weighted full-ring area on the same clipped support.
template< class NORMAL, class POSITION >
GEOS_HOST_DEVICE real64 lineArea( Support const & s, NORMAL const & normal, POSITION const & position )
{
  real64 const norm=::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
  if( !(norm>0 && s.hi[0]>s.lo[0] && s.hi[1]>s.lo[1]) ) return 0;
  real64 const n[2]={normal[0]/norm,normal[1]/norm};
  real64 const distance=n[0]*position[0]+n[1]*position[1];
  real64 const center[2]={distance*n[0],distance*n[1]}, tangent[2]={-n[1],n[0]};
  real64 lo=-DBL_MAX,hi=DBL_MAX;
  for( int d=0; d<2; ++d )
  {
    if( ::fabs(tangent[d])<DBL_EPSILON )
    { if( center[d]<s.lo[d] || center[d]>s.hi[d] ) return 0; }
    else
    {
      real64 a=(s.lo[d]-center[d])/tangent[d], b=(s.hi[d]-center[d])/tangent[d];
      if( a>b ) { real64 const tmp=a; a=b; b=tmp; }
      lo=::fmax(lo,a); hi=::fmin(hi,b);
    }
  }
  if( !(hi>lo) ) return 0;
  real64 cuts[4]={lo,hi,0,0}; int count=2;
  for( int d=0; d<2; ++d ) if( ::fabs(tangent[d])>=DBL_EPSILON )
  {
    real64 const t=-center[d]/tangent[d];
    if( t>lo && t<hi ) cuts[count++]=t;
  }
  for( int i=1; i<count; ++i ) for( int j=i; j>0 && cuts[j]<cuts[j-1]; --j )
  { real64 const tmp=cuts[j]; cuts[j]=cuts[j-1]; cuts[j-1]=tmp; }
  real64 area=0;
  for( int i=1; i<count; ++i )
  {
    real64 const mid=(cuts[i]+cuts[i-1])/2,half=(cuts[i]-cuts[i-1])/2;
    for( int sign=-1; sign<=1; sign+=2 )
    {
      real64 const t=mid+sign*half/::sqrt(3.0);
      area+=half*s.weight(center[0]+t*tangent[0],center[1]+t*tangent[1]);
    }
  }
  return ::fmax(0.0,area);
}

/** Dimensionless, area-weighted logistic normal fit in the meridional plane.
 * visit(callback) supplies (x/hx,y/hy,label=-1/+1,area/(hx*hy)), including
 * reflected axis samples. A failed/one-sided fit leaves the caller's normal.
 * Normalized coordinates and pivoting avoid dimensional determinant cutoffs.
 */
template< class VISIT >
GEOS_HOST_DEVICE bool logisticNormal( VISIT const & visit, real64 const hx, real64 const hy,
                                      int const maxIter, real64 const tolerance, real64 (&normal)[3] )
{
  real64 const length=::sqrt(hx*hy), penalty=1.e-7;
  real64 phi[3]={normal[0]*hx/length,normal[1]*hy/length,0};
  real64 previous[2]={normal[0],normal[1]};
  real64 classWeight[2]={};
  visit([&](real64,real64,real64 label,real64 weight){classWeight[label>0?1:0]+=weight;});
  if( !(classWeight[0]>0 && classWeight[1]>0) ) return false;
  auto objective=[&](real64 const (&p)[3])
  {
    real64 value=penalty*(p[0]*p[0]+p[1]*p[1]);
    visit([&](real64 x,real64 y,real64 label,real64 weight)
    { real64 const residual=label-::tanh((p[0]*x+p[1]*y+p[2])/2); value+=weight*residual*residual; });
    return value;
  };
  for( int iter=0; iter<maxIter; ++iter )
  {
    real64 system[3][4]={{penalty,0,0,-penalty*phi[0]},
                       {0,penalty,0,-penalty*phi[1]},{0,0,0,0}};
    visit([&](real64 x,real64 y,real64 label,real64 weight)
    {
      real64 const v[3]={x,y,1},f=::tanh((phi[0]*x+phi[1]*y+phi[2])/2),deriv=(1-f*f)/2;
      for( int i=0; i<3; ++i )
      {
        system[i][3]+=weight*deriv*(label-f)*v[i];
        for( int j=0; j<3; ++j ) system[i][j]+=weight*deriv*deriv*v[i]*v[j];
      }
    });
    real64 matrixScale=0;
    for( int i=0; i<3; ++i ) matrixScale=::fmax(matrixScale,::fabs(system[i][i]));
    for( int k=0; k<3; ++k )
    {
      int pivot=k;
      for( int i=k+1; i<3; ++i ) if( ::fabs(system[i][k])>::fabs(system[pivot][k]) ) pivot=i;
      if( ::fabs(system[pivot][k])<=64*DBL_EPSILON*matrixScale ) return false;
      for( int j=k; j<4; ++j ) { real64 const tmp=system[k][j]; system[k][j]=system[pivot][j]; system[pivot][j]=tmp; }
      real64 const diagonal=system[k][k];
      for( int j=k; j<4; ++j ) system[k][j]/=diagonal;
      for( int i=0; i<3; ++i ) if( i!=k )
      { real64 const factor=system[i][k]; for( int j=k; j<4; ++j ) system[i][j]-=factor*system[k][j]; }
    }
    real64 const oldObjective=objective(phi);
    real64 step=1,trial[3]={};
    for( int backtrack=0; backtrack<16; ++backtrack )
    {
      for( int i=0; i<3; ++i ) trial[i]=phi[i]+step*system[i][3];
      if( objective(trial)<=oldObjective ) break;
      step/=2;
    }
    if( objective(trial)>oldObjective ) break;
    for( int i=0; i<3; ++i ) phi[i]=trial[i];
    real64 const n[2]={phi[0]/hx,phi[1]/hy},norm=::sqrt(n[0]*n[0]+n[1]*n[1]);
    if( !(norm>0 && norm<DBL_MAX) ) return false;
    real64 const change=::hypot(n[0]/norm-previous[0],n[1]/norm-previous[1]);
    previous[0]=n[0]/norm; previous[1]=n[1]/norm;
    if( iter>2 && change<tolerance && step*::fabs(system[2][3])<tolerance ) break;
  }
  real64 const norm=::hypot(phi[0]/hx,phi[1]/hy);
  if( !(norm>0 && norm<DBL_MAX) ) return false;
  normal[0]=phi[0]/(hx*norm); normal[1]=phi[1]/(hy*norm); normal[2]=0;
  return true;
}

} // namespace axisymmetricSurface
} // namespace geos
#endif
