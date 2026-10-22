/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
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
/** Integral of the native Q1 nodal hat over an XY interface, times thickness.
 * Clip the infinite line to the nodal support rectangle, then split at x=0
 * and y=0. Two-point Gauss integration is exact on each quadratic segment.
 * This is physical crack/contact area; the corresponding geometric measure
 * is line length, independent of the domain thickness.
 */
template<bool AXISYMMETRIC, class NORMAL, class POSITION>
GEOS_HOST_DEVICE real64 interfaceLineIntegral( real64 const hx, real64 const hy, real64 const measureScale,
                                              NORMAL const & normal, POSITION const & position )
{
  real64 const norm=::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
  if(!(norm>0 && hx>0 && hy>0)) return 0;
  if constexpr(!AXISYMMETRIC) { if(!(measureScale>0)) return 0; }
  real64 const n[2]={normal[0]/norm,normal[1]/norm};
  real64 const distance=n[0]*position[0]+n[1]*position[1];
  real64 const center[2]={distance*n[0],distance*n[1]};
  real64 const tangent[2]={-n[1],n[0]}, h[2]={hx,hy};
  real64 lo=-DBL_MAX, hi=DBL_MAX;
  for(int d=0; d<2; ++d)
  {
    if(::fabs(tangent[d])<DBL_EPSILON)
    {
      if(::fabs(center[d])>=h[d]) return 0;
      continue;
    }
    real64 a=(-h[d]-center[d])/tangent[d], b=(h[d]-center[d])/tangent[d];
    if(a>b) { real64 const tmp=a; a=b; b=tmp; }
    if(a>lo) lo=a;
    if(b<hi) hi=b;
  }
  // measureScale is the reference node radius in axisymmetry. Clip the
  // meridional line at r=0 before integrating; a node on the axis still has
  // finite support at positive radius. There is no radius floor or division.
  if constexpr(AXISYMMETRIC)
  {
    if(::fabs(tangent[0])<DBL_EPSILON)
    {
      if(!(measureScale+center[0]>0)) return 0;
    }
    else
    {
      real64 const axis=-(measureScale+center[0])/tangent[0];
      if(tangent[0]>0 && axis>lo) lo=axis;
      if(tangent[0]<0 && axis<hi) hi=axis;
    }
  }
  if(!(hi>lo)) return 0;
  real64 breaks[4]={lo,hi,0,0}; int count=2;
  for(int d=0; d<2; ++d) if(::fabs(tangent[d])>=DBL_EPSILON)
  {
    real64 const t=-center[d]/tangent[d];
    if(t>lo && t<hi) breaks[count++]=t;
  }
  for(int a=1; a<count; ++a) for(int b=a; b>0 && breaks[b]<breaks[b-1]; --b)
  { real64 const tmp=breaks[b]; breaks[b]=breaks[b-1]; breaks[b-1]=tmp; }
  real64 integral=0;
  for(int a=1; a<count; ++a)
  {
    real64 const mid=0.5*(breaks[a]+breaks[a-1]), half=0.5*(breaks[a]-breaks[a-1]);
    for(int sign=-1; sign<=1; sign+=2)
    {
      real64 const t=mid+sign*half/::sqrt(3.0);
      real64 weight=1;
      for(int d=0; d<2; ++d) weight*=1-::fabs(center[d]+t*tangent[d])/h[d];
      // Q1 restricted to a line is piecewise quadratic. Multiplication by
      // radius makes it cubic, so the same two-point Gauss rule is exact.
      if constexpr(AXISYMMETRIC) weight*=measureScale+center[0]+t*tangent[0];
      integral+=half*weight;
    }
  }
  return integral*(AXISYMMETRIC ? 6.283185307179586476925286766559 : measureScale);
}

/// Plane-strain physical area: constant thickness times Q1-weighted length.
template<class NORMAL, class POSITION>
GEOS_HOST_DEVICE real64 lineArea( real64 const hx, real64 const hy, real64 const thickness,
                                 NORMAL const & normal, POSITION const & position )
{ return interfaceLineIntegral<false>(hx,hy,thickness,normal,position); }

/// Full-ring interface area integral: integral_Gamma 2*pi*r*N_I ds.
/// position is relative to the reference node, whose absolute radius is rI.
template<class NORMAL, class POSITION>
GEOS_HOST_DEVICE real64 axisymmetricLineArea( real64 const hx, real64 const hy, real64 const rI,
                                             NORMAL const & normal, POSITION const & position )
{ return interfaceLineIntegral<true>(hx,hy,rI,normal,position); }

/// Integral of the Q1 nodal hat over r>=0 with full-ring volume weighting.
/// At the axis this is 2*pi*hx^2*hy/6; away from it, 2*pi*rI*hx*hy.
GEOS_HOST_DEVICE inline real64 axisymmetricNodalVolume( real64 const hx, real64 const hy, real64 const rI )
{
  real64 const normal[2]={0,1}, position[2]={0,0};
  return hy*axisymmetricLineArea(hx,hy,rI,normal,position);
}
} // namespace mpmMeasures
} // namespace geos
#endif
