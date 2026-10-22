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
template<class NORMAL, class POSITION>
GEOS_HOST_DEVICE real64 lineArea( real64 const hx, real64 const hy, real64 const thickness,
                                 NORMAL const & normal, POSITION const & position )
{
  real64 const norm=::sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
  if(!(norm>0 && hx>0 && hy>0 && thickness>0)) return 0;
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
      integral+=half*weight;
    }
  }
  return integral*thickness;
}
} // namespace mpmMeasures
} // namespace geos
#endif
