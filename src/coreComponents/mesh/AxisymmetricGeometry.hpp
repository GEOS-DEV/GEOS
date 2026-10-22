/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MESH_AXISYMMETRICGEOMETRY_HPP
#define GEOS_MESH_AXISYMMETRICGEOMETRY_HPP
#include "mesh/ParticleGeometry.hpp"
namespace geos { namespace axisymmetric {
constexpr real64 twoPi = 6.283185307179586476925286766559;
/// Full-ring volume of an affine meridional domain (radius is its area centroid).
GEOS_HOST_DEVICE inline real64 volume(real64 area, real64 radius)
{ return twoPi * radius * area; }
/// Ring volume of an axis-aligned meridional rectangle, including a solid cylinder.
GEOS_HOST_DEVICE inline real64 boxVolume(real64 rmin, real64 rmax, real64 height)
{ return 0.5 * twoPi * (rmax-rmin) * (rmax+rmin) * height; }
/** Fit ONLY the temporary quadrature domain inside r>=0. The physical domain,
 * mass, reference/current volumes, and stored generators are never modified.
 * Return the existing projection status bits: changed=1, measure reduced=2,
 * invalid=4. There is no arbitrary radius floor or 1/r regularization.
 */
template<class MATRIX>
GEOS_HOST_DEVICE int fitIntegrationDomain(real64 radius, particleGeometry::Shape shape, MATRIX & g)
{
  if(!(radius>0 && radius<=DBL_MAX)) return 4;
  auto domain=particleGeometry::affine::Domain::fromGenerators(shape,g);
  real64 inward=0;
  for(int a=0;a<domain.verticesCount();++a)
  {
    real64 x[3]; domain.vertexOffset(a,x);
    if(-x[0]>inward) inward=-x[0];
  }
  if(inward<=radius) return 0;
  real64 const factor=(radius/inward)*(1-64*DBL_EPSILON);
  for(int a=0;a<2;++a) for(int d=0;d<2;++d) g[a][d]*=factor;
  return 3;
}
/// Convert geometric subdivision fraction to full-ring extensive fraction.
GEOS_HOST_DEVICE inline real64 splitFraction(real64 fraction, real64 radius, real64 radialOffset)
{ return fraction * (1 + radialOffset/radius); }
} }
#endif
