/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */
#ifndef GEOS_SOLIDMECHANICS_MPMSURFACEDETECTION_HPP
#define GEOS_SOLIDMECHANICS_MPMSURFACEDETECTION_HPP

#include "mesh/AxisymmetricGeometry.hpp"
#include <cmath>

namespace geos
{
namespace mpmSurface
{
/** Geometric centroid of the meridional neighborhood, relative to its query.
 *
 * Ring volume would bias this geometric test toward larger radii. Use area
 * V/(2*pi*r) and include an image at (-r,y), when within the kernel support,
 * so the axis itself is not classified as a material free surface. The image
 * needs no additional neighbor search: its distance to a positive-radius
 * query is never smaller than the distance of its real sample.
 *
 * The compact cubic kernel is the same polynomial as kernelDevice(). Its
 * common normalization cancels. Accumulating offsets avoids subtracting two
 * large absolute centroids. No images enter physical mass/force transfers.
 */
class AxisymmetricCentroid
{
public:
  GEOS_HOST_DEVICE
  AxisymmetricCentroid( real64 const radius, real64 const axial, real64 const supportRadius ):
    m_radius( radius ), m_axial( axial ), m_supportRadius( supportRadius )
  {}

  /// Input radii/supportRadius are positive; volume is the current full ring.
  GEOS_HOST_DEVICE
  void add( real64 const radius, real64 const axial, real64 const volume )
  {
    real64 const area = volume / ( axisymmetric::twoPi * radius );
    addImage( radius, axial, area );
    addImage( -radius, axial, area );
  }

  /// No supported samples is not evidence of an exposed surface.
  GEOS_HOST_DEVICE
  bool displacement( real64 (&result)[3] ) const
  {
    for( int d = 0; d < 3; ++d ) result[d] = 0.0;
    if( !( m_weight > 0.0 ) ) return false;
    for( int d = 0; d < 2; ++d ) result[d] = m_offset[d] / m_weight;
    return true;
  }

private:
  GEOS_HOST_DEVICE
  void addImage( real64 const radius, real64 const axial, real64 const area )
  {
    real64 const delta[2] = { radius - m_radius, axial - m_axial };
    real64 const distance = ::sqrt( delta[0] * delta[0] + delta[1] * delta[1] );
    if( distance >= m_supportRadius ) return;
    real64 const q = distance / m_supportRadius;
    real64 const weight = area * ( 1.0 - q ) * ( 1.0 - q ) * ( 1.0 + 2.0 * q );
    m_weight += weight;
    for( int d = 0; d < 2; ++d ) m_offset[d] += weight * delta[d];
  }

  real64 m_radius, m_axial, m_supportRadius;
  real64 m_weight = 0.0, m_offset[2] = {};
};
}
}
#endif
