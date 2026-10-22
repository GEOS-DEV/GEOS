/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */
#ifndef GEOS_SOLIDMECHANICS_MPMDAMAGEGRADIENT_HPP
#define GEOS_SOLIDMECHANICS_MPMDAMAGEGRADIENT_HPP

#include "common/DataTypes.hpp"
#include "mesh/AxisymmetricGeometry.hpp"
#include <cmath>

namespace geos
{
namespace mpmDamage
{
/** Meridional normalized-kernel damage gradient for no-swirl DFG.
 *
 * Damage describes a geometric partition in (r,y). Weight its reconstruction
 * with meridional area V/(2*pi*r), not full-ring volume. Reflect each sample
 * across r=0: a scalar axisymmetric damage field is even in r, so its radial
 * derivative vanishes at the axis. These images are quadrature samples only;
 * they create no particles, mass, momentum, contact fields or CZ interfaces.
 *
 * The compact cubic kernel is the existing MPM kernelDevice() polynomial.
 * Its common normalization factor cancels in grad(sum(A*d*K)/sum(A*K)).
 * The caller retains the existing neighbor set and surface-as-damage
 * convention. Physical P2G quantities retain ring weights.
 */
class AxisymmetricDamageGradient
{
public:
  GEOS_HOST_DEVICE
  AxisymmetricDamageGradient( real64 const radius, real64 const axial,
                              real64 const supportRadius ):
    m_radius( radius ), m_axial( axial ), m_supportRadius( supportRadius )
  {}

  /// Particle radii and supportRadius must be positive; volume is full-ring.
  GEOS_HOST_DEVICE
  void add( real64 const radius, real64 const axial, real64 const volume,
            real64 const damage )
  {
    real64 const area = volume / ( axisymmetric::twoPi * radius );
    addImage( radius, axial, area, damage );
    addImage( -radius, axial, area, damage );
  }

  GEOS_HOST_DEVICE
  void gradient( real64 (&result)[3] ) const
  {
    for( int d = 0; d < 3; ++d ) result[d] = 0.0;
    if( m_weight > 0.0 )
    {
      real64 const mean = m_damageWeight / m_weight;
      for( int d = 0; d < 2; ++d )
        result[d] = ( m_damageGradient[d] - mean * m_weightGradient[d] ) / m_weight;
    }
  }

private:
  GEOS_HOST_DEVICE
  void addImage( real64 const radius, real64 const axial, real64 const area,
                 real64 const damage )
  {
    real64 const delta[2] = { m_radius - radius, m_axial - axial };
    real64 const distance = ::sqrt( delta[0] * delta[0] + delta[1] * delta[1] );
    if( distance >= m_supportRadius ) return;
    real64 const q = distance / m_supportRadius;
    real64 const weight = area * ( 1.0 - q ) * ( 1.0 - q ) * ( 1.0 + 2.0 * q );
    real64 const scale = -6.0 * area * ( 1.0 - q ) / ( m_supportRadius * m_supportRadius );
    m_weight += weight;
    m_damageWeight += damage * weight;
    for( int d = 0; d < 2; ++d )
    {
      m_weightGradient[d] += scale * delta[d];
      m_damageGradient[d] += damage * scale * delta[d];
    }
  }

  real64 m_radius, m_axial, m_supportRadius;
  real64 m_weight = 0.0, m_damageWeight = 0.0;
  real64 m_weightGradient[2] = {}, m_damageGradient[2] = {};
};
}
}
#endif
