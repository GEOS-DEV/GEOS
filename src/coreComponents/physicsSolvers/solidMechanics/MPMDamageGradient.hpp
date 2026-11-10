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

/**
 * @brief Reconstructs the normalized damage gradient in the meridional plane
 *        for no-swirl, axisymmetric dual-field (DFG) calculations.
 *
 * Damage is a scalar geometric partition in (r, y), rather than a quantity
 * integrated over the full three-dimensional ring. Each physical particle
 * therefore contributes its meridional area A = V / (2*pi*r), where V is its
 * full-ring volume. The radial coordinate must be strictly positive.
 *
 * The scalar damage field is even across the symmetry axis. For every real
 * sample at (r, y), add() also evaluates the same damage and area at (-r, y).
 * These mirrored positions contribute only to the reconstruction quadrature:
 * they do not introduce particles, material mass, momentum, additional contact
 * fields, or cohesive interfaces. At a query point on the axis, their radial
 * kernel-gradient contributions cancel pairwise.
 *
 * For kernel K and samples indexed by p, the reconstructed field is
 *
 *   d_bar(x) = D(x) / W(x),
 *   D(x)     = sum_p A_p * d_p * K(x - x_p),
 *   W(x)     = sum_p A_p * K(x - x_p).
 *
 * gradient() evaluates the quotient-rule derivative of D/W. The compact
 * cubic kernel is the polynomial used by the existing MPM kernelDevice()
 * path; a common multiplicative normalization of K cancels from D/W.
 *
 * The caller supplies the existing particle-neighbor list and the current
 * surface-as-damage convention. Mechanical P2G quantities continue to use
 * physical ring-volume weights, not the meridional areas used here.
 */
class AxisymmetricDamageGradient
{
public:
  /**
   * @brief Initializes one damage-gradient reconstruction at a query location.
   * @param radius Radial coordinate of the query point.
   * @param axial Axial coordinate of the query point.
   * @param supportRadius Compact kernel support radius; must be positive.
   *
   * The constructor stores the query and support radius. Accumulators start
   * at zero so each instance represents one independent reconstruction.
   */
  GEOS_HOST_DEVICE
  AxisymmetricDamageGradient( real64 const radius, real64 const axial,
                              real64 const supportRadius ):
    m_radius( radius ), m_axial( axial ), m_supportRadius( supportRadius )
  {}

  /**
   * @brief Adds a physical particle's damage and its reflection across r = 0.
   * @param radius Physical particle radius, required to be positive.
   * @param axial Particle axial coordinate.
   * @param volume Full-ring particle volume.
   * @param damage Scalar damage value used by the caller's DFG partition.
   *
   * Dividing by the ring circumference converts the particle's physical
   * volume into a meridional area. Both images use this same area and damage,
   * implementing even reflection of the scalar field about the axis.
   *
   * This routine does not guard against radius = 0; callers must supply
   * positive particle radii and a positive support radius.
   */
  GEOS_HOST_DEVICE
  void add( real64 const radius, real64 const axial, real64 const volume,
            real64 const damage )
  {
    // Undo the circumferential ring weight for the meridional reconstruction.
    real64 const area = volume / ( axisymmetric::twoPi * radius );

    // The reflected image enforces even radial parity without creating state.
    // addImage() independently rejects either image outside kernel support.
    addImage( radius, axial, area, damage );
    addImage( -radius, axial, area, damage );
  }

  /**
   * @brief Returns the meridional gradient of the normalized damage field.
   * @param[out] result Three-component output: radial, axial, and zero hoop.
   *
   * Using the accumulated weighted sums W and D, and their spatial
   * derivatives, apply
   *
   *   grad(D/W) = (grad(D) - (D/W)*grad(W)) / W.
   *
   * An empty support leaves W = 0; in that case the result is identically
   * zero. The hoop component is also zero for the no-swirl scalar field.
   */
  GEOS_HOST_DEVICE
  void gradient( real64 ( &result )[3] ) const
  {
    // Initialize every component, including the inactive hoop direction.
    for( int d = 0; d < 3; ++d ) result[d] = 0.0;
    if( m_weight > 0.0 )
    {
      // Subtract the normalization derivative before dividing by total weight.
      real64 const mean = m_damageWeight / m_weight;
      for( int d = 0; d < 2; ++d )
        result[d] = ( m_damageGradient[d] - mean * m_weightGradient[d] ) / m_weight;
    }
  }

private:
  /**
   * @brief Accumulates the kernel moments for one real or reflected sample.
   * @param radius Radial coordinate of the sample, possibly negative.
   * @param axial Sample axial coordinate.
   * @param area Meridional area associated with the physical particle.
   * @param damage Scalar damage copied unchanged to the mirrored sample.
   *
   * The support-limited cubic kernel and its query-point gradient are
   *
   *   q       = |x_query - x_sample| / h,
   *   K(q)    = (1-q)^2 * (1+2*q),              0 <= q < 1,
   *   grad K  = -6*(1-q)*(x_query-x_sample)/h^2.
   *
   * Multiplication by sample area supplies the moments W, D, grad(W), and
   * grad(D). These four sums are sufficient for gradient() to differentiate
   * the normalized reconstruction without evaluating a separate field value.
   */
  GEOS_HOST_DEVICE
  void addImage( real64 const radius, real64 const axial, real64 const area,
                 real64 const damage )
  {
    // Kernel distances are measured in the meridional plane, never in hoop.
    real64 const delta[2] = { m_radius - radius, m_axial - axial };
    real64 const distance = ::sqrt( delta[0] * delta[0] + delta[1] * delta[1] );
    if( distance >= m_supportRadius ) return;

    // Compact support and the analytic derivative of its cubic polynomial.
    // The derivative is regular at coincident query and sample locations.
    real64 const q = distance / m_supportRadius;
    real64 const weight = area * ( 1.0 - q ) * ( 1.0 - q ) * ( 1.0 + 2.0 * q );
    real64 const scale = -6.0 * area * ( 1.0 - q ) / ( m_supportRadius * m_supportRadius );

    // Accumulate the unnormalized area-weighted damage and support weights.
    m_weight += weight;
    m_damageWeight += damage * weight;

    // grad K is proportional to the query-to-sample displacement; the
    // quotient rule in gradient() accounts for the varying normalization.
    for( int d = 0; d < 2; ++d )
    {
      m_weightGradient[d] += scale * delta[d];
      m_damageGradient[d] += damage * scale * delta[d];
    }
  }

  // Reconstruction location and compact support in meridional coordinates.
  real64 m_radius, m_axial, m_supportRadius;

  // Zeroth kernel moments: W = sum(A*K) and D = sum(A*damage*K).
  real64 m_weight = 0.0, m_damageWeight = 0.0;

  // First kernel moments: gradients of W and D in (r, y).
  real64 m_weightGradient[2] = {}, m_damageGradient[2] = {};
};

} // namespace mpmDamage
} // namespace geos

#endif
