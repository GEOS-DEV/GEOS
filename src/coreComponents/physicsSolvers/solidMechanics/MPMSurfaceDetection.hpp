/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS.
 */

/**
 * @file MPMSurfaceDetection.hpp
 * @brief Axisymmetric, meridional-neighborhood centroid reconstruction for MPM surface detection.
 */

#ifndef GEOS_SOLIDMECHANICS_MPMSURFACEDETECTION_HPP
#define GEOS_SOLIDMECHANICS_MPMSURFACEDETECTION_HPP

#include "mesh/AxisymmetricGeometry.hpp"
#include <cmath>

namespace geos
{
namespace mpmSurface
{

/**
 * @brief Accumulates the geometric displacement of an axisymmetric particle
 *        neighborhood's kernel-weighted centroid from a query point.
 *
 * Surface detection uses the spatial distribution of neighboring particles,
 * not their share of physical ring mass. A particle at radius r represents a
 * meridional area A = V / (2*pi*r), where V is its current full-ring volume.
 * Using V directly would give otherwise equivalent samples greater geometric
 * weight at larger radii and bias the centroid away from the symmetry axis.
 *
 * Reflect every physical sample (r, y) into (-r, y), carrying the same area.
 * This supplies the missing half of an even, axisymmetric material distribution
 * near r = 0. At an on-axis query point, equal reflected contributions cancel
 * in the radial centroid displacement; the axis is therefore not mistaken for
 * a free surface solely because no particles occupy negative physical radii.
 * The images are quadrature contributions only: no particles, additional
 * contact fields, momentum, mass, or cohesive interfaces are created.
 *
 * For a query x_q = (r_q, y_q), and real and reflected samples x_p, accumulate
 *
 *   W = sum_p A_p K(|x_p-x_q|/h),
 *   M = sum_p A_p K(|x_p-x_q|/h) (x_p-x_q).
 *
 * displacement() returns M/W, the vector from the query point toward the
 * reconstructed neighborhood centroid. The compact kernel is
 * K(q) = (1-q)^2 (1+2q) for 0 <= q < 1, and zero otherwise. Any common
 * multiplicative kernel normalization cancels in the ratio.
 *
 * The caller supplies the physical neighbor list and interprets the returned
 * displacement when applying its surface-detection criterion. A reflected
 * point cannot be closer to a nonnegative-radius query than its physical
 * counterpart, so an additional neighbor search is unnecessary.
 */
class AxisymmetricCentroid
{
public:
  /**
   * @brief Starts one independent centroid reconstruction at a query point.
   * @param radius Query-point radial coordinate, normally nonnegative.
   * @param axial Query-point axial coordinate.
   * @param supportRadius Compact-kernel support radius, required to be positive.
   *
   * The stored query coordinates are used by both physical and reflected
   * contributions. The total kernel weight and weighted offsets start at zero.
   */
  GEOS_HOST_DEVICE
  AxisymmetricCentroid( real64 const radius, real64 const axial, real64 const supportRadius ):
    m_radius( radius ), m_axial( axial ), m_supportRadius( supportRadius )
  {}

  /**
   * @brief Adds a physical particle and its even reflection across the axis.
   * @param radius Physical particle radius; must be strictly positive.
   * @param axial Particle axial coordinate.
   * @param volume Current physical volume of the complete axisymmetric ring.
   *
   * Dividing the ring volume by its circumference produces the meridional
   * area used for geometric surface detection. Both samples use that same
   * area; addImage() excludes either one when outside the kernel support.
   * This routine does not check for zero radius, so the caller must enforce
   * the stated precondition before division.
   */
  GEOS_HOST_DEVICE
  void add( real64 const radius, real64 const axial, real64 const volume )
  {
    // Convert physical ring volume into the particle's meridional area.
    // This weighting is local to the geometric reconstruction, not P2G.
    real64 const area = volume / ( axisymmetric::twoPi * radius );

    // The reflection represents scalar-even geometric sampling across r = 0.
    // Each image independently passes through the compact-support cutoff.
    addImage( radius, axial, area );
    addImage( -radius, axial, area );
  }

  /**
   * @brief Returns the kernel-weighted neighborhood-centroid displacement.
   * @param[out] result Displacement in (radial, axial, zero hoop) components.
   * @return True if positive kernel weight was accumulated; otherwise false.
   *
   * The first two components equal M/W, where M is the sum of weighted
   * sample-to-query offsets and W is the sum of sample weights. Working with
   * offsets avoids cancellation from subtracting large absolute coordinates.
   * If no supported sample contributed, the method returns false and leaves
   * all components zero; empty support alone is not evidence of an exposed
   * material surface. The third component is inactive in the no-swirl plane.
   */
  GEOS_HOST_DEVICE
  bool displacement( real64 ( &result )[3] ) const
  {
    // Always initialize the complete three-component output, including hoop.
    for( int d = 0; d < 3; ++d ) result[d] = 0.0;
    if( !( m_weight > 0.0 ) ) return false;

    // Normalize the first moment by the geometric kernel weight.
    for( int d = 0; d < 2; ++d ) result[d] = m_offset[d] / m_weight;
    return true;
  }

private:
  /**
   * @brief Accumulates one physical or reflected sample inside kernel support.
   * @param radius Sample radial coordinate; may be negative for an image.
   * @param axial Sample axial coordinate.
   * @param area Meridional area of the associated physical particle.
   *
   * Form the sample-to-query offset delta and its meridional distance. For
   * distance >= h the sample contributes nothing. Otherwise, with q = d/h,
   * accumulate A*K(q) in W and A*K(q)*delta in the first moment M.
   * No spatial derivatives, ring-radius factors, or third-direction moments
   * are needed for this centroid-based geometric test.
   */
  GEOS_HOST_DEVICE
  void addImage( real64 const radius, real64 const axial, real64 const area )
  {
    // Use only the meridional coordinates; the image may have negative r.
    real64 const delta[2] = { radius - m_radius, axial - m_axial };
    real64 const distance = ::sqrt( delta[0] * delta[0] + delta[1] * delta[1] );
    if( distance >= m_supportRadius ) return;

    // Compact cubic weight. Its constant normalization is irrelevant to M/W.
    real64 const q = distance / m_supportRadius;
    real64 const weight = area * ( 1.0 - q ) * ( 1.0 - q ) * ( 1.0 + 2.0 * q );

    // Maintain the zeroth (W) and first (M) weighted geometric moments.
    m_weight += weight;
    for( int d = 0; d < 2; ++d ) m_offset[d] += weight * delta[d];
  }

  // Query location in meridional coordinates and compact support radius.
  real64 m_radius, m_axial, m_supportRadius;

  // Zeroth geometric kernel moment and first offset moment in (r, y).
  real64 m_weight = 0.0, m_offset[2] = {};
};

} // namespace mpmSurface
} // namespace geos

#endif
