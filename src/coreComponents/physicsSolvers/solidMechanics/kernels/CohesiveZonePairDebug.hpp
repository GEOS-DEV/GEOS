/*
 * SPDX-License-Identifier: LGPL-2.1-only
 * Copyright (c) 2019- GEOS/GEOSX Contributors
 */
#ifndef GEOS_PHYSICSSOLVERS_COHESIVEZONEPAIRDEBUG_HPP_
#define GEOS_PHYSICSSOLVERS_COHESIVEZONEPAIRDEBUG_HPP_

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>

namespace geos
{
namespace solidMechanicsMPMKernels
{

/** Flat device-record layout. Vector/matrix offsets reserve 3/9 components. */
struct CohesivePairDebug
{
  enum Column : int
  {
    massA = 0, massB = 1, referenceAreaA = 2, referenceAreaB = 3,
    displacementA = 4, displacementB = 7,
    mappedNormalA = 10, mappedNormalB = 13,
    referenceNormalA = 16, referenceNormalB = 19,
    cofactorA = 22, cofactorB = 31,
    rawJump = 40, jump = 43, rawNormal = 46, normal = 49,
    normalMagnitude = 52, normalDisplacement = 53, tangentialDisplacement = 54,
    tangentialJump = 55, currentAreaA = 58, currentAreaB = 61,
    pairAreaVector = 64, pairArea = 67,
    returnedNormalStress = 68, appliedNormalStress = 69, shearStress = 70,
    tractionA = 71, tractionB = 74, pairForceA = 77, pairForceB = 80,
    slotForceA = 83, slotForceB = 86,
    temperature = 89, damageBefore = 90, damageAfter = 91,
    referencePositionA = 92, referencePositionB = 95
  };

  // LvArray requires integral dimension types; an enum does not satisfy std::is_integral.
  static constexpr int numColumns = referencePositionB + 3;

  enum Status : int
  {
    inactive = 0, degenerateNormal = 1, evaluated = 2,
    inactiveNormalA = 3, inactiveMassB = 4, inactiveNormalB = 5
  };
};

/** Time cadence takes precedence over cycle cadence; bounds include their endpoints. */
inline bool cohesivePairDumpDue( double const time, int const cycle,
                                int const cycleInterval, double const timeInterval,
                                double const startTime, double const endTime,
                                double const nextTime, int const lastCycle )
{
  if( cycle == lastCycle || time < startTime || time > endTime )
  {
    return false;
  }
  return timeInterval > 0.0
           ? time >= std::max( nextTime, startTime )
           : cycleInterval > 0 && cycle % cycleInterval == 0;
}

/** Advance past missed samples without an unbounded loop on adaptive/restarted steps. */
inline double nextCohesivePairDumpTime( double const time, double const interval,
                                       double const startTime )
{
  double const remainder = std::fmod( std::max( time - startTime, 0.0 ), interval );
  double const tolerance = std::min( 0.5 * interval, 4.0 * std::numeric_limits< double >::epsilon() *
                                    std::max( { std::abs( time ), std::abs( startTime ), interval } ) );
  // Decimal sample times can land one ulp below an interval multiple (e.g. 0.3/0.1).
  // Treat that as the sampled boundary rather than scheduling an adjacent-cycle duplicate.
  bool const atBoundary = remainder <= tolerance || interval - remainder <= tolerance;
  double const next = time + ( atBoundary ? interval : interval - remainder );
  return next > time ? next : std::nextafter( time, std::numeric_limits< double >::infinity() );
}

/** Write one captured kernel record. This accepts LvArray slices or standard arrays. */
template< typename ROW >
void writeCohesivePairDebugRecord( std::ostream & out, ROW const & row,
                                  int const status, double const smallMass )
{
  using D = CohesivePairDebug;
  auto scalar = [&]( char const * name, int const offset )
  {
    out << "  " << name << " = " << row[offset] << '\n';
  };
  auto vector = [&]( char const * name, int const offset )
  {
    out << "  " << name << " = [" << row[offset] << ", " << row[offset + 1]
        << ", " << row[offset + 2] << "]\n";
  };
  auto matrix = [&]( char const * name, int const offset )
  {
    out << "  " << name << " =\n";
    for( int i = 0; i < 3; ++i )
    {
      out << "    [" << row[offset + 3 * i] << ", " << row[offset + 3 * i + 1]
          << ", " << row[offset + 3 * i + 2] << "]\n";
    }
  };
  bool const inactive = status != D::evaluated && status != D::degenerateNormal;
  out << "  status = " << ( status == D::evaluated ? "evaluated" :
                              status == D::degenerateNormal ? "skipped_degenerate_normal" :
                              status == D::inactiveNormalA ? "skipped_mapped_normal_A" :
                              status == D::inactiveMassB ? "skipped_mass_B" :
                              status == D::inactiveNormalB ? "skipped_mapped_normal_B" :
                              "skipped_mass_A" ) << '\n';
  if( inactive )
  {
    out << "  first_failed_criterion: mass must exceed " << smallMass
        << "; mapped normal squared magnitude must exceed 1e-16\n";
  }
  scalar( "mass_A", D::massA );
  scalar( "mass_B", D::massB );
  scalar( "reference_area_A", D::referenceAreaA );
  scalar( "reference_area_B", D::referenceAreaB );
  vector( "reference_node_position_A", D::referencePositionA );
  vector( "reference_node_position_B", D::referencePositionB );
  vector( "displacement_A", D::displacementA );
  vector( "displacement_B", D::displacementB );
  vector( "mapped_normal_A", D::mappedNormalA );
  vector( "mapped_normal_B", D::mappedNormalB );
  vector( "reference_normal_A", D::referenceNormalA );
  vector( "reference_normal_B", D::referenceNormalB );
  matrix( "deformation_gradient_cofactor_A", D::cofactorA );
  matrix( "deformation_gradient_cofactor_B", D::cofactorB );
  scalar( "temperature", D::temperature );
  scalar( "damage_before", D::damageBefore );
  scalar( "damage_after", D::damageAfter );
  if( !inactive )
  {
    vector( "current_area_vector_A", D::currentAreaA );
    vector( "current_area_vector_B", D::currentAreaB );
    vector( "mass_weighted_normal_before_normalization", D::rawNormal );
    scalar( "normal_magnitude", D::normalMagnitude );
  }
  if( status == D::evaluated )
  {
    vector( "unit_normal_A_to_B", D::normal );
    vector( "jump_dA_minus_dB_before_periodic_reduction", D::rawJump );
    vector( "jump_after_periodic_reduction", D::jump );
    scalar( "normal_opening", D::normalDisplacement );
    vector( "tangential_jump", D::tangentialJump );
    scalar( "tangential_opening_magnitude", D::tangentialDisplacement );
    vector( "mass_weighted_pair_area_vector", D::pairAreaVector );
    scalar( "projected_pair_area", D::pairArea );
    scalar( "normal_stress_returned_by_model", D::returnedNormalStress );
    scalar( "normal_stress_after_compression_suppression", D::appliedNormalStress );
    scalar( "shear_stress_returned_by_model", D::shearStress );
    vector( "traction_A", D::tractionA );
    vector( "traction_B", D::tractionB );
    vector( "this_pair_force_A", D::pairForceA );
    vector( "this_pair_force_B", D::pairForceB );
  }
  else
  {
    out << "  constitutive_update = not_called; this_pair_force = [0, 0, 0]\n";
  }
  vector( "accumulated_slot_force_A_all_pairs", D::slotForceA );
  vector( "accumulated_slot_force_B_all_pairs", D::slotForceB );
}

} // namespace solidMechanicsMPMKernels
} // namespace geos
#endif
