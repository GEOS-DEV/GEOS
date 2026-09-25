/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProjectedGaussSeidelContact.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_PROJECTEDGAUSSSEIDELCONTACT_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_PROJECTEDGAUSSSEIDELCONTACT_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <vector>

namespace geos
{
namespace mpm
{

/** @brief Available nodal material-contact algorithms. */
enum struct ContactSolverOption : integer
{
  Pairwise,
  ProjectedGaussSeidel,
  NewtonRaphson
};

ENUM_STRINGS( ContactSolverOption,
              "Pairwise",
              "ProjectedGaussSeidel",
              "NewtonRaphson" );

/**
 * @brief One node-local constraint used by a coupled contact solve.
 *
 * The normal points out of field A toward field B. For a unilateral
 * constraint, the admissible relative normal velocity is
 *
 *   (vB-vA).n >= targetNormalVelocity.
 *
 * accumulatedNormalImpulse is the non-negative magnitude of the normal
 * impulse; the corresponding impulse on A is -lambda*n. The tangential
 * impulse is the impulse on A and is projected onto the Coulomb disk after
 * every PGS update.
 */
struct ProjectedGaussSeidelContactConstraint
{
  localIndex fieldA = -1;
  localIndex fieldB = -1;
  std::array< real64, 3 > normal = {{ 0.0, 0.0, 0.0 }};
  real64 targetNormalVelocity = 0.0;
  real64 frictionCoefficient = 0.0;
  real64 normalBiasImpulse = 0.0;
  real64 gap = 0.0;
  real64 gapActivationTolerance = 0.0;
  bool hasGap = false;
  bool bilateral = false;
  real64 accumulatedNormalImpulse = 0.0;
  std::array< real64, 3 > accumulatedTangentialImpulse = {{ 0.0, 0.0, 0.0 }};
};

/** @brief Convergence information returned by the node-local PGS solve. */
struct ProjectedGaussSeidelContactResult
{
  integer iterations = 0;
  integer numericalGuardActivations = 0;
  real64 residual = 0.0;
  bool converged = true;
  bool numericalFailure = false;
  bool usedSafeFallback = false;
};

namespace projectedGaussSeidelContact
{

/**
 * @brief Temporarily mask floating-point traps while a guarded node solve runs.
 *
 * Every arithmetic result is still checked before it is accepted.  Masking the
 * traps is the final containment layer: malformed nodal data can therefore be
 * rejected and rolled back instead of terminating the full MPM calculation
 * with SIGFPE.  The caller's floating-point environment is restored on exit.
 */
class FloatingPointEnvironmentGuard
{
public:
  FloatingPointEnvironmentGuard()
    : m_active( std::feholdexcept( &m_environment ) == 0 )
  {}

  ~FloatingPointEnvironmentGuard()
  {
    if( m_active )
    {
      std::fesetenv( &m_environment );
    }
  }

  FloatingPointEnvironmentGuard( FloatingPointEnvironmentGuard const & ) = delete;
  FloatingPointEnvironmentGuard & operator=( FloatingPointEnvironmentGuard const & ) = delete;

private:
  std::fenv_t m_environment;
  bool const m_active;
};

inline bool finiteValue( real64 const value )
{
  return std::isfinite( value );
}

inline bool finiteVector( std::array< real64, 3 > const & value )
{
  return finiteValue( value[0] ) &&
         finiteValue( value[1] ) &&
         finiteValue( value[2] );
}

/** @brief Add two finite values without executing an overflowing addition. */
inline bool safeAdd( real64 const a,
                     real64 const b,
                     real64 & result )
{
  real64 const maximum = std::numeric_limits< real64 >::max();
  if( !finiteValue( a ) || !finiteValue( b ) ||
      ( b > 0.0 && a > maximum - b ) ||
      ( b < 0.0 && a < -maximum - b ) )
  {
    return false;
  }
  result = a + b;
  return finiteValue( result );
}

/** @brief Subtract two finite values without executing an overflowing subtraction. */
inline bool safeSubtract( real64 const a,
                          real64 const b,
                          real64 & result )
{
  if( !finiteValue( b ) )
  {
    return false;
  }
  return safeAdd( a, -b, result );
}

/** @brief Multiply two finite values without executing an overflowing product. */
inline bool safeMultiply( real64 const a,
                          real64 const b,
                          real64 & result )
{
  if( !finiteValue( a ) || !finiteValue( b ) )
  {
    return false;
  }
  if( isZero(a) || isZero(b) )
  {
    result = 0.0;
    return true;
  }

  real64 const absA = std::abs( a );
  real64 const absB = std::abs( b );
  // Only divide the maximum by a value greater than one.  Dividing the
  // maximum by a fractional operand can itself overflow even when a*b is
  // perfectly representable (for example DBL_MAX * 0.5).
  if( absA > 1.0 &&
      absB > std::numeric_limits< real64 >::max() / absA )
  {
    return false;
  }
  result = a * b;
  return finiteValue( result );
}

/** @brief Divide two finite values without executing a zero or overflowing quotient. */
inline bool safeDivide( real64 const numerator,
                        real64 const denominator,
                        real64 & result )
{
  if( !finiteValue( numerator ) || !finiteValue( denominator ) ||
      isZero(denominator) )
  {
    return false;
  }

  real64 const absNumerator = std::abs( numerator );
  real64 const absDenominator = std::abs( denominator );
  if( absDenominator < 1.0 &&
      absNumerator > std::numeric_limits< real64 >::max() * absDenominator )
  {
    return false;
  }
  result = numerator / denominator;
  return finiteValue( result );
}

/** @brief Compute a three-vector norm without overflowing its sum of squares. */
inline bool robustNorm( std::array< real64, 3 > const & value,
                        real64 & result )
{
  if( !finiteVector( value ) )
  {
    return false;
  }
  real64 const scale = std::max( std::abs( value[0] ),
                                 std::max( std::abs( value[1] ),
                                           std::abs( value[2] ) ) );
  if( isZero(scale) )
  {
    result = 0.0;
    return true;
  }

  real64 const x = value[0] / scale;
  real64 const y = value[1] / scale;
  real64 const z = value[2] / scale;
  real64 const scaledNorm = std::sqrt( x * x + y * y + z * z );
  return safeMultiply( scale, scaledNorm, result );
}

/** @brief Compute a checked dot product. */
inline bool safeDot( std::array< real64, 3 > const & a,
                     std::array< real64, 3 > const & b,
                     real64 & result )
{
  result = 0.0;
  for( localIndex component = 0; component < 3; ++component )
  {
    real64 product = 0.0;
    real64 updated = 0.0;
    if( !safeMultiply( a[component], b[component], product ) ||
        !safeAdd( result, product, updated ) )
    {
      return false;
    }
    result = updated;
  }
  return true;
}

/** @brief Compute the two-body effective mass without reciprocal overflow. */
inline bool effectiveMass( real64 const massA,
                           real64 const massB,
                           real64 & result )
{
  if( !finiteValue( massA ) || !finiteValue( massB ) ||
      !( massA > 0.0 ) || !( massB > 0.0 ) )
  {
    return false;
  }
  real64 const smaller = std::min( massA, massB );
  real64 const larger = std::max( massA, massB );
  real64 const denominator = 1.0 + smaller / larger;
  return safeDivide( smaller, denominator, result ) && result > 0.0;
}

/**
 * @brief Compute the node-local mass cutoff used to construct contact topology.
 *
 * The global MPM small-mass threshold continues to govern ordinary grid-field
 * activity.  The two additional terms are contact-only conditioning controls:
 * an absolute floor and a fraction of the largest otherwise contact-capable
 * field on the node.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 contactMassCutoff( real64 const smallMass,
                          real64 const contactMinimumMass,
                          real64 const contactMinimumMassFraction,
                          real64 const nodeMaximumContactMass )
{
  real64 cutoff = smallMass > contactMinimumMass
    ? smallMass
    : contactMinimumMass;
  real64 const relativeCutoff =
    contactMinimumMassFraction * nodeMaximumContactMass;
  return cutoff > relativeCutoff ? cutoff : relativeCutoff;
}

/** @brief Select the optional pairwise LR normal policy for one node. */
inline bool useMultifieldLogisticRegression( int const enabled,
                                             localIndex const activeContactFieldCount,
                                             bool const rigidBodyNode )
{
  return enabled != 0 && activeContactFieldCount > 2 && !rigidBodyNode;
}

/**
 * @brief Select a pair for implicit coupled contact from its geometric gap.
 *
 * The initial normal velocity is deliberately not part of this candidate
 * test. At a multifield node, an impulse from another pair can reverse the
 * pair velocity. The unilateral projection inside either coupled solver
 * decides whether the final normal impulse is positive or zero.
 */
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
bool implicitGapConstraintCandidate(
  real64 const gap,
  real64 const gapActivationTolerance )
{
  return gap <= gapActivationTolerance;
}

inline real64 dot( std::array< real64, 3 > const & a,
                   std::array< real64, 3 > const & b )
{
  real64 result = 0.0;
  return safeDot( a, b, result ) ? result : 0.0;
}

inline real64 norm( std::array< real64, 3 > const & a )
{
  real64 result = 0.0;
  return robustNorm( a, result )
    ? result
    : std::numeric_limits< real64 >::max();
}

inline bool tryRelativeVelocity(
  std::vector< std::array< real64, 3 > > const & velocity,
  ProjectedGaussSeidelContactConstraint const & constraint,
  std::array< real64, 3 > & relative )
{
  if( constraint.fieldA < 0 || constraint.fieldB < 0 ||
      static_cast< std::size_t >( constraint.fieldA ) >= velocity.size() ||
      static_cast< std::size_t >( constraint.fieldB ) >= velocity.size() )
  {
    return false;
  }
  for( localIndex component = 0; component < 3; ++component )
  {
    if( !safeSubtract( velocity[constraint.fieldB][component],
                       velocity[constraint.fieldA][component],
                       relative[component] ) )
    {
      return false;
    }
  }
  return true;
}

inline std::array< real64, 3 > relativeVelocity(
  std::vector< std::array< real64, 3 > > const & velocity,
  ProjectedGaussSeidelContactConstraint const & constraint )
{
  std::array< real64, 3 > relative = {{ 0.0, 0.0, 0.0 }};
  tryRelativeVelocity( velocity, constraint, relative );
  return relative;
}

inline bool applyImpulseToA( std::array< real64, 3 > const & impulseOnA,
                             ProjectedGaussSeidelContactConstraint const & constraint,
                             std::vector< real64 > const & mass,
                             std::vector< std::array< real64, 3 > > & velocity )
{
  if( constraint.fieldA < 0 || constraint.fieldB < 0 ||
      static_cast< std::size_t >( constraint.fieldA ) >= mass.size() ||
      static_cast< std::size_t >( constraint.fieldB ) >= mass.size() ||
      static_cast< std::size_t >( constraint.fieldA ) >= velocity.size() ||
      static_cast< std::size_t >( constraint.fieldB ) >= velocity.size() ||
      !finiteVector( impulseOnA ) )
  {
    return false;
  }

  std::array< real64, 3 > candidateA = velocity[constraint.fieldA];
  std::array< real64, 3 > candidateB = velocity[constraint.fieldB];
  for( localIndex component = 0; component < 3; ++component )
  {
    real64 incrementA = 0.0;
    real64 incrementB = 0.0;
    if( !safeDivide( impulseOnA[component],
                     mass[constraint.fieldA],
                     incrementA ) ||
        !safeDivide( impulseOnA[component],
                     mass[constraint.fieldB],
                     incrementB ) ||
        !safeAdd( candidateA[component], incrementA,
                  candidateA[component] ) ||
        !safeSubtract( candidateB[component], incrementB,
                       candidateB[component] ) )
    {
      return false;
    }
  }
  velocity[constraint.fieldA] = candidateA;
  velocity[constraint.fieldB] = candidateB;
  return true;
}

/**
 * @brief Return the accumulated tangential impulse projected onto the Coulomb disk.
 *
 * The returned vector is the impulse on field A. It is kept in the plane
 * orthogonal to the pair normal and satisfies
 *
 *   ||J_t|| <= mu max( 0, J_n - J_n,bias ).
 */
inline bool projectedTangentialImpulse(
  ProjectedGaussSeidelContactConstraint const & constraint,
  std::array< real64, 3 > const & relativeVelocityAB,
  real64 const effectiveMass,
  real64 const normalImpulse,
  real64 const relaxation,
  std::array< real64, 3 > & projectedImpulse )
{
  real64 normalVelocity = 0.0;
  if( !safeDot( relativeVelocityAB, constraint.normal, normalVelocity ) )
  {
    return false;
  }
  std::array< real64, 3 > tangentVelocity = {{ 0.0, 0.0, 0.0 }};
  for( localIndex component = 0; component < 3; ++component )
  {
    real64 normalComponent = 0.0;
    if( !safeMultiply( normalVelocity,
                       constraint.normal[component],
                       normalComponent ) ||
        !safeSubtract( relativeVelocityAB[component],
                       normalComponent,
                       tangentVelocity[component] ) )
    {
      return false;
    }
  }

  std::array< real64, 3 > trial = constraint.accumulatedTangentialImpulse;
  real64 relaxedMass = 0.0;
  if( !safeMultiply( relaxation, effectiveMass, relaxedMass ) )
  {
    return false;
  }
  for( localIndex component = 0; component < 3; ++component )
  {
    // A positive impulse on A decreases vB-vA, hence the positive sign.
    real64 increment = 0.0;
    if( !safeMultiply( relaxedMass, tangentVelocity[component], increment ) ||
        !safeAdd( trial[component], increment, trial[component] ) )
    {
      return false;
    }
  }

  // Remove roundoff-level normal drift before evaluating the disk radius.
  real64 trialNormalComponent = 0.0;
  if( !safeDot( trial, constraint.normal, trialNormalComponent ) )
  {
    return false;
  }
  for( localIndex component = 0; component < 3; ++component )
  {
    real64 normalDrift = 0.0;
    if( !safeMultiply( trialNormalComponent,
                       constraint.normal[component],
                       normalDrift ) ||
        !safeSubtract( trial[component], normalDrift, trial[component] ) )
    {
      return false;
    }
  }

  // The legacy contact law excludes prescribed overlap/penetration bias from
  // the Coulomb bound. This also prevents a cohesive penetration-only
  // constraint from acquiring shear traction.
  real64 unbiasedNormalImpulse = 0.0;
  if( !safeSubtract( normalImpulse,
                     constraint.normalBiasImpulse,
                     unbiasedNormalImpulse ) )
  {
    return false;
  }
  real64 const frictionNormalImpulse = std::max( 0.0,
                                                  unbiasedNormalImpulse );
  real64 radius = 0.0;
  real64 trialNorm = 0.0;
  if( !safeMultiply( constraint.frictionCoefficient,
                     frictionNormalImpulse,
                     radius ) ||
      !robustNorm( trial, trialNorm ) )
  {
    return false;
  }
  if( trialNorm > radius && trialNorm > 0.0 )
  {
    real64 scale = 0.0;
    if( !safeDivide( radius, trialNorm, scale ) )
    {
      return false;
    }
    for( localIndex component = 0; component < 3; ++component )
    {
      if( !safeMultiply( trial[component], scale, trial[component] ) )
      {
        return false;
      }
    }
  }
  projectedImpulse = trial;
  return finiteVector( projectedImpulse );
}

inline bool projectedResidual(
  std::vector< real64 > const & mass,
  std::vector< std::array< real64, 3 > > const & velocity,
  std::vector< ProjectedGaussSeidelContactConstraint > const & constraints,
  real64 & residual )
{
  residual = 0.0;
  for( ProjectedGaussSeidelContactConstraint const & constraint : constraints )
  {
    real64 pairEffectiveMass = 0.0;
    std::array< real64, 3 > relative = {{ 0.0, 0.0, 0.0 }};
    if( constraint.fieldA < 0 || constraint.fieldB < 0 ||
        static_cast< std::size_t >( constraint.fieldA ) >= mass.size() ||
        static_cast< std::size_t >( constraint.fieldB ) >= mass.size() ||
        !effectiveMass( mass[constraint.fieldA],
                        mass[constraint.fieldB],
                        pairEffectiveMass ) ||
        !tryRelativeVelocity( velocity, constraint, relative ) )
    {
      return false;
    }

    if( constraint.bilateral )
    {
      real64 relativeNorm = 0.0;
      if( !robustNorm( relative, relativeNorm ) )
      {
        return false;
      }
      residual = std::max( residual, relativeNorm );
      continue;
    }

    real64 normalVelocity = 0.0;
    real64 targetError = 0.0;
    real64 normalCorrection = 0.0;
    real64 trialNormalImpulse = 0.0;
    if( !safeDot( relative, constraint.normal, normalVelocity ) ||
        !safeSubtract( constraint.targetNormalVelocity,
                       normalVelocity,
                       targetError ) ||
        !safeMultiply( pairEffectiveMass,
                       targetError,
                       normalCorrection ) ||
        !safeAdd( constraint.accumulatedNormalImpulse,
                  normalCorrection,
                  trialNormalImpulse ) )
    {
      return false;
    }
    real64 const projectedNormalImpulse = std::max( 0.0,
                                                     trialNormalImpulse );
    real64 normalImpulseDifference = 0.0;
    real64 normalVelocityResidual = 0.0;
    if( !safeSubtract( projectedNormalImpulse,
                       constraint.accumulatedNormalImpulse,
                       normalImpulseDifference ) ||
        !safeDivide( std::abs( normalImpulseDifference ),
                     pairEffectiveMass,
                     normalVelocityResidual ) )
    {
      return false;
    }
    residual = std::max( residual, normalVelocityResidual );

    std::array< real64, 3 > projectedTangent = {{ 0.0, 0.0, 0.0 }};
    std::array< real64, 3 > tangentDifference = {{ 0.0, 0.0, 0.0 }};
    if( !projectedTangentialImpulse( constraint,
                                    relative,
                                    pairEffectiveMass,
                                    projectedNormalImpulse,
                                    1.0,
                                    projectedTangent ) )
    {
      return false;
    }
    for( localIndex component = 0; component < 3; ++component )
    {
      if( !safeSubtract( projectedTangent[component],
                         constraint.accumulatedTangentialImpulse[component],
                         tangentDifference[component] ) )
      {
        return false;
      }
    }
    real64 tangentImpulseResidual = 0.0;
    real64 tangentVelocityResidual = 0.0;
    if( !robustNorm( tangentDifference, tangentImpulseResidual ) ||
        !safeDivide( tangentImpulseResidual,
                     pairEffectiveMass,
                     tangentVelocityResidual ) )
    {
      return false;
    }
    residual = std::max( residual, tangentVelocityResidual );
  }
  return finiteValue( residual );
}

/**
 * @brief Solve all active contact constraints for one grid node.
 *
 * Storage is supplied by the caller and is proportional to the number of
 * fields and constraints that are actually active at this node. In
 * particular, no node-by-global-contact-group or global pair matrix is built.
 */
inline ProjectedGaussSeidelContactResult solve(
  std::vector< real64 > const & mass,
  std::vector< std::array< real64, 3 > > & velocity,
  std::vector< ProjectedGaussSeidelContactConstraint > & constraints,
  integer const maximumIterations,
  real64 const velocityTolerance,
  real64 const relaxation )
{
  FloatingPointEnvironmentGuard floatingPointEnvironmentGuard;
  ProjectedGaussSeidelContactResult result;
  if( constraints.empty() )
  {
    return result;
  }

  std::vector< std::array< real64, 3 > > const initialVelocity = velocity;
  auto rejectNode = [&]()
  {
    velocity = initialVelocity;
    for( ProjectedGaussSeidelContactConstraint & constraint : constraints )
    {
      constraint.accumulatedNormalImpulse = 0.0;
      constraint.accumulatedTangentialImpulse = {{ 0.0, 0.0, 0.0 }};
    }
    result.converged = false;
    result.numericalFailure = true;
    result.usedSafeFallback = true;
    ++result.numericalGuardActivations;
    result.residual = std::numeric_limits< real64 >::max();
  };

  if( mass.size() != velocity.size() || maximumIterations <= 0 ||
      !finiteValue( velocityTolerance ) || velocityTolerance < 0.0 ||
      !finiteValue( relaxation ) || !( relaxation > 0.0 ) )
  {
    rejectNode();
    return result;
  }

  for( std::array< real64, 3 > const & fieldVelocity : velocity )
  {
    if( !finiteVector( fieldVelocity ) )
    {
      rejectNode();
      return result;
    }
  }

  for( ProjectedGaussSeidelContactConstraint & constraint : constraints )
  {
    if( constraint.fieldA < 0 || constraint.fieldB < 0 ||
        constraint.fieldA == constraint.fieldB ||
        static_cast< std::size_t >( constraint.fieldA ) >= mass.size() ||
        static_cast< std::size_t >( constraint.fieldB ) >= mass.size() ||
        !finiteValue( constraint.targetNormalVelocity ) ||
        !finiteValue( constraint.frictionCoefficient ) ||
        constraint.frictionCoefficient < 0.0 ||
        !finiteValue( constraint.normalBiasImpulse ) ||
        constraint.normalBiasImpulse < 0.0 )
    {
      rejectNode();
      return result;
    }

    real64 pairEffectiveMass = 0.0;
    real64 normalNorm = 0.0;
    if( !effectiveMass( mass[constraint.fieldA],
                        mass[constraint.fieldB],
                        pairEffectiveMass ) ||
        !robustNorm( constraint.normal, normalNorm ) ||
        !( normalNorm > 1.0e-20 ) )
    {
      rejectNode();
      return result;
    }
    for( localIndex component = 0; component < 3; ++component )
    {
      if( !safeDivide( constraint.normal[component],
                       normalNorm,
                       constraint.normal[component] ) )
      {
        rejectNode();
        return result;
      }
    }
    constraint.accumulatedNormalImpulse = 0.0;
    constraint.accumulatedTangentialImpulse = {{ 0.0, 0.0, 0.0 }};
  }

  result.converged = false;
  for( integer iteration = 0; iteration < maximumIterations; ++iteration )
  {
    std::vector< std::array< real64, 3 > > const previousVelocity = velocity;
    std::vector< ProjectedGaussSeidelContactConstraint > const previousConstraints = constraints;
    bool iterationSafe = true;
    for( ProjectedGaussSeidelContactConstraint & constraint : constraints )
    {
      real64 pairEffectiveMass = 0.0;
      std::array< real64, 3 > relative = {{ 0.0, 0.0, 0.0 }};
      if( !effectiveMass( mass[constraint.fieldA],
                          mass[constraint.fieldB],
                          pairEffectiveMass ) ||
          !tryRelativeVelocity( velocity, constraint, relative ) )
      {
        iterationSafe = false;
        break;
      }

      if( constraint.bilateral )
      {
        real64 relaxedMass = 0.0;
        std::array< real64, 3 > impulseOnA = {{ 0.0, 0.0, 0.0 }};
        if( !safeMultiply( relaxation, pairEffectiveMass, relaxedMass ) )
        {
          iterationSafe = false;
          break;
        }
        for( localIndex component = 0; component < 3; ++component )
        {
          if( !safeMultiply( relaxedMass,
                             relative[component],
                             impulseOnA[component] ) )
          {
            iterationSafe = false;
            break;
          }
        }
        if( !iterationSafe ||
            !applyImpulseToA( impulseOnA, constraint, mass, velocity ) )
        {
          iterationSafe = false;
          break;
        }
        continue;
      }

      real64 normalVelocity = 0.0;
      real64 const oldNormalImpulse = constraint.accumulatedNormalImpulse;
      real64 targetError = 0.0;
      real64 relaxedMass = 0.0;
      real64 normalCorrection = 0.0;
      real64 trialNormalImpulse = 0.0;
      if( !safeDot( relative, constraint.normal, normalVelocity ) ||
          !safeSubtract( constraint.targetNormalVelocity,
                         normalVelocity,
                         targetError ) ||
          !safeMultiply( relaxation, pairEffectiveMass, relaxedMass ) ||
          !safeMultiply( relaxedMass, targetError, normalCorrection ) ||
          !safeAdd( oldNormalImpulse,
                    normalCorrection,
                    trialNormalImpulse ) )
      {
        iterationSafe = false;
        break;
      }
      real64 const newNormalImpulse = std::max( 0.0,
                                                trialNormalImpulse );
      real64 normalImpulseIncrement = 0.0;
      if( !safeSubtract( newNormalImpulse,
                         oldNormalImpulse,
                         normalImpulseIncrement ) )
      {
        iterationSafe = false;
        break;
      }
      constraint.accumulatedNormalImpulse = newNormalImpulse;

      std::array< real64, 3 > normalImpulseOnA = {{ 0.0, 0.0, 0.0 }};
      for( localIndex component = 0; component < 3; ++component )
      {
        if( !safeMultiply( -normalImpulseIncrement,
                           constraint.normal[component],
                           normalImpulseOnA[component] ) )
        {
          iterationSafe = false;
          break;
        }
      }
      if( !iterationSafe ||
          !applyImpulseToA( normalImpulseOnA,
                            constraint,
                            mass,
                            velocity ) )
      {
        iterationSafe = false;
        break;
      }

      if( !tryRelativeVelocity( velocity, constraint, relative ) )
      {
        iterationSafe = false;
        break;
      }
      std::array< real64, 3 > newTangentialImpulse = {{ 0.0, 0.0, 0.0 }};
      std::array< real64, 3 > tangentialImpulseIncrement = {{ 0.0, 0.0, 0.0 }};
      if( !projectedTangentialImpulse( constraint,
                                      relative,
                                      pairEffectiveMass,
                                      newNormalImpulse,
                                      relaxation,
                                      newTangentialImpulse ) )
      {
        iterationSafe = false;
        break;
      }
      for( localIndex component = 0; component < 3; ++component )
      {
        if( !safeSubtract( newTangentialImpulse[component],
                           constraint.accumulatedTangentialImpulse[component],
                           tangentialImpulseIncrement[component] ) )
        {
          iterationSafe = false;
          break;
        }
      }
      if( !iterationSafe )
      {
        break;
      }
      constraint.accumulatedTangentialImpulse = newTangentialImpulse;
      if( !applyImpulseToA( tangentialImpulseIncrement,
                            constraint,
                            mass,
                            velocity ) )
      {
        iterationSafe = false;
        break;
      }
    }

    if( !iterationSafe )
    {
      velocity = previousVelocity;
      constraints = previousConstraints;
      result.numericalFailure = true;
      result.usedSafeFallback = true;
      ++result.numericalGuardActivations;
      result.residual = std::numeric_limits< real64 >::max();
      break;
    }

    result.iterations = iteration + 1;
    if( !projectedResidual( mass,
                            velocity,
                            constraints,
                            result.residual ) )
    {
      velocity = previousVelocity;
      constraints = previousConstraints;
      result.numericalFailure = true;
      result.usedSafeFallback = true;
      ++result.numericalGuardActivations;
      result.residual = std::numeric_limits< real64 >::max();
      break;
    }
    if( result.residual <= velocityTolerance )
    {
      result.converged = true;
      break;
    }
  }

  return result;
}

inline bool suppressContactForCohesivePair( int const cohesiveFieldA,
                                            int const cohesiveFieldB,
                                            int const preventInterpenetration )
{
  return cohesiveFieldA != 0 && cohesiveFieldB != 0 &&
         preventInterpenetration != 1;
}

} // namespace projectedGaussSeidelContact
} // namespace mpm
} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_PROJECTEDGAUSSSEIDELCONTACT_HPP_
