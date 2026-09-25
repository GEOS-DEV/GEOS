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
 * @file NewtonRaphsonContact.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_NEWTONRAPHSONCONTACT_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_NEWTONRAPHSONCONTACT_HPP_

#include "ProjectedGaussSeidelContact.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace geos
{
namespace mpm
{

/** @brief Convergence information returned by the node-local Newton solve. */
struct NewtonRaphsonContactResult
{
  integer iterations = 0;
  integer regularizedSteps = 0;
  integer lineSearchReductions = 0;
  integer numericalGuardActivations = 0;
  real64 residual = 0.0;
  bool converged = true;
  bool numericalFailure = false;
  bool usedSafeFallback = false;
};

namespace newtonRaphsonContact
{

using Constraint = ProjectedGaussSeidelContactConstraint;
using Vector3 = std::array< real64, 3 >;

/**
 * @brief Reconstruct field velocities from simultaneous pair impulses.
 *
 * Each three-component unknown is the total impulse on field A for one
 * constraint. Equal-and-opposite impulses are applied to field B.
 */
inline bool reconstructVelocity( std::vector< real64 > const & mass,
                                 std::vector< Vector3 > const & initialVelocity,
                                 std::vector< Constraint > const & constraints,
                                 std::vector< real64 > const & impulse,
                                 std::vector< Vector3 > & velocity )
{
  if( mass.size() != initialVelocity.size() ||
      impulse.size() != 3 * constraints.size() )
  {
    return false;
  }
  velocity = initialVelocity;
  for( std::size_t constraintIndex = 0;
       constraintIndex < constraints.size();
       ++constraintIndex )
  {
    std::size_t const offset = 3 * constraintIndex;
    Vector3 const impulseOnA = {{ impulse[offset],
                                  impulse[offset + 1],
                                  impulse[offset + 2] }};
    if( !projectedGaussSeidelContact::applyImpulseToA(
          impulseOnA,
          constraints[constraintIndex],
          mass,
          velocity ) )
    {
      return false;
    }
  }
  return true;
}

/**
 * @brief Evaluate the semismooth contact fixed-point residual.
 *
 * Bilateral constraints use the common-velocity residual. Unilateral
 * constraints use the normal positive-half-line projection and the
 * tangential Coulomb-disk projection. The residual has impulse units so all
 * three components of every constraint can be solved simultaneously.
 */
inline bool evaluateResidual( std::vector< real64 > const & mass,
                              std::vector< Vector3 > const & initialVelocity,
                              std::vector< Constraint > const & constraints,
                              std::vector< real64 > const & impulse,
                              std::vector< real64 > & residual,
                              std::vector< Vector3 > & velocity )
{
  if( !reconstructVelocity( mass,
                            initialVelocity,
                            constraints,
                            impulse,
                            velocity ) )
  {
    return false;
  }
  residual.assign( impulse.size(), 0.0 );

  for( std::size_t constraintIndex = 0;
       constraintIndex < constraints.size();
       ++constraintIndex )
  {
    Constraint const & constraint = constraints[constraintIndex];
    std::size_t const offset = 3 * constraintIndex;
    real64 pairEffectiveMass = 0.0;
    Vector3 relative = {{ 0.0, 0.0, 0.0 }};
    if( constraint.fieldA < 0 || constraint.fieldB < 0 ||
        static_cast< std::size_t >( constraint.fieldA ) >= mass.size() ||
        static_cast< std::size_t >( constraint.fieldB ) >= mass.size() ||
        !projectedGaussSeidelContact::effectiveMass(
          mass[constraint.fieldA],
          mass[constraint.fieldB],
          pairEffectiveMass ) ||
        !projectedGaussSeidelContact::tryRelativeVelocity(
          velocity,
          constraint,
          relative ) )
    {
      return false;
    }

    if( constraint.bilateral )
    {
      for( localIndex component = 0; component < 3; ++component )
      {
        if( !projectedGaussSeidelContact::safeMultiply(
              pairEffectiveMass,
              relative[component],
              residual[offset + component] ) )
        {
          return false;
        }
      }
      continue;
    }

    Vector3 const impulseOnA = {{ impulse[offset],
                                  impulse[offset + 1],
                                  impulse[offset + 2] }};
    real64 impulseNormalComponent = 0.0;
    if( !projectedGaussSeidelContact::safeDot(
          impulseOnA,
          constraint.normal,
          impulseNormalComponent ) )
    {
      return false;
    }
    real64 const normalImpulse = -impulseNormalComponent;
    Vector3 tangentialImpulse = {{ 0.0, 0.0, 0.0 }};
    for( localIndex component = 0; component < 3; ++component )
    {
      real64 normalComponent = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            impulseNormalComponent,
            constraint.normal[component],
            normalComponent ) ||
          !projectedGaussSeidelContact::safeSubtract(
            impulseOnA[component],
            normalComponent,
            tangentialImpulse[component] ) )
      {
        return false;
      }
    }

    real64 normalVelocity = 0.0;
    real64 targetError = 0.0;
    real64 normalCorrection = 0.0;
    real64 trialNormalImpulse = 0.0;
    if( !projectedGaussSeidelContact::safeDot(
          relative,
          constraint.normal,
          normalVelocity ) ||
        !projectedGaussSeidelContact::safeSubtract(
          constraint.targetNormalVelocity,
          normalVelocity,
          targetError ) ||
        !projectedGaussSeidelContact::safeMultiply(
          pairEffectiveMass,
          targetError,
          normalCorrection ) ||
        !projectedGaussSeidelContact::safeAdd(
          normalImpulse,
          normalCorrection,
          trialNormalImpulse ) )
    {
      return false;
    }
    real64 const projectedNormalImpulse = std::max( 0.0,
                                                     trialNormalImpulse );
    real64 normalResidual = 0.0;
    if( !projectedGaussSeidelContact::safeSubtract(
          normalImpulse,
          projectedNormalImpulse,
          normalResidual ) )
    {
      return false;
    }

    Vector3 tangentVelocity = {{ 0.0, 0.0, 0.0 }};
    for( localIndex component = 0; component < 3; ++component )
    {
      real64 normalComponent = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            normalVelocity,
            constraint.normal[component],
            normalComponent ) ||
          !projectedGaussSeidelContact::safeSubtract(
            relative[component],
            normalComponent,
            tangentVelocity[component] ) )
      {
        return false;
      }
    }
    Vector3 projectedTangentialImpulse = tangentialImpulse;
    for( localIndex component = 0; component < 3; ++component )
    {
      real64 correction = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            pairEffectiveMass,
            tangentVelocity[component],
            correction ) ||
          !projectedGaussSeidelContact::safeAdd(
            projectedTangentialImpulse[component],
            correction,
            projectedTangentialImpulse[component] ) )
      {
        return false;
      }
    }

    real64 projectedTangentNormal = 0.0;
    if( !projectedGaussSeidelContact::safeDot(
          projectedTangentialImpulse,
          constraint.normal,
          projectedTangentNormal ) )
    {
      return false;
    }
    for( localIndex component = 0; component < 3; ++component )
    {
      real64 normalDrift = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            projectedTangentNormal,
            constraint.normal[component],
            normalDrift ) ||
          !projectedGaussSeidelContact::safeSubtract(
            projectedTangentialImpulse[component],
            normalDrift,
            projectedTangentialImpulse[component] ) )
      {
        return false;
      }
    }

    real64 unbiasedNormalImpulse = 0.0;
    real64 frictionRadius = 0.0;
    real64 projectedTangentNorm = 0.0;
    if( !projectedGaussSeidelContact::safeSubtract(
          projectedNormalImpulse,
          constraint.normalBiasImpulse,
          unbiasedNormalImpulse ) ||
        !projectedGaussSeidelContact::safeMultiply(
          constraint.frictionCoefficient,
          std::max( 0.0, unbiasedNormalImpulse ),
          frictionRadius ) ||
        !projectedGaussSeidelContact::robustNorm(
          projectedTangentialImpulse,
          projectedTangentNorm ) )
    {
      return false;
    }
    if( projectedTangentNorm > frictionRadius &&
        projectedTangentNorm > 0.0 )
    {
      real64 scale = 0.0;
      if( !projectedGaussSeidelContact::safeDivide(
            frictionRadius,
            projectedTangentNorm,
            scale ) )
      {
        return false;
      }
      for( localIndex component = 0; component < 3; ++component )
      {
        if( !projectedGaussSeidelContact::safeMultiply(
              projectedTangentialImpulse[component],
              scale,
              projectedTangentialImpulse[component] ) )
        {
          return false;
        }
      }
    }

    for( localIndex component = 0; component < 3; ++component )
    {
      real64 normalPart = 0.0;
      real64 combined = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            -normalResidual,
            constraint.normal[component],
            normalPart ) ||
          !projectedGaussSeidelContact::safeAdd(
            normalPart,
            tangentialImpulse[component],
            combined ) ||
          !projectedGaussSeidelContact::safeSubtract(
            combined,
            projectedTangentialImpulse[component],
            residual[offset + component] ) )
      {
        return false;
      }
    }
  }
  return true;
}

/** @brief Return the maximum fixed-point residual in velocity units. */
inline bool velocityResidualNorm(
  std::vector< real64 > const & mass,
  std::vector< Constraint > const & constraints,
  std::vector< real64 > const & residual,
  real64 & result )
{
  result = 0.0;
  if( residual.size() != 3 * constraints.size() )
  {
    return false;
  }
  for( std::size_t constraintIndex = 0;
       constraintIndex < constraints.size();
       ++constraintIndex )
  {
    Constraint const & constraint = constraints[constraintIndex];
    std::size_t const offset = 3 * constraintIndex;
    Vector3 const block = {{ residual[offset],
                             residual[offset + 1],
                             residual[offset + 2] }};
    real64 pairEffectiveMass = 0.0;
    real64 blockNorm = 0.0;
    real64 blockVelocityNorm = 0.0;
    if( constraint.fieldA < 0 || constraint.fieldB < 0 ||
        static_cast< std::size_t >( constraint.fieldA ) >= mass.size() ||
        static_cast< std::size_t >( constraint.fieldB ) >= mass.size() ||
        !projectedGaussSeidelContact::effectiveMass(
          mass[constraint.fieldA],
          mass[constraint.fieldB],
          pairEffectiveMass ) ||
        !projectedGaussSeidelContact::robustNorm( block, blockNorm ) ||
        !projectedGaussSeidelContact::safeDivide(
          blockNorm,
          pairEffectiveMass,
          blockVelocityNorm ) )
    {
      return false;
    }
    result = std::max( result, blockVelocityNorm );
  }
  return projectedGaussSeidelContact::finiteValue( result );
}

/** @brief Return the Euclidean norm used by the Newton line search. */
inline bool meritNorm( std::vector< real64 > const & residual,
                       real64 & result )
{
  real64 scale = 0.0;
  for( real64 const component : residual )
  {
    if( !projectedGaussSeidelContact::finiteValue( component ) )
    {
      return false;
    }
    scale = std::max( scale, std::abs( component ) );
  }
  if( isZero(scale) )
  {
    result = 0.0;
    return true;
  }

  real64 sum = 0.0;
  for( real64 const component : residual )
  {
    real64 const scaled = component / scale;
    real64 square = 0.0;
    if( !projectedGaussSeidelContact::safeMultiply(
          scaled,
          scaled,
          square ) ||
        !projectedGaussSeidelContact::safeAdd( sum, square, sum ) )
    {
      return false;
    }
  }
  return projectedGaussSeidelContact::safeMultiply(
    scale,
    std::sqrt( sum ),
    result );
}

/**
 * @brief Solve a dense square system by Gaussian elimination with pivoting.
 */
inline bool solveDenseLinearSystem( std::vector< real64 > matrix,
                                    std::vector< real64 > rightHandSide,
                                    real64 const relativePivotTolerance,
                                    std::vector< real64 > & solution )
{
  std::size_t const size = rightHandSide.size();
  solution.assign( size, 0.0 );
  if( size == 0 )
  {
    return true;
  }
  if( matrix.size() != size * size ||
      !projectedGaussSeidelContact::finiteValue(
        relativePivotTolerance ) ||
      relativePivotTolerance < 0.0 )
  {
    return false;
  }

  real64 matrixScale = 0.0;
  for( real64 const value : matrix )
  {
    if( !projectedGaussSeidelContact::finiteValue( value ) )
    {
      return false;
    }
    matrixScale = std::max( matrixScale, std::abs( value ) );
  }
  if( isZero( matrixScale ) ||
      !projectedGaussSeidelContact::finiteValue( matrixScale ) )
  {
    return false;
  }
  for( real64 & value : matrix )
  {
    if( !projectedGaussSeidelContact::safeDivide(
          value,
          matrixScale,
          value ) )
    {
      return false;
    }
  }
  for( real64 & value : rightHandSide )
  {
    if( !projectedGaussSeidelContact::safeDivide(
          value,
          matrixScale,
          value ) )
    {
      return false;
    }
  }
  real64 const pivotTolerance = relativePivotTolerance;

  for( std::size_t column = 0; column < size; ++column )
  {
    std::size_t pivotRow = column;
    real64 pivotMagnitude = std::abs( matrix[column * size + column] );
    for( std::size_t row = column + 1; row < size; ++row )
    {
      real64 const candidate = std::abs( matrix[row * size + column] );
      if( candidate > pivotMagnitude )
      {
        pivotMagnitude = candidate;
        pivotRow = row;
      }
    }
    if( !projectedGaussSeidelContact::finiteValue( pivotMagnitude ) ||
        pivotMagnitude <= pivotTolerance )
    {
      return false;
    }

    if( pivotRow != column )
    {
      for( std::size_t entry = column; entry < size; ++entry )
      {
        std::swap( matrix[column * size + entry],
                   matrix[pivotRow * size + entry] );
      }
      std::swap( rightHandSide[column], rightHandSide[pivotRow] );
    }

    real64 const pivot = matrix[column * size + column];
    for( std::size_t row = column + 1; row < size; ++row )
    {
      real64 factor = 0.0;
      if( !projectedGaussSeidelContact::safeDivide(
            matrix[row * size + column],
            pivot,
            factor ) )
      {
        return false;
      }
      matrix[row * size + column] = 0.0;
      for( std::size_t entry = column + 1; entry < size; ++entry )
      {
        real64 product = 0.0;
        if( !projectedGaussSeidelContact::safeMultiply(
              factor,
              matrix[column * size + entry],
              product ) ||
            !projectedGaussSeidelContact::safeSubtract(
              matrix[row * size + entry],
              product,
              matrix[row * size + entry] ) )
        {
          return false;
        }
      }
      real64 rhsProduct = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            factor,
            rightHandSide[column],
            rhsProduct ) ||
          !projectedGaussSeidelContact::safeSubtract(
            rightHandSide[row],
            rhsProduct,
            rightHandSide[row] ) )
      {
        return false;
      }
    }
  }

  for( std::size_t reverseRow = size; reverseRow-- > 0; )
  {
    real64 value = rightHandSide[reverseRow];
    for( std::size_t column = reverseRow + 1; column < size; ++column )
    {
      real64 product = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            matrix[reverseRow * size + column],
            solution[column],
            product ) ||
          !projectedGaussSeidelContact::safeSubtract(
            value,
            product,
            value ) )
      {
        return false;
      }
    }
    real64 const diagonal = matrix[reverseRow * size + reverseRow];
    if( !projectedGaussSeidelContact::finiteValue( diagonal ) ||
        std::abs( diagonal ) <= pivotTolerance )
    {
      return false;
    }
    if( !projectedGaussSeidelContact::safeDivide(
          value,
          diagonal,
          solution[reverseRow] ) )
    {
      return false;
    }
  }
  return true;
}

/** @brief Compute a regularized least-squares Newton direction. */
inline bool regularizedNewtonStep( std::vector< real64 > const & jacobian,
                                   std::vector< real64 > const & residual,
                                   real64 const regularization,
                                   std::vector< real64 > & step )
{
  std::size_t const size = residual.size();
  if( jacobian.size() != size * size ||
      !projectedGaussSeidelContact::finiteValue( regularization ) ||
      regularization < 0.0 )
  {
    return false;
  }

  real64 commonScale = 0.0;
  for( real64 const value : jacobian )
  {
    if( !projectedGaussSeidelContact::finiteValue( value ) )
    {
      return false;
    }
    commonScale = std::max( commonScale, std::abs( value ) );
  }
  for( real64 const value : residual )
  {
    if( !projectedGaussSeidelContact::finiteValue( value ) )
    {
      return false;
    }
    commonScale = std::max( commonScale, std::abs( value ) );
  }
  if( isZero(commonScale) )
  {
    return false;
  }

  std::vector< real64 > normalMatrix( size * size, 0.0 );
  std::vector< real64 > rightHandSide( size, 0.0 );
  for( std::size_t row = 0; row < size; ++row )
  {
    for( std::size_t column = 0; column < size; ++column )
    {
      real64 value = 0.0;
      real64 scaledResidual = 0.0;
      if( !projectedGaussSeidelContact::safeDivide(
            jacobian[row * size + column],
            commonScale,
            value ) ||
          !projectedGaussSeidelContact::safeDivide(
            residual[row],
            commonScale,
            scaledResidual ) )
      {
        return false;
      }
      real64 rhsContribution = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            value,
            scaledResidual,
            rhsContribution ) ||
          !projectedGaussSeidelContact::safeSubtract(
            rightHandSide[column],
            rhsContribution,
            rightHandSide[column] ) )
      {
        return false;
      }
      for( std::size_t otherColumn = 0;
           otherColumn < size;
           ++otherColumn )
      {
        real64 otherValue = 0.0;
        real64 contribution = 0.0;
        if( !projectedGaussSeidelContact::safeDivide(
              jacobian[row * size + otherColumn],
              commonScale,
              otherValue ) ||
            !projectedGaussSeidelContact::safeMultiply(
              value,
              otherValue,
              contribution ) ||
            !projectedGaussSeidelContact::safeAdd(
              normalMatrix[column * size + otherColumn],
              contribution,
              normalMatrix[column * size + otherColumn] ) )
        {
          return false;
        }
      }
    }
  }

  real64 diagonalScale = 0.0;
  for( std::size_t row = 0; row < size; ++row )
  {
    diagonalScale = std::max( diagonalScale,
                              std::abs( normalMatrix[row * size + row] ) );
  }
  if( isZero( diagonalScale ) ||
      !projectedGaussSeidelContact::finiteValue( diagonalScale ) )
  {
    return false;
  }
  for( std::size_t row = 0; row < size; ++row )
  {
    real64 diagonalRegularization = 0.0;
    if( !projectedGaussSeidelContact::safeMultiply(
          regularization,
          diagonalScale,
          diagonalRegularization ) ||
        !projectedGaussSeidelContact::safeAdd(
          normalMatrix[row * size + row],
          diagonalRegularization,
          normalMatrix[row * size + row] ) )
    {
      return false;
    }
  }
  return solveDenseLinearSystem( normalMatrix,
                                 rightHandSide,
                                 1.0e-14,
                                 step );
}

/**
 * @brief Solve all active constraints at one node with damped semismooth Newton.
 *
 * The projection residual is differentiable inside each active-set region and
 * semismooth at contact/friction transitions. A forward-difference generalized
 * Jacobian is used. A regularized least-squares direction handles redundant
 * constraints, while backtracking prevents a Newton step from increasing the
 * residual.
 */
inline NewtonRaphsonContactResult solve(
  std::vector< real64 > const & mass,
  std::vector< Vector3 > & velocity,
  std::vector< Constraint > & constraints,
  integer const maximumIterations,
  real64 const velocityTolerance,
  real64 const finiteDifferenceRelativeStep,
  real64 const minimumLineSearchScale,
  real64 const regularization )
{
  projectedGaussSeidelContact::FloatingPointEnvironmentGuard floatingPointEnvironmentGuard;
  NewtonRaphsonContactResult result;
  if( constraints.empty() )
  {
    return result;
  }

  std::vector< Vector3 > const initialVelocity = velocity;
  auto rejectNode = [&]()
  {
    velocity = initialVelocity;
    for( Constraint & constraint : constraints )
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
      !projectedGaussSeidelContact::finiteValue( velocityTolerance ) ||
      velocityTolerance < 0.0 ||
      !projectedGaussSeidelContact::finiteValue(
        finiteDifferenceRelativeStep ) ||
      !( finiteDifferenceRelativeStep > 0.0 ) ||
      !projectedGaussSeidelContact::finiteValue(
        minimumLineSearchScale ) ||
      !( minimumLineSearchScale > 0.0 ) ||
      minimumLineSearchScale > 1.0 ||
      !projectedGaussSeidelContact::finiteValue( regularization ) ||
      regularization < 0.0 ||
      constraints.size() > std::numeric_limits< std::size_t >::max() / 3 )
  {
    rejectNode();
    return result;
  }

  for( Vector3 const & fieldVelocity : velocity )
  {
    if( !projectedGaussSeidelContact::finiteVector( fieldVelocity ) )
    {
      rejectNode();
      return result;
    }
  }

  for( Constraint & constraint : constraints )
  {
    if( constraint.fieldA < 0 || constraint.fieldB < 0 ||
        constraint.fieldA == constraint.fieldB ||
        static_cast< std::size_t >( constraint.fieldA ) >= mass.size() ||
        static_cast< std::size_t >( constraint.fieldB ) >= mass.size() ||
        !projectedGaussSeidelContact::finiteValue(
          constraint.targetNormalVelocity ) ||
        !projectedGaussSeidelContact::finiteValue(
          constraint.frictionCoefficient ) ||
        constraint.frictionCoefficient < 0.0 ||
        !projectedGaussSeidelContact::finiteValue(
          constraint.normalBiasImpulse ) ||
        constraint.normalBiasImpulse < 0.0 )
    {
      rejectNode();
      return result;
    }

    real64 pairEffectiveMass = 0.0;
    real64 normalNorm = 0.0;
    if( !projectedGaussSeidelContact::effectiveMass(
          mass[constraint.fieldA],
          mass[constraint.fieldB],
          pairEffectiveMass ) ||
        !projectedGaussSeidelContact::robustNorm(
          constraint.normal,
          normalNorm ) ||
        !( normalNorm > 1.0e-20 ) )
    {
      rejectNode();
      return result;
    }
    for( localIndex component = 0; component < 3; ++component )
    {
      if( !projectedGaussSeidelContact::safeDivide(
            constraint.normal[component],
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

  std::size_t const systemSize = 3 * constraints.size();
  if( systemSize > 0 &&
      systemSize > std::numeric_limits< std::size_t >::max() / systemSize )
  {
    rejectNode();
    return result;
  }
  std::vector< real64 > impulse( systemSize, 0.0 );
  std::vector< real64 > residual;
  if( !evaluateResidual( mass,
                         initialVelocity,
                         constraints,
                         impulse,
                         residual,
                         velocity ) ||
      !velocityResidualNorm( mass,
                             constraints,
                             residual,
                             result.residual ) )
  {
    rejectNode();
    return result;
  }
  if( result.residual <= velocityTolerance )
  {
    return result;
  }

  std::vector< real64 > impulseScale( systemSize, 0.0 );
  for( std::size_t constraintIndex = 0;
       constraintIndex < constraints.size();
       ++constraintIndex )
  {
    Constraint const & constraint = constraints[constraintIndex];
    real64 pairEffectiveMass = 0.0;
    Vector3 initialRelative = {{ 0.0, 0.0, 0.0 }};
    real64 initialRelativeNorm = 0.0;
    if( !projectedGaussSeidelContact::effectiveMass(
          mass[constraint.fieldA],
          mass[constraint.fieldB],
          pairEffectiveMass ) ||
        !projectedGaussSeidelContact::tryRelativeVelocity(
          initialVelocity,
          constraint,
          initialRelative ) ||
        !projectedGaussSeidelContact::robustNorm(
          initialRelative,
          initialRelativeNorm ) )
    {
      rejectNode();
      return result;
    }
    real64 const velocityScale = std::max(
      velocityTolerance,
      std::max( initialRelativeNorm,
                std::abs( constraint.targetNormalVelocity ) ) );
    for( localIndex component = 0; component < 3; ++component )
    {
      if( !projectedGaussSeidelContact::safeMultiply(
            pairEffectiveMass,
            velocityScale,
            impulseScale[3 * constraintIndex + component] ) ||
          !( impulseScale[3 * constraintIndex + component] > 0.0 ) )
      {
        rejectNode();
        return result;
      }
    }
  }

  result.converged = false;
  std::vector< real64 > jacobian( systemSize * systemSize, 0.0 );
  std::vector< real64 > perturbedImpulse = impulse;
  std::vector< real64 > perturbedResidual;
  std::vector< Vector3 > perturbedVelocity;
  std::vector< real64 > rightHandSide( systemSize, 0.0 );
  std::vector< real64 > step;
  std::vector< real64 > candidateImpulse( systemSize, 0.0 );
  std::vector< real64 > candidateResidual;
  std::vector< Vector3 > candidateVelocity;

  for( integer iteration = 0; iteration < maximumIterations; ++iteration )
  {
    bool jacobianSafe = true;
    for( std::size_t column = 0; column < systemSize; ++column )
    {
      perturbedImpulse = impulse;
      real64 perturbation = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            finiteDifferenceRelativeStep,
            std::max( std::abs( impulse[column] ), impulseScale[column] ),
            perturbation ) ||
          !( perturbation > 0.0 ) ||
          !projectedGaussSeidelContact::safeAdd(
            perturbedImpulse[column],
            perturbation,
            perturbedImpulse[column] ) ||
          !evaluateResidual( mass,
                             initialVelocity,
                             constraints,
                             perturbedImpulse,
                             perturbedResidual,
                             perturbedVelocity ) )
      {
        jacobianSafe = false;
        break;
      }
      for( std::size_t row = 0; row < systemSize; ++row )
      {
        real64 difference = 0.0;
        if( !projectedGaussSeidelContact::safeSubtract(
              perturbedResidual[row],
              residual[row],
              difference ) ||
            !projectedGaussSeidelContact::safeDivide(
              difference,
              perturbation,
              jacobian[row * systemSize + column] ) )
        {
          jacobianSafe = false;
          break;
        }
      }
      if( !jacobianSafe )
      {
        break;
      }
    }
    if( !jacobianSafe )
    {
      rejectNode();
      return result;
    }

    for( std::size_t row = 0; row < systemSize; ++row )
    {
      rightHandSide[row] = -residual[row];
    }
    bool const directStepAvailable = solveDenseLinearSystem(
      jacobian,
      rightHandSide,
      1.0e-10,
      step );

    real64 currentMerit = 0.0;
    if( !meritNorm( residual, currentMerit ) )
    {
      rejectNode();
      return result;
    }
    bool accepted = false;
    for( integer directionAttempt = 0;
         directionAttempt < 2 && !accepted;
         ++directionAttempt )
    {
      bool const regularizedDirection = directionAttempt == 1;
      if( !regularizedDirection && !directStepAvailable )
      {
        continue;
      }
      if( regularizedDirection )
      {
        if( !regularizedNewtonStep( jacobian,
                                    residual,
                                    regularization,
                                    step ) )
        {
          continue;
        }
      }

      real64 lineSearchScale = 1.0;
      while( true )
      {
        bool candidateSafe = true;
        for( std::size_t entry = 0; entry < systemSize; ++entry )
        {
          real64 scaledStep = 0.0;
          if( !projectedGaussSeidelContact::safeMultiply(
                lineSearchScale,
                step[entry],
                scaledStep ) ||
              !projectedGaussSeidelContact::safeAdd(
                impulse[entry],
                scaledStep,
                candidateImpulse[entry] ) )
          {
            candidateSafe = false;
            break;
          }
        }
        real64 candidateMerit = 0.0;
        real64 candidateVelocityResidual = 0.0;
        candidateSafe = candidateSafe &&
          evaluateResidual( mass,
                            initialVelocity,
                            constraints,
                            candidateImpulse,
                            candidateResidual,
                            candidateVelocity ) &&
          meritNorm( candidateResidual, candidateMerit ) &&
          velocityResidualNorm( mass,
                                constraints,
                                candidateResidual,
                                candidateVelocityResidual );
        if( candidateSafe &&
            ( candidateMerit < currentMerit ||
              candidateVelocityResidual <= velocityTolerance ) )
        {
          impulse = candidateImpulse;
          residual = candidateResidual;
          velocity = candidateVelocity;
          result.residual = candidateVelocityResidual;
          accepted = true;
          if( regularizedDirection )
          {
            ++result.regularizedSteps;
          }
          break;
        }

        if( !candidateSafe )
        {
          ++result.numericalGuardActivations;
        }

        if( lineSearchScale <= minimumLineSearchScale )
        {
          break;
        }
        lineSearchScale = std::max( minimumLineSearchScale,
                                    0.5 * lineSearchScale );
        ++result.lineSearchReductions;
      }
    }

    if( !accepted )
    {
      break;
    }

    result.iterations = iteration + 1;
    if( result.residual <= velocityTolerance )
    {
      result.converged = true;
      break;
    }
  }

  for( std::size_t constraintIndex = 0;
       constraintIndex < constraints.size();
       ++constraintIndex )
  {
    Constraint & constraint = constraints[constraintIndex];
    constraint.accumulatedNormalImpulse = 0.0;
    constraint.accumulatedTangentialImpulse = {{ 0.0, 0.0, 0.0 }};
    if( constraint.bilateral )
    {
      continue;
    }

    std::size_t const offset = 3 * constraintIndex;
    Vector3 const impulseOnA = {{ impulse[offset],
                                  impulse[offset + 1],
                                  impulse[offset + 2] }};
    real64 impulseNormalComponent = 0.0;
    if( !projectedGaussSeidelContact::safeDot(
          impulseOnA,
          constraint.normal,
          impulseNormalComponent ) )
    {
      rejectNode();
      return result;
    }
    constraint.accumulatedNormalImpulse =
      std::max( 0.0, -impulseNormalComponent );
    for( localIndex component = 0; component < 3; ++component )
    {
      real64 normalComponent = 0.0;
      if( !projectedGaussSeidelContact::safeMultiply(
            impulseNormalComponent,
            constraint.normal[component],
            normalComponent ) ||
          !projectedGaussSeidelContact::safeSubtract(
            impulseOnA[component],
            normalComponent,
            constraint.accumulatedTangentialImpulse[component] ) )
      {
        rejectNode();
        return result;
      }
    }
  }

  return result;
}

} // namespace newtonRaphsonContact
} // namespace mpm
} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_NEWTONRAPHSONCONTACT_HPP_
