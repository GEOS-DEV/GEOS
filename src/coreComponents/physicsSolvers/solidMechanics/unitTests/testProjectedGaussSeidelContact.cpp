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

#include "physicsSolvers/solidMechanics/ProjectedGaussSeidelContact.hpp"
#include "physicsSolvers/solidMechanics/NewtonRaphsonContact.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

namespace geos
{
namespace
{

using Constraint = mpm::ProjectedGaussSeidelContactConstraint;
using Velocity = std::array< real64, 3 >;

Constraint normalConstraint( localIndex const fieldA,
                             localIndex const fieldB,
                             Velocity const & normal,
                             real64 const frictionCoefficient = 0.0 )
{
  Constraint constraint;
  constraint.fieldA = fieldA;
  constraint.fieldB = fieldB;
  constraint.normal = normal;
  constraint.frictionCoefficient = frictionCoefficient;
  return constraint;
}

real64 momentum( std::vector< real64 > const & mass,
                 std::vector< Velocity > const & velocity,
                 localIndex const component )
{
  real64 value = 0.0;
  for( std::size_t field = 0; field < mass.size(); ++field )
  {
    value += mass[field] * velocity[field][component];
  }
  return value;
}

TEST( ProjectedGaussSeidelContact, MultifieldLogisticRegressionSelectionIsScoped )
{
  using mpm::projectedGaussSeidelContact::useMultifieldLogisticRegression;

  EXPECT_FALSE( useMultifieldLogisticRegression( 0, 3, false ) );
  EXPECT_FALSE( useMultifieldLogisticRegression( 1, 2, false ) );
  EXPECT_TRUE( useMultifieldLogisticRegression( 1, 3, false ) );
  EXPECT_FALSE( useMultifieldLogisticRegression( 1, 3, true ) );
}

TEST( ProjectedGaussSeidelContact, ContactMassCutoffCombinesGlobalAbsoluteAndRelativeFloors )
{
  using mpm::projectedGaussSeidelContact::contactMassCutoff;

  EXPECT_DOUBLE_EQ( contactMassCutoff( 1.0e-12, 0.0, 0.0, 1.0e-4 ),
                    1.0e-12 );
  EXPECT_DOUBLE_EQ( contactMassCutoff( 1.0e-12, 1.0e-8, 0.0, 1.0e-4 ),
                    1.0e-8 );
  EXPECT_NEAR( contactMassCutoff( 1.0e-12, 1.0e-8, 1.0e-3, 1.0e-4 ),
               1.0e-7,
               1.0e-22 );
  EXPECT_DOUBLE_EQ( contactMassCutoff( 1.0e-6, 1.0e-8, 1.0e-3, 1.0e-4 ),
                    1.0e-6 );

  // Representative masses from the failing Voronoi-pillar node: retain the
  // two resolved fields and exclude the trace-mass third field.
  real64 const voronoiCutoff =
    contactMassCutoff( 1.0e-20, 1.0e-8, 1.0e-3, 9.805e-6 );
  EXPECT_GT( 9.805e-6, voronoiCutoff );
  EXPECT_LE( 4.780e-9, voronoiCutoff );
  EXPECT_GT( 2.357e-6, voronoiCutoff );
}

TEST( ProjectedGaussSeidelContact, ImplicitGapActivationIncludesNumericalZero )
{
  using mpm::projectedGaussSeidelContact::implicitGapConstraintCandidate;

  real64 const tolerance = 1.0e-8;
  EXPECT_TRUE( implicitGapConstraintCandidate( -1.0e-3, tolerance ) );
  EXPECT_TRUE( implicitGapConstraintCandidate( 0.0, tolerance ) );
  EXPECT_TRUE( implicitGapConstraintCandidate( 0.5 * tolerance, tolerance ) );
  EXPECT_TRUE( implicitGapConstraintCandidate( tolerance, tolerance ) );
  EXPECT_FALSE( implicitGapConstraintCandidate( 2.0 * tolerance, tolerance ) );
}

TEST( ProjectedGaussSeidelContact, SeparatingTouchingPairIsNonadhesiveForBothSolvers )
{
  std::vector< real64 > const mass = { 1.0, 1.0 };
  std::vector< Velocity > const initialVelocity = {
    {{ -1.0, 0.0, 0.0 }},
    {{ 1.0, 0.0, 0.0 }}
  };
  Constraint constraint = normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} );
  constraint.gap = 0.0;
  constraint.gapActivationTolerance = 1.0e-8;
  constraint.hasGap = true;

  std::vector< Velocity > pgsVelocity = initialVelocity;
  std::vector< Constraint > pgsConstraints = { constraint };
  mpm::ProjectedGaussSeidelContactResult const pgsResult =
    mpm::projectedGaussSeidelContact::solve(
      mass, pgsVelocity, pgsConstraints, 20, 1.0e-13, 1.0 );

  EXPECT_TRUE( pgsResult.converged );
  EXPECT_DOUBLE_EQ( pgsConstraints[0].accumulatedNormalImpulse, 0.0 );
  EXPECT_EQ( pgsVelocity, initialVelocity );

  std::vector< Velocity > newtonVelocity = initialVelocity;
  std::vector< Constraint > newtonConstraints = { constraint };
  mpm::NewtonRaphsonContactResult const newtonResult =
    mpm::newtonRaphsonContact::solve(
      mass,
      newtonVelocity,
      newtonConstraints,
      20,
      1.0e-13,
      1.0e-6,
      1.0e-4,
      1.0e-12 );

  EXPECT_TRUE( newtonResult.converged );
  EXPECT_DOUBLE_EQ( newtonConstraints[0].accumulatedNormalImpulse, 0.0 );
  EXPECT_EQ( newtonVelocity, initialVelocity );
}

TEST( ProjectedGaussSeidelContact, IsolatedFrictionlessPairMatchesImpulseProjection )
{
  std::vector< real64 > const mass = { 2.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 1.0, 1.0, 0.0 }},
    {{ -2.0, -1.0, 0.0 }}
  };
  real64 const initialMomentum = momentum( mass, velocity, 0 );
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 20, 1.0e-13, 1.0 );

  EXPECT_TRUE( result.converged );
  EXPECT_NEAR( velocity[0][0], 0.0, 1.0e-13 );
  EXPECT_NEAR( velocity[1][0], 0.0, 1.0e-13 );
  EXPECT_NEAR( velocity[0][1], 1.0, 1.0e-13 );
  EXPECT_NEAR( velocity[1][1], -1.0, 1.0e-13 );
  EXPECT_NEAR( momentum( mass, velocity, 0 ), initialMomentum, 1.0e-13 );
}

TEST( ProjectedGaussSeidelContact, CoupledThreeFieldNodeConvergesSimultaneously )
{
  std::vector< real64 > const mass = { 1.0, 1.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 2.0, 0.0, 0.0 }},
    {{ 0.0, 0.0, 0.0 }},
    {{ -2.0, 0.0, 0.0 }}
  };
  real64 const initialMomentum = momentum( mass, velocity, 0 );
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} ),
    normalConstraint( 1, 2, {{ 1.0, 0.0, 0.0 }} )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 100, 1.0e-12, 1.0 );

  EXPECT_TRUE( result.converged );
  EXPECT_GT( result.iterations, 1 );
  EXPECT_GE( velocity[1][0] - velocity[0][0], -1.0e-12 );
  EXPECT_GE( velocity[2][0] - velocity[1][0], -1.0e-12 );
  EXPECT_NEAR( velocity[0][0], 0.0, 2.0e-12 );
  EXPECT_NEAR( velocity[1][0], 0.0, 2.0e-12 );
  EXPECT_NEAR( velocity[2][0], 0.0, 2.0e-12 );
  EXPECT_NEAR( momentum( mass, velocity, 0 ), initialMomentum, 1.0e-13 );
}

TEST( ProjectedGaussSeidelContact, IndependentCornerNormalsRemainSharp )
{
  std::vector< real64 > const mass = { 1.0, 1.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 1.0, 1.0, 0.0 }},
    {{ 0.0, 0.0, 0.0 }},
    {{ 0.0, 0.0, 0.0 }}
  };
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} ),
    normalConstraint( 0, 2, {{ 0.0, 1.0, 0.0 }} )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 50, 1.0e-13, 1.0 );

  EXPECT_TRUE( result.converged );
  EXPECT_GE( velocity[1][0] - velocity[0][0], -1.0e-13 );
  EXPECT_GE( velocity[2][1] - velocity[0][1], -1.0e-13 );
  EXPECT_NEAR( momentum( mass, velocity, 0 ), 1.0, 1.0e-13 );
  EXPECT_NEAR( momentum( mass, velocity, 1 ), 1.0, 1.0e-13 );
}

TEST( ProjectedGaussSeidelContact, BilateralConstraintEnforcesBondedVelocity )
{
  std::vector< real64 > const mass = { 1.0, 3.0 };
  std::vector< Velocity > velocity = {
    {{ 4.0, -2.0, 1.0 }},
    {{ 0.0, 2.0, -1.0 }}
  };
  Constraint constraint = normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} );
  constraint.bilateral = true;
  std::vector< Constraint > constraints = { constraint };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 20, 1.0e-13, 1.0 );

  EXPECT_TRUE( result.converged );
  EXPECT_NEAR( velocity[0][0], 1.0, 1.0e-13 );
  EXPECT_NEAR( velocity[0][1], 1.0, 1.0e-13 );
  EXPECT_NEAR( velocity[0][2], -0.5, 1.0e-13 );
  for( localIndex i = 0; i < 3; ++i )
  {
    EXPECT_NEAR( velocity[1][i], velocity[0][i], 1.0e-13 );
  }
}

TEST( ProjectedGaussSeidelContact, CoulombProjectionLimitsTangentialImpulse )
{
  std::vector< real64 > const mass = { 2.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 1.0, 1.0, 0.0 }},
    {{ -2.0, -1.0, 0.0 }}
  };
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }}, 0.5 )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 50, 1.0e-13, 1.0 );

  EXPECT_TRUE( result.converged );
  real64 const tangentialImpulse =
    mpm::projectedGaussSeidelContact::norm(
      constraints[0].accumulatedTangentialImpulse );
  EXPECT_NEAR( constraints[0].accumulatedNormalImpulse, 2.0, 1.0e-13 );
  EXPECT_NEAR( tangentialImpulse, 1.0, 1.0e-13 );
  EXPECT_LE( tangentialImpulse,
             0.5 * constraints[0].accumulatedNormalImpulse + 1.0e-13 );
}

TEST( ProjectedGaussSeidelContact, CoulombProjectionUsesRotatedPairFrame )
{
  real64 const angle = 0.37;
  Velocity const normal = {{ std::cos( angle ), std::sin( angle ), 0.0 }};
  Velocity const tangent = {{ -std::sin( angle ), std::cos( angle ), 0.0 }};
  std::vector< real64 > const mass = { 1.0, 1.0 };
  std::vector< Velocity > velocity( 2 );
  for( localIndex component = 0; component < 3; ++component )
  {
    velocity[0][component] = normal[component] + tangent[component];
    velocity[1][component] = -normal[component] - tangent[component];
  }
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, normal, 0.3 )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 100, 1.0e-13, 1.0 );

  ASSERT_TRUE( result.converged );
  real64 const normalImpulse = constraints[0].accumulatedNormalImpulse;
  Velocity const & tangentialImpulse = constraints[0].accumulatedTangentialImpulse;
  EXPECT_NEAR( mpm::projectedGaussSeidelContact::dot( tangentialImpulse, normal ),
               0.0,
               1.0e-13 );
  EXPECT_NEAR( mpm::projectedGaussSeidelContact::norm( tangentialImpulse ),
               0.3 * normalImpulse,
               1.0e-13 );
}

TEST( ProjectedGaussSeidelContact, CoupledFrictionConstraintsShareAField )
{
  std::vector< real64 > const mass = { 1.0, 1.0, 2.0 };
  std::vector< Velocity > velocity = {
    {{ 0.10, 0.20, 0.0 }},
    {{ 0.05, -0.20, 0.0 }},
    {{ -0.10, 0.0, 0.0 }}
  };
  Velocity const initialMomentum = {{
    momentum( mass, velocity, 0 ),
    momentum( mass, velocity, 1 ),
    momentum( mass, velocity, 2 )
  }};
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 2, {{ 1.0, 0.0, 0.0 }}, 0.3 ),
    normalConstraint( 1, 2, {{ 1.0, 0.0, 0.0 }}, 0.3 )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 200, 1.0e-12, 1.0 );

  ASSERT_TRUE( result.converged );
  for( Constraint const & constraint : constraints )
  {
    real64 const tangentialImpulse =
      mpm::projectedGaussSeidelContact::norm(
        constraint.accumulatedTangentialImpulse );
    EXPECT_LE( tangentialImpulse,
               constraint.frictionCoefficient * constraint.accumulatedNormalImpulse + 1.0e-12 );
  }
  EXPECT_GE( velocity[2][0] - velocity[0][0], -1.0e-12 );
  EXPECT_GE( velocity[2][0] - velocity[1][0], -1.0e-12 );
  for( localIndex component = 0; component < 3; ++component )
  {
    EXPECT_NEAR( momentum( mass, velocity, component ),
                 initialMomentum[component],
                 1.0e-13 );
  }
}

TEST( NewtonRaphsonContact, IsolatedFrictionlessPairMatchesImpulseProjection )
{
  std::vector< real64 > const mass = { 2.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 1.0, 1.0, 0.0 }},
    {{ -2.0, -1.0, 0.0 }}
  };
  real64 const initialMomentum = momentum( mass, velocity, 0 );
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} )
  };

  mpm::NewtonRaphsonContactResult const result =
    mpm::newtonRaphsonContact::solve(
      mass, velocity, constraints, 50, 1.0e-12, 1.0e-6, 1.0e-4, 1.0e-12 );

  EXPECT_TRUE( result.converged );
  EXPECT_NEAR( velocity[0][0], 0.0, 1.0e-10 );
  EXPECT_NEAR( velocity[1][0], 0.0, 1.0e-10 );
  EXPECT_NEAR( velocity[0][1], 1.0, 1.0e-10 );
  EXPECT_NEAR( velocity[1][1], -1.0, 1.0e-10 );
  EXPECT_NEAR( momentum( mass, velocity, 0 ), initialMomentum, 1.0e-12 );
}

TEST( NewtonRaphsonContact, RedundantThreeFieldConstraintsConverge )
{
  std::vector< real64 > const mass = { 1.0, 1.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 2.0, 0.0, 0.0 }},
    {{ 0.0, 0.0, 0.0 }},
    {{ -2.0, 0.0, 0.0 }}
  };
  real64 const initialMomentum = momentum( mass, velocity, 0 );
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} ),
    normalConstraint( 1, 2, {{ 1.0, 0.0, 0.0 }} ),
    normalConstraint( 0, 2, {{ 1.0, 0.0, 0.0 }} )
  };

  mpm::NewtonRaphsonContactResult const result =
    mpm::newtonRaphsonContact::solve(
      mass, velocity, constraints, 50, 1.0e-12, 1.0e-6, 1.0e-4, 1.0e-12 );

  EXPECT_TRUE( result.converged );
  EXPECT_LE( result.residual, 1.0e-12 );
  EXPECT_NEAR( velocity[0][0], 0.0, 1.0e-10 );
  EXPECT_NEAR( velocity[1][0], 0.0, 1.0e-10 );
  EXPECT_NEAR( velocity[2][0], 0.0, 1.0e-10 );
  EXPECT_NEAR( momentum( mass, velocity, 0 ), initialMomentum, 1.0e-12 );
}

TEST( NewtonRaphsonContact, IndependentCornerNormalsRemainSharp )
{
  std::vector< real64 > const mass = { 1.0, 1.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 1.0, 1.0, 0.0 }},
    {{ 0.0, 0.0, 0.0 }},
    {{ 0.0, 0.0, 0.0 }}
  };
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} ),
    normalConstraint( 0, 2, {{ 0.0, 1.0, 0.0 }} )
  };

  mpm::NewtonRaphsonContactResult const result =
    mpm::newtonRaphsonContact::solve(
      mass, velocity, constraints, 50, 1.0e-12, 1.0e-6, 1.0e-4, 1.0e-12 );

  ASSERT_TRUE( result.converged );
  EXPECT_GE( velocity[1][0] - velocity[0][0], -1.0e-10 );
  EXPECT_GE( velocity[2][1] - velocity[0][1], -1.0e-10 );
  EXPECT_NEAR( momentum( mass, velocity, 0 ), 1.0, 1.0e-12 );
  EXPECT_NEAR( momentum( mass, velocity, 1 ), 1.0, 1.0e-12 );
}

TEST( NewtonRaphsonContact, RegularizationHandlesDuplicateBilateralConstraints )
{
  std::vector< real64 > const mass = { 1.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 4.0, -2.0, 1.0 }},
    {{ 0.0, 2.0, -1.0 }}
  };
  Constraint constraint = normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} );
  constraint.bilateral = true;
  std::vector< Constraint > constraints = { constraint, constraint };

  mpm::NewtonRaphsonContactResult const result =
    mpm::newtonRaphsonContact::solve(
      mass, velocity, constraints, 50, 1.0e-12, 1.0e-6, 1.0e-4, 1.0e-12 );

  ASSERT_TRUE( result.converged );
  EXPECT_GT( result.regularizedSteps, 0 );
  for( localIndex component = 0; component < 3; ++component )
  {
    EXPECT_NEAR( velocity[0][component], velocity[1][component], 1.0e-10 );
  }
}

TEST( NewtonRaphsonContact, CoulombProjectionUsesRotatedPairFrame )
{
  real64 const angle = 0.37;
  Velocity const normal = {{ std::cos( angle ), std::sin( angle ), 0.0 }};
  Velocity const tangent = {{ -std::sin( angle ), std::cos( angle ), 0.0 }};
  std::vector< real64 > const mass = { 1.0, 1.0 };
  std::vector< Velocity > velocity( 2 );
  for( localIndex component = 0; component < 3; ++component )
  {
    velocity[0][component] = normal[component] + tangent[component];
    velocity[1][component] = -normal[component] - tangent[component];
  }
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, normal, 0.3 )
  };

  mpm::NewtonRaphsonContactResult const result =
    mpm::newtonRaphsonContact::solve(
      mass, velocity, constraints, 50, 1.0e-12, 1.0e-6, 1.0e-4, 1.0e-12 );

  ASSERT_TRUE( result.converged );
  real64 const normalImpulse = constraints[0].accumulatedNormalImpulse;
  Velocity const & tangentialImpulse = constraints[0].accumulatedTangentialImpulse;
  EXPECT_NEAR( mpm::projectedGaussSeidelContact::dot( tangentialImpulse, normal ),
               0.0,
               1.0e-10 );
  EXPECT_NEAR( mpm::projectedGaussSeidelContact::norm( tangentialImpulse ),
               0.3 * normalImpulse,
               1.0e-10 );
}

TEST( NewtonRaphsonContact, CoupledCoulombConstraintsShareAField )
{
  std::vector< real64 > const mass = { 1.0, 1.0, 2.0 };
  std::vector< Velocity > velocity = {
    {{ 0.10, 0.20, 0.0 }},
    {{ 0.05, -0.20, 0.0 }},
    {{ -0.10, 0.0, 0.0 }}
  };
  Velocity const initialMomentum = {{
    momentum( mass, velocity, 0 ),
    momentum( mass, velocity, 1 ),
    momentum( mass, velocity, 2 )
  }};
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 2, {{ 1.0, 0.0, 0.0 }}, 0.3 ),
    normalConstraint( 1, 2, {{ 1.0, 0.0, 0.0 }}, 0.3 )
  };

  mpm::NewtonRaphsonContactResult const result =
    mpm::newtonRaphsonContact::solve(
      mass, velocity, constraints, 50, 1.0e-12, 1.0e-6, 1.0e-4, 1.0e-12 );

  ASSERT_TRUE( result.converged );
  for( Constraint const & contactConstraint : constraints )
  {
    real64 const tangentialImpulse =
      mpm::projectedGaussSeidelContact::norm(
        contactConstraint.accumulatedTangentialImpulse );
    EXPECT_LE( tangentialImpulse,
               contactConstraint.frictionCoefficient *
               contactConstraint.accumulatedNormalImpulse + 1.0e-10 );
  }
  for( localIndex component = 0; component < 3; ++component )
  {
    EXPECT_NEAR( momentum( mass, velocity, component ),
                 initialMomentum[component],
                 1.0e-12 );
  }
}

TEST( ProjectedGaussSeidelContact, CohesivePairIsSuppressedUnlessPreventionIsEnabled )
{
  EXPECT_TRUE(
    mpm::projectedGaussSeidelContact::suppressContactForCohesivePair( 1, 1, 0 ) );
  EXPECT_FALSE(
    mpm::projectedGaussSeidelContact::suppressContactForCohesivePair( 1, 1, 1 ) );
  EXPECT_FALSE(
    mpm::projectedGaussSeidelContact::suppressContactForCohesivePair( 1, 0, 0 ) );
}

TEST( ProjectedGaussSeidelContact, SafeMultiplyDoesNotOverflowDuringItsGuardCheck )
{
  real64 const maximum = std::numeric_limits< real64 >::max();
  real64 product = 0.0;

  EXPECT_TRUE( mpm::projectedGaussSeidelContact::safeMultiply(
    maximum, 0.5, product ) );
  EXPECT_DOUBLE_EQ( product, maximum * 0.5 );
  EXPECT_FALSE( mpm::projectedGaussSeidelContact::safeMultiply(
    maximum, 2.0, product ) );
}

TEST( ProjectedGaussSeidelContact, ZeroMassIsRejectedAndVelocityIsRolledBack )
{
  std::vector< real64 > const mass = { 0.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 1.0, 0.0, 0.0 }},
    {{ -1.0, 0.0, 0.0 }}
  };
  std::vector< Velocity > const initialVelocity = velocity;
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 20, 1.0e-12, 1.0 );

  EXPECT_FALSE( result.converged );
  EXPECT_TRUE( result.numericalFailure );
  EXPECT_TRUE( result.usedSafeFallback );
  EXPECT_GT( result.numericalGuardActivations, 0 );
  EXPECT_EQ( velocity, initialVelocity );
}

TEST( ProjectedGaussSeidelContact, OverflowingRelativeVelocityUsesFiniteRollback )
{
  real64 const maximum = std::numeric_limits< real64 >::max();
  std::vector< real64 > const mass = { 1.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ maximum, 0.0, 0.0 }},
    {{ -maximum, 0.0, 0.0 }}
  };
  std::vector< Velocity > const initialVelocity = velocity;
  std::vector< Constraint > constraints = {
    normalConstraint( 0, 1, {{ 1.0, 0.0, 0.0 }} )
  };

  mpm::ProjectedGaussSeidelContactResult const result =
    mpm::projectedGaussSeidelContact::solve(
      mass, velocity, constraints, 20, 1.0e-12, 1.0 );

  EXPECT_FALSE( result.converged );
  EXPECT_TRUE( result.numericalFailure );
  EXPECT_TRUE( result.usedSafeFallback );
  EXPECT_EQ( velocity, initialVelocity );
}

TEST( NewtonRaphsonContact, NonfiniteNormalIsRejectedAndVelocityIsRolledBack )
{
  std::vector< real64 > const mass = { 1.0, 1.0 };
  std::vector< Velocity > velocity = {
    {{ 1.0, 0.0, 0.0 }},
    {{ -1.0, 0.0, 0.0 }}
  };
  std::vector< Velocity > const initialVelocity = velocity;
  std::vector< Constraint > constraints = {
    normalConstraint(
      0,
      1,
      {{ std::numeric_limits< real64 >::infinity(), 0.0, 0.0 }} )
  };

  mpm::NewtonRaphsonContactResult const result =
    mpm::newtonRaphsonContact::solve(
      mass, velocity, constraints, 20, 1.0e-12, 1.0e-6, 1.0e-4, 1.0e-12 );

  EXPECT_FALSE( result.converged );
  EXPECT_TRUE( result.numericalFailure );
  EXPECT_TRUE( result.usedSafeFallback );
  EXPECT_GT( result.numericalGuardActivations, 0 );
  EXPECT_EQ( velocity, initialVelocity );
}

} // namespace
} // namespace geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
