/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TestFluidUtilities.hpp
 */

#ifndef GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_
#define GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_

#include "codingUtilities/UnitTestUtilities.hpp"

namespace geos
{

namespace testing
{

namespace internal
{

static constexpr real64 relTol = 1.0e-4;
static constexpr real64 absTol = 1.0e-7;

/**
 * @brief Tests a function against a derivative
 * @details Will calculate the left-sided and the right-sided numerical derivatives of a function
 *          and compare this against a analytically calculated value provided.
 * @param x The value at which the function should be evaluated
 * @param dx The value to use to perturb @c x in the calculation of the numerical derivatives
 * @param derivative The value of the analytically calculated derivative to use for comparison
 * @param function The function which is being tested. This should be a function that returns
          a single real value.
 * @param absTolerance The absolute tolerance to use for the comparison
 * @param relTolerance The relative tolerance to use for the comparison
 */
template< typename FUNCTION >
void testNumericalDerivative( real64 const x, real64 const dx, real64 const derivative,
                              FUNCTION && function,
                              real64 const absTolerance = absTol, real64 const relTolerance = relTol )
{
  real64 const leftValue = function( x-dx );
  real64 const centreValue = function( x );
  real64 const rightValue = function( x+dx );
  real64 const leftDerivative = (centreValue - leftValue) / dx;
  real64 const rightDerivative = (rightValue - centreValue) / dx;
  checkRelativeError( derivative, leftDerivative, relTolerance, absTolerance, "Left derivative" );
  checkRelativeError( derivative, rightDerivative, relTolerance, absTolerance, "Right derivative" );
}

}// namespace internal

}// namespace testing

}// namespace geos

#endif //GEOS_CONSTITUTIVE_UNITTESTS_TESTFLUIDUTILITIES_HPP_
