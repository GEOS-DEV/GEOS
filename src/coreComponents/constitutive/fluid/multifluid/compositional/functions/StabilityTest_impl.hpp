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
 * @file StabilityTest_impl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_IMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_IMPL_HPP_

#include "StabilityTest.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
bool StabilityTest::compute( integer const numComps,
                             real64 const pressure,
                             real64 const temperature,
                             arraySlice1d< real64 const, USD1 > const & composition,
                             ComponentProperties::KernelWrapper const & componentProperties,
                             FlashData const & flashData,
                             arraySlice1d< real64 const, USD2 > const & kValues,
                             arraySlice1d< real64 const > const & continuousFlashParameters,
                             arraySlice1d< integer const > const & discreteFlashParameters,
                             bool & unstableMixture,
                             EquationOfStateType & incipientEquationOfState,
                             arraySlice1d< real64 > const & incipientComposition )
{
  integer constexpr maxDofs = maxNumComps + 2;

  integer const numDofs = 2 + numComps;

  stackArray2d< real64, 5*maxNumComps > workSpace( 5, numComps );
  stackArray2d< real64, maxNumComps *maxNumComps > jacobian( numComps, numComps );
  stackArray2d< real64, maxNumComps *maxDofs > derivatives( numComps, numDofs );
  arraySlice1d< real64 > logFugacity = workSpace[0];
  arraySlice1d< real64 > logTestComposition = workSpace[1];
  arraySlice1d< real64 > testComposition = workSpace[2];
  arraySlice1d< real64 > residual = workSpace[3];
  arraySlice1d< real64 > hyperplane = workSpace[4];   // h-parameter
  stackArray1d< integer, maxNumComps > availableComponents( numComps );
  stackArray1d< integer, maxNumComps > unavailableComponents( numComps );

  calculatePresentComponents( numComps, composition, availableComponents );
  calculateAbsentComponents( numComps, composition, unavailableComponents );
  auto const presentComponents = availableComponents.toSliceConst();
  auto const absentComponents = unavailableComponents.toSliceConst();

  LvArray::forValuesInSlice( workSpace.toSlice(), setZero );

  // Extract flash parameters
  integer const maxIterations = discreteFlashParameters[FlashParameters::STABILITY_MAX_ITERATIONS];
  real64 const stabilityThreshold = continuousFlashParameters[FlashParameters::STABILITY_THRESHOLD];
  real64 const stabilityTolerance = continuousFlashParameters[FlashParameters::STABILITY_TOLERANCE];
  real64 const stabilitySSITolerance = LvArray::math::max( stabilityTolerance, continuousFlashParameters[FlashParameters::SSI_TOLERANCE] );

  integer const numConfigSteps = flashData.liquidEos == flashData.vapourEos ? 1 : 2;

  // Initialise by assuming a stable mixture
  unstableMixture = false;

  // Flag to indicate all trial compositions converged
  bool allConverged = true;

  // Measure of distance to trivial solution
  real64 maxDistanceToTrivialSolution = 0.0;

  for( integer configStep = 0; configStep < numConfigSteps; ++configStep )
  {
    EquationOfStateType const sampleEos = (configStep == 0) ? flashData.liquidEos : flashData.vapourEos;
    EquationOfStateType const incipientEos = sampleEos;

    real64 tpd = 0.0;

    // Calculate the hyperplane parameter
    // h_i = log( z_i ) + log( phi_i )
    LvArray::forValuesInSlice( hyperplane, setZero );
    for( integer const ic : presentComponents )
    {
      logTestComposition[ic] = LvArray::math::log( composition[ic] );
    }
    calculateResidual( numComps,
                       presentComponents,
                       pressure,
                       temperature,
                       componentProperties,
                       sampleEos,
                       flashData,
                       logTestComposition.toSliceConst(),
                       tpd,
                       testComposition,
                       logFugacity,
                       hyperplane );

    for( real64 const alpha : { 1.0, -1.0, 0.0 } )
    {
      if( LvArray::math::abs( alpha ) < MultiFluidConstants::epsilon )
      {
        for( integer const ic : presentComponents )
        {
          logTestComposition[ic] = 1.0;
        }
      }
      else
      {
        for( integer const ic : presentComponents )
        {
          real64 const logK = LvArray::math::log( kValues[ic] );
          logTestComposition[ic] = LvArray::math::log( composition[ic] ) + alpha * logK;
        }
      }

      // Start iterations for this sample
      bool converged = false;

      localIndex iterationCount = 0;

      // Start with SSI iterations
      for(; iterationCount < maxIterations; ++iterationCount )
      {
        for( integer const ic : presentComponents )
        {
          residual[ic] = -hyperplane[ic];
        }
        real64 error = calculateResidual( numComps,
                                          presentComponents,
                                          pressure,
                                          temperature,
                                          componentProperties,
                                          incipientEos,
                                          flashData,
                                          logTestComposition.toSliceConst(),
                                          tpd,
                                          testComposition,
                                          logFugacity,
                                          residual );
        // Update to next step
        for( integer const ic : presentComponents )
        {
          logTestComposition[ic] -= residual[ic];
        }

        // Check stationarity
        if( error < stabilitySSITolerance )
        {
          converged = true;
          break;
        }
      }

      // Start with Newton iterations
      for(; iterationCount < maxIterations; ++iterationCount )
      {
        for( integer const ic : presentComponents )
        {
          residual[ic] = -hyperplane[ic];
        }
        real64 error = calculateResidualAndJacobian( numComps,
                                                     presentComponents,
                                                     pressure,
                                                     temperature,
                                                     componentProperties,
                                                     incipientEos,
                                                     flashData,
                                                     logTestComposition.toSliceConst(),
                                                     tpd,
                                                     testComposition,
                                                     logFugacity,
                                                     derivatives.toSlice(),
                                                     residual,
                                                     jacobian.toSlice() );
        for( integer const ic : absentComponents )
        {
          jacobian(ic, ic) = 1.0;
        }
        solveLinearSystem( jacobian.toSlice(), residual );
        // Update to next step
        for( integer const ic : presentComponents )
        {
          logTestComposition[ic] -= residual[ic];
        }

        // Check stationarity
        if( error < stabilityTolerance )
        {
          converged = true;
          break;
        }
      }
      allConverged = allConverged && converged;

      // Calculate the tangent-plane-distance (TPD) and distance to the trivial solution
      real64 distance = 0.0;
      for( integer const ic : presentComponents )
      {
        real64 const dZ = testComposition[ic] - composition[ic];
        distance += (dZ*dZ);
      }

      if( maxDistanceToTrivialSolution < distance )
      {
        maxDistanceToTrivialSolution = distance;
        incipientEquationOfState = incipientEos;
        for( integer ic = 0; ic < numComps; ++ic )
        {
          incipientComposition[ic] = testComposition[ic];
        }
      }
      // The mixture is unstable if either the TPD is negative or if the TPD is zero but the incipient composition
      // is different from the sample composition
      // stabilityThreshold is negative (default -1e-8)
      if( tpd < stabilityThreshold )
      {
        unstableMixture = true;
      }
      else if((tpd < -stabilityThreshold) && (-stabilityThreshold < distance))
      {
        unstableMixture = true;
      }
      if( unstableMixture )
      {
        return true;
      }
    }
  }

  // The test is successful if either we have an unstable mixture or all test compositions converged to stationarity
  return unstableMixture || allConverged;
}

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
real64 StabilityTest::calculateResidual( integer const numComps,
                                         arraySlice1d< integer const > const & presentComponents,
                                         real64 const pressure,
                                         real64 const temperature,
                                         ComponentProperties::KernelWrapper const & componentProperties,
                                         EquationOfStateType const equationOfState,
                                         FlashData const & flashData,
                                         arraySlice1d< real64 const, USD1 > const & logTestComposition,
                                         real64 & tangentPlaneDistance,
                                         arraySlice1d< real64, USD2 > const & testComposition,
                                         arraySlice1d< real64 > const & testLogFugacity,
                                         arraySlice1d< real64, USD3 > const & residual )
{
  for( integer const ic : presentComponents )
  {
    testComposition[ic] = LvArray::math::exp( logTestComposition[ic] );
  }
  real64 const sumY = normalizeComposition( numComps, testComposition );

  FugacityCalculator::computeLogFugacity( numComps,
                                          pressure,
                                          temperature,
                                          testComposition.toSliceConst(),
                                          componentProperties,
                                          equationOfState,
                                          flashData,
                                          testLogFugacity );

  real64 error = 0.0;
  tangentPlaneDistance = 1.0;
  for( integer const ic : presentComponents )
  {
    residual[ic] += logTestComposition[ic] + testLogFugacity[ic];
    tangentPlaneDistance += sumY * testComposition[ic] * (residual[ic] - 1.0);
    error += (residual[ic]*residual[ic]);
  }
  error = LvArray::math::sqrt( error );
  return error;
}

template< int USD1, int USD2, int USD3, int USD4 >
GEOS_HOST_DEVICE
real64 StabilityTest::calculateResidualAndJacobian( integer const numComps,
                                                    arraySlice1d< integer const > const & presentComponents,
                                                    real64 const pressure,
                                                    real64 const temperature,
                                                    ComponentProperties::KernelWrapper const & componentProperties,
                                                    EquationOfStateType const equationOfState,
                                                    FlashData const & flashData,
                                                    arraySlice1d< real64 const, USD1 > const & logTestComposition,
                                                    real64 & tangentPlaneDistance,
                                                    arraySlice1d< real64, USD2 > const & testComposition,
                                                    arraySlice1d< real64 > const & testLogFugacity,
                                                    arraySlice2d< real64 > const & logTestFugacityDerivs,
                                                    arraySlice1d< real64, USD3 > const & residual,
                                                    arraySlice2d< real64, USD4 > const & jacobian )
{
  for( integer const ic : presentComponents )
  {
    testComposition[ic] = LvArray::math::exp( logTestComposition[ic] );
  }
  real64 const sumY = normalizeComposition( numComps, testComposition );

  FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                     pressure,
                                                     temperature,
                                                     testComposition.toSliceConst(),
                                                     componentProperties,
                                                     equationOfState,
                                                     flashData,
                                                     testLogFugacity,
                                                     logTestFugacityDerivs );

  LvArray::forValuesInSlice( jacobian, setZero );

  tangentPlaneDistance = 1.0;
  real64 error = 0.0;
  for( integer const ic : presentComponents )
  {
    for( integer const jc : presentComponents )
    {
      real64 const yj = testComposition[jc];

      real64 const dphi_ic_dyjc = logTestFugacityDerivs( ic, Deriv::dC + jc );

      real64 sum_dphi_ic_dykc = dphi_ic_dyjc;

      for( integer const kc : presentComponents )
      {
        real64 const yk = testComposition[kc];
        real64 const dphi_ic_dykc = logTestFugacityDerivs( ic, Deriv::dC + kc );
        sum_dphi_ic_dykc -= (dphi_ic_dykc * yk);
      }
      jacobian( ic, jc ) = yj * sum_dphi_ic_dykc;
    }
    jacobian( ic, ic ) += 1.0;

    residual[ic] += logTestComposition[ic] + testLogFugacity[ic];
    tangentPlaneDistance += sumY * testComposition[ic] * (residual[ic] - 1.0);
    error += (residual[ic]*residual[ic]);
  }
  error = LvArray::math::sqrt( error );
  return error;
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_IMPL_HPP_
