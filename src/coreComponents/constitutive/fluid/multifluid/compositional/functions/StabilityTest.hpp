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
 * @file StabilityTest.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_

#include "KValueInitialization.hpp"
#include "FugacityCalculator.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/FlashData.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct StabilityTest
{
private:
  static constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
public:
  /**
   * @brief Perform a two-phase stability test
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] equationOfState The equation of state
   * @param[in] flashData The parameters required for the flash
   * @param[out] tangentPlaneDistance the minimum tangent plane distance (TPD)
   * @param[out] kValues the k-values estimated from the stationary points
   * @return a flag indicating that 2 stationary points have been found
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, USD1 > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       EquationOfStateType const & equationOfState,
                       FlashData const & flashData,
                       real64 & tangentPlaneDistance,
                       arraySlice1d< real64, USD2 > const & kValues )
  {
    stackArray2d< real64, 4*maxNumComps > workSpace( 4, numComps );
    arraySlice1d< real64 > logFugacity = workSpace[0];
    arraySlice1d< real64 > normalizedComposition = workSpace[1];
    arraySlice1d< real64 > logTrialComposition = workSpace[2];
    arraySlice1d< real64 > hyperplane = workSpace[3]; // h-parameter
    stackArray1d< integer, maxNumComps > availableComponents( numComps );

    calculatePresentComponents( numComps, composition, availableComponents );
    auto const presentComponents = availableComponents.toSliceConst();

    LvArray::forValuesInSlice( workSpace.toSlice(), []( real64 & a ){ a = 0.0; } );

    // Calculate the hyperplane parameter
    // h_i = log( z_i ) + log( phi_i )
    FugacityCalculator::computeLogFugacity( numComps,
                                            pressure,
                                            temperature,
                                            composition,
                                            componentProperties,
                                            equationOfState,
                                            flashData,
                                            logFugacity );
    for( integer const ic : presentComponents )
    {
      hyperplane[ic] = LvArray::math::log( composition[ic] ) + logFugacity[ic];
    }

    if( equationOfState == EquationOfStateType::SoreideWhitson )
    {
      // Initialise the trial compositions using SW uniform values
      KValueInitialization::computeSoreideWhitsonKvalue( numComps,
                                                         pressure,
                                                         temperature,
                                                         componentProperties,
                                                         presentComponents,
                                                         kValues );
    }
    else
    {
      // Initialise the trial compositions using Wilson k-Values
      KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                          pressure,
                                                          temperature,
                                                          componentProperties,
                                                          kValues );
    }

    // The mimimum TPD over all trial compositions
    tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;

    // Maximum error (convergence failure) over all trial compositions
    real64 maxError = -LvArray::NumericLimits< real64 >::max;

    for( real64 const & alpha : { 1.0, -1.0, 3.0, -3.0 } )
    {
      // Initialise next sample
      for( integer const ic : presentComponents )
      {
        logTrialComposition[ic] = LvArray::math::log( composition[ic] ) + alpha*LvArray::math::log( kValues[ic] );
        normalizedComposition[ic] = LvArray::math::exp( logTrialComposition[ic] );
      }

      real64 trialError = LvArray::NumericLimits< real64 >::max;
      for( localIndex iterationCount = 0; iterationCount < MultiFluidConstants::maxSSIIterations; ++iterationCount )
      {
        // Normalise the composition and calculate the fugacity
        real64 const totalMoles = normalizeComposition( numComps, normalizedComposition );
        FugacityCalculator::computeLogFugacity( numComps,
                                                pressure,
                                                temperature,
                                                normalizedComposition.toSliceConst(),
                                                componentProperties,
                                                equationOfState,
                                                flashData,
                                                logFugacity );

        // Calculate the TPD
        real64 tpd = 0.0;
        for( integer const ic : presentComponents )
        {
          tpd += composition[ic] + totalMoles * normalizedComposition[ic] * (logTrialComposition[ic] + logFugacity[ic] - hyperplane[ic] - 1.0);
        }
        if( tpd < tangentPlaneDistance )
        {
          tangentPlaneDistance = tpd;
        }
        if( tangentPlaneDistance < -MultiFluidConstants::fugacityTolerance )
        {
          break;
        }

        // Check stationarity
        real64 error = 0.0;
        for( integer const ic : presentComponents )
        {
          real64 const dG =  logTrialComposition[ic] + logFugacity[ic] - hyperplane[ic];
          error += (dG*dG);
        }
        error = LvArray::math::sqrt( error );
        trialError = LvArray::math::min( trialError, error );
        if( trialError < MultiFluidConstants::fugacityTolerance )
        {
          break;
        }

        // Update to next step
        for( integer const ic : presentComponents )
        {
          logTrialComposition[ic] = hyperplane[ic] - logFugacity[ic];
          normalizedComposition[ic] = LvArray::math::exp( logTrialComposition[ic] );
        }
      }
      maxError = LvArray::math::max( maxError, trialError );
      if( tangentPlaneDistance < -MultiFluidConstants::fugacityTolerance )
      {
        break;
      }
    }
    return ( tangentPlaneDistance < -MultiFluidConstants::fugacityTolerance ) || (maxError < MultiFluidConstants::fugacityTolerance);
  }

private:
  /**
   * @brief Calculate which components are present.
   * @details Creates a list of indices whose components have non-zero mole fraction.
   * @param[in] numComps number of components
   * @param[in] composition the composition of the fluid
   * @param[out] presentComponents the list of present components
   * @return the number of present components
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static integer calculatePresentComponents( integer const numComps,
                                             arraySlice1d< real64 const > const & composition,
                                             stackArray1d< integer, maxNumComps > & presentComponents )
  {
    // Check for machine-zero feed values
    integer presentCount = 0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      if( MultiFluidConstants::epsilon < composition[ic] )
      {
        presentComponents[presentCount++] = ic;
      }
    }
    presentComponents.resize( presentCount );
    return presentCount;
  }

  /**
   * @brief Normalise a composition in place to ensure that the components add up to unity
   * @param[in] numComps number of components
   * @param[in/out] composition composition to be normalized
   * @return the sum of the given values
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 normalizeComposition( integer const numComps,
                                      arraySlice1d< real64, USD > const & composition )
  {
    real64 totalMoles = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      totalMoles += composition[ic];
    }
    GEOS_ASSERT( MultiFluidConstants::epsilon < totalMoles );
    real64 const oneOverTotalMoles = 1.0 / totalMoles;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      composition[ic] *= oneOverTotalMoles;
    }
    return totalMoles;
  }
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
