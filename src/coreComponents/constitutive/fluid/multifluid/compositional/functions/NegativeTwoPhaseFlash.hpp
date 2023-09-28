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
 * @file NegativeTwoPhaseFlash.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_

#include "common/DataTypes.hpp"
#include "CubicEOSPhaseModel.hpp"
#include "RachfordRice.hpp"
#include "KValueInitialization.hpp"

namespace geos
{

namespace constitutive
{

struct NegativeTwoPhaseFlash
{
public:
  /// Max number of components alloweeed in the class for now
  static constexpr integer maxNumComps = 5;
  /// Max number of iterations
  static constexpr integer maxIterations = 200;
  /// Epsilon used in the calculations
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;
  /// Tolerance for checking fugacity ratio convergence
  static constexpr real64 fugacityTolerance = 1.0e-8;

  /**
   * @brief Perform negative two-phase EOS flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] criticalPressure critical pressures
   * @param[in] criticalTemperature critical temperatures
   * @param[in] acentricFactor acentric factors
   * @param[in] binaryInteractionCoefficients binary coefficients (currently not implemented)
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   * @return an indicator of success of the flash
   */
  template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arrayView1d< real64 const > const composition,
                       arrayView1d< real64 const > const criticalPressure,
                       arrayView1d< real64 const > const criticalTemperature,
                       arrayView1d< real64 const > const acentricFactor,
                       real64 const & binaryInteractionCoefficients,
                       real64 & vapourPhaseMoleFraction,
                       arrayView1d< real64 > const liquidComposition,
                       arrayView1d< real64 > const vapourComposition )
  {
    stackArray1d< real64, maxNumComps > logLiquidFugacity( numComps );
    stackArray1d< real64, maxNumComps > logVapourFugacity( numComps );
    stackArray1d< real64, maxNumComps > kVapourLiquid( numComps );
    stackArray1d< real64, maxNumComps > fugacityRatios( numComps );
    stackArray1d< integer, maxNumComps > presentComponentIds( numComps );

    // Initialise compositions to feed composition
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic];
      vapourComposition[ic] = composition[ic];
    }

    // Check for machine-zero feed values
    integer presentCount = 0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      if( epsilon < composition[ic] )
      {
        presentComponentIds[presentCount++] = ic;
      }
    }
    presentComponentIds.resize( presentCount );

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        criticalPressure,
                                                        criticalTemperature,
                                                        acentricFactor,
                                                        kVapourLiquid );

    bool converged = false;
    for( localIndex iterationCount = 0; iterationCount < maxIterations; ++iterationCount )
    {
      // Solve Rachford-Rice Equation
      vapourPhaseMoleFraction = RachfordRice::solve( kVapourLiquid, composition, presentComponentIds );

      // Assign phase compositions
      for( integer const ic : presentComponentIds )
      {
        liquidComposition[ic] = composition[ic] / ( 1.0 + vapourPhaseMoleFraction * ( kVapourLiquid[ic] - 1.0 ) );
        vapourComposition[ic] = kVapourLiquid[ic] * liquidComposition[ic];
      }

      normalizeComposition( numComps, liquidComposition );
      normalizeComposition( numComps, vapourComposition );

      // Compute the phase fugacities
      CubicEOSPhaseModel< EOS_TYPE_LIQUID >::compute( numComps,
                                                      pressure,
                                                      temperature,
                                                      liquidComposition,
                                                      criticalPressure,
                                                      criticalTemperature,
                                                      acentricFactor,
                                                      binaryInteractionCoefficients,
                                                      logLiquidFugacity );
      CubicEOSPhaseModel< EOS_TYPE_VAPOUR >::compute( numComps,
                                                      pressure,
                                                      temperature,
                                                      vapourComposition,
                                                      criticalPressure,
                                                      criticalTemperature,
                                                      acentricFactor,
                                                      binaryInteractionCoefficients,
                                                      logVapourFugacity );

      // Compute fugacity ratios and check convergence
      converged = true;

      for( integer const ic : presentComponentIds )
      {
        fugacityRatios[ic] = exp( logLiquidFugacity[ic] - logVapourFugacity[ic] ) * liquidComposition[ic] / vapourComposition[ic];
        if( fugacityTolerance < fabs( fugacityRatios[ic] - 1.0 ) )
        {
          converged = false;
        }
      }

      if( converged )
      {
        break;
      }

      // Update K-values
      for( integer const ic : presentComponentIds )
      {
        kVapourLiquid[ic] *= fugacityRatios[ic];
      }
    }

    // Retrieve physical bounds from negative flash values
    if( vapourPhaseMoleFraction <= 0.0 )
    {
      vapourPhaseMoleFraction = 0.0;
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidComposition[ic] = composition[ic];
      }
    }
    else if( 1.0 <= vapourPhaseMoleFraction )
    {
      vapourPhaseMoleFraction = 1.0;
      for( integer ic = 0; ic < numComps; ++ic )
      {
        vapourComposition[ic] = composition[ic];
      }
    }

    return converged;
  }

private:
  /**
   * @brief Normalise a composition in place to ensure that the components add up to unity
   * @param[in] numComps number of components
   * @param[in/out] composition composition to be normalized
   * @return the sum of the given values
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 normalizeComposition( integer const numComps,
                                      arraySlice1d< real64 > const composition )
  {
    real64 totalMoles = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      totalMoles += composition[ic];
    }
    real64 const oneOverTotalMoles = 1.0 / (totalMoles + epsilon);
    for( integer ic = 0; ic < numComps; ++ic )
    {
      composition[ic] *= oneOverTotalMoles;
    }
    return totalMoles;
  }
};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
