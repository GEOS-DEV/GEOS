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
 * @file NegativeTwoPhaseFlash_impl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_

#include "CubicEOSPhaseModel.hpp"
#include "RachfordRice.hpp"
//#include "KValueInitialization.hpp"

namespace geos
{

namespace constitutive
{

struct KValueInitialization
{
  static void
  computeWilsonGasLiquidKvalue( integer const numComps,
                                real64 const GEOS_UNUSED_PARAM( pressure ),
                                real64 const GEOS_UNUSED_PARAM( temperature ),
                                arrayView1d< real64 const > const GEOS_UNUSED_PARAM( criticalPressure ),
                                arrayView1d< real64 const > const GEOS_UNUSED_PARAM( criticalTemperature ),
                                arrayView1d< real64 const > const GEOS_UNUSED_PARAM( acentricFactor ),
                                arraySlice1d< real64 > const kValues )
  {
    for( integer ic = 0; ic <numComps; ++ic )
    {
      kValues[ic] = 1.0;
    }
  }
};

template< typename EOS_TYPE >
bool NegativeTwoPhaseFlash< EOS_TYPE >::compute(
  integer const numComps,
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
  stackArray1d< integer, maxNumComps > presentComponentIds;

  // Initialise compositions to feed composition
  for( integer ic = 0; ic<numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic];
    vapourComposition[ic] = composition[ic];
  }

  // Check for machine-zero feed values
  for( integer ic = 0; ic<numComps; ++ic )
  {
    if( epsilon < composition[ic] )
    {
      presentComponentIds.emplace_back( ic );
    }
  }

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
    CubicEOSPhaseModel< EOS_TYPE >::compute( numComps,
                                             pressure,
                                             temperature,
                                             liquidComposition,
                                             criticalPressure,
                                             criticalTemperature,
                                             acentricFactor,
                                             binaryInteractionCoefficients,
                                             logLiquidFugacity );
    CubicEOSPhaseModel< EOS_TYPE >::compute( numComps,
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
      fugacityRatios[ic]= exp( logLiquidFugacity[ic] - logVapourFugacity[ic] ) * liquidComposition[ic] / vapourComposition[ic];
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
    for( auto ic : presentComponentIds )
    {
      kVapourLiquid[ic] *= fugacityRatios[ic];
    }
  }

  // Retrieve physical bounds from negative flash values
  if( vapourPhaseMoleFraction <= 0.0 )
  {
    vapourPhaseMoleFraction = 0.0;
    for( integer ic = 0; ic<numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic];
    }
  }
  else if( 1.0 <= vapourPhaseMoleFraction )
  {
    vapourPhaseMoleFraction = 1.0;
    for( integer ic = 0; ic<numComps; ++ic )
    {
      vapourComposition[ic] = composition[ic];
    }
  }

  return converged;
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_
