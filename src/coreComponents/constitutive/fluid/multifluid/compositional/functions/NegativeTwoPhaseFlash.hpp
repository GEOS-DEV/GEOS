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
#include "RachfordRice.hpp"
#include "KValueInitialization.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ComponentProperties.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct NegativeTwoPhaseFlash
{
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Perform negative two-phase EOS flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
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
                       arraySlice1d< real64 const > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       real64 & vapourPhaseMoleFraction,
                       arraySlice1d< real64 > const & liquidComposition,
                       arraySlice1d< real64 > const & vapourComposition )
  {
    constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
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

    calculatePresentComponents( numComps, composition, presentComponentIds );

    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        pressure,
                                                        temperature,
                                                        componentProperties,
                                                        kVapourLiquid );

    bool converged = false;
    for( localIndex iterationCount = 0; iterationCount < MultiFluidConstants::maxSSIIterations; ++iterationCount )
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
      EOS_TYPE_LIQUID::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       logLiquidFugacity );
      EOS_TYPE_VAPOUR::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       logVapourFugacity );

      // Compute fugacity ratios and check convergence
      converged = true;

      for( integer const ic : presentComponentIds )
      {
        fugacityRatios[ic] = exp( logLiquidFugacity[ic] - logVapourFugacity[ic] ) * liquidComposition[ic] / vapourComposition[ic];
        if( MultiFluidConstants::fugacityTolerance < LvArray::math::abs( fugacityRatios[ic] - 1.0 ) )
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

  /**
   * @brief Calculate derivatives from the two-phase negative flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] vapourFraction the calculated vapour (gas) mole fraction
   * @param[in] liquidComposition the calculated liquid phase composition
   * @param[in] vapourComposition the calculated vapour phase composition
   * @param[out] vapourFractionDerivs derivatives of the calculated vapour (gas) mole fraction
   * @param[out] liquidCompositionDerivs derivatives of the calculated liquid phase composition
   * @param[out] vapourCompositionDerivs derivatives of the calculated vapour phase composition
   */
  template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
  GEOS_HOST_DEVICE
  static void computeDerivatives( integer const numComps,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  real64 const vapourFraction,
                                  arraySlice1d< real64 const > const & liquidComposition,
                                  arraySlice1d< real64 const > const & vapourComposition,
                                  arraySlice1d< real64 > const & vapourFractionDerivs,
                                  arraySlice2d< real64 > const & liquidCompositionDerivs,
                                  arraySlice2d< real64 > const & vapourCompositionDerivs )
  {
    constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
    constexpr integer maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;

    integer const numDofs = numComps + 2;

    auto const setZero = []( real64 & val ) { val = 0.0; };
    LvArray::forValuesInSlice( vapourFractionDerivs, setZero );
    LvArray::forValuesInSlice( liquidCompositionDerivs, setZero );
    LvArray::forValuesInSlice( vapourCompositionDerivs, setZero );

    // Check if we are single or 2-phase
    if( vapourFraction < MultiFluidConstants::epsilon )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
      }
    }
    if( 1.0 - vapourFraction < MultiFluidConstants::epsilon )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        vapourCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
      }
    }
    else
    {
      stackArray1d< integer, maxNumComps > presentComponents( numComps );
      calculatePresentComponents( numComps, composition, presentComponents );

      // Calculate the liquid and vapour fugacities and derivatives
      stackArray1d< real64, maxNumComps > logLiquidFugacity( numComps );
      stackArray1d< real64, maxNumComps > logVapourFugacity( numComps );
      stackArray2d< real64, maxNumComps * maxNumDofs > logLiquidFugacityDerivs( numComps, numDofs );
      stackArray2d< real64, maxNumComps * maxNumDofs > logVapourFugacityDerivs( numComps, numDofs );
      EOS_TYPE_LIQUID::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       logLiquidFugacity );
      EOS_TYPE_VAPOUR::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       logVapourFugacity );
      EOS_TYPE_LIQUID::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       logLiquidFugacity,
                                                       logLiquidFugacityDerivs );

      EOS_TYPE_VAPOUR::computeLogFugacityCoefficients( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       logVapourFugacity,
                                                       logVapourFugacityDerivs );

      // Calculate kValues and derivatives
      stackArray1d< real64, maxNumComps > kValues( numComps );
      stackArray2d< real64, maxNumComps * maxNumDofs > kValueDerivs( numComps, numDofs );
      kValues.zero();
      kValueDerivs.zero();
      for( integer const ic : presentComponents )
      {
        kValues[ic] = vapourComposition[ic] / liquidComposition[ic];
        std::cout << "K******* " << exp( logLiquidFugacity[ic] - logVapourFugacity[ic] ) << " " << kValues[ic] << std::endl;
      }
      for( integer const kc : {Deriv::dP, Deriv::dT} )
      {
        for( integer const ic : presentComponents )
        {
          kValueDerivs( ic, kc ) = (logLiquidFugacityDerivs( ic, kc )-logVapourFugacityDerivs( ic, kc )) * kValues[ic];
          std::cout << "<<<****>>> " << kValueDerivs( ic, kc ) << std::endl;
        }
      }
      // Calculate vapour fraction derivatives
      RachfordRice::computeDerivatives( kValues,
                                        kValueDerivs,
                                        composition,
                                        presentComponents,
                                        vapourFraction,
                                        vapourFractionDerivs );
      std::cout << "<<<****>>> dP = " << vapourFractionDerivs( Deriv::dP ) << std::endl;
      std::cout << "<<<****>>> dT = " << vapourFractionDerivs( Deriv::dT ) << std::endl;
    }


    real64 displacedVapourFraction = -1.0;
    stackArray1d< real64, maxNumComps > displacedLiquidComposition( numComps );
    stackArray1d< real64, maxNumComps > displacedVapourComposition( numComps );

    // Pressure derivatives
    real64 const dp = 1.0e-4 * pressure;
    compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >( numComps,
                                                 pressure + dp,
                                                 temperature,
                                                 composition,
                                                 componentProperties,
                                                 displacedVapourFraction,
                                                 displacedLiquidComposition,
                                                 displacedVapourComposition );

    //vapourFractionDerivs[Deriv::dP] = (displacedVapourFraction - vapourFraction) / dp;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidCompositionDerivs( ic, Deriv::dP ) = (displacedLiquidComposition[ic] - liquidComposition[ic]) / dp;
      vapourCompositionDerivs( ic, Deriv::dP ) = (displacedVapourComposition[ic] - vapourComposition[ic]) / dp;
    }

    // Temperature derivatives
    real64 const dT = 1.0e-6 * temperature;
    compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >( numComps,
                                                 pressure,
                                                 temperature + dT,
                                                 composition,
                                                 componentProperties,
                                                 displacedVapourFraction,
                                                 displacedLiquidComposition,
                                                 displacedVapourComposition );

    //vapourFractionDerivs[Deriv::dT] = (displacedVapourFraction - vapourFraction) / dT;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidCompositionDerivs( ic, Deriv::dT ) = (displacedLiquidComposition[ic] - liquidComposition[ic]) / dT;
      vapourCompositionDerivs( ic, Deriv::dT ) = (displacedVapourComposition[ic] - vapourComposition[ic]) / dT;
    }

    // Composition derivatives
    real64 constexpr dz = 1.0e-7;
    stackArray1d< real64, maxNumComps > displacedComposition( numComps );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      displacedComposition[ic] = composition[ic];
    }

    for( integer jc = 0; jc < numComps; ++jc )
    {
      displacedComposition[jc] += dz;
      compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >( numComps,
                                                   pressure,
                                                   temperature,
                                                   displacedComposition,
                                                   componentProperties,
                                                   displacedVapourFraction,
                                                   displacedLiquidComposition,
                                                   displacedVapourComposition );
      displacedComposition[jc] = composition[jc];

      integer const kc = Deriv::dC + jc;
      vapourFractionDerivs[kc] = (displacedVapourFraction - vapourFraction) / dz;
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidCompositionDerivs( ic, kc ) = (displacedLiquidComposition[ic] - liquidComposition[ic]) / dz;
        vapourCompositionDerivs( ic, kc ) = (displacedVapourComposition[ic] - vapourComposition[ic]) / dz;
      }
    }
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
  template< typename ARRAY >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static integer calculatePresentComponents( integer const numComps,
                                             arraySlice1d< real64 const > const & composition,
                                             ARRAY & presentComponents )
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
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 normalizeComposition( integer const numComps,
                                      arraySlice1d< real64 > const & composition )
  {
    real64 totalMoles = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      totalMoles += composition[ic];
    }
    real64 const oneOverTotalMoles = 1.0 / (totalMoles + MultiFluidConstants::epsilon);
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

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
