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
 * @file NegativeTwoPhaseFlash_impl.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_

#include "RachfordRice.hpp"
#include "KValueInitialization.hpp"
#include "FugacityCalculator.hpp"
#include "Utilities.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/FlashParameters.hpp"

namespace geos
{

namespace constitutive
{
namespace compositional
{

namespace negative_two_phase_flash_internal
{

/**
 * @brief Run a short successive-substitution flash using K-values initialized at @p wilsonPressure,
 *        while evaluating all fugacities at @p fugacityPressure (the true cell pressure).
 * @return Equilibrium residual norm after the last (or converged) iteration.
 */
template< int USD2 >
GEOS_HOST_DEVICE
real64 probeEquilibriumResidual( integer const numComps,
                                 real64 const wilsonPressure,
                                 real64 const fugacityPressure,
                                 real64 const temperature,
                                 arraySlice1d< real64 const > const & composition,
                                 ComponentProperties::KernelWrapper const & componentProperties,
                                 FlashData const & flashData,
                                 arraySlice1d< integer const > const & presentComponents,
                                 real64 const flashTolerance,
                                 integer const maxProbeIterations,
                                 arraySlice1d< real64 > const & logLiquidFugacity,
                                 arraySlice1d< real64 > const & logVapourFugacity,
                                 arraySlice1d< real64 > const & fugacityRatios,
                                 arraySlice1d< real64, USD2 > const & liquidComposition,
                                 arraySlice1d< real64, USD2 > const & vapourComposition )
{
  constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  StackArray< real64, 1, maxNumComps > kProbe( numComps );
  if( flashData.liquidEos == EquationOfStateType::SoreideWhitson )
  {
    KValueInitialization::computeSoreideWhitsonKvalue( numComps,
                                                       wilsonPressure,
                                                       temperature,
                                                       componentProperties,
                                                       composition,
                                                       kProbe.toSlice() );
  }
  else
  {
    KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                        wilsonPressure,
                                                        temperature,
                                                        componentProperties,
                                                        kProbe.toSlice() );
  }
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic];
    vapourComposition[ic] = composition[ic];
  }

  real64 error = LvArray::NumericLimits< real64 >::max;
  for( localIndex iterationCount = 0; iterationCount < maxProbeIterations; ++iterationCount )
  {
    real64 const vapourFrac = RachfordRice::solve( kProbe.toSliceConst(), composition, presentComponents );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic] / ( 1.0 + vapourFrac * ( kProbe[ic] - 1.0 ) );
      vapourComposition[ic] = kProbe[ic] * liquidComposition[ic];
    }

    normalizeComposition( numComps, liquidComposition );
    normalizeComposition( numComps, vapourComposition );

    FugacityCalculator::computeLogFugacity( numComps,
                                            fugacityPressure,
                                            temperature,
                                            liquidComposition.toSliceConst(),
                                            componentProperties,
                                            flashData.liquidEos,
                                            flashData,
                                            logLiquidFugacity );
    FugacityCalculator::computeLogFugacity( numComps,
                                            fugacityPressure,
                                            temperature,
                                            vapourComposition.toSliceConst(),
                                            componentProperties,
                                            flashData.vapourEos,
                                            flashData,
                                            logVapourFugacity );

    error = 0.0;
    for( integer const ic : presentComponents )
    {
      fugacityRatios[ic] = ( logLiquidFugacity[ic] - logVapourFugacity[ic] ) + log( liquidComposition[ic] ) - log( vapourComposition[ic] );
      error += (fugacityRatios[ic]*fugacityRatios[ic]);
    }
    error = LvArray::math::sqrt( error );

    if( error < flashTolerance )
    {
      break;
    }

    for( integer ic = 0; ic < numComps; ++ic )
    {
      kProbe[ic] *= exp( fugacityRatios[ic] );
    }
  }
  return error;
}

} // namespace negative_two_phase_flash_internal

template< int USD1, int USD2 >
GEOS_HOST_DEVICE
bool NegativeTwoPhaseFlash::compute( integer const numComps,
                                     real64 const pressure,
                                     real64 const temperature,
                                     arraySlice1d< real64 const > const & composition,
                                     ComponentProperties::KernelWrapper const & componentProperties,
                                     FlashData const & flashData,
                                     arraySlice1d< real64 const > const & continuousFlashParameters,
                                     arraySlice1d< integer const > const & discreteFlashParameters,
                                     arraySlice2d< real64, USD1 > const & kValues,
                                     real64 & vapourPhaseMoleFraction,
                                     arraySlice1d< real64, USD2 > const & liquidComposition,
                                     arraySlice1d< real64, USD2 > const & vapourComposition )
{
  constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  StackArray< real64, 2, 3*maxNumComps > workSpace( 3, maxNumComps );
  arraySlice1d< real64 > logLiquidFugacity = workSpace[0];
  arraySlice1d< real64 > logVapourFugacity = workSpace[1];
  arraySlice1d< real64 > fugacityRatios = workSpace[2];
  stackArray1d< integer, maxNumComps > componentIndices( numComps );
  auto const & kVapourLiquid = kValues[0];

  calculatePresentComponents( numComps, composition, componentIndices );

  auto const presentComponents = componentIndices.toSliceConst();

  // Extract flash parameters
  integer const maxIterations = discreteFlashParameters[FlashParameters::FLASH_MAX_ITERATIONS];
  real64 const flashTolerance = continuousFlashParameters[FlashParameters::FLASH_TOLERANCE];

  // Initialise compositions to feed composition
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic];
    vapourComposition[ic] = composition[ic];
  }

  bool converged = false;
  real64 lastInnerError = 0.0;
  real64 errAfterWilsonAtPressure = 0.0;
  bool haveErrAfterWilsonAtPressure = false;

  integer constexpr numRecoveryPasses = 7;
  integer const maxProbeIterations = LvArray::math::min( maxIterations, integer{ 25 } );

  for( integer recoveryPass = 0; recoveryPass < numRecoveryPasses && !converged; ++recoveryPass )
  {
    if( recoveryPass > 0 )
    {
      real64 wilsonReferencePressure = pressure;

      if( recoveryPass == 6 )
      {
        // Secant / Newton-style update of the Wilson reference pressure using finite differences.
        // Fugacity evaluations remain at the true cell pressure.
        real64 const relPerturbation = 5.0e-4;
        real64 const wPlus = pressure * ( 1.0 + relPerturbation );
        real64 const wMinus = pressure * ( 1.0 - relPerturbation );

        real64 const errPlus = negative_two_phase_flash_internal::probeEquilibriumResidual(
          numComps, wPlus, pressure, temperature, composition, componentProperties, flashData,
          presentComponents, flashTolerance, maxProbeIterations,
          logLiquidFugacity, logVapourFugacity, fugacityRatios,
          liquidComposition, vapourComposition );

        real64 const errMinus = negative_two_phase_flash_internal::probeEquilibriumResidual(
          numComps, wMinus, pressure, temperature, composition, componentProperties, flashData,
          presentComponents, flashTolerance, maxProbeIterations,
          logLiquidFugacity, logVapourFugacity, fugacityRatios,
          liquidComposition, vapourComposition );

        real64 const denom = errPlus - errMinus;
        wilsonReferencePressure = pressure;
        if( haveErrAfterWilsonAtPressure &&
            LvArray::math::abs( denom ) > flashTolerance * flashTolerance + pressure * LvArray::NumericLimits< real64 >::epsilon )
        {
          wilsonReferencePressure = pressure - errAfterWilsonAtPressure * ( wPlus - wMinus ) / denom;
        }
        wilsonReferencePressure = LvArray::math::min( LvArray::math::max( wilsonReferencePressure, 0.85 * pressure ), 1.15 * pressure );
      }
      else if( recoveryPass == 2 )
      {
        wilsonReferencePressure = pressure * 0.999;
      }
      else if( recoveryPass == 3 )
      {
        wilsonReferencePressure = pressure * 1.001;
      }
      else if( recoveryPass == 4 )
      {
        wilsonReferencePressure = pressure * 0.99;
      }
      else if( recoveryPass == 5 )
      {
        wilsonReferencePressure = pressure * 1.01;
      }

      if( flashData.liquidEos == EquationOfStateType::SoreideWhitson )
      {
        KValueInitialization::computeSoreideWhitsonKvalue( numComps,
                                                           wilsonReferencePressure,
                                                           temperature,
                                                           componentProperties,
                                                           composition,
                                                           kVapourLiquid );
      }
      else
      {
        KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                            wilsonReferencePressure,
                                                            temperature,
                                                            componentProperties,
                                                            kVapourLiquid );
      }

      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidComposition[ic] = composition[ic];
        vapourComposition[ic] = composition[ic];
      }
    }

    converged = false;
    for( localIndex iterationCount = 0; iterationCount < maxIterations; ++iterationCount )
    {
      // Solve Rachford-Rice Equation
      vapourPhaseMoleFraction = RachfordRice::solve( kVapourLiquid.toSliceConst(), composition, presentComponents );

      // Assign phase compositions
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidComposition[ic] = composition[ic] / ( 1.0 + vapourPhaseMoleFraction * ( kVapourLiquid[ic] - 1.0 ) );
        vapourComposition[ic] = kVapourLiquid[ic] * liquidComposition[ic];
      }

      normalizeComposition( numComps, liquidComposition );
      normalizeComposition( numComps, vapourComposition );

      FugacityCalculator::computeLogFugacity( numComps,
                                              pressure,
                                              temperature,
                                              liquidComposition.toSliceConst(),
                                              componentProperties,
                                              flashData.liquidEos,
                                              flashData,
                                              logLiquidFugacity );
      FugacityCalculator::computeLogFugacity( numComps,
                                              pressure,
                                              temperature,
                                              vapourComposition.toSliceConst(),
                                              componentProperties,
                                              flashData.vapourEos,
                                              flashData,
                                              logVapourFugacity );

      // Compute fugacity ratios and calculate the error
      real64 error = 0.0;
      for( integer const ic : presentComponents )
      {
        fugacityRatios[ic] = ( logLiquidFugacity[ic] - logVapourFugacity[ic] ) + log( liquidComposition[ic] ) - log( vapourComposition[ic] );
        error += (fugacityRatios[ic]*fugacityRatios[ic]);
      }
      error = LvArray::math::sqrt( error );
      lastInnerError = error;

      // Compute fugacity ratios and check convergence
      converged = (error < flashTolerance);

      if( converged )
      {
        break;
      }

      // Update K-values
      for( integer ic = 0; ic < numComps; ++ic )
      {
        kVapourLiquid[ic] *= exp( fugacityRatios[ic] );
      }
    }

    if( recoveryPass == 1 && !converged )
    {
      errAfterWilsonAtPressure = lastInnerError;
      haveErrAfterWilsonAtPressure = true;
    }
  }

  // Test if we have converged to a null or trivial solution
  bool const testNegativeFlash = truncateCompositions( numComps,
                                                       composition,
                                                       vapourPhaseMoleFraction,
                                                       liquidComposition.toSliceConst(),
                                                       vapourComposition.toSliceConst());
  if( testNegativeFlash )
  {
    vapourPhaseMoleFraction = LvArray::math::min( 1.0, LvArray::math::max( 0.0, vapourPhaseMoleFraction ));
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic];
      vapourComposition[ic] = composition[ic];
    }
  }

  return converged;
}

template< integer USD1, integer USD2, integer USD3 >
GEOS_HOST_DEVICE
void NegativeTwoPhaseFlash::computeDerivatives(
  integer const numComps,
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const > const & composition,
  ComponentProperties::KernelWrapper const & componentProperties,
  FlashData const & flashData,
  real64 const & vapourFraction,
  arraySlice1d< real64 const, USD1 > const & liquidComposition,
  arraySlice1d< real64 const, USD1 > const & vapourComposition,
  arraySlice1d< real64, USD2 > const & vapourFractionDerivs,
  arraySlice2d< real64, USD3 > const & liquidCompositionDerivs,
  arraySlice2d< real64, USD3 > const & vapourCompositionDerivs )
{
  constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  constexpr integer maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;

  integer const numDofs = numComps + 2;

  LvArray::forValuesInSlice( vapourFractionDerivs, setZero );
  LvArray::forValuesInSlice( liquidCompositionDerivs, setZero );
  LvArray::forValuesInSlice( vapourCompositionDerivs, setZero );

  // Check for a trivial solution
  real64 diffXY = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const dxy = liquidComposition[ic] - vapourComposition[ic];
    diffXY += (dxy*dxy);
  }
  diffXY = LvArray::math::sqrt( diffXY );

  if( diffXY <  MultiFluidConstants::minForSpeciesPresence )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
      vapourCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
    }
  }
  else
  {
    // Calculate the liquid and vapour fugacities and derivatives
    StackArray< real64, 2, 2*maxNumComps > logFugacity( 2, numComps );
    StackArray< real64, 3, 2*maxNumComps * maxNumDofs > logFugacityDerivs( 2, numComps, numDofs );
    arraySlice1d< real64 > const logLiquidFugacity = logFugacity[0];
    arraySlice1d< real64 > const logVapourFugacity = logFugacity[1];
    arraySlice2d< real64 > const logLiquidFugacityDerivs = logFugacityDerivs[0];
    arraySlice2d< real64 > const logVapourFugacityDerivs = logFugacityDerivs[1];

    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       flashData.liquidEos,
                                                       flashData,
                                                       logLiquidFugacity,
                                                       logLiquidFugacityDerivs );
    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       flashData.vapourEos,
                                                       flashData,
                                                       logVapourFugacity,
                                                       logVapourFugacityDerivs );

    arraySlice1d< real64 > const & liquidFugacity = logLiquidFugacity;
    arraySlice1d< real64 > const & vapourFugacity = logVapourFugacity;
    for( integer i = 0; i < numComps; ++i )
    {
      liquidFugacity[i] = LvArray::math::exp( logLiquidFugacity[i] );
      vapourFugacity[i] = LvArray::math::exp( logVapourFugacity[i] );
    }

    if( vapourFraction < 0.5 )
    {
      computeDerivatives( numComps,
                          composition,
                          vapourFraction,
                          vapourComposition,
                          liquidComposition,
                          vapourFugacity,
                          liquidFugacity,
                          logVapourFugacityDerivs,
                          logLiquidFugacityDerivs,
                          vapourFractionDerivs,
                          vapourCompositionDerivs,
                          liquidCompositionDerivs );
    }
    else
    {
      computeDerivatives( numComps,
                          composition,
                          1.0 - vapourFraction,
                          liquidComposition,
                          vapourComposition,
                          liquidFugacity,
                          vapourFugacity,
                          logLiquidFugacityDerivs,
                          logVapourFugacityDerivs,
                          vapourFractionDerivs,
                          liquidCompositionDerivs,
                          vapourCompositionDerivs );
      LvArray::forValuesInSlice( vapourFractionDerivs, []( real64 & v ){ v *= -1.0; } );
    }
  }
}

template< integer USD1, integer USD2, integer USD3 >
GEOS_HOST_DEVICE
void NegativeTwoPhaseFlash::computeDerivatives(
  integer const numComps,
  arraySlice1d< real64 const > const & totalComposition,
  real64 const & phase1Fraction,
  arraySlice1d< real64 const, USD1 > const & phase1Composition,
  arraySlice1d< real64 const, USD1 > const & phase2Composition,
  arraySlice1d< real64 const > const & phase1Fugacity,
  arraySlice1d< real64 const > const & phase2Fugacity,
  arraySlice2d< real64 const > const & phase1LogFugacityDerivs,
  arraySlice2d< real64 const > const & phase2LogFugacityDerivs,
  arraySlice1d< real64, USD2 > const & phase1FractionDerivs,
  arraySlice2d< real64, USD3 > const & phase1CompositionDerivs,
  arraySlice2d< real64, USD3 > const & phase2CompositionDerivs )
{
  constexpr integer maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;
  constexpr integer maxNumVals = MultiFluidConstants::MAX_NUM_COMPONENTS+1;

  integer const numDofs = numComps + 2;

  StackArray< real64, 2, maxNumVals * maxNumVals, MatrixLayout::COL_MAJOR_PERM > A( numComps + 1, numComps + 1 );
  StackArray< real64, 2, maxNumVals * maxNumDofs, MatrixLayout::COL_MAJOR_PERM > X( numComps + 1, numDofs );

  real64 const VL = phase1Fraction;
  real64 const factor1 = VL / (1.0 - VL);
  real64 const factor2 = 1.0 / (1.0 - VL);
  for( integer i = 0; i < numComps; ++i )
  {
    real64 const phi_2_i = phase2Fugacity[i];
    real64 const phi_1_i = phase1Fugacity[i];

    real64 const xi = phase2Composition[i];
    real64 const yi = phase1Composition[i];

    real64 col_N = 0.0;
    for( integer j = 0; j < numComps; ++j )
    {
      real64 const xj = phase2Composition[j];
      real64 const yj = phase1Composition[j];

      real64 const c1_ij = xi * phi_2_i * phase2LogFugacityDerivs( i, Deriv::dC+j );
      real64 const c2_ij = -yi * phi_1_i * phase1LogFugacityDerivs( i, Deriv::dC+j );

      A( i, j ) = c2_ij - factor1 * c1_ij;

      col_N += c1_ij * (yj - xj);
    }

    A( i, i ) -= (phi_1_i + factor1 * phi_2_i);

    col_N += phi_2_i * (yi - xi);

    A( i, numComps ) = -factor2 * col_N;
    A( numComps, i ) = factor2;
  }
  A( numComps, numComps ) = 0.0;

  // Pressure and temperature derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const phi2 = -phase2Composition[ic] * phase2Fugacity[ic];
    real64 const phi1 = phase1Composition[ic] * phase1Fugacity[ ic ];
    X( ic, Deriv::dP ) = phi2 * phase2LogFugacityDerivs( ic, Deriv::dP ) +
                         phi1 * phase1LogFugacityDerivs( ic, Deriv::dP );
    X( ic, Deriv::dT ) = phi2 * phase2LogFugacityDerivs( ic, Deriv::dT ) +
                         phi1 * phase1LogFugacityDerivs( ic, Deriv::dT );
  }
  X( numComps, Deriv::dP ) = 0.0;
  X( numComps, Deriv::dT ) = 0.0;

  // Composition derivatives
  for( integer k = 0; k < numComps; ++k )
  {
    integer const idof = Deriv::dC + k;

    for( integer i = 0; i < numComps; ++i )
    {
      real64 const phi_2_i = phase2Fugacity[i];
      real64 const xi = phase2Composition[i];
      real64 const zi = totalComposition[i];

      real64 row_i = 0.0;
      for( integer j = 0; j < numComps; ++j )
      {
        row_i -= xi * phi_2_i * phase2LogFugacityDerivs( i, Deriv::dC + j ) * totalComposition[j];
      }
      row_i += xi * phi_2_i * phase2LogFugacityDerivs( i, Deriv::dC + k );
      row_i -= phi_2_i * zi;

      X( i, idof ) = -factor2 * row_i;
    }
    X( k, idof ) -= factor2 * phase2Fugacity[k];
    X( numComps, idof ) = 0.0;
  }

  // Solve linear system
  solveLinearSystem( A.toSlice(), X.toSlice() );

  // Fill in the derivatives
  for( integer i = 0; i < numComps; ++i )
  {
    real64 const xi = phase2Composition[i];
    real64 const yi = phase1Composition[i];
    real64 const zi = totalComposition[i];

    for( integer const idof : {Deriv::dP, Deriv::dT} )
    {
      real64 const dVL = X( numComps, idof );
      real64 const dy = X( i, idof );
      real64 const dx = factor2*((xi - yi)*dVL - VL*dy);

      phase1FractionDerivs( idof ) = dVL;
      phase1CompositionDerivs( i, idof ) = dy;
      phase2CompositionDerivs( i, idof ) = dx;
    }
    for( integer j = 0; j < numComps; ++j )
    {
      integer const idof = Deriv::dC+j;
      real64 const dVL = X( numComps, idof );
      real64 const dy = X( i, idof );
      real64 const dx = factor2*(-zi + (xi - yi)*dVL - VL*dy);

      phase1FractionDerivs( idof ) = dVL;
      phase1CompositionDerivs( i, idof ) = dy;
      phase2CompositionDerivs( i, idof ) = dx;
    }
    phase2CompositionDerivs( i, Deriv::dC+i ) += factor2;
  }
}

template< integer USD1, integer USD2 >
GEOS_FORCE_INLINE
GEOS_HOST_DEVICE
bool NegativeTwoPhaseFlash::truncateCompositions( integer const numComps,
                                                  arraySlice1d< real64 const, USD1 > const & totalComposition,
                                                  real64 const & vapourPhaseMoleFraction,
                                                  arraySlice1d< real64 const, USD2 > const & liquidComposition,
                                                  arraySlice1d< real64 const, USD2 > const & vapourComposition )
{
  GEOS_UNUSED_VAR( numComps );
  GEOS_UNUSED_VAR( totalComposition );
  GEOS_UNUSED_VAR( liquidComposition );
  GEOS_UNUSED_VAR( vapourComposition );
  real64 const V = vapourPhaseMoleFraction;
  real64 const L = 1.0 - vapourPhaseMoleFraction;
  return (V < MultiFluidConstants::epsilon|| L < MultiFluidConstants::epsilon);
}

} // namespace compositional
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_
