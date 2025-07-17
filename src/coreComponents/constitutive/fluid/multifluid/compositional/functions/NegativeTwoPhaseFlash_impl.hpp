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
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

namespace geos
{

namespace constitutive
{
namespace compositional
{

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
  for( localIndex iterationCount = 0; iterationCount < maxIterations; ++iterationCount )
  {
    real64 const error = computeFugacityRatio( numComps,
                                               pressure,
                                               temperature,
                                               composition,
                                               componentProperties,
                                               flashData,
                                               kVapourLiquid.toSliceConst(),
                                               presentComponents,
                                               vapourPhaseMoleFraction,
                                               liquidComposition,
                                               vapourComposition,
                                               logLiquidFugacity,
                                               logVapourFugacity,
                                               fugacityRatios );

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

  // Retrieve physical bounds from negative flash values
  if( vapourPhaseMoleFraction < MultiFluidConstants::epsilon )
  {
    vapourPhaseMoleFraction = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic];
      vapourComposition[ic] = composition[ic];
    }
  }
  else if( 1.0 - vapourPhaseMoleFraction < MultiFluidConstants::epsilon )
  {
    vapourPhaseMoleFraction = 1.0;
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
      vapourCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
    }
  }
  else if( 1.0 - vapourFraction < MultiFluidConstants::epsilon )
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
    StackArray< real64, 1, maxNumComps > logLiquidFugacity( numComps );
    StackArray< real64, 1, maxNumComps > logVapourFugacity( numComps );
    StackArray< real64, 2, maxNumComps * maxNumDofs, MatrixLayout::ROW_MAJOR_PERM > logLiquidFugacityDerivs( numComps, numDofs );
    StackArray< real64, 2, maxNumComps * maxNumDofs, MatrixLayout::ROW_MAJOR_PERM > logVapourFugacityDerivs( numComps, numDofs );

    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       flashData.liquidEos,
                                                       flashData,
                                                       logLiquidFugacity.toSlice(),
                                                       logLiquidFugacityDerivs.toSlice() );
    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       flashData.vapourEos,
                                                       flashData,
                                                       logVapourFugacity.toSlice(),
                                                       logVapourFugacityDerivs.toSlice() );

    constexpr integer maxNumVals = 2*MultiFluidConstants::MAX_NUM_COMPONENTS+1;
    integer const numVals = 2*numComps;
    StackArray< real64, 2, maxNumVals * maxNumVals, MatrixLayout::COL_MAJOR_PERM > A( numVals + 1, numVals + 1 );
    StackArray< real64, 2, maxNumVals * maxNumVals, MatrixLayout::COL_MAJOR_PERM > X( numVals + 1, numVals + 1 );

    LvArray::forValuesInSlice( A.toSlice(), setZero );
    LvArray::forValuesInSlice( X.toSlice(), setZero );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      integer const xi = ic;
      integer const yi = ic + numComps;
      integer const vi = numVals;

      integer e = ic;
      A( e, xi ) = 1.0 - vapourFraction;
      A( e, yi ) = vapourFraction;
      A( e, vi ) = vapourComposition[ic] - liquidComposition[ic];

      e = ic + numComps;
      real64 const phiL = exp( logLiquidFugacity( ic ) );
      real64 const phiV = exp( logVapourFugacity( ic ) );
      for( integer jc = 0; jc < numComps; ++jc )
      {
        integer const xj = jc;
        integer const yj = jc + numComps;
        real64 const dPhiLdx = logLiquidFugacityDerivs( ic, Deriv::dC+jc );
        real64 const dPhiVdy = logVapourFugacityDerivs( ic, Deriv::dC+jc );
        A( e, xj ) =  liquidComposition[ic] * phiL * dPhiLdx;
        A( e, yj ) = -vapourComposition[ic] * phiV * dPhiVdy;
      }
      A( e, xi ) += phiL;
      A( e, yi ) -= phiV;

      e = numVals;
      A( e, xi ) = -1.0;
      A( e, yi ) =  1.0;
    }

    // Pressure and temperature derivatives
    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const phiL = -liquidComposition[ic] * exp( logLiquidFugacity( ic ) );
      real64 const phiV = vapourComposition[ic] * exp( logVapourFugacity( ic ) );
      X( ic + numComps, Deriv::dP ) = phiL * logLiquidFugacityDerivs( ic, Deriv::dP ) +
                                      phiV * logVapourFugacityDerivs( ic, Deriv::dP );
      X( ic + numComps, Deriv::dT ) = phiL * logLiquidFugacityDerivs( ic, Deriv::dT ) +
                                      phiV * logVapourFugacityDerivs( ic, Deriv::dT );
    }

    // Composition derivatives
    for( integer kc = 0; kc < numComps; ++kc )
    {
      integer const idof = Deriv::dC + kc;

      for( integer ic = 0; ic < numComps; ++ic )
      {
        X( ic, idof ) = -composition[ic];
      }
      X( kc, idof ) += 1.0;
    }

    // Solve linear system
    solveLinearSystem( A.toSlice(), X.toSlice() );

    // Fill in the derivatives
    for( integer idof = 0; idof < numDofs; ++idof )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidCompositionDerivs( ic, idof ) = X( ic, idof );
        vapourCompositionDerivs( ic, idof ) = X( ic + numComps, idof );
      }
      vapourFractionDerivs( idof ) = X( numVals, idof );
    }
  }
}

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
real64 NegativeTwoPhaseFlash::computeFugacityRatio(
  integer const numComps,
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const > const & composition,
  ComponentProperties::KernelWrapper const & componentProperties,
  FlashData const & flashData,
  arraySlice1d< real64 const, USD1 > const & kValues,
  arraySlice1d< integer const > const & presentComponents,
  real64 & vapourPhaseMoleFraction,
  arraySlice1d< real64, USD2 > const & liquidComposition,
  arraySlice1d< real64, USD2 > const & vapourComposition,
  arraySlice1d< real64 > const & logLiquidFugacity,
  arraySlice1d< real64 > const & logVapourFugacity,
  arraySlice1d< real64 > const & fugacityRatios )
{
  // Solve Rachford-Rice Equation
  vapourPhaseMoleFraction = RachfordRice::solve( kValues, composition, presentComponents );

  // Assign phase compositions
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic] / ( 1.0 + vapourPhaseMoleFraction * ( kValues[ic] - 1.0 ) );
    vapourComposition[ic] = kValues[ic] * liquidComposition[ic];
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
  return LvArray::math::sqrt( error );
}

template< int USD >
GEOS_HOST_DEVICE
bool NegativeTwoPhaseFlash::solveLinearSystem( arraySlice2d< real64, USD > const & A,
                                               arraySlice2d< real64, USD > const & X )
{
#if defined(GEOS_DEVICE_COMPILE)
  GEOS_UNUSED_VAR( A );
  GEOS_UNUSED_VAR( X );
  return false;
#else
  BlasLapackLA::solveLinearSystem( A, X );
  return true;
#endif
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_
