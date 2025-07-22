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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

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

    constexpr integer maxNumVals = MultiFluidConstants::MAX_NUM_COMPONENTS+1;
    StackArray< real64, 2, maxNumVals * maxNumVals, MatrixLayout::COL_MAJOR_PERM > A( numComps + 1, numComps + 1 );
    StackArray< real64, 2, maxNumVals * maxNumDofs, MatrixLayout::COL_MAJOR_PERM > X( numComps + 1, numDofs );

    arraySlice1d< real64 > const & liquidFugacity = logLiquidFugacity;
    arraySlice1d< real64 > const & vapourFugacity = logVapourFugacity;
    for( integer i = 0; i < numComps; ++i )
    {
      liquidFugacity[i] = LvArray::math::exp( logLiquidFugacity[i] );
      vapourFugacity[i] = LvArray::math::exp( logVapourFugacity[i] );
    }

    // Precalculate some factors
    real64 const V = vapourFraction;
    real64 const factor1 = V / (1.0 - V);
    real64 const factor2 = 1.0 / (1.0 - V);

    for( integer i = 0; i < numComps; ++i )
    {
      real64 const phi_L_i = liquidFugacity[i];
      real64 const phi_V_i = vapourFugacity[i];

      real64 const xi = liquidComposition[i];
      real64 const yi = vapourComposition[i];

      real64 col_N = 0.0;
      for( integer j = 0; j < numComps; ++j )
      {
        real64 const xj = liquidComposition[j];
        real64 const yj = vapourComposition[j];

        real64 const c1_ij = xi * phi_L_i * logLiquidFugacityDerivs( i, Deriv::dC+j );
        real64 const c2_ij = -yi * phi_V_i * logVapourFugacityDerivs( i, Deriv::dC+j );

        A( i, j ) = c2_ij - factor1 * c1_ij;

        col_N += c1_ij * (yj - xj);
      }

      A( i, i ) -= (phi_V_i + factor1 * phi_L_i);

      col_N += phi_L_i * (yi - xi);

      A( i, numComps ) = -factor2 * col_N;
      A( numComps, i ) = factor2;
    }
    A( numComps, numComps ) = 0.0;

    // Pressure and temperature derivatives
    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const phiL = -liquidComposition[ic] * liquidFugacity[ic];
      real64 const phiV = vapourComposition[ic] * vapourFugacity[ ic ];
      X( ic, Deriv::dP ) = phiL * logLiquidFugacityDerivs( ic, Deriv::dP ) +
                           phiV * logVapourFugacityDerivs( ic, Deriv::dP );
      X( ic, Deriv::dT ) = phiL * logLiquidFugacityDerivs( ic, Deriv::dT ) +
                           phiV * logVapourFugacityDerivs( ic, Deriv::dT );
    }
    X( numComps, Deriv::dP ) = 0.0;
    X( numComps, Deriv::dT ) = 0.0;

    // Composition derivatives
    for( integer k = 0; k < numComps; ++k )
    {
      integer const idof = Deriv::dC + k;

      for( integer i = 0; i < numComps; ++i )
      {
        real64 const phi_L_i = liquidFugacity[i];
        real64 const xi = liquidComposition[i];
        real64 const zi = composition[i];

        real64 row_i = 0.0;
        for( integer j = 0; j < numComps; ++j )
        {
          row_i -= xi * phi_L_i * logLiquidFugacityDerivs( i, Deriv::dC + j ) * composition[j];
        }
        row_i += xi * phi_L_i * logLiquidFugacityDerivs( i, Deriv::dC + k );
        row_i -= phi_L_i * zi;

        X( i, idof ) = -factor2 * row_i;
      }
      X( k, idof ) -= factor2 * liquidFugacity[k];
      X( numComps, idof ) = 0.0;
    }

    // Solve linear system
    solveLinearSystem( A.toSlice(), X.toSlice() );

    // Fill in the derivatives
    for( integer i = 0; i < numComps; ++i )
    {
      real64 const xi = liquidComposition[i];
      real64 const yi = vapourComposition[i];
      real64 const zi = composition[i];

      for( integer const idof : {Deriv::dP, Deriv::dT} )
      {
        real64 const dV = X( numComps, idof );
        real64 const dy = X( i, idof );
        real64 const dx = factor2*((xi - yi)*dV - V*dy);

        vapourFractionDerivs( idof ) = dV;
        vapourCompositionDerivs( i, idof ) = dy;
        liquidCompositionDerivs( i, idof ) = dx;
      }
      for( integer j = 0; j < numComps; ++j )
      {
        integer const idof = Deriv::dC+j;
        real64 const dV = X( numComps, idof );
        real64 const dy = X( i, idof );
        real64 const dx = factor2*(-zi + (xi - yi)*dV - V*dy);

        vapourFractionDerivs( idof ) = dV;
        vapourCompositionDerivs( i, idof ) = dy;
        liquidCompositionDerivs( i, idof ) = dx;
      }
      liquidCompositionDerivs( i, Deriv::dC+i ) += factor2;
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

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_
