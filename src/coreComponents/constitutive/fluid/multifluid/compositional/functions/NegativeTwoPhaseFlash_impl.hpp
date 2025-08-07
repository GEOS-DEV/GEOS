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
  constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  constexpr integer maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;

  integer const numDofs = numComps + 2;

  StackArray< real64, 2, 3*maxNumComps > workSpace( 3, numComps );
  StackArray< real64, 1, maxNumComps+1 > residualSpace( numComps );
  arraySlice1d< real64 > logLiquidFugacity = workSpace[0];
  arraySlice1d< real64 > logVapourFugacity = workSpace[1];
  arraySlice1d< real64 > residual = residualSpace.toSlice();
  auto const & kVapourLiquid = kValues[0];
  stackArray1d< integer, maxNumComps > availableComponents( numComps );
  stackArray1d< integer, maxNumComps > unavailableComponents( numComps );

  calculatePresentComponents( numComps, composition, availableComponents );
  calculateAbsentComponents( numComps, composition, unavailableComponents );
  auto const presentComponents = availableComponents.toSliceConst();
  auto const absentComponents = unavailableComponents.toSliceConst();

  // Extract flash parameters
  integer const maxIterations = discreteFlashParameters[FlashParameters::FLASH_MAX_ITERATIONS];
  real64 const flashTolerance = continuousFlashParameters[FlashParameters::FLASH_TOLERANCE];
  real64 const flashSSITolerance = LvArray::math::max( flashTolerance, continuousFlashParameters[FlashParameters::SSI_TOLERANCE] );

  // Initialise compositions to feed composition
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic];
    vapourComposition[ic] = composition[ic];
  }

  bool converged = false;
  localIndex iterationCount = 0; 

  // Start SSI iterations
  for( ; ( !converged ) && (iterationCount < maxIterations); ++iterationCount )
  {
    // Solve Rachford-Rice Equation
    vapourPhaseMoleFraction = RachfordRice::solve( kVapourLiquid.toSliceConst(), composition, presentComponents );

     real64 const error = calculateResidual(  numComps,
                        pressure,
                        temperature,
                        composition,
                        componentProperties,
                        flashData,
                        kVapourLiquid.toSliceConst(),
                        vapourPhaseMoleFraction,
                        liquidComposition,
                        vapourComposition,
                       logLiquidFugacity,
                        logVapourFugacity,
                        residual );

    // Update K-values
    for( integer const ic : presentComponents )
    {
      kVapourLiquid[ic] = exp( logLiquidFugacity[ic] - logVapourFugacity[ic] );
    }

    // Check convergence
    converged = ( error < flashTolerance );

    if (error < flashSSITolerance)
    {
      break;
    }
  }

  // Start Newton iterations
    StackArray< real64, 3, 2*maxNumComps * maxNumDofs > logFugacityDerivs( 2, numComps, numDofs );
    StackArray< real64, 2, (maxNumComps+1) * (maxNumComps+1) > jacobian( numComps+1, numComps+1 );
    arraySlice2d< real64 > const logLiquidFugacityDerivs = logFugacityDerivs[0];
    arraySlice2d< real64 > const logVapourFugacityDerivs = logFugacityDerivs[1];

for( ; ( !converged ) && (iterationCount < maxIterations); ++iterationCount )
  {
 real64 const error = calculateResidualAndJacobian( numComps,
                        pressure,
                        temperature,
                        composition,
                        componentProperties,
                        flashData,
                        kVapourLiquid.toSliceConst(),
                        vapourPhaseMoleFraction,
                        liquidComposition,
                        vapourComposition,
                       logLiquidFugacity,
                        logVapourFugacity,
                                            logLiquidFugacityDerivs,
                                            logVapourFugacityDerivs,
                                            residual,
                                            jacobian.toSlice() );

    // Account for missing components in Jacobian
                                            for( integer const ic : absentComponents )
        {
          jacobian(ic, ic) = 1.0;
        }

        // Solve for next step
        solveLinearSystem( jacobian.toSlice(), residual );
            // Update K-values
    for( integer const ic : presentComponents )
    {
      kVapourLiquid[ic] *= exp( -residual[ic] );
    }       
    vapourPhaseMoleFraction -= residual[numComps];       
        // Check convergence
        converged = ( error < flashTolerance );
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

template< int USD >
void calculateResidualAndJacobian_gen( int const numComps,
                                       double const pressure,
                                       double const temperature,
                                       ComponentProperties::KernelWrapper const & componentProperties,
                                       FlashData const & flashData,
                                       std::vector< double > const & totalComposition,
                                       double const V,
                                       std::vector< double > const & kValues,
                                       std::vector< double > & residual,
                                       std::vector< std::vector< double > > & jacobian );
template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
real64 NegativeTwoPhaseFlash::calculateResidual( integer const numComps,
                                                 real64 const pressure,
                                                 real64 const temperature,
                                                 arraySlice1d< real64 const > const & composition,
                                                 ComponentProperties::KernelWrapper const & componentProperties,
                                                 FlashData const & flashData,
                                                 arraySlice1d< real64 const, USD1 > const & kValues,
                                                 real64 const & vapourPhaseMoleFraction,
                                                 arraySlice1d< real64, USD2 > const & liquidComposition,
                                                 arraySlice1d< real64, USD2 > const & vapourComposition,
                                                 arraySlice1d< real64 > const & logLiquidFugacity,
                                                 arraySlice1d< real64 > const & logVapourFugacity,
                                                 arraySlice1d< real64, USD3 > const & residual )
{
  // Calculate phase and normalize compositions
  real64 massBalance = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic] / ( 1.0 + vapourPhaseMoleFraction * ( kValues[ic] - 1.0 ) );
    vapourComposition[ic] = kValues[ic] * liquidComposition[ic];
    massBalance += (vapourComposition[ic] - liquidComposition[ic]);
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
  for( integer ic = 0; ic < numComps; ++ic )
  {
    residual[ic] = liquidComposition[ic]*LvArray::math::exp( logLiquidFugacity[ic] ) - vapourComposition[ic]*LvArray::math::exp( logVapourFugacity[ic] );
    error += (residual[ic]*residual[ic]);
  }
  residual[numComps] = massBalance;
  error += (residual[numComps]*residual[numComps]);
  error = LvArray::math::sqrt( error );

  return error;
}

template< int USD1, int USD2, int USD3, int USD4 >
GEOS_HOST_DEVICE
real64 NegativeTwoPhaseFlash::calculateResidualAndJacobian( integer const numComps,
                                                          real64 const pressure,
                                                          real64 const temperature,
                                                          arraySlice1d< real64 const > const & composition,
                                                          ComponentProperties::KernelWrapper const & componentProperties,
                                                          FlashData const & flashData,
                                                          arraySlice1d< real64 const, USD1 > const & kValues,
                                                          real64 const & vapourPhaseMoleFraction,
                                                          arraySlice1d< real64, USD2 > const & liquidComposition,
                                                          arraySlice1d< real64, USD2 > const & vapourComposition,
                                                          arraySlice1d< real64 > const & logLiquidFugacity,
                                                          arraySlice1d< real64 > const & logVapourFugacity,
                                                          arraySlice2d< real64 > const & logLiquidFugacityDerivs,
                                                          arraySlice2d< real64 > const & logVapourFugacityDerivs,
                                                          arraySlice1d< real64, USD3 > const & residual,
                                                          arraySlice2d< real64, USD4 > const & jacobian )
{
  /*
     GEOS_UNUSED_VAR( numComps );
     GEOS_UNUSED_VAR( pressure );
     GEOS_UNUSED_VAR( temperature );
     GEOS_UNUSED_VAR( composition );
     GEOS_UNUSED_VAR( componentProperties );
     GEOS_UNUSED_VAR( flashData );
     GEOS_UNUSED_VAR( kValues );
     GEOS_UNUSED_VAR( vapourPhaseMoleFraction );
     GEOS_UNUSED_VAR( liquidComposition );
     GEOS_UNUSED_VAR( vapourComposition );
     GEOS_UNUSED_VAR( logLiquidFugacity );
     GEOS_UNUSED_VAR( logVapourFugacity );
     GEOS_UNUSED_VAR( logLiquidFugacityDerivs );
     GEOS_UNUSED_VAR( logVapourFugacityDerivs );
     GEOS_UNUSED_VAR( residual );
     GEOS_UNUSED_VAR( jacobian );

     std::vector< double > compositionVector( numComps );
     std::vector< double > kValuesVector( numComps );
     std::vector< double > residualVector( numComps+1, 0.0 );
     std::vector< std::vector< double > > jacobianVector( numComps+1, residualVector );
     for( integer ic = 0; ic < numComps; ++ic )
     {
     compositionVector[ic] = composition[ic];
     kValuesVector[ic] = kValues[ic];
     }

     calculateResidualAndJacobian_gen< 0 >( numComps,
                                         pressure,
                                         temperature,
                                         componentProperties,
                                         flashData,
                                         compositionVector,
                                         vapourPhaseMoleFraction,
                                         kValuesVector,
                                         residualVector,
                                         jacobianVector );
     for( integer ic = 0; ic <= numComps; ++ic )
     {
     residual[ic] = residualVector[ic];
     for( integer jc = 0; jc <= numComps; ++jc )
     {
      jacobian( ic, jc ) = jacobianVector[ic][jc];
     }
     }

     for( integer ic = 0; ic <= numComps; ++ic )
     {
     residual[ic] = 0.0;
     //   for( integer jc = 0; jc <= numComps; ++jc )
     //   {
     //     jacobian( ic, jc ) = 0.0;
     //   }
     }
   */
  real64 const & V = vapourPhaseMoleFraction;

  // Calculate phase compositions ---
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const denominator = 1.0 + V * (kValues[ic] - 1.0);
    liquidComposition[ic] = composition[ic] / denominator;
    vapourComposition[ic] = kValues[ic] * liquidComposition[ic];
  }

  real64 const sumX = normalizeComposition( numComps, liquidComposition );
  real64 const sumY = normalizeComposition( numComps, vapourComposition );

  // Calculate the fugacity and derivatives
  FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                     pressure,
                                                     temperature,
                                                     liquidComposition.toSliceConst(),
                                                     componentProperties,
                                                     flashData.liquidEos,
                                                     flashData,
                                                     logLiquidFugacity,
                                                     logLiquidFugacityDerivs );
  FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                     pressure,
                                                     temperature,
                                                     vapourComposition.toSliceConst(),
                                                     componentProperties,
                                                     flashData.vapourEos,
                                                     flashData,
                                                     logVapourFugacity,
                                                     logVapourFugacityDerivs );

  // Convert from log (use same memory space: log fugacity no longer required)
  arraySlice1d< real64 > const & liquidFugacity = logLiquidFugacity;
  arraySlice1d< real64 > const & vapourFugacity = logVapourFugacity;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidFugacity[ic] = LvArray::math::exp( logLiquidFugacity[ic] );
    vapourFugacity[ic] = LvArray::math::exp( logVapourFugacity[ic] );
  }

  // Calculate the residual vector
  real64 error = 0.0;
  residual[numComps] = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const zi = composition[ic];
    real64 const xi = liquidComposition[ic];
    real64 const yi = vapourComposition[ic];
    real64 const Ki = kValues[ic];
    residual[ic] = xi * liquidFugacity[ic] - yi * vapourFugacity[ic];
    residual[numComps] += zi * (Ki - 1.0) / (1.0 + V * (Ki - 1.0));
    error += (residual[ic]*residual[ic]);
  }
  error += (residual[numComps]*residual[numComps]);
  error = LvArray::math::sqrt( error );

  real64 const invSumX = 1.0/sumX;
  real64 const invSumY = 1.0/sumY;

  real64 sum_d_xs_d_V = 0.0;
  real64 sum_d_ys_d_V = 0.0;
  for( int ic = 0; ic < numComps; ++ic )
  {
    real64 const zi = composition[ic];
    real64 const Ki = kValues[ic];

    real64 const invKi = 1.0 / ( 1.0 + V * (Ki - 1.0) );
    real64 const d_xsi_d_V = -zi * (Ki - 1.0) * invKi * invKi;
    real64 const d_ysi_d_V = Ki * d_xsi_d_V;

    sum_d_xs_d_V += d_xsi_d_V;
    sum_d_ys_d_V += d_ysi_d_V;
  }

  // Start calculating the Jacobian
  real64 sum_j_N_V = 0.0;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    real64 const zi = composition[ic];
    real64 const xi = liquidComposition[ic];
    real64 const yi = vapourComposition[ic];
    real64 const Ki = kValues[ic];

    real64 const invKi = 1.0 / ( 1.0 + V * (Ki - 1.0) );
    real64 const d_xsi_d_V = -zi * (Ki - 1.0) * invKi * invKi;
    real64 const d_ysi_d_V = Ki * d_xsi_d_V;

    real64 const d_xi_d_V = (d_xsi_d_V - xi * sum_d_xs_d_V) * invSumX;
    real64 const d_yi_d_V = (d_ysi_d_V - yi * sum_d_ys_d_V) * invSumY;

    real64 sum_dXid_xj_dxj_d_V = 0.0;
    real64 sum_dYid_yj_dyj_d_V = 0.0;

    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 const zj = composition[jc];
      real64 const xj = liquidComposition[jc];
      real64 const yj = vapourComposition[jc];
      real64 const Kj = kValues[jc];

      real64 const dphiLi_dxj = logLiquidFugacityDerivs( ic, Deriv::dC + jc );
      real64 const dphiVi_dyj = logVapourFugacityDerivs( ic, Deriv::dC + jc );

      real64 const invKj = 1.0 / ( 1.0 + V * (Kj - 1.0) );
      real64 const d_xsj_d_V = -zj * (Kj - 1.0) * invKj * invKj;
      real64 const d_ysj_d_V = Kj * d_xsj_d_V;

      real64 const d_xsj_d_kj = -zj * V * Kj * invKj * invKj;
      real64 const d_ysj_d_kj = Kj * (xj * sumX + d_xsj_d_kj);

      real64 const d_xj_d_V = (d_xsj_d_V - xj * sum_d_xs_d_V) * invSumX;
      real64 const d_yj_d_V = (d_ysj_d_V - yj * sum_d_ys_d_V) * invSumY;

      sum_dXid_xj_dxj_d_V += dphiLi_dxj * d_xj_d_V;
      sum_dYid_yj_dyj_d_V += dphiVi_dyj * d_yj_d_V;

      real64 const dij = (ic == jc) ? 1.0 : 0.0;
      real64 const d_xi_d_kj = d_xsj_d_kj * (dij - xi) * invSumX;
      real64 const d_yi_d_kj = d_ysj_d_kj * (dij - yi) * invSumY;

      real64 sum_dXid_xk_dxk_d_kj = 0.0;
      real64 sum_dYid_yk_dyk_d_kj = 0.0;
      for( integer kc = 0; kc < numComps; ++kc )
      {
        real64 const djk = (jc == kc) ? 1.0 : 0.0;
        real64 const d_xk_d_kj = d_xsj_d_kj * (djk - liquidComposition[kc]) * invSumX;
        real64 const d_yk_d_kj = d_ysj_d_kj * (djk - vapourComposition[kc]) * invSumY;

        real64 const dphiLi_dxk = logLiquidFugacityDerivs( ic, Deriv::dC + kc );
        real64 const dphiVi_dyk = logVapourFugacityDerivs( ic, Deriv::dC + kc );

        sum_dXid_xk_dxk_d_kj += dphiLi_dxk * d_xk_d_kj;
        sum_dYid_yk_dyk_d_kj += dphiVi_dyk * d_yk_d_kj;
      }
      jacobian( ic, jc ) = liquidFugacity[ic] * (d_xi_d_kj + xi * sum_dXid_xk_dxk_d_kj) -
                           vapourFugacity[ic] * (d_yi_d_kj + yi * sum_dYid_yk_dyk_d_kj);
    }
    jacobian( ic, numComps ) = liquidFugacity[ic] * (d_xi_d_V + xi * sum_dXid_xj_dxj_d_V) -
                               vapourFugacity[ic] * (d_yi_d_V + yi * sum_dYid_yj_dyj_d_V);
    jacobian( numComps, ic ) = zi * Ki * invKi* invKi;

    sum_j_N_V -= zi * (Ki - 1.0) * (Ki - 1.0) * invKi * invKi;
  }
  jacobian( numComps, numComps ) = sum_j_N_V;

  return error;
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

// Assume calculateFugacity is defined elsewhere.
// This function calculates the log fugacities and their derivatives for a given phase.
template< int USD >
void calculateFugacity( int const numComps,
                        double const pressure,
                        double const temperature,
                        ComponentProperties::KernelWrapper const & componentProperties,
                        FlashData const & flashData,
                        EquationOfStateType const equationOfState,
                        std::vector< double > const & phaseComposition,
                        std::vector< double > & logFugacity,
                        std::vector< std::vector< double > > & logFugacityDerivatives )
{
  StackArray< real64, 2, 2*9 > workSpace( 2, numComps );
  arraySlice1d< real64, USD > phaseCompositionSlice = workSpace[0];
  arraySlice1d< real64, USD > logFugacitySlice = workSpace[1];
  StackArray< real64, 2, 9*11 > derivatives( numComps, numComps+2 );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    phaseCompositionSlice[ic] = phaseComposition[ic];
  }
  FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                     pressure,
                                                     temperature,
                                                     phaseCompositionSlice.toSliceConst(),
                                                     componentProperties,
                                                     equationOfState,
                                                     flashData,
                                                     logFugacitySlice,
                                                     derivatives.toSlice() );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    logFugacity[ic] = logFugacitySlice[ic];
    for( integer jc = 0; jc < numComps; ++jc )
    {
      logFugacityDerivatives[ic][jc] = derivatives( ic, 2 + jc );
    }
  }
}

/**
 * @brief Calculates the residual vector and Jacobian matrix for a flash calculation.
 * * @param numComps Number of components (N).
 * @param pressure System pressure.
 * @param temperature System temperature.
 * @param totalComposition The total mole fractions (z_i).
 * @param V The vapor fraction guess.
 * @param kValues The K-values (K_i).
 * @param residual Output vector for the residual F.
 * @param jacobian Output matrix for the Jacobian of F.
 */
template< int USD >
void calculateResidualAndJacobian_gen( int const numComps,
                                       double const pressure,
                                       double const temperature,
                                       ComponentProperties::KernelWrapper const & componentProperties,
                                       FlashData const & flashData,
                                       std::vector< double > const & totalComposition,
                                       double const V,
                                       std::vector< double > const & kValues,
                                       std::vector< double > & residual,
                                       std::vector< std::vector< double > > & jacobian )
{

  // Helper lambda for Kronecker delta
  auto kronecker_delta = []( int i, int j ) { return (i == j) ? 1.0 : 0.0; };

  // --- Part 1: Calculate intermediate values and phase compositions ---
  std::vector< double > x_s( numComps );
  std::vector< double > y_s( numComps );
  double x_sum = 0.0;
  double y_sum = 0.0;

  for( int i = 0; i < numComps; ++i )
  {
    double denominator = 1.0 + V * (kValues[i] - 1.0);
    x_s[i] = totalComposition[i] / denominator;
    y_s[i] = kValues[i] * x_s[i];
    x_sum += x_s[i];
    y_sum += y_s[i];
  }

  std::vector< double > x_i( numComps );
  std::vector< double > y_i( numComps );
  for( int i = 0; i < numComps; ++i )
  {
    x_i[i] = x_s[i] / x_sum;
    y_i[i] = y_s[i] / y_sum;
  }

  // --- Part 2: Calculate log fugacities and their derivatives ---
  std::vector< double > X_i( numComps );
  std::vector< std::vector< double > > dXid_xj( numComps, std::vector< double >( numComps ));
  calculateFugacity< USD >( numComps,
                            pressure,
                            temperature,
                            componentProperties,
                            flashData,
                            flashData.liquidEos,
                            x_i,
                            X_i,
                            dXid_xj );

  std::vector< double > Y_i( numComps );
  std::vector< std::vector< double > > dYid_yj( numComps, std::vector< double >( numComps ));
  calculateFugacity< USD >( numComps,
                            pressure,
                            temperature,
                            componentProperties,
                            flashData,
                            flashData.vapourEos,
                            y_i,
                            Y_i,
                            dYid_yj );

  // --- Part 3: Calculate the residual vector ---
  residual.resize( numComps + 1 );
  for( int i = 0; i < numComps; ++i )
  {
    residual[i] = x_i[i] * std::exp( X_i[i] ) - y_i[i] * std::exp( Y_i[i] );
  }

  double fV_sum = 0.0;
  for( int i = 0; i < numComps; ++i )
  {
    fV_sum += totalComposition[i] * (kValues[i] - 1.0) / (1.0 + V * (kValues[i] - 1.0));
  }
  residual[numComps] = fV_sum;


  // --- Part 4: Calculate the Jacobian matrix ---
  jacobian.assign( numComps + 1, std::vector< double >( numComps + 1, 0.0 ));

  // Pre-calculate intermediate derivatives to avoid redundant calculations
  std::vector< double > d_xs_d_V( numComps );
  std::vector< double > d_ys_d_V( numComps );
  double sum_d_xs_d_V = 0.0;
  double sum_d_ys_d_V = 0.0;
  for( int i = 0; i < numComps; ++i )
  {
    double denominator = 1.0 + V * (kValues[i] - 1.0);
    d_xs_d_V[i] = -totalComposition[i] * (kValues[i] - 1.0) / (denominator * denominator);
    d_ys_d_V[i] = kValues[i] * d_xs_d_V[i];
    sum_d_xs_d_V += d_xs_d_V[i];
    sum_d_ys_d_V += d_ys_d_V[i];
  }

  std::vector< double > d_xsj_d_kj( numComps );
  std::vector< double > d_ysj_d_kj( numComps );
  for( int j = 0; j < numComps; ++j )
  {
    double denominator = 1.0 + V * (kValues[j] - 1.0);
    d_xsj_d_kj[j] = -totalComposition[j] * V * kValues[j] / (denominator * denominator);
    d_ysj_d_kj[j] = kValues[j] * x_s[j] + kValues[j] * d_xsj_d_kj[j];
  }

  // Fill the top-left block J_ij (i=1..N, j=1..N)
  for( int i = 0; i < numComps; ++i )
  {
    for( int j = 0; j < numComps; ++j )
    {
      double d_xi_d_kj = d_xsj_d_kj[j] * (kronecker_delta( i, j ) / x_sum - x_s[i] / (x_sum * x_sum));
      double d_yi_d_kj = d_ysj_d_kj[j] * (kronecker_delta( i, j ) / y_sum - y_s[i] / (y_sum * y_sum));

      double sum_dXid_xk_dxk_d_kj = 0.0;
      double sum_dYid_yk_dyk_d_kj = 0.0;
      for( int k = 0; k < numComps; ++k )
      {
        double d_xk_d_kj = d_xsj_d_kj[j] * (kronecker_delta( k, j ) / x_sum - x_s[k] / (x_sum * x_sum));
        double d_yk_d_kj = d_ysj_d_kj[j] * (kronecker_delta( k, j ) / y_sum - y_s[k] / (y_sum * y_sum));
        sum_dXid_xk_dxk_d_kj += dXid_xj[i][k] * d_xk_d_kj;
        sum_dYid_yk_dyk_d_kj += dYid_yj[i][k] * d_yk_d_kj;
      }

      jacobian[i][j] = std::exp( X_i[i] ) * d_xi_d_kj + x_i[i] * std::exp( X_i[i] ) * sum_dXid_xk_dxk_d_kj -
                       std::exp( Y_i[i] ) * d_yi_d_kj - y_i[i] * std::exp( Y_i[i] ) * sum_dYid_yk_dyk_d_kj;
    }
  }

  // Fill the top-right block J_i,N+1 (i=1..N)
  for( int i = 0; i < numComps; ++i )
  {
    double d_xi_d_V = (d_xs_d_V[i] * x_sum - x_s[i] * sum_d_xs_d_V) / (x_sum * x_sum);
    double d_yi_d_V = (d_ys_d_V[i] * y_sum - y_s[i] * sum_d_ys_d_V) / (y_sum * y_sum);

    double sum_dXid_xk_dxk_d_V = 0.0;
    double sum_dYid_yk_dyk_d_V = 0.0;
    for( int k = 0; k < numComps; ++k )
    {
      double d_xk_d_V = (d_xs_d_V[k] * x_sum - x_s[k] * sum_d_xs_d_V) / (x_sum * x_sum);
      double d_yk_d_V = (d_ys_d_V[k] * y_sum - y_s[k] * sum_d_ys_d_V) / (y_sum * y_sum);
      sum_dXid_xk_dxk_d_V += dXid_xj[i][k] * d_xk_d_V;
      sum_dYid_yk_dyk_d_V += dYid_yj[i][k] * d_yk_d_V;
    }

    jacobian[i][numComps] = std::exp( X_i[i] ) * d_xi_d_V + x_i[i] * std::exp( X_i[i] ) * sum_dXid_xk_dxk_d_V -
                            std::exp( Y_i[i] ) * d_yi_d_V - y_i[i] * std::exp( Y_i[i] ) * sum_dYid_yk_dyk_d_V;
  }

  // Fill the bottom-left block J_N+1,j (j=1..N)
  for( int j = 0; j < numComps; ++j )
  {
    double denominator = 1.0 + V * (kValues[j] - 1.0);
    jacobian[numComps][j] = totalComposition[j] * kValues[j] / (denominator * denominator);
  }

  // Fill the bottom-right entry J_N+1,N+1
  double sum_j_N_V = 0.0;
  for( int i = 0; i < numComps; ++i )
  {
    double denominator = 1.0 + V * (kValues[i] - 1.0);
    sum_j_N_V += -totalComposition[i] * std::pow( kValues[i] - 1.0, 2 ) / (denominator * denominator);
  }
  jacobian[numComps][numComps] = sum_j_N_V;
}

} // namespace compositional
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_IMPL_HPP_
