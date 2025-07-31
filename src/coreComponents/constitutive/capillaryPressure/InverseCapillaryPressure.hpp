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
 * @file InverseCapillaryPressure.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_INVERSECAPILLARYPRESSURE_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_INVERSECAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

namespace geos
{
namespace constitutive
{

template< typename CAP_PRESSURE >
struct CapillaryPressureEvaluate
{
  GEOS_HOST_DEVICE
  static void compute( typename CAP_PRESSURE::KernelWrapper const & capPressureWrapper,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolumeFraction,
                       arraySlice1d< real64 const > const & jFuncMultiplier,
                       arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                       arraySlice2d< real64 > const & dPhaseCapPres_dSaturation )
  {
    GEOS_UNUSED_VAR( jFuncMultiplier );
    integer constexpr MAX_NUM_PHASES = CapillaryPressureBase::MAX_NUM_PHASES;
    integer const numPhases = phaseVolumeFraction.size();
    StackArray< real64, 4, MAX_NUM_PHASES *MAX_NUM_PHASES, constitutive::cappres::LAYOUT_CAPPRES_DS > dPhaseCapPres_dPhaseVolFrac( 1, 1, numPhases, numPhases );
    capPressureWrapper.compute( phaseVolumeFraction, phaseCapPres, dPhaseCapPres_dPhaseVolFrac[0][0] );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      for( integer jp = 0; jp < numPhases; ++jp )
      {
        dPhaseCapPres_dSaturation[ip][jp] = dPhaseCapPres_dPhaseVolFrac[0][0][ip][jp];
      }
    }
  }
};

template< typename CAP_PRESSURE >
class InverseCapillaryPressureUpdate
{
private:
  static constexpr integer MAX_NUM_PHASES = CapillaryPressureBase::MAX_NUM_PHASES;
  static constexpr integer MAX_ITERATIONS = 20;
  static constexpr real64 STEP_TOLERANCE = 1.0e-6;
  static constexpr integer USD_SAT = compflow::USD_PHASE - 1;
  static constexpr integer USD_PC = cappres::USD_CAPPRES - 2;

public:
  InverseCapillaryPressureUpdate( CAP_PRESSURE & capPressure,
                                  arrayView1d< real64 const > const & phaseMinVolumeFraction );

  GEOS_HOST_DEVICE
  bool compute( arraySlice1d< real64 const, USD_PC > const & phaseCapillapryPressure,
                arraySlice1d< real64 const > const & jFunctionMultiplier,
                arraySlice1d< real64, USD_SAT > const & phaseVolumeFraction ) const;

private:
  GEOS_HOST_DEVICE
  static void normalizeSaturations( integer numPhases,
                                    integer dependentPhase,
                                    arraySlice1d< real64, USD_SAT > const & phaseVolumeFraction );

private:
  typename CAP_PRESSURE::KernelWrapper m_capPressureWrapper;
  arrayView1d< real64 const > const m_phaseMinVolumeFraction;
  real64 m_sumMinVolumeFraction{0.0};
};

template< typename CAP_PRESSURE >
class InverseCapillaryPressure
{
public:
  InverseCapillaryPressure( CAP_PRESSURE & capPressure );

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = InverseCapillaryPressureUpdate< CAP_PRESSURE >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

private:
  CAP_PRESSURE & m_capPressure;
};

template< typename CAP_PRESSURE >
GEOS_HOST_DEVICE
bool
InverseCapillaryPressureUpdate< CAP_PRESSURE >::compute(
  arraySlice1d< real64 const, USD_PC > const & phaseCapillapryPressure,
  arraySlice1d< real64 const > const & jFunctionMultiplier,
  arraySlice1d< real64, USD_SAT > const & phaseVolumeFraction ) const
{
  constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;

  integer const numPhases = phaseVolumeFraction.size();
  StackArray< real64, 2, MAX_NUM_PHASES, compflow::LAYOUT_PHASE > sat( 1, numPhases );
  StackArray< real64, 3, 7*MAX_NUM_PHASES, constitutive::cappres::LAYOUT_CAPPRES > workSpace( 1, 7, numPhases );
  StackArray< real64, 2, MAX_NUM_PHASES *MAX_NUM_PHASES > dPhaseCapPres_dSaturation( numPhases, numPhases );
  StackArray< integer, 1, MAX_NUM_PHASES > freePhases( numPhases-1 );

  auto const minCapPres = workSpace[0][0];
  auto const maxCapPres = workSpace[0][1];
  auto const capPres = workSpace[0][2];
  auto const targetCapPres = workSpace[0][3];
  auto const saturationStep = workSpace[0][4];
  auto const minSat = workSpace[0][5];
  auto const maxSat = workSpace[0][6];
  auto const saturation = sat[0];
  auto const jacobian = dPhaseCapPres_dSaturation.toSlice();

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    minCapPres[ip] = LvArray::NumericLimits< real64 >::max;
    maxCapPres[ip] = -LvArray::NumericLimits< real64 >::max;
  }

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    for( integer jp = 0; jp < numPhases; ++jp )
    {
      saturation[jp] = m_phaseMinVolumeFraction[jp];
    }
    saturation[ip] = 1.0 - m_sumMinVolumeFraction + m_phaseMinVolumeFraction[ip];
    CapillaryPressureEvaluate< CAP_PRESSURE >::compute( m_capPressureWrapper,
                                                        saturation.toSliceConst(),
                                                        jFunctionMultiplier,
                                                        capPres,
                                                        jacobian );
    for( integer jp = 0; jp < numPhases; ++jp )
    {
      if( capPres[jp] < minCapPres[jp] )
      {
        minSat[jp] = saturation[jp];
        minCapPres[jp] = capPres[jp];
      }
      if( maxCapPres[jp] < capPres[jp] )
      {
        maxSat[jp] = saturation[jp];
        maxCapPres[jp] = capPres[jp];
      }
    }
  }

  // Choose one of the phases to be dependent
  real64 minDP = LvArray::NumericLimits< real64 >::max;
  integer dependentPhase = -1;
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    real64 const dp = maxCapPres[ip] - minCapPres[ip];
    if( dp < minDP )
    {
      minDP = dp;
      dependentPhase = ip;
    }
  }

  integer k = 0;
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    if( ip != dependentPhase )
    {
      freePhases[k++] = ip;
    }
  }

  // Use the limits to solve
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    targetCapPres[ip] = phaseCapillapryPressure[ip];
    targetCapPres[ip] = LvArray::math::max( minCapPres[ip], targetCapPres[ip] );
    targetCapPres[ip] = LvArray::math::min( maxCapPres[ip], targetCapPres[ip] );
  }

  // Initial solution
  real64 sumSaturations = 0.0;
  for( integer const ip : freePhases )
  {
    real64 S = phaseVolumeFraction[ip];
    if( phaseCapillapryPressure[ip] - STEP_TOLERANCE < minCapPres[ip] )
    {
      S = minSat[ip];
    }
    else if( maxCapPres[ip] < phaseCapillapryPressure[ip] + STEP_TOLERANCE )
    {
      S = maxSat[ip];
    }
    else
    {
      S = LvArray::math::max( S, LvArray::math::min( maxSat[ip], minSat[ip] ));
      S = LvArray::math::min( S, LvArray::math::max( maxSat[ip], minSat[ip] ));
    }
    saturation[ip] = S;
    sumSaturations += S;
  }
  saturation[dependentPhase] = 1.0 - sumSaturations;
  normalizeSaturations( numPhases, dependentPhase, saturation );

  bool converged = false;
  for( integer iterationCount = 0; (iterationCount < MAX_ITERATIONS) && !converged; ++iterationCount )
  {
    CapillaryPressureEvaluate< CAP_PRESSURE >::compute( m_capPressureWrapper,
                                                        saturation.toSliceConst(),
                                                        jFunctionMultiplier,
                                                        capPres,
                                                        jacobian );

    // Calculate Newton step
    // Assume diagonal Jacobian
    real64 stepSize = 1.0;
    for( integer const ip : freePhases )
    {
      real64 const dp = targetCapPres[ip] - capPres[ip];
      real64 const Aii = jacobian[ip][ip];
      real64 const dS = ( epsilon < LvArray::math::abs( Aii )) ? dp/Aii : 0.0;
      if( epsilon < LvArray::math::abs( dS ))
      {
        // Calculate a step size that does not take us outside the bounds
        real64 phaseStepSize = LvArray::math::max((maxSat[ip] - saturation[ip])/dS, (minSat[ip] - saturation[ip])/dS );
        if( phaseStepSize < 1.0 )
        {
          phaseStepSize *= (1.0 - STEP_TOLERANCE);
        }
        stepSize = LvArray::math::min( stepSize, phaseStepSize );
      }
      saturationStep[ip] = dS;
    }
    real64 saturationChange = 0.0;
    sumSaturations = 0.0;
    for( integer const ip : freePhases )
    {
      real64 const dS = stepSize* saturationStep[ip];
      saturation[ip] += dS;
      sumSaturations += saturation[ip];
      saturationChange += dS*dS;
    }
    saturation[dependentPhase] = 1.0 - sumSaturations;

    // Check for convergence: when saturations stop changing
    saturationChange = LvArray::math::sqrt( saturationChange );
    if( saturationChange < STEP_TOLERANCE )
    {
      converged = true;
    }
  }

  // Copy back the solution
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    phaseVolumeFraction[ip] = saturation[ip];
  }

  return converged;
}

template< typename CAP_PRESSURE >
GEOS_HOST_DEVICE
void
InverseCapillaryPressureUpdate< CAP_PRESSURE >::normalizeSaturations(
  integer numPhases,
  integer dependentPhase,
  arraySlice1d< real64, USD_SAT > const & phaseVolumeFraction )
{
  constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;
  real64 const freeSaturation = phaseVolumeFraction[dependentPhase];
  if( freeSaturation < epsilon )
  {
    real64 const scale = 1.0 / (1.0 - freeSaturation);
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseVolumeFraction[ip] *= scale;
    }
    phaseVolumeFraction[dependentPhase] = 0.0;
  }
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_INVERSECAPILLARYPRESSURE_HPP
