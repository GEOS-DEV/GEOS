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

class JFunctionCapillaryPressure;


struct NoOpCapillaryPressure
{
  struct KernelWrapper
  {
    GEOS_HOST_DEVICE
    void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                  arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                  arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dPhaseVolFrac ) const
    {
      GEOS_UNUSED_VAR( phaseVolFraction );
      LvArray::forValuesInSlice( phaseCapPres, setZero );
      LvArray::forValuesInSlice( dPhaseCapPres_dPhaseVolFrac, setZero );
    }

    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    static void setZero( real64 & value )
    {
      value = 0.0;
    }
  };

  NoOpCapillaryPressure( integer const numPhases, arrayView1d< integer const > const & phaseOrder );
  KernelWrapper createKernelWrapper();

  integer numFluidPhases() const;
  arrayView1d< integer const > const & phaseOrder() const;

private:
  integer const m_numPhases;
  arrayView1d< integer const > const & m_phaseOrder;
};

template< typename CAP_PRESSURE >
struct CapillaryPressureEvaluate
{
  GEOS_HOST_DEVICE
  static void compute( typename CAP_PRESSURE::KernelWrapper const & capPressureWrapper,
                       arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolumeFraction,
                       arraySlice1d< real64 const > const & jFuncMultiplier,
                       arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPres,
                       arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPres_dSaturation )
  {
    constexpr bool isJFunction = std::is_same_v< CAP_PRESSURE, JFunctionCapillaryPressure >;
    if constexpr ( isJFunction )
    {
      capPressureWrapper.compute( phaseVolumeFraction, jFuncMultiplier, phaseCapPres, dPhaseCapPres_dSaturation );
    }
    else
    {
      GEOS_UNUSED_VAR( jFuncMultiplier );
      capPressureWrapper.compute( phaseVolumeFraction, phaseCapPres, dPhaseCapPres_dSaturation );
    }
  }

  GEOS_HOST_DEVICE
  static real64 applyScale( real64 const capPressure,
                            real64 const jFuncMultiplier )
  {
    constexpr bool isJFunction = std::is_same_v< CAP_PRESSURE, JFunctionCapillaryPressure >;
    if constexpr ( isJFunction )
    {
      return capPressure * jFuncMultiplier;
    }
    else
    {
      GEOS_UNUSED_VAR( jFuncMultiplier );
      return capPressure;
    }
  }
};

template< typename CAP_PRESSURE >
class InverseCapillaryPressureUpdate
{
public:
  static constexpr integer MAX_ITERATIONS = 20;
  static constexpr real64 STEP_TOLERANCE = 1.0e-6;
  static constexpr integer USD_SAT = compflow::USD_PHASE - 1;
  static constexpr integer USD_PC = cappres::USD_CAPPRES - 2;

public:
  // Property limits. Used to index `m_propertyLimits`
  // Minimum phase volume fraction for phase
  static constexpr integer MIN_PORE_VOLUME = 0;
  // End points for the capillary pressure
  static constexpr integer MIN_CAP_PRESSURE = 1;
  static constexpr integer MAX_CAP_PRESSURE = 2;
  // End points for the saturation dimension
  // Min saturation is saturation at which PC is minimum
  // Max saturation is saturation at which PC is maximum
  static constexpr integer MIN_SATURATION = 3;
  static constexpr integer MAX_SATURATION = 4;
public:
  InverseCapillaryPressureUpdate( CAP_PRESSURE & capPressure,
                                  arrayView2d< real64 const > const & propertyLimits,
                                  arrayView1d< integer const > const & independentPhases,
                                  integer dependentPhase,
                                  arrayView1d< integer const > const & jFunctionIndex );

  GEOS_HOST_DEVICE
  bool compute( arraySlice1d< real64 const, USD_PC > const & phaseCapillaryPressure,
                arraySlice1d< real64 const > const & jFunctionMultiplier,
                arraySlice1d< real64, USD_SAT > const & phaseVolumeFraction ) const
  {
    constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;
    constexpr integer MAX_NUM_PHASES = CapillaryPressureBase::MAX_NUM_PHASES;

    integer const numPhases = phaseVolumeFraction.size();
    StackArray< real64, 3, 3*MAX_NUM_PHASES, cappres::LAYOUT_CAPPRES > workSpace( 1, 3, numPhases );
    StackArray< real64, 4, MAX_NUM_PHASES *MAX_NUM_PHASES, cappres::LAYOUT_CAPPRES_DS > dPhaseCapPres_dSaturation( 1, 1, numPhases, numPhases );

    auto const & minSaturation = m_propertyLimits[MIN_SATURATION];
    auto const & maxSaturation = m_propertyLimits[MAX_SATURATION];
    auto const & minCapPressure = m_propertyLimits[MIN_CAP_PRESSURE];
    auto const & maxCapPressure = m_propertyLimits[MAX_CAP_PRESSURE];

    auto const capPres = workSpace[0][0];
    auto const targetCapPres = workSpace[0][1];
    auto const saturationStep = workSpace[0][2];
    auto const & saturation = phaseVolumeFraction;
    auto const jacobian = dPhaseCapPres_dSaturation[0][0];

    // Initial solution
    real64 sumSaturations = 0.0;
    for( integer const ip : m_independentPhases )
    {
      real64 S = phaseVolumeFraction[ip];
      real64 const minPc = CapillaryPressureEvaluate< CAP_PRESSURE >
                           ::applyScale( minCapPressure[ip],
                                         jFunctionMultiplier[m_jFunctionIndex[ip]] );
      real64 const maxPc = CapillaryPressureEvaluate< CAP_PRESSURE >
                           ::applyScale( maxCapPressure[ip],
                                         jFunctionMultiplier[m_jFunctionIndex[ip]] );

      targetCapPres[ip] = phaseCapillaryPressure[ip];
      targetCapPres[ip] = LvArray::math::max( minPc, targetCapPres[ip] );
      targetCapPres[ip] = LvArray::math::min( maxPc, targetCapPres[ip] );

      if( phaseCapillaryPressure[ip] - STEP_TOLERANCE < minPc )
      {
        S = minSaturation[ip];
      }
      else if( maxPc < phaseCapillaryPressure[ip] + STEP_TOLERANCE )
      {
        S = maxSaturation[ip];
      }
      else
      {
        S = LvArray::math::max( S, LvArray::math::min( maxSaturation[ip], minSaturation[ip] ));
        S = LvArray::math::min( S, LvArray::math::max( maxSaturation[ip], minSaturation[ip] ));
      }
      saturation[ip] = S;
      sumSaturations += S;
    }
    saturation[m_dependentPhase] = 1.0 - sumSaturations;
    normalizeSaturations( numPhases, m_dependentPhase, saturation );

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
      for( integer const ip : m_independentPhases )
      {
        real64 const dp = targetCapPres[ip] - capPres[ip];
        real64 const Aii = jacobian[ip][ip];
        real64 const dS = ( epsilon < LvArray::math::abs( Aii )) ? dp/Aii : 0.0;
        if( epsilon < LvArray::math::abs( dS ))
        {
          // Calculate a step size that does not take us outside the bounds
          real64 phaseStepSize = LvArray::math::max((minSaturation[ip] - saturation[ip])/dS, (maxSaturation[ip] - saturation[ip])/dS );
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
      for( integer const ip : m_independentPhases )
      {
        real64 const dS = stepSize* saturationStep[ip];
        saturation[ip] += dS;
        sumSaturations += saturation[ip];
        saturationChange += dS*dS;
      }
      saturation[m_dependentPhase] = 1.0 - sumSaturations;

      // Check for convergence: when saturations stop changing
      saturationChange = LvArray::math::sqrt( saturationChange );
      if( saturationChange < STEP_TOLERANCE )
      {
        converged = true;
      }
    }

    // For 3-phase cases the oil phase might end up with a negative saturation. This is an
    // indication that we have reached the gas-oil contact while still in the oil transision
    // zone. In this case, we will solve the problem as a 2-phase gas-water system.
    if( 3 <= numPhases && saturation[m_dependentPhase] < -LvArray::NumericLimits< real64 >::epsilon )
    {
      saturation[m_dependentPhase] = 0.0;

      integer const phase0 = m_independentPhases[0];
      integer const phase1 = m_independentPhases[1];

      real64 const minPhase0Pc = CapillaryPressureEvaluate< CAP_PRESSURE >::applyScale( minCapPressure[phase0],
                                                                                        jFunctionMultiplier[m_jFunctionIndex[phase0]] )
                                 - CapillaryPressureEvaluate< CAP_PRESSURE >::applyScale( maxCapPressure[phase1],
                                                                                          jFunctionMultiplier[m_jFunctionIndex[phase1]] );

      real64 const maxPhase0Pc = CapillaryPressureEvaluate< CAP_PRESSURE >::applyScale( maxCapPressure[phase0],
                                                                                        jFunctionMultiplier[m_jFunctionIndex[phase0]] )
                                 - CapillaryPressureEvaluate< CAP_PRESSURE >::applyScale( minCapPressure[phase1],
                                                                                          jFunctionMultiplier[m_jFunctionIndex[phase1]] );

      real64 const targetPc = phaseCapillaryPressure[phase0] - phaseCapillaryPressure[phase1];

      converged = false;
      if( targetPc - STEP_TOLERANCE < minPhase0Pc )
      {
        saturation[phase0] = minSaturation[phase0];
        converged = true;
      }
      else if( maxPhase0Pc < targetPc + STEP_TOLERANCE )
      {
        saturation[phase0] = maxSaturation[phase0];
        converged = true;
      }
      else
      {
        saturation[phase0] = 0.5*(minSaturation[phase0] + maxSaturation[phase0]);
      }
      saturation[phase1] = 1.0 - saturation[phase0];

      for( integer iterationCount = 0; (iterationCount < MAX_ITERATIONS) && !converged; ++iterationCount )
      {
        CapillaryPressureEvaluate< CAP_PRESSURE >::compute( m_capPressureWrapper,
                                                            saturation.toSliceConst(),
                                                            jFunctionMultiplier,
                                                            capPres,
                                                            jacobian );

        real64 const currentPc = capPres[phase0] - capPres[phase1];
        real64 const dp = targetPc - currentPc;
        real64 const Aii = jacobian[phase0][phase0] + jacobian[phase1][phase1];
        real64 const dS = ( epsilon < LvArray::math::abs( Aii )) ? dp/Aii : 0.0;
        real64 stepSize = 1.0;
        if( epsilon < LvArray::math::abs( dS ))
        {
          real64 phaseStepSize = LvArray::math::max((minSaturation[phase0] - saturation[phase0])/dS, (maxSaturation[phase0] - saturation[phase0])/dS );
          if( phaseStepSize < 1.0 )
          {
            phaseStepSize *= (1.0 - STEP_TOLERANCE);
          }
          stepSize = LvArray::math::min( stepSize, phaseStepSize );
        }

        saturation[phase0] += stepSize*dS;
        saturation[phase1] = 1.0 - saturation[phase0];

        converged = LvArray::math::abs( stepSize*dS ) < STEP_TOLERANCE;
      }
    }

    return converged;
  }

private:
  GEOS_HOST_DEVICE
  static void normalizeSaturations( integer numPhases,
                                    integer dependentPhase,
                                    arraySlice1d< real64, USD_SAT > const & phaseVolumeFraction );

private:
  typename CAP_PRESSURE::KernelWrapper m_capPressureWrapper;
  arrayView2d< real64 const > const m_propertyLimits;
  arrayView1d< integer const > const m_independentPhases;
  integer const m_dependentPhase{-1};
  arrayView1d< integer const > const m_jFunctionIndex;
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
  // Calculate the model limits
  void calculatePropertyLimits( integer numPhases,
                                arrayView1d< integer const > const & phaseOrder,
                                typename CAP_PRESSURE::KernelWrapper capPressureWrapper,
                                arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                arrayView2d< real64 > const & propertyLimits ) const;

  // Determine which phases are independent and which is dependent
  void calculateIndependentPhases( integer numPhases,
                                   arrayView1d< integer const > const & phaseOrder,
                                   integer & dependentPhase,
                                   array1d< integer > & independentPhases ) const;

  // Find the indices for the JFunction multiplier
  void calculateJFunctionIndex( integer numPhases,
                                arrayView1d< integer const > const & phaseOrder,
                                array1d< integer > & jFunctionIndex ) const;

private:
  CAP_PRESSURE & m_capPressure;

  /// Array for saturation and capillary pressure
  array2d< real64 > m_propertyLimits;

  // List of free phases
  array1d< integer > m_independentPhases;

  // The non-free phase
  integer m_dependentPhase{-1};

  // Indices for the JFunction multiplier for each phase
  array1d< integer > m_jFunctionIndex;
};

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
