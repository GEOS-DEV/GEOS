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
 * @file InverseCapillaryPressure.cpp
 */

#include "InverseCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/BrooksCoreyCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/VanGenuchtenCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/JFunctionCapillaryPressure.hpp"

namespace geos
{

namespace constitutive
{

template< typename CAP_PRESSURE >
InverseCapillaryPressureUpdate< CAP_PRESSURE >::InverseCapillaryPressureUpdate( CAP_PRESSURE & capPressure,
                                                                                arrayView2d< real64 const > const & propertyLimits,
                                                                                arrayView1d< integer const > const & independentPhases,
                                                                                integer dependentPhase,
                                                                                arrayView1d< integer const > const & jFunctionIndex )
  : m_capPressureWrapper( capPressure.createKernelWrapper() ),
  m_propertyLimits( propertyLimits ),
  m_independentPhases( independentPhases ),
  m_dependentPhase( dependentPhase ),
  m_jFunctionIndex( jFunctionIndex )
{
  auto const minSaturation = propertyLimits[MIN_PORE_VOLUME];
  for( real64 const saturation : minSaturation )
  {
    m_sumMinVolumeFraction += saturation;
  }
}

template< typename CAP_PRESSURE >
InverseCapillaryPressure< CAP_PRESSURE >::InverseCapillaryPressure( CAP_PRESSURE & capPressure )
  : m_capPressure( capPressure ),
  m_propertyLimits( 5, capPressure.numFluidPhases())
{
  string const mivVolumeKey = CapillaryPressureBase::viewKeyStruct::phaseMinVolumeFractionString();
  array1d< real64 > & phaseMinVolumeFraction = capPressure.template getReference< array1d< real64 > >( mivVolumeKey );

  string const phaseOrderKey = CapillaryPressureBase::viewKeyStruct::phaseOrderString();
  array1d< integer > const & phaseOrder = capPressure.template getReference< array1d< integer > >( phaseOrderKey );;

  calculatePropertyLimits( capPressure.numFluidPhases(),
                           capPressure.createKernelWrapper(),
                           phaseMinVolumeFraction,
                           m_propertyLimits );

  calculateIndependentPhases( capPressure.numFluidPhases(),
                              m_propertyLimits.toViewConst(),
                              m_dependentPhase,
                              m_independentPhases );

  calculateJFunctionIndex( capPressure.numFluidPhases(),
                           phaseOrder.toViewConst(),
                           m_jFunctionIndex );
}

template< typename CAP_PRESSURE >
typename InverseCapillaryPressure< CAP_PRESSURE >::KernelWrapper
InverseCapillaryPressure< CAP_PRESSURE >::createKernelWrapper()
{
  return KernelWrapper( m_capPressure,
                        m_propertyLimits.toViewConst(),
                        m_independentPhases.toViewConst(),
                        m_dependentPhase,
                        m_jFunctionIndex.toViewConst() );
}

template< typename CAP_PRESSURE >
void InverseCapillaryPressure< CAP_PRESSURE >::calculatePropertyLimits( integer numPhases,
                                                                        typename CAP_PRESSURE::KernelWrapper capPressureWrapper,
                                                                        arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                                                        arrayView2d< real64 > const & propertyLimits ) const
{
  static constexpr integer MAX_NUM_PHASES = CapillaryPressureBase::MAX_NUM_PHASES;

  auto const minPoreVolume = propertyLimits[KernelWrapper::MIN_PORE_VOLUME];
  auto const minSaturation = propertyLimits[KernelWrapper::MIN_SATURATION];
  auto const maxSaturation = propertyLimits[KernelWrapper::MAX_SATURATION];
  auto const minCapPressure = propertyLimits[KernelWrapper::MIN_CAP_PRESSURE];
  auto const maxCapPressure = propertyLimits[KernelWrapper::MAX_CAP_PRESSURE];

  real64 sumMinVolumeFraction = 0.0;
  stackArray1d< real64, MAX_NUM_PHASES > jFunctionMultiplier( numPhases );
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    minCapPressure[ip] = LvArray::NumericLimits< real64 >::max;
    maxCapPressure[ip] = -LvArray::NumericLimits< real64 >::max;
    minPoreVolume[ip] = phaseMinVolumeFraction[ip];
    sumMinVolumeFraction += phaseMinVolumeFraction[ip];
    jFunctionMultiplier[ip] = 1.0;  // Call with multipliers of 1
  }

  StackArray< real64, 2, MAX_NUM_PHASES, compflow::LAYOUT_PHASE > saturation( 1, numPhases );
  StackArray< real64, 3, MAX_NUM_PHASES, cappres::LAYOUT_CAPPRES > workSpace( 1, 1, numPhases );
  StackArray< real64, 4, MAX_NUM_PHASES *MAX_NUM_PHASES, cappres::LAYOUT_CAPPRES_DS > dPhaseCapPres_dSaturation( 1, 1, numPhases, numPhases );

  auto const capPressure = workSpace[0][0];

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    for( integer jp = 0; jp < numPhases; ++jp )
    {
      saturation[0][jp] = phaseMinVolumeFraction[jp];
    }
    saturation[0][ip] = 1.0 - sumMinVolumeFraction + phaseMinVolumeFraction[ip];
    CapillaryPressureEvaluate< CAP_PRESSURE >::compute( capPressureWrapper,
                                                        saturation[0].toSliceConst(),
                                                        jFunctionMultiplier.toSliceConst(),
                                                        capPressure,
                                                        dPhaseCapPres_dSaturation[0][0] );
    for( integer jp = 0; jp < numPhases; ++jp )
    {
      if( capPressure[jp] < minCapPressure[jp] )
      {
        minSaturation[jp] = saturation[0][jp];
        minCapPressure[jp] = capPressure[jp];
      }
      if( maxCapPressure[jp] < capPressure[jp] )
      {
        maxSaturation[jp] = saturation[0][jp];
        maxCapPressure[jp] = capPressure[jp];
      }
    }
  }
}

template< typename CAP_PRESSURE >
void InverseCapillaryPressure< CAP_PRESSURE >::calculateIndependentPhases( integer numPhases,
                                                                           arrayView2d< real64 const > const & propertyLimits,
                                                                           integer & dependentPhase,
                                                                           array1d< integer > & independentPhases ) const
{
  // The dependent phase is the one with the least variation in capillary pressure
  auto const minCapPressure = propertyLimits[KernelWrapper::MIN_CAP_PRESSURE];
  auto const maxCapPressure = propertyLimits[KernelWrapper::MAX_CAP_PRESSURE];

  // Choose one of the phases to be dependent
  real64 minDP = LvArray::NumericLimits< real64 >::max;
  dependentPhase = -1;
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    real64 const dp = maxCapPressure[ip] - minCapPressure[ip];
    if( dp < minDP )
    {
      minDP = dp;
      dependentPhase = ip;
    }
  }

  independentPhases.clear();
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    if( ip != dependentPhase )
    {
      independentPhases.emplace_back( ip );
    }
  }
}

template< typename CAP_PRESSURE >
void InverseCapillaryPressure< CAP_PRESSURE >::calculateJFunctionIndex( integer numPhases,
                                                                        arrayView1d< integer const > const & phaseOrder,
                                                                        array1d< integer > & jFunctionIndex ) const
{
  jFunctionIndex.resize( numPhases );
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    jFunctionIndex[ip] = 0;
  }
  if constexpr (std::is_same_v< CAP_PRESSURE, JFunctionCapillaryPressure >)
  {
    if( numPhases == 3 )
    {
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        if( phaseOrder[ip] == CapillaryPressureBase::PhaseType::WATER )
        {
          jFunctionIndex[ip] = JFunctionCapillaryPressure::ThreePhasePairPhaseType::INTERMEDIATE_WETTING;
        }
        if( phaseOrder[ip] == CapillaryPressureBase::PhaseType::GAS )
        {
          jFunctionIndex[ip] = JFunctionCapillaryPressure::ThreePhasePairPhaseType::INTERMEDIATE_NONWETTING;
        }
      }
    }
  }
  else
  {
    GEOS_UNUSED_VAR( phaseOrder );
  }
}

template class InverseCapillaryPressure< BrooksCoreyCapillaryPressure >;
template class InverseCapillaryPressure< TableCapillaryPressure >;
template class InverseCapillaryPressure< JFunctionCapillaryPressure >;
template class InverseCapillaryPressure< VanGenuchtenCapillaryPressure >;

} // namespace constitutive

} // namespace geos
