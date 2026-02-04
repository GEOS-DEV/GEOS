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

NoOpCapillaryPressure::NoOpCapillaryPressure( integer const numPhases,
                                              arrayView1d< integer const > const & phaseOrder )
  : m_numPhases( numPhases ),
  m_phaseOrder( phaseOrder )
{}

integer NoOpCapillaryPressure::numFluidPhases() const
{
  return m_numPhases;
}

arrayView1d< integer const > const & NoOpCapillaryPressure::phaseOrder() const
{
  return m_phaseOrder;
}

NoOpCapillaryPressure::KernelWrapper NoOpCapillaryPressure::createKernelWrapper()
{
  return KernelWrapper();
}

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
  m_propertyLimits( 5, capPressure.numFluidPhases() )
{
  string const mivVolumeKey = CapillaryPressureBase::viewKeyStruct::phaseMinVolumeFractionString();
  array1d< real64 > & phaseMinVolumeFraction = capPressure.template getReference< array1d< real64 > >( mivVolumeKey );

  string const phaseOrderKey = CapillaryPressureBase::viewKeyStruct::phaseOrderString();
  array1d< integer > const & phaseOrder = capPressure.template getReference< array1d< integer > >( phaseOrderKey );

  calculatePropertyLimits( capPressure.numFluidPhases(),
                           phaseOrder,
                           capPressure.createKernelWrapper(),
                           phaseMinVolumeFraction,
                           m_propertyLimits );

  calculateIndependentPhases( capPressure.numFluidPhases(),
                              phaseOrder.toViewConst(),
                              m_dependentPhase,
                              m_independentPhases );

  calculateJFunctionIndex( capPressure.numFluidPhases(),
                           phaseOrder.toViewConst(),
                           m_jFunctionIndex );
}

template< >
InverseCapillaryPressure< NoOpCapillaryPressure >::InverseCapillaryPressure( NoOpCapillaryPressure & capPressure )
  : m_capPressure( capPressure ),
  m_propertyLimits( 5, capPressure.numFluidPhases() ),
  m_jFunctionIndex( capPressure.numFluidPhases() )
{
  integer const numPhases = capPressure.numFluidPhases();
  arrayView1d< integer const > const & phaseOrder = capPressure.phaseOrder();

  auto const minPoreVolume = m_propertyLimits[KernelWrapper::MIN_PORE_VOLUME];
  auto const minSaturation = m_propertyLimits[KernelWrapper::MIN_SATURATION];
  auto const maxSaturation = m_propertyLimits[KernelWrapper::MAX_SATURATION];
  auto const minCapPressure = m_propertyLimits[KernelWrapper::MIN_CAP_PRESSURE];
  auto const maxCapPressure = m_propertyLimits[KernelWrapper::MAX_CAP_PRESSURE];

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    minCapPressure[ip] = 0.0;
    maxCapPressure[ip] = 0.0;
    minPoreVolume[ip] = 0.0;
    m_jFunctionIndex[ip] = 0;
  }

  for( integer phaseType = 0; phaseType < phaseOrder.size(); ++phaseType )
  {
    integer const phaseIndex = phaseOrder[phaseType];
    if( phaseIndex < 0 )
    {
      continue;
    }
    // The water indexed capillary pressure is supposed to be decreasing
    bool const isDecreasing = (phaseType == CapillaryPressureBase::PhaseType::WATER);
    if( isDecreasing )
    {
      minSaturation[phaseIndex] = 1.0;
      maxSaturation[phaseIndex] = 0.0;
    }
    else
    {
      minSaturation[phaseIndex] = 0.0;
      maxSaturation[phaseIndex] = 1.0;
    }
  }

  calculateIndependentPhases( capPressure.numFluidPhases(),
                              phaseOrder.toViewConst(),
                              m_dependentPhase,
                              m_independentPhases );
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
                                                                        arrayView1d< integer const > const & phaseOrder,
                                                                        typename CAP_PRESSURE::KernelWrapper capPressureWrapper,
                                                                        arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                                                        arrayView2d< real64 > const & propertyLimits ) const
{
  constexpr integer MAX_NUM_PHASES = CapillaryPressureBase::MAX_NUM_PHASES;

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

  // If the capillary pressure variation in a phase is zero then we need to change the mimimum and maximum saturations
  auto const getPhase = [&]( integer const phaseIndex ) -> integer {
    for( integer ip = 0; ip < phaseOrder.size(); ++ip )
    {
      if( phaseOrder[ip] == phaseIndex )
      {
        return ip;
      }
    }
    return -1;
  };
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    real64 const dp = maxCapPressure[ip] - minCapPressure[ip];
    if( dp < LvArray::NumericLimits< real64 >::epsilon )
    {
      bool const isDecreasing = (getPhase( ip ) != CapillaryPressureBase::PhaseType::OIL);
      minSaturation[ip] = phaseMinVolumeFraction[ip];
      maxSaturation[ip] = 1.0 - sumMinVolumeFraction + phaseMinVolumeFraction[ip];
      if( isDecreasing )
      {
        std::swap( minSaturation[ip], maxSaturation[ip] );
      }
    }
  }
}

template< typename CAP_PRESSURE >
void InverseCapillaryPressure< CAP_PRESSURE >::calculateIndependentPhases( integer numPhases,
                                                                           arrayView1d< integer const > const & phaseOrder,
                                                                           integer & dependentPhase,
                                                                           array1d< integer > & independentPhases ) const
{
  // Precedence of phases on independence
  // Always favour water as the independent phase, then gas, then oil
  if( 0 <= phaseOrder[CapillaryPressureBase::PhaseType::OIL] )
  {
    dependentPhase = phaseOrder[CapillaryPressureBase::PhaseType::OIL];
  }
  else if( 0 <= phaseOrder[CapillaryPressureBase::PhaseType::GAS] )
  {
    dependentPhase = phaseOrder[CapillaryPressureBase::PhaseType::GAS];
  }
  else
  {
    dependentPhase = 0;
  }

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    if( ip != dependentPhase )
    {
      independentPhases.emplace_back( ip );
    }
  }

  // Ensure water is the first independent phase
  if( numPhases == 3 )
  {
    if( phaseOrder[CapillaryPressureBase::PhaseType::GAS] == independentPhases[0] )
    {
      std::swap( independentPhases[0], independentPhases[1] );
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
template class InverseCapillaryPressure< NoOpCapillaryPressure >;

} // namespace constitutive

} // namespace geos
