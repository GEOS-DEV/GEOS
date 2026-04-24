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

#ifndef GEOS_CONSTITUTIVEDRIVERS_RELATIVEPERMEABILITY_RELPERMDRIVERRUNTEST_HPP
#define GEOS_CONSTITUTIVEDRIVERS_RELATIVEPERMEABILITY_RELPERMDRIVERRUNTEST_HPP

#include "constitutiveDrivers/relativePermeability/RelpermDriver.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "constitutive/relativePermeability/KilloughHysteresis.hpp"
#include "common/DataLayouts.hpp"

namespace geos
{

// Hysteresis traits
template< typename RELPERM_TYPE >
struct HasHysteresis : std::false_type
{};

template< typename RELPERM_TYPE >
void
RelpermDriver::runTest( RELPERM_TYPE & relperm, const arrayView2d< real64 > & table )
{
  // Get the number of phases
  integer const numPhases = relperm.numFluidPhases();

  // Create the kernel wrapper
  typename RELPERM_TYPE::KernelWrapper const kernelWrapper = relperm.createKernelWrapper();

  // Offset for saturations in table
  constexpr integer SATURATION = 1;

  // Offset for relative permeability data
  integer RELPERM = SATURATION + numPhases;

  // Number of "cells"
  integer const numRows = m_table.size( 0 );

  // If we have hysteresis, we need to populate the historical saturations
  if constexpr (HasHysteresis< RELPERM_TYPE >::value)
  {
    // Shift the relative permeability offset
    RELPERM += numPhases;

    // Offset for the historical saturations
    integer const OFFSET = SATURATION + numPhases;

    using CurveType = constitutive::KilloughHysteresis::HysteresisCurve;
    using KeyStruct = typename RELPERM_TYPE::viewKeyStruct;

    auto [ipWetting, ipNonWetting] = relperm.wettingAndNonWettingPhaseIndices();

    stackArray1d< real64, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES > historicalSaturations( numPhases );
    historicalSaturations.zero();

    if( m_historicalSaturations.empty())
    {
      auto const & phaseHasHysteresis = relperm.template getReference< array1d< integer > >( KeyStruct::phaseHasHysteresisString());
      if( phaseHasHysteresis[ipNonWetting] )
      {
        CurveType const nonWettingCurve = relperm.template getReference< CurveType >( KeyStruct::nonWettingCurveString() );
        historicalSaturations[ipNonWetting] = nonWettingCurve.m_extremumPhaseVolFraction;
      }
      if( phaseHasHysteresis[ipWetting] )
      {
        CurveType const wettingCurve = relperm.template getReference< CurveType >( KeyStruct::wettingCurveString());
        historicalSaturations[ipWetting] = wettingCurve.m_extremumPhaseVolFraction;
      }
    }
    else
    {
      for( integer p = 0; p < numPhases; ++p )
      {
        historicalSaturations[p] = m_historicalSaturations[p];
      }
    }
    auto historicalView = historicalSaturations.toView();

    arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = relperm.template getField< fields::relperm::phaseMaxHistoricalVolFraction >().reference();
    arrayView2d< real64, compflow::USD_PHASE > phaseMinHistoricalVolFraction = relperm.template getField< fields::relperm::phaseMinHistoricalVolFraction >().reference();
    forAll< parallelDevicePolicy<> >( numRows,
                                      [numPhases, table,
                                       ipWetting, ipNonWetting,
                                       OFFSET,
                                       phaseMaxHistoricalVolFraction,
                                       phaseMinHistoricalVolFraction,
                                       historicalView] GEOS_HOST_DEVICE ( integer const n )
    {
      phaseMaxHistoricalVolFraction( n, ipNonWetting ) = historicalView[ipNonWetting];
      phaseMinHistoricalVolFraction( n, ipWetting ) = historicalView[ipWetting];
      for( integer p = 0; p < numPhases; ++p )
      {
        table( n, p + OFFSET ) = historicalView[p];
      }
    } );
  }

  forAll< parallelDevicePolicy<> >( numRows,
                                    [numPhases, kernelWrapper, table,
                                     RELPERM] GEOS_HOST_DEVICE ( integer const n )
  {
    StackArray< real64, 2, constitutive::RelativePermeabilityBase::MAX_NUM_PHASES, compflow::LAYOUT_PHASE > saturation( 1, numPhases );

    for( integer p = 0; p < numPhases; ++p )
    {
      saturation[0][p] = table( n, SATURATION + p );
    }
    kernelWrapper.update( n, 0, saturation[0] );
    for( integer p = 0; p < numPhases; ++p )
    {
      table( n, RELPERM + p ) = kernelWrapper.relperm()( n, 0, p );
    }
  } );
}

}

#endif //GEOS_CONSTITUTIVEDRIVERS_RELATIVEPERMEABILITY_RELPERMDRIVERRUNTEST_HPP
