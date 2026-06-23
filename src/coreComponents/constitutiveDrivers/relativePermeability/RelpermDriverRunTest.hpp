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
struct HasHysteresis : std::false_type {};

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

    integer ipWetting, ipNonWetting;
    std::tie( ipWetting, ipNonWetting ) = relperm.wettingAndNonWettingPhaseIndices();

    // For the wetting phase, we need the minimum saturation
    // For the non-wetting phase, we need the maximum saturation
    real64 maxNonWetting = 0.0;
    real64 minWetting = 1.0;

    arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = relperm.template getField< fields::relperm::phaseMaxHistoricalVolFraction >().reference();
    arrayView2d< real64, compflow::USD_PHASE > phaseMinHistoricalVolFraction = relperm.template getField< fields::relperm::phaseMinHistoricalVolFraction >().reference();

    for( integer step = 0; step < numRows; ++step )
    {
      real64 const sw = table( step, SATURATION + ipWetting );
      real64 const snw = table( step, SATURATION + ipNonWetting );

      minWetting = LvArray::math::min( minWetting, sw );
      maxNonWetting = LvArray::math::max( maxNonWetting, snw );

      phaseMinHistoricalVolFraction( step, ipWetting ) = minWetting;
      phaseMaxHistoricalVolFraction( step, ipNonWetting ) = maxNonWetting;

      table( step, OFFSET + ipWetting ) = minWetting;
      table( step, OFFSET + ipNonWetting ) = maxNonWetting;
    }
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
