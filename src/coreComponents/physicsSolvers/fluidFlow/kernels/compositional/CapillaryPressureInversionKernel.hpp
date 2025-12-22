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
 * @file CapillaryPressureInversionKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_CAPILLARYPRESSUREINVERSIONKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_CAPILLARYPRESSUREINVERSIONKERNEL_HPP

#include "constitutive/capillaryPressure/InverseCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** CapillaryPressureInversionKernel ********************************/
/**
 * @brief Kernel for the inversion of capillary pressure to saturation when we have multiphase initialisation
 */
template< typename CAP_PRESSURE = NoOpFunc >
struct CapillaryPressureInversionKernel
{
  template< integer numComps, integer numPhases >
  static void launch( SortedArrayView< localIndex const > const & targetSet,
                      CAP_PRESSURE & capPressure,
                      arrayView1d< integer const > const & phaseOrder,
                      arrayView2d< real64 const > const & elementCenter,
                      TableFunction::KernelWrapper const & elevationIndexTable,
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & pressureValues,
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDensityValues,
                      arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > const & phaseComponentFractions,
                      arrayView2d< real64, compflow::USD_COMP > const globalComponentFractions )
  {
    constitutive::InverseCapillaryPressure< CAP_PRESSURE > inverseCapPressureType( capPressure );
    auto capPressureWrapper = inverseCapPressureType.createKernelWrapper();

    array2d< real64 > jFuncMultiplierArray( 1, numPhases - 1 );
    arrayView2d< real64 const > jFuncMultiplier = jFuncMultiplierArray.toViewConst();
    constexpr bool isJFunction = std::is_same_v< CAP_PRESSURE, constitutive::JFunctionCapillaryPressure >;
    if constexpr ( isJFunction )
    {
      jFuncMultiplier = capPressure.template getField< geos::fields::cappres::jFuncMultiplier >().reference().toViewConst();
    }

    integer const ipGas = phaseOrder[constitutive::CapillaryPressureBase::PhaseType::GAS];
    integer const ipWater = phaseOrder[constitutive::CapillaryPressureBase::PhaseType::WATER];
    integer const ipOil = phaseOrder[constitutive::CapillaryPressureBase::PhaseType::OIL];
    integer phases[numPhases]{};
    if constexpr ( numPhases == 2  )
    {
      if( 0 <= ipGas && 0 <= ipWater )
      {
        phases[0] = ipGas;
        phases[1] = ipWater;
      }
      if( 0 <= ipOil && 0 <= ipWater )
      {
        phases[0] = ipOil;
        phases[1] = ipWater;
      }
      if( 0 <= ipGas && 0 <= ipOil )
      {
        phases[0] = ipGas;
        phases[1] = ipOil;
      }
    }
    if constexpr ( numPhases == 3  )
    {
      phases[0] = ipGas;
      phases[1] = ipOil;
      phases[2] = ipWater;
    }

    localIndex const numPoints = pressureValues.size( 0 );

    forAll< parallelDevicePolicy<> >( targetSet.size(), [targetSet,
                                                         elementCenter,
                                                         capPressureWrapper,
                                                         elevationIndexTable,
                                                         phases,
                                                         numPoints,
                                                         jFuncMultiplier,
                                                         pressureValues,
                                                         phaseDensityValues,
                                                         phaseComponentFractions,
                                                         globalComponentFractions] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      real64 const elevation = elementCenter[k][2];

      real64 ea = elevationIndexTable.compute( &elevation );
      integer const en = LvArray::math::min( static_cast< integer >(ea), numPoints - 2 );
      ea -= en;

      StackArray< real64, 3, numPhases, geos::constitutive::cappres::LAYOUT_CAPPRES > targetPhaseCapPressure( 1, 1, numPhases );
      StackArray< real64, 2, numPhases, compflow::LAYOUT_PHASE > targetPhaseVolumeFraction( 1, numPhases );
      calculateCapillaryPressure< numPhases >( ea,
                                               pressureValues[en][0],
                                               pressureValues[en+1][0],
                                               phases,
                                               targetPhaseCapPressure[0][0] );


      real64 constexpr initialPhaseVolumeFractionGuess = 1.0 / numPhases;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        targetPhaseVolumeFraction[0][ip] = initialPhaseVolumeFractionGuess;
      }

      localIndex const jFunctionIndex = isJFunction ? k : 0;
      capPressureWrapper.compute( targetPhaseCapPressure[0][0],
                                  jFuncMultiplier[jFunctionIndex],
                                  targetPhaseVolumeFraction[0] );

      for( integer ic = 0; ic < numComps; ++ic )
      {
        globalComponentFractions[k][ic] = 0.0;
      }
      real64 totalMass = 0.0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        real64 const density0 = phaseDensityValues[en][0][ip];
        real64 const density1 = phaseDensityValues[en+1][0][ip];
        real64 const phaseMass = ((1.0-ea)*density0 + ea*density1) * targetPhaseVolumeFraction[0][ip];
        for( integer ic = 0; ic < numComps; ++ic )
        {
          real64 const fraction0 = phaseComponentFractions[en][0][ip][ic];
          real64 const fraction1 = phaseComponentFractions[en+1][0][ip][ic];
          real64 const componentMass = phaseMass * ((1.0-ea)*fraction0 + ea*fraction1);
          globalComponentFractions[k][ic] += componentMass;
          totalMass += componentMass;
        }
      }
      for( integer ic = 0; ic < numComps; ++ic )
      {
        globalComponentFractions[k][ic] /= totalMass;
      }
    } );
  }

  template< integer numPhases >
  static void GEOS_HOST_DEVICE calculateCapillaryPressure( real64 const alpha,
                                                           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE-2 > const & phasePressures0,
                                                           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE-2 > const & phasePressures1,
                                                           integer const phaseOrder[numPhases],
                                                           arraySlice1d< real64, geos::constitutive::cappres::USD_CAPPRES-2 > const & phaseCapPressure )
  {
    if constexpr (numPhases == 2)
    {
      real64 const capPressure0 = phasePressures0[phaseOrder[0]] - phasePressures0[phaseOrder[1]];
      real64 const capPressure1 = phasePressures1[phaseOrder[0]] - phasePressures1[phaseOrder[1]];
      phaseCapPressure[phaseOrder[1]] = (1.0-alpha)*capPressure0 + alpha*capPressure1;
    }
    if constexpr (numPhases == 3)
    {
      real64 const capPressure00 = phasePressures0[phaseOrder[1]] - phasePressures0[phaseOrder[0]];
      real64 const capPressure01 = phasePressures1[phaseOrder[1]] - phasePressures1[phaseOrder[0]];
      phaseCapPressure[phaseOrder[0]] = (1.0-alpha)*capPressure00 + alpha*capPressure01;
      real64 const capPressure10 = phasePressures0[phaseOrder[1]] - phasePressures0[phaseOrder[2]];
      real64 const capPressure11 = phasePressures1[phaseOrder[1]] - phasePressures1[phaseOrder[2]];
      phaseCapPressure[phaseOrder[2]] = (1.0-alpha)*capPressure10 + alpha*capPressure11;
    }
  }
};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_KERNELS_COMPOSITIONAL_CAPILLARYPRESSUREINVERSIONKERNEL_HPP
