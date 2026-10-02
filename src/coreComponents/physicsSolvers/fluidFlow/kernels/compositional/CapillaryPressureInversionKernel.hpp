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
#include "HydrostaticMobility.hpp"
#include "physicsSolvers/fluidFlow/kernels/HydrostaticCoordinate.hpp"

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
                      HydrostaticCoordinate const & coordinate,
                      arrayView1d< real64 const > const & phaseContacts,
                      TableFunction::KernelWrapper elevationIndexTable,
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & pressureValues,
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDensityValues,
                      arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > const & phaseComponentFractions,
                      arrayView2d< real64, compflow::USD_COMP > const globalComponentFractions,
                      arrayView1d< real64 > const pressure,
                      HydrostaticMobility const & evaluateRelativePermeability,
                      real64 const equilTolerance )
  {
    // Actual RP evaluation is a host callback, not a device-callable surrogate.
    targetSet.move( hostMemorySpace, false );
    elementCenter.move( hostMemorySpace, false );
    phaseOrder.move( hostMemorySpace, false );
    phaseContacts.move( hostMemorySpace, false );
    elevationIndexTable.move( hostMemorySpace, false );
    pressureValues.move( hostMemorySpace, false );
    phaseDensityValues.move( hostMemorySpace, false );
    phaseComponentFractions.move( hostMemorySpace, false );
    globalComponentFractions.move( hostMemorySpace, true );
    pressure.move( hostMemorySpace, true );
    if constexpr ( !std::is_same_v< CAP_PRESSURE, constitutive::NoOpCapillaryPressure > )
    {
      GEOS_THROW_IF( !evaluateRelativePermeability, "Missing actual relative-permeability evaluator for hydrostatic capillary initialization", InputError );
      capPressure.forWrappers( []( dataRepository::WrapperBase const & wrapper ) { wrapper.move( hostMemorySpace, false ); } );
    }
    constitutive::InverseCapillaryPressure< CAP_PRESSURE > inverseCapPressureType( capPressure );
    auto capPressureWrapper = inverseCapPressureType.createKernelWrapper();
    auto forwardCapPressureWrapper = capPressure.createKernelWrapper();
    RAJA::ReduceMax< serialReduce, integer > initializationFailure( 0 );

    array2d< real64 > jFuncMultiplierArray( 1, numPhases - 1 );
    arrayView2d< real64 const > jFuncMultiplier = jFuncMultiplierArray.toViewConst();
    constexpr bool isJFunction = std::is_same_v< CAP_PRESSURE, constitutive::JFunctionCapillaryPressure >;
    if constexpr ( isJFunction )
    {
      jFuncMultiplier = capPressure.template getField< geos::fields::cappres::jFuncMultiplier >().reference().toViewConst();
    }

    jFuncMultiplier.move( hostMemorySpace, false );
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

    forAll< serialPolicy >( targetSet.size(), [targetSet,
                                                         elementCenter,
                                                         coordinate,
                                                         phaseContacts,
                                                         capPressureWrapper,
                                                         forwardCapPressureWrapper,
                                                         initializationFailure,
                                                         ipGas,
                                                         ipOil,
                                                         ipWater,
                                                         pressure,
                                                         &evaluateRelativePermeability,
                                                         equilTolerance,
                                                         elevationIndexTable,
                                                         phases,
                                                         numPoints,
                                                         jFuncMultiplier,
                                                         pressureValues,
                                                         phaseDensityValues,
                                                         phaseComponentFractions,
                                                         globalComponentFractions] ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      real64 const elevation = coordinate.elevation( elementCenter[k] );

      real64 ea = elevationIndexTable.compute( &elevation );
      integer const en = LvArray::math::max( 0, LvArray::math::min( static_cast< integer >(ea), numPoints - 2 ) );
      integer const next = LvArray::math::min( en + 1, numPoints - 1 );
      ea = numPoints == 1 ? 0.0 : ea - en;

      StackArray< real64, 3, numPhases, geos::constitutive::cappres::LAYOUT_CAPPRES > targetPhaseCapPressure( 1, 1, numPhases );
      StackArray< real64, 2, numPhases, compflow::LAYOUT_PHASE > targetPhaseVolumeFraction( 1, numPhases );
      calculateCapillaryPressure< numPhases >( ea,
                                               pressureValues[en][0],
                                               pressureValues[next][0],
                                               phases,
                                               targetPhaseCapPressure[0][0] );
      if constexpr ( numPhases == 2 )
      {
        if( ( coordinate.gravityAligned || !std::is_same_v< CAP_PRESSURE, constitutive::NoOpCapillaryPressure > ) &&
            ipWater < 0 && ipGas >= 0 && ipOil >= 0 )
        {
          // Oil is the primary-pressure phase for gas/oil: p_g = P - Pc_g.
          targetPhaseCapPressure[0][0][ipGas] = -targetPhaseCapPressure[0][0][ipOil];
          targetPhaseCapPressure[0][0][ipOil] = 0.0;
        }
      }


      real64 constexpr initialPhaseVolumeFractionGuess = 1.0 / numPhases;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        targetPhaseVolumeFraction[0][ip] = initialPhaseVolumeFractionGuess;
      }

      localIndex const jFunctionIndex = isJFunction ? k : 0;
      if constexpr ( std::is_same_v< CAP_PRESSURE, constitutive::NoOpCapillaryPressure > )
      {
        if( coordinate.gravityAligned )
        {
          // Without capillarity the authored potential-distance contacts define
          // a sharp partition, including the zero-gravity degenerate case.
          integer selected = phases[numPhases-1];
          for( integer contact = 0; contact < phaseContacts.size(); ++contact )
          {
            if( coordinate.isBelow( elevation, phaseContacts[contact], coordinate.projectionErrorBound( elementCenter[k] ) ) ) break;
            selected = phases[numPhases-2-contact];
          }
          for( integer ip = 0; ip < numPhases; ++ip ) targetPhaseVolumeFraction[0][ip] = ip == selected ? 1.0 : 0.0;
        }
        else
        {
          capPressureWrapper.compute( targetPhaseCapPressure[0][0], jFuncMultiplier[jFunctionIndex], targetPhaseVolumeFraction[0] );
        }
      }
      else
      {
        bool const inverseConverged = capPressureWrapper.compute( targetPhaseCapPressure[0][0], jFuncMultiplier[jFunctionIndex], targetPhaseVolumeFraction[0] );
        {
          if( !inverseConverged ) { initializationFailure.max( 1 ); return; }
          real64 constexpr saturationTolerance = 128.0 * LvArray::NumericLimits< real64 >::epsilon;
          real64 saturationSum = 0.0;
          for( integer ip = 0; ip < numPhases; ++ip )
          {
            real64 & saturation = targetPhaseVolumeFraction[0][ip];
            if( !isFinite( saturation ) || saturation < -saturationTolerance || saturation > 1.0 + saturationTolerance )
            { initializationFailure.max( 2 ); return; }
            // Remove only out-of-domain arithmetic roundoff; never snap a
            // small positive phase into absence or a residual-saturation band.
            saturation = LvArray::math::max( 0.0, LvArray::math::min( 1.0, saturation ) );
            saturationSum += saturation;
          }
          if( LvArray::math::abs( saturationSum - 1.0 ) > saturationTolerance )
          { initializationFailure.max( 2 ); return; }
          StackArray< real64, 3, numPhases, geos::constitutive::cappres::LAYOUT_CAPPRES > actualPc( 1, 1, numPhases );
          StackArray< real64, 4, numPhases * numPhases, geos::constitutive::cappres::LAYOUT_CAPPRES_DS > dPc_dS( 1, 1, numPhases, numPhases );
          constitutive::CapillaryPressureEvaluate< CAP_PRESSURE >::compute( forwardCapPressureWrapper,
                                                                            targetPhaseVolumeFraction[0].toSliceConst(),
                                                                            jFuncMultiplier[jFunctionIndex],
                                                                            actualPc[0][0], dPc_dS[0][0] );
          HydrostaticPhaseVector saturationValues{}, relativePermeability{};
          for( integer ip = 0; ip < numPhases; ++ip ) saturationValues[ip] = targetPhaseVolumeFraction[0][ip];
          evaluateRelativePermeability( saturationValues, relativePermeability );
          // The geometric contact phase can be absent after endpoint clipping
          // (for example before a positive entry pressure is reached). Preserve
          // an actually present phase's integrated pressure instead.
          // Use the actual constitutive law: neither a saturation epsilon nor
          // a reported model minimum proves zero mobility for blended laws.
          integer referencePhase = -1;
          for( integer ip = 0; ip < numPhases; ++ip )
          {
            real64 const saturation = targetPhaseVolumeFraction[0][ip];
            if( !( saturation >= -saturationTolerance && saturation <= 1.0 + saturationTolerance ) ||
                !isFinite( actualPc[0][0][ip] ) || !isFinite( relativePermeability[ip] ) || relativePermeability[ip] < 0.0 )
            { initializationFailure.max( 2 ); return; }
            if( saturation > 0.0 && relativePermeability[ip] > 0.0 &&
                ( referencePhase < 0 || saturation > targetPhaseVolumeFraction[0][referencePhase] ) )
            { referencePhase = ip; }
          }
          if( referencePhase < 0 ) { initializationFailure.max( 2 ); return; }
          real64 const referencePressure = (1.0-ea)*pressureValues[en][0][referencePhase] + ea*pressureValues[next][0][referencePhase];
          real64 const primaryPressure = referencePressure + actualPc[0][0][referencePhase];
          if( !isFinite( primaryPressure ) || primaryPressure < 0.0 )
          { initializationFailure.max( 2 ); return; }
          for( integer ip = 0; ip < numPhases; ++ip )
          {
            bool const mobile = targetPhaseVolumeFraction[0][ip] > 0.0 && relativePermeability[ip] > 0.0;
            real64 const phasePressure = (1.0-ea)*pressureValues[en][0][ip] + ea*pressureValues[next][0][ip];
            real64 const residual = (phasePressure-referencePressure) + (actualPc[0][0][ip]-actualPc[0][0][referencePhase]);
            // Subtracting O(P) phase pressures incurs pressure-scaled roundoff.
            // No saturation-based pressure error allowance is substituted.
            real64 const scale = LvArray::math::max( 1.0, LvArray::math::max( LvArray::math::abs( phasePressure ), LvArray::math::abs( referencePressure ) ) );
            real64 const tolerance = equilTolerance + 128.0 * LvArray::NumericLimits< real64 >::epsilon * scale;
            if( mobile && !( LvArray::math::abs( residual ) <= tolerance ) )
            { initializationFailure.max( 3 ); return; }
            // An excluded continued phase must not have excess pressure that
            // would require its appearance beyond a capillary endpoint.
            if( !mobile && !( residual <= tolerance ) )
            { initializationFailure.max( 5 ); return; }
          }
          pressure[k] = primaryPressure;
        }
      }

      for( integer ic = 0; ic < numComps; ++ic )
      {
        globalComponentFractions[k][ic] = 0.0;
      }
      real64 totalMass = 0.0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        real64 const density0 = phaseDensityValues[en][0][ip];
        real64 const density1 = phaseDensityValues[next][0][ip];
        real64 const phaseMass = ((1.0-ea)*density0 + ea*density1) * targetPhaseVolumeFraction[0][ip];
        for( integer ic = 0; ic < numComps; ++ic )
        {
          real64 const fraction0 = phaseComponentFractions[en][0][ip][ic];
          real64 const fraction1 = phaseComponentFractions[next][0][ip][ic];
          real64 const componentMass = phaseMass * ((1.0-ea)*fraction0 + ea*fraction1);
          globalComponentFractions[k][ic] += componentMass;
          totalMass += componentMass;
        }
      }
      if( ( coordinate.gravityAligned || !std::is_same_v< CAP_PRESSURE, constitutive::NoOpCapillaryPressure > ) &&
          ( !isFinite( totalMass ) || totalMass <= 0.0 ) )
      { initializationFailure.max( 4 ); return; }
      for( integer ic = 0; ic < numComps; ++ic )
      {
        globalComponentFractions[k][ic] /= totalMass;
      }
    } );
    GEOS_THROW_IF( initializationFailure.get() != 0,
                   GEOS_FMT( "Hydrostatic capillary initialization failed: {}. "
                             "The active capillary law must reproduce the integrated pressure differences for every phase with positive actual relative permeability and saturation; "
                             "check phase roles, capillary endpoints and equilibrationTolerance.",
                             initializationFailure.get() == 1 ? "capillary inversion did not converge" :
                             initializationFailure.get() == 2 ? "invalid saturation or reference pressure" :
                             initializationFailure.get() == 3 ? "inconsistent mobile-phase capillary pressure" :
                             initializationFailure.get() == 5 ? "inconsistent capillary endpoint complementarity" : "invalid reconstructed component mass" ),
                   InputError );
  }

  GEOS_HOST_DEVICE
  static bool isFinite( real64 const value )
  {
    return LvArray::math::abs( value ) <= LvArray::NumericLimits< real64 >::max;
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
