/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseWellUpdateVolRatesForConstraint.cpp
 */

#include "CompositionalMultiphaseWell.hpp"

#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace constitutive;

void CompositionalMultiphaseWell::updateVolRatesForConstraint( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  integer constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const & temp = subRegion.getField< fields::well::temperature >();
  arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::mixtureConnectionRate >();

  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();
  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens = subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >();

  // fluid data

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac = fluid.phaseFraction();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac = fluid.dPhaseFraction();

  arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens = fluid.totalDensity();
  arrayView3d< real64 const, multifluid::USD_FLUID_DC > const & dTotalDens = fluid.dTotalDensity();

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens = fluid.phaseDensity();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens = fluid.dPhaseDensity();

  // control data

  WellControls & wellControls = getWellControls( subRegion );
  string const wellControlsName = wellControls.getName();
  integer const logLevel = wellControls.getLogLevel();
  string const massUnit = m_useMass ? "kg" : "mol";

  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 const & surfacePres = wellControls.getSurfacePressure();
  real64 const & surfaceTemp = wellControls.getSurfaceTemperature();

  arrayView1d< real64 > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dPres =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dPresString() );
  arrayView2d< real64 > const & dCurrentPhaseVolRate_dCompDens =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dRateString() );

  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
  real64 & dCurrentTotalVolRate_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dPresString() );
  arrayView1d< real64 > const & dCurrentTotalVolRate_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dCompDensString() );
  real64 & dCurrentTotalVolRate_dRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dRateString() );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [&numComp,
                                &numPhase,
                                fluidWrapper,
                                pres,
                                temp,
                                compFrac,
                                dCompFrac_dCompDens,
                                connRate,
                                totalDens,
                                dTotalDens,
                                phaseDens,
                                dPhaseDens,
                                phaseFrac,
                                dPhaseFrac,
                                &useSurfaceConditions,
                                &surfacePres,
                                &surfaceTemp,
                                &currentTotalVolRate,
                                &dCurrentTotalVolRate_dPres,
                                dCurrentTotalVolRate_dCompDens,
                                &dCurrentTotalVolRate_dRate,
                                currentPhaseVolRate,
                                dCurrentPhaseVolRate_dPres,
                                dCurrentPhaseVolRate_dCompDens,
                                dCurrentPhaseVolRate_dRate,
                                &iwelemRef,
                                &logLevel,
                                &wellControlsName,
                                &massUnit] ( localIndex const )
    {
      GEOS_UNUSED_VAR( massUnit );
      using Deriv = multifluid::DerivativeOffset;

      stackArray1d< real64, maxNumComp > work( numComp );

      // Step 1: evaluate the phase and total density in the reference element

      //    We need to evaluate the density as follows:
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element
      if( useSurfaceConditions )
      {
        // we need to compute the surface density
        fluidWrapper.update( iwelemRef, 0, surfacePres, surfaceTemp, compFrac[iwelemRef] );
        if( logLevel >= 2 )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: surface density computed with P_surface = {} Pa and T_surface = {} K",
                                   wellControlsName, surfacePres, surfaceTemp ) );
#ifdef GEOS_USE_HIP
          GEOS_UNUSED_VAR( wellControlsName );
#endif
        }
      }
      else
      {
        real64 const refPres = pres[iwelemRef];
        fluidWrapper.update( iwelemRef, 0, refPres, temp[iwelemRef], compFrac[iwelemRef] );
      }

      // Step 2: update the total volume rate

      real64 const currentTotalRate = connRate[iwelemRef];

      // Step 2.1: compute the inverse of the total density and derivatives

      real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];
      real64 const dTotalDensInv_dPres = -dTotalDens[iwelemRef][0][Deriv::dP] * totalDensInv * totalDensInv;
      stackArray1d< real64, maxNumComp > dTotalDensInv_dCompDens( numComp );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalDensInv_dCompDens[ic] = -dTotalDens[iwelemRef][0][Deriv::dC+ic] * totalDensInv * totalDensInv;
      }
      applyChainRuleInPlace( numComp, dCompFrac_dCompDens[iwelemRef], dTotalDensInv_dCompDens, work.data() );

      // Step 2.2: divide the total mass/molar rate by the total density to get the total volumetric rate
      currentTotalVolRate = currentTotalRate * totalDensInv;
      dCurrentTotalVolRate_dPres = ( useSurfaceConditions ==  0 ) * currentTotalRate * dTotalDensInv_dPres;
      dCurrentTotalVolRate_dRate = totalDensInv;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dCurrentTotalVolRate_dCompDens[ic] = currentTotalRate * dTotalDensInv_dCompDens[ic];
      }

      if( logLevel >= 2 && useSurfaceConditions )
      {
        GEOS_LOG_RANK( GEOS_FMT( "{}: The total fluid density at surface conditions is {} {}/sm3. \n"
                                 "The total rate is {} {}/s, which corresponds to a total surface volumetric rate of {} sm3/s",
                                 wellControlsName, totalDens[iwelemRef][0], massUnit,
                                 currentTotalRate, massUnit, currentTotalVolRate ) );
      }

      // Step 3: update the phase volume rate
      for( integer ip = 0; ip < numPhase; ++ip )
      {

        // Step 3.1: compute the inverse of the (phase density * phase fraction) and derivatives

        // skip the rest of this function if phase ip is absent
        bool const phaseExists = (phaseFrac[iwelemRef][0][ip] > 0);
        if( !phaseExists )
        {
          continue;
        }

        real64 const phaseDensInv =  1.0 / phaseDens[iwelemRef][0][ip];
        real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;
        real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac[iwelemRef][0][ip][Deriv::dP] * phaseDensInv
                                                         - dPhaseDens[iwelemRef][0][ip][Deriv::dP] * phaseFracTimesPhaseDensInv * phaseDensInv;


        // Step 3.2: divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
        currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
        dCurrentPhaseVolRate_dPres[ip] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
        dCurrentPhaseVolRate_dRate[ip] = phaseFracTimesPhaseDensInv;
        for( integer ic = 0; ic < numComp; ++ic )
        {
          dCurrentPhaseVolRate_dCompDens[ip][ic] = -phaseFracTimesPhaseDensInv * dPhaseDens[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] += dPhaseFrac[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] *= currentTotalRate;
        }
        applyChainRuleInPlace( numComp, dCompFrac_dCompDens[iwelemRef], dCurrentPhaseVolRate_dCompDens[ip], work.data() );

        if( logLevel >= 2 && useSurfaceConditions )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: The density of phase {} at surface conditions is {} {}/sm3. \n"
                                   "The phase surface volumetric rate is {} sm3/s",
                                   wellControlsName, ip, phaseDens[iwelemRef][0][ip], massUnit, currentPhaseVolRate[ip] ) );
        }
      }
    } );
  } );
}

} // namespace geos
