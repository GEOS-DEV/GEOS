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

// Source includes
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "functions/FunctionManager.hpp"
#include "mainInterface/GeosxState.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::constitutive::multifluid;
using namespace geosx::dataRepository;
using namespace geosx::constitutive::PVTProps;

/// Black-oil tables written into temporary files during testing

static const char * pvtoTableContent = "# Rs[sm3/sm3]\tPbub[Pa]\tBo[m3/sm3]\tVisc(Pa.s)\n"
                                       "\n"
                                       "  2\t            2000000\t    1.02\t    0.000975\n"
                                       "  5\t            5000000\t    1.03\t    0.00091\n"
                                       " 10\t            10000000\t1.04\t    0.00083\n"
                                       " 15\t            20000000\t1.05\t    0.000695\n"
                                       "                90000000\t1.03\t    0.000985  -- some line comment\n"
                                       " 30\t            30000000\t1.07\t    0.000594\n"
                                       " 40\t            40000000\t1.08\t    0.00051\n"
                                       "                50000000\t1.07\t    0.000549  -- another one\n"
                                       "                90000000\t1.06\t    0.00074\n"
                                       " 50\t            50000000.7\t1.09\t    0.000449\n"
                                       "                90000000.7\t1.08\t    0.000605";

static const char * pvtwTableContent = "#\tPref[bar]\tBw[m3/sm3]\tCp[1/bar]\t    Visc[cP]\n"
                                       "\t30600000.1\t1.03\t\t0.00000000041\t0.0003";

/// Dead-oil tables written into temporary files during testing

static const char * pvdgTableContent = "# Pg(Pa) Bg(m3/sm3) Visc(Pa.s)\n"
                                       "3000000  0.04234  0.00001344\n"
                                       "6000000  0.02046  0.0000142\n"
                                       "9000000  0.01328  0.00001526\n"
                                       "12000000 0.00977  0.0000166\n"
                                       "15000000 0.00773  0.00001818\n"
                                       "18000000 0.006426 0.00001994\n"
                                       "21000000 0.005541 0.00002181\n"
                                       "24000000 0.004919 0.0000237\n"
                                       "27000000 0.004471 0.00002559\n"
                                       "29500000 0.004194 0.00002714\n"
                                       "31000000 0.004031 0.00002806\n"
                                       "33000000 0.00391  0.00002832\n"
                                       "53000000 0.003868 0.00002935";

static const char * pvdoTableContent = "#P[Pa] Bo[m3/sm3] Visc(Pa.s)\n"
                                       "2000000  1.02 0.000975\n"
                                       "5000000  1.03 0.00091\n"
                                       "10000000 1.04 0.00083\n"
                                       "20000000 1.05 0.000695\n"
                                       "30000000 1.07 0.000594\n"
                                       "40000000 1.08 0.00051\n"
                                       "50000000.7 1.09 0.000449";

static const char * pvdwTableContent = "# Pref[bar] Bw[m3/sm3] Cp[1/bar]     Visc[cP]\n"
                                       " 30600000.1 1.03  0.00000000041 0.0003";

// CO2-brine model

static const char * pvtLiquidPhillipsTableContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                                    "ViscosityFun PhillipsBrineViscosity 0.1";

// the last are set relatively high (1e-4) to increase derivative value and check it properly
static const char * pvtLiquidEzrokhiTableContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                                   "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4";

static const char * pvtGasTableContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                         "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1";

static const char * co2FlashTableContent = "FlashModel CO2Solubility 1e6 1.5e7 5e4 367.15 369.15 1 0.15";

void testNumericalDerivatives( MultiFluidBase & fluid,
                               Group & parent,
                               real64 const P,
                               real64 const T,
                               arraySlice1d< real64 > const & compositionInput,
                               real64 const perturbParameter,
                               bool usePVTPackage,
                               real64 const relTol,
                               real64 const absTol = std::numeric_limits< real64 >::max() )
{
  localIndex const NC = fluid.numFluidComponents();
  localIndex const NP = fluid.numFluidPhases();

  // Copy input values into an array with expected layout
  array2d< real64, compflow::LAYOUT_COMP > compositionValues( 1, NC );
  for( localIndex i = 0; i < NC; ++i )
  {
    compositionValues[0][i] = compositionInput[i];
  }
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition = compositionValues[0];

  auto const & components = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  auto const & phases     = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );

  // create a clone of the fluid to run updates on
  std::unique_ptr< ConstitutiveBase > fluidCopyPtr = fluid.deliverClone( "fluidCopy", &parent );
  MultiFluidBase & fluidCopy = dynamicCast< MultiFluidBase & >( *fluidCopyPtr );

  fluid.allocateConstitutiveData( fluid.getParent(), 1 );
  fluidCopy.allocateConstitutiveData( fluid.getParent(), 1 );

  // extract data views from both fluids
  #define GET_FLUID_DATA( FLUID, DIM, PERM, KEY ) \
    FLUID.getReference< Array< real64, DIM, PERM > >( MultiFluidBase::viewKeyStruct::KEY() )[0][0]

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseFrac {
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, phaseFractionString ),
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, dPhaseFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, dPhaseFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, LAYOUT_PHASE_DC, dPhaseFraction_dGlobalCompFractionString )
  };

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseDens {
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, phaseDensityString ),
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, dPhaseDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, dPhaseDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, LAYOUT_PHASE_DC, dPhaseDensity_dGlobalCompFractionString )
  };

  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 > phaseVisc {
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, phaseViscosityString ),
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, dPhaseViscosity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, LAYOUT_PHASE, dPhaseViscosity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, LAYOUT_PHASE_DC, dPhaseViscosity_dGlobalCompFractionString )
  };

  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 > phaseCompFrac {
    GET_FLUID_DATA( fluid, 4, LAYOUT_PHASE_COMP, phaseCompFractionString ),
    GET_FLUID_DATA( fluid, 4, LAYOUT_PHASE_COMP, dPhaseCompFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 4, LAYOUT_PHASE_COMP, dPhaseCompFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 5, LAYOUT_PHASE_COMP_DC, dPhaseCompFraction_dGlobalCompFractionString )
  };

  MultiFluidVarSlice< real64, 0, USD_FLUID - 2, USD_FLUID_DC - 2 > totalDens {
    GET_FLUID_DATA( fluid, 2, LAYOUT_FLUID, totalDensityString ),
    GET_FLUID_DATA( fluid, 2, LAYOUT_FLUID, dTotalDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 2, LAYOUT_FLUID, dTotalDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 3, LAYOUT_FLUID_DC, dTotalDensity_dGlobalCompFractionString )
  };

  auto const & phaseFracCopy     = GET_FLUID_DATA( fluidCopy, 3, LAYOUT_PHASE, phaseFractionString );
  auto const & phaseDensCopy     = GET_FLUID_DATA( fluidCopy, 3, LAYOUT_PHASE, phaseDensityString );
  auto const & phaseViscCopy     = GET_FLUID_DATA( fluidCopy, 3, LAYOUT_PHASE, phaseViscosityString );
  auto const & phaseCompFracCopy = GET_FLUID_DATA( fluidCopy, 4, LAYOUT_PHASE_COMP, phaseCompFractionString );
  auto const & totalDensCopy     = GET_FLUID_DATA( fluidCopy, 2, LAYOUT_FLUID, totalDensityString );

#undef GET_FLUID_DATA

  // set the original fluid state to current
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    fluidWrapper.update( 0, 0, P, T, composition );
  } );

  // now perturb variables and update the copied fluid's state
  constitutive::constitutiveUpdatePassThru( fluidCopy, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // update pressure and check derivatives
    {
      real64 const dP = perturbParameter * (P + perturbParameter);
      fluidWrapper.update( 0, 0, P + dP, T, composition );

      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), phaseFrac.dPres.toSliceConst(),
                       dP, relTol, absTol, "phaseFrac", "Pres", phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), phaseDens.dPres.toSliceConst(),
                       dP, relTol, absTol, "phaseDens", "Pres", phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), phaseVisc.dPres.toSliceConst(),
                       dP, relTol, absTol, "phaseVisc", "Pres", phases );
      checkDerivative( totalDensCopy, totalDens.value, totalDens.dPres,
                       dP, relTol, absTol, "totalDens", "Pres" );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), phaseCompFrac.dPres.toSliceConst(),
                       dP, relTol, absTol, "phaseCompFrac", "Pres", phases, components );
    }

    // update temperature and check derivatives
    {
      real64 const dT = perturbParameter * (T + perturbParameter);
      fluidWrapper.update( 0, 0, P, T + dT, composition );

      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), phaseFrac.dTemp.toSliceConst(),
                       dT, relTol, absTol, "phaseFrac", "Temp", phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), phaseDens.dTemp.toSliceConst(),
                       dT, relTol, absTol, "phaseDens", "Temp", phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), phaseVisc.dTemp.toSliceConst(),
                       dT, relTol, absTol, "phaseVisc", "Temp", phases );
      checkDerivative( totalDensCopy, totalDens.value, totalDens.dTemp,
                       dT, relTol, absTol, "totalDens", "Temp" );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), phaseCompFrac.dTemp.toSliceConst(),
                       dT, relTol, absTol, "phaseCompFrac", "Temp", phases, components );
    }

    // update composition and check derivatives
    auto dPhaseFrac_dC     = invertLayout( phaseFrac.dComp.toSliceConst(), NP, NC );
    auto dPhaseDens_dC     = invertLayout( phaseDens.dComp.toSliceConst(), NP, NC );
    auto dPhaseVisc_dC     = invertLayout( phaseVisc.dComp.toSliceConst(), NP, NC );
    auto dTotalDens_dC     = invertLayout( totalDens.dComp.toSliceConst(), NC );
    auto dPhaseCompFrac_dC = invertLayout( phaseCompFrac.dComp.toSliceConst(), NP, NC, NC );

    array2d< real64, compflow::LAYOUT_COMP > compNew( 1, NC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      real64 const dC = perturbParameter * ( composition[jc] + perturbParameter );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compNew[0][ic] = composition[ic];
      }
      compNew[0][jc] += dC;

      // Note: in PVTPackage, derivatives are obtained with finite-difference approx **with normalization of the comp fraction**
      //       The component fraction is perturbed (just as above), and then all the component fractions are normalized (as below)
      //       But, in the native DO model and in CO2BrinePhillips, derivatives are computed analytically, which results in different
      //       derivatives wrt component fractions--although the derivatives wrt component densities obtained with the chain rule
      //       in the solver will be very similar (see discussion on PR #1325 on GitHub).
      //
      //       Since both approaches--FD approximation of derivatives with normalization, and analytical derivatives--are correct,
      //       we have to support both when we check the intermediate derivatives wrt component fractions below. Therefore, if the
      //       PVTPackage is used, then we normalize the perturbed component fractions before taking the FD approx. If the native
      //       DO or CO2-brine models are used, we skip the normalization below.
      if( usePVTPackage )
      {
        // renormalize
        real64 sum = 0.0;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          sum += compNew[0][ic];
        }
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          compNew[0][ic] /= sum;
        }
      }

      fluidWrapper.update( 0, 0, P, T, compNew[0] );

      string const var = "compFrac[" + components[jc] + "]";
      checkDerivative( phaseFracCopy.toSliceConst(), phaseFrac.value.toSliceConst(), dPhaseFrac_dC[jc].toSliceConst(),
                       dC, relTol, absTol, "phaseFrac", var, phases );
      checkDerivative( phaseDensCopy.toSliceConst(), phaseDens.value.toSliceConst(), dPhaseDens_dC[jc].toSliceConst(),
                       dC, relTol, absTol, "phaseDens", var, phases );
      checkDerivative( phaseViscCopy.toSliceConst(), phaseVisc.value.toSliceConst(), dPhaseVisc_dC[jc].toSliceConst(),
                       dC, relTol, absTol, "phaseVisc", var, phases );
      checkDerivative( totalDensCopy, totalDens.value, dTotalDens_dC[jc],
                       dC, relTol, absTol, "totalDens", var );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), dPhaseCompFrac_dC[jc].toSliceConst(),
                       dC, relTol, absTol, "phaseCompFrac", var, phases, components );
    }
  } );
}

void testValuesAgainstPreviousImplementation( CO2BrinePhillipsFluid::KernelWrapper const & wrapper,
                                              real64 const P,
                                              real64 const T,
                                              arraySlice1d< real64 const > const & compositionInput,
                                              real64 const & savedTotalDensity,
                                              real64 const & savedGasPhaseFrac,
                                              real64 const & savedWaterDens,
                                              real64 const & savedGasDens,
                                              real64 const & savedWaterMassDens,
                                              real64 const & savedGasMassDens,
                                              real64 const & savedWaterVisc,
                                              real64 const & savedGasVisc,
                                              real64 const & savedWaterPhaseGasComp,
                                              real64 const & savedWaterPhaseWaterComp,
                                              real64 const relTol )
{
  // Copy input values into an array with expected layout
  array2d< real64, compflow::LAYOUT_COMP > compositionValues( 1, 2 );
  for( localIndex i = 0; i < 2; ++i )
  {
    compositionValues[0][i] = compositionInput[i];
  }
  arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition = compositionValues[0];

  StackArray< real64, 3, 2, LAYOUT_PHASE > phaseFraction( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseFraction_dPressure( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseFraction_dTemperature( 1, 1, 2 );
  StackArray< real64, 4, 4, LAYOUT_PHASE_DC > dPhaseFraction_dGlobalCompFraction( 1, 1, 2, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > phaseDensity( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseDensity_dPressure( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseDensity_dTemperature( 1, 1, 2 );
  StackArray< real64, 4, 4, LAYOUT_PHASE_DC > dPhaseDensity_dGlobalCompFraction( 1, 1, 2, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > phaseMassDensity( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseMassDensity_dPressure( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseMassDensity_dTemperature( 1, 1, 2 );
  StackArray< real64, 4, 4, LAYOUT_PHASE_DC > dPhaseMassDensity_dGlobalCompFraction( 1, 1, 2, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > phaseViscosity( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseViscosity_dPressure( 1, 1, 2 );
  StackArray< real64, 3, 2, LAYOUT_PHASE > dPhaseViscosity_dTemperature( 1, 1, 2 );
  StackArray< real64, 4, 4, LAYOUT_PHASE_DC > dPhaseViscosity_dGlobalCompFraction( 1, 1, 2, 2 );
  StackArray< real64, 4, 4, LAYOUT_PHASE_COMP > phaseCompFraction( 1, 1, 2, 2 );
  StackArray< real64, 4, 4, LAYOUT_PHASE_COMP > dPhaseCompFraction_dPressure( 1, 1, 2, 2 );
  StackArray< real64, 4, 4, LAYOUT_PHASE_COMP > dPhaseCompFraction_dTemperature( 1, 1, 2, 2 );
  StackArray< real64, 5, 8, LAYOUT_PHASE_COMP_DC > dPhaseCompFraction_dGlobalCompFraction( 1, 1, 2, 2, 2 );
  StackArray< real64, 2, 1, LAYOUT_FLUID > totalDensity( 1, 1 );
  StackArray< real64, 2, 1, LAYOUT_FLUID > dTotalDensity_dPressure( 1, 1 );
  StackArray< real64, 2, 1, LAYOUT_FLUID > dTotalDensity_dTemperature( 1, 1 );
  StackArray< real64, 3, 2, LAYOUT_FLUID_DC >  dTotalDensity_dGlobalCompFraction( 1, 1, 2 );

  wrapper.compute( P, T, composition,
  {
    phaseFraction[0][0],
    dPhaseFraction_dPressure[0][0],
    dPhaseFraction_dTemperature[0][0],
    dPhaseFraction_dGlobalCompFraction[0][0]
  },
  {
    phaseDensity[0][0],
    dPhaseDensity_dPressure[0][0],
    dPhaseDensity_dTemperature[0][0],
    dPhaseDensity_dGlobalCompFraction[0][0]
  },
  {
    phaseMassDensity[0][0],
    dPhaseMassDensity_dPressure[0][0],
    dPhaseMassDensity_dTemperature[0][0],
    dPhaseMassDensity_dGlobalCompFraction[0][0]
  },
  {
    phaseViscosity[0][0],
    dPhaseViscosity_dPressure[0][0],
    dPhaseViscosity_dTemperature[0][0],
    dPhaseViscosity_dGlobalCompFraction[0][0]
  },
  {
    phaseCompFraction[0][0],
    dPhaseCompFraction_dPressure[0][0],
    dPhaseCompFraction_dTemperature[0][0],
    dPhaseCompFraction_dGlobalCompFraction[0][0]
  },
  {
    totalDensity[0][0],
    dTotalDensity_dPressure[0][0],
    dTotalDensity_dTemperature[0][0],
    dTotalDensity_dGlobalCompFraction[0][0]
  } );

  checkRelativeError( totalDensity[0][0], savedTotalDensity, relTol );
  for( localIndex ip = 0; ip < 2; ++ip )
  {
    real64 const savedPhaseFrac = ( ip == 0 ) ? savedGasPhaseFrac : 1 - savedGasPhaseFrac;
    checkRelativeError( phaseFraction[0][0][ip], savedPhaseFrac, relTol );
    real64 const savedPhaseDens = ( ip == 0 ) ? savedGasDens : savedWaterDens;
    checkRelativeError( phaseDensity[0][0][ip], savedPhaseDens, relTol );
    real64 const savedPhaseMassDens = ( ip == 0 ) ? savedGasMassDens : savedWaterMassDens;
    checkRelativeError( phaseMassDensity[0][0][ip], savedPhaseMassDens, relTol );
    real64 const savedPhaseVisc = ( ip == 0 ) ? savedGasVisc : savedWaterVisc;
    checkRelativeError( phaseViscosity[0][0][ip], savedPhaseVisc, relTol );
    for( localIndex ic = 0; ic < 2; ++ic )
    {
      real64 savedCompFrac = 0.0;
      if( ip == 0 )
      {
        savedCompFrac = ( ic == 0 ) ? 1 : 0;
      }
      else
      {
        savedCompFrac = ( ic == 0 ) ? savedWaterPhaseGasComp : savedWaterPhaseWaterComp;
      }
      checkRelativeError( phaseCompFraction[0][0][ip][ic], savedCompFrac, relTol );
    }
  }
}

MultiFluidBase & makeCompositionalFluid( string const & name, Group & parent )
{
  CompositionalMultiphaseFluid & fluid = parent.registerGroup< CompositionalMultiphaseFluid >( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 4 );
  compNames[0] = "N2"; compNames[1] = "C10"; compNames[2] = "C20"; compNames[3] = "H20";

  auto & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 4 );
  molarWgt[0] = 28e-3; molarWgt[1] = 134e-3; molarWgt[2] = 275e-3; molarWgt[3] = 18e-3;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas";

  auto & eqnOfState = fluid.getReference< string_array >( CompositionalMultiphaseFluid::viewKeyStruct::equationsOfStateString() );
  eqnOfState.resize( 2 );
  eqnOfState[0] = "PR"; eqnOfState[1] = "PR";

  auto & critPres = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalPressureString() );
  critPres.resize( 4 );
  critPres[0] = 34e5; critPres[1] = 25.3e5; critPres[2] = 14.6e5; critPres[3] = 220.5e5;

  auto & critTemp = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalTemperatureString() );
  critTemp.resize( 4 );
  critTemp[0] = 126.2; critTemp[1] = 622.0; critTemp[2] = 782.0; critTemp[3] = 647.0;

  auto & acFactor = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluid::viewKeyStruct::componentAcentricFactorString() );
  acFactor.resize( 4 );
  acFactor[0] = 0.04; acFactor[1] = 0.443; acFactor[2] = 0.816; acFactor[3] = 0.344;

  fluid.postProcessInputRecursive();
  return fluid;
}

class CompositionalFluidTestBase : public ::testing::Test
{
public:
  CompositionalFluidTestBase():
    node(),
    parent( "parent", node )
  {}

protected:
  conduit::Node node;
  Group parent;
  MultiFluidBase * fluid;
};

class CompositionalFluidTest : public CompositionalFluidTestBase
{
public:
  CompositionalFluidTest()
  {
    parent.resize( 1 );
    fluid = &makeCompositionalFluid( "fluid", parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }
};

TEST_F( CompositionalFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 4 );
  comp[0] = 0.099; comp[1] = 0.3; comp[2] = 0.6; comp[3] = 0.001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, true, relTol );
}

TEST_F( CompositionalFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 4 );
  comp[0] = 0.099; comp[1] = 0.3; comp[2] = 0.6; comp[3] = 0.001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-2;

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, true, relTol );
}

MultiFluidBase & makeLiveOilFluid( string const & name, Group * parent )
{
  BlackOilFluid & fluid = parent->registerGroup< BlackOilFluid >( name );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 3 );
  compNames[0] = "oil"; compNames[1] = "gas"; compNames[2] = "water";

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 3 );
  molarWgt[0] = 114e-3; molarWgt[1] = 16e-3; molarWgt[2] = 18e-3;

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluidBase::viewKeyStruct::surfacePhaseMassDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0; surfaceDens[1] = 0.9907; surfaceDens[2] = 1022.0;

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluidBase::viewKeyStruct::tableFilesString() );
  tableNames.resize( 3 );
  tableNames[0] = "pvto.txt"; tableNames[1] = "pvdg.txt"; tableNames[2] = "pvtw.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

MultiFluidBase & makeDeadOilFluid( string const & name, Group * parent )
{
  DeadOilFluid & fluid = parent->registerGroup< DeadOilFluid >( name );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 3 );
  compNames[0] = "oil"; compNames[1] = "water"; compNames[2] = "gas";

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 3 );
  molarWgt[0] = 114e-3; molarWgt[1] = 16e-3; molarWgt[2] = 18e-3;

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "water"; phaseNames[2] = "gas";

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluidBase::viewKeyStruct::surfacePhaseMassDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0; surfaceDens[1] = 1022.0; surfaceDens[2] = 0.9907;

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluidBase::viewKeyStruct::tableFilesString() );
  tableNames.resize( 3 );
  tableNames[0] = "pvdo.txt"; tableNames[1] = "pvdw.txt"; tableNames[2] = "pvdg.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

MultiFluidBase & makeDeadOilFluidFromTable( string const & name, Group * parent )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // 1) First, define the tables (PVDO, PVDG)

  // 1D table with linear interpolation
  localIndex const NaxisPVDO = 7;
  localIndex const NaxisPVDG = 13;

  array1d< real64_array > coordinatesPVDO;
  real64_array valuesPVDO_Bo( NaxisPVDO );
  real64_array valuesPVDO_visc( NaxisPVDO );
  coordinatesPVDO.resize( 1 );
  coordinatesPVDO[0].resize( NaxisPVDO );
  coordinatesPVDO[0][0] =  2000000; valuesPVDO_Bo[0] = 1.02; valuesPVDO_visc[0] = 0.000975;
  coordinatesPVDO[0][1] =  5000000; valuesPVDO_Bo[1] = 1.03; valuesPVDO_visc[1] = 0.00091;
  coordinatesPVDO[0][2] = 10000000; valuesPVDO_Bo[2] = 1.04; valuesPVDO_visc[2] = 0.00083;
  coordinatesPVDO[0][3] = 20000000; valuesPVDO_Bo[3] = 1.05; valuesPVDO_visc[3] = 0.000695;
  coordinatesPVDO[0][4] = 30000000; valuesPVDO_Bo[4] = 1.07; valuesPVDO_visc[4] = 0.000594;
  coordinatesPVDO[0][5] = 40000000; valuesPVDO_Bo[5] = 1.08; valuesPVDO_visc[5] = 0.00051;
  coordinatesPVDO[0][6] = 50000000; valuesPVDO_Bo[6] = 1.09; valuesPVDO_visc[6] = 0.000449;

  array1d< real64_array > coordinatesPVDG;
  real64_array valuesPVDG_Bg( NaxisPVDG );
  real64_array valuesPVDG_visc( NaxisPVDG );
  coordinatesPVDG.resize( 1 );
  coordinatesPVDG[0].resize( NaxisPVDG );
  coordinatesPVDG[0][0]  =  3000000; valuesPVDG_Bg[0]  = 0.04234;  valuesPVDG_visc[0] = 0.00001344;
  coordinatesPVDG[0][1]  =  6000000; valuesPVDG_Bg[1]  = 0.02046;  valuesPVDG_visc[1] = 0.0000142;
  coordinatesPVDG[0][2]  =  9000000; valuesPVDG_Bg[2]  = 0.01328;  valuesPVDG_visc[2] = 0.00001526;
  coordinatesPVDG[0][3]  = 12000000; valuesPVDG_Bg[3]  = 0.00977;  valuesPVDG_visc[3] = 0.0000166;
  coordinatesPVDG[0][4]  = 15000000; valuesPVDG_Bg[4]  = 0.00773;  valuesPVDG_visc[4] = 0.00001818;
  coordinatesPVDG[0][5]  = 18000000; valuesPVDG_Bg[5]  = 0.006426; valuesPVDG_visc[5] = 0.00001994;
  coordinatesPVDG[0][6]  = 21000000; valuesPVDG_Bg[6]  = 0.005541; valuesPVDG_visc[6] = 0.00002181;
  coordinatesPVDG[0][7]  = 24000000; valuesPVDG_Bg[7]  = 0.004919; valuesPVDG_visc[7] = 0.0000237;
  coordinatesPVDG[0][8]  = 27000000; valuesPVDG_Bg[8]  = 0.004471; valuesPVDG_visc[8] = 0.00002559;
  coordinatesPVDG[0][9]  = 29500000; valuesPVDG_Bg[9]  = 0.004194; valuesPVDG_visc[9] = 0.00002714;
  coordinatesPVDG[0][10] = 31000000; valuesPVDG_Bg[10] = 0.004031; valuesPVDG_visc[10] = 0.00002806;
  coordinatesPVDG[0][11] = 33000000; valuesPVDG_Bg[11] = 0.00391;  valuesPVDG_visc[11] = 0.00002832;
  coordinatesPVDG[0][12] = 53000000; valuesPVDG_Bg[12] = 0.003868; valuesPVDG_visc[12] = 0.00002935;

  TableFunction & tablePVDO_Bo = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDO_Bo" ) );
  tablePVDO_Bo.setTableCoordinates( coordinatesPVDO );
  tablePVDO_Bo.setTableValues( valuesPVDO_Bo );
  tablePVDO_Bo.reInitializeFunction();
  tablePVDO_Bo.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDO_visc = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDO_visc" ) );
  tablePVDO_visc.setTableCoordinates( coordinatesPVDO );
  tablePVDO_visc.setTableValues( valuesPVDO_visc );
  tablePVDO_visc.reInitializeFunction();
  tablePVDO_visc.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDG_Bg = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDG_Bg" ) );
  tablePVDG_Bg.setTableCoordinates( coordinatesPVDG );
  tablePVDG_Bg.setTableValues( valuesPVDG_Bg );
  tablePVDG_Bg.reInitializeFunction();
  tablePVDG_Bg.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  TableFunction & tablePVDG_visc = dynamicCast< TableFunction & >( *functionManager.createChild( "TableFunction", "PVDG_visc" ) );
  tablePVDG_visc.setTableCoordinates( coordinatesPVDG );
  tablePVDG_visc.setTableValues( valuesPVDG_visc );
  tablePVDG_visc.reInitializeFunction();
  tablePVDG_visc.setInterpolationMethod( TableFunction::InterpolationType::Linear );

  // 2) Then, define the Dead-Oil constitutive model

  DeadOilFluid & fluid = parent->registerGroup< DeadOilFluid >( name );

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 3 );
  compNames[0] = "gas"; compNames[1] = "water"; compNames[2] = "oil";

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 3 );
  molarWgt[0] = 18e-3; molarWgt[1] = 16e-3; molarWgt[2] = 114e-3;

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "gas"; phaseNames[1] = "water"; phaseNames[2] = "oil";

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( DeadOilFluid::viewKeyStruct::surfacePhaseMassDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 0.9907; surfaceDens[1] = 1022.0; surfaceDens[2] = 800.0;

  string_array & FVFTableNames = fluid.getReference< string_array >( DeadOilFluid::viewKeyStruct::formationVolumeFactorTableNamesString() );
  FVFTableNames.resize( 2 );
  FVFTableNames[0] = "PVDG_Bg"; FVFTableNames[1] = "PVDO_Bo";

  string_array & viscosityTableNames = fluid.getReference< string_array >( DeadOilFluid::viewKeyStruct::viscosityTableNamesString() );
  viscosityTableNames.resize( 2 );
  viscosityTableNames[0] = "PVDG_visc"; viscosityTableNames[1] = "PVDO_visc";

  real64 & waterRefPressure = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterRefPressureString() );
  waterRefPressure = 30600000.1;
  real64 & waterFormationVolumeFactor = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterFormationVolumeFactorString() );
  waterFormationVolumeFactor = 1.03;
  real64 & waterCompressibility = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterCompressibilityString() );
  waterCompressibility = 0.00000000041;
  real64 & waterViscosity = fluid.getReference< real64 >( DeadOilFluid::viewKeyStruct::waterViscosityString() );
  waterViscosity = 0.0003;

  fluid.postProcessInputRecursive();
  return fluid;
}

void writeTableToFile( string const & filename, char const * str )
{
  std::ofstream os( filename );
  ASSERT_TRUE( os.is_open() );
  os << str;
  os.close();
}

void removeFile( string const & filename )
{
  int const ret = std::remove( filename.c_str() );
  ASSERT_TRUE( ret == 0 );
}

class LiveOilFluidTest : public CompositionalFluidTestBase
{
public:
  LiveOilFluidTest()
  {
    writeTableToFile( "pvto.txt", pvtoTableContent );
    writeTableToFile( "pvdg.txt", pvdgTableContent );
    writeTableToFile( "pvtw.txt", pvtwTableContent );

    parent.resize( 1 );
    fluid = &makeLiveOilFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~LiveOilFluidTest()
  {
    removeFile( "pvto.txt" );
    removeFile( "pvdg.txt" );
    removeFile( "pvtw.txt" );
  }
};

TEST_F( LiveOilFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.79999; comp[1] = 0.2; comp[2] = 0.00001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-12;

  for( localIndex i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

TEST_F( LiveOilFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.79999; comp[1] = 0.2; comp[2] = 0.00001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-12;

  for( localIndex i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

class DeadOilFluidTest : public CompositionalFluidTestBase
{
public:

  DeadOilFluidTest()
  {
    writeTableToFile( "pvdo.txt", pvdoTableContent );
    writeTableToFile( "pvdg.txt", pvdgTableContent );
    writeTableToFile( "pvdw.txt", pvdwTableContent );

    parent.resize( 1 );
    fluid = &makeDeadOilFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~DeadOilFluidTest()
  {
    removeFile( "pvdo.txt" );
    removeFile( "pvdg.txt" );
    removeFile( "pvdw.txt" );
  }
};

TEST_F( DeadOilFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( localIndex i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

TEST_F( DeadOilFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;
  real64 const absTol = 1e-14;

  for( localIndex i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol, absTol );
  }
}

class DeadOilFluidFromTableTest : public CompositionalFluidTestBase
{
public:

  DeadOilFluidFromTableTest()
  {
    parent.resize( 1 );
    fluid = &makeDeadOilFluidFromTable( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }
};

TEST_F( DeadOilFluidFromTableTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 5.4e6, 1.24e7, 3.21e7 };
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( localIndex i = 0; i < 3; ++i )
  {
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, false, relTol );
  }
}

MultiFluidBase & makeCO2BrinePhillipsFluid( string const & name, Group * parent )
{
  CO2BrinePhillipsFluid & fluid = parent->registerGroup< CO2BrinePhillipsFluid >( name );

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 2 );
  compNames[0] = "co2"; compNames[1] = "water";

  auto & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 2 );
  molarWgt[0] = 44e-3; molarWgt[1] = 18e-3;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "gas"; phaseNames[1] = "liquid";

  auto & phasePVTParaFileNames = fluid.getReference< path_array >( CO2BrinePhillipsFluid::viewKeyStruct::phasePVTParaFilesString() );
  phasePVTParaFileNames.resize( 2 );
  phasePVTParaFileNames[0] = "pvtgas.txt"; phasePVTParaFileNames[1] = "pvtliquid.txt";

  auto & flashModelParaFileName = fluid.getReference< Path >( CO2BrinePhillipsFluid::viewKeyStruct::flashModelParaFileString() );
  flashModelParaFileName = "co2flash.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

class CO2BrinePhillipsFluidTest : public CompositionalFluidTestBase
{
protected:

  CO2BrinePhillipsFluidTest()
  {
    writeTableToFile( "pvtliquid.txt", pvtLiquidPhillipsTableContent );
    writeTableToFile( "pvtgas.txt", pvtGasTableContent );
    writeTableToFile( "co2flash.txt", co2FlashTableContent );

    parent.resize( 1 );
    fluid = &makeCO2BrinePhillipsFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~CO2BrinePhillipsFluidTest()
  {
    removeFile( "pvtliquid.txt" );
    removeFile( "pvtgas.txt" );
    removeFile( "co2flash.txt" );
  }

};


TEST_F( CO2BrinePhillipsFluidTest, checkAgainstPreviousImplementationMolar )
{
  fluid->setMassFlag( false );

  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const relTol = 1e-10;

  fluid->allocateConstitutiveData( fluid->getParent(), 1 );

  CO2BrinePhillipsFluid::KernelWrapper wrapper =
    dynamicCast< CO2BrinePhillipsFluid * >( fluid )->createKernelWrapper();

  real64 const savedTotalDens[] =
  { 5881.8128183956969224, 5869.522096458530541, 5854.9469601674582009, 9180.9455320478591602, 9157.2045503913905122, 9129.1751063784995495, 15755.475565136142905, 15696.691553847707837,
    15627.990771463533747 };
  real64 const savedGasPhaseFrac[] =
  { 0.29413690046142371148, 0.29415754810481165027, 0.29418169867697463449, 0.29194010802017489326, 0.29196434961986583723, 0.29199266189550621142, 0.2890641335638892695, 0.28908718137828937067,
    0.28911404840933618843 };
  real64 const savedWaterDens[] =
  { 53286.457784368176362, 53264.389103437584708, 53237.751306267287873, 53229.257940878436784, 53207.597127679167897, 53181.436584967217641, 53197.49848403003125, 53176.033397316634364,
    53150.105086882285832 };
  real64 const savedGasDens[] =
  { 1876.2436091302606656, 1872.184636376355229, 1867.3711104617746059, 3053.1548401973859654, 3044.5748249030266379, 3034.4507978134674886, 5769.0622621289458039, 5742.8476745352018042,
    5712.2837704249559465 };
  real64 const savedWaterMassDens[] =
  { 970.85108546544745423, 970.4075834766143771, 969.87385780866463847, 974.23383396044232541, 973.78856424100911227, 973.25280170872576946, 979.48333010951580491, 979.04147229150635212,
    978.50977403260912979 };
  real64 const savedGasMassDens[] =
  { 82.554718801731468147, 82.376124000559627802, 82.164328860318079251, 134.33881296868497657, 133.96129229573315911, 133.51583510379256836, 253.83873953367358922, 252.68529767954885301,
    251.34048589869803436 };
  real64 const savedWaterVisc[] =
  { 0.00090094759910161340347, 0.00090096652240945261734, 0.00090098923037885969567, 0.00090094759910161340347, 0.00090096652240945261734, 0.00090098923037885969567, 0.00090094759910161340347,
    0.00090096652240945261734, 0.00090098923037885969567 };
  real64 const savedGasVisc[] =
  { 1.9042384704865343673e-05, 1.9062615947696152414e-05, 1.9086923154230274463e-05, 2.0061713844617985449e-05, 2.0075955757102255573e-05, 2.0093249989250199265e-05, 2.3889596884008691474e-05,
    2.3865756080512667728e-05, 2.3839170076324036522e-05  };
  real64 const savedWaterPhaseGasComp[] =
  { 0.0083062842389820552153, 0.008277274736736653718, 0.0082433415400525456018, 0.011383065290266058955, 0.011349217198060387521, 0.011309682362800700661, 0.015382352969377973903,
    0.015350431636424789056, 0.015313218057419366105  };
  real64 const savedWaterPhaseWaterComp[] =
  { 0.99169371576101794652, 0.9917227252632633272, 0.99175665845994742664, 0.98861693470973388553, 0.98865078280193963156, 0.98869031763719927852, 0.98461764703062204518, 0.98464956836357520054,
    0.98468678194258063563 };

  localIndex counter = 0;
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      testValuesAgainstPreviousImplementation( wrapper,
                                               P[i], T[j], comp,
                                               savedTotalDens[counter], savedGasPhaseFrac[counter],
                                               savedWaterDens[counter], savedGasDens[counter],
                                               savedWaterMassDens[counter], savedGasMassDens[counter],
                                               savedWaterVisc[counter], savedGasVisc[counter],
                                               savedWaterPhaseGasComp[counter], savedWaterPhaseWaterComp[counter],
                                               relTol );
      counter++;
    }
  }
}

TEST_F( CO2BrinePhillipsFluidTest, checkAgainstPreviousImplementationMass )
{
  fluid->setMassFlag( true );

  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const relTol = 1e-10;

  fluid->allocateConstitutiveData( fluid->getParent(), 1 );

  CO2BrinePhillipsFluid::KernelWrapper wrapper =
    dynamicCast< CO2BrinePhillipsFluid * >( fluid )->createKernelWrapper();

  real64 const savedTotalDens[] =
  { 238.33977561940088208, 237.86350488026934613, 237.29874890241927687, 354.01144731214282046, 353.18618684355078585, 352.21120673560858449, 550.02182875764299297, 548.3889751707506548,
    546.47580480217254717 };
  real64 const savedGasPhaseFrac[] =
  { 0.28562868803317220667, 0.28567941665326646028, 0.285738749802139258, 0.28022484140718162404, 0.2802844989853667812, 0.28035417162546172332, 0.2731355646393489045, 0.27319238868618361815,
    0.2732586251114847431 };
  real64 const savedWaterDens[] =
  { 970.85108546544745423, 970.4075834766143771, 969.87385780866463847, 974.23383396044232541, 973.78856424100911227, 973.25280170872576946, 979.48333010951580491, 979.04147229150635212,
    978.50977403260912979 };
  real64 const savedGasDens[] =
  { 82.554718801731468147, 82.376124000559627802, 82.164328860318079251, 134.33881296868497657, 133.96129229573315911, 133.51583510379256836, 253.83873953367358922, 252.68529767954885301,
    251.34048589869803436 };
  real64 const savedWaterMassDens[] =
  { 970.85108546544745423, 970.4075834766143771, 969.87385780866463847, 974.23383396044232541, 973.78856424100911227, 973.25280170872576946, 979.48333010951580491, 979.04147229150635212,
    978.50977403260912979 };
  real64 const savedGasMassDens[] =
  { 82.554718801731468147, 82.376124000559627802, 82.164328860318079251, 134.33881296868497657, 133.96129229573315911, 133.51583510379256836, 253.83873953367358922, 252.68529767954885301,
    251.34048589869803436 };
  real64 const savedWaterVisc[] =
  { 0.00090094759910161340347, 0.00090096652240945261734, 0.00090098923037885969567, 0.00090094759910161340347, 0.00090096652240945261734, 0.00090098923037885969567, 0.00090094759910161340347,
    0.00090096652240945261734, 0.00090098923037885969567 };
  real64 const savedGasVisc[] =
  { 1.9042384704865343673e-05, 1.9062615947696152414e-05, 1.9086923154230274463e-05, 2.0061713844617985449e-05, 2.0075955757102255573e-05, 2.0093249989250199265e-05, 2.3889596884008691474e-05,
    2.3865756080512667728e-05, 2.3839170076324036522e-05  };
  real64 const savedWaterPhaseGasComp[] =
  { 0.02005966592318779787, 0.019990461277537684842, 0.019909503061226688919, 0.027365230280837819082, 0.027285226317914228894, 0.027191770514265831832, 0.036759501299346700187,
    0.036684965747010883641, 0.036598063202886929601 };
  real64 const savedWaterPhaseWaterComp[] =
  { 0.9797478006656266114, 0.97981828292617156873, 0.97990072673671935188, 0.97227194517798976037, 0.97235397979974458327, 0.97244979317916002692, 0.9625743441996873484, 0.9626514061444874093,
    0.962741233539301966 };

  localIndex counter = 0;
  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      testValuesAgainstPreviousImplementation( wrapper,
                                               P[i], T[j], comp,
                                               savedTotalDens[counter], savedGasPhaseFrac[counter],
                                               savedWaterDens[counter], savedGasDens[counter],
                                               savedWaterMassDens[counter], savedGasMassDens[counter],
                                               savedWaterVisc[counter], savedGasVisc[counter],
                                               savedWaterPhaseGasComp[counter], savedWaterPhaseWaterComp[counter],
                                               relTol );
      counter++;
    }
  }
}

TEST_F( CO2BrinePhillipsFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

TEST_F( CO2BrinePhillipsFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-8;

  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

MultiFluidBase & makeCO2BrineEzrokhiFluid( string const & name, Group * parent )
{
  CO2BrineEzrokhiFluid & fluid = parent->registerGroup< CO2BrineEzrokhiFluid >( name );

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 2 );
  compNames[0] = "co2"; compNames[1] = "water";

  auto & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 2 );
  molarWgt[0] = 44e-3; molarWgt[1] = 18e-3;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 2 );
  phaseNames[0] = "gas"; phaseNames[1] = "liquid";

  auto & phasePVTParaFileNames = fluid.getReference< path_array >( CO2BrineEzrokhiFluid::viewKeyStruct::phasePVTParaFilesString() );
  phasePVTParaFileNames.resize( 2 );
  phasePVTParaFileNames[0] = "pvtgas.txt"; phasePVTParaFileNames[1] = "pvtliquid.txt";

  auto & flashModelParaFileName = fluid.getReference< Path >( CO2BrineEzrokhiFluid::viewKeyStruct::flashModelParaFileString() );
  flashModelParaFileName = "co2flash.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

class CO2BrineEzrokhiFluidTest : public CompositionalFluidTestBase
{
protected:

  CO2BrineEzrokhiFluidTest()
  {
    writeTableToFile( "pvtliquid.txt", pvtLiquidEzrokhiTableContent );
    writeTableToFile( "pvtgas.txt", pvtGasTableContent );
    writeTableToFile( "co2flash.txt", co2FlashTableContent );

    parent.resize( 1 );
    fluid = &makeCO2BrineEzrokhiFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~CO2BrineEzrokhiFluidTest()
  {
    removeFile( "pvtliquid.txt" );
    removeFile( "pvtgas.txt" );
    removeFile( "co2flash.txt" );
  }

};

TEST_F( CO2BrineEzrokhiFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

TEST_F( CO2BrineEzrokhiFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const T[3] = { 367.65, 368.15, 368.75 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-8;

  for( localIndex i = 0; i < 3; ++i )
  {
    for( localIndex j = 0; j < 3; ++j )
    {
      testNumericalDerivatives( *fluid, parent, P[i], T[j], comp, eps, false, relTol );
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
