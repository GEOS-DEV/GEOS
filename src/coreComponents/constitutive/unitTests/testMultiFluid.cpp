/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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
#include "physicsSolvers/fluidFlow/unitTests/testCompFlowUtils.hpp"
#include "managers/initialization.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/GeosxState.hpp"


// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;

/// Black-oil tables written into temporary files during testing

static const char * pvtg_str = "#\tPg(Pa)\t\tRv(sm3/sm3)\tBg(m3/sm3)\tVisc(Pa.s)\n"
                               "\n"
                               "\t3000000\t\t0.000132\t0.04234\t    0.00001344\n"
                               "\t\t\t\t0\t\t\t0.04231\t    0.00001389\n"
                               "\t6000000\t\t0.000124\t0.02046\t    0.0000142\n"
                               "\t\t\t\t0\t\t\t0.02043\t    0.0000145\n"
                               "\t9000000\t\t0.000126\t0.01328\t    0.00001526\n"
                               "\t\t\t\t0\t\t\t0.01325\t    0.00001532\n"
                               "   12000000\t\t0.000135\t0.00977\t    0.0000166\n"
                               "\t\t\t\t0\t\t\t0.00973\t    0.00001634\n"
                               "   15000000\t\t0.000149\t0.00773\t    0.00001818\n"
                               "\t\t\t\t0\t\t\t0.00769\t    0.00001752\n"
                               "   18000000\t\t0.000163\t0.006426\t0.00001994\n"
                               "\t\t\t\t0\t\t\t0.006405\t0.00001883\n"
                               "   21000000\t\t0.000191\t0.005541\t0.00002181\n"
                               "\t\t\t\t0\t\t\t0.005553\t0.00002021\n"
                               "   24000000\t\t0.000225\t0.004919\t0.0000237\n"
                               "\t\t\t\t0\t\t\t0.004952\t0.00002163\n"
                               "   27000000\t\t0.000272\t0.004471\t0.00002559\n"
                               "\t\t\t\t0\t\t\t0.004511\t0.00002305\n"
                               "   29500000\t\t0.000354\t0.004194\t0.00002714\n"
                               "\t\t\t\t0\t\t\t0.004225\t0.00002423\n"
                               "   31000000\t\t0.000403\t0.004031\t0.00002806\n"
                               "\t\t\t\t0.000354\t0.004059\t0.00002768\n"
                               "   33000000\t\t0.000354\t0.00391\t    0.00002832\n"
                               "\t\t\t\t0\t\t\t0.003913\t0.00002583\n"
                               "   53000000\t\t0.000479\t0.003868\t0.00002935\n"
                               "\t\t\t\t0.000354\t0.0039\t\t0.00002842\n"
                               "\t\t\t\t0\t\t\t0.003903\t0.00002593";

static const char * pvto_str = "# Rs[sm3/sm3]\tPbub[Pa]\tBo[m3/sm3]\tVisc(Pa.s)\n"
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

static const char * pvtw_str = "#\tPref[bar]\tBw[m3/sm3]\tCp[1/bar]\t    Visc[cP]\n"
                               "\t30600000.1\t1.03\t\t0.00000000041\t0.0003";

/// Dead-oil tables written into temporary files during testing

static const char * pvdg_str = "# Pg(Pa) Bg(m3/sm3) Visc(Pa.s)\n"
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

static const char * pvdo_str = "#P[Pa] Bo[m3/sm3] Visc(Pa.s)\n"
                               "2000000  1.02 0.000975\n"
                               "5000000  1.03 0.00091\n"
                               "10000000 1.04 0.00083\n"
                               "20000000 1.05 0.000695\n"
                               "30000000 1.07 0.000594\n"
                               "40000000 1.08 0.00051\n"
                               "50000000.7 1.09 0.000449";

static const char * pvdw_str = "# Pref[bar] Bw[m3/sm3] Cp[1/bar]     Visc[cP]\n"
                               " 30600000.1 1.03  0.00000000041 0.0003";

void testNumericalDerivatives( MultiFluidBase & fluid,
                               Group & parent,
                               real64 const P,
                               real64 const T,
                               arraySlice1d< real64 > const & composition,
                               real64 const perturbParameter,
                               real64 const relTol,
                               real64 const absTol = std::numeric_limits< real64 >::max() )
{
  localIndex const NC = fluid.numFluidComponents();
  localIndex const NP = fluid.numFluidPhases();

  auto const & components = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  auto const & phases     = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );

  // create a clone of the fluid to run updates on
  std::unique_ptr< ConstitutiveBase > fluidCopyPtr = fluid.deliverClone( "fluidCopy", &parent );
  MultiFluidBase & fluidCopy = dynamicCast< MultiFluidBase & >( *fluidCopyPtr );

  fluid.allocateConstitutiveData( fluid.getParent(), 1 );
  fluidCopy.allocateConstitutiveData( fluid.getParent(), 1 );

  // extract data views from both fluids
  #define GET_FLUID_DATA( FLUID, DIM, KEY ) \
    FLUID.getReference< Array< real64, DIM > >( MultiFluidBase::viewKeyStruct::KEY() )[0][0]

  CompositionalVarContainer< 1 > phaseFrac {
    GET_FLUID_DATA( fluid, 3, phaseFractionString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseFraction_dGlobalCompFractionString )
  };

  CompositionalVarContainer< 1 > phaseDens {
    GET_FLUID_DATA( fluid, 3, phaseDensityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseDensity_dGlobalCompFractionString )
  };

  CompositionalVarContainer< 1 > phaseVisc {
    GET_FLUID_DATA( fluid, 3, phaseViscosityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseViscosity_dGlobalCompFractionString )
  };

  CompositionalVarContainer< 2 > phaseCompFrac {
    GET_FLUID_DATA( fluid, 4, phaseCompFractionString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 5, dPhaseCompFraction_dGlobalCompFractionString )
  };

  CompositionalVarContainer< 0 > totalDens {
    GET_FLUID_DATA( fluid, 2, totalDensityString ),
    GET_FLUID_DATA( fluid, 2, dTotalDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 2, dTotalDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 3, dTotalDensity_dGlobalCompFractionString )
  };

  auto const & phaseFracCopy     = GET_FLUID_DATA( fluidCopy, 3, phaseFractionString );
  auto const & phaseDensCopy     = GET_FLUID_DATA( fluidCopy, 3, phaseDensityString );
  auto const & phaseViscCopy     = GET_FLUID_DATA( fluidCopy, 3, phaseViscosityString );
  auto const & phaseCompFracCopy = GET_FLUID_DATA( fluidCopy, 4, phaseCompFractionString );
  auto const & totalDensCopy     = GET_FLUID_DATA( fluidCopy, 2, totalDensityString );

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

      checkDerivative( phaseFracCopy, phaseFrac.value, phaseFrac.dPres, dP, relTol, absTol, "phaseFrac", "Pres", phases );
      checkDerivative( phaseDensCopy, phaseDens.value, phaseDens.dPres, dP, relTol, absTol, "phaseDens", "Pres", phases );
      checkDerivative( phaseViscCopy, phaseVisc.value, phaseVisc.dPres, dP, relTol, absTol, "phaseVisc", "Pres", phases );
      checkDerivative( totalDensCopy, totalDens.value, totalDens.dPres, dP, relTol, absTol, "totalDens", "Pres" );
      checkDerivative( phaseCompFracCopy.toSliceConst(),
                       phaseCompFrac.value.toSliceConst(),
                       phaseCompFrac.dPres.toSliceConst(),
                       dP,
                       relTol,
                       absTol,
                       "phaseCompFrac",
                       "Pres",
                       phases,
                       components );
      std::cout << phaseDens.value << std::endl;
    }

    // update temperature and check derivatives
    {
      real64 const dT = perturbParameter * (T + perturbParameter);
      fluidWrapper.update( 0, 0, P, T + dT, composition );

      checkDerivative( phaseFracCopy, phaseFrac.value, phaseFrac.dTemp, dT, relTol, absTol, "phaseFrac", "Temp", phases );
      checkDerivative( phaseDensCopy, phaseDens.value, phaseDens.dTemp, dT, relTol, absTol, "phaseDens", "Temp", phases );
      checkDerivative( phaseViscCopy, phaseVisc.value, phaseVisc.dTemp, dT, relTol, absTol, "phaseVisc", "Temp", phases );
      checkDerivative( totalDensCopy, totalDens.value, totalDens.dTemp, dT, relTol, absTol, "totalDens", "Temp" );
      checkDerivative( phaseCompFracCopy.toSliceConst(),
                       phaseCompFrac.value.toSliceConst(),
                       phaseCompFrac.dTemp.toSliceConst(),
                       dT,
                       relTol,
                       absTol,
                       "phaseCompFrac",
                       "Temp",
                       phases,
                       components );
    }

    // update composition and check derivatives
    auto dPhaseFrac_dC     = invertLayout( phaseFrac.dComp, NP, NC );
    auto dPhaseDens_dC     = invertLayout( phaseDens.dComp, NP, NC );
    auto dPhaseVisc_dC     = invertLayout( phaseVisc.dComp, NP, NC );
    auto dTotalDens_dC     = invertLayout( totalDens.dComp, NC );
    auto dPhaseCompFrac_dC = invertLayout( phaseCompFrac.dComp, NP, NC, NC );

    array1d< real64 > compNew( NC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      real64 const dC = perturbParameter * ( composition[jc] + perturbParameter );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compNew[ic] = composition[ic];
      }
      compNew[jc] += dC;

      // renormalize
      real64 sum = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic )
        sum += compNew[ic];
      for( localIndex ic = 0; ic < NC; ++ic )
        compNew[ic] /= sum;

      fluidWrapper.update( 0, 0, P, T, compNew );

      string const var = "compFrac[" + components[jc] + "]";
      checkDerivative( phaseFracCopy, phaseFrac.value, dPhaseFrac_dC[jc], dC, relTol, absTol, "phaseFrac", var, phases );
      checkDerivative( phaseDensCopy, phaseDens.value, dPhaseDens_dC[jc], dC, relTol, absTol, "phaseDens", var, phases );
      checkDerivative( phaseViscCopy, phaseVisc.value, dPhaseVisc_dC[jc], dC, relTol, absTol, "phaseVisc", var, phases );
      checkDerivative( totalDensCopy, totalDens.value, dTotalDens_dC[jc], dC, relTol, absTol, "totalDens", var );
      checkDerivative( phaseCompFracCopy.toSliceConst(), phaseCompFrac.value.toSliceConst(), dPhaseCompFrac_dC[jc].toSliceConst(), dC, relTol, absTol,
                       "phaseCompFrac", var, phases, components );
    }
  } );
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

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, relTol );
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

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, relTol );
}

MultiFluidBase & makeLiveOilFluid( string const & name, Group * parent )
{
  BlackOilFluid & fluid = parent->registerGroup< BlackOilFluid >( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  string_array & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames.resize( 3 );
  compNames[0] = "oil"; compNames[1] = "gas"; compNames[2] = "water";

  array1d< real64 > & molarWgt = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  molarWgt.resize( 3 );
  molarWgt[0] = 114e-3; molarWgt[1] = 16e-3; molarWgt[2] = 18e-3;

  string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil"; phaseNames[1] = "gas"; phaseNames[2] = "water";

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluid::viewKeyStruct::surfaceDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0; surfaceDens[1] = 0.9907; surfaceDens[2] = 1022.0;

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluid::viewKeyStruct::tableFilesString() );
  tableNames.resize( 3 );
  tableNames[0] = "pvto.txt"; tableNames[1] = "pvtg.txt"; tableNames[2] = "pvtw.txt";

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

  array1d< real64 > & surfaceDens = fluid.getReference< array1d< real64 > >( BlackOilFluid::viewKeyStruct::surfaceDensitiesString() );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0; surfaceDens[1] = 1022.0; surfaceDens[2] = 0.9907;

  path_array & tableNames = fluid.getReference< path_array >( BlackOilFluid::viewKeyStruct::tableFilesString() );
  tableNames.resize( 3 );
  tableNames[0] = "pvdo.txt"; tableNames[1] = "pvdw.txt"; tableNames[2] = "pvdg.txt";

  fluid.postProcessInputRecursive();
  return fluid;
}

MultiFluidBase & makeDeadOilFluidFromTable( string const & name, Group * parent )
{
  FunctionManager & functionManager = getGlobalState().getFunctionManager();

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
    writeTableToFile( "pvto.txt", pvto_str );
    writeTableToFile( "pvtg.txt", pvtg_str );
    writeTableToFile( "pvtw.txt", pvtw_str );

    parent.resize( 1 );
    fluid = &makeLiveOilFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~LiveOilFluidTest()
  {
    removeFile( "pvto.txt" );
    removeFile( "pvtg.txt" );
    removeFile( "pvtw.txt" );
  }
};

TEST_F( LiveOilFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-4;

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, relTol );
}

TEST_F( LiveOilFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1; comp[1] = 0.3; comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 1e-2;
  real64 const absTol = 1e-14;

  testNumericalDerivatives( *fluid, parent, P, T, comp, eps, relTol, absTol );
}

class DeadOilFluidTest : public CompositionalFluidTestBase
{
public:

  DeadOilFluidTest()
  {
    writeTableToFile( "pvdo.txt", pvdo_str );
    writeTableToFile( "pvdg.txt", pvdg_str );
    writeTableToFile( "pvdw.txt", pvdw_str );

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
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, relTol );
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
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, relTol, absTol );
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
    testNumericalDerivatives( *fluid, parent, P[i], T, comp, eps, relTol );
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
