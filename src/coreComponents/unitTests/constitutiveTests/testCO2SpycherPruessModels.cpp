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
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Solubility.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "functions/FunctionManager.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::testing;
using namespace geos::stringutilities;
using namespace geos::constitutive;
using namespace geos::constitutive::multifluid;
using namespace geos::constitutive::PVTProps;

// Test parameter
using TestParam = std::tuple<
  real64 const,     // Pressure [Pa]
  real64 const,     // Temperature [C]
  real64 const,     // Total co2 mole fraction (z_co2 = 1-z_wat)
  real64 const,     // Expected vapour fraction (V)
  real64 const,     // Expected dissolved CO2 fraction (x_co2)
  real64 const      // Expected vapourised water fraction (y_wat)
  >;

class CO2SolubilitySpycherPruessTestFixture : public ::testing::TestWithParam< TestParam >
{
protected:
  static integer constexpr numPhase = 2;
  static integer constexpr numComp  = 2;
  static integer constexpr numDof   = numComp + 2;
  static real64 constexpr relTol = 1.0e-5;
  static real64 constexpr absTol = 1.0e-7;
  static real64 constexpr pertubation = 1.0e-6;

public:
  CO2SolubilitySpycherPruessTestFixture() = default;
  ~CO2SolubilitySpycherPruessTestFixture() override = default;

protected:
  static std::unique_ptr< CO2Solubility > makeFlashModel( string const & fileContent );
};

std::unique_ptr< CO2Solubility >
CO2SolubilitySpycherPruessTestFixture::makeFlashModel( string const & fileContent )
{
  // Define phase names
  string_array phaseNames;
  phaseNames.resize( 2 );
  phaseNames[0] = "gas";
  phaseNames[1] = "liquid";

  // Define component names and molar weight
  string_array componentNames;
  componentNames.resize( 2 );
  componentNames[0] = "co2";
  componentNames[1] = "water";

  array1d< real64 > componentMolarWeight;
  componentMolarWeight.resize( 2 );
  componentMolarWeight[0] = 44.0e-3;
  componentMolarWeight[1] = 18.0e-3;

  // Read file parameters
  array1d< string > const strs = stringutilities::tokenizeBySpaces< array1d >( fileContent );

  return std::make_unique< CO2Solubility >( strs[1],
                                            strs,
                                            phaseNames,
                                            componentNames,
                                            componentMolarWeight,
                                            false );
}

TEST_P( CO2SolubilitySpycherPruessTestFixture, testExpectedValues )
{
  constexpr char const * fileContent = "FlashModel CO2Solubility 1e5 1e7 9.9e5 283.15 383.15 10.0 0.15 1.0e-8 SpycherPruess";
  auto flashModel = makeFlashModel( fileContent );

  auto [pressure, temperature, z_co2, expected_V, expected_x_co2, expected_y_wat] = GetParam();

  StackArray< real64, 1, numComp > composition( 2 );
  composition[0] = z_co2;
  composition[1] = 1.0 - z_co2;

  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseFrac( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseFrac( 1, 1, numPhase, numDof );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac[0][0], dPhaseFrac[0][0] };

  StackArray< real64, 4, numComp *numPhase, LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhase, numComp );
  StackArray< real64, 5, numDof *numComp *numPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFrac( 1, 1, numPhase, numComp, numDof );
  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 >
  phaseCompFracAndDeriv { phaseCompFrac[0][0], dPhaseCompFrac[0][0] };

  auto flashModelWrapper = flashModel->createKernelWrapper();

  flashModelWrapper.compute( pressure,
                             temperature,
                             composition.toSliceConst(),
                             phaseFracAndDeriv,
                             phaseCompFracAndDeriv );

  real64 const V = phaseFracAndDeriv.value[0];
  real64 const x_co2 = phaseCompFracAndDeriv.value[1][0];
  real64 const y_wat = phaseCompFracAndDeriv.value[0][1];

  checkRelativeError( V, expected_V, relTol, absTol );
  checkRelativeError( x_co2, expected_x_co2, relTol, absTol );
  checkRelativeError( y_wat, expected_y_wat, relTol, absTol );
}

TEST_P( CO2SolubilitySpycherPruessTestFixture, testNumericalDerivatives )
{
  using Deriv = multifluid::DerivativeOffset;

  constexpr char const * fileContent = "FlashModel CO2Solubility 1e5 1e7 9.9e5 283.15 383.15 10.0 0.15 1.0e-8 SpycherPruess";
  auto flashModel = makeFlashModel( fileContent );

  auto [pressure, temperature, z_co2, expected_V, expected_x_co2, expected_y_wat] = GetParam();
  GEOS_UNUSED_VAR( expected_V, expected_x_co2, expected_y_wat );

  StackArray< real64, 1, numComp > composition( numComp );
  composition[0] = z_co2;
  composition[1] = 1.0 - z_co2;
  StackArray< real64, 1, numComp > perturbedComposition( numComp );

  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseFrac( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseFrac( 1, 1, numPhase, numDof );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac[0][0], dPhaseFrac[0][0] };

  StackArray< real64, 4, numComp *numPhase, LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhase, numComp );
  StackArray< real64, 5, numDof *numComp *numPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFrac( 1, 1, numPhase, numComp, numDof );
  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 >
  phaseCompFracAndDeriv { phaseCompFrac[0][0], dPhaseCompFrac[0][0] };

  // 1) First compute the unperturbed parameters

  auto flashModelWrapper = flashModel->createKernelWrapper();

  flashModelWrapper.compute( pressure,
                             temperature,
                             composition.toSliceConst(),
                             phaseFracAndDeriv,
                             phaseCompFracAndDeriv );

  StackArray< real64, 3, numPhase, LAYOUT_PHASE > perturbedPhaseFrac( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPerturbedPhaseFrac( 1, 1, numPhase, numDof );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  perturbedPhaseFracAndDeriv { perturbedPhaseFrac[0][0], dPerturbedPhaseFrac[0][0] };

  StackArray< real64, 4, numComp *numPhase, LAYOUT_PHASE_COMP > perturbedPhaseCompFrac( 1, 1, numPhase, numComp );
  StackArray< real64, 5, numDof *numComp *numPhase, LAYOUT_PHASE_COMP_DC > dPerturbedPhaseCompFrac( 1, 1, numPhase, numComp, numDof );
  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 >
  perturbedPhaseCompFracAndDeriv { perturbedPhaseCompFrac[0][0], dPerturbedPhaseCompFrac[0][0] };

  real64 numericalDerivative = 0.0;
  real64 analyticalDerivative = 0.0;

  // 2) Check derivative with respect to pressure
  real64 const dP = pertubation * pressure;
  flashModelWrapper.compute( pressure + dP,
                             temperature,
                             composition.toSliceConst(),
                             perturbedPhaseFracAndDeriv,
                             perturbedPhaseCompFracAndDeriv );

  for( integer i = 0; i < numPhase; ++i )
  {
    numericalDerivative = (perturbedPhaseFracAndDeriv.value[i]-phaseFracAndDeriv.value[i])/dP;
    analyticalDerivative = perturbedPhaseFracAndDeriv.derivs[i][Deriv::dP];
    checkRelativeError( numericalDerivative, analyticalDerivative, relTol, absTol );
    for( integer j = 0; j < numComp; ++j )
    {
      numericalDerivative = (perturbedPhaseCompFracAndDeriv.value[i][j]-phaseCompFracAndDeriv.value[i][j])/dP;
      analyticalDerivative = perturbedPhaseCompFracAndDeriv.derivs[i][j][Deriv::dP];
      checkRelativeError( numericalDerivative, analyticalDerivative, relTol, absTol );
    }
  }

  // 3) Check derivative with respect to temperature
  real64 const dT = pertubation * temperature;
  flashModelWrapper.compute( pressure,
                             temperature + dT,
                             composition.toSliceConst(),
                             perturbedPhaseFracAndDeriv,
                             perturbedPhaseCompFracAndDeriv );

  for( integer i = 0; i < numPhase; ++i )
  {
    numericalDerivative = (perturbedPhaseFracAndDeriv.value[i]-phaseFracAndDeriv.value[i])/dT;
    analyticalDerivative = perturbedPhaseFracAndDeriv.derivs[i][Deriv::dT];
    checkRelativeError( numericalDerivative, analyticalDerivative, relTol, absTol );

    for( integer j = 0; j < numComp; ++j )
    {
      numericalDerivative = (perturbedPhaseCompFracAndDeriv.value[i][j]-phaseCompFracAndDeriv.value[i][j])/dT;
      analyticalDerivative = perturbedPhaseCompFracAndDeriv.derivs[i][j][Deriv::dT];
      checkRelativeError( numericalDerivative, analyticalDerivative, relTol, absTol );
    }
  }

  // 4) Check derivative with respect to composition
  for( integer ic=0; ic < numComp; ++ic )
  {
    real64 const dC = pertubation;
    for( integer k=0; k < numComp; ++k )
      perturbedComposition[k] = composition[k];
    perturbedComposition[ic] += dC;

    flashModelWrapper.compute( pressure,
                               temperature,
                               perturbedComposition.toSliceConst(),
                               perturbedPhaseFracAndDeriv,
                               perturbedPhaseCompFracAndDeriv );

    for( integer i = 0; i < numPhase; ++i )
    {
      numericalDerivative = (perturbedPhaseFracAndDeriv.value[i]-phaseFracAndDeriv.value[i])/dC;
      analyticalDerivative = perturbedPhaseFracAndDeriv.derivs[i][Deriv::dC+ic];
      checkRelativeError( numericalDerivative, analyticalDerivative, relTol, absTol );

      for( integer j = 0; j < numComp; ++j )
      {
        numericalDerivative = (perturbedPhaseCompFracAndDeriv.value[i][j]-phaseCompFracAndDeriv.value[i][j])/dC;
        analyticalDerivative = perturbedPhaseCompFracAndDeriv.derivs[i][j][Deriv::dC+ic];
        checkRelativeError( numericalDerivative, analyticalDerivative, relTol, absTol );
      }
    }
  }
}

// Test data
std::vector< TestParam > generateTestData()
{
  return {
    {1.0000000000e+05, 1.0000000000e+01, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+06, 2.0000000000e+01, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {5.0000000000e+06, 1.1000000000e+02, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+07, 1.1000000000e+02, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+05, 1.0000000000e+01, 1.0000000000e-05, 0.0000000000e+00, 1.0000000000e-05, 0.0000000000e+00},
    {1.0000000000e+06, 1.0000000000e+01, 1.0000000000e-05, 0.0000000000e+00, 1.0000000000e-05, 0.0000000000e+00},
    {1.0000000000e+07, 5.0000000000e+01, 1.0000000000e-05, 0.0000000000e+00, 1.0000000000e-05, 0.0000000000e+00},
    {1.0000000000e+05, 1.0000000000e+02, 1.0000000000e-05, 1.1768795168e-03, 1.4667486313e-06, 9.9274778987e-01},
    {1.0000000000e+07, 1.1000000000e+02, 1.0000000000e-05, 0.0000000000e+00, 1.0000000000e-05, 0.0000000000e+00},
    {1.0000000000e+05, 1.0000000000e+01, 2.0000000000e-01, 2.0181739006e-01, 8.5835590596e-04, 1.2399896845e-02},
    {1.0000000000e+06, 1.0000000000e+01, 2.0000000000e-01, 1.9380702561e-01, 8.2261526034e-03, 2.2645827445e-03},
    {2.0000000000e+06, 1.0000000000e+01, 2.0000000000e-01, 1.8748831090e-01, 1.5573461467e-02, 7.5700923149e-04},
    {3.0000000000e+06, 1.0000000000e+01, 2.0000000000e-01, 1.8210136343e-01, 2.2006648400e-02, 5.5228118318e-04},
    {5.0000000000e+06, 1.0000000000e+01, 2.0000000000e-01, 1.7480487929e-01, 3.0942122380e-02, 1.9345438433e-03},
    {1.0000000000e+07, 1.0000000000e+01, 2.0000000000e-01, 1.7167871541e-01, 3.4616909097e-02, 2.0540579107e-03},
    {3.0000000000e+06, 1.5000000000e+01, 2.0000000000e-01, 1.8414138671e-01, 1.9615654890e-02, 7.8737167031e-04},
    {1.0000000000e+05, 2.0000000000e+01, 2.0000000000e-01, 2.0428469777e-01, 6.5298948908e-04, 2.3517627865e-02},
    {3.0000000000e+06, 2.0000000000e+01, 2.0000000000e-01, 1.8618214427e-01, 1.7212941749e-02, 1.0223515873e-03},
    {1.0000000000e+06, 5.0000000000e+01, 2.0000000000e-01, 2.0205285848e-01, 3.3678845233e-03, 2.3460456550e-02},
    {2.0000000000e+06, 5.0000000000e+01, 2.0000000000e-01, 1.9616766138e-01, 6.5333931378e-03, 7.2357189752e-03},
    {1.0000000000e+07, 5.0000000000e+01, 2.0000000000e-01, 1.8234424244e-01, 2.2438771918e-02, 3.7919129223e-03},
    {1.0000000000e+05, 1.0000000000e+02, 2.0000000000e-01, 1.0000000000e+00, 0.0000000000e+00, 8.0000000000e-01},
    {1.0000000000e+06, 1.0000000000e+02, 2.0000000000e-01, 1.0000000000e+00, 0.0000000000e+00, 8.0000000000e-01},
    {2.0000000000e+06, 1.0000000000e+02, 2.0000000000e-01, 2.0878320780e-01, 3.6560536960e-03, 5.5923745016e-02},
    {3.0000000000e+06, 1.0000000000e+02, 2.0000000000e-01, 2.0344428398e-01, 5.4491512731e-03, 3.8265201781e-02},
    {1.0000000000e+07, 1.0000000000e+02, 2.0000000000e-01, 1.9054484917e-01, 1.5654545984e-02, 1.6880551049e-02},
    {1.0000000000e+05, 1.1000000000e+02, 2.0000000000e-01, 2.0000000000e-01, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+06, 1.1000000000e+02, 2.0000000000e-01, 1.9863195029e-01, 1.7071428219e-03, 0.0000000000e+00},
    {2.0000000000e+06, 1.1000000000e+02, 2.0000000000e-01, 2.1398644726e-01, 3.4584564637e-03, 7.8064948154e-02},
    {1.0000000000e+05, 1.0000000000e+01, 5.0000000000e-01, 5.0584830883e-01, 8.5835590596e-04, 1.2399896845e-02},
    {1.0000000000e+06, 1.0000000000e+01, 5.0000000000e-01, 4.9698761291e-01, 8.2261526034e-03, 2.2645827445e-03},
    {2.0000000000e+06, 1.0000000000e+01, 5.0000000000e-01, 4.9246878561e-01, 1.5573461467e-02, 7.5700923149e-04},
    {3.0000000000e+06, 1.0000000000e+01, 5.0000000000e-01, 4.8902523750e-01, 2.2006648400e-02, 5.5228118318e-04},
    {1.0000000000e+06, 1.5000000000e+01, 5.0000000000e-01, 4.9797507513e-01, 7.2906770529e-03, 3.2836517466e-03},
    {1.0000000000e+05, 2.0000000000e+01, 5.0000000000e-01, 5.1171548981e-01, 6.5298948908e-04, 2.3517627865e-02},
    {5.0000000000e+06, 2.0000000000e+01, 5.0000000000e-01, 4.8732521146e-01, 2.5470263012e-02, 7.8627814230e-04},
    {1.0000000000e+07, 2.0000000000e+01, 5.0000000000e-01, 4.8519098453e-01, 3.1323342486e-02, 2.7134132315e-03},
    {1.0000000000e+07, 5.0000000000e+01, 5.0000000000e-01, 4.9042542278e-01, 2.2438771918e-02, 3.7919129223e-03},
    {1.0000000000e+05, 1.0000000000e+02, 5.0000000000e-01, 1.0000000000e+00, 0.0000000000e+00, 5.0000000000e-01},
    {1.0000000000e+06, 1.0000000000e+02, 5.0000000000e-01, 1.0000000000e+00, 0.0000000000e+00, 5.0000000000e-01},
    {2.0000000000e+06, 1.0000000000e+02, 5.0000000000e-01, 5.2778954091e-01, 3.6560536960e-03, 5.5923745016e-02},
    {5.0000000000e+06, 1.0000000000e+02, 5.0000000000e-01, 5.0852379246e-01, 8.7725862968e-03, 2.5240333089e-02},
    {1.0000000000e+07, 1.0000000000e+02, 5.0000000000e-01, 5.0063361733e-01, 1.5654545984e-02, 1.6880551049e-02},
    {1.0000000000e+05, 1.1000000000e+02, 5.0000000000e-01, 5.0000000000e-01, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+06, 1.1000000000e+02, 5.0000000000e-01, 4.9903272484e-01, 1.9308150658e-03, 0.0000000000e+00},
    {1.0000000000e+07, 1.1000000000e+02, 5.0000000000e-01, 5.0375053647e-01, 1.5364362528e-02, 2.2580805990e-02},
    {1.0000000000e+05, 1.0000000000e+01, 9.9999000000e-01, 1.0000000000e+00, 0.0000000000e+00, 1.0000000000e-05},
    {1.0000000000e+06, 1.0000000000e+01, 9.9999000000e-01, 1.0000000000e+00, 0.0000000000e+00, 1.0000000000e-05},
    {2.0000000000e+06, 5.0000000000e+01, 9.9999000000e-01, 1.0000000000e+00, 0.0000000000e+00, 1.0000000000e-05},
    {1.0000000000e+05, 1.1000000000e+02, 9.9999000000e-01, 9.9999000000e-01, 0.0000000000e+00, 0.0000000000e+00},
    {2.0000000000e+06, 1.1000000000e+02, 9.9999000000e-01, 1.0000000000e+00, 0.0000000000e+00, 1.0000000000e-05},
    {3.0000000000e+06, 1.1000000000e+02, 9.9999000000e-01, 1.0000000000e+00, 0.0000000000e+00, 1.0000000000e-05},
    {1.0000000000e+05, 1.0000000000e+01, 1.0000000000e+00, 1.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+07, 1.0000000000e+01, 1.0000000000e+00, 1.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+06, 2.0000000000e+01, 1.0000000000e+00, 1.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+07, 5.0000000000e+01, 1.0000000000e+00, 1.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00},
    {1.0000000000e+07, 1.1000000000e+02, 1.0000000000e+00, 1.0000000000e+00, 0.0000000000e+00, 0.0000000000e+00}
  };
}

INSTANTIATE_TEST_SUITE_P(
  CO2SolubilitySpycherPruessTest,
  CO2SolubilitySpycherPruessTestFixture,
  ::testing::ValuesIn( generateTestData() )
  );

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
