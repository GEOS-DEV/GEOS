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

// Source includes
#include "common/initializeEnvironment.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2Solubility.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "functions/FunctionManager.hpp"

// TPL includes
#include <conduit.hpp>
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
  static std::unique_ptr< CO2Solubility > makeFlashModel();
};

std::unique_ptr< CO2Solubility >
CO2SolubilitySpycherPruessTestFixture::makeFlashModel()
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

  // Fluid model parameters
  BrineFluidParameters parameters;
  parameters.m_solubilityModel = BrineFluidParameters::SolubilityModel::SpycherPruess;
  parameters.m_pressureCoordinates.emplace_back( 1.0e5 );
  parameters.m_pressureCoordinates.emplace_back( 1.0e7 );
  parameters.m_pressureInterval = 9.9e5;
  parameters.m_temperatureCoordinates.emplace_back( 283.15 );
  parameters.m_temperatureCoordinates.emplace_back( 383.15 );
  parameters.m_temperatureInterval = 10.0;
  parameters.m_salinity = 0.15;
  parameters.m_tolerance = 1.0e-8;

  TableFunction::OutputOptions const flashOutputOpts = {
    false,  // writeCSV
    false,   // writeInLog
  };

  return std::make_unique< CO2Solubility >( "FlashModel",
                                            parameters,
                                            phaseNames,
                                            componentNames,
                                            componentMolarWeight,
                                            flashOutputOpts );
}

TEST_P( CO2SolubilitySpycherPruessTestFixture, testExpectedValues )
{
  auto flashModel = makeFlashModel();

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

  auto flashModel = makeFlashModel();

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
stdVector< TestParam > generateTestData()
{
  return {
    {1.00000000e+05, 1.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+05, 5.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+05, 8.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+05, 1.10000000e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 1.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 5.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 8.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 1.10000000e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 1.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 5.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 8.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 1.10000000e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 1.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 5.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 8.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 1.10000000e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 1.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 5.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 8.00000000e+01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 1.10000000e+02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+05, 1.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+05, 5.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+05, 8.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+05, 1.10000000e+02, 1.00000000e-05, 1.00000000e-05, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 1.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+06, 5.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+06, 8.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+06, 1.10000000e+02, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {5.00000000e+06, 1.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {5.00000000e+06, 5.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {5.00000000e+06, 8.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {5.00000000e+06, 1.10000000e+02, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+07, 1.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+07, 5.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+07, 8.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+07, 1.10000000e+02, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+08, 1.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+08, 5.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+08, 8.00000000e+01, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+08, 1.10000000e+02, 1.00000000e-05, 0.00000000e+00, 1.00000000e-05, 0.00000000e+00},
    {1.00000000e+05, 1.00000000e+01, 3.00000000e-01, 0.30316629266172906, 0.00085090329021506083, 1.23998968e-02},
    {1.00000000e+05, 5.00000000e+01, 3.00000000e-01, 0.34221891417081973, 0.00030737250473219488, 1.23958952e-01},
    {1.00000000e+05, 8.00000000e+01, 3.00000000e-01, 0.56694554161147426, 0.00012308123367675458, 4.70942662e-01},
    {1.00000000e+05, 1.10000000e+02, 3.00000000e-01, 3.00000000e-01, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 1.00000000e+01, 3.00000000e-01, 0.29499882319684639, 0.0080450957564913543, 2.27331474e-03},
    {1.00000000e+06, 5.00000000e+01, 3.00000000e-01, 0.30488166349930129, 0.0032989000731056487, 2.35330303e-02},
    {1.00000000e+06, 8.00000000e+01, 3.00000000e-01, 0.33577292210903398, 0.002127595591701263, 1.10747849e-01},
    {1.00000000e+06, 1.10000000e+02, 3.00000000e-01, 0.34096859661779633, 0.0016079594339768342, 1.23261476e-01},
    {5.00000000e+06, 1.00000000e+01, 3.00000000e-01, 0.2799234586625956, 0.02866391855130497, 2.01367187e-03},
    {5.00000000e+06, 5.00000000e+01, 3.00000000e-01, 0.29146281934391888, 0.013570676580556669, 3.69909366e-03},
    {5.00000000e+06, 8.00000000e+01, 3.00000000e-01, 0.29719329085522239, 0.0094347679322459715, 1.28674138e-02},
    {5.00000000e+06, 1.10000000e+02, 3.00000000e-01, 0.30538788828878705, 0.0079712840925638476, 3.57736478e-02},
    {1.00000000e+07, 1.00000000e+01, 3.00000000e-01, 0.27904865261135298, 0.029922554402427436, 2.22670319e-03},
    {1.00000000e+07, 5.00000000e+01, 3.00000000e-01, 0.28710559375554839, 0.019724098763057561, 4.06398711e-03},
    {1.00000000e+07, 8.00000000e+01, 3.00000000e-01, 0.2919777311217493, 0.015313414145722523, 9.65816584e-03},
    {1.00000000e+07, 1.10000000e+02, 3.00000000e-01, 0.2974591103618851, 0.013746139793011094, 2.39237441e-02},
    {1.00000000e+08, 1.00000000e+01, 3.00000000e-01, 0.27904865261135298, 0.029922554402427436, 2.22670319e-03},
    {1.00000000e+08, 5.00000000e+01, 3.00000000e-01, 0.28710559375554839, 0.019724098763057561, 4.06398711e-03},
    {1.00000000e+08, 8.00000000e+01, 3.00000000e-01, 0.2919777311217493, 0.015313414145722523, 9.65816584e-03},
    {1.00000000e+08, 1.10000000e+02, 3.00000000e-01, 0.2974591103618851, 0.013746139793011094, 2.39237441e-02},
    {1.00000000e+05, 1.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+05, 5.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+05, 8.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+05, 1.10000000e+02, 9.99990000e-01, 9.99990000e-01, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 1.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+06, 5.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+06, 8.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+06, 1.10000000e+02, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {5.00000000e+06, 1.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {5.00000000e+06, 5.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {5.00000000e+06, 8.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {5.00000000e+06, 1.10000000e+02, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+07, 1.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+07, 5.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+07, 8.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+07, 1.10000000e+02, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+08, 1.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+08, 5.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+08, 8.00000000e+01, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+08, 1.10000000e+02, 9.99990000e-01, 1.00000000e+00, 0.00000000e+00, 1.00000000e-05},
    {1.00000000e+05, 1.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+05, 5.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+05, 8.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+05, 1.10000000e+02, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 1.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 5.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 8.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+06, 1.10000000e+02, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 1.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 5.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 8.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {5.00000000e+06, 1.10000000e+02, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 1.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 5.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 8.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+07, 1.10000000e+02, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 1.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 5.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 8.00000000e+01, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00},
    {1.00000000e+08, 1.10000000e+02, 1.00000000e+00, 1.00000000e+00, 0.00000000e+00, 0.00000000e+00}
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

  geos::setupEnvironment( argc, argv );

  conduit::Node node;
  dataRepository::Group parent( "parent", node );
  auto functionManager = std::make_unique< FunctionManager >( FunctionManager::catalogName(), &parent );

  int const result = RUN_ALL_TESTS();

  geos::cleanupEnvironment();

  return result;
}
