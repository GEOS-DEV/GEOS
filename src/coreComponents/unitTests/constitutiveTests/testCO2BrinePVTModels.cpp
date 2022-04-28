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
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineViscosity.hpp"
#include "constitutive/fluid/PVTFunctions/EzrokhiBrineViscosity.hpp"
#include "constitutive/fluid/PVTFunctions/FenghourCO2Viscosity.hpp"
#include "constitutive/fluid/PVTFunctions/PhillipsBrineDensity.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2Density.hpp"
#include "constitutive/fluid/PVTFunctions/CO2Solubility.hpp"
#include "constitutive/fluid/PVTFunctions/BrineEnthalpy.hpp"
#include "constitutive/fluid/PVTFunctions/CO2Enthalpy.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;
using namespace geosx::stringutilities;
using namespace geosx::constitutive::PVTProps;
using namespace geosx::constitutive::multifluid;

/// Input tables written into temporary files during testing

static const char * pvtLiquidPhillipsTableContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                                    "ViscosityFun PhillipsBrineViscosity 0.1";

// Used also for gas phase
static const char * pvtLiquidEnthalpyTableContent = "EnthalpyFun BrineEnthalpy 1e6 1.5e7 5e4 367.15 369.15 1 0.2";

// the last are set relatively high (1e-4) to increase derivative value and check it properly
static const char * pvtLiquidEzrokhiTableContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                                   "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4";

static const char * pvtGasTableContent = "DensityFun SpanWagnerCO2Density 1e5 7.5e7 5e4 285.15 369.15 4.0\n" // we want to test the full
                                                                                                             // (pres, temp) range here
                                         "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1.0";

static const char * co2FlashTableContent = "FlashModel CO2Solubility 1e5 7.5e7 5e4 285.15 369.15 4.0 0.15"; // we want to test the full
                                                                                                            // (pres, temp) range here

template< typename PVT_WRAPPER >
void testValuesAgainstPreviousImplementation( PVT_WRAPPER const & pvtFunctionWrapper,
                                              real64 const pressure,
                                              real64 const temperature,
                                              arraySlice1d< real64 const > const & phaseComposition,
                                              real64 const & oldImplValue,
                                              bool const useMass,
                                              real64 const relTol )
{
  integer constexpr numPhase = 2;
  integer constexpr numComp  = 2;
  integer constexpr numDof   = numComp + 2;

  real64 value = 0.0;
  stackArray1d< real64, numDof > dValue( numDof );
  stackArray2d< real64, numDof *numPhase > dPhaseComposition( numPhase, numDof );

  pvtFunctionWrapper.compute( pressure,
                              temperature,
                              phaseComposition,
                              dPhaseComposition.toSliceConst(),
                              value,
                              dValue.toSlice(),
                              useMass );

  checkRelativeError( value, oldImplValue, relTol );
}

template< typename FLASH_WRAPPER >
void testValuesAgainstPreviousImplementation( FLASH_WRAPPER const & flashModelWrapper,
                                              real64 const & pressure,
                                              real64 const & temperature,
                                              arraySlice1d< real64 const > const & compFraction,
                                              real64 const & savedGasPhaseFrac,
                                              real64 const & savedWaterPhaseGasComp,
                                              real64 const relTol )
{
  integer constexpr numPhase = 2;
  integer constexpr numComp  = 2;
  integer constexpr numDof   = numComp + 2;

  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseFrac( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseFrac( 1, 1, numPhase, numDof );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac[0][0], dPhaseFrac[0][0] };

  StackArray< real64, 4, numComp *numPhase, LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhase, numComp );
  StackArray< real64, 5, numDof *numComp *numPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFrac( 1, 1, numPhase, numComp, numDof );
  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 >
  phaseCompFracAndDeriv { phaseCompFrac[0][0], dPhaseCompFrac[0][0] };

  flashModelWrapper.compute( pressure,
                             temperature,
                             compFraction,
                             phaseFracAndDeriv,
                             phaseCompFracAndDeriv );

  for( integer i = 0; i < numPhase; ++i )
  {
    real64 const savedPhaseFrac = (i == 0) ? savedGasPhaseFrac : 1.0 - savedGasPhaseFrac;
    checkRelativeError( phaseFracAndDeriv.value[i], savedPhaseFrac, relTol );

    for( integer j = 0; j < numComp; ++j )
    {
      real64 savedCompFrac = 0.0;
      if( i == 0 )
      {
        savedCompFrac = ( j == 0 ) ? 1 : 0;
      }
      else
      {
        savedCompFrac = ( j == 0 ) ? savedWaterPhaseGasComp : 1 - savedWaterPhaseGasComp;
      }
      checkRelativeError( phaseCompFracAndDeriv.value[i][j], savedCompFrac, relTol );
    }
  }
}

template< typename PVT_WRAPPER >
void testNumericalDerivatives( PVT_WRAPPER const & pvtFunctionWrapper,
                               real64 const pressure,
                               real64 const temperature,
                               arraySlice1d< real64 const > const & phaseComposition,
                               bool const useMass,
                               real64 const perturbParameter,
                               real64 const relTol )
{
  using Deriv = multifluid::DerivativeOffset;

  integer constexpr numComp  = 2;
  integer constexpr numDof   = numComp + 2;

  // 1) First compute the unperturbed pressure
  real64 value = 0.0;
  stackArray1d< real64, numDof > dValue( numDof );
  stackArray2d< real64, numDof *numComp > dPhaseComposition( numComp, numDof );

  dPhaseComposition[0][Deriv::dC]   = 1.0;
  dPhaseComposition[1][Deriv::dC+1] = 1.0;
  pvtFunctionWrapper.compute( pressure,
                              temperature,
                              phaseComposition,
                              dPhaseComposition.toSliceConst(),
                              value,
                              dValue.toSlice(),
                              useMass );
  real64 perturbedValue = 0.0;
  stackArray1d< real64, numDof > dPerturbedValue( numDof );

  // 2) Check derivative with respect to pressure
  real64 const dP = perturbParameter * (pressure + perturbParameter);
  pvtFunctionWrapper.compute( pressure + dP,
                              temperature,
                              phaseComposition,
                              dPhaseComposition.toSliceConst(),
                              perturbedValue,
                              dPerturbedValue.toSlice(),
                              useMass );
  checkRelativeError( (perturbedValue-value)/dP, dValue[Deriv::dP], relTol );

  // 3) Check derivative with respect to temperature
  real64 const dT = perturbParameter * (temperature + perturbParameter);
  pvtFunctionWrapper.compute( pressure,
                              temperature + dT,
                              phaseComposition,
                              dPhaseComposition.toSliceConst(),
                              perturbedValue,
                              dPerturbedValue.toSlice(),
                              useMass );
  checkRelativeError( (perturbedValue-value)/dT, dValue[Deriv::dT], relTol );

  // 4) Check derivatives with respect to phaseComposition
  for( integer i = 0; i < numComp; ++i )
  {
    real64 const dC = perturbParameter * (phaseComposition[i] + perturbParameter);
    stackArray1d< real64, numComp > perturbedPhaseComposition( numComp );
    for( integer j = 0; j < numComp; ++j )
    {
      perturbedPhaseComposition[j] = phaseComposition[j] + ( (i == j) ? dC : 0.0 );
    }
    pvtFunctionWrapper.compute( pressure,
                                temperature,
                                perturbedPhaseComposition.toSliceConst(),
                                dPhaseComposition.toSliceConst(),
                                perturbedValue,
                                dPerturbedValue.toSlice(),
                                useMass );
    checkRelativeError( (perturbedValue-value)/dC, dValue[Deriv::dC+i], relTol );
  }
}

template< typename FLASH_WRAPPER >
void testNumericalDerivatives( FLASH_WRAPPER const & flashModelWrapper,
                               real64 const pressure,
                               real64 const temperature,
                               arraySlice1d< real64 const > const & compFraction,
                               real64 const perturbParameter,
                               real64 const relTol )
{
  using namespace multifluid;
  using Deriv = multifluid::DerivativeOffset;

  integer constexpr numPhase = 2;
  integer constexpr numComp  = 2;
  integer constexpr numDof   = numComp + 2;

  // 1) First compute the unperturbed pressure

  StackArray< real64, 3, numPhase, LAYOUT_PHASE > phaseFrac( 1, 1, numPhase );
  StackArray< real64, 4, numDof *numPhase, LAYOUT_PHASE_DC > dPhaseFrac( 1, 1, numPhase, numDof );
  MultiFluidVarSlice< real64, 1, USD_PHASE - 2, USD_PHASE_DC - 2 >
  phaseFracAndDeriv { phaseFrac[0][0], dPhaseFrac[0][0] };

  StackArray< real64, 4, numComp *numPhase, LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhase, numComp );
  StackArray< real64, 5, numDof *numComp *numPhase, LAYOUT_PHASE_COMP_DC > dPhaseCompFrac( 1, 1, numPhase, numComp, numDof );
  MultiFluidVarSlice< real64, 2, USD_PHASE_COMP - 2, USD_PHASE_COMP_DC - 2 >
  phaseCompFracAndDeriv { phaseCompFrac[0][0], dPhaseCompFrac[0][0] };

  flashModelWrapper.compute( pressure,
                             temperature,
                             compFraction,
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

  // 2) Check derivative with respect to pressure
  real64 const dP = perturbParameter * (pressure + perturbParameter);
  flashModelWrapper.compute( pressure + dP,
                             temperature,
                             compFraction,
                             perturbedPhaseFracAndDeriv,
                             perturbedPhaseCompFracAndDeriv );

  for( integer i = 0; i < numPhase; ++i )
  {
    checkRelativeError( (perturbedPhaseFracAndDeriv.value[i]-phaseFracAndDeriv.value[i])/dP,
                        phaseFracAndDeriv.derivs[i][Deriv::dP],
                        relTol );
    for( integer j = 0; j < numComp; ++j )
    {
      checkRelativeError( (perturbedPhaseCompFracAndDeriv.value[i][j]-phaseCompFracAndDeriv.value[i][j])/dP,
                          phaseCompFracAndDeriv.derivs[i][j][Deriv::dP],
                          relTol );
    }
  }

  // 3) Check derivative with respect to temperature
  real64 const dT = perturbParameter * (temperature + perturbParameter);
  flashModelWrapper.compute( pressure,
                             temperature + dT,
                             compFraction,
                             perturbedPhaseFracAndDeriv,
                             perturbedPhaseCompFracAndDeriv );
  for( integer i = 0; i < numPhase; ++i )
  {
    checkRelativeError( (perturbedPhaseFracAndDeriv.value[i]-phaseFracAndDeriv.value[i])/dT,
                        phaseFracAndDeriv.derivs[i][Deriv::dT],
                        relTol );
    for( integer j = 0; j < numComp; ++j )
    {
      checkRelativeError( (perturbedPhaseCompFracAndDeriv.value[i][j]-phaseCompFracAndDeriv.value[i][j])/dT,
                          phaseCompFracAndDeriv.derivs[i][j][Deriv::dT],
                          relTol );
    }
  }

  // 4) Check derivatives with respect to phaseComposition
  for( integer i = 0; i < numComp; ++i )
  {
    real64 const dC = perturbParameter * (compFraction[i] + perturbParameter);
    stackArray1d< real64, numComp > perturbedCompFraction( numComp );
    for( integer j = 0; j < numComp; ++j )
    {
      perturbedCompFraction[j] = compFraction[j] + ( (i == j) ? dC : 0.0 );
    }
    flashModelWrapper.compute( pressure,
                               temperature,
                               perturbedCompFraction.toSliceConst(),
                               perturbedPhaseFracAndDeriv,
                               perturbedPhaseCompFracAndDeriv );
    for( integer j = 0; j < numPhase; ++j )
    {
      checkRelativeError( (perturbedPhaseFracAndDeriv.value[j]-phaseFracAndDeriv.value[j])/dC,
                          phaseFracAndDeriv.derivs[j][Deriv::dC+i],
                          relTol );
      for( integer k = 0; k < numComp; ++k )
      {
        checkRelativeError( (perturbedPhaseCompFracAndDeriv.value[j][k]-phaseCompFracAndDeriv.value[j][k])/dC,
                            phaseCompFracAndDeriv.derivs[j][k][Deriv::dC+i],
                            relTol );
      }
    }
  }

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

template< typename MODEL >
std::unique_ptr< MODEL > makePVTFunction( string const & filename,
                                          string const & key )
{
  // define component names and molar weight
  string_array componentNames;
  componentNames.resize( 2 );
  componentNames[0] = "co2"; componentNames[1] = "water";

  array1d< real64 > componentMolarWeight;
  componentMolarWeight.resize( 2 );
  componentMolarWeight[0] = 44e-3; componentMolarWeight[1] = 18e-3;

  // read parameters from file
  std::ifstream is( filename );
  std::unique_ptr< MODEL > pvtFunction = nullptr;
  string str;
  while( std::getline( is, str ) )
  {
    string_array const strs = stringutilities::tokenize( str, " " );

    if( strs[0] == key )
    {
      pvtFunction = std::make_unique< MODEL >( strs[1],
                                               strs,
                                               componentNames,
                                               componentMolarWeight );
    }
  }
  GEOSX_ERROR_IF( pvtFunction == nullptr,
                  "Could not find " << key << " in " << filename );

  return pvtFunction;
}

template< typename MODEL >
std::unique_ptr< MODEL > makeFlashModel( string const & filename,
                                         string const & key )
{
  // define phase names
  string_array phaseNames;
  phaseNames.resize( 2 );
  phaseNames[0] = "gas"; phaseNames[1] = "liquid";

  // define component names and molar weight
  string_array componentNames;
  componentNames.resize( 2 );
  componentNames[0] = "co2"; componentNames[1] = "water";

  array1d< real64 > componentMolarWeight;
  componentMolarWeight.resize( 2 );
  componentMolarWeight[0] = 44e-3; componentMolarWeight[1] = 18e-3;

  // read parameters from file
  std::ifstream is( filename );
  std::unique_ptr< MODEL > flashModel;
  string str;
  while( std::getline( is, str ) )
  {
    string_array const strs = stringutilities::tokenize( str, " " );

    if( strs[0] == key )
    {
      flashModel = std::make_unique< MODEL >( strs[1],
                                              strs,
                                              phaseNames,
                                              componentNames,
                                              componentMolarWeight );
    }
  }
  GEOSX_ERROR_IF( flashModel == nullptr,
                  "Could not find " << key << " in " << filename );

  return flashModel;
}

class PhillipsBrineViscosityTest : public ::testing::Test
{
public:
  PhillipsBrineViscosityTest()
  {
    writeTableToFile( filename, pvtLiquidPhillipsTableContent );
    pvtFunction = makePVTFunction< PhillipsBrineViscosity >( filename, key );
  }

  ~PhillipsBrineViscosityTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "ViscosityFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< PhillipsBrineViscosity > pvtFunction;
};

TEST_F( PhillipsBrineViscosityTest, brineViscosityValuesAndDeriv )
{
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const TC[3] = { 94.5, 95, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;


  real64 const savedValues[] = { 0.000303214, 0.000301571, 0.000299598, 0.000303214, 0.000301571,
                                 0.000299598, 0.000303214, 0.000301571, 0.000299598, 0.000303214,
                                 0.000301571, 0.000299598, 0.000303214, 0.000301571, 0.000299598,
                                 0.000303214, 0.000301571, 0.000299598, 0.000303214, 0.000301571,
                                 0.000299598, 0.000303214, 0.000301571, 0.000299598, 0.000303214,
                                 0.000301571, 0.000299598 };

  PhillipsBrineViscosity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

class EzrokhiBrineViscosityTest : public ::testing::Test
{
public:
  EzrokhiBrineViscosityTest()
  {
    writeTableToFile( filename, pvtLiquidEzrokhiTableContent );
    pvtFunction = makePVTFunction< EzrokhiBrineViscosity >( filename, key );
  }

  ~EzrokhiBrineViscosityTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "ViscosityFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< EzrokhiBrineViscosity > pvtFunction;
};

TEST_F( EzrokhiBrineViscosityTest, brineViscosityValuesAndDeriv )
{
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const TC[3] = { 94.5, 95, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  EzrokhiBrineViscosity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

class FenghourCO2ViscosityTest : public ::testing::Test
{
public:
  FenghourCO2ViscosityTest()
  {
    writeTableToFile( filename, pvtGasTableContent );
    pvtFunction = makePVTFunction< FenghourCO2Viscosity >( filename, key );
  }

  ~FenghourCO2ViscosityTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "ViscosityFun";
  string const filename = "pvtgas.txt";
  std::unique_ptr< FenghourCO2Viscosity > pvtFunction;
};

TEST_F( FenghourCO2ViscosityTest, fenghourCO2ViscosityValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  real64 const savedValues[] = { 1.904605302e-05, 1.907031467e-05, 1.909055363e-05, 2.008611673e-05, 2.010303796e-05, 2.011728461e-05,
                                 2.508888392e-05, 2.504245053e-05, 2.500575965e-05, 1.904605302e-05, 1.907031467e-05, 1.909055363e-05,
                                 2.008611673e-05, 2.010303796e-05, 2.011728461e-05, 2.508888392e-05, 2.504245053e-05, 2.500575965e-05,
                                 1.904605302e-05, 1.907031467e-05, 1.909055363e-05, 2.008611673e-05, 2.010303796e-05, 2.011728461e-05,
                                 2.508888392e-05, 2.504245053e-05, 2.500575965e-05 };

  FenghourCO2Viscosity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

class PhillipsBrineDensityTest : public ::testing::Test
{
public:
  PhillipsBrineDensityTest()
  {
    writeTableToFile( filename, pvtLiquidPhillipsTableContent );
    pvtFunction = makePVTFunction< PhillipsBrineDensity >( filename, key );
  }

  ~PhillipsBrineDensityTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "DensityFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< PhillipsBrineDensity > pvtFunction;
};

TEST_F( PhillipsBrineDensityTest, brineCO2DensityMassValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  real64 const savedValues[] = { 1186.281618, 1185.321099, 1184.511961, 1186.997612, 1186.037426, 1185.228539, 1188.470081, 1187.510437,
                                 1186.70195, 1477.603717, 1476.040357, 1474.719028, 1476.803422, 1475.235146, 1473.909633, 1475.091089,
                                 1473.51237, 1472.177974, 2162.60433, 2159.623476, 2157.097807, 2158.238706, 2155.240595, 2152.700314,
                                 2149.037782, 2146.003403, 2143.432409 };

  PhillipsBrineDensity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

TEST_F( PhillipsBrineDensityTest, brineCO2DensityMolarValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  PhillipsBrineDensity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

class EzrokhiBrineDensityTest : public ::testing::Test
{
public:
  EzrokhiBrineDensityTest()
  {
    writeTableToFile( filename, pvtLiquidEzrokhiTableContent );
    pvtFunction = makePVTFunction< EzrokhiBrineDensity >( filename, key );
  }

  ~EzrokhiBrineDensityTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "DensityFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< EzrokhiBrineDensity > pvtFunction;
};

TEST_F( EzrokhiBrineDensityTest, brineCO2DensityMassValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  EzrokhiBrineDensity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

TEST_F( EzrokhiBrineDensityTest, brineCO2DensityMolarValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  EzrokhiBrineDensity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}


class SpanWagnerCO2DensityTest : public ::testing::Test
{
public:
  SpanWagnerCO2DensityTest()
  {
    writeTableToFile( filename, pvtGasTableContent );
    pvtFunction = makePVTFunction< SpanWagnerCO2Density >( filename, key );
  }

  ~SpanWagnerCO2DensityTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "DensityFun";
  string const filename = "pvtgas.txt";
  std::unique_ptr< SpanWagnerCO2Density > pvtFunction;
};

TEST_F( SpanWagnerCO2DensityTest, spanWagnerCO2DensityMassValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTolPrevImpl = 1e-3; // the saved values have been generated with a different pressure spacing, so we loosen the tol a
                                      // little bit
  real64 const relTolDeriv = 5e-5;

  real64 const savedValues[] = { 82.78363562, 82.56888654, 82.39168811, 135.3774839, 134.9199659, 134.5440568, 281.9140962, 280.2559694,
                                 278.9092508, 82.78363562, 82.56888654, 82.39168811, 135.3774839, 134.9199659, 134.5440568, 281.9140962,
                                 280.2559694, 278.9092508, 82.78363562, 82.56888654, 82.39168811, 135.3774839, 134.9199659, 134.5440568,
                                 281.9140962, 280.2559694, 278.9092508 };

  SpanWagnerCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTolPrevImpl );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTolDeriv );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

TEST_F( SpanWagnerCO2DensityTest, spanWagnerCO2DensityMolarValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTolPrevImpl = 1e-3;
  real64 const relTolDeriv = 5e-5;

  real64 const savedValues[] = { 1881.446264, 1876.565603, 1872.538366, 3076.760999, 3066.362862, 3057.819473, 6407.138549, 6369.45385,
                                 6338.846609, 1881.446264, 1876.565603, 1872.538366, 3076.760999, 3066.362862, 3057.819473, 6407.138549,
                                 6369.45385, 6338.846609, 1881.446264, 1876.565603, 1872.538366, 3076.760999, 3066.362862, 3057.819473,
                                 6407.138549, 6369.45385, 6338.846609 };

  SpanWagnerCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], false, relTolPrevImpl );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTolDeriv );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}


class CO2SolubilityTest : public ::testing::Test
{
public:
  CO2SolubilityTest()
  {
    writeTableToFile( filename, co2FlashTableContent );
    flashModel = makeFlashModel< CO2Solubility >( filename, key );
  }

  ~CO2SolubilityTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "FlashModel";
  string const filename = "co2flash.txt";
  std::unique_ptr< CO2Solubility > flashModel;
};

TEST_F( CO2SolubilityTest, co2SolubilityValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTolPrevImpl = 5e-4;
  real64 const relTolDeriv = 5e-5;

  real64 const savedGasPhaseFrac[] = { 0.298158785, 0.298183347, 0.2982033821, 0.295950309, 0.2959791448, 0.2960026365, 0.2926988393,
                                       0.292724834, 0.2927459702, 0.499837295, 0.499854799, 0.4998690769, 0.4982634386, 0.4982839883,
                                       0.4983007295, 0.4959462993, 0.4959648242, 0.4959798868, 0.7015158051, 0.701526251, 0.7015347717,
                                       0.7005765682, 0.7005888317, 0.7005988224, 0.6991937592, 0.6992048145, 0.6992138034 };
  real64 const savedWaterPhaseGasComp[] = { 0.008322701666, 0.008287995083, 0.008259683449, 0.01143341315, 0.0113929227, 0.01135993384,
                                            0.01597786252, 0.01594169644, 0.015912288, 0.008322701666, 0.008287995083, 0.008259683449,
                                            0.01143341315, 0.0113929227, 0.01135993384, 0.01597786252, 0.01594169644, 0.015912288,
                                            0.008322701666, 0.008287995083, 0.008259683449, 0.01143341315, 0.0113929227, 0.01135993384,
                                            0.01597786252, 0.01594169644, 0.015912288 };

  CO2Solubility::KernelWrapper flashModelWrapper = flashModel->createKernelWrapper();

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( flashModelWrapper,
                                                 P[iPres], TC[iTemp], comp,
                                                 savedGasPhaseFrac[counter], savedWaterPhaseGasComp[counter], relTolPrevImpl );
        testNumericalDerivatives( flashModelWrapper, P[iPres], TC[iTemp], comp, eps, relTolDeriv );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

class BrineEnthalpyTest : public ::testing::Test
{
public:
  BrineEnthalpyTest()
  {
    writeTableToFile( filename, pvtLiquidEnthalpyTableContent );
    pvtFunction = makePVTFunction< BrineEnthalpy >( filename, key );
  }

  ~BrineEnthalpyTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "EnthalpyFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< BrineEnthalpy > pvtFunction;
};

TEST_F( BrineEnthalpyTest, BrineEnthalpyMassValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;


  BrineEnthalpy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] =   {    279338.177769, 281327.917595, 282984.982790, 273523.601187, 275547.730627,
                                      277232.967486, 259305.594861, 261448.359022, 263229.980673, 207243.035969,
                                      208857.558651, 210202.000531, 197603.080058, 199274.617099, 200665.764632,
                                      174031.122202, 175899.343122, 177450.286494, 135147.894170, 136387.199707,
                                      137419.018273, 121682.558928, 123001.503570, 124098.561778, 88756.649542,
                                      90350.327221, 91670.592316 };
  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

TEST_F( BrineEnthalpyTest, BrineEnthalpyMolarValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;



  BrineEnthalpy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] =  {  15234891.499346, 15338606.577025, 15424985.891636, 15102742.031582, 15207238.691387, 15294258.271084,
                                   14779605.524173, 14886798.427634, 14976008.570783, 11042832.057960, 11121210.933610, 11186485.534383,
                                   10823742.150877, 10903416.807419, 10969752.900309, 10288015.835964, 10372160.580672, 10442128.397178,
                                   6850772.616574, 6903815.290194, 6947985.177129, 6544742.270173, 6599594.923452, 6645247.529535,
                                   5796426.147754, 5857522.733710, 5908248.223573};
  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], false, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}


class CO2EnthalpyTest : public ::testing::Test
{
public:
  CO2EnthalpyTest()
  {
    // gas enthalpy model repeats liquid parameters (except m), use them here
    writeTableToFile( filename, pvtLiquidEnthalpyTableContent );
    pvtFunction = makePVTFunction< CO2Enthalpy >( filename, key );

  }

  ~CO2EnthalpyTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "EnthalpyFun";
  string const filename = "pvtgas.txt";
  std::unique_ptr< CO2Enthalpy > pvtFunction;
};

TEST_F( CO2EnthalpyTest, CO2EnthalpyMassValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;



  CO2Enthalpy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] =   {     28447.084306, 29131.068469, 29700.204529, 9320.187656, 10117.295548, 10779.101555, -37449.569995,
                                       -36262.216311, -35283.355068, 28447.084306, 29131.068469, 29700.204529, 9320.187656, 10117.295548,
                                       10779.101555, -37449.569995, -36262.216311, -35283.355068, 28447.084306, 29131.068469, 29700.204529,
                                       9320.187656, 10117.295548, 10779.101555, -37449.569995, -36262.216311, -35283.355068};

  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
}

TEST_F( CO2EnthalpyTest, CO2EnthalpyMolarValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.304; comp[1] = 0.696;
  real64 const deltaComp = 0.2;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;



  CO2Enthalpy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] =   {    646524.643323, 662069.737939, 675004.648394, 211822.446731, 229938.535180, 244979.580788,
                                      -851126.590796, -824141.279795, -801894.433361, 646524.643323, 662069.737939, 675004.648394,
                                      211822.446731, 229938.535180, 244979.580788, -851126.590796, -824141.279795, -801894.433361,
                                      646524.643323, 662069.737939, 675004.648394, 211822.446731, 229938.535180, 244979.580788,
                                      -851126.590796, -824141.279795, -801894.433361 };
  integer counter = 0;
  for( integer iComp = 0; iComp < 3; ++iComp )
  {
    for( integer iPres = 0; iPres < 3; ++iPres )
    {
      for( integer iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], false, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
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
