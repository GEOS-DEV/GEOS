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
#include "constitutive/fluid/PVTFunctions/BrineInternalEnergy.hpp"
#include "constitutive/fluid/PVTFunctions/CO2InternalEnergy.hpp"
#include "constitutive/fluid/PVTFunctions/old/BrineEnthalpyFunction.hpp"
#include "constitutive/fluid/PVTFunctions/old/CO2EnthalpyFunction.hpp"
#include "constitutive/fluid/PVTFunctions/old/BrineInternalEnergyFunction.hpp"
#include "constitutive/fluid/PVTFunctions/old/CO2InternalEnergyFunction.hpp"
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
using namespace geosx::PVTProps;

/// Input tables written into temporary files during testing

static const char * pvtLiquidPhillipsTableContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                                    "ViscosityFun PhillipsBrineViscosity 0.1";

// Temperature is specified in Celcius!!!
static const char * pvtLiquidOldEnthalpyTableContent = "EnthalpyFun BrineEnthalpy 1e6 1.5e7 5e4 94 96 1 0.2";

//Temperature is specified in Kelvin, should match the range specified in pvtLiquidOldEnthalpyTableContent
static const char * pvtLiquidEnthalpyTableContent = "EnthalpyFun BrineEnthalpy 1e6 1.5e7 5e4 367.15 369.15 1 0.2";

// Model has no parameters
static const char * pvtLiquidInternalEnergyTableContent = "InternalEnergyFun BrineInternalEnergy";

// the last are set relatively high (1e-4) to increase derivative value and check it properly
static const char * pvtLiquidEzrokhiTableContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                                   "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4";

static const char * pvtGasTableContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                         "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1";

static const char * co2FlashTableContent = "FlashModel CO2Solubility 1e6 1.5e7 5e4 367.15 369.15 1 0.15";

template< typename PVT_WRAPPER >
void testValuesAgainstPreviousImplementation( PVT_WRAPPER const & pvtFunctionWrapper,
                                              real64 const pressure,
                                              real64 const temperature,
                                              arraySlice1d< real64 const > const & phaseComposition,
                                              real64 const & oldImplValue,
                                              bool const useMass,
                                              real64 const relTol )
{
  real64 value = 0.0;
  real64 dValue_dPressure = 0.0;
  real64 dValue_dTemperature = 0.0;
  stackArray1d< real64, 2 > dPhaseComposition_dPressure( 2 );
  stackArray1d< real64, 2 > dPhaseComposition_dTemperature( 2 );
  stackArray2d< real64, 4 > dPhaseComposition_dGlobalCompFraction( 2, 2 );
  stackArray1d< real64, 2 > dValue_dGlobalCompFraction( 2 );
  pvtFunctionWrapper.compute( pressure,
                              temperature,
                              phaseComposition,
                              dPhaseComposition_dPressure.toSliceConst(),
                              dPhaseComposition_dTemperature.toSliceConst(),
                              dPhaseComposition_dGlobalCompFraction.toSliceConst(),
                              value,
                              dValue_dPressure,
                              dValue_dTemperature,
                              dValue_dGlobalCompFraction.toSlice(),
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
  stackArray1d< real64, 2 > phaseFraction( 2 );
  stackArray1d< real64, 2 > dPhaseFraction_dPres( 2 );
  stackArray1d< real64, 2 > dPhaseFraction_dTemp( 2 );
  stackArray2d< real64, 4 > dPhaseFraction_dCompFraction( 2, 2 );
  stackArray2d< real64, 4 > phaseCompFraction( 2, 2 );
  stackArray2d< real64, 4 > dPhaseCompFraction_dPres( 2, 2 );
  stackArray2d< real64, 4 > dPhaseCompFraction_dTemp( 2, 2 );
  stackArray3d< real64, 8 > dPhaseCompFraction_dCompFraction( 2, 2, 2 );
  flashModelWrapper.compute( pressure,
                             temperature,
                             compFraction,
                             phaseFraction.toSlice(),
                             dPhaseFraction_dPres.toSlice(),
                             dPhaseFraction_dTemp.toSlice(),
                             dPhaseFraction_dCompFraction.toSlice(),
                             phaseCompFraction.toSlice(),
                             dPhaseCompFraction_dPres.toSlice(),
                             dPhaseCompFraction_dTemp.toSlice(),
                             dPhaseCompFraction_dCompFraction.toSlice() );

  for( localIndex i = 0; i < 2; ++i )
  {
    real64 const savedPhaseFrac = (i == 0) ? savedGasPhaseFrac : 1.0 - savedGasPhaseFrac;
    checkRelativeError( phaseFraction[i], savedPhaseFrac, relTol );

    for( localIndex j = 0; j < 2; ++j )
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
      checkRelativeError( phaseCompFraction[i][j], savedCompFrac, relTol );
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
  // 1) First compute the unperturbed pressure
  real64 value = 0.0;
  real64 dValue_dPressure = 0.0;
  real64 dValue_dTemperature = 0.0;
  stackArray1d< real64, 2 > dPhaseComposition_dPressure( 2 );
  stackArray1d< real64, 2 > dPhaseComposition_dTemperature( 2 );
  stackArray2d< real64, 4 > dPhaseComposition_dGlobalCompFraction( 2, 2 );
  stackArray1d< real64, 2 > dValue_dGlobalCompFraction( 2 );
  dPhaseComposition_dGlobalCompFraction[0][0] = 1.0;
  dPhaseComposition_dGlobalCompFraction[1][1] = 1.0;
  pvtFunctionWrapper.compute( pressure,
                              temperature,
                              phaseComposition,
                              dPhaseComposition_dPressure.toSliceConst(),
                              dPhaseComposition_dTemperature.toSliceConst(),
                              dPhaseComposition_dGlobalCompFraction.toSliceConst(),
                              value,
                              dValue_dPressure,
                              dValue_dTemperature,
                              dValue_dGlobalCompFraction.toSlice(),
                              useMass );
  real64 perturbedValue = 0.0;
  real64 dPerturbedValue_dPressure = 0.0;
  real64 dPerturbedValue_dTemperature = 0.0;
  stackArray1d< real64, 2 > dPerturbedValue_dGlobalCompFraction( 2 );

  // 2) Check derivative with respect to pressure
  real64 const dP = perturbParameter * (pressure + perturbParameter);
  pvtFunctionWrapper.compute( pressure + dP,
                              temperature,
                              phaseComposition,
                              dPhaseComposition_dPressure.toSliceConst(),
                              dPhaseComposition_dTemperature.toSliceConst(),
                              dPhaseComposition_dGlobalCompFraction.toSliceConst(),
                              perturbedValue,
                              dPerturbedValue_dPressure,
                              dPerturbedValue_dTemperature,
                              dPerturbedValue_dGlobalCompFraction.toSlice(),
                              useMass );
  checkRelativeError( (perturbedValue-value)/dP, dValue_dPressure, relTol );

  // 3) Check derivative with respect to temperature
  real64 const dT = perturbParameter * (temperature + perturbParameter);
  pvtFunctionWrapper.compute( pressure,
                              temperature + dT,
                              phaseComposition,
                              dPhaseComposition_dPressure.toSliceConst(),
                              dPhaseComposition_dTemperature.toSliceConst(),
                              dPhaseComposition_dGlobalCompFraction.toSliceConst(),
                              perturbedValue,
                              dPerturbedValue_dPressure,
                              dPerturbedValue_dTemperature,
                              dPerturbedValue_dGlobalCompFraction.toSlice(),
                              useMass );
  checkRelativeError( (perturbedValue-value)/dT, dValue_dTemperature, relTol );

  // 4) Check derivatives with respect to phaseComposition
  for( localIndex i = 0; i < 2; ++i )
  {
    real64 const dC = perturbParameter * (phaseComposition[i] + perturbParameter);
    stackArray1d< real64, 2 > perturbedPhaseComposition( 2 );
    for( localIndex j = 0; j < 2; ++j )
    {
      perturbedPhaseComposition[j] = phaseComposition[j] + ( (i == j) ? dC : 0.0 );
    }
    pvtFunctionWrapper.compute( pressure,
                                temperature,
                                perturbedPhaseComposition.toSliceConst(),
                                dPhaseComposition_dPressure.toSliceConst(),
                                dPhaseComposition_dTemperature.toSliceConst(),
                                dPhaseComposition_dGlobalCompFraction.toSliceConst(),
                                perturbedValue,
                                dPerturbedValue_dPressure,
                                dPerturbedValue_dTemperature,
                                dPerturbedValue_dGlobalCompFraction.toSlice(),
                                useMass );
    checkRelativeError( (perturbedValue-value)/dC, dValue_dGlobalCompFraction[i], relTol );
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
  // 1) First compute the unperturbed pressure
  stackArray1d< real64, 2 > phaseFraction( 2 );
  stackArray1d< real64, 2 > dPhaseFraction_dPres( 2 );
  stackArray1d< real64, 2 > dPhaseFraction_dTemp( 2 );
  stackArray2d< real64, 4 > dPhaseFraction_dCompFraction( 2, 2 );
  stackArray2d< real64, 4 > phaseCompFraction( 2, 2 );
  stackArray2d< real64, 4 > dPhaseCompFraction_dPres( 2, 2 );
  stackArray2d< real64, 4 > dPhaseCompFraction_dTemp( 2, 2 );
  stackArray3d< real64, 8 > dPhaseCompFraction_dCompFraction( 2, 2, 2 );
  flashModelWrapper.compute( pressure,
                             temperature,
                             compFraction,
                             phaseFraction.toSlice(),
                             dPhaseFraction_dPres.toSlice(),
                             dPhaseFraction_dTemp.toSlice(),
                             dPhaseFraction_dCompFraction.toSlice(),
                             phaseCompFraction.toSlice(),
                             dPhaseCompFraction_dPres.toSlice(),
                             dPhaseCompFraction_dTemp.toSlice(),
                             dPhaseCompFraction_dCompFraction.toSlice() );
  stackArray1d< real64, 2 > perturbedPhaseFraction( 2 );
  stackArray1d< real64, 2 > dPerturbedPhaseFraction_dPres( 2 );
  stackArray1d< real64, 2 > dPerturbedPhaseFraction_dTemp( 2 );
  stackArray2d< real64, 4 > dPerturbedPhaseFraction_dCompFraction( 2, 2 );
  stackArray2d< real64, 4 > perturbedPhaseCompFraction( 2, 2 );
  stackArray2d< real64, 4 > dPerturbedPhaseCompFraction_dPres( 2, 2 );
  stackArray2d< real64, 4 > dPerturbedPhaseCompFraction_dTemp( 2, 2 );
  stackArray3d< real64, 8 > dPerturbedPhaseCompFraction_dCompFraction( 2, 2, 2 );

  // 2) Check derivative with respect to pressure
  real64 const dP = perturbParameter * (pressure + perturbParameter);
  flashModelWrapper.compute( pressure + dP,
                             temperature,
                             compFraction,
                             perturbedPhaseFraction.toSlice(),
                             dPerturbedPhaseFraction_dPres.toSlice(),
                             dPerturbedPhaseFraction_dTemp.toSlice(),
                             dPerturbedPhaseFraction_dCompFraction.toSlice(),
                             perturbedPhaseCompFraction.toSlice(),
                             dPerturbedPhaseCompFraction_dPres.toSlice(),
                             dPerturbedPhaseCompFraction_dTemp.toSlice(),
                             dPerturbedPhaseCompFraction_dCompFraction.toSlice() );
  for( localIndex i = 0; i < 2; ++i )
  {
    checkRelativeError( (perturbedPhaseFraction[i]-phaseFraction[i])/dP, dPhaseFraction_dPres[i], relTol );
    for( localIndex j = 0; j < 2; ++j )
    {
      checkRelativeError( (perturbedPhaseCompFraction[i][j]-phaseCompFraction[i][j])/dP, dPhaseCompFraction_dPres[i][j], relTol );
    }
  }

  // 3) Check derivative with respect to temperature
  real64 const dT = perturbParameter * (temperature + perturbParameter);
  flashModelWrapper.compute( pressure,
                             temperature + dT,
                             compFraction,
                             perturbedPhaseFraction.toSlice(),
                             dPerturbedPhaseFraction_dPres.toSlice(),
                             dPerturbedPhaseFraction_dTemp.toSlice(),
                             dPerturbedPhaseFraction_dCompFraction.toSlice(),
                             perturbedPhaseCompFraction.toSlice(),
                             dPerturbedPhaseCompFraction_dPres.toSlice(),
                             dPerturbedPhaseCompFraction_dTemp.toSlice(),
                             dPerturbedPhaseCompFraction_dCompFraction.toSlice() );
  for( localIndex i = 0; i < 2; ++i )
  {
    checkRelativeError( (perturbedPhaseFraction[i]-phaseFraction[i])/dT, dPhaseFraction_dTemp[i], relTol );
    for( localIndex j = 0; j < 2; ++j )
    {
      checkRelativeError( (perturbedPhaseCompFraction[i][j]-phaseCompFraction[i][j])/dT, dPhaseCompFraction_dTemp[i][j], relTol );
    }
  }

  // 4) Check derivatives with respect to phaseComposition
  for( localIndex i = 0; i < 2; ++i )
  {
    real64 const dC = perturbParameter * (compFraction[i] + perturbParameter);
    stackArray1d< real64, 2 > perturbedCompFraction( 2 );
    for( localIndex j = 0; j < 2; ++j )
    {
      perturbedCompFraction[j] = compFraction[j] + ( (i == j) ? dC : 0.0 );
    }
    flashModelWrapper.compute( pressure,
                               temperature,
                               perturbedCompFraction.toSliceConst(),
                               perturbedPhaseFraction.toSlice(),
                               dPerturbedPhaseFraction_dPres.toSlice(),
                               dPerturbedPhaseFraction_dTemp.toSlice(),
                               dPerturbedPhaseFraction_dCompFraction.toSlice(),
                               perturbedPhaseCompFraction.toSlice(),
                               dPerturbedPhaseCompFraction_dPres.toSlice(),
                               dPerturbedPhaseCompFraction_dTemp.toSlice(),
                               dPerturbedPhaseCompFraction_dCompFraction.toSlice() );
    for( localIndex j = 0; j < 2; ++j )
    {
      checkRelativeError( (perturbedPhaseFraction[j]-phaseFraction[j])/dC, dPhaseFraction_dCompFraction[j][i], relTol );
      for( localIndex k = 0; k < 2; ++k )
      {
        checkRelativeError( (perturbedPhaseCompFraction[j][k]-phaseCompFraction[j][k])/dC, dPhaseCompFraction_dCompFraction[j][k][i], relTol );
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

template< typename MODEL >
std::unique_ptr< MODEL > makeOldPVTFunction( string const & filename,
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
      pvtFunction = std::make_unique< MODEL >( strs,
                                               componentNames,
                                               componentMolarWeight );
    }
  }
  GEOSX_ERROR_IF( pvtFunction == nullptr,
                  "Could not find " << key << " in " << filename );

  return pvtFunction;
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

  real64 const savedValues[] = { 0.0009009475991, 0.0009009665224, 0.0009009892304, 0.0009009475991, 0.0009009665224,
                                 0.0009009892304, 0.0009009475991, 0.0009009665224, 0.0009009892304, 0.0009009475991,
                                 0.0009009665224, 0.0009009892304, 0.0009009475991, 0.0009009665224, 0.0009009892304,
                                 0.0009009475991, 0.0009009665224, 0.0009009892304, 0.0009009475991, 0.0009009665224,
                                 0.0009009892304, 0.0009009475991, 0.0009009665224, 0.0009009892304, 0.0009009475991,
                                 0.0009009665224, 0.0009009892304 };

  PhillipsBrineViscosity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        //testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
        //                                         P[iPres], TC[iTemp], comp, savedValues[counter], false, relTol );
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

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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
  real64 const relTol = 5e-5;

  real64 const savedValues[] = { 82.78363562, 82.56888654, 82.39168811, 135.3774839, 134.9199659, 134.5440568, 281.9140962, 280.2559694,
                                 278.9092508, 82.78363562, 82.56888654, 82.39168811, 135.3774839, 134.9199659, 134.5440568, 281.9140962,
                                 280.2559694, 278.9092508, 82.78363562, 82.56888654, 82.39168811, 135.3774839, 134.9199659, 134.5440568,
                                 281.9140962, 280.2559694, 278.9092508 };

  SpanWagnerCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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
  real64 const relTol = 5e-5;

  real64 const savedValues[] = { 1881.446264, 1876.565603, 1872.538366, 3076.760999, 3066.362862, 3057.819473, 6407.138549, 6369.45385,
                                 6338.846609, 1881.446264, 1876.565603, 1872.538366, 3076.760999, 3066.362862, 3057.819473, 6407.138549,
                                 6369.45385, 6338.846609, 1881.446264, 1876.565603, 1872.538366, 3076.760999, 3066.362862, 3057.819473,
                                 6407.138549, 6369.45385, 6338.846609 };

  SpanWagnerCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
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
  real64 const relTol = 5e-5;

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

  localIndex counter = 0;
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        testValuesAgainstPreviousImplementation( flashModelWrapper,
                                                 P[iPres], TC[iTemp], comp,
                                                 savedGasPhaseFrac[counter], savedWaterPhaseGasComp[counter], relTol );
        testNumericalDerivatives( flashModelWrapper, P[iPres], TC[iTemp], comp, eps, relTol );
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
    writeTableToFile( filename, pvtLiquidOldEnthalpyTableContent );
    pvtFunctionOld = makeOldPVTFunction< BrineEnthalpyFunction >( filename, key );
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
  std::unique_ptr< BrineEnthalpyFunction > pvtFunctionOld;
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
  localIndex counter = 0;
  printf ( " Mass old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, true );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, true );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
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
  localIndex counter = 0;
  printf ( " Molar old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, false );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, false );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], false, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
}


class CO2EnthalpyTest : public ::testing::Test
{
public:
  CO2EnthalpyTest()
  {
    // gas enthalpy model repeats liquid parameters (except m), use them here
    writeTableToFile( filename, pvtLiquidOldEnthalpyTableContent );
    pvtFunctionOld = makeOldPVTFunction< CO2EnthalpyFunction >( filename, key );
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
  std::unique_ptr< CO2EnthalpyFunction > pvtFunctionOld;
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

  localIndex counter = 0;
  printf ( " CO2Enthalpy Mass old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, true );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, true );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
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
  localIndex counter = 0;
  printf ( " CO2Enthalpy Molar old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, false );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, false );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], false, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
}

class BrineInternalEnergyTest : public ::testing::Test
{
public:
  BrineInternalEnergyTest()
  {
    writeTableToFile( filename, pvtLiquidInternalEnergyTableContent );

    pvtFunctionOld = makeOldPVTFunction< BrineInternalEnergyFunction >( filename, key );
    writeTableToFile( filename, pvtLiquidInternalEnergyTableContent );
    pvtFunction = makePVTFunction< BrineInternalEnergy >( filename, key );

  }

  ~BrineInternalEnergyTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "InternalEnergyFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< BrineInternalEnergyFunction > pvtFunctionOld;
  std::unique_ptr< BrineInternalEnergy > pvtFunction;
};

TEST_F( BrineInternalEnergyTest, BrineInternalEnergyMassValuesAndDeriv )
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



  BrineInternalEnergy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] = {      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000, 5106.500000, 5107.100000, 5107.600000,
                                      7640.500000, 7641.100000, 7641.600000, 12984.500000, 12985.100000, 12985.600000,
                                      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000};


  localIndex counter = 0;
  printf ( " BrineInternalEnergy Mass old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, true );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, true );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
}

TEST_F( BrineInternalEnergyTest, BrineInternalEnergyMolarValuesAndDeriv )
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



  BrineInternalEnergy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] = {      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000, 5106.500000, 5107.100000, 5107.600000,
                                      7640.500000, 7641.100000, 7641.600000, 12984.500000, 12985.100000, 12985.600000,
                                      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000};
  localIndex counter = 0;
  printf ( " BrineInternalEnergy Molar old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, false );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, false );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], false, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
}

class CO2InternalEnergyTest : public ::testing::Test
{
public:
  CO2InternalEnergyTest()
  {
    writeTableToFile( filename, pvtLiquidInternalEnergyTableContent );

    pvtFunctionOld = makeOldPVTFunction< CO2InternalEnergyFunction >( filename, key );
    writeTableToFile( filename, pvtLiquidInternalEnergyTableContent );
    pvtFunction = makePVTFunction< CO2InternalEnergy >( filename, key );

  }

  ~CO2InternalEnergyTest() override
  {
    removeFile( filename );
  }

protected:
  string const key = "InternalEnergyFun";
  string const filename = "pvtgas.txt";
  std::unique_ptr< CO2InternalEnergyFunction > pvtFunctionOld;
  std::unique_ptr< CO2InternalEnergy > pvtFunction;
};

TEST_F( CO2InternalEnergyTest, CO2InternalEnergyMassValuesAndDeriv )
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



  CO2InternalEnergy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] = {      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000, 5106.500000, 5107.100000, 5107.600000,
                                      7640.500000, 7641.100000, 7641.600000, 12984.500000, 12985.100000, 12985.600000,
                                      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000};


  localIndex counter = 0;
  printf ( " CO2InternalEnergy Mass old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, true );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, true );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], true, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
}

TEST_F( CO2InternalEnergyTest, CO2InternalEnergyMolarValuesAndDeriv )
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



  CO2InternalEnergy::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();
  real64 const savedValues[] = {      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000, 5106.500000, 5107.100000, 5107.600000,
                                      7640.500000, 7641.100000, 7641.600000, 12984.500000, 12985.100000, 12985.600000,
                                      5106.500000, 5107.100000, 5107.600000, 7640.500000, 7641.100000, 7641.600000,
                                      12984.500000, 12985.100000, 12985.600000};
  localIndex counter = 0;
  printf ( " CO2InternalEnergy Molar old values: \n { " );
  for( localIndex iComp = 0; iComp < 3; ++iComp )
  {
    for( localIndex iPres = 0; iPres < 3; ++iPres )
    {
      for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
      {
        EvalVarArgs const pressure_old = P[iPres];
        EvalVarArgs const temp_old = TC[iTemp];
        array1d< EvalVarArgs > comp_old ( 2 );
        comp_old[0] = comp[0];
        comp_old[1] = comp[1];
        EvalVarArgs value_old;
        real64 value;

        pvtFunctionOld->evaluation( pressure_old, temp_old, comp_old, value_old, false );
        pvtFunctionWrapper.compute( pressure_old.m_var, temp_old.m_var, comp.toSliceConst(), value, false );
        checkRelativeError ( value, value_old.m_var, relTol );
        printf ( " %15lf, ", value_old.m_var );

        testValuesAgainstPreviousImplementation( pvtFunctionWrapper,
                                                 P[iPres], TC[iTemp], comp, savedValues[counter], false, relTol );
        testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
        counter++;
      }
    }
    comp[0] += deltaComp;
    comp[1] = 1 - comp[0];
  }
  printf ( " } \n " );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
