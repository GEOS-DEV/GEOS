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
#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/fluid/PVTFunctions/BrineViscosity.hpp"
#include "constitutive/fluid/PVTFunctions/FenghourCO2Viscosity.hpp"
#include "constitutive/fluid/PVTFunctions/BrineCO2Density.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2Density.hpp"
#include "constitutive/fluid/PVTFunctions/CO2Solubility.hpp"
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

/// Input tables written into temporary files during testing

static const char * pvtLiquidTableContent = "DensityFun BrineCO2Density 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                            "ViscosityFun BrineViscosity 0.1";

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

  constexpr std::streamsize buf_size = 256;
  char buf[buf_size];

  std::unique_ptr< MODEL > pvtFunction = nullptr;

  while( is.getline( buf, buf_size ))
  {
    string const str( buf );
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

  constexpr std::streamsize buf_size = 256;
  char buf[buf_size];

  std::unique_ptr< MODEL > flashModel = nullptr;

  while( is.getline( buf, buf_size ))
  {
    string const str( buf );
    string_array const strs = stringutilities::tokenize( str, " " );

    if( strs[0] == key )
    {
      flashModel = std::make_unique< MODEL >( strs,
                                              phaseNames,
                                              componentNames,
                                              componentMolarWeight );
    }
  }
  GEOSX_ERROR_IF( flashModel == nullptr,
                  "Could not find " << key << " in " << filename );

  return flashModel;
}


class BrineViscosityTest : public ::testing::Test
{
public:
  BrineViscosityTest()
  {
    writeTableToFile( filename, pvtLiquidTableContent );
    pvtFunction = makePVTFunction< BrineViscosity >( filename, key );
  }

  ~BrineViscosityTest()
  {
    removeFile( filename );
  }

protected:
  string const key = "ViscosityFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< BrineViscosity > pvtFunction;
};

TEST_F( BrineViscosityTest, brineViscosityValuesAndDeriv )
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

  BrineViscosity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

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

class FenghourCO2ViscosityTest : public ::testing::Test
{
public:
  FenghourCO2ViscosityTest()
  {
    writeTableToFile( filename, pvtGasTableContent );
    pvtFunction = makePVTFunction< FenghourCO2Viscosity >( filename, key );
  }

  ~FenghourCO2ViscosityTest()
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

class BrineCO2DensityTest : public ::testing::Test
{
public:
  BrineCO2DensityTest()
  {
    writeTableToFile( filename, pvtLiquidTableContent );
    pvtFunction = makePVTFunction< BrineCO2Density >( filename, key );
  }

  ~BrineCO2DensityTest()
  {
    removeFile( filename );
  }

protected:
  string const key = "DensityFun";
  string const filename = "pvtliquid.txt";
  std::unique_ptr< BrineCO2Density > pvtFunction;
};

TEST_F( BrineCO2DensityTest, brineCO2DensityMassValuesAndDeriv )
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

  BrineCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

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

TEST_F( BrineCO2DensityTest, brineCO2DensityMolarValuesAndDeriv )
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

  BrineCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

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


class SpanWagnerCO2DensityTest : public ::testing::Test
{
public:
  SpanWagnerCO2DensityTest()
  {
    writeTableToFile( filename, pvtGasTableContent );
    pvtFunction = makePVTFunction< SpanWagnerCO2Density >( filename, key );
  }

  ~SpanWagnerCO2DensityTest()
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

  ~CO2SolubilityTest()
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

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
