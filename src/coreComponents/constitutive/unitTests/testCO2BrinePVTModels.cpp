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
#include "constitutive/fluid/PVTFunctions/NewBrineViscosityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/NewFenghourCO2ViscosityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/NewBrineCO2DensityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/NewSpanWagnerCO2DensityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/NewCO2SolubilityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/BrineViscosityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/FenghourCO2ViscosityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/BrineCO2DensityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/CO2SolubilityFunction.hpp"
#include "managers/GeosxState.hpp"
#include "managers/initialization.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;
using namespace geosx::stringutilities;
using namespace geosx::PVTProps;

/// Input tables written into temporary files during testing

static const char * pvtliquid_str = "DensityFun BrineCO2Density 1e6 1.5e7 5e4 94 96 1 0.2\n"
                                    "ViscosityFun BrineViscosity 0.1";

static const char * pvtgas_str = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 94 96 1\n"
                                 "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 94 96 1";

static const char * co2flash_str = "FlashModel CO2Solubility 1e6 1.5e7 5e4 94 96 1 0.15";

// TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
void testValuesAgainstPreviousImplementation( PVTFunctionBaseUpdate const & pvtFunctionWrapper,
                                              PVTFunction const & oldPvtFunction,
                                              real64 const pressure,
                                              real64 const temperature,
                                              arraySlice1d< real64 > const & phaseComposition,
                                              bool const useMass,
                                              real64 const relTol )
{
  real64 value = 0.0;
  real64 dValue_dPressure = 0.0;
  real64 dValue_dTemperature = 0.0;
  stackArray1d< real64, 2 > dValue_dPhaseComposition( 2 );
  pvtFunctionWrapper.compute( pressure,
                              temperature,
                              phaseComposition,
                              value,
                              dValue_dPressure,
                              dValue_dTemperature,
                              dValue_dPhaseComposition,
                              useMass );

  EvalVarArgs oldPres( pressure );
  EvalVarArgs oldTemp( temperature );
  EvalVarArgs oldValue( 0.0 );
  stackArray1d< EvalVarArgs, 2 > oldPhaseComposition( 2 );
  for( localIndex ic = 0; ic < 2; ic++ )
  {
    oldPhaseComposition[ic].m_var = phaseComposition[ic];
  }
  oldPvtFunction.evaluation( oldPres,
                             oldTemp,
                             oldPhaseComposition,
                             oldValue,
                             useMass );

  checkRelativeError( value, oldValue.m_var, relTol );
}

// TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
void testValuesAgainstPreviousImplementation( FlashModelBaseUpdate const & flashModelWrapper,
                                              FlashModel const & oldFlashModel,
                                              real64 const pressure,
                                              real64 const temperature,
                                              arraySlice1d< real64 > const & compFraction,
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
                             phaseFraction,
                             dPhaseFraction_dPres,
                             dPhaseFraction_dTemp,
                             dPhaseFraction_dCompFraction,
                             phaseCompFraction,
                             dPhaseCompFraction_dPres,
                             dPhaseCompFraction_dTemp,
                             dPhaseCompFraction_dCompFraction );

  EvalVarArgs oldPres( pressure );
  EvalVarArgs oldTemp( temperature );
  stackArray1d< EvalVarArgs, 2 > oldCompFraction( 2 );
  stackArray1d< EvalVarArgs, 2 > oldPhaseFraction( 2 );
  stackArray2d< EvalVarArgs, 4 > oldPhaseCompFraction( 2, 2 );
  for( localIndex ic = 0; ic < 2; ic++ )
  {
    oldCompFraction[ic].m_var = compFraction[ic];
  }
  oldFlashModel.partition( oldPres,
                           oldTemp,
                           oldCompFraction,
                           oldPhaseFraction,
                           oldPhaseCompFraction );

  for( localIndex i = 0; i < 2; ++i )
  {
    checkRelativeError( phaseFraction[i], oldPhaseFraction[i].m_var, relTol );
    for( localIndex j = 0; j < 2; ++j )
    {
      checkRelativeError( phaseCompFraction[i][j], oldPhaseCompFraction[i][j].m_var, relTol );
    }
  }
}

void testNumericalDerivatives( PVTFunctionBaseUpdate const & pvtFunctionWrapper,
                               real64 const pressure,
                               real64 const temperature,
                               arraySlice1d< real64 > const & phaseComposition,
                               bool const useMass,
                               real64 const perturbParameter,
                               real64 const relTol )
{
  // 1) First compute the unperturbed pressure
  real64 value = 0.0;
  real64 dValue_dPressure = 0.0;
  real64 dValue_dTemperature = 0.0;
  stackArray1d< real64, 2 > dValue_dPhaseComposition( 2 );
  pvtFunctionWrapper.compute( pressure,
                              temperature,
                              phaseComposition,
                              value,
                              dValue_dPressure,
                              dValue_dTemperature,
                              dValue_dPhaseComposition,
                              useMass );
  real64 perturbedValue = 0.0;
  real64 dPerturbedValue_dPressure = 0.0;
  real64 dPerturbedValue_dTemperature = 0.0;
  stackArray1d< real64, 2 > dPerturbedValue_dPhaseComposition( 2 );

  // 2) Check derivative with respect to pressure
  real64 const dP = perturbParameter * (pressure + perturbParameter);
  pvtFunctionWrapper.compute( pressure + dP,
                              temperature,
                              phaseComposition,
                              perturbedValue,
                              dPerturbedValue_dPressure,
                              dPerturbedValue_dTemperature,
                              dPerturbedValue_dPhaseComposition,
                              useMass );
  checkRelativeError( (perturbedValue-value)/dP, dValue_dPressure, relTol );

  // 3) Check derivative with respect to temperature
  real64 const dT = perturbParameter * (temperature + perturbParameter);
  pvtFunctionWrapper.compute( pressure,
                              temperature + dT,
                              phaseComposition,
                              perturbedValue,
                              dPerturbedValue_dPressure,
                              dPerturbedValue_dTemperature,
                              dPerturbedValue_dPhaseComposition,
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
      perturbedPhaseComposition[j] /= ( 1 + dC );
    }
    pvtFunctionWrapper.compute( pressure,
                                temperature,
                                perturbedPhaseComposition,
                                perturbedValue,
                                dPerturbedValue_dPressure,
                                dPerturbedValue_dTemperature,
                                dPerturbedValue_dPhaseComposition,
                                useMass );
    checkRelativeError( (perturbedValue-value)/dC, dValue_dPhaseComposition[i], relTol );
  }
}

void testNumericalDerivatives( FlashModelBaseUpdate const & flashModelWrapper,
                               real64 const pressure,
                               real64 const temperature,
                               arraySlice1d< real64 > const & compFraction,
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
                             phaseFraction,
                             dPhaseFraction_dPres,
                             dPhaseFraction_dTemp,
                             dPhaseFraction_dCompFraction,
                             phaseCompFraction,
                             dPhaseCompFraction_dPres,
                             dPhaseCompFraction_dTemp,
                             dPhaseCompFraction_dCompFraction );
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
                             perturbedPhaseFraction,
                             dPerturbedPhaseFraction_dPres,
                             dPerturbedPhaseFraction_dTemp,
                             dPerturbedPhaseFraction_dCompFraction,
                             perturbedPhaseCompFraction,
                             dPerturbedPhaseCompFraction_dPres,
                             dPerturbedPhaseCompFraction_dTemp,
                             dPerturbedPhaseCompFraction_dCompFraction );
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
                             perturbedPhaseFraction,
                             dPerturbedPhaseFraction_dPres,
                             dPerturbedPhaseFraction_dTemp,
                             dPerturbedPhaseFraction_dCompFraction,
                             perturbedPhaseCompFraction,
                             dPerturbedPhaseCompFraction_dPres,
                             dPerturbedPhaseCompFraction_dTemp,
                             dPerturbedPhaseCompFraction_dCompFraction );
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
      perturbedCompFraction[j] /= ( 1 + dC );
    }
    flashModelWrapper.compute( pressure,
                               temperature,
                               perturbedCompFraction,
                               perturbedPhaseFraction,
                               dPerturbedPhaseFraction_dPres,
                               dPerturbedPhaseFraction_dTemp,
                               dPerturbedPhaseFraction_dCompFraction,
                               perturbedPhaseCompFraction,
                               dPerturbedPhaseCompFraction_dPres,
                               dPerturbedPhaseCompFraction_dTemp,
                               dPerturbedPhaseCompFraction_dCompFraction );
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
  array1d< string > componentNames;
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
    string_array const strs = Tokenize( str, " " );

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
  array1d< string > phaseNames;
  phaseNames.resize( 2 );
  phaseNames[0] = "gas"; phaseNames[1] = "liquid";

  // define component names and molar weight
  array1d< string > componentNames;
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
    string_array const strs = Tokenize( str, " " );

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
protected:

  virtual void SetUp() override
  {
    writeTableToFile( filename, pvtliquid_str );

    pvtFunction = makePVTFunction< BrineViscosity >( filename, key );
    oldPvtFunction = makePVTFunction< BrineViscosityFunction >( filename, key ); // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
  }

  virtual void TearDown() override
  {
    removeFile( filename );
  }

  string const key = "ViscosityFun";
  string const filename = "pvtliquid.txt";

  std::unique_ptr< BrineViscosity > pvtFunction;
  std::unique_ptr< BrineViscosityFunction > oldPvtFunction; // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
};

TEST_F( BrineViscosityTest, brineViscosityValuesAndDeriv )
{
  real64 const P[3] = { 5e6, 7.5e6, 1.2e7 };
  real64 const TC[3] = { 94.5, 95, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  BrineViscosity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  for( localIndex iPres = 0; iPres < 3; ++iPres )
  {
    for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
    {
      testValuesAgainstPreviousImplementation( pvtFunctionWrapper, // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
                                               *oldPvtFunction, P[iPres], TC[iTemp], comp, true, relTol );
      testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
    }
  }
}

class FenghourCO2ViscosityTest : public ::testing::Test
{
protected:

  virtual void SetUp() override
  {
    writeTableToFile( filename, pvtgas_str );

    pvtFunction = makePVTFunction< FenghourCO2Viscosity >( filename, key );
    oldPvtFunction = makePVTFunction< FenghourCO2ViscosityFunction >( filename, key ); // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
  }

  virtual void TearDown() override
  {
    removeFile( filename );
  }

  string const key = "ViscosityFun";
  string const filename = "pvtgas.txt";

  std::unique_ptr< FenghourCO2Viscosity > pvtFunction;
  std::unique_ptr< FenghourCO2ViscosityFunction > oldPvtFunction; // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
};

TEST_F( FenghourCO2ViscosityTest, fenghourCO2ViscosityValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  FenghourCO2Viscosity::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  for( localIndex iPres = 0; iPres < 3; ++iPres )
  {
    for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
    {
      testValuesAgainstPreviousImplementation( pvtFunctionWrapper, // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
                                               *oldPvtFunction, P[iPres], TC[iTemp], comp, true, relTol );
      testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
    }
  }
}

class BrineCO2DensityTest : public ::testing::Test
{
protected:

  virtual void SetUp() override
  {
    writeTableToFile( filename, pvtliquid_str );

    pvtFunction = makePVTFunction< BrineCO2Density >( filename, key );
    oldPvtFunction = makePVTFunction< BrineCO2DensityFunction >( filename, key ); // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
  }

  virtual void TearDown() override
  {
    removeFile( filename );
  }

  string const key = "DensityFun";
  string const filename = "pvtliquid.txt";

  std::unique_ptr< BrineCO2Density > pvtFunction;
  std::unique_ptr< BrineCO2DensityFunction > oldPvtFunction; // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
};

TEST_F( BrineCO2DensityTest, brineCO2DensityMassValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  BrineCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  for( localIndex iPres = 0; iPres < 3; ++iPres )
  {
    for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
    {
      testValuesAgainstPreviousImplementation( pvtFunctionWrapper, // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
                                               *oldPvtFunction, P[iPres], TC[iTemp], comp, true, relTol );
      testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
    }
  }
}

TEST_F( BrineCO2DensityTest, brineCO2DensityMolarValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  BrineCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  for( localIndex iPres = 0; iPres < 3; ++iPres )
  {
    for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
    {
      testValuesAgainstPreviousImplementation( pvtFunctionWrapper, // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
                                               *oldPvtFunction, P[iPres], TC[iTemp], comp, false, relTol );
      testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
    }
  }
}


class SpanWagnerCO2DensityTest : public ::testing::Test
{
protected:

  virtual void SetUp() override
  {
    writeTableToFile( filename, pvtgas_str );

    pvtFunction = makePVTFunction< SpanWagnerCO2Density >( filename, key );
    oldPvtFunction = makePVTFunction< SpanWagnerCO2DensityFunction >( filename, key ); // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
  }

  virtual void TearDown() override
  {
    removeFile( filename );
  }

  string const key = "DensityFun";
  string const filename = "pvtgas.txt";

  std::unique_ptr< SpanWagnerCO2Density > pvtFunction;
  std::unique_ptr< SpanWagnerCO2DensityFunction > oldPvtFunction; // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
};

TEST_F( SpanWagnerCO2DensityTest, spanWagnerCO2DensityMassValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  SpanWagnerCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  for( localIndex iPres = 0; iPres < 3; ++iPres )
  {
    for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
    {
      testValuesAgainstPreviousImplementation( pvtFunctionWrapper, // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
                                               *oldPvtFunction, P[iPres], TC[iTemp], comp, true, relTol );
      testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, true, eps, relTol );
    }
  }
}

TEST_F( SpanWagnerCO2DensityTest, spanWagnerCO2DensityMolarValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  SpanWagnerCO2Density::KernelWrapper pvtFunctionWrapper = pvtFunction->createKernelWrapper();

  for( localIndex iPres = 0; iPres < 3; ++iPres )
  {
    for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
    {
      testValuesAgainstPreviousImplementation( pvtFunctionWrapper, // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
                                               *oldPvtFunction, P[iPres], TC[iTemp], comp, false, relTol );
      testNumericalDerivatives( pvtFunctionWrapper, P[iPres], TC[iTemp], comp, false, eps, relTol );
    }
  }
}


class CO2SolubilityTest : public ::testing::Test
{
protected:

  virtual void SetUp() override
  {
    writeTableToFile( filename, co2flash_str );

    flashModel = makeFlashModel< CO2Solubility >( filename, key );
    oldFlashModel = makeFlashModel< CO2SolubilityFunction >( filename, key ); // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
  }

  virtual void TearDown() override
  {
    removeFile( filename );
  }

  string const key = "FlashModel";
  string const filename = "co2flash.txt";

  std::unique_ptr< CO2Solubility > flashModel;
  std::unique_ptr< CO2SolubilityFunction > oldFlashModel; // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
};

TEST_F( CO2SolubilityTest, co2SolubilityValuesAndDeriv )
{
  // when checking numerical derivatives, do not fall on the coordinate points of the tables!!
  // (see the txt file defined at the top of the file for the definition of the coordinates)
  real64 const P[3] = { 5.012e6, 7.546e6, 1.289e7 };
  real64 const TC[3] = { 94.5, 95.1, 95.6 };
  array1d< real64 > comp( 2 );
  comp[0] = 0.3; comp[1] = 0.7;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const relTol = 5e-5;

  CO2Solubility::KernelWrapper flashModelWrapper = flashModel->createKernelWrapper();

  for( localIndex iPres = 0; iPres < 3; ++iPres )
  {
    for( localIndex iTemp = 0; iTemp < 3; ++iTemp )
    {
      testValuesAgainstPreviousImplementation( flashModelWrapper, // TEMPORARY CODE TO VALIDATE NEW IMPLEMENTATION
                                               *oldFlashModel, P[iPres], TC[iTemp], comp, relTol );
      testNumericalDerivatives( flashModelWrapper, P[iPres], TC[iTemp], comp, eps, relTol );
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
