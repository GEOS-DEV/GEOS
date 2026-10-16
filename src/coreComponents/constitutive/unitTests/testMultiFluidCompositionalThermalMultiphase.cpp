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
 * @file testMultiFluidCompositionalThermalMultiphase.cpp
 */

#include "FluidModelTest.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/HeatCapacityCoefficients.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/KValueFlashParameters.hpp"
#include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

enum class FLASH_TYPE : int { NEGATIVE, KVALUE };
ENUM_STRINGS( FLASH_TYPE, "Negative", "KValue" );

enum class VISCOSITY_TYPE : int { LBC, PHILLIPS };
ENUM_STRINGS( VISCOSITY_TYPE, "LohrenzBrayClark", "Phillips" );

template< FLASH_TYPE FLASH, VISCOSITY_TYPE VISCOSITY >
struct FluidModel {};

template<>
struct FluidModel< FLASH_TYPE::NEGATIVE, VISCOSITY_TYPE::LBC >
{
  using type = CompositionalThermalTwoPhaseLohrenzBrayClarkViscosity;
};
template<>
struct FluidModel< FLASH_TYPE::NEGATIVE, VISCOSITY_TYPE::PHILLIPS >
{
  using type = CompositionalThermalTwoPhasePhillipsBrine;
};
template<>
struct FluidModel< FLASH_TYPE::KVALUE, VISCOSITY_TYPE::PHILLIPS >
{
  using type = CompositionalThermalKValuePhillipsBrine;
};

template< typename FluidModel, integer NUM_COMP >
struct FluidData
{};

template< typename TEST_TYPE >
class MultiThermalFluidCompositionalMultiphaseTestFixture : public FluidModelTest<
    typename FluidModel< std::tuple_element_t< 1, TEST_TYPE >::value,
                         std::tuple_element_t< 2, TEST_TYPE >::value >::type,
    std::tuple_element_t< 3, TEST_TYPE >::value >
{
public:
  static constexpr EquationOfStateType EquationOfState = std::tuple_element_t< 0, TEST_TYPE >::value;
  static constexpr FLASH_TYPE FlashModel = std::tuple_element_t< 1, TEST_TYPE >::value;
  static constexpr VISCOSITY_TYPE ViscosityModel = std::tuple_element_t< 2, TEST_TYPE >::value;
  using CompositionalMultiphaseFluid = typename FluidModel< FlashModel, ViscosityModel >::type;
  using Base = FluidModelTest< CompositionalMultiphaseFluid, std::tuple_element_t< 3, TEST_TYPE >::value >;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-3;

  static constexpr integer pointCount = 10;

  static constexpr real64 minPressure = 0.99;
  static constexpr real64 maxPressure = 600.1;

  static constexpr real64 minTemperature = 9.9;
  static constexpr real64 maxTemperature = 151.0;

public:
  MultiThermalFluidCompositionalMultiphaseTestFixture()
  {
    Base::createFluid( getFluidName(), [&]( CompositionalMultiphaseFluid & fluid ){
      fillPhysicalProperties( fluid );
      if constexpr (ViscosityModel == VISCOSITY_TYPE::PHILLIPS)
      {
        generateInterpolationPoints( fluid );
      }
      if constexpr (FlashModel == FLASH_TYPE::KVALUE)
      {
        generateTables( fluid );
      }
    } );
  }

  ~MultiThermalFluidCompositionalMultiphaseTestFixture() override
  {
    for( string const & fileName : m_fileNames )
    {
      removeFile( fileName );
    }
  }

  void testNumericalDerivatives( const bool useMass )
  {
    CompositionalMultiphaseFluid * fluid = this->getFluid( this->getFluidName() );

    fluid->setMassFlag( useMass );

    array2d< real64 > samples;
    FluidData< CompositionalMultiphaseFluid, Base::numComp >::getSamples( samples );
    integer const sampleCount = samples.size( 0 );
    Feed< Base::numComp > sample;

    real64 constexpr eps = 1.0e-6;

    constexpr real64 pressures[] = { 10.0e5, 50.0e5, 100.0e5, 600.0e5 };
    constexpr real64 temperatures[] = { 15.5, 24.0, 40.0, 80.0 };

    for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
    {
      for( integer ic = 0; ic < Base::numComp; ++ic )
      {
        sample[ic] = samples( sampleIndex, ic );
      }
      for( real64 const pressure : pressures )
      {
        for( real64 const temperature : temperatures )
        {
          typename Base::TestPoint const data ( pressure, units::convertCToK( temperature ), sample );
          Base::testNumericalDerivatives( fluid, data, eps, relTol, absTol );
        }
      }
    }
  }

  static string getFluidName();

private:
  static void fillPhysicalProperties( CompositionalMultiphaseFluid & fluid );
  void generateInterpolationPoints( CompositionalMultiphaseFluid & fluid );
  void generateTables( CompositionalMultiphaseFluid & fluid );

  static void writeToFile( string const & fileName, string const & content );

  static void removeFile( string const & fileName );
private:
  string_array m_fileNames;
};

template< typename TEST_TYPE >
string MultiThermalFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::getFluidName()
{
  return GEOS_FMT( "fluid_{}_{}_{}_{}",
                   Base::numComp,
                   EnumStrings< EquationOfStateType >::toString( EquationOfState ),
                   EnumStrings< FLASH_TYPE >::toString( FlashModel ),
                   EnumStrings< VISCOSITY_TYPE >::toString( ViscosityModel ) );
}

template< integer NUM_COMP >
static void fillBinaryCoeffs( array2d< real64 > & binaryCoeff, std::array< real64 const, NUM_COMP *(NUM_COMP-1)/2 > const data )
{
  auto bic = data.begin();
  binaryCoeff.resize( NUM_COMP, NUM_COMP );
  for( integer i = 0; i < NUM_COMP; ++i )
  {
    binaryCoeff( i, i ) = 0.0;
    for( integer j = i+1; j < NUM_COMP; ++j )
    {
      binaryCoeff( i, j ) = *bic++;
      binaryCoeff( j, i ) = binaryCoeff( i, j );
    }
  }
}

template< integer NUM_COMP >
static void populateArray( arraySlice1d< real64 > array, std::array< real64 const, NUM_COMP > const data )
{
  for( integer i = 0; i < NUM_COMP; ++i )
  {
    array[i] = data[i];
  }
}

template< typename FluidModel >
struct FluidData< FluidModel, 4 >
{
  static constexpr integer numComps = 4;
  static constexpr integer numCoeffs = 5;

  static void fillProperties( dataRepository::Group & fluid )
  {
    using Keys = typename FluidModel::viewKeyStruct;

    string_array & componentNames = fluid.getReference< string_array >( Keys::componentNamesString() );
    componentNames = {"N2", "C5", "C10", "H2O"};

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( Keys::componentMolarWeightString() );
    TestFluid< numComps >::createArray( molarWeight, Feed< numComps >{28e-3, 134e-3, 142e-3, 18e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( Keys::componentCriticalPressureString() );
    TestFluid< numComps >::createArray( criticalPressure, Feed< numComps >{34e5, 33.68e5, 21.03e5, 220.5e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( Keys::componentCriticalTemperatureString() );
    TestFluid< numComps >::createArray( criticalTemperature, Feed< numComps >{126.2, 469.7, 617.7, 647.0} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( Keys::componentAcentricFactorString() );
    TestFluid< numComps >::createArray( acentricFactor, Feed< numComps >{0.04, 0.2510, 0.4884, 0.344} );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( Keys::componentBinaryCoeffString() );
    fillBinaryCoeffs< numComps >( binaryCoeff, {0.0, 0.1, 0.0, 0.0, 0.0, 0.0} );

    array2d< real64 > & referenceEnthalpy = fluid.getReference< array2d< real64 > >( HeatCapacityCoefficients::viewKeyStruct::referenceEnthalpyString() );
    referenceEnthalpy.resize( 2, numComps );
    referenceEnthalpy.zero();

    array2d< real64 > & coefficients = fluid.getReference< array2d< real64 > >( HeatCapacityCoefficients::viewKeyStruct::componentHeatCapacityCoefficientsString() );
    coefficients.resize( numComps, numCoeffs );
    std::array< real64, numComps * numCoeffs > coefficientsData {
      2.89614125e+01, 2.03166820e-03, 2.47038198e-05, -9.54674193e-08, 1.15955044e-10, /* N2 */
      1.21004526e+02, 2.47982404e-01, 7.27115529e-04, -2.64999262e-06, 2.52476069e-09, /* C5 */
      3.11952000e+02, 5.26890000e-01, 8.16470000e-04, -4.03480000e-06, 9.60000000e-09, /* C10 */
      2.69032982e+01, 5.51993739e-02, 6.31405927e-04, -3.66166336e-06, 4.56291585e-09, /* H2O */
    };
    for( int ic = 0; ic < numComps; ++ic )
    {
      for( int j = 0; j < numCoeffs; ++j )
      {
        coefficients( ic, j ) = coefficientsData[ic*numCoeffs + j];
      }
    }
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 3, numComps );
    populateArray< numComps >( samples[0], {0.099, 0.300, 0.600, 0.001} );
    populateArray< numComps >( samples[1], {0.100, 0.250, 0.350, 0.300} );
    populateArray< numComps >( samples[2], {0.050, 0.050, 0.100, 0.800} );
  }
};

template< typename FluidModel >
struct FluidData< FluidModel, 5 >
{
  static constexpr integer numComps = 5;
  static constexpr integer numCoeffs = 5;

  static void fillProperties( dataRepository::Group & fluid )
  {
    using Keys = typename FluidModel::viewKeyStruct;

    string_array & componentNames = fluid.getReference< string_array >( Keys::componentNamesString() );
    componentNames = {"CO2", "N2", "C1", "H2O", "C4"};

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( Keys::componentMolarWeightString() );
    TestFluid< 5 >::createArray( molarWeight, Feed< 5 >{44.0098e-3, 28.0135e-3, 16.0428e-3, 18e-3, 82.4191e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( Keys::componentCriticalPressureString() );
    TestFluid< 5 >::createArray( criticalPressure, Feed< 5 >{73.77300e5, 33.95800e5, 45.99200e5, 48.71800e5, 33.20710e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( Keys::componentCriticalTemperatureString() );
    TestFluid< 5 >::createArray( criticalTemperature, Feed< 5 >{304.1280, 126.1920, 190.5640, 305.3300, 504.2160} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( Keys::componentAcentricFactorString() );
    TestFluid< 5 >::createArray( acentricFactor, Feed< 5 >{0.223000, 0.037200, 0.010400, 0.099100, 0.250274} );
    array1d< real64 > & volumeShift = fluid.getReference< array1d< real64 > >( Keys::componentVolumeShiftString() );
    TestFluid< 5 >::createArray( volumeShift, Feed< 5 >{1.845465e-01, -1.283880e-01, 9.225800e-02, 6.458060e-02, 0.000000e+00} );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( Keys::componentBinaryCoeffString() );
    fillBinaryCoeffs< 5 >( binaryCoeff, {0.0, 0.1, 0.03, 0.139, 0.032, 0.0, 0.12, 0.03, 0.0, 0.0} );

    array2d< real64 > & referenceEnthalpy = fluid.getReference< array2d< real64 > >( HeatCapacityCoefficients::viewKeyStruct::referenceEnthalpyString() );
    referenceEnthalpy.resize( 2, numComps );
    referenceEnthalpy.zero();

    array2d< real64 > & coefficients = fluid.getReference< array2d< real64 > >( HeatCapacityCoefficients::viewKeyStruct::componentHeatCapacityCoefficientsString() );
    coefficients.resize( numComps, numCoeffs );
    std::array< real64, numComps * numCoeffs > coefficientsData {
      3.40933377e+01, 7.52742132e-02, 2.00271751e-04, -1.61900936e-06, 2.12199565e-09, /* CO2 */
      2.89614125e+01, 2.03166820e-03, 2.47038198e-05, -9.54674193e-08, 1.15955044e-10, /* N2 */
      3.54137115e+01, 4.34392612e-02, 1.42852598e-04, -4.70393972e-07, 4.80614694e-10, /* CH4 */
      2.69032982e+01, 5.51993739e-02, 6.31405927e-04, -3.66166336e-06, 4.56291585e-09, /* H2O */
      5.28630051e+01, 1.16281994e-01, 1.40990713e-04, -6.38108433e-07, 7.14166667e-10, /* C4H10 */
    };
    for( int ic = 0; ic < numComps; ++ic )
    {
      for( int j = 0; j < numCoeffs; ++j )
      {
        coefficients( ic, j ) = coefficientsData[ic*numCoeffs + j];
      }
    }
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 1, 5 );
    populateArray< 5 >( samples[0], {0.050, 0.150, 0.550, 0.150, 0.100} );
  }
};

template< typename TEST_TYPE >
void MultiThermalFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::fillPhysicalProperties( CompositionalMultiphaseFluid & fluid )
{
  string_array & phaseNames = fluid.template getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames = {"oil", "gas"};

  string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfState );
  string_array & equationOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  equationOfState = {eosName, eosName};
  if constexpr (EquationOfState == EquationOfStateType::SoreideWhitson)
  {
    phaseNames[0] = "water";
    equationOfState[1] = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
  }

  dataRepository::Group & group = fluid;
  FluidData< CompositionalMultiphaseFluid, Base::numComp >::fillProperties( group );
}

template< typename TEST_TYPE >
void MultiThermalFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::generateInterpolationPoints( CompositionalMultiphaseFluid & fluid )
{
  integer constexpr N = pointCount;

  array1d< real64 > & pressureCoordinates = fluid.template getReference< array1d< real64 > >( PressureTemperatureCoordinates::viewKeyStruct::pressureCoordinatesString() );
  array1d< real64 > & temperatureCoordinates = fluid.template getReference< array1d< real64 > >( PressureTemperatureCoordinates::viewKeyStruct::temperatureCoordinatesString() );
  pressureCoordinates.resize( N+1 );
  temperatureCoordinates.resize( N+1 );

  real64 const dp = pow( maxPressure/minPressure, 1.0/N );
  real64 constexpr BAR_2_PA = 1.0e5;
  pressureCoordinates[0] = BAR_2_PA * minPressure;
  for( integer i = 1; i <= N; ++i )
  {
    pressureCoordinates[i] = pressureCoordinates[i-1] * dp;
  }

  real64 const dT = (maxTemperature-minTemperature)/N;
  for( integer i = 0; i <= N; ++i )
  {
    temperatureCoordinates[i] = minTemperature + i*dT + 273.15;
  }
}

template< typename TEST_TYPE >
void MultiThermalFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::generateTables( CompositionalMultiphaseFluid & fluid )
{
  string_array & kValueTables = fluid.template getReference< string_array >( KValueFlashParameters< 2 >::viewKeyStruct::kValueTablesString() );

  string_array const & compNames = fluid.componentNames();
  string const fluidName = fluid.getName();

  // Crookston correlations pressure (bar), temperature (K)
  mapBase< string, std::array< real64 const, 5 >, std::false_type > crookstonCoefficients = {
    {"H2O", {3.0620e+00, 8.9414e+02, 1.1912e-02, 5.3659e+02, 1.1951e+02}},
    {"CO2", {-1.8141e+00, 6.2655e+02, 6.7489e-03, 5.0732e+00, 2.5249e+02}},
    {"N2", {-3.5742e-01, 4.8660e+02, 5.7887e-03, 9.1910e+01, 1.8027e+02}},
    {"C1", {3.6577e+00, 7.1889e+02, 1.4154e-02, 5.7323e+02, 1.1088e+02}},
    {"C2", {8.2845e+00, 1.0144e+03, 5.2968e-02, 9.5800e+02, 9.3548e+01}},
    {"C4", {2.5068e+01, 2.2625e+03, 5.3888e-01, 1.7981e+03, 7.1948e+01}},
    {"C5", {3.6243e+01, 4.1985e+03, 1.9055e+00, 2.3575e+03, 5.8025e+01}},
    {"C10", {-2.3346e-01, 7.9356e-01, 6.8406e-03, 7.1217e+02, 1.3827e+02}}
  };
  std::ostringstream content;

  // Create arrays of pressure and temperature
  integer constexpr N = pointCount;
  array1d< real64 > temperature( N+1 );
  array1d< real64 > pressure( N+1 );

  real64 const dp = pow( maxPressure/minPressure, 1.0/N );
  content.str( "" );
  real64 constexpr BAR_2_PA = 1.0e5;
  pressure[0] = minPressure;
  content << BAR_2_PA * pressure[0] << "\n";
  for( integer i = 1; i <= N; ++i )
  {
    pressure[i] = pressure[i-1] * dp;
    content << BAR_2_PA * pressure[i] << "\n";
  }
  string const pressureFileName = GEOS_FMT( "{}_PRESSURE.txt", fluidName );
  writeToFile( pressureFileName, content.str());
  m_fileNames.emplace_back( pressureFileName );

  real64 const dT = (maxTemperature-minTemperature)/N;
  content.str( "" );
  for( integer i = 0; i <= N; ++i )
  {
    temperature[i] = minTemperature + i*dT + 273.15;
    content << temperature[i] << "\n";
  }
  string const temperatureFileName = GEOS_FMT( "{}_TEMP.txt", fluidName );
  writeToFile( temperatureFileName, content.str());
  m_fileNames.emplace_back( temperatureFileName );

  FunctionManager & functionManager = FunctionManager::getInstance();

  kValueTables.resize( Base::numComp );
  for( integer ic = 0; ic < Base::numComp; ++ic )
  {
    auto const it = crookstonCoefficients.find( compNames[ic] );
    ASSERT_FALSE( it == crookstonCoefficients.end() );
    auto const [A, B, C, D, E] = std::apply( []( auto &&... args ) { return std::make_tuple( args ... ); }, it->second );
    string const tableName = GEOS_FMT( "{}_KVALUE_{}", fluidName, ic+1 );
    content.str( "" );
    for( integer pi = 0; pi <= N; ++pi )
    {
      for( integer ti = 0; ti <= N; ++ti )
      {
        real64 const kValue = (A + B/pressure[pi] + C*pressure[pi])*LvArray::math::exp( -D/(temperature[ti]-E));
        content << LvArray::math::max( kValue, 1.0e-8 ) << "\n";
      }
    }
    string const kValueFileName = GEOS_FMT( "{}.txt", tableName );
    writeToFile( kValueFileName, content.str());
    m_fileNames.emplace_back( kValueFileName );

    TableFunction * tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );

    string_array & inputVarNames = tableFunction->getWrapper< string_array >( dataRepository::keys::inputVarNames ).reference();
    inputVarNames.emplace_back( "pressure" );
    inputVarNames.emplace_back( "temperature" );

    path_array & coordinateFiles = tableFunction->getWrapper< path_array >( TableFunction::viewKeyStruct::coordinateFilesString())
                                     .reference();
    coordinateFiles.emplace_back( Path() );
    coordinateFiles.emplace_back( Path() );
    dynamicCast< string & >( coordinateFiles[0] ) = pressureFileName;
    dynamicCast< string & >( coordinateFiles[1] ) = temperatureFileName;
    Path & voxelFile = tableFunction->getWrapper< Path >( TableFunction::viewKeyStruct::voxelFileString())
                         .reference();
    dynamicCast< string & >( voxelFile ) = kValueFileName;
    tableFunction->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    tableFunction->initializeFunction();
    kValueTables[ic] = tableName;
  }
}

template< typename TEST_TYPE >
void MultiThermalFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::writeToFile( string const & fileName, string const & content )
{
  std::ofstream os( fileName );
  ASSERT_TRUE( os.is_open() );
  os << content;
  os.close();
}

template< typename TEST_TYPE >
void MultiThermalFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::removeFile( string const & fileName )
{
  int const ret = std::remove( fileName.c_str() );
  ASSERT_TRUE( ret == 0 );
}

template< EquationOfStateType EOS, FLASH_TYPE FLASH, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
using TestType = std::tuple<
  std::integral_constant< EquationOfStateType, EOS >,
  std::integral_constant< FLASH_TYPE, FLASH >,
  std::integral_constant< VISCOSITY_TYPE, VISCOSITY >,
  std::integral_constant< integer, NUM_COMP >
  >;

using TestTypes = ::testing::Types<
  TestType< EquationOfStateType::PengRobinson, FLASH_TYPE::NEGATIVE, VISCOSITY_TYPE::LBC, 4 >,
  TestType< EquationOfStateType::SoreideWhitson, FLASH_TYPE::NEGATIVE, VISCOSITY_TYPE::PHILLIPS, 4 >,
  TestType< EquationOfStateType::PengRobinson, FLASH_TYPE::KVALUE, VISCOSITY_TYPE::PHILLIPS, 4 >,
  TestType< EquationOfStateType::SoaveRedlichKwong, FLASH_TYPE::NEGATIVE, VISCOSITY_TYPE::LBC, 5 >,
  TestType< EquationOfStateType::PengRobinson, FLASH_TYPE::NEGATIVE, VISCOSITY_TYPE::PHILLIPS, 5 >,
  TestType< EquationOfStateType::PengRobinson, FLASH_TYPE::KVALUE, VISCOSITY_TYPE::PHILLIPS, 5 >
  >;

class NameGenerator
{
public:
  template< typename T >
  static std::string GetName( int index )
  {
    return GEOS_FMT( "CompositionalFluid_{}_{}", MultiThermalFluidCompositionalMultiphaseTestFixture< T >::getFluidName(), index );
  }
};

TYPED_TEST_SUITE( MultiThermalFluidCompositionalMultiphaseTestFixture, TestTypes, NameGenerator );

TYPED_TEST( MultiThermalFluidCompositionalMultiphaseTestFixture, numericalDerivativesMolar )
{
  this->testNumericalDerivatives( false );
}

TYPED_TEST( MultiThermalFluidCompositionalMultiphaseTestFixture, numericalDerivativesMass )
{
  this->testNumericalDerivatives( true );
}

} // testing
} // geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::setupEnvironment( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::cleanupEnvironment();

  return result;
}
