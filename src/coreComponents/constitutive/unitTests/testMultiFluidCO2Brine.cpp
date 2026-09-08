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
 * @file testMultiFluidCO2Brine.cpp
 */

 #include "FluidModelTest.hpp"
 #include "constitutive/fluid/multifluid/CO2Brine/CO2BrineFluid.hpp"
 #include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::PVTProps;

namespace geos
{
namespace testing
{

enum class BrineModelType : int {Phillips, Ezrokhi};
enum class FlashType : int {DuanSun, SpycherPruess};

ENUM_STRINGS( BrineModelType, "Phillips", "Ezrokhi" );
ENUM_STRINGS( FlashType, "DuanSun", "SpycherPruess" );

static void populateFluidParameters( BrineFluidParameters & parameters )
{
  parameters.m_pressureCoordinates.emplace_back( 1.0e6 );
  parameters.m_pressureCoordinates.emplace_back( 1.5e7 );
  parameters.m_pressureInterval = 5.0e4;
  parameters.m_temperatureCoordinates.emplace_back( 367.15 );
  parameters.m_temperatureCoordinates.emplace_back( 369.15 );
  parameters.m_temperatureInterval = 1.0;
  parameters.m_salinity = 0.2;
  parameters.m_ezrokhiDensityCoefficients.emplace_back( 2.01e-6 );
  parameters.m_ezrokhiDensityCoefficients.emplace_back( -6.34e-7 );
  parameters.m_ezrokhiDensityCoefficients.emplace_back( 1e-4 );
  parameters.m_ezrokhiViscosityCoefficients.emplace_back( 2.42e-7 );
  parameters.m_ezrokhiViscosityCoefficients.emplace_back( 0.0 );
  parameters.m_ezrokhiViscosityCoefficients.emplace_back( 1e-4 );
}

template< BrineModelType BRINE, bool THERMAL >
struct FluidType {};

template<>
struct FluidType< BrineModelType::Phillips, false >
{
  using type = CO2BrinePhillipsFluid;
};
template<>
struct FluidType< BrineModelType::Phillips, true >
{
  using type = CO2BrinePhillipsThermalFluid;
};
template<>
struct FluidType< BrineModelType::Ezrokhi, false >
{
  using type = CO2BrineEzrokhiFluid;
};
template<>
struct FluidType< BrineModelType::Ezrokhi, true >
{
  using type = CO2BrineEzrokhiThermalFluid;
};

template< FlashType FLASH >
struct FlashModel {};

template<>
struct FlashModel< FlashType::DuanSun >
{
  static void populate( BrineFluidParameters & GEOS_UNUSED_PARAM( parameters ) )
  {}
};

template<>
struct FlashModel< FlashType::SpycherPruess >
{
  static void populate( BrineFluidParameters & parameters )
  {
    parameters.m_solubilityModel = BrineFluidParameters::SolubilityModel::SpycherPruess;
    parameters.m_tolerance = 1.0e-10;
  }
};

template< BrineModelType BRINE, FlashType FLASH, bool THERMAL >
class MultiFluidCO2BrineTestFixture : public FluidModelTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >,
  public ::testing::WithParamInterface<
    std::tuple< typename FluidModelTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >::TestPoint,
                typename FluidModelTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >::TestResult > >
{
public:
  using CO2BrineFluid = typename FluidType< BRINE, THERMAL >::type;
  using Base = FluidModelTest< CO2BrineFluid, 2, 2 >;

public:
  MultiFluidCO2BrineTestFixture()
  {
    Base::createFluid( getFluidName(), []( CO2BrineFluid & fluid ){
      fillPhysicalProperties( fluid );
    } );
  }

  ~MultiFluidCO2BrineTestFixture() override = default;

  // Test numerical derivatives at selected data points
  void testNumericalDerivatives( bool const useMass )
  {
    CO2BrineFluid * fluid = this->getFluid( this->getFluidName() );

    fluid->setMassFlag( useMass );

    real64 constexpr eps = 1.0e-6;

    // Some of the functions are simply table lookups. We need to keep the test points away from
    // the table nodes because the kink in the linear interpolation might cause numerical derivative
    // mismatches. Some of these values have been manually inspected and the differences, although
    // not meeting the tolerance here, are small as expected.
    constexpr real64 temperatures[] = { 367.65, 368.00, 368.75 };
    constexpr real64 pressures[] = { 20.01e5, 75.01e5, 120.1e5 };
    auto const samples = { Feed< 2 >{0.7, 0.3}, Feed< 2 >{0.01, 0.99}, Feed< 2 >{0.99, 0.01} };

    for( auto const & sample : samples )
    {
      for( real64 const pressure : pressures )
      {
        for( real64 const temperature : temperatures )
        {
          typename Base::TestPoint const data ( pressure, temperature, sample );
          Base::testNumericalDerivatives( fluid, data, eps );
        }
      }
    }
  }

  static string getFluidName();

private:
  static void fillPhysicalProperties( CO2BrineFluid & fluid );
};

template< BrineModelType BRINE, FlashType FLASH, bool THERMAL >
string MultiFluidCO2BrineTestFixture< BRINE, FLASH, THERMAL >::getFluidName()
{
  return GEOS_FMT( "fluid{}{}{}",
                   EnumStrings< BrineModelType >::toString( BRINE ),
                   EnumStrings< FlashType >::toString( FLASH ),
                   (THERMAL ? "Thermal" : ""));
}

template< BrineModelType BRINE, FlashType FLASH, bool THERMAL >
void MultiFluidCO2BrineTestFixture< BRINE, FLASH, THERMAL >::fillPhysicalProperties( CO2BrineFluid & fluid )
{
  dataRepository::Group & group = fluid;

  auto & phaseNames = group.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames = {"gas", "liquid"};

  auto & compNames = group.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  compNames = {"co2", "water"};

  auto & logLevel = group.getReference< integer >( dataRepository::Group::viewKeyStruct::logLevelString() );
  logLevel = 0;

  auto & molarWeight = group.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  Base::fill( molarWeight, Feed< 2 >{44e-3, 18e-3} );

  BrineFluidParameters parameters;
  populateFluidParameters( parameters );
  FlashModel< FLASH >::populate( parameters );

  using Keys = BrineFluidParameters::viewKeyStruct;
  auto setValue = [&]( string const key, auto const value )
  {
    using T = typename std::remove_const< decltype(value) >::type;
    group.getReference< T >( key ) = value;
  };
  setValue( Keys::solubilityModelString(), parameters.m_solubilityModel );
  setValue( Keys::pressureCoordinatesString(), parameters.m_pressureCoordinates );
  setValue( Keys::pressureIntervalString(), parameters.m_pressureInterval );
  setValue( Keys::temperatureCoordinatesString(), parameters.m_temperatureCoordinates );
  setValue( Keys::temperatureIntervalString(), parameters.m_temperatureInterval );
  setValue( Keys::salinityString(), parameters.m_salinity );
  setValue( Keys::toleranceString(), parameters.m_tolerance );
  setValue( Keys::waterCompressibilityString(), parameters.m_waterCompressibility );
  setValue( Keys::solubilityTablesString(), parameters.m_solubilityTables );
  if constexpr (BRINE == BrineModelType::Ezrokhi)
  {
    setValue( Keys::ezrokhiDensityCoefficientsString(), parameters.m_ezrokhiDensityCoefficients );
    setValue( Keys::ezrokhiViscosityCoefficientsString(), parameters.m_ezrokhiViscosityCoefficients );
  }
}

using CO2BrinePhillipsTest = MultiFluidCO2BrineTestFixture< BrineModelType::Phillips,
                                                            FlashType::DuanSun,
                                                            false >;
using CO2BrineEzrokhiTest = MultiFluidCO2BrineTestFixture< BrineModelType::Ezrokhi,
                                                           FlashType::DuanSun,
                                                           false >;
using CO2BrinePhillipsThermalTest = MultiFluidCO2BrineTestFixture< BrineModelType::Phillips,
                                                                   FlashType::DuanSun,
                                                                   true >;
#if !defined(GEOS_DEVICE_COMPILE)
using CO2BrineEzrokhiThermalTest = MultiFluidCO2BrineTestFixture< BrineModelType::Ezrokhi,
                                                                  FlashType::DuanSun,
                                                                  true >;
#endif
using CO2BrinePhillipsSpycherPruessTest = MultiFluidCO2BrineTestFixture< BrineModelType::Phillips,
                                                                         FlashType::SpycherPruess,
                                                                         false >;

TEST_F( CO2BrinePhillipsTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( CO2BrinePhillipsTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}
TEST_F( CO2BrineEzrokhiTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( CO2BrineEzrokhiTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}
TEST_F( CO2BrinePhillipsThermalTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( CO2BrinePhillipsThermalTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}
#if !defined(GEOS_DEVICE_COMPILE)
TEST_F( CO2BrineEzrokhiThermalTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( CO2BrineEzrokhiThermalTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}
#endif
TEST_F( CO2BrinePhillipsSpycherPruessTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( CO2BrinePhillipsSpycherPruessTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

TEST_P( CO2BrinePhillipsTest, testFluidValues )
{
  auto const [testPoint, testResult] = GetParam();
  CO2BrineFluid * fluid = this->getFluid( this->getFluidName() );
  testValuesAgainstPreviousImplementation( fluid, testPoint, testResult );
}

TEST_P( CO2BrineEzrokhiTest, testFluidValues )
{
  auto const [testPoint, testResult] = GetParam();
  CO2BrineFluid * fluid = this->getFluid( this->getFluidName() );
  testValuesAgainstPreviousImplementation( fluid, testPoint, testResult );
}

TEST_P( CO2BrinePhillipsThermalTest, testFluidValues )
{
  auto const [testPoint, testResult] = GetParam();
  CO2BrineFluid * fluid = this->getFluid( this->getFluidName() );
  testValuesAgainstPreviousImplementation( fluid, testPoint, testResult );
}

#if !defined(GEOS_DEVICE_COMPILE)
TEST_P( CO2BrineEzrokhiThermalTest, testFluidValues )
{
  auto const [testPoint, testResult] = GetParam();
  CO2BrineFluid * fluid = this->getFluid( this->getFluidName() );
  testValuesAgainstPreviousImplementation( fluid, testPoint, testResult );
}
#endif

TEST_P( CO2BrinePhillipsSpycherPruessTest, testFluidValues )
{
  auto const [testPoint, testResult] = GetParam();
  CO2BrineFluid * fluid = this->getFluid( this->getFluidName() );
  testValuesAgainstPreviousImplementation( fluid, testPoint, testResult );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrinePhillipsTest,
  ::testing::ValuesIn<CO2BrinePhillipsTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity           | phase enthalpy            | phase internal energy     | density    
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.94192e-01, 7.05808e-01}, {1.87668e+03, 5.33005e+04}, {8.25738e+01, 9.70812e+02}, {1.90427e-05, 3.06891e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.88220e+03}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.94206e-01, 7.05794e-01}, {1.87383e+03, 5.32850e+04}, {8.24487e+01, 9.70502e+02}, {1.90569e-05, 3.05730e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.87360e+03}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.94236e-01, 7.05764e-01}, {1.86780e+03, 5.32516e+04}, {8.21833e+01, 9.69835e+02}, {1.90872e-05, 3.03243e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.85534e+03}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.92018e-01, 7.07982e-01}, {3.05367e+03, 5.32538e+04}, {1.34361e+02, 9.74179e+02}, {2.00622e-05, 3.06891e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.18078e+03}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.92035e-01, 7.07965e-01}, {3.04766e+03, 5.32386e+04}, {1.34097e+02, 9.73867e+02}, {2.00722e-05, 3.05730e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.16417e+03}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.92070e-01, 7.07930e-01}, {3.03496e+03, 5.32057e+04}, {1.33538e+02, 9.73198e+02}, {2.00938e-05, 3.03243e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.12901e+03}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{2.89169e-01, 7.10831e-01}, {5.77607e+03, 5.32402e+04}, {2.54147e+02, 9.79416e+02}, {2.39022e-05, 3.06891e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.57692e+04}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{2.89185e-01, 7.10815e-01}, {5.75768e+03, 5.32251e+04}, {2.53338e+02, 9.79107e+02}, {2.38854e-05, 3.05730e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.57280e+04}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{2.89218e-01, 7.10782e-01}, {5.71917e+03, 5.31925e+04}, {2.51643e+02, 9.78443e+02}, {2.38514e-05, 3.03243e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.56415e+04}},
  })
);

INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrineEzrokhiTest,
  ::testing::ValuesIn<CO2BrineEzrokhiTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy     | density    
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.94192e-01, 7.05808e-01}, {1.87668e+03, 5.51522e+04}, {8.25738e+01, 1.00454e+03}, {1.90427e-05, 3.12027e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.89763e+03}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.94206e-01, 7.05794e-01}, {1.87383e+03, 5.51511e+04}, {8.24487e+01, 1.00449e+03}, {1.90569e-05, 3.10902e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.88910e+03}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.94236e-01, 7.05764e-01}, {1.86780e+03, 5.51487e+04}, {8.21833e+01, 1.00439e+03}, {1.90872e-05, 3.08490e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.87101e+03}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.92018e-01, 7.07982e-01}, {3.05367e+03, 5.57997e+04}, {1.34361e+02, 1.02075e+03}, {2.00622e-05, 3.16706e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.23220e+03}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.92035e-01, 7.07965e-01}, {3.04766e+03, 5.58041e+04}, {1.34097e+02, 1.02080e+03}, {2.00722e-05, 3.15594e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.21580e+03}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.92070e-01, 7.07930e-01}, {3.03496e+03, 5.58136e+04}, {1.33538e+02, 1.02090e+03}, {2.00938e-05, 3.13210e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.18112e+03}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{2.89169e-01, 7.10831e-01}, {5.77607e+03, 5.66765e+04}, {2.54147e+02, 1.04263e+03}, {2.39022e-05, 3.22839e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.59731e+04}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{2.89185e-01, 7.10815e-01}, {5.75768e+03, 5.66893e+04}, {2.53338e+02, 1.04283e+03}, {2.38854e-05, 3.21754e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.59325e+04}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{2.89218e-01, 7.10782e-01}, {5.71917e+03, 5.67171e+04}, {2.51643e+02, 1.04328e+03}, {2.38514e-05, 3.19427e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.58474e+04}},
  })
);

INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrinePhillipsThermalTest,
  ::testing::ValuesIn<CO2BrinePhillipsThermalTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy     | density    
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.94192e-01, 7.05808e-01}, {1.87668e+03, 5.33005e+04}, {8.25738e+01, 9.70812e+02}, {1.90427e-05, 3.06891e-04}, {1.21447e+07, 2.15290e+07}, {1.20841e+07, 2.15238e+07}, 5.88220e+03}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.94206e-01, 7.05794e-01}, {1.87383e+03, 5.32850e+04}, {8.24487e+01, 9.70502e+02}, {1.90569e-05, 3.05730e-04}, {1.21537e+07, 2.16115e+07}, {1.20931e+07, 2.16064e+07}, 5.87360e+03}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.94236e-01, 7.05764e-01}, {1.86780e+03, 5.32516e+04}, {8.21833e+01, 9.69835e+02}, {1.90872e-05, 3.03243e-04}, {1.21731e+07, 2.17883e+07}, {1.21123e+07, 2.17832e+07}, 5.85534e+03}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.92018e-01, 7.07982e-01}, {3.05367e+03, 5.32538e+04}, {1.34361e+02, 9.74179e+02}, {2.00622e-05, 3.06891e-04}, {1.17163e+07, 2.14953e+07}, {1.16605e+07, 2.14876e+07}, 9.18078e+03}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.92035e-01, 7.07965e-01}, {3.04766e+03, 5.32386e+04}, {1.34097e+02, 9.73867e+02}, {2.00722e-05, 3.05730e-04}, {1.17269e+07, 2.15777e+07}, {1.16709e+07, 2.15700e+07}, 9.16417e+03}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.92070e-01, 7.07930e-01}, {3.03496e+03, 5.32057e+04}, {1.33538e+02, 9.73198e+02}, {2.00938e-05, 3.03243e-04}, {1.17494e+07, 2.17542e+07}, {1.16932e+07, 2.17465e+07}, 9.12901e+03}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{2.89169e-01, 7.10831e-01}, {5.77607e+03, 5.32402e+04}, {2.54147e+02, 9.79416e+02}, {2.39022e-05, 3.06891e-04}, {1.08306e+07, 2.14426e+07}, {1.07833e+07, 2.14304e+07}, 1.57692e+04}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{2.89185e-01, 7.10815e-01}, {5.75768e+03, 5.32251e+04}, {2.53338e+02, 9.79107e+02}, {2.38854e-05, 3.05730e-04}, {1.08453e+07, 2.15248e+07}, {1.07979e+07, 2.15125e+07}, 1.57280e+04}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{2.89218e-01, 7.10782e-01}, {5.71917e+03, 5.31925e+04}, {2.51643e+02, 9.78443e+02}, {2.38514e-05, 3.03243e-04}, {1.08767e+07, 2.17008e+07}, {1.08290e+07, 2.16885e+07}, 1.56415e+04}},
  })
);

#if !defined(GEOS_DEVICE_COMPILE)
INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrineEzrokhiThermalTest,
  ::testing::ValuesIn<CO2BrineEzrokhiThermalTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    | density      |
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.94192e-01, 7.05808e-01}, {1.87668e+03, 5.51522e+04}, {8.25738e+01, 1.00454e+03}, {1.90427e-05, 3.12027e-04}, {1.21447e+07, 2.15290e+07}, {1.20841e+07, 2.15240e+07}, 5.89763e+03}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.94206e-01, 7.05794e-01}, {1.87383e+03, 5.51511e+04}, {8.24487e+01, 1.00449e+03}, {1.90569e-05, 3.10902e-04}, {1.21537e+07, 2.16115e+07}, {1.20931e+07, 2.16065e+07}, 5.88910e+03}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.94236e-01, 7.05764e-01}, {1.86780e+03, 5.51487e+04}, {8.21833e+01, 1.00439e+03}, {1.90872e-05, 3.08490e-04}, {1.21731e+07, 2.17883e+07}, {1.21123e+07, 2.17834e+07}, 5.87101e+03}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.92018e-01, 7.07982e-01}, {3.05367e+03, 5.57997e+04}, {1.34361e+02, 1.02075e+03}, {2.00622e-05, 3.16706e-04}, {1.17163e+07, 2.14953e+07}, {1.16605e+07, 2.14880e+07}, 9.23220e+03}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.92035e-01, 7.07965e-01}, {3.04766e+03, 5.58041e+04}, {1.34097e+02, 1.02080e+03}, {2.00722e-05, 3.15594e-04}, {1.17269e+07, 2.15777e+07}, {1.16709e+07, 2.15704e+07}, 9.21580e+03}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.92070e-01, 7.07930e-01}, {3.03496e+03, 5.58136e+04}, {1.33538e+02, 1.02090e+03}, {2.00938e-05, 3.13210e-04}, {1.17494e+07, 2.17542e+07}, {1.16932e+07, 2.17468e+07}, 9.18112e+03}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{2.89169e-01, 7.10831e-01}, {5.77607e+03, 5.66765e+04}, {2.54147e+02, 1.04263e+03}, {2.39022e-05, 3.22839e-04}, {1.08306e+07, 2.14426e+07}, {1.07833e+07, 2.14311e+07}, 1.59731e+04}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{2.89185e-01, 7.10815e-01}, {5.75768e+03, 5.66893e+04}, {2.53338e+02, 1.04283e+03}, {2.38854e-05, 3.21754e-04}, {1.08453e+07, 2.15248e+07}, {1.07979e+07, 2.15133e+07}, 1.59325e+04}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{2.89218e-01, 7.10782e-01}, {5.71917e+03, 5.67171e+04}, {2.51643e+02, 1.04328e+03}, {2.38514e-05, 3.19427e-04}, {1.08767e+07, 2.17008e+07}, {1.08290e+07, 2.16893e+07}, 1.58474e+04}},
  })
);
#endif

INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrinePhillipsSpycherPruessTest,
  ::testing::ValuesIn<CO2BrinePhillipsSpycherPruessTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    | density      |
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{3.00488e-01, 6.99512e-01}, {1.87668e+03, 5.32842e+04}, {8.25738e+01, 9.70985e+02}, {1.90427e-05, 3.06891e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.77217e+03}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{3.00580e-01, 6.99420e-01}, {1.87383e+03, 5.32686e+04}, {8.24487e+01, 9.70677e+02}, {1.90569e-05, 3.05730e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.76239e+03}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{3.00781e-01, 6.99219e-01}, {1.86780e+03, 5.32348e+04}, {8.21833e+01, 9.70013e+02}, {1.90872e-05, 3.03243e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 5.74154e+03}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{2.96731e-01, 7.03269e-01}, {3.05367e+03, 5.32295e+04}, {1.34361e+02, 9.74440e+02}, {2.00622e-05, 3.06891e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.05928e+03}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{2.96804e-01, 7.03196e-01}, {3.04766e+03, 5.32141e+04}, {1.34097e+02, 9.74130e+02}, {2.00722e-05, 3.05730e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.04143e+03}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{2.96962e-01, 7.03038e-01}, {3.03496e+03, 5.31808e+04}, {1.33538e+02, 9.73465e+02}, {2.00938e-05, 3.03243e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 9.00358e+03}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{2.93050e-01, 7.06950e-01}, {5.77607e+03, 5.32048e+04}, {2.54147e+02, 9.79798e+02}, {2.39022e-05, 3.06891e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.56195e+04}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{2.93105e-01, 7.06895e-01}, {5.75768e+03, 5.31894e+04}, {2.53338e+02, 9.79492e+02}, {2.38854e-05, 3.05730e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.55771e+04}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{2.93225e-01, 7.06775e-01}, {5.71917e+03, 5.31561e+04}, {2.51643e+02, 9.78834e+02}, {2.38514e-05, 3.03243e-04}, {0.00000e+00, 0.00000e+00}, {0.00000e+00, 0.00000e+00}, 1.54878e+04}},
  })
);

/* UNCRUSTIFY-ON */

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
