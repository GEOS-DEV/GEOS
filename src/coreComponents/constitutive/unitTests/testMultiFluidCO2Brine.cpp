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
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29425923996094233, 0.70574076003905761}, {1876.6771637012978, 53305.013964111968}, {82.573795202857099, 970.76392873772056}, {1.9042690398106135e-05, 0.00030689078075780307}, {0, 0}, {0, 0}, 5881.0495330397516}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29427338244869544, 0.7057266175513045}, {1873.8349354054078, 53289.497799926183}, {82.44873715783794, 970.45381664230172}, {1.9056851253889581e-05, 0.00030573032322826445}, {0, 0}, {0, 0}, 5872.451138710232}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29430307582076554, 0.70569692417923446}, {1867.8017150981846, 53256.094960688722}, {82.183275464320118, 969.78772926333068}, {1.908722569364096e-05, 0.00030324340751602518}, {0, 0}, {0, 0}, 5854.1990797320386}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29211045005945158, 0.70788954994054842}, {3053.6680173832319, 53259.968646826892}, {134.36139276486219, 974.11281327244239}, {2.0062244341948444e-05, 0.00030689078075780307}, {0, 0}, {0, 0}, 9178.5148241844036}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29212705231288177, 0.70787294768711817}, {3047.6603410545627, 53244.705654864556}, {134.09705500640075, 973.8015500126096}, {2.0072211093265288e-05, 0.00030573032322826445}, {0, 0}, {0, 0}, 9161.9054628817212}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29216186940719618, 0.70783813059280376}, {3034.9588123747385, 53211.834244050187}, {133.53818774448848, 973.13306413318696}, {2.0093772480867798e-05, 0.00030324340751602518}, {0, 0}, {0, 0}, 9126.7704334118753}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{0.28929387549866292, 0.71070612450133697}, {5776.0682830082569, 53248.450886578539}, {254.14700445236329, 979.32770964345616}, {2.3902169768739735e-05, 0.00030689078075780307}, {0, 0}, {0, 0}, 15764.939745650383}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{0.28930965444324924, 0.7106903455567507}, {5757.678026500017, 53233.305520700458}, {253.33783316600073, 979.0188946299935}, {2.3885372962202469e-05, 0.00030573032322826445}, {0, 0}, {0, 0}, 15723.737825049826}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{0.28934269592211376, 0.71065730407788619}, {5719.1666421929785, 53200.683841331775}, {251.64333225649102, 978.35560186955672}, {2.3851411358173899e-05, 0.00030324340751602518}, {0, 0}, {0, 0}, 15637.263889592703}},
  })
);

INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrineEzrokhiTest,
  ::testing::ValuesIn<CO2BrineEzrokhiTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy     | density    
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29425923996094233, 0.70574076003905761}, {1876.6771637012978, 55134.109523922023}, {82.573795202857099, 1004.0744910963268}, {1.9042690398106135e-05, 0.00031188253134058777}, {0, 0}, {0, 0}, 5896.2804531365746}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29427338244869544, 0.7057266175513045}, {1873.8349354054078, 55132.800406255315}, {82.44873715783794, 1004.022157936396}, {1.9056851253889581e-05, 0.00031075711115315639}, {0, 0}, {0, 0}, 5887.760285482098}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29430307582076554, 0.70569692417923446}, {1867.8017150981846, 55130.087877585065}, {82.183275464320118, 1003.9129375962731}, {1.908722569364096e-05, 0.0003083448869739105}, {0, 0}, {0, 0}, 5869.6769017971983}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29211045005945158, 0.70788954994054842}, {3053.6680173832319, 55774.915348820658}, {134.36139276486219, 1020.1106211824431}, {2.0062244341948444e-05, 0.00031650765699059208}, {0, 0}, {0, 0}, 9229.2833853882985}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29212705231288177, 0.70787294768711817}, {3047.6603410545627, 55779.080038991713}, {134.09705500640075, 1020.153157617934}, {2.0072211093265288e-05, 0.00031539513670803205}, {0, 0}, {0, 0}, 9212.8925404158108}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29216186940719618, 0.70783813059280376}, {3034.9588123747385, 55788.169647549992}, {133.53818774448848, 1020.2488458208462}, {2.0093772480867798e-05, 0.00031301031595972021}, {0, 0}, {0, 0}, 9178.2295176963744}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{0.28929387549866292, 0.71070612450133697}, {5776.0682830082569, 56643.071805783446}, {254.14700445236329, 1041.7604428884836}, {2.3902169768739735e-05, 0.00032257012724654506}, {0, 0}, {0, 0}, 15966.276717799727}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{0.28930965444324924, 0.7106903455567507}, {5757.678026500017, 56655.596406194534}, {253.33783316600073, 1041.9585788567383}, {2.3885372962202469e-05, 0.00032148402900293452}, {0, 0}, {0, 0}, 15925.678980727182}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{0.28934269592211376, 0.71065730407788619}, {5719.1666421929785, 56682.73097440914}, {251.64333225649102, 1042.3901230945153}, {2.3851411358173899e-05, 0.00031915537643445436}, {0, 0}, {0, 0}, 15840.526389642862}},
  })
);

INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrinePhillipsThermalTest,
  ::testing::ValuesIn<CO2BrinePhillipsThermalTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy     | density    
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29425923996094233, 0.70574076003905761}, {1876.6771637012978, 53305.013964111968}, {82.573795202857099, 970.76392873772056}, {1.9042690398106135e-05, 0.00030689078075780307}, {23512.098196924297, 7134.9209150422157}, {20847.281841330816, 7041.1023484929792}, 5881.0495330397516}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29427338244869544, 0.7057266175513045}, {1873.8349354054078, 53289.497799926183}, {82.44873715783794, 970.45381664230172}, {1.9056851253889581e-05, 0.00030573032322826445}, {23529.648484157646, 7161.4005915002635}, {20860.790141564077, 7067.5547080435254}, 5872.451138710232}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29430307582076554, 0.70569692417923446}, {1867.8017150981846, 53256.094960688722}, {82.183275464320118, 969.78772926333068}, {1.908722569364096e-05, 0.00030324340751602518}, {23567.202639605734, 7218.1284222828845}, {20889.723569116115, 7124.2236776017071}, 5854.1990797320386}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29211045005945158, 0.70788954994054842}, {3053.6680173832319, 53259.968646826892}, {134.36139276486219, 974.11281327244239}, {2.0062244341948444e-05, 0.00030689078075780307}, {22682.777979740778, 7175.3907169914946}, {20226.387842600787, 7034.5532326561206}, 9178.5148241844036}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29212705231288177, 0.70787294768711817}, {3047.6603410545627, 53244.705654864556}, {134.09705500640075, 973.8015500126096}, {2.0072211093265288e-05, 0.00030573032322826445}, {22703.191217798885, 7201.8386517373838}, {20241.958941040488, 7060.9607952867927}, 9161.9054628817212}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29216186940719618, 0.70783813059280376}, {3034.9588123747385, 53211.834244050187}, {133.53818774448848, 973.13306413318696}, {2.0093772480867798e-05, 0.00030324340751602518}, {22746.780966371905, 7258.4994748145491}, {20275.24824921758, 7117.5345915970775}, 9126.7704334118753}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{0.28929387549866292, 0.71070612450133697}, {5776.0682830082569, 53248.450886578539}, {254.14700445236329, 979.32770964345616}, {2.3902169768739735e-05, 0.00030689078075780307}, {20967.805250887584, 7211.0130368847103}, {18888.536202886473, 6985.4665693342868}, 15764.939745650383}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{0.28930965444324924, 0.7106903455567507}, {5757.678026500017, 53233.305520700458}, {253.33783316600073, 979.0188946299935}, {2.3885372962202469e-05, 0.00030573032322826445}, {20996.406587066806, 7237.6170471546557}, {18910.496269622283, 7012.0064095561738}, 15723.737825049826}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{0.28934269592211376, 0.71065730407788619}, {5719.1666421929785, 53200.683841331775}, {251.64333225649102, 978.35560186955672}, {2.3851411358173899e-05, 0.00030324340751602518}, {21057.120003813936, 7294.6025609247663}, {18957.163707489821, 7068.8535830390838}, 15637.263889592703}},
  })
);

#if !defined(GEOS_DEVICE_COMPILE)
INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrineEzrokhiThermalTest,
  ::testing::ValuesIn<CO2BrineEzrokhiThermalTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    | density      |
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29425923996094233, 0.70574076003905761}, {1876.6771637012978, 55134.109523922023}, {82.573795202857099, 1004.0744910963268}, {1.9042690398106135e-05, 0.00031188253134058777}, {23512.098196924297, 7134.9209150422157}, {20847.281841330816, 7044.2148159833359}, 5896.2804531365746}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29427338244869544, 0.7057266175513045}, {1873.8349354054078, 55132.800406255315}, {82.44873715783794, 1004.022157936396}, {1.9056851253889581e-05, 0.00031075711115315639}, {23529.648484157646, 7161.4005915002635}, {20860.790141564077, 7070.6923386426306}, 5887.760285482098}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29430307582076554, 0.70569692417923446}, {1867.8017150981846, 55130.087877585065}, {82.183275464320118, 1003.9129375962731}, {1.908722569364096e-05, 0.0003083448869739105}, {23567.202639605734, 7218.1284222828845}, {20889.723569116115, 7127.4157063680359}, 5869.6769017971983}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29211045005945158, 0.70788954994054842}, {3053.6680173832319, 55774.915348820658}, {134.36139276486219, 1020.1106211824431}, {2.0062244341948444e-05, 0.00031650765699059208}, {22682.777979740778, 7175.3907169914946}, {20226.387842600787, 7040.9037356471363}, 9229.2833853882985}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29212705231288177, 0.70787294768711817}, {3047.6603410545627, 55779.080038991713}, {134.09705500640075, 1020.153157617934}, {2.0072211093265288e-05, 0.00031539513670803205}, {22703.191217798885, 7201.8386517373838}, {20241.958941040488, 7067.3617117312651}, 9212.8925404158108}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29216186940719618, 0.70783813059280376}, {3034.9588123747385, 55788.169647549992}, {133.53818774448848, 1020.2488458208462}, {2.0093772480867798e-05, 0.00031301031595972021}, {22746.780966371905, 7258.4994748145491}, {20275.24824921758, 7124.0444452376933}, 9178.2295176963744}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{0.28929387549866292, 0.71070612450133697}, {5776.0682830082569, 56643.071805783446}, {254.14700445236329, 1041.7604428884836}, {2.3902169768739735e-05, 0.00032257012724654506}, {20967.805250887584, 7211.0130368847103}, {18888.536202886473, 6998.9835756086759}, 15966.276717799727}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{0.28930965444324924, 0.7106903455567507}, {5757.678026500017, 56655.596406194534}, {253.33783316600073, 1041.9585788567383}, {2.3885372962202469e-05, 0.00032148402900293452}, {20996.406587066806, 7237.6170471546557}, {18910.496269622283, 7025.6344582874617}, 15925.678980727182}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{0.28934269592211376, 0.71065730407788619}, {5719.1666421929785, 56682.73097440914}, {251.64333225649102, 1042.3901230945153}, {2.3851411358173899e-05, 0.00031915537643445436}, {21057.120003813936, 7294.6025609247663}, {18957.163707489821, 7082.7214501606741}, 15840.526389642862}},
  })
);
#endif

INSTANTIATE_TEST_SUITE_P(
  FluidValueTest, CO2BrinePhillipsSpycherPruessTest,
  ::testing::ValuesIn<CO2BrinePhillipsSpycherPruessTest::ParamType>({
     //| pressure  | temp       | composition    | phase fraction             |  phase density            |  phase mass density       | phase viscosity          | phase enthalpy             | phase internal energy    | density      |
      {{5.00100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.30055901953547864, 0.6994409804645213}, {1876.6771637012978, 53288.926502975992}, {82.573795202857099, 970.93528779977464}, {1.9042690398106135e-05, 0.00030689078075780307}, {0, 0}, {0, 0}, 5770.995331923642}},
      {{5.00100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.30065099332750278, 0.69934900667249711}, {1873.8349354054078, 53273.265089499575}, {82.44873715783794, 970.62645912075811}, {1.9056851253889581e-05, 0.00030573032322826445}, {0, 0}, {0, 0}, 5761.2158486916578}},
      {{5.00100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.30085124145379533, 0.69914875854620462}, {1867.8017150981846, 53239.531571829604}, {82.183275464320118, 969.96330012860449}, {1.908722569364096e-05, 0.00030324340751602518}, {0, 0}, {0, 0}, 5740.3797818912908}},
      {{7.50100e+06, 3.67650e+02, {0.300, 0.700}}, {{0.29682857182081501, 0.70317142817918499}, {3053.6680175979145, 53235.944853242792}, {134.36139277430823, 974.37066442194998}, {2.0062244342169086e-05, 0.00030689078075780307}, {0, 0}, {0, 0}, 9056.9433658296184}},
      {{7.50100e+06, 3.68000e+02, {0.300, 0.700}}, {{0.29690136340438578, 0.70309863659561422}, {3047.6603414195229, 53220.487452567322}, {134.09705502245899, 974.06107893123135}, {2.0072211093640384e-05, 0.00030573032322826445}, {0, 0}, {0, 0}, 9039.0992701421255}},
      {{7.50100e+06, 3.68750e+02, {0.300, 0.700}}, {{0.29705903023548158, 0.70294096976451825}, {3034.9588125464843, 53187.17198747053}, {133.5381877520453, 973.39644199990892}, {2.0093772481044313e-05, 0.00030324340751602518}, {0, 0}, {0, 0}, 9001.2678623312495}},
      {{1.20100e+07, 3.67650e+02, {0.300, 0.700}}, {{0.29318209547310786, 0.70681790452689197}, {5776.0682830082569, 53213.417559934132}, {254.14700445236329, 979.7053949450401}, {2.3902169768739735e-05, 0.00030689078075780307}, {0, 0}, {0, 0}, 15615.051798198872}},
      {{1.20100e+07, 3.68000e+02, {0.300, 0.700}}, {{0.29323762048186658, 0.70676237951813337}, {5757.678026500017, 53197.977653174094}, {253.33783316600073, 979.39914170660757}, {2.3885372962202469e-05, 0.00030573032322826445}, {0, 0}, {0, 0}, 15572.600113462746}},
      {{1.20100e+07, 3.68750e+02, {0.300, 0.700}}, {{0.29335735073442271, 0.70664264926557729}, {5719.1666421929785, 53164.689255500016}, {251.64333225649102, 978.74165946166147}, {2.3851411358173899e-05, 0.00030324340751602518}, {0, 0}, {0, 0}, 15483.399024623406}},
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
