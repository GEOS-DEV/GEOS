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

namespace geos
{
namespace testing
{

enum class BrineModelType : int {Phillips, Ezrokhi};
enum class FlashType : int {DuanSun, SpycherPruess};

ENUM_STRINGS( BrineModelType, "Phillips", "Ezrokhi" );
ENUM_STRINGS( FlashType, "DuanSun", "SpycherPruess" );

template< BrineModelType BRINE, bool THERMAL >
struct FluidType {};

template<>
struct FluidType< BrineModelType::Phillips, false >
{
  using type = CO2BrinePhillipsFluid;
  static constexpr const char * brineContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                               "ViscosityFun PhillipsBrineViscosity 0.1";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1";
};
template<>
struct FluidType< BrineModelType::Phillips, true >
{
  using type = CO2BrinePhillipsThermalFluid;
  static constexpr const char * brineContent = "DensityFun PhillipsBrineDensity 1e6 1.5e7 5e4 367.15 369.15 1 0.2\n"
                                               "ViscosityFun PhillipsBrineViscosity 0.1\n"
                                               "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 299.15 369.15 10 0";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "EnthalpyFun CO2Enthalpy 1e6 1.5e7 5e4 367.15 369.15 1";
};
template<>
struct FluidType< BrineModelType::Ezrokhi, false >
{
  using type = CO2BrineEzrokhiFluid;
  static constexpr const char * brineContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                               "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1";
};
template<>
struct FluidType< BrineModelType::Ezrokhi, true >
{
  using type = CO2BrineEzrokhiThermalFluid;
  static constexpr const char * brineContent = "DensityFun EzrokhiBrineDensity 2.01e-6 -6.34e-7 1e-4\n"
                                               "ViscosityFun EzrokhiBrineViscosity 2.42e-7 0 1e-4\n"
                                               "EnthalpyFun BrineEnthalpy 1e6 7.5e7 5e5 299.15 369.15 10 0";
  static constexpr const char * gasContent = "DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 367.15 369.15 1\n"
                                             "EnthalpyFun CO2Enthalpy 1e6 1.5e7 5e4 367.15 369.15 1";
};

template< FlashType FLASH >
struct FlashModel {};

template<>
struct FlashModel< FlashType::DuanSun >
{
  static constexpr const char * flashContent = "FlashModel CO2Solubility 1e6 1.5e7 5e4 367.15 369.15 1 0.15";
};
template<>
struct FlashModel< FlashType::SpycherPruess >
{
  static constexpr const char * flashContent = "FlashModel CO2Solubility 1e6 1.5e7 5e4 367.15 369.15 1 0.15 1.0e-10 SpycherPruess";
};

template< BrineModelType BRINE, FlashType FLASH, bool THERMAL >
class MultiFluidCO2BrineTestFixture : public FluidModelTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >
{
public:
  using CO2BrineFluid = typename FluidType< BRINE, THERMAL >::type;
  using Base = FluidModelTest< CO2BrineFluid, 2, 2 >;

public:
  MultiFluidCO2BrineTestFixture()
  {
    Base::writeTableToFile( pvtGasFileName, FluidType< BRINE, THERMAL >::gasContent );
    Base::writeTableToFile( pvtLiquidFileName, FluidType< BRINE, THERMAL >::brineContent );
    Base::writeTableToFile( pvtFlashFileName, FlashModel< FLASH >::flashContent );

    Base::createFluid( getFluidName(), []( CO2BrineFluid & fluid ){
      fillPhysicalProperties( fluid );
    } );
  }

  ~MultiFluidCO2BrineTestFixture() override
  {
    Base::removeFile( pvtGasFileName );
    Base::removeFile( pvtLiquidFileName );
    Base::removeFile( pvtFlashFileName );
  }

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
  static constexpr const char * pvtGasFileName = "pvtgas.txt";
  static constexpr const char * pvtLiquidFileName = "pvtliquid.txt";
  static constexpr const char * pvtFlashFileName = "co2flash.txt";
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

  auto & phasePVTParaFileNames = group.getReference< path_array >( CO2BrineFluid::viewKeyStruct::phasePVTParaFilesString() );
  phasePVTParaFileNames.emplace_back( Path( pvtGasFileName ) );
  phasePVTParaFileNames.emplace_back( Path( pvtLiquidFileName ) );

  auto & flashModelParaFileName = group.getReference< Path >( CO2BrineFluid::viewKeyStruct::flashModelParaFileString() );
  flashModelParaFileName = pvtFlashFileName;
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
using CO2BrinePhillipsSpycherPruessTest = MultiFluidCO2BrineTestFixture< BrineModelType::Phillips,
                                                                         FlashType::SpycherPruess,
                                                                         false >;
#if !defined(GEOS_DEVICE_COMPILE)
using CO2BrineEzrokhiThermalTest = MultiFluidCO2BrineTestFixture< BrineModelType::Ezrokhi,
                                                                  FlashType::DuanSun,
                                                                  true >;
#endif

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
TEST_F( CO2BrinePhillipsSpycherPruessTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( CO2BrinePhillipsSpycherPruessTest, numericalDerivativesMass )
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

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

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
