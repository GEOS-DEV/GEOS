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

/**
 * @file testMultiFluidCO2Brine.cpp
 */

#include "MultiFluidTest.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

enum class BrineModelType : int {Phillips, Ezrokhi};
enum class FlashType : int {DuanSun, SpycherPruess};

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
class MultiFluidCO2BrineTest : public MultiFluidTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >
{
public:
  using CO2BrineFluid = typename FluidType< BRINE, THERMAL >::type;
  using Base = MultiFluidTest< typename FluidType< BRINE, THERMAL >::type, 2, 2 >;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-4;
public:
  MultiFluidCO2BrineTest()
  {
    Base::writeTableToFile( pvtGasFileName, FluidType< BRINE, THERMAL >::gasContent );
    Base::writeTableToFile( pvtLiquidFileName, FluidType< BRINE, THERMAL >::brineContent );
    Base::writeTableToFile( pvtFlashFileName, FlashModel< FLASH >::flashContent );

    auto & parent = this->m_parent;
    parent.resize( 1 );
    this->m_model = makeCO2BrineFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~MultiFluidCO2BrineTest() override
  {
    Base::removeFile( pvtGasFileName );
    Base::removeFile( pvtLiquidFileName );
    Base::removeFile( pvtFlashFileName );
  }

  // Test numerical derivatives at selected data points
  void testNumericalDerivatives( bool const useMass )
  {
    auto * parent = &(this->getParent());
    auto & fluid = this->getFluid();
    fluid.setMassFlag( useMass );

    real64 constexpr eps = 1.0e-5;

    // Some of the functions are simply table lookups. We need to keep the test points away from
    // the table nodes because the kink in the linear interpolation might cause numerical derivative
    // mismatches. Some of these values have been manually inspected and the differences, although
    // not meeting the tolerance here, are small as expected.
    constexpr real64 temperatures[] = { 367.65, 368.00, 368.75 };
    constexpr real64 pressures[] = { 5.001e6, 7.501e6, 1.201e7 };

    for( real64 const pressure : pressures )
    {
      for( real64 const temperature : temperatures )
      {
        typename Base::TestData data ( pressure, temperature, { 0.3, 0.7 } );
        Base::testNumericalDerivatives( fluid, parent, data, eps, relTol, absTol );
      }
    }
  }

private:
  static CO2BrineFluid * makeCO2BrineFluid( string const & name, Group * parent );
  static constexpr const char * pvtGasFileName = "pvtgas.txt";
  static constexpr const char * pvtLiquidFileName = "pvtliquid.txt";
  static constexpr const char * pvtFlashFileName = "co2flash.txt";
};

template< BrineModelType BRINE, FlashType FLASH, bool THERMAL >
typename MultiFluidCO2BrineTest< BRINE, FLASH, THERMAL >::CO2BrineFluid *
MultiFluidCO2BrineTest< BRINE, FLASH, THERMAL >::makeCO2BrineFluid( string const & name, Group * parent )
{
  CO2BrineFluid & co2BrineFluid = parent->registerGroup< CO2BrineFluid >( name );

  Group & fluid = co2BrineFluid;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  fill< 2 >( phaseNames, {"gas", "liquid"} );

  auto & compNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
  fill< 2 >( compNames, {"co2", "water"} );

  auto & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
  fill< 2 >( molarWeight, {44e-3, 18e-3} );

  auto & phasePVTParaFileNames = fluid.getReference< path_array >( CO2BrinePhillipsFluid::viewKeyStruct::phasePVTParaFilesString() );
  fill< 2 >( phasePVTParaFileNames, {pvtGasFileName, pvtLiquidFileName} );

  auto & flashModelParaFileName = fluid.getReference< Path >( CO2BrinePhillipsFluid::viewKeyStruct::flashModelParaFileString() );
  flashModelParaFileName = pvtFlashFileName;

  co2BrineFluid.postProcessInputRecursive();
  return &co2BrineFluid;
}

using MultiFluidCO2BrineBrinePhillipsTest = MultiFluidCO2BrineTest< BrineModelType::Phillips, FlashType::DuanSun, false >;
TEST_F( MultiFluidCO2BrineBrinePhillipsTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrinePhillipsTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

using MultiFluidCO2BrineBrineEzrokhiTest = MultiFluidCO2BrineTest< BrineModelType::Ezrokhi, FlashType::DuanSun, false >;
TEST_F( MultiFluidCO2BrineBrineEzrokhiTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrineEzrokhiTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

using MultiFluidCO2BrineBrinePhillipsThermalTest = MultiFluidCO2BrineTest< BrineModelType::Phillips, FlashType::DuanSun, true >;
TEST_F( MultiFluidCO2BrineBrinePhillipsThermalTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}

TEST_F( MultiFluidCO2BrineBrinePhillipsThermalTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

using MultiFluidCO2BrineBrinePhillipsSpycherPruessTest = MultiFluidCO2BrineTest< BrineModelType::Phillips, FlashType::SpycherPruess, false >;
TEST_F( MultiFluidCO2BrineBrinePhillipsSpycherPruessTest, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( MultiFluidCO2BrineBrinePhillipsSpycherPruessTest, numericalDerivativesMass )
{
  testNumericalDerivatives( true );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
