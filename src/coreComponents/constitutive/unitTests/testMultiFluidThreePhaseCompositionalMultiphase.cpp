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
 * @file testMultiFluidTwoPhaseCompositionalMultiphase.cpp
 */

#include "FluidModelTest.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ImmiscibleWaterParameters.hpp"
#include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< typename TEST_TYPE >
class MultiFluidCompositionalMultiphaseTestFixture : public FluidModelTest< CompositionalThreePhaseLohrenzBrayClarkViscosity,
                                                                               getNumberOfComponents< std::tuple_element_t< 1, TEST_TYPE >::value >, 3 >
{
public:
  static constexpr EquationOfStateType EquationOfState = std::tuple_element_t< 0, TEST_TYPE >::value;
  static constexpr TestFluidName FluidName = std::tuple_element_t< 1, TEST_TYPE >::value;
  using Base = FluidModelTest< CompositionalThreePhaseLohrenzBrayClarkViscosity, getNumberOfComponents< FluidName >, 3 >;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-3;

public:
  MultiFluidCompositionalMultiphaseTestFixture()
  {
    Base::createFluid( getFluidName(), []( CompositionalThreePhaseLohrenzBrayClarkViscosity & fluid ){
      fillPhysicalProperties( fluid );
    } );
  }

  ~MultiFluidCompositionalMultiphaseTestFixture() override = default;

  void testNumericalDerivatives( const bool useMass )
  {
    CompositionalThreePhaseLohrenzBrayClarkViscosity * fluid = this->getFluid( this->getFluidName() );

    fluid->setMassFlag( useMass );

    real64 constexpr eps = 1.0e-7;

    constexpr real64 pressures[] = { 10.0e5, 50.0e5, 100.0e5, 600.0e5 };
    constexpr real64 temperatures[] = { 15.5, 24.0, 40.0, 80.0 };

    for( auto const & sample : getTestSamples< FluidName >() )
    {
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
  static void fillPhysicalProperties( CompositionalThreePhaseLohrenzBrayClarkViscosity & fluid );
};

template< typename TEST_TYPE >
string MultiFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::getFluidName()
{
  return GEOS_FMT( "fluid_{}_{}",
                   Base::numComp,
                   EnumStrings< EquationOfStateType >::toString( EquationOfState ) );
}

template< typename TEST_TYPE >
void MultiFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::fillPhysicalProperties( CompositionalThreePhaseLohrenzBrayClarkViscosity & fluid )
{
  string_array & phaseNames = fluid.template getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames = {"oil", "gas", "water"};

  string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfState );
  string_array & equationOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  equationOfState = {eosName, eosName, eosName};

  dataRepository::Group & group = fluid;
  fillTestFluidComponentProperties< CompositionalThreePhaseLohrenzBrayClarkViscosity, FluidName >( group );

  // Water properties
  using WaterKeys = typename ImmiscibleWaterParameters::viewKeyStruct;
  real64 & waterDensity = fluid.template getReference< real64 >( WaterKeys::waterDensityString() );
  waterDensity = 1020.0;
  real64 & waterCompressibility = fluid.template getReference< real64 >( WaterKeys::waterCompressibilityString() );
  waterCompressibility = 4.4e-10;
  real64 & waterViscosity = fluid.template getReference< real64 >( WaterKeys::waterViscosityString() );
  waterViscosity = 5.4650e-04;
}

template< EquationOfStateType EOS, TestFluidName FLUID >
using TestType = std::tuple<
  std::integral_constant< EquationOfStateType, EOS >,
  std::integral_constant< TestFluidName, FLUID >
  >;

using TestTypes = ::testing::Types<
  TestType< EquationOfStateType::PengRobinson, TestFluidName::AquiferGasStorage >,
  TestType< EquationOfStateType::SoaveRedlichKwong, TestFluidName::AquiferGasStorage >
  >;

class NameGenerator
{
public:
  template< typename T >
  static std::string GetName( int index )
  {
    return GEOS_FMT( "CompositionalFluid_{}_{}", MultiFluidCompositionalMultiphaseTestFixture< T >::getFluidName(), index );
  }
};

TYPED_TEST_SUITE( MultiFluidCompositionalMultiphaseTestFixture, TestTypes, NameGenerator );

TYPED_TEST( MultiFluidCompositionalMultiphaseTestFixture, numericalDerivativesMolar )
{
  this->testNumericalDerivatives( false );
}

TYPED_TEST( MultiFluidCompositionalMultiphaseTestFixture, numericalDerivativesMass )
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
