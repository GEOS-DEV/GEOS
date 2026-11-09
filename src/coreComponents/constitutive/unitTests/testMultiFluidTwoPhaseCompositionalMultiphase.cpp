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
#include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

enum class VISCOSITY_TYPE : int { CONSTANT, LBC };
ENUM_STRINGS( VISCOSITY_TYPE, "Constant", "LohrenzBrayClark" );

template< VISCOSITY_TYPE VISCOSITY >
struct Viscosity {};

template<>
struct Viscosity< VISCOSITY_TYPE::CONSTANT >
{
  using FluidType = CompositionalTwoPhaseConstantViscosity;
};
template<>
struct Viscosity< VISCOSITY_TYPE::LBC >
{
  using FluidType = CompositionalTwoPhaseLohrenzBrayClarkViscosity;
};

template< typename TEST_TYPE >
class MultiFluidCompositionalMultiphaseTestFixture : public FluidModelTest<
    typename Viscosity< std::tuple_element_t< 1, TEST_TYPE >::value >::FluidType,
    getNumberOfComponents< std::tuple_element_t< 2, TEST_TYPE >::value > >
{
public:
  static constexpr EquationOfStateType EquationOfState = std::tuple_element_t< 0, TEST_TYPE >::value;
  static constexpr VISCOSITY_TYPE ViscosityModel = std::tuple_element_t< 1, TEST_TYPE >::value;
  static constexpr TestFluidName FluidName = std::tuple_element_t< 2, TEST_TYPE >::value;
  using CompositionalMultiphaseFluid = typename Viscosity< ViscosityModel >::FluidType;
  using Base = FluidModelTest< CompositionalMultiphaseFluid, getNumberOfComponents< FluidName > >;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-3;

public:
  MultiFluidCompositionalMultiphaseTestFixture()
  {
    Base::createFluid( getFluidName(), []( CompositionalMultiphaseFluid & fluid ){
      fillPhysicalProperties( fluid );
    } );
  }

  ~MultiFluidCompositionalMultiphaseTestFixture() override = default;

  void testNumericalDerivatives( const bool useMass )
  {
    CompositionalMultiphaseFluid * fluid = this->getFluid( this->getFluidName() );

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
  static void fillPhysicalProperties( CompositionalMultiphaseFluid & fluid );
};

template< typename TEST_TYPE >
string MultiFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::getFluidName()
{
  return GEOS_FMT( "fluid_{}_{}_{}",
                   Base::numComp,
                   EnumStrings< EquationOfStateType >::toString( EquationOfState ),
                   EnumStrings< VISCOSITY_TYPE >::toString( ViscosityModel ) );
}

template< typename TEST_TYPE >
void MultiFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::fillPhysicalProperties( CompositionalMultiphaseFluid & fluid )
{
  string_array & phaseNames = fluid.template getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames = {"oil", "gas"};

  string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfState );
  string_array & equationOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  equationOfState = {eosName, eosName};

  dataRepository::Group & group = fluid;
  fillTestFluidComponentProperties< CompositionalMultiphaseFluid, FluidName >( group );
}

template< EquationOfStateType EOS, VISCOSITY_TYPE VISCOSITY, TestFluidName FLUID >
using TestType = std::tuple<
  std::integral_constant< EquationOfStateType, EOS >,
  std::integral_constant< VISCOSITY_TYPE, VISCOSITY >,
  std::integral_constant< TestFluidName, FLUID >
  >;

using TestTypes = ::testing::Types<
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::CONSTANT, TestFluidName::EnhancedOilRecovery >,
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::CONSTANT, TestFluidName::SourNaturalGas >,
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::LBC, TestFluidName::EnhancedOilRecovery >,
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::LBC, TestFluidName::SourNaturalGas >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::CONSTANT, TestFluidName::EnhancedOilRecovery >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::CONSTANT, TestFluidName::SourNaturalGas >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::LBC, TestFluidName::EnhancedOilRecovery >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::LBC, TestFluidName::SourNaturalGas >
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
