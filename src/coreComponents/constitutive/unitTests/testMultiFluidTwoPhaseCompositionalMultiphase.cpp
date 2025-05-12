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

template< typename FluidModel, integer NUM_COMP >
struct FluidData
{};

template< typename TEST_TYPE >
class MultiFluidCompositionalMultiphaseTestFixture : public FluidModelTest<
    typename Viscosity< std::tuple_element_t< 1, TEST_TYPE >::value >::FluidType,
    std::tuple_element_t< 2, TEST_TYPE >::value >
{
public:
  static constexpr EquationOfStateType EquationOfState = std::tuple_element_t< 0, TEST_TYPE >::value;
  static constexpr VISCOSITY_TYPE ViscosityModel = std::tuple_element_t< 1, TEST_TYPE >::value;
  using CompositionalMultiphaseFluid = typename Viscosity< ViscosityModel >::FluidType;
  using Base = FluidModelTest< CompositionalMultiphaseFluid, std::tuple_element_t< 2, TEST_TYPE >::value >;
  static constexpr real64 relTol = 1.0e-4;
  static constexpr real64 absTol = 1.0e-4;

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

    array2d< real64 > samples;
    FluidData< CompositionalMultiphaseFluid, Base::numComp >::getSamples( samples );
    integer const sampleCount = samples.size( 0 );
    Feed< Base::numComp > sample;

    real64 constexpr eps = 1.0e-7;

    constexpr real64 pressures[] = { 1.0e5, 50.0e5, 100.0e5, 600.0e5 };
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
};

template< typename TEST_TYPE >
string MultiFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::getFluidName()
{
  return GEOS_FMT( "fluid_{}_{}_{}",
                   Base::numComp,
                   EnumStrings< EquationOfStateType >::toString( EquationOfState ),
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
  static void fillProperties( dataRepository::Group & fluid )
  {
    using Keys = typename FluidModel::viewKeyStruct;

    string_array & componentNames = fluid.getReference< string_array >( Keys::componentNamesString() );
    componentNames = {"N2", "C10", "C20", "H20"};

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( Keys::componentMolarWeightString() );
    TestFluid< 4 >::createArray( molarWeight, Feed< 4 >{28e-3, 134e-3, 275e-3, 18e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( Keys::componentCriticalPressureString() );
    TestFluid< 4 >::createArray( criticalPressure, Feed< 4 >{34e5, 25.3e5, 14.6e5, 220.5e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( Keys::componentCriticalTemperatureString() );
    TestFluid< 4 >::createArray( criticalTemperature, Feed< 4 >{126.2, 622.0, 782.0, 647.0} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( Keys::componentAcentricFactorString() );
    TestFluid< 4 >::createArray( acentricFactor, Feed< 4 >{0.04, 0.443, 0.816, 0.344} );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( Keys::componentBinaryCoeffString() );
    fillBinaryCoeffs< 4 >( binaryCoeff, {0.0, 0.1, 0.0, 0.0, 0.0, 0.0} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 3, 4 );
    populateArray< 4 >( samples[0], {0.099, 0.300, 0.600, 0.001} );
    populateArray< 4 >( samples[1], {0.350, 0.350, 0.200, 0.100} );
    populateArray< 4 >( samples[2], {0.050, 0.050, 0.100, 0.800} );
  }
};

template< typename FluidModel >
struct FluidData< FluidModel, 5 >
{
  static void fillProperties( dataRepository::Group & fluid )
  {
    using Keys = typename FluidModel::viewKeyStruct;

    string_array & componentNames = fluid.getReference< string_array >( Keys::componentNamesString() );
    componentNames = {"CO2", "N2", "C1", "C2", "C4"};

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( Keys::componentMolarWeightString() );
    TestFluid< 5 >::createArray( molarWeight, Feed< 5 >{44.0098e-3, 28.0135e-3, 16.0428e-3, 30.0700e-3, 82.4191e-3} );

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
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 1, 5 );
    populateArray< 5 >( samples[0], {0.050, 0.150, 0.550, 0.150, 0.100} );
  }
};

template< typename TEST_TYPE >
void MultiFluidCompositionalMultiphaseTestFixture< TEST_TYPE >::fillPhysicalProperties( CompositionalMultiphaseFluid & fluid )
{
  string_array & phaseNames = fluid.template getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames = {"oil", "gas"};

  string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfState );
  string_array & equationOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  equationOfState = {eosName, eosName};

  dataRepository::Group & group = fluid;
  FluidData< CompositionalMultiphaseFluid, Base::numComp >::fillProperties( group );
}

template< EquationOfStateType EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
using TestType = std::tuple<
  std::integral_constant< EquationOfStateType, EOS >,
  std::integral_constant< VISCOSITY_TYPE, VISCOSITY >,
  std::integral_constant< integer, NUM_COMP >
  >;

using TestTypes = ::testing::Types<
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::CONSTANT, 4 >,
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::CONSTANT, 5 >,
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::LBC, 4 >,
  TestType< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::LBC, 5 >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::CONSTANT, 4 >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::CONSTANT, 5 >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::LBC, 4 >,
  TestType< EquationOfStateType::SoaveRedlichKwong, VISCOSITY_TYPE::LBC, 5 >
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

TYPED_TEST( MultiFluidCompositionalMultiphaseTestFixture, numericalDerivatives )
{
  this->testNumericalDerivatives( false );
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
