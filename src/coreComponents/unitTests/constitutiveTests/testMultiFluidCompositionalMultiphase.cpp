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
 * @file testMultiFluidCompositionalMultiphase.cpp
 */

#include "MultiFluidTest.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

enum class EOS_TYPE : int { PR, SRK };
enum class VISCOSITY_TYPE : int { CONSTANT, LBC };

template< EOS_TYPE eos, VISCOSITY_TYPE viscosity >
struct FluidType {};

template<>
struct FluidType< EOS_TYPE::PR, VISCOSITY_TYPE::CONSTANT >
{
  using type = CompositionalTwoPhasePengRobinsonConstantViscosity;
};
template<>
struct FluidType< EOS_TYPE::PR, VISCOSITY_TYPE::LBC >
{
  using type = CompositionalTwoPhasePengRobinsonLBCViscosity;
};
template<>
struct FluidType< EOS_TYPE::SRK, VISCOSITY_TYPE::CONSTANT >
{
  using type = CompositionalTwoPhasePengRobinsonConstantViscosity;
};
template<>
struct FluidType< EOS_TYPE::SRK, VISCOSITY_TYPE::LBC >
{
  using type = CompositionalTwoPhaseSoaveRedlichKwongConstantViscosity;
};

template< EOS_TYPE EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
class MultiFluidCompositionalMultiphaseTest : public MultiFluidTest< typename FluidType< EOS, VISCOSITY >::type, 2, NUM_COMP >
{
public:
  using CompositionalMultiphaseFluid = typename FluidType< EOS, VISCOSITY >::type;
  static constexpr real64 relTol = 1e-4;
  static constexpr real64 absTol = 1e-4;
public:
  MultiFluidCompositionalMultiphaseTest()
  {
    auto & parent = this->m_parent;
    parent.resize( 1 );
    this->m_model = makeFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~MultiFluidCompositionalMultiphaseTest() override = default;

private:
  static CompositionalMultiphaseFluid * makeFluid( string const & name, Group * parent );
};

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

template< typename FLUID_TYPE, integer NUM_COMP >
struct Fluid
{};

template< typename FLUID_TYPE >
struct Fluid< FLUID_TYPE, 2 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & phaseNames = fluid.getReference< string_array >( string( MultiFluidBase::viewKeyStruct::phaseNamesString()) );
    fill< 2 >( phaseNames, {"oil", "gas"} );

    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 2 >( componentNames, {"C1", "CO2"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 2 >( molarWeight, {16.043e-3, 44.01e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalPressureString() );
    fill< 2 >( criticalPressure, {46.04208e5, 73.865925e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 2 >( criticalTemperature, {190.6, 304.7} );
    array1d< real64 > & criticalVolume = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalVolumeString() );
    fill< 2 >( criticalVolume, {0.098e-3, 0.094e-3} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentAcentricFactorString() );
    fill< 2 >( acentricFactor, {0.013, 0.225 } );
    array1d< real64 > & volumeShift = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentVolumeShiftString() );
    fill< 2 >( volumeShift, {-0.1486264, -0.04958 } );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( FLUID_TYPE::viewKeyStruct::componentBinaryCoeffString() );
    fillBinaryCoeffs< 2 >( binaryCoeff, {0.1896} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 1, 2 );
    fill< 2 >( samples[0], {0.500000, 0.500000} );
    //fill< 2 >( samples[1], {0.100000, 0.900000} );
    //fill< 2 >( samples[2], {0.000010, 0.999990} );
    //fill< 2 >( samples[3], {0.999990, 0.000010} );
    //fill< 2 >( samples[4], {1.000000, 0.000000} );
  }
};

template< EOS_TYPE EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
typename MultiFluidCompositionalMultiphaseTest< EOS, VISCOSITY, NUM_COMP >::CompositionalMultiphaseFluid *
MultiFluidCompositionalMultiphaseTest< EOS, VISCOSITY, NUM_COMP >::
makeFluid( string const & name, Group * parent )
{
  CompositionalMultiphaseFluid & fluid = parent->registerGroup< CompositionalMultiphaseFluid >( name );

  Fluid< CompositionalMultiphaseFluid, NUM_COMP >::fillProperties( fluid );

  fluid.postProcessInputRecursive();
  return &fluid;
}

using PengRobinsonConstantViscosity2Test = MultiFluidCompositionalMultiphaseTest< EOS_TYPE::PR, VISCOSITY_TYPE::CONSTANT, 2 >;

TEST_F( PengRobinsonConstantViscosity2Test, numericalDerivatives )
{
  auto & fluid = getFluid();
  array2d< real64 > samples;
  Fluid< CompositionalMultiphaseFluid, numComp >::getSamples( samples );
  integer const sampleCount = samples.size( 1 );

  real64 constexpr eps = 1.0e-7;

  constexpr real64 pressures[] = { 1.0e5 };
  constexpr real64 temperatures[] = { 15.5 };
  //constexpr real64 pressures[] = { 1.0e5, 10.0e5, 100.0e5, 600.0e5 };
  //constexpr real64 temperatures[] = { 15.5, 25.0, 40.0, 80.0 };

  fluid.setMassFlag( false );
  for( real64 const pressure : pressures )
  {
    for( real64 const temperature : temperatures )
    {
      for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
      {
        TestData data ( pressure, units::convertCToK( temperature ), samples[sampleIndex].toSliceConst() );
        testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
      }
    }
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
