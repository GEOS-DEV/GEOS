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
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
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

protected:
  void resetFluid( MultiFluidBase & fluid ) const override;

private:
  static CompositionalMultiphaseFluid * makeFluid( string const & name, Group * parent );
};

template< EOS_TYPE EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
void MultiFluidCompositionalMultiphaseTest< EOS, VISCOSITY, NUM_COMP >::
resetFluid( MultiFluidBase & fluid ) const
{
  auto & kValues = fluid.getReference< fields::multifluid::array4dLayoutPhaseComp >( "kValues" );
  kValues.zero();
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

template< typename FLUID_TYPE >
struct Fluid< FLUID_TYPE, 5 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & phaseNames = fluid.getReference< string_array >( string( MultiFluidBase::viewKeyStruct::phaseNamesString()) );
    fill< 2 >( phaseNames, {"oil", "gas"} );

    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 5 >( componentNames, {"CO2", "N2", "C1", "C2", "C3"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 5 >( molarWeight, {44.0098e-3, 28.0135e-3, 16.0428e-3, 30.0700e-3, 44.1000e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalPressureString() );
    fill< 5 >( criticalPressure, {73.77300e5, 33.95800e5, 45.99200e5, 48.71800e5, 42.46000e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 5 >( criticalTemperature, {304.1280, 126.1920, 190.5640, 305.3300, 369.8000} );
    array1d< real64 > & criticalVolume = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalVolumeString() );
    fill< 5 >( criticalVolume, {0.094472e-3, 0.089800e-3, 0.099200e-3, 0.148300e-3, 0.203000e-3} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentAcentricFactorString() );
    fill< 5 >( acentricFactor, {0.223000, 0.037200, 0.010400, 0.099100, 0.145400} );
    array1d< real64 > & volumeShift = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentVolumeShiftString() );
    fill< 5 >( volumeShift, {1.845465e-01, -1.283880e-01, 9.225800e-02, 6.458060e-02, -5.000000e-02} );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( FLUID_TYPE::viewKeyStruct::componentBinaryCoeffString() );
    fillBinaryCoeffs< 5 >( binaryCoeff, {0.000000, 0.100000, 0.139000, 0.120000, 0.030000, 0.032000, 0.030000, 0.000000, 0.000000, 0.000000} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 1, 5 );
    fill< 5 >( samples[0], {0.13683, 0.02680, 0.50376, 0.02848, 0.30413} );
    //fill< 5 >( samples[0], {0.03683, 0.02680, 0.90376, 0.02848, 0.00413} );
  }
};

template< typename FLUID_TYPE >
struct Fluid< FLUID_TYPE, 9 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & phaseNames = fluid.getReference< string_array >( string( MultiFluidBase::viewKeyStruct::phaseNamesString()) );
    fill< 2 >( phaseNames, {"oil", "gas"} );

    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 9 >( componentNames, {"CO2", "N2", "C1", "C2", "C3", "C4", "C5", "C6", "C7+"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 9 >( molarWeight, {44.01e-3, 28.01e-3, 16.04e-3, 30.07e-3, 44.1e-3, 58.12e-3, 72.15e-3, 84e-3, 173e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalPressureString() );
    fill< 9 >( criticalPressure, {73.8659e5, 33.9439e5, 46.0421e5, 48.8387e5, 42.4552e5, 37.47e5, 33.5892e5, 30.1037e5, 20.549e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 9 >( criticalTemperature, {304.7, 126.2, 190.6, 305.43, 369.8, 419.5, 465.9, 507.5, 678.8} );
    array1d< real64 > & criticalVolume = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentCriticalVolumeString() );
    fill< 9 >( criticalVolume, {9.3999e-05, 9.0001e-05, 9.7999e-05, 1.4800e-04, 2.0000e-04, 2.5800e-04, 3.1000e-04, 3.5100e-04, 6.8243e-04} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentAcentricFactorString() );
    fill< 9 >( acentricFactor, {00.225, 0.04, 0.013, 0.0986, 0.1524, 0.1956, 0.2413, 0.299, 0.5618} );
    array1d< real64 > & volumeShift = fluid.getReference< array1d< real64 > >( FLUID_TYPE::viewKeyStruct::componentVolumeShiftString() );
    fill< 9 >( volumeShift, {-0.04958, -0.136012, -0.1486264, -0.10863408, -0.08349872, -0.06331568, -0.04196464, -0.0150072, 0.0000} );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( FLUID_TYPE::viewKeyStruct::componentBinaryCoeffString() );
    fillBinaryCoeffs< 9 >( binaryCoeff, {
      1.0000e-02,
      0.0000e+00, 3.7320e-03,
      0.0000e+00, 1.0000e-02, 0.0000e+00,
      0.0000e+00, 1.0000e-02, 0.0000e+00, 0.0000e+00,
      0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
      0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
      1.0000e-02, 0.0000e+00, 2.8000e-02, 1.0000e-02, 1.0000e-02, 0.0000e+00, 0.0000e+00,
      1.0000e-02, 0.0000e+00, 4.5320e-02, 1.0000e-02, 1.0000e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 1, 9 );
    fill< 9 >( samples[0], {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920} );
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

/*
   using PengRobinsonConstantViscosity2Test = MultiFluidCompositionalMultiphaseTest< EOS_TYPE::PR, VISCOSITY_TYPE::CONSTANT, 2 >;
   TEST_F( PengRobinsonConstantViscosity2Test, numericalDerivatives )
   {
   auto & fluid = getFluid();
   array2d< real64 > samples;
   Fluid< CompositionalMultiphaseFluid, numComp >::getSamples( samples );
   integer const sampleCount = samples.size( 0 );

   real64 constexpr eps = 1.0e-5;

   constexpr real64 pressures[] = { 1.0e7 };
   //constexpr real64 temperatures[] = { 15.5 };
   //constexpr real64 pressures[] = { 1.0e5, 10.0e5, 100.0e5, 600.0e5 };
   constexpr real64 temperatures[] = { 15.5, 25.0, 40.0, 80.0 };

   fluid.setMassFlag( false );
   for( real64 const pressure : pressures )
   {
    for( real64 const temperature : temperatures )
    {
      for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
      {
        std::cout << "P: " << pressure << " T: " << temperature << " S: " << sampleIndex << "\n";
        TestData data ( pressure, units::convertCToK( temperature ), samples[sampleIndex].toSliceConst() );
        testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
      }
    }
   }
   }

   using SoaveRedlichKwongConstantViscosity5Test = MultiFluidCompositionalMultiphaseTest< EOS_TYPE::SRK, VISCOSITY_TYPE::CONSTANT, 5 >;
   TEST_F( SoaveRedlichKwongConstantViscosity5Test, numericalDerivatives )
   {
   auto & fluid = getFluid();
   array2d< real64 > samples;
   Fluid< CompositionalMultiphaseFluid, numComp >::getSamples( samples );
   integer const sampleCount = samples.size( 0 );

   real64 constexpr eps = 1.0e-5;

   constexpr real64 pressures[] = { 206e5 };
   constexpr real64 temperatures[] = { 80.0 };
   //constexpr real64 pressures[] = { 1.0e5, 10.0e5, 100.0e5, 600.0e5 };
   //constexpr real64 temperatures[] = { 15.5, 25.0, 40.0, 80.0 };

   fluid.setMassFlag( false );
   for( real64 const pressure : pressures )
   {
    for( real64 const temperature : temperatures )
    {
      for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
      {
        std::cout << "P: " << pressure << " T: " << temperature << " S: " << sampleIndex << "\n";
        TestData data ( pressure, units::convertCToK( temperature ), samples[sampleIndex].toSliceConst() );
        testNumericalDerivatives( fluid, &getParent(), data, eps, relTol, absTol );
      }
    }
   }
   }
 */

using PengRobinsonConstantViscosity9Test = MultiFluidCompositionalMultiphaseTest< EOS_TYPE::PR, VISCOSITY_TYPE::CONSTANT, 9 >;
TEST_F( PengRobinsonConstantViscosity9Test, numericalDerivatives )
{
  auto & fluid = getFluid();
  array2d< real64 > samples;
  Fluid< CompositionalMultiphaseFluid, numComp >::getSamples( samples );
  integer const sampleCount = samples.size( 0 );

  real64 constexpr eps = 1.0e-5;

  constexpr real64 pressures[] = { 1.0e5, 10.0e5, 100.0e5, 600.0e5 };
  constexpr real64 temperatures[] = { 15.5, 25.0, 40.0, 80.0 };

  fluid.setMassFlag( false );
  for( real64 const pressure : pressures )
  {
    for( real64 const temperature : temperatures )
    {
      for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
      {
        std::cout << "P: " << pressure << " T: " << temperature << " S: " << sampleIndex << "\n";
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
