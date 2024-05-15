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
 * @file testMultiFluidCompositionalMultiphasePVTPackage.cpp
 */

#include "MultiFluidTest.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;

enum class EOS_TYPE : int { PR, SRK };

template< EOS_TYPE EOS >
struct EquationOfStateName {};
template<> struct EquationOfStateName< EOS_TYPE::PR > { static constexpr char const * name = "PR"; };
template<> struct EquationOfStateName< EOS_TYPE::SRK > { static constexpr char const * name = "SRK"; };


namespace geos
{
namespace testing
{

template< >
struct UsePVTPackage< CompositionalMultiphaseFluidPVTPackage >
{
  static constexpr bool value = true;
};

} // namespace testing

} // namespace geos

template< EOS_TYPE EOS, integer NUM_COMP >
class MultiFluidCompositionalMultiphasePVTPackageTest : public MultiFluidTest< CompositionalMultiphaseFluidPVTPackage, 2, NUM_COMP >
{
public:
  static constexpr real64 relTol = 1e-4;
  static constexpr real64 absTol = 1e-4;
public:
  MultiFluidCompositionalMultiphasePVTPackageTest()
  {
    auto & parent = this->m_parent;
    parent.resize( 1 );
    this->m_model = makeFluid( "fluid", &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~MultiFluidCompositionalMultiphasePVTPackageTest() override = default;

private:
  static CompositionalMultiphaseFluidPVTPackage * makeFluid( string const & name, Group * parent );
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

template< integer NUM_COMP >
struct Fluid
{};

template< >
struct Fluid< 2 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & phaseNames = fluid.getReference< string_array >( string( MultiFluidBase::viewKeyStruct::phaseNamesString()) );
    fill< 2 >( phaseNames, {"oil", "gas"} );

    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 2 >( componentNames, {"C1", "CO2"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 2 >( molarWeight, {16.043e-3, 44.01e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentCriticalPressureString() );
    fill< 2 >( criticalPressure, {46.04208e5, 73.865925e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 2 >( criticalTemperature, {190.6, 304.7} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentAcentricFactorString() );
    fill< 2 >( acentricFactor, {0.013, 0.225 } );
    array1d< real64 > & volumeShift = fluid.getReference< array1d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentVolumeShiftString() );
    fill< 2 >( volumeShift, {-0.1486264, -0.04958 } );
    array2d< real64 > & binaryCoeff = fluid.getReference< array2d< real64 > >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::componentBinaryCoeffString() );
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

template< EOS_TYPE EOS, integer NUM_COMP >
CompositionalMultiphaseFluidPVTPackage *
MultiFluidCompositionalMultiphasePVTPackageTest< EOS, NUM_COMP >::
makeFluid( string const & name, Group * parent )
{
  CompositionalMultiphaseFluidPVTPackage & fluid = parent->registerGroup< CompositionalMultiphaseFluidPVTPackage >( name );

  Fluid< NUM_COMP >::fillProperties( fluid );

  string_array & equationOfState = fluid.getReference< string_array >( CompositionalMultiphaseFluidPVTPackage::viewKeyStruct::equationsOfStateString() );
  fill< 2 >( equationOfState, {EquationOfStateName< EOS >::name, EquationOfStateName< EOS >::name} );

  fluid.postProcessInputRecursive();
  return &fluid;
}

using PengRobinsonConstantViscosity2Test = MultiFluidCompositionalMultiphasePVTPackageTest< EOS_TYPE::PR, 2 >;

TEST_F( PengRobinsonConstantViscosity2Test, numericalDerivatives )
{
  auto & fluid = getFluid();
  array2d< real64 > samples;
  Fluid< numComp >::getSamples( samples );
  integer const sampleCount = samples.size( 0 );

  real64 constexpr eps = 1.0e-7;

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

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
