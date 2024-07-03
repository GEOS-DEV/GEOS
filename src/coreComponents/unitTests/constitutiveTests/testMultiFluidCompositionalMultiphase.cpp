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
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

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
struct Fluid
{};

template< EquationOfStateType EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
class MultiFluidCompositionalMultiphaseTest : public MultiFluidTest< typename Viscosity< VISCOSITY >::FluidType, 2, NUM_COMP >
{
public:
  using FluidModel = typename Viscosity< VISCOSITY >::FluidType;
  using Base = MultiFluidTest< FluidModel, 2, NUM_COMP >;
  static constexpr real64 relTol = 1.0e-4;
public:
  MultiFluidCompositionalMultiphaseTest()
  {
    auto & parent = this->m_parent;
    parent.resize( 1 );

    string fluidName = GEOS_FMT( "fluid{}{}{}",
                                 NUM_COMP,
                                 EnumStrings< EquationOfStateType >::toString( EOS ),
                                 EnumStrings< VISCOSITY_TYPE >::toString( VISCOSITY ) );
    this->m_model = makeFluid( fluidName, &parent );

    parent.initialize();
    parent.initializePostInitialConditions();
  }

  ~MultiFluidCompositionalMultiphaseTest() override = default;

  void testNumericalDerivatives( const bool useMass )
  {
    auto & fluid = this->getFluid();
    fluid.setMassFlag( useMass );

    auto * parent = &(this->getParent());

    array2d< real64 > samples;
    Fluid< FluidModel, Base::numComp >::getSamples( samples );
    integer const sampleCount = samples.size( 0 );

    real64 constexpr eps = 1.0e-7;

    constexpr real64 pressures[] = { 1.0e5, 50.0e5, 100.0e5, 600.0e5 };
    constexpr real64 temperatures[] = { 15.5, 24.0, 40.0, 80.0 };

    for( integer sampleIndex = 0; sampleIndex < sampleCount; ++sampleIndex )
    {
      for( real64 const pressure : pressures )
      {
        for( real64 const temperature : temperatures )
        {
          typename Base::TestData data ( pressure, units::convertCToK( temperature ), samples[sampleIndex].toSliceConst() );
          Base::testNumericalDerivatives( fluid, parent, data, eps, relTol );
        }
      }
    }
  }

private:
  static FluidModel * makeFluid( string const & name, Group * parent );
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

template< typename FluidModel >
struct Fluid< FluidModel, 4 >
{
  static void fillProperties( Group & fluid )
  {
    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    fill< 4 >( componentNames, {"N2", "C10", "C20", "H20"} );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    fill< 4 >( molarWeight, {28e-3, 134e-3, 275e-3, 18e-3} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalPressureString() );
    fill< 4 >( criticalPressure, {34e5, 25.3e5, 14.6e5, 220.5e5} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalTemperatureString() );
    fill< 4 >( criticalTemperature, {126.2, 622.0, 782.0, 647.0} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentAcentricFactorString() );
    fill< 4 >( acentricFactor, {0.04, 0.443, 0.816, 0.344} );
  }

  static void getSamples( array2d< real64 > & samples )
  {
    samples.resize( 2, 4 );
    fill< 4 >( samples[0], {0.099, 0.300, 0.600, 0.001} );
    fill< 4 >( samples[1], {0.000, 0.000, 0.000, 1.000} );
  }
};

template< EquationOfStateType EOS, VISCOSITY_TYPE VISCOSITY, integer NUM_COMP >
typename MultiFluidCompositionalMultiphaseTest< EOS, VISCOSITY, NUM_COMP >::FluidModel *
MultiFluidCompositionalMultiphaseTest< EOS, VISCOSITY, NUM_COMP >::
makeFluid( string const & name, Group * parent )
{
  FluidModel & fluid = parent->registerGroup< FluidModel >( name );

  string_array & phaseNames = fluid.template getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  fill< 2 >( phaseNames, {"oil", "gas"} );

  string const eosName = EnumStrings< EquationOfStateType >::toString( EOS );
  string_array & equationOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  fill< 2 >( equationOfState, {eosName, eosName} );

  Fluid< FluidModel, NUM_COMP >::fillProperties( fluid );

  fluid.postInputInitializationRecursive();
  return &fluid;
}

using PengRobinson4Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::CONSTANT, 4 >;
using PengRobinsonLBC4Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::PengRobinson, VISCOSITY_TYPE::LBC, 4 >;
//using SoaveRedlichKwong4Test = MultiFluidCompositionalMultiphaseTest< EquationOfStateType::So, 4 >;

TEST_F( PengRobinson4Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}
TEST_F( PengRobinsonLBC4Test, numericalDerivativesMolar )
{
  testNumericalDerivatives( false );
}

/**
   TEST_F( PengRobinson4Test, numericalDerivativesMass )
   {
   testNumericalDerivatives( true );
   }

   TEST_F( SoaveRedlichKwong4Test, numericalDerivativesMolar )
   {
   testNumericalDerivatives( false );
   }
   TEST_F( SoaveRedlichKwong4Test, numericalDerivativesMass )
   {
   testNumericalDerivatives( true );
   }
 */

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
