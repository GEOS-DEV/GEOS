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

#include "FluidModelTest.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PressureTemperatureCoordinates.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NUM_COMP >
struct FluidData;

template<>
struct FluidData< 3 >
{
  static void populateFluid( dataRepository::Group & fluid )
  {
    using Keys = CompositionalTwoPhasePhillipsBrine::viewKeyStruct;

    string_array & phaseNames = fluid.getReference< string_array >( Keys::phaseNamesString() );
    phaseNames = {"liquid", "gas"};

    string_array & equationsOfState = fluid.getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
    equationsOfState.resize( 2 );
    equationsOfState[0] = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::SoreideWhitson );
    equationsOfState[1] = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );

    string_array & componentNames = fluid.getReference< string_array >( Keys::componentNamesString() );
    componentNames = {"CO2", "C1", "H2O"};

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( Keys::componentMolarWeightString() );
    TestFluid< 3 >::createArray( molarWeight, Feed< 3 >{4.40095e-02, 1.60425e-02, 1.80153e-02} );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( Keys::componentCriticalPressureString() );
    TestFluid< 3 >::createArray( criticalPressure, Feed< 3 >{7.37730e+06, 4.59920e+06, 2.20640e+07} );
    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( Keys::componentCriticalTemperatureString() );
    TestFluid< 3 >::createArray( criticalTemperature, Feed< 3 >{3.04128e+02, 1.90564e+02, 6.47096e+02} );
    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( Keys::componentAcentricFactorString() );
    TestFluid< 3 >::createArray( acentricFactor, Feed< 3 >{2.23940e-01, 1.14200e-02, 3.44300e-01} );
    array1d< real64 > & volumeShift = fluid.getReference< array1d< real64 > >( CriticalVolume::viewKeyStruct::componentCriticalVolumeString() );
    TestFluid< 3 >::createArray( volumeShift, Feed< 3 >{9.41185e-05, 9.86278e-05, 5.59480e-05} );
  }
};

template< integer NUM_COMP >
class SoreideWhitsonTestFixture : public FluidModelTest< CompositionalTwoPhasePhillipsBrine, NUM_COMP >
{
public:
  using Base = FluidModelTest< CompositionalTwoPhasePhillipsBrine, NUM_COMP >;
  using FluidModel = typename Base::FluidModel;
  using FluidWrapper = typename Base::FluidWrapper;
  using TestPoint = typename Base::TestPoint;

public:
  SoreideWhitsonTestFixture()
  {
    Base::createFluid( "testFluid", []( FluidModel & model ){
      FluidData< Base::numComp >::populateFluid( model );
      populatePressureCoordinates( model );
      populateTemperatureCoordinates( model );
    } );
  }
  ~SoreideWhitsonTestFixture() override = default;

protected:
  static void populatePressureCoordinates( dataRepository::Group & fluid )
  {
    array1d< real64 > & pressureCoordinates = fluid.getReference< array1d< real64 > >( PressureTemperatureCoordinates::viewKeyStruct::pressureCoordinatesString() );
    Base::populateLogScale( pressureCoordinates, 1.0e5, 600.0e5, 20 );
  }

  static void populateTemperatureCoordinates( dataRepository::Group & fluid )
  {
    array1d< real64 > & temperatureCoordinates = fluid.getReference< array1d< real64 > >( PressureTemperatureCoordinates::viewKeyStruct::temperatureCoordinatesString() );
    Base::populateLinearScale( temperatureCoordinates, 288.15, 393.15, 21 );
  }
};

using SoreideWhitsonTestFixture3 = SoreideWhitsonTestFixture< 3 >;

TEST_F( SoreideWhitsonTestFixture3, testSolubility )
{
  FluidModel * fluid = getFluid( "testFluid" );
  EXPECT_TRUE( fluid != nullptr );

  m_parent.initialize();
  m_parent.initializePostInitialConditions();

  //TestPoint const data{ 50.e5, 50.1+273.15, { 0.01, 0.02, 0.97 } };
  TestPoint const data{ 50.e5, 50.1+273.15, { 0.00294844, 0.000553996, 0.9964981 } };
  Base::testNumericalDerivatives( fluid, data );
}

} // testing

} // geos
