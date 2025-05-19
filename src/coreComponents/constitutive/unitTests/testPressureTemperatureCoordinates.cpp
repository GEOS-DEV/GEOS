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

// Source includes
#include "constitutive/fluid/multifluid/compositional/parameters/PressureTemperatureCoordinates.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "TestFluid.hpp"

#include <conduit.hpp>
#include <gtest/gtest.h>

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

class MockFluid : public MultiFluidBase
{
public:
  MockFluid( string const & name,
             Group * const parent ):
    MultiFluidBase( name, parent )
  {}
  ~MockFluid() override = default;

  string getCatalogName() const override { return "FluidModel"; }
  integer getWaterPhaseIndex() const override { return 1; }
  void checkTablesParameters( real64, real64 ) const override {}
};

class PressureTemperatureCoordinatesTestFixture : public ::testing::Test
{
public:
  using Keys = PressureTemperatureCoordinates::viewKeyStruct;
  static constexpr integer numComps = 4;
  static constexpr real64 absTol = 1.0e-10;

public:
  PressureTemperatureCoordinatesTestFixture();
  ~PressureTemperatureCoordinatesTestFixture() override = default;

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  std::unique_ptr< MockFluid > m_fluid{};
  std::unique_ptr< TestFluid< numComps > > m_testFluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
  PressureTemperatureCoordinates * m_coordinates{};
};

PressureTemperatureCoordinatesTestFixture::PressureTemperatureCoordinatesTestFixture():
  m_parent( "parent", m_node ),
  m_fluid( std::make_unique< MockFluid >( "fluid", &m_parent )),
  m_testFluid( TestFluid< numComps >::create( {Fluid::CO2, Fluid::H2O, Fluid::C1, Fluid::C5} ))
{
  m_parameters = PressureTemperatureCoordinates::create( std::move( m_parameters ));
  m_coordinates = dynamic_cast< PressureTemperatureCoordinates * >(m_parameters.get());
  GEOS_ERROR_IF( m_coordinates == nullptr, "Cannot create PressureTemperatureCoordinates" );
}

TEST_F( PressureTemperatureCoordinatesTestFixture, testPressureCoordinates )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get());

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::pressureCoordinatesString()));

  auto & pressureCoordinates = m_fluid->getWrapper< array1d< real64 > >( Keys::pressureCoordinatesString() ).reference();

  // Can be empty
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  // Should have at least 2 elements
  pressureCoordinates.resize( 1 );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );

  // Should be increasing
  pressureCoordinates.resize( 3 );
  TestFluid< 3 >::populateArray( pressureCoordinates, Feed< 3 >{ 100e5, 99e5, 98e5 } );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );

  // Should be strictly increasing
  pressureCoordinates.resize( 3 );
  TestFluid< 3 >::populateArray( pressureCoordinates, Feed< 3 >{ 100e5, 101e5, 101e5 } );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
}

TEST_F( PressureTemperatureCoordinatesTestFixture, testTemperatureCoordinates )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get() );

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::temperatureCoordinatesString()));

  auto & temperatureCoordinates = m_fluid->getWrapper< array1d< real64 > >( Keys::temperatureCoordinatesString() ).reference();

  // Can be empty
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  // Should have at least 2 elements
  temperatureCoordinates.resize( 1 );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );

  // Should be increasing
  temperatureCoordinates.resize( 3 );
  TestFluid< 3 >::populateArray( temperatureCoordinates, Feed< 3 >{ 303.15, 313.15, 298.15 } );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );

  // Should be strictly increasing
  temperatureCoordinates.resize( 3 );
  TestFluid< 3 >::populateArray( temperatureCoordinates, Feed< 3 >{ 303.15, 313.15, 313.15 } );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
}

} // namespace testing

} //namespace geos
