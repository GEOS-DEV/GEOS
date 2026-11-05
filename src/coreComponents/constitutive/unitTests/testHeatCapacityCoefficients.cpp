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
#include "constitutive/fluid/multifluid/compositional/parameters/HeatCapacityCoefficients.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
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
  MockFluid( string const & name, Group * const parent ):
    MultiFluidBase( name, parent )
  {}
  ~MockFluid() override = default;

  string getCatalogName() const override { return "FluidModel"; }
  integer getWaterPhaseIndex() const override { return 1; }
  void checkTablesParameters( real64, real64 ) const override {}
};

class HeatCapacityCoefficientsTestFixture : public ::testing::Test
{
public:
  using Keys = HeatCapacityCoefficients::viewKeyStruct;
  static constexpr integer numComps = 4;
  static constexpr integer numCoeffs = 5;
  static constexpr real64 absTol = 1.0e-10;

public:
  HeatCapacityCoefficientsTestFixture();
  ~HeatCapacityCoefficientsTestFixture() override = default;

protected:
  void createFluid();
protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  std::unique_ptr< MockFluid > m_fluid{};
  std::unique_ptr< TestFluid< numComps > > m_testFluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
  HeatCapacityCoefficients * m_heatCapacityCoefficients{};
};

HeatCapacityCoefficientsTestFixture::HeatCapacityCoefficientsTestFixture():
  m_parent( "parent", m_node ),
  m_fluid( std::make_unique< MockFluid >( "fluid", &m_parent )),
  m_testFluid( TestFluid< numComps >::create( {Fluid::CO2, Fluid::H2O, Fluid::CH4, Fluid::C5H12} ))
{
  m_parameters = HeatCapacityCoefficients::create( std::move( m_parameters ));
  m_heatCapacityCoefficients = dynamic_cast< HeatCapacityCoefficients * >(m_parameters.get());
  GEOS_ERROR_IF( m_heatCapacityCoefficients == nullptr, "Cannot create HeatCapacityCoefficients" );
}

void HeatCapacityCoefficientsTestFixture::createFluid()
{
  using FluidKeys = MultiFluidBase::viewKeyStruct;
  auto & phaseNames = m_fluid->getWrapper< string_array >( FluidKeys::phaseNamesString() ).reference();
  phaseNames.emplace_back( "gas" );
  phaseNames.emplace_back( "water" );
  phaseNames.emplace_back( "oil" );
  integer const numPhases = static_cast< integer >(phaseNames.size());

  auto & componentNames = m_fluid->getWrapper< string_array >( FluidKeys::componentNamesString() ).reference();
  for( integer ic = 0; ic < numComps; ++ic )
  {
    componentNames.emplace_back( m_testFluid->componentNames[ic] );
  }

  auto & equationsOfState = m_fluid->getWrapper< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() ).reference();
  string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    equationsOfState.emplace_back( eosName );
  }

  m_heatCapacityCoefficients->m_referenceEnthalpy.resize( numPhases, numComps );
  m_heatCapacityCoefficients->m_referenceEnthalpy.zero();
  m_heatCapacityCoefficients->m_coefficients.resize( numComps, numCoeffs );
  m_heatCapacityCoefficients->m_coefficients.zero();
  for( integer ic = 0; ic < numComps; ic++ )
  {
    m_heatCapacityCoefficients->m_coefficients( ic, 0 ) = 1.0;
  }
}

TEST_F( HeatCapacityCoefficientsTestFixture, testReferenceTemperature )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get());
  createFluid();

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::enthalpyReferenceTemperatureString()) );

  real64 & referenceTemperature = m_fluid->getWrapper< real64 >( Keys::enthalpyReferenceTemperatureString() ).reference();
  referenceTemperature = 288.706;
  EXPECT_NEAR( m_heatCapacityCoefficients->m_referenceTemperature, referenceTemperature, absTol );

  referenceTemperature = -1.0;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
}

TEST_F( HeatCapacityCoefficientsTestFixture, testReferenceEnthalpy )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get());
  createFluid();

  integer const numPhases = m_fluid->numFluidPhases();

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::referenceEnthalpyString()) );

  array2d< real64 > & referenceEnthalpy = m_fluid->getWrapper< array2d< real64 > >( Keys::referenceEnthalpyString() ).reference();
  referenceEnthalpy.resize( numPhases - 1, numComps );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  referenceEnthalpy.resize( numPhases + 1, numComps );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  referenceEnthalpy.resize( numPhases, numComps - 1 );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  referenceEnthalpy.resize( numPhases, numComps + 1 );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  referenceEnthalpy.resize( numPhases, numComps );
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  constexpr integer vapour = 0;
  constexpr integer liquid = 2;
  constexpr integer aqueous = 1;

  referenceEnthalpy( vapour, 2 ) = 100.0;
  referenceEnthalpy( liquid, 2 ) = 100.0;
  referenceEnthalpy( aqueous, 2 ) = 100.0;

  // Reference enthalpy of gas must be greater than oil
  referenceEnthalpy( liquid, 2 ) = 110.0;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  referenceEnthalpy( liquid, 2 ) = 100.0;
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  // Reference enthalpy of gas must be greater than water
  referenceEnthalpy( aqueous, 2 ) = 120.0;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  referenceEnthalpy( aqueous, 2 ) = 100.0;
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  // Reference enthalpy of water must be greater than oil
  referenceEnthalpy( aqueous, 2 ) = 90.0;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  referenceEnthalpy( aqueous, 2 ) = 100.0;
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );
}

TEST_F( HeatCapacityCoefficientsTestFixture, testHeatCapacityCoefficients )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get());
  createFluid();

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::componentHeatCapacityCoefficientsString()) );

  array2d< real64 > & coefficients = m_fluid->getWrapper< array2d< real64 > >( Keys::componentHeatCapacityCoefficientsString() ).reference();

  auto const setValues = [&]( std::array< real64 const, numCoeffs > const & data )
  {
    for( integer ic = 0; ic < numCoeffs; ++ic )
    {
      coefficients( 0, ic ) = data[ic];
    }
  };

  setValues( {0.0, 0.0, 0.0, 0.0, 0.0} );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  setValues( {1.0, 0.0, 0.0, 0.0, 0.0} );
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  setValues( {9.7529018563e+00, -1.7379753331e-01, 1.3828200497e-03, -4.4470070713e-06, 4.3687436501e-09} );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  setValues( {1.0, 0.0, 0.0, 0.0, 0.0} );
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  setValues( {8.440483444351720e-02, -8.156515860472481e-04, 2.921178619177309e-06, -4.649781358999949e-09, 2.775470116151693e-12} );
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
  setValues( {1.0, 0.0, 0.0, 0.0, 0.0} );
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );

  setValues( {1.469119322002621e-02, -1.822650742125867e-04, 8.537827827668615e-07, -1.777519183416209e-09, 1.387740262103405e-12} );
  EXPECT_NO_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ) );
}

} // namespace testing

} //namespace geos
