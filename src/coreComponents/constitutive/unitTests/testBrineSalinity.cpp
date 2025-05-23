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
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"
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

class BrineSalinityTestFixture : public ::testing::Test
{
public:
  using Keys = BrineSalinity::viewKeyStruct;
  static constexpr integer numComps = 4;
  static constexpr real64 absTol = 1.0e-10;

public:
  BrineSalinityTestFixture();
  ~BrineSalinityTestFixture() override = default;

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  std::unique_ptr< MockFluid > m_fluid{};
  std::unique_ptr< TestFluid< numComps > > m_testFluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
  BrineSalinity * m_salinity{};
};

BrineSalinityTestFixture::BrineSalinityTestFixture():
  m_parent( "parent", m_node ),
  m_fluid( std::make_unique< MockFluid >( "fluid", &m_parent )),
  m_testFluid( TestFluid< numComps >::create( {Fluid::CO2, Fluid::H2O, Fluid::C1, Fluid::C5} ))
{
  m_parameters = BrineSalinity::create( std::move( m_parameters ));
  m_salinity = dynamic_cast< BrineSalinity * >(m_parameters.get());
  GEOS_ERROR_IF( m_salinity == nullptr, "Cannot create BrineSalinity" );
}

TEST_F( BrineSalinityTestFixture, testSalinity )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get());

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::salinityString()));

  real64 & salinity = m_fluid->getWrapper< real64 >( Keys::salinityString()).reference();
  salinity = 0.45;
  EXPECT_NEAR( m_salinity->m_salinity, salinity, absTol );

  salinity = -1.0;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
}

TEST_F( BrineSalinityTestFixture, testWaterCompressibility )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get());

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::waterCompressibilityString()));

  real64 & waterCompressibility = m_fluid->getWrapper< real64 >( Keys::waterCompressibilityString()).reference();
  waterCompressibility = 2.0e-10;
  EXPECT_NEAR( m_salinity->m_waterCompressibility, waterCompressibility, absTol );

  waterCompressibility = 0.0;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );

  waterCompressibility = -1.0e-10;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
}

TEST_F( BrineSalinityTestFixture, testSaltMolarWeight )
{
  auto componentProperties = m_testFluid->getComponentProperties();
  m_parameters->registerParameters( m_fluid.get());

  EXPECT_TRUE( m_fluid->hasWrapper( Keys::saltMolarWeightString()));

  real64 & saltMolarWeight = m_fluid->getWrapper< real64 >( Keys::saltMolarWeightString()).reference();
  saltMolarWeight = 43.22;
  EXPECT_NEAR( m_salinity->m_saltMolarWeight, saltMolarWeight, absTol );

  saltMolarWeight = 0.0;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );

  saltMolarWeight = -18.22;
  EXPECT_THROW( m_parameters->postInputInitialization( m_fluid.get(), componentProperties ), InputError );
}

} // namespace testing

} //namespace geos
