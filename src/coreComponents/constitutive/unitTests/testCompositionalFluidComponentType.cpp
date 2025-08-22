/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentType.hpp"

// TPL includes
#include <gtest/gtest.h>

namespace geos
{

using namespace constitutive;
using namespace constitutive::compositional;

namespace testing
{

TEST( ComponentTypeTest, testValidComponents )
{
  ComponentType const water = ComponentType::Water;
  integer const comp1 = static_cast< integer >(water);
  EXPECT_TRUE( isComponentType( comp1, water ));

  ComponentType const hydrogen = ComponentType::Hydrogen;
  integer const comp2 = static_cast< integer >(hydrogen);
  EXPECT_TRUE( isComponentType( comp2, hydrogen ));

  ComponentType const hydrogenSulphide = ComponentType::HydrogenSulphide;
  integer const comp3 = static_cast< integer >(hydrogenSulphide);
  EXPECT_TRUE( isComponentType( comp3, hydrogenSulphide ));

  ComponentType const carbonDioxide = ComponentType::CarbonDioxide;
  integer const comp4 = static_cast< integer >(carbonDioxide);
  EXPECT_TRUE( isComponentType( comp4, carbonDioxide ));

  ComponentType const nitrogen = ComponentType::Nitrogen;
  integer const comp5 = static_cast< integer >(nitrogen);
  EXPECT_TRUE( isComponentType( comp5, nitrogen ));
}

TEST( ComponentTypeTest, testComponentsNames )
{
  EXPECT_EQ( getComponentTypeFromName( "c1" ), ComponentType::HydroCarbon );
  EXPECT_EQ( getComponentTypeFromName( "ch4" ), ComponentType::HydroCarbon );
  EXPECT_EQ( getComponentTypeFromName( "CH4" ), ComponentType::HydroCarbon );

  EXPECT_EQ( getComponentTypeFromName( "h2" ), ComponentType::Hydrogen );
  EXPECT_EQ( getComponentTypeFromName( "H2" ), ComponentType::Hydrogen );

  EXPECT_EQ( getComponentTypeFromName( "h2o" ), ComponentType::Water );
  EXPECT_EQ( getComponentTypeFromName( "H2O" ), ComponentType::Water );

  EXPECT_EQ( getComponentTypeFromName( "CO2" ), ComponentType::CarbonDioxide );
  EXPECT_EQ( getComponentTypeFromName( "co2" ), ComponentType::CarbonDioxide );

  EXPECT_EQ( getComponentTypeFromName( "h2s" ), ComponentType::HydrogenSulphide );
  EXPECT_EQ( getComponentTypeFromName( "H2S" ), ComponentType::HydrogenSulphide );

  EXPECT_EQ( getComponentTypeFromName( "N2" ), ComponentType::Nitrogen );
  EXPECT_EQ( getComponentTypeFromName( "n2" ), ComponentType::Nitrogen );
}

} // namespace testing

} // namespace geos
