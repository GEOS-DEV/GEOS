/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "dataRepository/DynamicCasts.hpp"
#include "dataRepository/Group.hpp"

#include <gtest/gtest.h>

namespace geosx
{
namespace dataRepository
{
namespace tests
{

class Dog : private GroupDownCastHelper< Dog >
{
  using GroupDownCastHelper::GroupDownCastHelper;
  virtual ~Dog(){};
};

class Cat : public Group
{
  using Group::Group;
  virtual ~Cat(){};
};

TEST( testGeosxTraits, test_dynamic_casts_behavior )
{
  using details::DynamicCast;

  Dog * dog = new Dog( "dog", nullptr );
  Cat * cat = new Cat( "cat", nullptr );

  Group * dogAsGroup = reinterpret_cast< Group * >( dog );
  Group const * dogAsConstGroup = reinterpret_cast< Group const * >( dog );

  Group * catAsGroup = cat;
  Group const * catAsConstGroup = cat;

  // dogs to dogs
  EXPECT_EQ( dog, DynamicCast< Dog * >( dogAsGroup ) );
  EXPECT_EQ( dog, DynamicCast< Dog const * >( dogAsGroup ) );
  EXPECT_EQ( dog, DynamicCast< Dog const * >( dogAsConstGroup ) );

  // cats to cats
  EXPECT_EQ( cat, DynamicCast< Cat * >( catAsGroup ) );
  EXPECT_EQ( cat, DynamicCast< Cat const * >( catAsGroup ) );
  EXPECT_EQ( cat, DynamicCast< Cat const * >( catAsConstGroup ) );

  // interweave
  EXPECT_EQ( nullptr, DynamicCast< Cat * >( dogAsGroup ) );
  EXPECT_EQ( nullptr, DynamicCast< Dog * >( catAsGroup ) );
}

} // end of namespace tests
} // end of namespace dataRepository
} // end of namespace geosx
