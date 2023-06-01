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

#include "common/TypeDispatch.hpp"

#include <gtest/gtest.h>

using namespace geos;

namespace typeDispatchTest
{
class A
{
public:
  A() = default;
  virtual ~A() = default;
  virtual string getTypeName() const {return "A";}
};

class B : public A
{
public:
  virtual string getTypeName() const { return "B";}
};

class C : public A
{
public:
  virtual string getTypeName() const { return "C";}
};

class D : public A
{
public:
  virtual string getTypeName() const { return "D";}
};

}

template< typename ... Ts >
void testDispatchSimple( types::TypeList< Ts... > const list )
{
  bool const result = types::dispatchCombinations( list, []( auto tupleOfTypes )
  {
    EXPECT_TRUE( ( std::is_same< camp::first< decltype(tupleOfTypes) >, float >::value ) ) << "Dispatch matched the wrong type";
    EXPECT_TRUE( ( std::is_same< camp::second< decltype(tupleOfTypes) >, int >::value ) ) << "Dispatch matched the wrong type";
  }, 0.0f, 0 );
  EXPECT_TRUE( result ) << "Dispatch failed to match the type";
}

template< typename TRUE_TYPES_LIST, typename ... Ts >
void testDispatchVirtual( types::TypeList< Ts... > const list,
                          typeDispatchTest::A & t1,
                          typeDispatchTest::A & t2 )
{
  bool const result = types::dispatchCombinations( list, []( auto tupleOfTypes )
  {
    EXPECT_TRUE( ( std::is_same< camp::first< decltype(tupleOfTypes) >, camp::first< TRUE_TYPES_LIST > >::value ) ) << "Dispatch matched the wrong type";
    EXPECT_TRUE( ( std::is_same< camp::second< decltype(tupleOfTypes) >, camp::second< TRUE_TYPES_LIST > >::value ) ) << "Dispatch matched the wrong type";
  }, t1, t2 );
  EXPECT_TRUE( result ) << "Dispatch failed to match the type";
}


TEST( testDispatchSimple, DispatchSimpleTypes )
{
  using Types = types::TypeList< types::TypeList< int, int >,
                                 types::TypeList< float, int > >;

  testDispatchSimple( Types{} );
}

TEST( testDispatchVirtual, DispatchVirtualTypes )
{

  using namespace typeDispatchTest;
  using Types = types::TypeList< types::TypeList< B, B >,
                                 types::TypeList< B, C >,
                                 types::TypeList< B, D > >;

  std::unique_ptr< A > b = std::make_unique< B >();
  std::unique_ptr< A > c = std::make_unique< C >();

  // std::cout << "I am " << b->getTypeName() << std::endl;
  // std::cout << "I am " << c->getTypeName() << std::endl;

  testDispatchVirtual< types::TypeList< B, C > >( Types{}, *b, *c );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
