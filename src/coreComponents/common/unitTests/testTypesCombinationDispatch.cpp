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

template < typename TRUE_TYPES_LIST, size_t I = 0, typename... Ts >
typename std::enable_if<I == sizeof...(Ts), void>::type
checkType(types::TypeList<Ts...> )
{
    return;
}
 
template< typename TRUE_TYPES_LIST, size_t I = 0, typename... Ts>
typename std::enable_if< (I < sizeof...(Ts)), void>::type
checkType( types::TypeList<Ts...> list )
{
    EXPECT_TRUE( ( std::is_same< camp::at_t< decltype(list), camp::num<I> >, 
                                 camp::at_t< TRUE_TYPES_LIST, camp::num<I> > >::value ) ) << "Dispatch matched the wrong type";
    // check next element
    checkType< TRUE_TYPES_LIST, I + 1>( list );
}

template< typename TRUE_TYPES_LIST, typename ... Ts, typename ...Vs >
void testDispatch( types::TypeList< Ts... > const list,
                          Vs && ... objects )
{
  bool const result = types::dispatchCombinations( list, []( auto tupleOfTypes )
  { 
    checkType< TRUE_TYPES_LIST >( tupleOfTypes );
  }, std::forward< Vs >( objects ) ... );
  EXPECT_TRUE( result ) << "Dispatch failed to match the type";
}

TEST( testDispatch, DispatchSimpleTypesSingles )
{
  using Types = types::ListofTypeList< types::TypeList< int, float > >;

  int a = 0;                               
  testDispatch< types::TypeList< int > >( Types{}, a );                             
}


TEST( testDispatch, DispatchSimpleTypes )
{
  using Types = types::TypeList< types::TypeList< int, int >,
                                 types::TypeList< int, float > >;

  int a = 0;
  float b = 0.0;
  
  testDispatch< types::TypeList< int, float > >( Types{}, a, b );   

}

TEST( testDispatch, DispatchVirtualTypePairs )
{

  using namespace typeDispatchTest;
  using Types = types::TypeList< types::TypeList< B, B >,
                                 types::TypeList< B, C >,
                                 types::TypeList< B, D > >;

  std::unique_ptr< A > b = std::make_unique< B >();
  std::unique_ptr< A > c = std::make_unique< C >();

  testDispatch< types::TypeList< B, C > >( Types{}, *b, *c );  
}

TEST( testDispatch, DispatchVirtualTypeTriplets )
{

  using namespace typeDispatchTest;
  using Types = types::TypeList< types::TypeList< B, B, B>,
                                 types::TypeList< B, C, C >,
                                 types::TypeList< B, C, D > >;


  std::cout << typeid(types::TypeList< C, B, B>).name() << std::endl; 
  std::cout << typeid(types::TypeList< B, B, B>).name() << std::endl;                                                     

  std::unique_ptr< A > b = std::make_unique< B >();
  std::unique_ptr< A > c = std::make_unique< C >();
  std::unique_ptr< A > d = std::make_unique< D >();

  testDispatch< types::TypeList< B, C, D > >( Types{}, *b, *c, *d);
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
