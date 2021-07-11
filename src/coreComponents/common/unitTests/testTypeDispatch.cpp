/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "common/TypeDispatch.hpp"

#include <gtest/gtest.h>

using namespace geosx;

template< typename TYPES >
class TypeDispatchTest : public ::testing::Test
{
protected:

  using TypeList = TYPES;

  template< typename ... Ts >
  void matchAllTypes( types::TypeList< Ts... >,
                      types::TypeList<> )
  {}

  template< typename ... Ts, typename T, typename ... Rest >
  void matchAllTypes( types::TypeList< Ts... > const list,
                      types::TypeList< T, Rest... > )
  {
    // Test both that dispatch occurred and that type given to lambda matches expected
    bool const success = types::dispatch( list, typeid( T ), false, []( auto v )
    {
      EXPECT_TRUE( ( std::is_same< decltype( v ), T >::value ) ) << "Dispatch matched the wrong type";
    } );
    EXPECT_TRUE( success ) << "Dispatch failed to match the type";
    matchAllTypes( list, types::TypeList< Rest... >{} );
  }

  template< typename ... Ts >
  void checkTypeNotInList( types::TypeList< Ts... > const list )
  {
    struct TypeNotInList {};
    bool const success = types::dispatch( list, typeid( TypeNotInList ), false, []( auto ){} );
    EXPECT_FALSE( success ) << "Dispatch matched a type not in list";
  }

};

struct X {};

using TypeLists = ::testing::Types<
  types::TypeList<>,
  types::TypeList< int, double, X >,
  types::ArrayTypes< types::TypeList< int, double, X >, types::DimsUpTo< 3 > >
  >;

TYPED_TEST_SUITE( TypeDispatchTest, TypeLists, );

TYPED_TEST( TypeDispatchTest, MatchAllTypes )
{
  using Types = typename TestFixture::TypeList;
  this->matchAllTypes( Types{}, Types{} );
}

TYPED_TEST( TypeDispatchTest, CheckTypeNotInList )
{
  using Types = typename TestFixture::TypeList;
  this->checkTypeNotInList( Types{} );
}
