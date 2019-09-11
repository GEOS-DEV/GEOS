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

#include <gtest/gtest.h>

#include "common/DataTypes.hpp"
#include <string>
#include <typeindex>

using namespace geosx;

template< typename T >
std::string Func1()
{
  return rtTypes::typeNames( std::type_index( typeid(T)));
}

TEST( testDataTypes, applyTypeLambda )
{
  std::string funcReturn_integer  = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::integer_id, []( auto a ) -> std::string {
    return Func1< decltype(a) >();
  } );
  std::string funcReturn_real32 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::real32_id, []( auto a ) -> std::string {
    return Func1< decltype(a) >();
  } );
  std::string funcReturn_real64 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::real64_id, []( auto a ) -> std::string {
    return Func1< decltype(a) >();
  } );

  EXPECT_TRUE( funcReturn_integer.compare( rtTypes::typeNames( std::type_index( typeid(integer) ) ) ) == 0 );
  EXPECT_TRUE( funcReturn_real32.compare( rtTypes::typeNames( std::type_index( typeid(real32) ) ) ) == 0 );
  EXPECT_TRUE( funcReturn_real64.compare( rtTypes::typeNames( std::type_index( typeid(real64) ) ) ) == 0 );

}
