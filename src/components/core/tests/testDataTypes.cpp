/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include <gtest/gtest.h>

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "common/DataTypes.hpp"
#include <string>
#include <typeindex>

using namespace geosx;

template< typename T >
std::string Func1()
{
  return rtTypes::typeNames(std::type_index(typeid(T)));
}

TEST(testDataTypes,applyTypeLambda)
{
  std::string funcReturn_integer  = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::integer_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });
  std::string funcReturn_real32 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::real32_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });
  std::string funcReturn_real64 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::real64_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });

  EXPECT_TRUE( funcReturn_integer.compare(rtTypes::typeNames( std::type_index( typeid(integer) ) ) ) == 0);
  EXPECT_TRUE( funcReturn_real32.compare(rtTypes::typeNames( std::type_index( typeid(real32) ) ) ) == 0);
  EXPECT_TRUE( funcReturn_real64.compare(rtTypes::typeNames( std::type_index( typeid(real64) ) ) ) == 0);

}
