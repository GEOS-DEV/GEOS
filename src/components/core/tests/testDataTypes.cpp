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
#endif

#include "common/DataTypes.hpp"
#include <string>
#include <typeindex>


template< typename T >
std::string Func1()
{
  return rtTypes::typeNames(std::type_index(typeid(T)));
}

TEST(testDataTypes,applyTypeLambda)
{
  std::string funcReturn_int32  = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::int32_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });
  std::string funcReturn_int64  = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::int64_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });
  std::string funcReturn_uint32 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::uint32_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });
  std::string funcReturn_uint64 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::uint64_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });
  std::string funcReturn_real32 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::real32_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });
  std::string funcReturn_real64 = rtTypes::ApplyTypeLambda1( rtTypes::TypeIDs::real64_id, []( auto a ) -> std::string {
    return Func1<decltype(a)>();
  });

  EXPECT_TRUE( funcReturn_int32.compare(rtTypes::typeNames( std::type_index( typeid(int32) ) ) ) == 0);
  EXPECT_TRUE( funcReturn_int64.compare(rtTypes::typeNames( std::type_index( typeid(int64) ) ) ) == 0);
  EXPECT_TRUE( funcReturn_uint32.compare(rtTypes::typeNames( std::type_index( typeid(uint32) ) ) ) == 0);
  EXPECT_TRUE( funcReturn_uint64.compare(rtTypes::typeNames( std::type_index( typeid(uint64) ) ) ) == 0);
  EXPECT_TRUE( funcReturn_real32.compare(rtTypes::typeNames( std::type_index( typeid(real32) ) ) ) == 0);
  EXPECT_TRUE( funcReturn_real64.compare(rtTypes::typeNames( std::type_index( typeid(real64) ) ) ) == 0);

}
