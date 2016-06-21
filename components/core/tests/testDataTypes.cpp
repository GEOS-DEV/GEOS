/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"
#include "DataTypes.hpp"
#include <unordered_map>
#include <string>
#include <typeindex>
// API coverage tests
// Each test should be documented with the interface functions being tested

  template< typename T >
  void Func1()
  {

    std::cout << "called " << rtTypes::typeNames(std::type_index(typeid(T)))<<" version of Func1" << std::endl;
//    std::cout<<"called "<<typeid(T).name()<<std::endl;
  }

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(testDataTypes,applyTypeLambda)
{

  rtTypes::ApplyTypeLambda( rtTypes::TypeIDs::int32_id, []( auto a ) -> void { return Func1<decltype(a)>(); });
  rtTypes::ApplyTypeLambda( rtTypes::TypeIDs::int64_id, []( auto a ) -> void { return Func1<decltype(a)>(); });
  rtTypes::ApplyTypeLambda( rtTypes::TypeIDs::uint32_id, []( auto a ) -> void { return Func1<decltype(a)>(); });
  rtTypes::ApplyTypeLambda( rtTypes::TypeIDs::uint64_id, []( auto a ) -> void { return Func1<decltype(a)>(); });
  rtTypes::ApplyTypeLambda( rtTypes::TypeIDs::real32_id, []( auto a ) -> void { return Func1<decltype(a)>(); });
  rtTypes::ApplyTypeLambda( rtTypes::TypeIDs::real64_id, []( auto a ) -> void { return Func1<decltype(a)>(); });

}



