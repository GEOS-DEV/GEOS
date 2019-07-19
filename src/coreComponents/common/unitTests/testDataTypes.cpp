/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
