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

#include "Array.hpp"
#include "dataRepository/DefaultValue.hpp"

#include <functional>
#include <string>
#include <typeindex>
#include <vector>

using namespace geosx;
using namespace LvArray;
using namespace dataRepository;


TEST(testDefaultValue,testScalar)
{
  EXPECT_TRUE( DefaultValue<int>::has_default_value == true );
  EXPECT_TRUE( DefaultValue<long int>::has_default_value == true );
  EXPECT_TRUE( DefaultValue<long long int>::has_default_value == true );
  EXPECT_TRUE( DefaultValue<double>::has_default_value == true );
}

TEST(testDefaultValue,testArray)
{
  using array1 = Array<double,1,int>;
  using array2 = Array<double,2,int>;
  using array3 = Array<double,3,int>;
  using array4 = Array<int,1,int>;
  using array5 = Array<long int,1,long int>;
  using array6 = Array<long long int,1,long long int>;
  EXPECT_TRUE( DefaultValue<array1>::has_default_value==true );
  EXPECT_TRUE( DefaultValue<array2>::has_default_value==true );
  EXPECT_TRUE( DefaultValue<array3>::has_default_value==true );
  EXPECT_TRUE( DefaultValue<array4>::has_default_value==true );
  EXPECT_TRUE( DefaultValue<array5>::has_default_value==true );
  EXPECT_TRUE( DefaultValue<array6>::has_default_value==true );
}

