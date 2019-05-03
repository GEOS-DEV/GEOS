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

#include "codingUtilities/GeosxTraits.hpp"

using namespace geosx;
using namespace geosx::traits;
using namespace bufferOps;

TEST(testGeosxTraits,test_is_tensorT)
{
  static_assert(is_tensorT<R1Tensor>, "Should be true");
  static_assert(is_tensorT<R2Tensor>, "Should be true");
  static_assert(is_tensorT<R2SymTensor>, "Should be true");

  static_assert(!is_tensorT<int>, "Should be false");
  static_assert(!is_tensorT<double>, "Should be false");
  static_assert(!is_tensorT<void>, "Should be false");
}


TEST(testGeosxTraits,test_is_string)
{
  static_assert(is_string<string>, "Should be true");

  static_assert(!is_string<int>, "Should be false");
  static_assert(!is_string<double>, "Should be false");
  static_assert(!is_string<void>, "Should be false");
}

TEST(testGeosxTraits,test_is_array)
{
  static_assert(is_array<LvArray::Array<int,1,int>>, "Should be true");
  
  static_assert(!is_array<int>, "Should be false");
  static_assert(!is_array<double>, "Should be false");
  static_assert(!is_array<void>, "Should be false");
}

TEST(testGeosxTraits,test_is_map)
{
  static_assert(is_map<map<string,int>>, "Should be true");
  static_assert(is_map<unordered_map<string,int>>, "Should be true");

  static_assert(!is_map<int>, "Should be false");
  static_assert(!is_map<double>, "Should be false");
  static_assert(!is_map<void>, "Should be false");
  SUCCEED();
}

TEST(testGeosxTraits,test_is_set)
{
  static_assert(is_set<set<string>>, "Should be true");

  static_assert(!is_set<int>, "Should be false");
  static_assert(!is_set<double>, "Should be false");
  static_assert(!is_set<void>, "Should be false");
  SUCCEED();
}

TEST(testGeosxTraits,test_is_pair)
{
  static_assert(is_pair<std::pair<string,int>>, "Should be true");

  static_assert(!is_pair<int>, "Should be false");
  static_assert(!is_pair<double>, "Should be false");
  static_assert(!is_pair<void>, "Should be false");
  SUCCEED();
}

TEST(testGeosxTraits,test_is_noncontainer_type_packable)
{
  EXPECT_TRUE( is_noncontainer_type_packable<int>::value );
  EXPECT_TRUE( is_noncontainer_type_packable<double>::value );
  EXPECT_FALSE( is_noncontainer_type_packable<void>::value );
  EXPECT_TRUE( is_noncontainer_type_packable<R1Tensor>::value );
  EXPECT_FALSE( is_noncontainer_type_packable< array1d<double> >::value );
  EXPECT_FALSE( is_noncontainer_type_packable< set<double> >::value );
  using mapType = map<string,int>;
  EXPECT_FALSE( is_noncontainer_type_packable< mapType >::value );
  using pairType = std::pair<string,int>;
  EXPECT_FALSE( is_noncontainer_type_packable< pairType >::value );
}


TEST(testGeosxTraits,test_is_array_packable)
{
  using arrayType = LvArray::Array<double,2,long int>;
  EXPECT_TRUE( is_packable_array<arrayType>::value );

  using arrayType0 = LvArray::Array<void,1,int>;
  EXPECT_FALSE( is_packable_array<arrayType0>::value );

  EXPECT_FALSE( is_packable_array<int>::value );
  EXPECT_FALSE( is_packable_array<double>::value );
  EXPECT_FALSE( is_packable_array<void>::value );
}


TEST(testGeosxTraits,test_is_packable_map)
{
  using mapType0 = map<string,int>;
  EXPECT_TRUE( is_packable_map< mapType0 >::value );

  using mapType1 = map<string,std::pair<int,int> >;
  EXPECT_FALSE( is_packable_map< mapType1 >::value );

  using arrayType = LvArray::Array<int,1,int>;
  using mapType2 = map<string,arrayType>;
  EXPECT_TRUE( is_packable_map< mapType2 >::value );

}
