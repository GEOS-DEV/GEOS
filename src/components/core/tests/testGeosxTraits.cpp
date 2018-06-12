// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
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

#include "../src/codingUtilities/GeosxTraits.hpp"
#include<type_traits>

using namespace geosx;
using namespace geosx::traits;
using namespace bufferOps;

TEST(testGeosxTraits,test_is_tensorT)
{
  EXPECT_TRUE( is_tensorT<R1Tensor>::value );
  EXPECT_TRUE( is_tensorT<R2Tensor>::value );
  EXPECT_TRUE( is_tensorT<R2SymTensor>::value );

  EXPECT_FALSE( is_tensorT<int>::value );
  EXPECT_FALSE( is_tensorT<double>::value );
  EXPECT_FALSE( is_tensorT<void>::value );
}


TEST(testGeosxTraits,test_is_string)
{
  EXPECT_TRUE( is_string<string>::value );

  EXPECT_FALSE( is_string<int>::value );
  EXPECT_FALSE( is_string<double>::value );
  EXPECT_FALSE( is_string<void>::value );
}

TEST(testGeosxTraits,test_is_array)
{
  using arrayType = multidimensionalArray::ManagedArray<int,1,int>;
  EXPECT_TRUE( is_array<arrayType >::value );

  EXPECT_FALSE( is_array<int>::value );
  EXPECT_FALSE( is_array<double>::value );
  EXPECT_FALSE( is_array<void>::value );
}

TEST(testGeosxTraits,test_is_map)
{
  using mapType = map<string,int>;
  EXPECT_TRUE( is_map< mapType >::value );

  EXPECT_FALSE( is_map<int>::value );
  EXPECT_FALSE( is_map<double>::value );
  EXPECT_FALSE( is_map<void>::value );
}

TEST(testGeosxTraits,test_is_pair)
{
  using pairType = std::pair<string,int>;
  EXPECT_TRUE( is_pair< pairType >::value );

  EXPECT_FALSE( is_pair<int>::value );
  EXPECT_FALSE( is_pair<double>::value );
  EXPECT_FALSE( is_pair<void>::value );

}


TEST(testGeosxTraits,test_is_noncontainer_type_packable)
{

  EXPECT_TRUE( is_noncontainer_type_packable<int>::value );
  EXPECT_TRUE( is_noncontainer_type_packable<double>::value );
  EXPECT_FALSE( is_noncontainer_type_packable<void>::value );
  EXPECT_TRUE( is_noncontainer_type_packable<R1Tensor>::value );
  EXPECT_FALSE( is_noncontainer_type_packable< array<double> >::value );
  EXPECT_FALSE( is_noncontainer_type_packable< set<double> >::value );
  using mapType = map<string,int>;
  EXPECT_FALSE( is_noncontainer_type_packable< mapType >::value );
  using pairType = std::pair<string,int>;
  EXPECT_FALSE( is_noncontainer_type_packable< pairType >::value );
}


TEST(testGeosxTraits,test_is_array_packable)
{
  using arrayType = multidimensionalArray::ManagedArray<double,2,long int>;
  EXPECT_TRUE( is_packable_array<arrayType>::value );

  using arrayType0 = multidimensionalArray::ManagedArray<void,1,int>;
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

  using arrayType = multidimensionalArray::ManagedArray<int,1,int>;
  using mapType2 = map<string,arrayType>;
  EXPECT_TRUE( is_packable_map< mapType2 >::value );

}
