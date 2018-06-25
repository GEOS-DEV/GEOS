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

#include "dataRepository/ReferenceWrapper.hpp"
#include "ManagedArray.hpp"

#include <functional>
#include <string>
#include <typeindex>
#include <vector>

using namespace geosx;
using namespace multidimensionalArray;

TEST(testReferenceWrapper,testIntWrapper)
{

  int var = 0;

  ReferenceWrapper<int> wrappedVar( var );

  wrappedVar = 5;

  int var2 = wrappedVar;

  EXPECT_TRUE( var == 5 );
  EXPECT_TRUE( var2 == var );

  EXPECT_TRUE( &var == &(wrappedVar.get()) );
  EXPECT_TRUE( &var2 != &(wrappedVar.get()) );
  EXPECT_TRUE( &var != &var2 );

}

TEST(testReferenceWrapper,testArrayWrapper)
{
  ManagedArray< int, 1 , int > arr;
  arr.resize( 4 );

  ReferenceWrapper< ManagedArray< int, 1 , int > > wrappedArr( arr );

  for( int i=0 ; i<4 ; ++i )
  {
    wrappedArr[i] = 2 * i;
  }

  EXPECT_TRUE( arr[0] == 0 );
  EXPECT_TRUE( arr[1] == 2 );
  EXPECT_TRUE( arr[2] == 4 );
  EXPECT_TRUE( arr[3] == 6 );

}

TEST(testReferenceWrapper,testArrayOfWrappedInts)
{
  int val0 = 0;
  int val1 = 1;

  std::vector< ReferenceWrapper<int> > arr2;
  arr2.resize(2);
  arr2[0].set( val0 );
  arr2[1].set( val1 );

  arr2[0] = 10;
  arr2[1] = 11;

  EXPECT_TRUE( val0 == 10 );
  EXPECT_TRUE( val1 == 11 );
}

TEST(testReferenceWrapper,testOperatorParen)
{
  using array2d = ManagedArray< int, 2 , int >;
  array2d arr;
  arr.resize( 2,3 );

  ReferenceWrapper< array2d > wrappedArr( arr );

  for( int i=0 ; i<2 ; ++i )
  {
    for( int j=0 ; j<3 ; ++j )
    {
      wrappedArr(i,j) = 3*i+j;
    }
  }

  EXPECT_TRUE( arr[0][0] == 0 );
  EXPECT_TRUE( arr[0][1] == 1 );
  EXPECT_TRUE( arr[0][2] == 2 );
  EXPECT_TRUE( arr[1][0] == 3 );
  EXPECT_TRUE( arr[1][1] == 4 );
  EXPECT_TRUE( arr[1][2] == 5 );
}

TEST(testReferenceWrapper,testNestedOperatorSquare)
{
  using array2d = ManagedArray< int, 2 , int >;
  array2d arr;
  arr.resize( 2,3 );

  ReferenceWrapper< array2d > wrappedArr( arr );

  for( int i=0 ; i<2 ; ++i )
  {
    for( int j=0 ; j<3 ; ++j )
    {
      wrappedArr[i][j] = 3*i+j;
    }
  }

  EXPECT_TRUE( arr[0][0] == 0 );
  EXPECT_TRUE( arr[0][1] == 1 );
  EXPECT_TRUE( arr[0][2] == 2 );
  EXPECT_TRUE( arr[1][0] == 3 );
  EXPECT_TRUE( arr[1][1] == 4 );
  EXPECT_TRUE( arr[1][2] == 5 );
}
