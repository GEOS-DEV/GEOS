/*
* ------------------------------------------------------------------------------------------------------------
* SPDX-License-Identifier: LGPL-2.1-only
*
* Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
* Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
* Copyright (c) 2018-2020 TotalEnergies
* Copyright (c) 2019-     GEOSX Contributors
* All rights reserved
*
* See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
* ------------------------------------------------------------------------------------------------------------
*/

#include "../ConstexprConstructs.hpp"

#include <gtest/gtest.h>
#include <memory>

using namespace geos::compileTime;

TEST(ConstexprConstructs, StaticFor_SumOfIndices)
{
  constexpr int size = 10;
  int array[size] = {};
  auto fillArray = [&array](auto index)
  {
    array[index] = index;
  };

  static_for<0, size>(fillArray);

  int expectedSum = 0;
  for (int i = 0; i < size; ++i)
  {
    expectedSum += i;
  }

  int actualSum = 0;
  for (const auto & value : array)
  {
    actualSum += value;
  }

  EXPECT_EQ(expectedSum, actualSum);
}

TEST(ConstexprConstructs, ToString_Base10_Positive)
{
  constexpr auto str = to_string<12345>;
  constexpr auto buf = "12345";
  static_for<0,5>( [=]( auto idx )
  {
    constexpr auto ii = decltype(idx)::value;
    EXPECT_EQ( str[ii], buf[ii] );
  } );
}

TEST(ConstexprConstructs, ToString_Base10_Negative)
{
  constexpr auto str = to_string<-12345>;
  constexpr auto buf = "-12345";
  static_for<0,6>( [=]( auto idx )
  {
    constexpr auto ii = decltype(idx)::value;
    EXPECT_EQ( str[ii], buf[ii] );
  } );
}

TEST(ConstexprConstructs, ToString_Base2_Positive)
{
  constexpr auto str = to_string<13, 2>;
  constexpr auto buf = "1101";
  static_for<0,4>( [=]( auto idx )
  {
    constexpr auto ii = decltype(idx)::value;
    EXPECT_EQ( str[ii], buf[ii] );
  } );
}

TEST(ConstexprConstructs, ToString_Base16_Positive)
{
  constexpr auto str = to_string<255, 16>;
  constexpr auto buf = "FF";
  static_for<0,2>( [=]( auto idx )
  {
    constexpr auto ii = decltype(idx)::value;
    EXPECT_EQ( str[ii], buf[ii] );
  } );
}

TEST(ConstexprConstructs, ToString_NegativeZero)
{
  constexpr auto str = to_string<0>;
  constexpr auto buf = "0";
  static_for<0,1>( [=]( auto idx )
  {
    constexpr auto ii = decltype(idx)::value;
    EXPECT_EQ( str[ii], buf[ii] );
  } );
}

TEST(ConstexprConstructs, StrCat_SimpleConcatenation)
{
  constexpr auto result = strcat("Hello, ", "World!");
  EXPECT_STREQ(result, "Hello, World!");
}

TEST(ConstexprConstructs, StrCat_EmptyFirstString)
{
  constexpr auto result = strcat("", "Non-empty");
  EXPECT_STREQ(result, "Non-empty");
}

TEST(ConstexprConstructs, StrCat_EmptySecondString)
{
  constexpr auto result = strcat("Non-empty", "");
  EXPECT_STREQ(result, "Non-empty");
}