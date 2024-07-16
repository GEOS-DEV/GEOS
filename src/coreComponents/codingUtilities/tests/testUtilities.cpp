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

#include "codingUtilities/Utilities.hpp"

#include <gtest/gtest.h>

#include <map>

using namespace geos;


TEST(Utilities, IsEqualWithTolerance)
{
    EXPECT_TRUE(isEqual(1.0, 1.0));
    EXPECT_TRUE(isEqual(1.0, 1.01, 0.1));
    EXPECT_FALSE(isEqual(1.0, 2.0, 0.1));
}

TEST(Utilities, IsZeroWithTolerance)
{
    EXPECT_TRUE(isZero(0.0));
    EXPECT_TRUE(isZero(0.00001, 0.0001));
    EXPECT_FALSE(isZero(0.1, 0.001));
}

TEST(Utilities, IsOdd)
{
    EXPECT_TRUE(isOdd(1));
    EXPECT_TRUE(isOdd(-3));
    EXPECT_FALSE(isOdd(2));
}

TEST(Utilities, IsEven)
{
    EXPECT_TRUE(isEven(2));
    EXPECT_TRUE(isEven(-4));
    EXPECT_FALSE(isEven(3));
}

TEST(Utilities, ForEqualRanges)
{
    std::vector<int> values{1, 1, 2, 2, 3, 3, 3};
    int count = 0;
    auto lambda = [&count](auto start, auto end)
    {
        count += std::distance(start, end);
    };
    forEqualRanges(values.begin(), values.end(), lambda);
    EXPECT_EQ(count, values.size());
}

TEST(Utilities, ForUniqueValues)
{
    std::vector<int> values{1, 1, 2, 2, 2, 3};
    int numUnique = 0;
    auto lambda = [&numUnique](int value, int size)
    {
        numUnique++;
    };
    forUniqueValues(values.begin(), values.end(), lambda);
    EXPECT_EQ(numUnique, 3);  // Expects 3 unique values: 1, 2, and 3
}

TEST(Utilities, ForEachArgInTuple)
{
    std::tuple<int, double, std::string> myTuple{1, 2.0, "three"};
    std::vector<std::string> results;
    auto func = [&results](const auto& item, auto index)
    {
      if constexpr ( std::is_same_v< std::string const &, decltype( item ) > )
      {
        results.push_back( item + "-" + std::to_string(index));
      }
      else
      {
        results.push_back(std::to_string(item) + "-" + std::to_string(index));
      }
    };
    forEachArgInTuple(myTuple, func);
    std::vector<std::string> expected{"1-0", "2.000000-1", "three-2"};
    EXPECT_EQ(results, expected);
}

enum class Color : uint8_t { Red, Green, Blue };

TEST(Utilities, ToUnderlying)
{
    Color c = Color::Green;
    EXPECT_EQ(toUnderlying(c), static_cast<std::underlying_type_t<Color>>(1));
}

TEST(Utilities, ToUnderlyingPtr)
{
    Color c = Color::Blue;
    auto ptr = toUnderlyingPtr(&c);
    EXPECT_EQ(*ptr, static_cast<std::underlying_type_t<Color>>(2));
}


// Test for the copy function
TEST(Utilities, CopyVector)
{
    std::vector<int> src = {1, 2, 3, 4, 5};
    std::vector<int> dest(7, 0); // larger destination with an initial offset
    copy(5, src, dest, 2);
    EXPECT_EQ(dest, (std::vector<int>{0, 0, 1, 2, 3, 4, 5}));
}

// Test for applying the chain rule
TEST(Utilities, ApplyChainRule)
{
    const int N = 2;
    std::vector<std::vector<double>> dy_dx = {{1, 0}, {0, 1}}; // Identity matrix
    std::vector<double> df_dy = {1, 2}; // Derivatives of f with respect to y
    std::vector<double> df_dx(N, 0.0);
    applyChainRule(N, dy_dx, df_dy, df_dx, 0);
    EXPECT_EQ(df_dx, (std::vector<double>{1, 2}));
}

// Test for in-place chain rule application
TEST(Utilities, ApplyChainRuleInPlace)
{
    const int N = 2;
    std::vector<std::vector<double>> dy_dx = {{1, 0}, {0, 1}}; // Identity matrix
    std::vector<double> df_dy = {3, 4}; // Derivatives of f with respect to y
    std::vector<double> work(N, 0.0); // Work vector
    applyChainRuleInPlace(N, dy_dx, df_dy, work, 0);
    EXPECT_EQ(df_dy, (std::vector<double>{3, 4}));
}

// NoOpFunc test
TEST(Utilities, NoOpFuncTest)
{
    NoOpFunc noop;
    EXPECT_NO_THROW(noop(1, "test", 3.14)); // Should do nothing and not throw
}

// BitFlags test
TEST(Utilities, BitFlagsTest)
{
    BitFlags<Color> flags;
    flags.set(Color::Red);
    EXPECT_TRUE(flags.isSet(Color::Red));
    EXPECT_FALSE(flags.isSet(Color::Blue));
}

TEST( Utilities, MapExtraction )
{
  std::map< string, int > const m{
    { "k0", 0 },
    { "k1", 1 },
    { "k2", 2 }
  };

  EXPECT_EQ( mapKeys( m ), std::vector< string >( { "k0", "k1", "k2" } ) );
  EXPECT_EQ( mapKeys< std::set >( m ), std::set< string >( { "k0", "k1", "k2" } ) );
  EXPECT_EQ( mapValues( m ), std::vector< int >( { 0, 1, 2 } ) );
  EXPECT_EQ( mapValues< std::set >( m ), std::set< int >( { 0, 1, 2 } ) );
}
