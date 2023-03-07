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

// Source includes
#include "../StringUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace stringutilities;

TEST( testStringUtilities, tokenize )
{
  map< string, std::pair< std::vector< string >, std::vector< string > > >
  entries =
  {
    { "//entry0//entry1//entry2", { { "", "entry0", "entry1", "entry2" },
        { "", "", "entry0", "", "entry1", "", "entry2" } } },
    { "entry0//entry1/entry2", { { "entry0", "entry1", "entry2" },
        { "entry0", "", "entry1", "entry2" } } }
  };

  for( auto const & entry : entries )
  {
    string const & key = entry.first;
    std::vector< string > const & values0 = entry.second.first;
    std::vector< string > const & values1 = entry.second.second;

    std::vector< string > tokens0 = stringutilities::tokenize< std::vector< string > >( key, "/", true );
    std::vector< string > tokens1 = stringutilities::tokenize< std::vector< string > >( key, "/", false );


    EXPECT_TRUE( tokens0==values0 );
    EXPECT_TRUE( tokens1==values1 );
  }
}


TEST( testStringUtilities, toMetricPrefixString )
{
  double const values[36] = { 1.234567890e-15,
                              1.234567890e-14,
                              1.234567890e-13,
                              1.234567890e-12,
                              1.234567890e-11,
                              1.234567890e-10,
                              1.234567890e-9,
                              1.234567890e-8,
                              1.234567890e-7,
                              1.234567890e-6,
                              1.234567890e-5,
                              1.234567890e-4,
                              1.234567890e-3,
                              1.234567890e-2,
                              1.234567890e-1,
                              1.234567890e0,
                              1.234567890e1,
                              1.234567890e2,
                              1.234567890e3,
                              1.234567890e4,
                              1.234567890e5,
                              1.234567890e6,
                              1.234567890e7,
                              1.234567890e8,
                              1.234567890e9,
                              1.234567890e10,
                              1.234567890e11,
                              1.234567890e12,
                              1.234567890e13,
                              1.234567890e14,
                              1.234567890e15,
                              1.234567890e16,
                              1.234567890e17,
                              1.234567890e18,
                              1.234567890e19,
                              1.234567890e20 };

  string const answer[36] = { " 1.23 f",
                              " 12.3 f",
                              "  123 f",
                              " 1.23 p",
                              " 12.3 p",
                              "  123 p",
                              " 1.23 n",
                              " 12.3 n",
                              "  123 n",
                              " 1.23 u",
                              " 12.3 u",
                              "  123 u",
                              " 1.23 m",
                              " 12.3 m",
                              "  123 m",
                              " 1.23  ",
                              " 12.3  ",
                              "  123  ",
                              " 1.23 K",
                              " 12.3 K",
                              "  123 K",
                              " 1.23 M",
                              " 12.3 M",
                              "  123 M",
                              " 1.23 G",
                              " 12.3 G",
                              "  123 G",
                              " 1.23 T",
                              " 12.3 T",
                              "  123 T",
                              " 1.23 P",
                              " 12.3 P",
                              "  123 P",
                              " 1.23 E",
                              " 12.3 E",
                              "  123 E" };


  for( int a=0; a<36; ++a )
  {
    std::string const result = toMetricPrefixString( values[a] );
    std::string const negResult = toMetricPrefixString( -values[a] );

    EXPECT_STRCASEEQ( result.c_str(), answer[a].c_str() );
    std::string negAnswer = answer[a];
    int const sign = negAnswer.find_first_not_of( ' ' );
    negAnswer[sign-1] = '-';
    EXPECT_STRCASEEQ( negResult.c_str(), negAnswer.c_str() );
  }
}
