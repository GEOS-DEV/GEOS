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
  double const values[21] = { 1.234567890e0,
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

  string const answer[21] = { " 1.23  ",
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

  string const negAnswer[21] = { "-1.23  ",
                                 "-12.3  ",
                                 " -123  ",
                                 "-1.23 K",
                                 "-12.3 K",
                                 " -123 K",
                                 "-1.23 M",
                                 "-12.3 M",
                                 " -123 M",
                                 "-1.23 G",
                                 "-12.3 G",
                                 " -123 G",
                                 "-1.23 T",
                                 "-12.3 T",
                                 " -123 T",
                                 "-1.23 P",
                                 "-12.3 P",
                                 " -123 P",
                                 "-1.23 E",
                                 "-12.3 E",
                                 " -123 E" };


  for( int a=0; a<21; ++a )
  {
    std::string const result = toMetricPrefixString( values[a] );
    std::string const negResult = toMetricPrefixString( -values[a] );

    EXPECT_STRCASEEQ( result.c_str(), answer[a].c_str() );
    EXPECT_STRCASEEQ( negResult.c_str(), negAnswer[a].c_str() );
  }
}
