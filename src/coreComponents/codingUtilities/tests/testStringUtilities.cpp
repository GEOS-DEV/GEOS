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

    for( auto const & value : values0 )
    {
      std::cout<<value<<",";
    }
    std::cout<<std::endl;
    for( auto const & value : tokens0 )
    {
      std::cout<<value<<",";
    }
    std::cout<<std::endl;

    for( auto const & value : values1 )
    {
      std::cout<<value<<",";
    }
    std::cout<<std::endl;
    for( auto const & value : tokens1 )
    {
      std::cout<<value<<",";
    }
    std::cout<<std::endl;

    EXPECT_TRUE( tokens0==values0 );
    EXPECT_TRUE( tokens1==values1 );
  }
}
