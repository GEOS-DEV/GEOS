/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "dataRepository/Group.hpp"
// TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

using namespace geos;

TEST( testDataTypes, testBoundChecking )
{
  internal::StdVectorWrapper< std::string,
                              std::allocator< std::string >,
                              true > vectorBoundsChecking = {"test"};
  EXPECT_THROW( {
    try
    {
      std::string crash = vectorBoundsChecking[1];
    }
    catch( const std::out_of_range & e )
    {
      throw;
    }
  }, std::out_of_range );

  internal::StdMapWrapper< std::map< integer, integer >,
                           true > mapBoundsChecking{{0, 1}};
  EXPECT_THROW( {
    try
    {
      integer crash = mapBoundsChecking[1];
    }
    catch( const std::out_of_range & e )
    {
      throw;
    }
  }, std::out_of_range );

  internal::StdMapWrapper< std::unordered_map< integer, integer >,
                           true > unorderedMapBoundsChecking{{0, 1}};
  EXPECT_THROW( {
    try
    {
      integer crash = unorderedMapBoundsChecking[1];
    }
    catch( const std::out_of_range & e )
    {
      throw;
    }
  }, std::out_of_range );

}

TEST( testDataTypes, testNoBoundChecking )
{
  internal::StdVectorWrapper< std::string,
                              std::allocator< std::string >,
                              false > boundChecking = {"test"};

  EXPECT_NO_THROW( {
    std::string crash = boundChecking[1];
  } );

  internal::StdMapWrapper< std::map< integer, integer >,
                           false > mapBoundsChecking{{0, 1}};

  EXPECT_NO_THROW( {
    integer crash = mapBoundsChecking[1];
  } );

  internal::StdMapWrapper< std::unordered_map< integer, integer >,
                           false > unorderedMapBoundsChecking{{0, 1}};

  EXPECT_NO_THROW( {
    integer crash = unorderedMapBoundsChecking[1];
  } );

}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
