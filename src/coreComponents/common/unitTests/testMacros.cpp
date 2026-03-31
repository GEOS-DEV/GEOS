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
#include "common/GeosxMacros.hpp"
// TPL includes
#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>
// test dependancies
#include <array>
#include "common/StdContainerWrappers.hpp"

TEST( testMacros, testArgumentCount )
{
  EXPECT_EQ( 0, int( GEOS_DETAIL_MORE_THAN_ONE_ARG( 'a' ) ) );
  EXPECT_EQ( 1, int( GEOS_DETAIL_MORE_THAN_ONE_ARG( 'a', 'b' ) ) );
  EXPECT_EQ( 1, int( GEOS_DETAIL_MORE_THAN_ONE_ARG( 'a', 'b', 'c' ) ) );

  EXPECT_EQ( 1, int( GEOS_DETAIL_MORE_THAN_ONE_ARG( 'a', 'b', 'c', 'd',
                                                    'e', 'f', 'g', 'h',
                                                    'i', 'j', 'k', 'l',
                                                    'w', 'x', 'y', 'z' ) ) );

  // Expected out of bound (>16 params): wrongly cast the last '!' to integer type
  EXPECT_EQ( 33, int( GEOS_DETAIL_MORE_THAN_ONE_ARG( 'a', 'b', 'c', 'd',
                                                     'e', 'f', 'g', 'h',
                                                     'i', 'j', 'k', 'l',
                                                     'w', 'x', 'y', 'z', '!' ) ) );
}


TEST( testMacros, testFirstArgument )
{
  EXPECT_EQ( 'a', GEOS_DETAIL_FIRST_ARG( 'a' ) );
  EXPECT_EQ( 'a', GEOS_DETAIL_FIRST_ARG( 'a', 'b' ) );
  EXPECT_EQ( 'a', GEOS_DETAIL_FIRST_ARG( 'a', 'b', 'c' ) );

  EXPECT_EQ( 'a', GEOS_DETAIL_FIRST_ARG( 'a', 'b', 'c', 'd',
                                         'e', 'f', 'g', 'h',
                                         'i', 'j', 'k', 'l',
                                         'w', 'x', 'y', 'z' ) );

  // Out of bound (>16 params): Cannot compile here, not testable.
  // EXPECT_EXIT( GEOS_DETAIL_FIRST_ARG( 'a', 'b', 'c', 'd',
  //                                     'e', 'f', 'g', 'h',
  //                                     'i', 'j', 'k', 'l',
  //                                     'w', 'x', 'y', 'z', '!' ) );
}

TEST( testMacros, testRestArguments )
{
  // no param after the first -> double(void) called -> 0.0 value
  EXPECT_EQ( 0.0, double( GEOS_DETAIL_REST_ARGS( 1.0 ) ) );

  EXPECT_EQ( 2.0, double( GEOS_DETAIL_REST_ARGS( 1.0, 2.0 ) ) );

  {
    auto const generatedArray = stdArray< double, 3 >{ GEOS_DETAIL_REST_ARGS( 1.0, 2.0, 3.0, 4.0 ) };
    auto const expectedArray = stdArray< double, 3 >{ 2.0, 3.0, 4.0 };
    EXPECT_EQ( generatedArray, expectedArray );
  }

  {
    auto const generatedArray = stdArray< double, 16 >{ GEOS_DETAIL_REST_ARGS( 01.0, 02.0, 03.0, 04.0, 05.0, 06.0, 07.0, 08.0,
                                                                                 09.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 ) };
    auto const expectedArray = stdArray< double, 16 >{ 02.0, 03.0, 04.0, 05.0, 06.0, 07.0, 08.0, 09.0,
                                                         10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
    EXPECT_EQ( generatedArray, expectedArray );
  }

  // { // Out of bound (>16 params): Cannot compile here, not testable.
  //   auto const generatedArray = stdArray< double, 3 >{ GEOS_DETAIL_REST_ARGS( 01.0, 02.0, 03.0, 04.0, 05.0, 06.0, 07.0, 08.0,
  //                                                                               09.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0 ) };
  //   auto const expectedArray = stdArray< double, 3 >{ 02.0, 03.0, 04.0, 05.0, 06.0, 07.0, 08.0,
  //                                                       09.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0 };
  //   EXPECT_EQ( generatedArray, expectedArray );
  // }
}

int main( int argc, char * * argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
