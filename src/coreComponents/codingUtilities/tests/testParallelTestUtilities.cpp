/*	
 * ------------------------------------------------------------------------------------------------------------	
 * SPDX-License-Identifier: LGPL-2.1-only	
 *	
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC	
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University	
 * Copyright (c) 2018-2020 Total, S.A	
 * Copyright (c) 2020-     GEOSX Contributors	
 * All right reserved	
 *	
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.	
 * ------------------------------------------------------------------------------------------------------------	
 */

#include "mainInterface/initialization.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

#include <gtest/gtest.h>

TEST( UnitTestUtilities, skipInSerial )
{
  SKIP_TEST_IN_SERIAL( "This test should not be run in serial." );

  if( geosx::MpiWrapper::commSize() == 1 )
  {
    GTEST_FAIL();
  }
  else
  {
    GTEST_SUCCEED();
  }
}

TEST( UnitTestUtilities, skipInParallel )
{
  SKIP_TEST_IN_PARALLEL( "This test should not be run in parallel." );

  if( geosx::MpiWrapper::commSize() != 1 )
  {
    GTEST_FAIL();
  }
  else
  {
    GTEST_SUCCEED();
  }
}

TEST( UnitTestUtilities, expected )
{
  using namespace geosx;
  using namespace geosx::testing;

  if( MpiWrapper::commSize() == 1 )
  {
    ASSERT_EQ( 1, expected( 1, {} ) );
    ASSERT_EQ( 1, expected( 1, { 2 } ) );
    ASSERT_EQ( 1, expected( 1, { 2, 3 } ) );
  }
  else
  {
    ASSERT_EQ( MpiWrapper::commRank() == 0 ? 2 : 3, expected( 1, { 2, 3 } ) );
  }
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geosx::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
