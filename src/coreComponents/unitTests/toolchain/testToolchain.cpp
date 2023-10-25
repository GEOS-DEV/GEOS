/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOS Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include <gtest/gtest.h>

TEST( Toolchain, NDEBUGfromTPls )
{
  /*
   * This test guards against spurious propagation of -DNDEBUG preprocessor flag (HDF5 in our case),
   * which has the effect of disabling LvArray assertions:
   * instead, we always check that this test fails in CMake Debug mode, whilst it should always
   * pass in RelWithDebInfo or Release builds.
   */
#ifdef NDEBUG
  SUCCEED();  // RelWithDebInfo or Release builds
#else
  FAIL();  // Debug builds
#endif
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  int const result = RUN_ALL_TESTS();

  return result;
}
