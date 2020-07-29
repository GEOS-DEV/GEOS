/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "managers/ProblemManager.hpp"
#include "managers/initialization.hpp"
#include "managers/GeosxState.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

TEST( testXML, testXML )
{
  geosx::getGlobalState().getProblemManager().ParseCommandLineInput();
  geosx::getGlobalState().getProblemManager().ParseInputFile();
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv, true ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
