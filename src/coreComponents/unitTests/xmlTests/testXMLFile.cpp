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
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/FieldSpecificationBase.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;

TEST( testXML, testXMLFile )
{
  geosx::ProblemManager & problemManager = geosx::getGlobalState().getProblemManager();
  problemManager.parseCommandLineInput();
  problemManager.parseInputFile();

  // Check that we've read the full XML with all nested includes by inspecting boundary conditions
  EXPECT_TRUE( problemManager.getFieldSpecificationManager().hasGroup< FieldSpecificationBase >( "v0" ) );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  GeosxState state( basicSetup( argc, argv, true ) );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
