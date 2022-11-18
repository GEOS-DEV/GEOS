/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2022 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2022 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2022 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include <gtest/gtest.h>

#include <cstdlib>
#include <string>

#include "common/TimingMacros.hpp"
#include <adiak.hpp>

TEST( CaliperSmoke, SmokeTest )
{

  bool preset = false;

  // Check for CALI_CONFIG environment variable
  if( getenv( "CALI_CONFIG" ) != NULL )
  {
    std::string config( getenv( "CALI_CONFIG" ));
    if( config.find( "runtime-report" ) != std::string::npos )
    {
      preset = true;
    }
    else
    {
      FAIL() << "testCaliperSmoke failed - environment variable "
             << " CALI_CONFIG must contain 'runtime-report' or be unset";
    }
  }
  else
  {
    setenv( "CALI_CONFIG", "runtime-report", 0 );
  }

  GEOSX_CALIPER_MARK_FUNCTION_BEGIN;
  EXPECT_STRNE( cali_caliper_version(), NULL );
  GEOSX_CALIPER_MARK_FUNCTION_END;

  adiak::init( nullptr );
  adiak::fini();

  if( !preset )
  {
    unsetenv( "CALI_CONFIG" );
  }
}
