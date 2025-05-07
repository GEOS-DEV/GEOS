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

#include "common/logger/ErrorHandling.hpp"
#include "dataRepository/DataContext.hpp"

#include <gtest/gtest.h>

using namespace geos;

TEST( DataContext, testCompleteYaml )
{
  geos::ErrorLogger errorLogger; 
  int x = 5;
  geos::dataRepository::DataFileContext dataContext( "targetName",
                                                     "test1_file.xml",
                                                     42 );
  GEOS_THROW_CTX_IF( dataContext, x==5, "Here is the error message", std::runtime_error );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}
