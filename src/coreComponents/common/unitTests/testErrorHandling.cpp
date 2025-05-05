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

#include <gtest/gtest.h>

using namespace geos;

TEST( ErrorHandling, testYaml )
{
  ErrorLogger logger; 
  ErrorLogger::ErrorMsg msgStruct = logger.errorMsgFormatter( ErrorLogger::MsgType::Error, "msg content", "dev file name", 24, "input file name", 42 );
  logger.errorMsgWritter( msgStruct );
}

int main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  int const result = RUN_ALL_TESTS();
  return result;
}