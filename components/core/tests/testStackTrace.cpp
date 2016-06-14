/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"
#include <cstring>
#include "slic/slic.hpp"
#include "../src/codingUtilities/stackTrace.hpp"
#include "../src/codingUtilities/SetSignalHandling.hpp"
#include <fenv.h>
#include <xmmintrin.h>
// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(testStackTrace,stackTrace)
{
  std::cout<<"calling handler0"<<std::endl;
  geosx::stacktrace::handler0(0);
}



