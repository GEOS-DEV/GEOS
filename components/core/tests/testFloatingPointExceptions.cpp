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
#include "../src/codingUtilities/SetSignalHandling.hpp"
#include "../src/codingUtilities/stackTrace.hpp"
#include <fenv.h>
#include <xmmintrin.h>
// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(testStackTrace,stackTrace)
{

  geosx::setSignalHandling(geosx::stacktrace::handler1);

  std::cout<<"Perform 1.0/0.0"<<std::endl;
  double a = 1.0/0.0;
  std::cout<<"1.0/0.0 didn't kill program, result is "<<a<<std::endl<<std::endl;

//  exit(0);
}



