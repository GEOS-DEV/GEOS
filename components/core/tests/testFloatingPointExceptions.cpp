/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wexit-time-destructors"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#endif

#include "../src/codingUtilities/SetSignalHandling.hpp"
#include "../src/codingUtilities/stackTrace.hpp"
#include <fenv.h>
#include <xmmintrin.h>
// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------

void func3()
{
  std::cout<<"Perform 1.0/0.0"<<std::endl;
  double a = 1.0/0.0;
  std::cout<<"1.0/0.0 didn't kill program, result is "<<a<<std::endl<<std::endl;
}

void func2()
{
  func3();
}

void func1()
{
  func2();
}

void func0()
{
  func1();
}
TEST(testStackTrace,stackTrace)
{

  geosx::setSignalHandling(geosx::stacktrace::handler1);
  func0();

//  exit(0);
}



