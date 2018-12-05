/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "SetFPE.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include <fenv.h>
#include <xmmintrin.h>
#include <cmath>
#include <float.h>
// API coverage tests
// Each test should be documented with the interface functions being tested

const char IGNORE_OUTPUT[] = ".*";

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------

void func3(double divisor)
{
  double a = 1.0 / divisor;
  EXPECT_TRUE(false) << "1.0/0.0 didn't kill program, result is " << a;
}

void func2(double divisor)
{
  func3(divisor);
}

void func1(double divisor)
{
  func2(divisor);
}

void func0(double divisor)
{
  func1(divisor);
}

void testStackTrace(double divisor)
{
  cxx_utilities::setSignalHandling(cxx_utilities::handler1);
  func0(divisor);
}

//TEST(testStackTrace_DeathTest, stackTrace)
//{
//   EXPECT_DEATH_IF_SUPPORTED(testStackTrace(0), IGNORE_OUTPUT);
//}

double uf_test(double x, double denominator)
{
  return x/denominator;
}


TEST( TestFloatingPointEnvironment, test_FE_UNDERFLOW )
{
  cxx_utilities::UnsetUnderflowFlush();
  int temp = fegetexcept();
  feenableexcept( FE_ALL_EXCEPT );

  EXPECT_DEATH_IF_SUPPORTED( uf_test(DBL_MIN, 2), "");

  double normalNum = DBL_MIN*2;
  EXPECT_TRUE( std::fpclassify( normalNum ) == FP_NORMAL );

  fedisableexcept(FE_ALL_EXCEPT);
  feenableexcept( temp );
}


TEST( TestFloatingPointEnvironment, test_FE_UNDERFLOW_flush )
{
  cxx_utilities::SetFPE();

  double fpnum = uf_test(DBL_MIN, 2);
  int fpclassification = std::fpclassify( fpnum );
  EXPECT_TRUE( fpclassification != FP_SUBNORMAL );

}

TEST( TestFloatingPointEnvironment, test_FE_DIVBYZERO )
{
  cxx_utilities::SetFPE();
  EXPECT_DEATH_IF_SUPPORTED( func3(0.0) , "");
}

double of_test( double x, double y )
{
  return x*y;
}
TEST( TestFloatingPointEnvironment, test_FE_OVERFLOW )
{
  cxx_utilities::SetFPE();
  EXPECT_DEATH_IF_SUPPORTED( of_test(100,DBL_MAX) , "");
}

TEST( TestFloatingPointEnvironment, test_FE_INVALID )
{
  cxx_utilities::SetFPE();
  EXPECT_DEATH_IF_SUPPORTED( double junk0 = std::acos(2.0); , "");
}
