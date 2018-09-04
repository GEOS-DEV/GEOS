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
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include <fenv.h>
#include <xmmintrin.h>
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

TEST(testStackTrace_DeathTest, stackTrace)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wused-but-marked-unused"
  EXPECT_DEATH_IF_SUPPORTED(testStackTrace(0), IGNORE_OUTPUT);
#pragma GCC diagnostic pop
}
