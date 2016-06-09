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

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(testStackTrace,stackTrace)
{
//	signal(SIGSEGV, stacktrace::handler);   // install our handler
//	stacktrace::foo(); // this will call foo, bar, and baz.  baz segfaults.
  geosx::stacktrace::handler( SIGSEGV, 0 );
}

