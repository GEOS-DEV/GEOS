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
#include <fenv.h>
#include <xmmintrin.h>
// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(testStackTrace,stackTrace)
{
  signal(SIGHUP, geosx::stacktrace::handler0);
  signal(SIGINT, geosx::stacktrace::handler0);
  signal(SIGQUIT, geosx::stacktrace::handler0);
  signal(SIGILL, geosx::stacktrace::handler0);
  signal(SIGTRAP, geosx::stacktrace::handler0);
  signal(SIGABRT, geosx::stacktrace::handler0);
#if  (defined(_POSIX_C_SOURCE) && !defined(_DARWIN_C_SOURCE))
  signal(SIGPOLL, geosx::stacktrace::handler0);
#else
  signal(SIGIOT, geosx::stacktrace::handler0);
  signal(SIGEMT, geosx::stacktrace::handler0);
#endif
  signal(SIGFPE, geosx::stacktrace::handler0);
  signal(SIGKILL, geosx::stacktrace::handler0);
  signal(SIGBUS, geosx::stacktrace::handler0);
  signal(SIGSEGV, geosx::stacktrace::handler0);
  signal(SIGSYS, geosx::stacktrace::handler0);
  signal(SIGPIPE, geosx::stacktrace::handler0);
  signal(SIGALRM, geosx::stacktrace::handler0);
  signal(SIGTERM, geosx::stacktrace::handler0);
  signal(SIGURG, geosx::stacktrace::handler0);
  signal(SIGSTOP, geosx::stacktrace::handler0);
  signal(SIGTSTP, geosx::stacktrace::handler0);
  signal(SIGCONT, geosx::stacktrace::handler0);
  signal(SIGCHLD, geosx::stacktrace::handler0);

#ifdef __APPLE__// && __MACH__
//  _MM_SET_EXCEPTION_MASK(  _MM_EXCEPT_OVERFLOW);
  _MM_SET_EXCEPTION_MASK(_MM_EXCEPT_INVALID | _MM_EXCEPT_DIV_ZERO | _MM_EXCEPT_OVERFLOW);
  //  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
  //fesetenv(FE_ALL_EXCEPT);
//#endif

//#ifdef __INTEL_COMPILER
#else
  feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#endif
//	signal(SIGSEGV, stacktrace::handler);   // install our handler
//	stacktrace::foo(); // this will call foo, bar, and baz.  baz segfaults.

//  geosx::stacktrace::handler( SIGSEGV, 0 );

  double temp;
  temp = 1.0/0.0;

  exit(1);
}

