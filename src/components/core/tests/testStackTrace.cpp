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

#include <cstring>
#include "slic/slic.hpp"
#include "stackTrace.hpp"
#include "SetSignalHandling.hpp"
// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
//TEST(testStackTrace,stackTrace)
//{
//  std::cout<<"calling handler0"<<std::endl;
//  geosx::stacktrace::handler0(0);
//}
#include <execinfo.h>

void my_terminate(void);
namespace
{
// invoke set_terminate as part of global constant initialization
static const bool SET_TERMINATE = std::set_terminate(my_terminate);
}
void my_terminate()
{
  static bool tried_throw = false;

  try
  {
    // try once to re-throw currently active exception
    if( !tried_throw )
    {
      tried_throw = true;
      throw;
    }
  }
  catch( const std::exception &e )
  {
    std::cerr << __FUNCTION__ << " caught unhandled exception. what(): "
              << e.what()
              << std::endl;
  }
  catch( ... )
  {
    std::cerr << __FUNCTION__ << " caught unknown/unhandled exception."
              << std::endl;
  }

  void * array[50];
  int size = backtrace( array, 50 );

  std::cerr << __FUNCTION__ << " backtrace returned "
            << size
            << " frames\n\n";

  char ** messages = backtrace_symbols( array, size );

  for( int i = 0 ; i < size && messages != NULL ; ++i )
  {
    std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
  }
  std::cerr << std::endl;

  free( messages );

//  abort();
}

[[ noreturn ]] void throw_exception()
{
  // throw an unhandled runtime error
  throw std::runtime_error("RUNTIME ERROR!");
}

[[ noreturn ]] void foo2() {
  throw_exception();
}

[[ noreturn ]] void foo1() {
  foo2();
}
TEST(testStackTrace,uncaughtException)
{
  try
  {
    foo1();
  }
  catch(...)
  {
    std::cerr << __FUNCTION__ << " caught unknown/unhandled exception."
              << std::endl;
    void * array[50];
    int size = backtrace( array, 50 );

    std::cerr << __FUNCTION__ << " backtrace returned "
              << size
              << " frames\n\n";

    char ** messages = backtrace_symbols( array, size );

    for( int i = 0 ; i < size && messages != NULL ; ++i )
    {
      std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
    }
    std::cerr << std::endl;

    free( messages );
  }
}
