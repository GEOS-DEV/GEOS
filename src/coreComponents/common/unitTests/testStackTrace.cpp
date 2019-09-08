/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "gtest/gtest.h"

#include <cstring>
#include "stackTrace.hpp"
#include "SetSignalHandling.hpp"
// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
//TEST(testStackTrace,stackTrace)
//{
//  geosx::stacktrace::handler0(0);
//}
#include <execinfo.h>

void my_terminate( void );
namespace
{
// invoke set_terminate as part of global constant initialization
static const bool SET_TERMINATE = std::set_terminate( my_terminate );
}
void my_terminate()
{

  cxx_utilities::handler( 1, 1, 1 );

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
  catch( const std::exception & e )
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

  char * * messages = backtrace_symbols( array, size );

  for( int i = 0 ; i < size && messages != nullptr ; ++i )
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
  throw std::runtime_error( "RUNTIME ERROR!" );
}

[[ noreturn ]] void foo2() {
  throw_exception();
}

[[ noreturn ]] void foo1() {
  foo2();
}
TEST( testStackTrace, uncaughtException )
{

  cxx_utilities::setSignalHandling( cxx_utilities::handler1 );
}
