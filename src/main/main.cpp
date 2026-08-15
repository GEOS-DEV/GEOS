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

// Source includes
#include "common/logger/ErrorHandling.hpp"
#include "common/logger/Logger.hpp"
#include "common/MemoryInfos.hpp"
#include "common/TimingMacros.hpp"
#include "common/Units.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/version.hpp"

#include <fenv.h>
#include <csignal>      // for signal()
#include <iostream>

using namespace geos;

static void
fpe_signal_handler( int sig, siginfo_t *sip, void *scp )
{
  GEOS_UNUSED_VAR( sig );
  GEOS_UNUSED_VAR( scp );
  int fe_code = sip->si_code;
  printf( "In signal handler : " );
  if( fe_code == ILL_ILLTRP )
    printf( "Illegal trap detected\n" );
  else
    printf( "Code detected : %d\n", fe_code );
  abort();
}
void enable_floating_point_exceptions()
{
  fenv_t env;
  fegetenv( &env );
  env.__fpcr = env.__fpcr | __fpcr_trap_invalid; 
 fesetenv( &env );
  struct sigaction act;
  act.sa_sigaction = fpe_signal_handler;
  sigemptyset ( &act.sa_mask );
  act.sa_flags = SA_SIGINFO;
  sigaction( SIGILL, &act, NULL );
}

int main( int argc, char *argv[] )
{
  try
  {
    std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();

    std::unique_ptr< CommandLineOptions > commandLineOptions = basicSetup( argc, argv, true );

    outputVersionInfo();

    GEOS_LOG_RANK_0( GEOS_FMT( "Started at {:%Y-%m-%d %H:%M:%S}", startTime ) );

    std::chrono::system_clock::duration initTime;
    std::chrono::system_clock::duration runTime;
    {
      GeosxState state( std::move( commandLineOptions ) );

      bool const problemToRun = state.initializeDataRepository();
      if( problemToRun )
      {
        state.applyInitialConditions();
        state.run();
        GEOS_WARNING_IF( state.getState() != State::COMPLETED, "Simulation exited early." );
      }

      initTime = state.getInitTime();
      runTime = state.getRunTime();
    }

    MemoryLogging::getInstance().memoryStatsReport();

    basicCleanup( false );

    std::chrono::system_clock::time_point endTime = std::chrono::system_clock::now();
    std::chrono::system_clock::duration totalTime = endTime - startTime;

    GEOS_LOG_RANK_0( GEOS_FMT( "Finished at {:%Y-%m-%d %H:%M:%S}", endTime ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "total time            {}", units::TimeFormatInfo::fromDuration( totalTime ) ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "initialization time   {}", units::TimeFormatInfo::fromDuration( initTime ) ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "run time              {}", units::TimeFormatInfo::fromDuration( runTime ) ) );

    return 0;
  }
  // A NotAnError is thrown if "-h" or "--help" option is used.
  catch( NotAnError const & )
  {
    basicCleanup( false );
    return 0;
  }
  catch( geos::Exception & e )
  { // GEOS generated exceptions management
    ErrorLogger::global().flushCurrentExceptionMessage();
    basicCleanup( true );
    LvArray::system::callErrorHandler();
  }
  catch( std::exception const & e )
  { // native exceptions management
    ErrorLogger::global().flushErrorMsg( ErrorLogger::global().initCurrentExceptionMessage(
                                           MsgType::Exception, e.what(),
                                           ::geos::logger::internal::g_rank )
                                           .addCallStackInfo( LvArray::system::stackTrace( true ) )
                                           .getDiagnosticMsg());
    basicCleanup( true );
    LvArray::system::callErrorHandler();
  }
  return 0;
}
