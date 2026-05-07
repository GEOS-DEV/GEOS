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
#include "common/logger/StackTrace.hpp"
#include "common/GeosxConfig.hpp"
#include "common/MemoryInfos.hpp"
#include "common/TimingMacros.hpp"
#include "common/Units.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/version.hpp"
#ifdef GEOS_USE_CPPTRACE
#include <cpptrace/from_current.hpp>
#endif


using namespace geos;


int main( int argc, char *argv[] )
{
  auto runMain = [&]
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
  };

  // A NotAnError is thrown if "-h" or "--help" option is used.
  auto onNotAnError = [&]( NotAnError const & )
  {
    basicCleanup( false );
  };

  auto onGeosException = [&]( geos::Exception & )
  { // GEOS generated exceptions management
#ifdef GEOS_USE_CPPTRACE
    std::string const stacktrace = StackTrace::formatStackTrace( cpptrace::from_current_exception() );
#else
    std::string const stacktrace = LvArray::system::stackTrace( true );
#endif
    if( !stacktrace.empty() )
    {
      ErrorLogger::global().modifyCurrentExceptionMessage().addCallStackInfo( stacktrace );
    }
    ErrorLogger::global().flushCurrentExceptionMessage();
    basicCleanup( true );
    LvArray::system::callErrorHandler();
  };

  auto onStdException = [&]( std::exception const & e )
  { // native exceptions management
#ifdef GEOS_USE_CPPTRACE
    std::string const stacktrace = StackTrace::formatStackTrace( cpptrace::from_current_exception() );
#else
    std::string const stacktrace = LvArray::system::stackTrace( true );
#endif
    ErrorLogger::global().flushErrorMsg( ErrorLogger::global().initCurrentExceptionMessage(
                                           MsgType::Exception, e.what(),
                                           ::geos::logger::internal::g_rank )
                                           .addCallStackInfo( stacktrace )
                                           .getDiagnosticMsg());
    basicCleanup( true );
    LvArray::system::callErrorHandler();
  };

#ifdef GEOS_USE_CPPTRACE
  cpptrace::try_catch( runMain, onNotAnError, onGeosException, onStdException );
#else
  try
  {
    runMain();
  }
  catch( NotAnError const & )
  {
    onNotAnError();
  }
  catch( geos::Exception & )
  {
    onGeosException();
  }
  catch( std::exception const & e )
  {
    onStdException( e );
  }
#endif

  return 0;
}
