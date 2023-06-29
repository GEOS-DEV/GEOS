/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "common/DataTypes.hpp"
#include "common/Format.hpp"
#include "common/TimingMacros.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"


using namespace geos;


int main( int argc, char *argv[] )
{
  try
  {
    std::chrono::system_clock::time_point const startTime = std::chrono::system_clock::now();

    std::unique_ptr< CommandLineOptions > commandLineOptions = basicSetup( argc, argv, true );

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
        LVARRAY_WARNING_IF( state.getState() != State::COMPLETED, "Simulation exited early." );
      }

      initTime = state.getInitTime();
      runTime = state.getRunTime();
    }

    basicCleanup();

    std::chrono::system_clock::time_point const endTime = std::chrono::system_clock::now();
    std::chrono::system_clock::duration const totalTime = endTime - startTime;

    GEOS_LOG_RANK_0( GEOS_FMT( "Finished at {:%Y-%m-%d %H:%M:%S}", endTime ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "total time            {:%H:%M:%S}", totalTime ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "initialization time   {:%H:%M:%S}", initTime ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "run time              {:%H:%M:%S}", runTime ) );

    return 0;
  }
  // A NotAnError is thrown if "-h" or "--help" option is used.
  catch( NotAnError const & )
  {
    return 0;
  }
  catch( std::exception const & e )
  {
    GEOS_LOG( e.what() );
    LvArray::system::callErrorHandler();
    std::abort();
  }
  return 0;
}
