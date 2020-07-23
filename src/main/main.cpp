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

// Source includes
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/GeosxState.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

// TPL includes
#ifdef GEOSX_USE_CALIPER
#include <caliper/cali-manager.h>
#endif

// System includes
#include <chrono>

using namespace geosx;


int main( int argc, char *argv[] )
{
  std::chrono::system_clock::time_point const startTime = std::chrono::system_clock::now();

  basicSetup( argc, argv, true );

  std::chrono::system_clock::duration initTime;
  std::chrono::system_clock::duration runTime;
  {
    GeosxState state;

    bool const problemToRun = state.initializeDataRepository();
    if ( problemToRun )
    { LVARRAY_WARNING_IF( state.run(), "Simulation exited early." ); }

    initTime = state.getInitTime();
    runTime = state.getRunTime();
  }

  basicCleanup();

  std::chrono::system_clock::duration const totalTime = std::chrono::system_clock::now() - startTime;

  GEOSX_LOG_RANK_0( "total time          " << durationToString( totalTime ) );
  GEOSX_LOG_RANK_0( "initialization time " << durationToString( initTime ) );
  GEOSX_LOG_RANK_0( "run time            " << durationToString( runTime ) );

  return 0;
}
