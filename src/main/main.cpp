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

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include "SetFPE.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/ProblemManager.hpp"


#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace geosx;


int main( int argc, char *argv[] )
{
  timeval tim;
  gettimeofday(&tim, nullptr);
  real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);

#ifdef GEOSX_USE_MPI
  int rank = 0;
  int nranks = 1;

  MPI_Init(&argc,&argv);

  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );

  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);

  MPI_Comm_size(MPI_COMM_GEOSX, &nranks);

  logger::InitializeLogger(MPI_COMM_GEOSX);
#else
  logger::InitializeLogger():
#endif

#ifdef GEOSX_USE_OPENMP
  {
    GEOS_LOG_RANK_0("Number of threads: " << omp_get_max_threads());
  }
#endif

  // cxx_utilities::setSignalHandling(cxx_utilities::handler1);
  cxx_utilities::SetFPE();

  std::string restartFileName;
  bool restart = ProblemManager::ParseRestart( argc, argv, restartFileName );
  if (restart) {
    GEOS_LOG_RANK_0("Loading restart file " << restartFileName);
    dataRepository::SidreWrapper::reconstructTree( restartFileName, "sidre_hdf5", MPI_COMM_GEOSX );
  }

  ProblemManager problemManager( "Problem", nullptr );

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( argc, argv );

  if ( !problemManager.getSchemaFileName().empty() )
  {
    problemManager.GenerateDocumentation();
  }
  else
  {
    problemManager.ParseInputFile();

    problemManager.ProblemSetup();

    if (restart) {
      problemManager.ReadRestartOverwrite( restartFileName );
    }

    MPI_Barrier(MPI_COMM_GEOSX);
    GEOS_LOG_RANK_0("Running simulation");

    gettimeofday(&tim, nullptr);
    const real64 t_initialize = tim.tv_sec + (tim.tv_usec / 1000000.0);

    problemManager.RunSimulation();

    gettimeofday(&tim, nullptr);
    const real64 t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

    GEOS_LOG_RANK_0("\ninit time = " << std::setprecision(5) << t_initialize-t_start <<
                    "s, run time = " << t_run-t_initialize << "s");
  }
  

  problemManager.ClosePythonInterpreter();

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return 0;
}
