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



#include "common/Logger.hpp"
#include "common/TimingMacros.hpp"
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include "../dataRepository/SidreWrapper.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/ProblemManager.hpp"

#include <functional>
#include "coupling/TribolCoupling.hpp"


#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace geosx;

#ifdef USE_ATK
using namespace axom;
#endif



int main( int argc, char *argv[] )
{
  timeval tim;
  gettimeofday(&tim, nullptr);
  real64 t_start = tim.tv_sec + (tim.tv_usec / 1000000.0);
  real64 t_initialize, t_run;

#ifdef USE_MPI
  int rank;
  MPI_Init(&argc,&argv);

#if 1
  std::string cmdline ;
  // Combine the command line arguments into a single string.
  for (int i = 0 ; i < argc ; ++i) {
      cmdline += argv[i] ;
  }
  // Hash the command line.
  int codeID = std::hash<std::string>{}(cmdline) ;

  MPI_Comm MPI_OTHER_COMM ;
  TribolCoupling::InitCommSubset(MPI_COMM_WORLD, &MPI_COMM_GEOSX, &MPI_OTHER_COMM, codeID) ;
#else
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
#endif

  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
#endif

  std::cout<<"starting main"<<std::endl;

#ifdef USE_OPENMP
  {
    int noThreads = omp_get_max_threads();
    std::cout<<"No of threads: "<<noThreads<<std::endl;
  }
#endif




#ifdef USE_ATK
  slic::initialize();
  std::string format =  std::string( "***********************************\n" )+
                       std::string( "* <TIMESTAMP>\n\n" ) +
                       std::string( "* LEVEL=<LEVEL>\n" ) +
                       std::string( "* MESSAGE=<MESSAGE>\n" ) +
                       std::string( "* FILE=<FILE>\n" ) +
                       std::string( "* LINE=<LINE>\n" ) +
                       std::string( "***********************************\n" );
  slic::setLoggingMsgLevel( slic::message::Debug );
  slic::GenericOutputStream * const stream = new slic::GenericOutputStream(&std::cout, format );
  slic::addStreamToAllMsgLevels( stream );

#endif
  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  // Mark begin of "initialization" phase
  GEOS_MARK_BEGIN("Initialization");

  std::string restartFileName;
  bool restart = ProblemManager::ParseRestart( argc, argv, restartFileName );
  if (restart) {
    std::cout << "Loading restart file " << restartFileName << std::endl;
    dataRepository::SidreWrapper::reconstructTree( restartFileName, "sidre_hdf5", MPI_COMM_GEOSX );
  }

  ProblemManager problemManager( "ProblemManager", nullptr );
  problemManager.SetDocumentationNodes();
  problemManager.RegisterDocumentationNodes();  

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( argc, argv );

  problemManager.ParseInputFile();


  problemManager.Initialize( &problemManager );


  GEOS_MARK_END("Initialization");

  problemManager.ApplyInitialConditions();

  problemManager.FinalInitializationRecursive( &problemManager );

  if (restart) {
    problemManager.ReadRestartOverwrite( restartFileName );
  }

  std::cout << std::endl << "Running simulation:" << std::endl;

  GEOS_MARK_BEGIN("RunSimulation");
  gettimeofday(&tim, nullptr);
  t_initialize = tim.tv_sec + (tim.tv_usec / 1000000.0);

  problemManager.RunSimulation();
  gettimeofday(&tim, nullptr);
  t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

  GEOS_MARK_END("RunSimulation");

  gettimeofday(&tim, nullptr);
  t_run = tim.tv_sec + (tim.tv_usec / 1000000.0);

  printf("Done!\n\nScaling Data: initTime = %1.2fs, runTime = %1.2fs\n", t_initialize - t_start,  t_run - t_initialize );

  problemManager.ClosePythonInterpreter();

#ifdef USE_ATK
  slic::finalize();
#endif


#ifdef USE_MPI
  MPI_Finalize();
#endif


  return 0;
}
