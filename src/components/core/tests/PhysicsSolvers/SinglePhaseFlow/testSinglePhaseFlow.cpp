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

#include "gtest/gtest.h"

#include "common/Logger.hpp"
#include "common/TimingMacros.hpp"
#include <cmath>
#include <mpi.h>
#include <iostream>
#include <sys/time.h>
#include "dataRepository/SidreWrapper.hpp"
#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "managers/ProblemManager.hpp"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace geosx;
using namespace dataRepository;
#ifdef USE_ATK
using namespace axom;
#endif

namespace
{
int global_argc;
char** global_argv;
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
  {
    global_argv[i] = argv[i];
    std::cout<<argv[i]<<std::endl;
  }

  return RUN_ALL_TESTS();
}

TEST(singlePhaseFlow,analyticalTest)
{
#ifdef USE_MPI
  int rank;
  MPI_Init(&global_argc,&global_argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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

  ProblemManager problemManager( "ProblemManager", nullptr );
  problemManager.SetDocumentationNodes();
  problemManager.RegisterDocumentationNodes();  

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( global_argc, global_argv );
  problemManager.ParseInputFile();

  problemManager.Initialize( &problemManager );
  problemManager.FinalInitializationRecursive( &problemManager );
  problemManager.ApplyInitialConditions();

  std::cout << std::endl << "Running simulation:" << std::endl;
  problemManager.RunSimulation();
  std::cout << "Done!";

  problemManager.ClosePythonInterpreter();

#ifdef USE_ATK
  slic::finalize();
#endif

#ifdef USE_MPI
  MPI_Finalize();
#endif

  // TODO: check solution vs analytical
  EXPECT_TRUE(true);
}
