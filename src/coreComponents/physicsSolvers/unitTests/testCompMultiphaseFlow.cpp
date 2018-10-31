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

#include <constitutive/Fluid/CompositionalMultiphaseFluid.hpp>
#include "gtest/gtest.h"

#ifdef __clang__
#pragma clang diagnostic push
#define __null nullptr
#endif

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/EventManager.hpp"
#include "managers/DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/FiniteVolume/CompositionalMultiphaseFlow.hpp"

using namespace geosx;

namespace
{
int global_argc;
char** global_argv;
}

TEST(testCompMultiphaseFlow, numericalJacobian)
{
  ASSERT_GE(global_argc, 2);

  char buf[2][1024];

  char const * workdir  = global_argv[1];
  char const * filename = "testCompMultiphaseFlowNumJacobian.xml";

  strcpy(buf[0], "-i");
  sprintf(buf[1], "%s/%s", workdir, filename);

  constexpr int argc = 3;
  char * argv[argc] = {
    global_argv[0],
    buf[0],
    buf[1]
  };


  ProblemManager problemManager( "ProblemManager", nullptr );
  problemManager.SetDocumentationNodes();
  problemManager.RegisterDocumentationNodes();

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput( argc, argv );
  problemManager.ParseInputFile();

  problemManager.Initialize( &problemManager );
  problemManager.IntermediateInitializationRecursive( &problemManager );
  problemManager.ApplyInitialConditions();
  problemManager.FinalInitializationRecursive( &problemManager );

  CompositionalMultiphaseFlow * solver =
    problemManager.GetPhysicsSolverManager().GetGroup<CompositionalMultiphaseFlow>( "compflow" );

  //EventManager const * eventManager = problemManager.GetGroup<EventManager>(problemManager.groupKeys.eventManager.Key());
  //real64 const dt = eventManager->getReference<real64>(eventManager->viewKeys.dt);
  real64 const dt = 1;

  solver->ImplicitStepSetup( 0, dt, problemManager.getDomainPartition(), solver->getLinearSystemRepository() );

  bool res = solver->TestNumericalJacobian( problemManager.getDomainPartition(),
                                            solver->getLinearSystemRepository(),
                                            0.0, dt,
                                            1e-3, 1e-6 );

  EXPECT_TRUE( res );
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

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

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return result;
}
