/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#endif

#include "gtest/gtest.h"

#ifdef __clang__
#define __null nullptr
#endif

#include "SetSignalHandling.hpp"
#include "stackTrace.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/EventManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/SimpleSolvers/LaplaceFEM.hpp"

using namespace geosx;
using namespace geosx::systemSolverInterface;

namespace
{
int global_argc;
char** global_argv;
int mpiRank = 0;
int mpiSize = 1;
}

class LaplaceFEMTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    problemManager.InitializePythonInterpreter();
    problemManager.ParseCommandLineInput( global_argc, global_argv );
    problemManager.ParseInputFile();
    problemManager.ProblemSetup();

    solver = problemManager.GetPhysicsSolverManager().GetGroup<LaplaceFEM>( "laplace" );
  }

  static void TearDownTestCase()
  {

  }

  static ProblemManager problemManager;
  static LaplaceFEM * solver;
};

ProblemManager LaplaceFEMTest::problemManager( "Problem", nullptr );
LaplaceFEM * LaplaceFEMTest::solver = nullptr;

TEST_F(LaplaceFEMTest, laplaceSolverCheckSolution)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());

  real64 const time = 1.0;
  real64 const dt = 1.0;
  real64 const scalingFactor = 1.0;
  int const cycleNumber = 0;

  DomainPartition * domain = problemManager.getDomainPartition();

  // Create and solve the problem
  LaplaceFEM laplaceFEM( "Temperature", domain );
  EpetraBlockSystem * system = solver->getLinearSystemRepository();
  solver->ImplicitStepSetup( time, dt, domain, system );
  solver->AssembleSystem( domain, system, time, dt );
  solver->ApplyBoundaryConditions( domain, system, time, dt );
  solver->SolverStep( time, dt, cycleNumber, domain );
  solver->ApplySystemSolution( system, scalingFactor, domain );

  // Get matrix and matrix size
  Epetra_FECrsMatrix const * const matrix = system->GetMatrix( BlockIDs::dummyScalarBlock,
                                                               BlockIDs::dummyScalarBlock );
  real64 const matrixSize3 = std::pow( static_cast<real64>( matrix->NumGlobalRows64() ), 1.0/3.0 );
  real64 const tol = 4.0 * std::pow( matrixSize3, 2 ) * eps;

  // Get solution
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  real64_array & fieldVar = nodeManager->getReference<real64_array>(string("Temperature"));

  // Compute relative error
  real64 xMin = 0.0;
  real64 xMax = 0.0;
  // TODO read them on input
  real64 const vMin = 1000.0;
  real64 const vMax = 0.0;

  // Compute domain bounds (x direction)
  r1_array const & referencePosition = nodeManager->getReference<r1_array>(dataRepository::keys::referencePositionString);
  localIndex const numNodes = nodeManager->size();
  for( localIndex a = 0 ; a < numNodes ; ++a )
  {
    R1Tensor nodePosition;
    nodePosition = referencePosition[a];
    if( xMin > nodePosition[0] )
    {
      xMin = nodePosition[0];
    }
    if( xMax < nodePosition[0] )
    {
      xMax = nodePosition[0];
    }
  }

  // Compute xMax and xMin across ranks
  real64_array gather;
  CommunicationTools::allGather( xMax, gather );
  xMax = *std::max_element( gather.begin(), gather.end() );
  CommunicationTools::allGather( xMin, gather );
  xMin = *std::min_element( gather.begin(), gather.end() );

  // Compute true solution and error
  real64 const slope = ( vMax - vMin ) / ( xMax - xMin );
  real64 error = 0.0;
  real64 normSol = 0.0;
  for( localIndex a = 0 ; a < numNodes ; ++a )
  {
    R1Tensor nodePosition;
    nodePosition = referencePosition[a];
    real64 refVal = slope * ( nodePosition[0] - xMin ) + vMin;
    error += std::pow( fieldVar[a] - refVal, 2 );
    normSol += std::pow( refVal, 2 );
  }

  // Gather errors across ranks
  CommunicationTools::allGather( error, gather );
  error = 0.0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    error += gather[p];
  }
  error = std::sqrt( error );

  // Gather solution norms across ranks
  CommunicationTools::allGather( normSol, gather );
  normSol = 0.0;
  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    normSol += gather[p];
  }
  normSol = std::sqrt( normSol );

  // Compute and check relative error
  error /= normSol;
  EXPECT_NEAR( error, 0.0, tol );
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef GEOSX_USE_MPI
  MPI_Init(&argc,&argv);

  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );

  MPI_Comm_rank(MPI_COMM_GEOSX, &mpiRank);

  MPI_Comm_size(MPI_COMM_GEOSX, &mpiSize);

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

#ifdef __clang__
#pragma clang diagnostic pop
#endif
