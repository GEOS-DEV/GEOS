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
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/FiniteVolume/SinglePhaseFlow.hpp"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"

using namespace geosx;
using namespace dataRepository;

#ifdef GEOSX_USE_ATK
using namespace axom;
#include "slic/GenericOutputStream.hpp"
#endif

namespace
{
int global_argc;
char** global_argv;
}

// compute approximate discrete L2 norm of vector with element center values returned by a lambda
template <typename LAMBDA>
real64 computeL2Norm(DomainPartition const *domain,
                     LAMBDA &&lambda)
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<real64_array const>
  volume = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);

  // compute local norm
  real64 localNorm = sumOverElemsInMesh(mesh, [&] (localIndex const er,
                                                   localIndex const esr,
                                                   localIndex const k) -> real64
  {
    real64 const val = lambda(er, esr, k);
    return val * val * volume[er][esr][k]; // TODO how to compute error norm correctly?
  });

  // compute global norm
  real64 globalNorm;
  MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalNorm);
}


// compute relative discrete L2 norm of the solution error
template <typename LAMBDA>
real64 computeErrorNorm(DomainPartition const *domain,
                        std::string const &solutionFieldName,
                        LAMBDA &&solutionFunc)
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<real64_array const>
  solutionField = elemManager->ConstructViewAccessor<real64_array>(solutionFieldName);

  ElementRegionManager::ElementViewAccessor<r1_array const>
  center = elemManager->ConstructViewAccessor<r1_array>(CellBlock::viewKeyStruct::elementCenterString);

  real64 const errorNorm = computeL2Norm(domain, [&] (localIndex const er,
                                                      localIndex const esr,
                                                      localIndex const k) -> real64
  {
    return solutionField[er][esr][k] - solutionFunc(center[er][esr][k]);
  });

  real64 const solutionNorm = computeL2Norm(domain, [&] (localIndex const er,
                                                         localIndex const esr,
                                                         localIndex const k) -> real64
  {
    return solutionField[er][esr][k];
  });

  return errorNorm / solutionNorm;
}

// extracts the solution function if one is defined in the input
FunctionBase const * getSolutionFunction()
{
  FunctionBase const * fnFound = nullptr;
  NewFunctionManager const * fnMgr = NewFunctionManager::Instance();
  for (auto & subGroup : fnMgr->GetSubGroups())
  {
    auto fn = subGroup.second->group_cast<FunctionBase const *>();

    if (fn->getName() == "solutionFunction")
    {
      fnFound = fn;
      break;
    }
  }
  return fnFound;
}

void runProblem(ProblemManager & problemManager, int argc, char** argv)
{
  problemManager.SetDocumentationNodes();
  problemManager.RegisterDocumentationNodes();

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput(argc, argv);
  problemManager.ParseInputFile();

  problemManager.Initialize(&problemManager);
  problemManager.ApplyInitialConditions();
  problemManager.FinalInitializationRecursive(&problemManager);

  std::cout << std::endl << "Running simulation:" << std::endl;
  problemManager.RunSimulation();
  std::cout << "Done!" << std::endl;

  problemManager.ClosePythonInterpreter();
}

TEST(singlePhaseFlow,analyticalTest)
{
  ProblemManager problemManager("ProblemManager", nullptr);
  runProblem(problemManager, global_argc, global_argv);

  FunctionBase const * fn = getSolutionFunction();
  ASSERT_TRUE(fn != nullptr);

  real64 const err = computeErrorNorm(problemManager.getDomainPartition(),
                                      SinglePhaseFlow::viewKeyStruct::pressureString,
                                      [&](R1Tensor const &pt) -> real64
                                      {
                                        return fn->Evaluate(pt.Data());
                                      });

  std::cout << "Computed error norm: " << err << std::endl;
  EXPECT_LT(err, 1e-6); // TODO what is an appropriate error bound/estimate for FV?
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef GEOSX_USE_MPI
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
#endif

#ifdef GEOSX_USE_ATK
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

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
  {
    global_argv[i] = argv[i];
    std::cout<<argv[i]<<std::endl;
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

#ifdef GEOSX_USE_ATK
  slic::finalize();
#endif

#ifdef GEOSX_USE_MPI
  MPI_Finalize();
#endif

  return result;
}

