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
#include "physicsSolvers/FiniteVolume/SinglePhaseFlow_TPFA.hpp"
#include "physicsSolvers/BoundaryConditions/BoundaryConditionManager.hpp"

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

template <typename LAMBDA>
real64 computeErrorNorm(DomainPartition const *domain, std::string const &solutionFieldName,
                        R1Tensor const h, LAMBDA &&solutionFunc)
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<real64_array const>
  solutionField = elemManager->ConstructViewAccessor<real64_array>(solutionFieldName);

  ElementRegionManager::ElementViewAccessor<r1_array const>
  center = elemManager->ConstructViewAccessor<r1_array>(CellBlock::viewKeyStruct::elementCenterString);

  // compute local error norm
  real64 localErrorNorm = sumOverElemsInMesh(mesh, [&] (localIndex const er,
                                                        localIndex const esr,
                                                        localIndex const k) -> real64
  {
    real64 const val = solutionField[er][esr][k] - solutionFunc(center[er][esr][k]);
    return val * val * h[0] * h[1] * h[2]; // TODO how to compute error norm correctly?
  });

  // compute global error norm
  real64 globalErrorNorm;
  MPI_Allreduce(&localErrorNorm, &globalErrorNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  globalErrorNorm = sqrt(globalErrorNorm);

  return globalErrorNorm;
}

// extracts the function used in boundary conditions, assuming it's unique
std::pair<FunctionBase const *, real64>
getBoundaryFunctionAndScale(std::string const &fieldName)
{
  real64 scale = 1.0;
  BoundaryConditionBase const * bcFound = nullptr;
  BoundaryConditionManager const * bcMgr = BoundaryConditionManager::get();
  for (auto & subGroup : bcMgr->GetSubGroups())
  {
    BoundaryConditionBase const * bc = subGroup.second->group_cast<BoundaryConditionBase const *>();

    if (bc->initialCondition())
      continue;

    if (bc->getData<std::string>(BoundaryConditionBase::viewKeyStruct::fieldNameString) == fieldName)
    {
      bcFound = bc;
      scale = *bc->getData<real64>(BoundaryConditionBase::viewKeyStruct::scaleString);
      break;
    }
  }

  if (bcFound == nullptr)
    return { nullptr, scale };

  std::string funcName = bcFound->getData<std::string>(BoundaryConditionBase::viewKeyStruct::functionNameString);
  if (funcName.empty())
    return { nullptr, scale };

  FunctionBase const * fnFound = nullptr;
  NewFunctionManager const * fnMgr = NewFunctionManager::Instance();
  for (auto & subGroup : fnMgr->GetSubGroups())
  {
    FunctionBase const * fn = subGroup.second->group_cast<FunctionBase const *>();

    if (fn->getName() == funcName)
    {
      fnFound = fn;
      break;
    }
  }

  return { fnFound, scale };
}

void runProblem(ProblemManager & problemManager, int argc, char** argv)
{
  problemManager.SetDocumentationNodes();
  problemManager.RegisterDocumentationNodes();

  problemManager.InitializePythonInterpreter();
  problemManager.ParseCommandLineInput(argc, argv);
  problemManager.ParseInputFile();

  problemManager.Initialize(&problemManager);
  problemManager.FinalInitializationRecursive(&problemManager);
  problemManager.ApplyInitialConditions();

  std::cout << std::endl << "Running simulation:" << std::endl;
  problemManager.RunSimulation();
  std::cout << "Done!" << std::endl;

  problemManager.ClosePythonInterpreter();
}

TEST(singlePhaseFlow,analyticalTest)
{
  ProblemManager problemManager("ProblemManager", nullptr);
  runProblem(problemManager, global_argc, global_argv);

  auto const pair = getBoundaryFunctionAndScale(SinglePhaseFlow_TPFA::viewKeyStruct::fluidPressureString);
  FunctionBase const * fn = pair.first;
  real64 const scale = pair.second;

  R1Tensor h(0.1, 0.1, 1.0); // TODO how to get h from input?

  real64 const err = computeErrorNorm(problemManager.getDomainPartition(),
                                      SinglePhaseFlow_TPFA::viewKeyStruct::fluidPressureString,
                                      h,
                                      [&] (R1Tensor const & pt) -> real64
                                      {
                                        return scale * fn->Evaluate(pt.Data());
                                      });

  std::cout << "Computed error norm: " << err << std::endl;
  EXPECT_TRUE(err < 1.0); // TODO what is an appropriate error bound/estimate for FV?
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef USE_MPI
  int rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
  {
    global_argv[i] = argv[i];
    std::cout<<argv[i]<<std::endl;
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

#ifdef USE_ATK
  slic::finalize();
#endif

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return result;
}