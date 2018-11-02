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
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/FiniteVolume/CompositionalMultiphaseFlow.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::systemSolverInterface;

namespace
{
int global_argc;
char** global_argv;
}

bool compareMatrices( Epetra_FECrsMatrix const * matrix1,
                      Epetra_FECrsMatrix const * matrix2,
                      real64 relTol )
{
#if 1
  matrix1->Print(std::cout);
  matrix2->Print(std::cout);
#endif

  int numLocalRows = matrix1->NumMyRows();
  int numLocalRowsFD = matrix2->NumMyRows();

  if (numLocalRows != numLocalRowsFD)
  {
    GEOS_LOG("Mismatch in number of local rows: " << numLocalRows << " != " << numLocalRowsFD);
    return false;
  }

  bool result = true;
  double * row1 = nullptr;
  double * row2 = nullptr;
  int numEntries1, numEntries2;
  int * indices1 = nullptr;
  int* indices2 = nullptr;

  // check the accuracy across local rows
  for (int i = 0; i < numLocalRows; ++i)
  {
    matrix1->ExtractMyRowView( i, numEntries1, row1, indices1 );
    matrix2->ExtractMyRowView( i, numEntries2, row2, indices2 );

    if (numEntries1 != numEntries2)
    {
      GEOS_LOG( "Mismatch in number of entries in local row " << i << ": " << numEntries1 << " != " << numEntries2);
      result = false;
    }
    for (int j1 = 0, j2 = 0; j1 < numEntries1 && j2 < numEntries2; ++j1, ++j2)
    {
      while (j1 < numEntries1 && j2 < numEntries2 && indices1[j1] != indices1[j2])
      {
        while (j1 < numEntries1 && indices1[j1] < indices2[j2])
        {
          GEOS_LOG( "Entry (" << i << ", " << indices1[j1] << ") in matrix 1 does not have a match" );
          result = false;
          j1++;
        }
        while (j2 < numEntries2 && indices2[j2] < indices1[j1])
        {
          GEOS_LOG( "Entry (" << i << ", " << indices2[j2] << ") in matrix 2 does not have a match" );
          result = false;
          j2++;
        }
      }
      if (j1 < numEntries1 && j2 < numEntries2)
      {
        double const delta = std::fabs(row1[j1] - row2[j1]);
        double const value = std::fmax(std::fabs(row1[j1]), std::fabs(row2[j2]));
        if (value > 0 && delta / value > relTol)
        {
          GEOS_LOG( "Entry (" << i << ", " << indices1[j1] << ") relative error: " << delta/value );
          result = false;
        }
      }
    }
  }

  return result;
}

bool testNumericalJacobian( CompositionalMultiphaseFlow * solver,
                            DomainPartition * domain,
                            EpetraBlockSystem * blockSystem,
                            real64 const time_n,
                            real64 const dt,
                            double perturbParameter,
                            double relTol)
{
  localIndex const NC   = solver->numFluidComponents();
  localIndex const NDOF = solver->numDofPerCell();

  Epetra_FECrsMatrix const * jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock, BlockIDs::compositionalBlock );
  Epetra_FEVector const * residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );
  Epetra_Map      const * rowMap   = blockSystem->GetRowMap( BlockIDs::compositionalBlock );

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto elemGhostRank =
    elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  auto pres      = elemManager->ConstructViewAccessor<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );
  auto dPres     = elemManager->ConstructViewAccessor<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );
  auto compDens  = elemManager->ConstructViewAccessor<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );
  auto dCompDens = elemManager->ConstructViewAccessor<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

  auto blockLocalDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>>( CompositionalMultiphaseFlow::viewKeyStruct::blockLocalDofNumberString );

  // assemble the analytical residual
  solver->ResetStateToBeginningOfStep( domain );
  solver->AssembleSystem( domain, blockSystem, time_n, dt );

  // copy the analytical residual
  auto residualOrig = std::make_unique<Epetra_FEVector>( *residual );
  double* localResidualOrig = nullptr;
  residualOrig->ExtractView(&localResidualOrig, &localSizeInt);

  // create the numerical jacobian
  auto jacobianFD = std::make_unique<Epetra_FECrsMatrix>( Copy, jacobian->Graph() );
  jacobianFD->Scale( 0.0 );

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei ) -> void
  {
    if (elemGhostRank[er][esr][ei] >= 0)
      return;

    globalIndex offset = blockLocalDofNumber[er][esr][ei] * NDOF;

    real64 totalDensity = 0.0;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      totalDensity += compDens[er][esr][ei][ic];
    }

    {
      solver->ResetStateToBeginningOfStep(domain);
      real64 const dP = perturbParameter * (pres[er][esr][ei] + perturbParameter);
      dPres[er][esr][ei] = dP;
      solver->UpdateState( domain );
      solver->AssembleSystem(domain, blockSystem, time_n, dt);
      long long const dofIndex = integer_conversion<long long>(offset);

      for (int lid = 0; lid < localSizeInt; ++lid)
      {
        real64 dRdP = (localResidual[lid] - localResidualOrig[lid]) / dP;
        if (std::fabs(dRdP) > 0.0)
        {
          long long gid = rowMap->GID64(lid);
          jacobianFD->ReplaceGlobalValues(gid, 1, &dRdP, &dofIndex);
        }
      }
    }

    for (localIndex ic = 0; ic < NC; ++ic)
    {
      solver->ResetStateToBeginningOfStep(domain);
      real64 const dRho = perturbParameter * totalDensity;
      dCompDens[er][esr][ei][ic] = dRho;
      solver->UpdateState( domain );
      solver->AssembleSystem(domain, blockSystem, time_n, dt);
      long long const dofIndex = integer_conversion<long long>(offset + ic + 1);

      for (int lid = 0; lid < localSizeInt; ++lid)
      {
        real64 dRdRho = (localResidual[lid] - localResidualOrig[lid]) / dRho;
        if (std::fabs(dRdRho) > 0.0)
        {
          long long gid = rowMap->GID64(lid);
          jacobianFD->ReplaceGlobalValues(gid, 1, &dRdRho, &dofIndex);
        }
      }
    }
  });

  jacobianFD->GlobalAssemble(true);

  // assemble the analytical jacobian
  solver->ResetStateToBeginningOfStep( domain );
  solver->AssembleSystem( domain, blockSystem, time_n, dt );

  return compareMatrices( jacobian, jacobianFD.get(), relTol );
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

  real64 const dt = 1;

  solver->ImplicitStepSetup( 0, dt, problemManager.getDomainPartition(), solver->getLinearSystemRepository() );

  auto eps = sqrt(std::numeric_limits<real64>::epsilon());

  bool res = testNumericalJacobian( solver,
                                    problemManager.getDomainPartition(),
                                    solver->getLinearSystemRepository(),
                                    0.0, dt,
                                    eps, 1e-2 );

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
