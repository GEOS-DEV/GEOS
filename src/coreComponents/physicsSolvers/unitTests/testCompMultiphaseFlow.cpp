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
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/EventManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/FiniteVolume/CompositionalMultiphaseFlow.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::systemSolverInterface;
using namespace geosx::constitutive;

namespace
{
int global_argc;
char** global_argv;
}

template<typename T, int NDIM>
using array = multidimensionalArray::ManagedArray<T,NDIM,localIndex>;

// helper struct to represent a var and its derivatives (always with array views, not pointers)
template<int DIM>
struct TestCompositionalVarContainer
{
  array_view<real64,DIM>   value; // variable value
  array_view<real64,DIM>   dPres; // derivative w.r.t. pressure
  array_view<real64,DIM+1> dComp; // derivative w.r.t. composition
};

::testing::AssertionResult checkRelativeErrorFormat( const char *, const char *, const char *,
                                                     real64 v1, real64 v2, real64 relTol )
{
  real64 const delta = std::fabs( v1 - v2 );
  real64 const value = std::fmax( std::fabs(v1), std::fabs(v2) );
  if (delta > relTol * value)
  {
    return ::testing::AssertionFailure() << std::scientific << std::setprecision(5)
                                         << " relative error: " << delta / value
                                         << " (" << v1 << " vs " << v2 << "),"
                                         << " exceeds " << relTol << std::endl;
  }
  return ::testing::AssertionSuccess();
}

void checkRelativeError( real64 v1, real64 v2, real64 relTol )
{
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

void checkRelativeError( real64 v1, real64 v2, real64 relTol, string const & name )
{
  SCOPED_TRACE(name);
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

void checkDerivative( real64 const & valueEps, real64 const & value, real64 const & deriv,
                      real64 eps, real64 relTol, string const & name, string const & var )
{
  real64 const numDeriv = (valueEps - value) / eps;
  checkRelativeError( numDeriv, deriv, relTol, "d(" + name + ")/d(" + var + ")" );
}

template<int DIM, typename ... Args>
typename std::enable_if<DIM != 0, void>::type
checkDerivative( array_view<real64,DIM> const & valueEps,
                 array_view<real64,DIM> const & value,
                 array_view<real64,DIM> const & deriv,
                 real64 eps, real64 relTol,
                 string const & name, string const & var,
                 string_array const & labels,
                 Args ... label_lists)
{
  const auto size = valueEps.size(0);

  for (localIndex i = 0; i < size; ++i)
  {
    checkDerivative( valueEps.slice(i), value.slice(i), deriv.slice(i), eps, relTol,
                     name + "[" + labels[i] + "]", var, label_lists... );
  }
}

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
template<int DIM>
array<real64,DIM> invertLayout( array_view<real64,DIM> const & input )
{
  array<real64,DIM> output(input);
  return output;
}

template<>
array<real64,2> invertLayout( array_view<real64,2> const & input )
{
  array<real64,2> output(input.size(1), input.size(0));

  for (int i = 0; i < input.size(0); ++i)
    for (int j = 0; j < input.size(1); ++j)
      output[j][i] = input[i][j];

  return output;
}

template<>
array<real64,3> invertLayout( array_view<real64,3> const & input )
{
  array<real64,3> output(input.size(2), input.size(0), input.size(1));

  for (int i = 0; i < input.size(0); ++i)
    for (int j = 0; j < input.size(1); ++j)
      for (int k = 0; k < input.size(2); ++k)
        output[k][i][j] = input[i][j][k];

  return output;
}

void compareMatrixRows( int rowNumber, double relTol,
                        int numRowEntries1, int * indices1, double * values1,
                        int numRowEntries2, int * indices2, double * values2 )
{
  SCOPED_TRACE("Row " + std::to_string(rowNumber));

  EXPECT_EQ( numRowEntries1, numRowEntries2 );

  for (int j1 = 0, j2 = 0; j1 < numRowEntries1 && j2 < numRowEntries2; ++j1, ++j2)
  {
    while (j1 < numRowEntries1 && j2 < numRowEntries2 && indices1[j1] != indices1[j2])
    {
      while (j1 < numRowEntries1 && indices1[j1] < indices2[j2])
      {
        ADD_FAILURE() << "column " << indices1[j1] << ") in matrix 1 does not have a match";
      }
      while (j2 < numRowEntries2 && indices2[j2] < indices1[j1])
      {
        ADD_FAILURE() << "column " << indices2[j2] << ") in matrix 2 does not have a match";
      }
    }
    if (j1 < numRowEntries1 && j2 < numRowEntries2)
    {
      SCOPED_TRACE("Column " + std::to_string(indices1[j1]) );

      checkRelativeError( values1[j1], values2[j1], relTol );
    }
  }
}

void compareMatrices( Epetra_FECrsMatrix const * matrix1,
                      Epetra_FECrsMatrix const * matrix2,
                      real64 relTol )
{
  int numLocalRows1 = matrix1->NumMyRows();
  int numLocalRows2 = matrix2->NumMyRows();

  ASSERT_EQ(numLocalRows1, numLocalRows2);

  int nnz1, nnz2;
  int * indices1 = nullptr;
  int * indices2 = nullptr;
  double * row1 = nullptr;
  double * row2 = nullptr;

  // check the accuracy across local rows
  for (int i = 0; i < numLocalRows1; ++i)
  {
    matrix1->ExtractMyRowView( i, nnz1, row1, indices1 );
    matrix2->ExtractMyRowView( i, nnz2, row2, indices2 );

    compareMatrixRows( i, relTol, nnz1, indices1, row1, nnz2, indices2, row2 );
  }
}

void testNumericalJacobian( CompositionalMultiphaseFlow * solver,
                            DomainPartition * domain,
                            EpetraBlockSystem * blockSystem,
                            real64 const time_n,
                            real64 const dt,
                            double perturbParameter,
                            double relTol )
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

  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    for (localIndex esr = 0; esr < elemRegion->numSubRegions(); ++esr)
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(esr);

      for (localIndex ei = 0; ei < subRegion->size(); ++ei)
      {
        if (elemGhostRank[er][esr][ei] >= 0)
          continue;

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
          solver->UpdateStateAll(domain);
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
          solver->UpdateStateAll(domain);
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
      }
    }
  }

  jacobianFD->GlobalAssemble(true);

  // assemble the analytical jacobian
  solver->ResetStateToBeginningOfStep( domain );
  solver->AssembleSystem( domain, blockSystem, time_n, dt );

  compareMatrices( jacobian, jacobianFD.get(), relTol );

  if (::testing::Test::HasFatalFailure() || ::testing::Test::HasNonfatalFailure())
  {
    jacobian->Print(std::cout);
    jacobianFD->Print(std::cout);
  }
}

void testCompositionNumericalDerivatives( CompositionalMultiphaseFlow * solver,
                                          DomainPartition * domain,
                                          real64 perturbParameter,
                                          real64 relTol )
{
  localIndex const NC = solver->numFluidComponents();

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * elemManager = mesh->getElemManager();

  ConstitutiveManager * constitutiveManager = domain->getConstitutiveManager();
  MultiFluidBase * fluid = constitutiveManager->GetGroup<MultiFluidBase>( solver->fluidIndex() );
  ASSERT_NE( fluid, nullptr );

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );

  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegion * elemRegion = elemManager->GetRegion(er);
    SCOPED_TRACE( "Region " + std::to_string(er) + " (" + elemRegion->getName() + ")" );

    for (localIndex esr = 0; esr < elemRegion->numSubRegions(); ++esr)
    {
      CellBlockSubRegion * subRegion = elemRegion->GetSubRegion(esr);
      SCOPED_TRACE( "Subregion " + std::to_string(esr) + " (" + subRegion->getName() + ")" );

      auto & compDens  = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );
      auto & dCompDens = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );
      auto & compFrac  = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString );

      auto & dCompFrac_dCompDens = subRegion->getReference<array3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::
                                                                             dGlobalCompFraction_dGlobalCompDensityString );

      // reset the solver state to zero out variable updates
      solver->ResetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d<real64> compFracOrig( subRegion->size(), NC );
      compFracOrig = compFrac;

      // update component density and check derivatives
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver->ResetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        for (localIndex ei = 0; ei < subRegion->size(); ++ei)
        {
          real64 const dRho = perturbParameter * (compDens[ei][jc] + perturbParameter);
          dCompDens[ei][jc] = dRho;
        }

        // recompute component fractions
        solver->UpdateComponentFraction( subRegion );

        // check values in each cell
        for (localIndex ei = 0; ei < subRegion->size(); ++ei)
        {
          SCOPED_TRACE( "Element " + std::to_string(ei) );

          auto dZ_dRho = invertLayout( dCompFrac_dCompDens[ei] );
          string var = "compDens[" + components[jc] + "]";

          checkDerivative( compFrac[ei], compFracOrig[ei], dZ_dRho.slice(jc), dCompDens[ei][jc], relTol,
                           "compFrac", var, components );
        }

      }
    }
  }
}


void testphaseVolumeFractionNumericalDerivatives(CompositionalMultiphaseFlow * solver,
                                                 DomainPartition * domain,
                                                 real64 perturbParameter,
                                                 real64 relTol)
{
  localIndex const NC = solver->numFluidComponents();
  localIndex const NP = solver->numFluidPhases();

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * elemManager = mesh->getElemManager();

  ConstitutiveManager * constitutiveManager = domain->getConstitutiveManager();
  MultiFluidBase * fluid = constitutiveManager->GetGroup<MultiFluidBase>( solver->fluidIndex() );
  ASSERT_NE( fluid, nullptr );

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  auto const & phases     = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );

  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegion * elemRegion = elemManager->GetRegion(er);
    SCOPED_TRACE( "Region " + std::to_string(er) + " (" + elemRegion->getName() + ")" );

    for (localIndex esr = 0; esr < elemRegion->numSubRegions(); ++esr)
    {
      CellBlockSubRegion * subRegion = elemRegion->GetSubRegion(esr);
      SCOPED_TRACE( "Subregion " + std::to_string(esr) + " (" + subRegion->getName() + ")" );

      auto & pres  = subRegion->getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );
      auto & dPres = subRegion->getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

      auto & compDens  = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );
      auto & dCompDens = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

      auto & phaseVolFrac  = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::
                                                                       phaseVolumeFractionString );
      auto & dPhaseVolFrac_dPres  = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::
                                                                              dPhaseVolumeFraction_dPressureString );
      auto & dPhaseVolFrac_dCompDens = subRegion->getReference<array3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::
                                                                                 dPhaseVolumeFraction_dGlobalCompDensityString );

      // reset the solver state to zero out variable updates
      solver->ResetStateToBeginningOfStep( domain );

      // make a copy of unperturbed values of component fractions
      array2d<real64> phaseVolFracOrig( subRegion->size(), NP );
      phaseVolFracOrig = phaseVolFrac;

      // update pressure and check derivatives
      {
        // perturb pressure in each cell
        for (localIndex ei = 0; ei < subRegion->size(); ++ei)
        {
          real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
          dPres[ei] = dP;
        }

        // recompute component fractions
        solver->UpdateState( subRegion );

        // check values in each cell
        for (localIndex ei = 0; ei < subRegion->size(); ++ei)
        {
          SCOPED_TRACE( "Element " + std::to_string(ei) );

          checkDerivative( phaseVolFrac[ei], phaseVolFracOrig[ei], dPhaseVolFrac_dPres[ei], dPres[ei], relTol,
                           "phaseVolFrac", "Pres", phases );
        }
      }

      // update component density and check derivatives
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver->ResetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        for (localIndex ei = 0; ei < subRegion->size(); ++ei)
        {
          real64 const dRho = perturbParameter * (compDens[ei][jc] + perturbParameter);
          dCompDens[ei][jc] = dRho;
        }

        // recompute component fractions
        solver->UpdateState( subRegion );

        // check values in each cell
        for (localIndex ei = 0; ei < subRegion->size(); ++ei)
        {
          SCOPED_TRACE( "Element " + std::to_string(ei) );

          auto dS_dRho = invertLayout( dPhaseVolFrac_dCompDens[ei] );
          string var = "compDens[" + components[jc] + "]";

          checkDerivative( phaseVolFrac[ei], phaseVolFracOrig[ei], dS_dRho.slice(jc), dCompDens[ei][jc], relTol,
                           "phaseVolFrac", var, phases );
        }
      }
    }
  }
}


class CompositionalMultiphaseFlowTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    char buf[2][1024];

    char const * workdir  = global_argv[1];
    char const * filename = "testCompMultiphaseFlow.xml";

    strcpy(buf[0], "-i");
    sprintf(buf[1], "%s/%s", workdir, filename);

    constexpr int argc = 3;
    char * argv[argc] = {
      global_argv[0],
      buf[0],
      buf[1]
    };

    problemManager.SetDocumentationNodes();
    problemManager.RegisterDocumentationNodes();

    problemManager.InitializePythonInterpreter();
    problemManager.ParseCommandLineInput( argc, argv );
    problemManager.ParseInputFile();

    problemManager.Initialize( &problemManager );
    problemManager.IntermediateInitializationRecursive( &problemManager );
    problemManager.ApplyInitialConditions();
    problemManager.FinalInitializationRecursive( &problemManager );

    solver = problemManager.GetPhysicsSolverManager().GetGroup<CompositionalMultiphaseFlow>( "compflow" );
  }

  static void TearDownTestCase()
  {

  }

  static ProblemManager problemManager;
  static CompositionalMultiphaseFlow * solver;

};

ProblemManager CompositionalMultiphaseFlowTest::problemManager("ProblemManager", nullptr);
CompositionalMultiphaseFlow * CompositionalMultiphaseFlowTest::solver = nullptr;


TEST_F(CompositionalMultiphaseFlowTest, jacobianNumericalCheck)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-2;

  real64 const time = 0.0;
  real64 const dt = 1.0;

  DomainPartition   * domain = problemManager.getDomainPartition();
  EpetraBlockSystem * system = solver->getLinearSystemRepository();

  solver->ImplicitStepSetup( time, dt, domain, system );

  testNumericalJacobian( solver, domain, system, time, dt, eps, tol );
}

TEST_F(CompositionalMultiphaseFlowTest, compositionDerivativesNumericalCheck)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-5;

  DomainPartition * domain = problemManager.getDomainPartition();

  testCompositionNumericalDerivatives( solver, domain, eps, tol );
}

TEST_F(CompositionalMultiphaseFlowTest, phaseVolumeFractionDerivativesNumericalCheck)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-2;

  DomainPartition * domain = problemManager.getDomainPartition();

  testphaseVolumeFractionNumericalDerivatives(solver, domain, eps, tol);
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  if (argc < 2)
  {
    std::cerr << "Usage: testCompMultiphaseFlow <path/to/xml/dir>";
    return 1;
  }

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
