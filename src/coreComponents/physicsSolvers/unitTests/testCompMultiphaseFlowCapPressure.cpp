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
using array = LvArray::Array<T,NDIM,localIndex>;

// helper struct to represent a var and its derivatives (always with array views, not pointers)
template<int DIM>
struct TestCompositionalVarContainer
{
  array_slice<real64,DIM>   value; // variable value
  array_slice<real64,DIM>   dPres; // derivative w.r.t. pressure
  array_slice<real64,DIM+1> dComp; // derivative w.r.t. composition
};

template<typename T>
::testing::AssertionResult checkRelativeErrorFormat( const char *, const char *, const char *,
                                                     T v1, T v2, T relTol )
{
  T const delta = std::abs( v1 - v2 );
  T const value = std::max( std::abs(v1), std::abs(v2) );
  if (delta > relTol * value)
  {
    return ::testing::AssertionFailure() << std::scientific << std::setprecision(5)
                                         << " relative error: " << delta / value
                                         << " (" << v1 << " vs " << v2 << "),"
                                         << " exceeds " << relTol << std::endl;
  }
  return ::testing::AssertionSuccess();
}

template<typename T>
void checkRelativeError( T v1, T v2, T relTol )
{
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

template<typename T>
void checkRelativeError( T v1, T v2, T relTol, string const & name )
{
  SCOPED_TRACE(name);
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

template<typename T>
void checkDerivative( T valueEps, T value, T deriv, real64 eps, real64 relTol, string const & name, string const & var )
{
  T numDeriv = (valueEps - value) / eps;
  checkRelativeError( deriv, numDeriv, relTol, "d(" + name + ")/d(" + var + ")" );
}

template<typename T, typename ... Args>
void
checkDerivative( arraySlice1d<T> const & valueEps,
                 arraySlice1d<T> const & value,
                 arraySlice1d<T> const & deriv,
                 real64 eps, real64 relTol,
                 string const & name, string const & var,
                 string_array const & labels,
                 Args ... label_lists )
{
  localIndex const size = labels.size(0);

  for (localIndex i = 0; i < size; ++i)
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol,
                     name + "[" + labels[i] + "]", var, label_lists... );
  }
}

template<typename T, int DIM, typename ... Args>
typename std::enable_if<(DIM > 1), void>::type
checkDerivative( array_slice<T,DIM> const & valueEps,
                 array_slice<T,DIM> const & value,
                 array_slice<T,DIM> const & deriv,
                 real64 eps, real64 relTol,
                 string const & name, string const & var,
                 string_array const & labels,
                 Args ... label_lists )
{
  const auto size = labels.size(0);

  for (localIndex i = 0; i < size; ++i)
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol,
                     name + "[" + labels[i] + "]", var, label_lists... );
  }
}

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
array1d<real64> invertLayout( arraySlice1d<real64 const> const & input, localIndex N )
{
  array<real64,1> output( N );
  for (int i = 0; i < N; ++i)
    output[i] = input[i];

  return output;
}

array2d<real64> invertLayout( arraySlice2d<real64 const> const & input, localIndex N1, localIndex N2 )
{
  array<real64,2> output( N2, N1 );

  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j)
      output[j][i] = input[i][j];

  return output;
}

array3d<real64> invertLayout( arraySlice3d<real64 const> const & input, localIndex N1, localIndex N2, localIndex N3 )
{
  array<real64,3> output( N3, N1, N2 );

  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < N2; ++j)
      for (int k = 0; k < N3; ++k)
        output[k][i][j] = input[i][j][k];

  return output;
}

void compareMatrixRow( int rowNumber, double relTol,
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

    compareMatrixRow( i, relTol, nnz1, indices1, row1, nnz2, indices2, row2 );
  }
}

template<typename LAMBDA>
void testNumericalJacobian( CompositionalMultiphaseFlow * solver,
                            DomainPartition * domain,
                            EpetraBlockSystem * blockSystem,
                            double perturbParameter,
                            double relTol,
                            LAMBDA && assembleFunction )
{
  localIndex const NC   = solver->numFluidComponents();
  localIndex const NDOF = solver->numDofPerCell();

  Epetra_FECrsMatrix * jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock, BlockIDs::compositionalBlock );
  Epetra_FEVector    * residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );
  Epetra_Map const   * rowMap   = blockSystem->GetRowMap( BlockIDs::compositionalBlock );

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  // assemble the analytical residual
  solver->ResetStateToBeginningOfStep( domain );
  residual->Scale( 0.0 );
  assembleFunction( solver, domain, jacobian, residual );

  // copy the analytical residual
  auto residualOrig = std::make_unique<Epetra_FEVector>( *residual );
  double* localResidualOrig = nullptr;
  residualOrig->ExtractView(&localResidualOrig, &localSizeInt);

  // create the numerical jacobian
  auto jacobianFD = std::make_unique<Epetra_FECrsMatrix>( Copy, jacobian->Graph() );
  jacobianFD->Scale( 0.0 );

  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegion * const elemRegion = elemManager->GetRegion(er);
    elemRegion->forElementSubRegionsIndex([&]( localIndex const esr, auto * const subRegion )
    {
      arrayView1d<integer> & elemGhostRank =
        subRegion-> template getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

      arrayView1d<globalIndex> & dofNumber =
        subRegion-> template getReference<array1d<globalIndex >>( CompositionalMultiphaseFlow::viewKeyStruct::blockLocalDofNumberString );

      arrayView1d<real64> & pres =
        subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

      arrayView1d<real64> & dPres =
        subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

      arrayView2d<real64> & compDens =
        subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

      arrayView2d<real64> & dCompDens =
        subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

      for (localIndex ei = 0; ei < subRegion->size(); ++ei)
      {
        if (elemGhostRank[ei] >= 0)
          continue;

        globalIndex offset = dofNumber[ei] * NDOF;

        real64 totalDensity = 0.0;
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          totalDensity += compDens[ei][ic];
        }

        {
          solver->ResetStateToBeginningOfStep(domain);

          real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
          dPres[ei] = dP;

          applyToSubRegions( domain, [&] ( ElementSubRegionBase * subRegion2 )
          {
            solver->UpdateState( subRegion2 );
          });

          residual->Scale( 0.0 );
          assembleFunction( solver, domain, jacobian, residual );

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

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          solver->ResetStateToBeginningOfStep(domain);

          real64 const dRho = perturbParameter * totalDensity;
          dCompDens[ei][jc] = dRho;

          applyToSubRegions( domain, [&] ( ElementSubRegionBase * subRegion2 )
          {
            solver->UpdateState( subRegion2 );
          });

          residual->Scale( 0.0 );
          assembleFunction( solver, domain, jacobian, residual );

          long long const dofIndex = integer_conversion<long long>(offset + jc + 1);

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
    });
  }

  jacobianFD->GlobalAssemble(true);

  // assemble the analytical jacobian
  solver->ResetStateToBeginningOfStep( domain );
  jacobian->Scale( 0.0 );
  assembleFunction( solver, domain, jacobian, residual );

  compareMatrices( jacobian, jacobianFD.get(), relTol );

#if 0
  if (::testing::Test::HasFatalFailure() || ::testing::Test::HasNonfatalFailure())
  {
    jacobian->Print(std::cout);
    jacobianFD->Print(std::cout);
  }
#endif
}


class CompositionalMultiphaseFlowTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    char buf[2][1024];

    char const * workdir  = global_argv[1];
    char const * filename = "testCompMultiphaseFlowBrooksCoreyCapPressure.xml";

    strcpy(buf[0], "-i");
    sprintf(buf[1], "%s/%s", workdir, filename);

    constexpr int argc = 3;
    char * argv[argc] = {
      global_argv[0],
      buf[0],
      buf[1]
    };

    problemManager.InitializePythonInterpreter();
    problemManager.ParseCommandLineInput( argc, argv );
    problemManager.ParseInputFile();

    problemManager.ProblemSetup();

    solver = problemManager.GetPhysicsSolverManager().GetGroup<CompositionalMultiphaseFlow>( "compflow" );

  }

  static void TearDownTestCase()
  {

  }

  static ProblemManager problemManager;
  static CompositionalMultiphaseFlow * solver;

};

ProblemManager CompositionalMultiphaseFlowTest::problemManager("Problem", nullptr);
CompositionalMultiphaseFlow * CompositionalMultiphaseFlowTest::solver = nullptr;

TEST_F(CompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-1; // 10% error margin

  real64 const time = 0.0;
  real64 const dt = 1e4;

  DomainPartition   * domain = problemManager.getDomainPartition();
  EpetraBlockSystem * system = solver->getLinearSystemRepository();

  solver->ImplicitStepSetup( time, dt, domain, system );

  testNumericalJacobian( solver, domain, system, eps, tol,
                         [&] ( CompositionalMultiphaseFlow * const targetSolver,
                               DomainPartition * const targetDomain,
                               Epetra_FECrsMatrix * const targetJacobian,
                               Epetra_FEVector * const targetResidual ) -> void
                         {
                           targetSolver->AssembleFluxTerms( targetDomain, targetJacobian, targetResidual, time, dt );
                         });
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

#ifdef __clang__
#pragma clang diagnostic pop
#endif
