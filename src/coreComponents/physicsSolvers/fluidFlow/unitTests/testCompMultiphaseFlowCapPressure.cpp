/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlow.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

// helper struct to represent a var and its derivatives (always with array views, not pointers)
template< int DIM >
struct TestCompositionalVarContainer
{
  ArraySlice< real64, DIM >   value; // variable value
  ArraySlice< real64, DIM >   dPres; // derivative w.r.t. pressure
  ArraySlice< real64, DIM+1 > dComp; // derivative w.r.t. composition
};

template< typename LAMBDA >
void testNumericalJacobian( CompositionalMultiphaseFlow * solver,
                            DomainPartition * domain,
                            double perturbParameter,
                            double relTol,
                            LAMBDA && assembleFunction )
{
  localIndex const NC = solver->numFluidComponents();

  ParallelMatrix & jacobian = solver->getSystemMatrix();
  ParallelVector & residual = solver->getSystemRhs();
  DofManager const & dofManager = solver->getDofManager();

  // get a view into local residual vector
  real64 const * localResidual = residual.extractLocalVector();

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // assemble the analytical residual
  solver->ResetStateToBeginningOfStep( domain );
  residual.zero();
  jacobian.zero();
  residual.open();
  jacobian.open();
  assembleFunction( solver, domain, &jacobian, &residual, &dofManager );
  residual.close();
  jacobian.close();

  // copy the analytical residual
  ParallelVector residualOrig( residual );

  real64 const * localResidualOrig = residualOrig.extractLocalVector();

  // create the numerical jacobian
  ParallelMatrix jacobianFD( jacobian );
  jacobianFD.zero();
  jacobianFD.open();

  string const dofKey = dofManager.getKey( CompositionalMultiphaseFlow::viewKeyStruct::dofFieldString );

  solver->forTargetSubRegions( mesh, [&]( localIndex const,
                                          ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer > & elemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d< globalIndex > & dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 > & pres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d< real64 > & dPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d< real64 > & compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d< real64 > & dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      if( elemGhostRank[ei] >= 0 )
        continue;

      globalIndex const dofIndex = dofNumber[ei];

      real64 totalDensity = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        totalDensity += compDens[ei][ic];
      }

      {
        solver->ResetStateToBeginningOfStep( domain );

        real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
        dPres[ei] = dP;

        solver->forTargetSubRegions( mesh, [&]( localIndex const targetIndex2,
                                                ElementSubRegionBase & subRegion2 )
        {
          solver->UpdateState( subRegion2, targetIndex2 );
        } );

        residual.zero();
        jacobian.zero();
        residual.open();
        jacobian.open();
        assembleFunction( solver, domain, &jacobian, &residual, &dofManager );
        residual.close();
        jacobian.close();

        for( localIndex lid = 0; lid < residual.localSize(); ++lid )
        {
          real64 dRdP = ( localResidual[lid] - localResidualOrig[lid] ) / dP;
          if( std::fabs( dRdP ) > 0.0 )
          {
            globalIndex gid = residual.getGlobalRowID( lid );
            jacobianFD.set( gid, dofIndex, dRdP );
          }
        }
      }

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        solver->ResetStateToBeginningOfStep( domain );

        real64 const dRho = perturbParameter * totalDensity;
        dCompDens[ei][jc] = dRho;

        solver->forTargetSubRegions( mesh, [&]( localIndex const targetIndex2,
                                                ElementSubRegionBase & subRegion2 )
        {
          solver->UpdateState( subRegion2, targetIndex2 );
        } );

        residual.zero();
        jacobian.zero();
        residual.open();
        jacobian.open();
        assembleFunction( solver, domain, &jacobian, &residual, &dofManager );
        residual.close();
        jacobian.close();

        for( localIndex lid = 0; lid < residual.localSize(); ++lid )
        {
          real64 dRdRho = ( localResidual[lid] - localResidualOrig[lid] ) / dRho;
          if( std::fabs( dRdRho ) > 0.0 )
          {
            globalIndex gid = residual.getGlobalRowID( lid );
            jacobianFD.set( gid, dofIndex + jc + 1, dRdRho );
          }
        }
      }
    }
  } );

  jacobianFD.close();

  // assemble the analytical jacobian
  solver->ResetStateToBeginningOfStep( domain );
  residual.zero();
  jacobian.zero();
  residual.open();
  jacobian.open();
  assembleFunction( solver, domain, &jacobian, &residual, &dofManager );
  residual.close();
  jacobian.close();

  compareMatrices( jacobian, jacobianFD, relTol );

#if 0
  if( ::testing::Test::HasFatalFailure() || ::testing::Test::HasNonfatalFailure())
  {
    jacobian->Print( std::cout );
    jacobianFD->Print( std::cout );
  }
#endif
}


class CompositionalMultiphaseFlowTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    problemManager = new ProblemManager( "Problem", nullptr );

    problemManager->InitializePythonInterpreter();
    problemManager->ParseCommandLineInput();
    problemManager->ParseInputFile();

    problemManager->ProblemSetup();

    solver = problemManager->GetPhysicsSolverManager().GetGroup< CompositionalMultiphaseFlow >( "compflow" );
  }

  static void TearDownTestCase()
  {
    delete problemManager;
    problemManager = nullptr;
    solver = nullptr;
  }

  static ProblemManager * problemManager;
  static CompositionalMultiphaseFlow * solver;

};

ProblemManager * CompositionalMultiphaseFlowTest::problemManager = nullptr;
CompositionalMultiphaseFlow * CompositionalMultiphaseFlowTest::solver = nullptr;

TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 1e-1; // 10% error margin

  real64 const time = 0.0;
  real64 const dt = 1e4;

  DomainPartition * domain = problemManager->getDomainPartition();

  solver->ImplicitStepSetup( time,
                             dt,
                             domain,
                             solver->getDofManager(),
                             solver->getSystemMatrix(),
                             solver->getSystemRhs(),
                             solver->getSystemSolution() );

  testNumericalJacobian( solver, domain, eps, tol,
                         [&] ( CompositionalMultiphaseFlow * const targetSolver,
                               DomainPartition * const targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager )
  {
    targetSolver->AssembleFluxTerms( time, dt, targetDomain, targetDofManager, targetJacobian, targetResidual );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  GEOSX_ERROR_IF_NE( argc, 2 );

  std::string inputFileName = argv[ 1 ];
  inputFileName += "/testCompMultiphaseFlowBrooksCoreyCapPressure.xml";
  geosx::overrideInputFileName( inputFileName );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
