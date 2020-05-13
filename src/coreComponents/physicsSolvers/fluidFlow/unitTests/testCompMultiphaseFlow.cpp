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
#include "constitutive/fluid/MultiFluidBase.hpp"
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

void testCompositionNumericalDerivatives( CompositionalMultiphaseFlow * solver,
                                          DomainPartition * domain,
                                          real64 perturbParameter,
                                          real64 relTol )
{
  localIndex const NC = solver->numFluidComponents();

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  solver->forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                          ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent()->getParent()->getName() + "/" + subRegion.getName() );

    string const & fluidName = solver->fluidModelNames()[targetIndex];
    Group const * const constitutiveGroup = subRegion.GetConstitutiveModels();
    MultiFluidBase const & fluid = *constitutiveGroup->GetGroup< MultiFluidBase >( fluidName );
    arrayView1d< string const > const & components = fluid.componentNames();

    arrayView2d< real64 > & compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d< real64 > & dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d< real64 > & compFrac =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString );

    arrayView3d< real64 > & dCompFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver->ResetStateToBeginningOfStep( domain );

    // make a copy of unperturbed values of component fractions
    array2d< real64 > compFracOrig( subRegion.size(), NC );
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compFracOrig[ei][ic] = compFrac[ei][ic];
      }
    }

    // update component density and check derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
      solver->ResetStateToBeginningOfStep( domain );

      // perturb a single component density in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
        dCompDens[ei][jc] = dRho;
      }

      // recompute component fractions
      solver->UpdateComponentFraction( subRegion );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        auto dZ_dRho = invertLayout( dCompFrac_dCompDens[ei], NC, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( compFrac[ei].toSliceConst(),
                         compFracOrig[ei].toSliceConst(),
                         dZ_dRho[jc].toSliceConst(),
                         dCompDens[ei][jc],
                         relTol,
                         "compFrac",
                         var,
                         components );
      }
    }
  } );
}


void testPhaseVolumeFractionNumericalDerivatives( CompositionalMultiphaseFlow * solver,
                                                  DomainPartition * domain,
                                                  real64 perturbParameter,
                                                  real64 relTol )
{
  localIndex const NC = solver->numFluidComponents();
  localIndex const NP = solver->numFluidPhases();

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  solver->forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                          ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent()->getParent()->getName() + "/" + subRegion.getName() );

    string const & fluidName = solver->fluidModelNames()[targetIndex];
    Group const * const constitutiveGroup = subRegion.GetConstitutiveModels();
    MultiFluidBase const & fluid = *constitutiveGroup->GetGroup< MultiFluidBase >( fluidName );
    arrayView1d< string const > const & components = fluid.componentNames();
    arrayView1d< string const > const & phases = fluid.phaseNames();

    arrayView1d< real64 > & pres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d< real64 > & dPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d< real64 > & compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d< real64 > & dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d< real64 > & phaseVolFrac =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::phaseVolumeFractionString );

    arrayView2d< real64 > & dPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    arrayView3d< real64 > & dPhaseVolFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver->ResetStateToBeginningOfStep( domain );

    // make a copy of unperturbed values of component fractions
    array2d< real64 > phaseVolFracOrig( subRegion.size(), NP );
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        phaseVolFracOrig[ei][ip] = phaseVolFrac[ei][ip];
      }
    }

    // update pressure and check derivatives
    {
      // perturb pressure in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
        dPres[ei] = dP;
      }

      // recompute component fractions
      solver->UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        checkDerivative( phaseVolFrac[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dPhaseVolFrac_dPres[ei].toSliceConst(),
                         dPres[ei],
                         relTol,
                         "phaseVolFrac",
                         "Pres",
                         phases );
      }
    }

    // update component density and check derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
      solver->ResetStateToBeginningOfStep( domain );

      // perturb a single component density in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
        dCompDens[ei][jc] = dRho;
      }

      // recompute component fractions
      solver->UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        auto dS_dRho = invertLayout( dPhaseVolFrac_dCompDens[ei], NP, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( phaseVolFrac[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dS_dRho[jc].toSliceConst(),
                         dCompDens[ei][jc],
                         relTol,
                         "phaseVolFrac",
                         var,
                         phases );
      }
    }
  } );
}

void testPhaseMobilityNumericalDerivatives( CompositionalMultiphaseFlow * solver,
                                            DomainPartition * domain,
                                            real64 perturbParameter,
                                            real64 relTol )
{
  localIndex const NC = solver->numFluidComponents();
  localIndex const NP = solver->numFluidPhases();

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  solver->forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                          ElementSubRegionBase & subRegion )
  {
    SCOPED_TRACE( subRegion.getParent()->getName() + "/" + subRegion.getName() );

    string const & fluidName = solver->fluidModelNames()[targetIndex];
    Group const * const constitutiveGroup = subRegion.GetConstitutiveModels();
    MultiFluidBase const & fluid = *constitutiveGroup->GetGroup< MultiFluidBase >( fluidName );
    arrayView1d< string const > const & components = fluid.componentNames();
    arrayView1d< string const > const & phases = fluid.phaseNames();

    arrayView1d< real64 > & pres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d< real64 > & dPres =
      subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d< real64 > & compDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d< real64 > & dCompDens =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d< real64 > & phaseMob =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::phaseMobilityString );

    arrayView2d< real64 > & dPhaseMob_dPres =
      subRegion.getReference< array2d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dPressureString );

    arrayView3d< real64 > & dPhaseMob_dCompDens =
      subRegion.getReference< array3d< real64 > >( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver->ResetStateToBeginningOfStep( domain );

    // make a copy of unperturbed values of component fractions
    array2d< real64 > phaseVolFracOrig( subRegion.size(), NP );
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        phaseVolFracOrig[ei][ip] = phaseMob[ei][ip];
      }
    }

    // update pressure and check derivatives
    {
      // perturb pressure in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dP = perturbParameter * ( pres[ei] + perturbParameter );
        dPres[ei] = dP;
      }

      // recompute component fractions
      solver->UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        checkDerivative( phaseMob[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dPhaseMob_dPres[ei].toSliceConst(),
                         dPres[ei],
                         relTol,
                         "phaseVolFrac",
                         "Pres",
                         phases );
      }
    }

    // update component density and check derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
      solver->ResetStateToBeginningOfStep( domain );

      // perturb a single component density in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        real64 const dRho = perturbParameter * ( compDens[ei][jc] + perturbParameter );
        dCompDens[ei][jc] = dRho;
      }

      // recompute component fractions
      solver->UpdateState( subRegion, targetIndex );

      // check values in each cell
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        SCOPED_TRACE( "Element " + std::to_string( ei ) );

        auto dS_dRho = invertLayout( dPhaseMob_dCompDens[ei], NP, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( phaseMob[ei].toSliceConst(),
                         phaseVolFracOrig[ei].toSliceConst(),
                         dS_dRho[jc].toSliceConst(),
                         dCompDens[ei][jc],
                         relTol,
                         "phaseMob",
                         var,
                         phases );
      }
    }
  } );
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

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_composition )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 1e-4;

  DomainPartition * domain = problemManager->getDomainPartition();

  testCompositionNumericalDerivatives( solver, domain, eps, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseVolumeFraction )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition * domain = problemManager->getDomainPartition();

  testPhaseVolumeFractionNumericalDerivatives( solver, domain, eps, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseMobility )
{
  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon());
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition * domain = problemManager->getDomainPartition();

  testPhaseMobilityNumericalDerivatives( solver, domain, eps, tol );
}

TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_accumulation )
{
//  real64 const eps = sqrt( std::numeric_limits<real64>::epsilon() );
//  real64 const tol = 1e-1; // 10% error margin

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

//  testNumericalJacobian( solver, domain, eps, tol,
//                         [&] ( CompositionalMultiphaseFlow * targetSolver,
//                               DomainPartition * targetDomain,
//                               ParallelMatrix * targetJacobian,
//                               ParallelVector * targetResidual,
//                               DofManager const * targetDofManager )
//  {
//    targetSolver->AssembleAccumulationTerms( targetDomain, targetJacobian, targetResidual, targetDofManager, time, dt
// );
//  });
}

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


TEST_F( CompositionalMultiphaseFlowTest, jacobianNumericalCheck_volumeBalance )
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
                         [&] ( CompositionalMultiphaseFlow * targetSolver,
                               DomainPartition * targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager )
  {
    targetSolver->AssembleVolumeBalanceTerms( time, dt, targetDomain, targetDofManager, targetJacobian, targetResidual );
  } );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  GEOSX_ERROR_IF_NE( argc, 2 );

  std::string inputFileName = argv[ 1 ];
  inputFileName += "/testCompMultiphaseFlow.xml";
  geosx::overrideInputFileName( inputFileName );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
