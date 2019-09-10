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

#include "testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/FiniteVolume/CompositionalMultiphaseFlow.hpp"

using namespace geosx;
using namespace geosx::dataRepository;
using namespace geosx::constitutive;
using namespace geosx::testing;

namespace
{
int global_argc;
char** global_argv;
}

// helper struct to represent a var and its derivatives (always with array views, not pointers)
template<int DIM>
struct TestCompositionalVarContainer
{
  array_slice<real64,DIM>   value; // variable value
  array_slice<real64,DIM>   dPres; // derivative w.r.t. pressure
  array_slice<real64,DIM+1> dComp; // derivative w.r.t. composition
};

template<typename LAMBDA>
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
  real64 * localResidual = nullptr;
  residual.extractLocalVector( &localResidual );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  // assemble the analytical residual
  solver->ResetStateToBeginningOfStep( domain );
  residual.zero();
  assembleFunction( solver, domain, &jacobian, &residual, &dofManager );

  // copy the analytical residual
  ParallelVector residualOrig( residual );

  real64 * localResidualOrig = nullptr;
  residualOrig.extractLocalVector( &localResidualOrig );

  // create the numerical jacobian
  ParallelMatrix jacobianFD( jacobian );
  jacobianFD.zero();

  string const dofKey = dofManager.getKey( CompositionalMultiphaseFlow::viewKeyStruct::dofFieldString );

  solver->applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                         ElementRegionBase * const region,
                                         ElementSubRegionBase * const subRegion )
  {
    arrayView1d<integer> & elemGhostRank =
      subRegion-> template getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d<globalIndex> & dofNumber =
      subRegion-> template getReference<array1d<globalIndex >>( dofKey );

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

      globalIndex const dofIndex = dofNumber[ei];

      real64 totalDensity = 0.0;
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        totalDensity += compDens[ei][ic];
      }

      {
        solver->ResetStateToBeginningOfStep(domain);

        real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
        dPres[ei] = dP;

        solver->applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion2 )
        {
          solver->UpdateState( subRegion2 );
        });

        residual.zero();
        assembleFunction( solver, domain, &jacobian, &residual, &dofManager );

        for (localIndex lid = 0; lid < residual.localSize(); ++lid)
        {
          real64 dRdP = (localResidual[lid] - localResidualOrig[lid]) / dP;
          if (std::fabs(dRdP) > 0.0)
          {
            globalIndex gid = residual.getGlobalRowID( lid );
            jacobianFD.set( gid, dofIndex, dRdP );
          }
        }
      }

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        solver->ResetStateToBeginningOfStep(domain);

        real64 const dRho = perturbParameter * totalDensity;
        dCompDens[ei][jc] = dRho;

        solver->applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion2 )
        {
          solver->UpdateState( subRegion2 );
        });

        residual.zero();
        assembleFunction( solver, domain, &jacobian, &residual, &dofManager );

        for (localIndex lid = 0; lid < residual.localSize(); ++lid)
        {
          real64 dRdRho = (localResidual[lid] - localResidualOrig[lid]) / dRho;
          if (std::fabs(dRdRho) > 0.0)
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
  jacobian.zero();
  assembleFunction( solver, domain, &jacobian, &residual, &dofManager );

  compareMatrices( jacobian, jacobianFD, relTol );

#if 0
  if (::testing::Test::HasFatalFailure() || ::testing::Test::HasNonfatalFailure())
  {
    jacobian->Print(std::cout);
    jacobianFD->Print(std::cout);
  }
#endif
}

void testCompositionNumericalDerivatives( CompositionalMultiphaseFlow * solver,
                                          DomainPartition * domain,
                                          real64 perturbParameter,
                                          real64 relTol )
{
  localIndex const NC = solver->numFluidComponents();

  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();
  MultiFluidBase * const fluid = constitutiveManager->GetGroup<MultiFluidBase>( solver->fluidIndex() );
  ASSERT_NE( fluid, nullptr );

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );

  solver->applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                         ElementRegionBase * const region,
                                         ElementSubRegionBase * const subRegion )
  {
    SCOPED_TRACE( region->getName() + "/" + subRegion->getName() );

    arrayView2d<real64> & compDens  =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d<real64> & dCompDens =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d<real64> & compFrac  =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString );

    arrayView3d<real64> & dCompFrac_dCompDens =
      subRegion-> template getReference<array3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver->ResetStateToBeginningOfStep( domain );

    // make a copy of unperturbed values of component fractions
    array2d<real64> compFracOrig( subRegion->size(), NC );
    for (localIndex ei = 0; ei < subRegion->size(); ++ei)
    {
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        compFracOrig[ei][ic] = compFrac[ei][ic];
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
      solver->UpdateComponentFraction( subRegion );

      // check values in each cell
      for (localIndex ei = 0; ei < subRegion->size(); ++ei)
      {
        SCOPED_TRACE( "Element " + std::to_string(ei) );

        auto dZ_dRho = invertLayout( dCompFrac_dCompDens[ei], NC, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( compFrac[ei], compFracOrig[ei], dZ_dRho[jc], dCompDens[ei][jc], relTol,
                         "compFrac", var, components );
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

  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();
  MultiFluidBase * const fluid = constitutiveManager->GetGroup<MultiFluidBase>( solver->fluidIndex() );
  ASSERT_NE( fluid, nullptr );

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  auto const & phases     = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );

  solver->applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                         ElementRegionBase * const region,
                                         ElementSubRegionBase * const subRegion )
  {
    SCOPED_TRACE( region->getName() + "/" + subRegion->getName() );

    arrayView1d<real64> & pres =
      subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d<real64> & dPres =
      subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d<real64> & compDens =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d<real64> & dCompDens =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d<real64> & phaseVolFrac =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::phaseVolumeFractionString );

    arrayView2d<real64> & dPhaseVolFrac_dPres =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    arrayView3d<real64> & dPhaseVolFrac_dCompDens =
      subRegion-> template getReference<array3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver->ResetStateToBeginningOfStep( domain );

    // make a copy of unperturbed values of component fractions
    array2d<real64> phaseVolFracOrig( subRegion->size(), NP );
    for (localIndex ei = 0; ei < subRegion->size(); ++ei)
    {
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        phaseVolFracOrig[ei][ip] = phaseVolFrac[ei][ip];
      }
    }

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

        auto dS_dRho = invertLayout( dPhaseVolFrac_dCompDens[ei], NP, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( phaseVolFrac[ei], phaseVolFracOrig[ei], dS_dRho[jc], dCompDens[ei][jc], relTol,
                         "phaseVolFrac", var, phases );
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

  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();
  MultiFluidBase * const fluid = constitutiveManager->GetGroup<MultiFluidBase>( solver->fluidIndex() );
  ASSERT_NE( fluid, nullptr );

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );
  auto const & phases     = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::phaseNamesString );

  solver->applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                         ElementRegionBase * const region,
                                         ElementSubRegionBase * const subRegion )
  {
    SCOPED_TRACE( region->getName() + "/" + subRegion->getName() );

    arrayView1d<real64> & pres =
      subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

    arrayView1d<real64> & dPres =
      subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

    arrayView2d<real64> & compDens =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

    arrayView2d<real64> & dCompDens =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d<real64> & phaseMob =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::phaseMobilityString );

    arrayView2d<real64> & dPhaseMob_dPres =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dPressureString );

    arrayView3d<real64> & dPhaseMob_dCompDens =
      subRegion-> template getReference<array3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver->ResetStateToBeginningOfStep( domain );

    // make a copy of unperturbed values of component fractions
    array2d<real64> phaseVolFracOrig( subRegion->size(), NP );
    for (localIndex ei = 0; ei < subRegion->size(); ++ei)
    {
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        phaseVolFracOrig[ei][ip] = phaseMob[ei][ip];
      }
    }

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

        checkDerivative( phaseMob[ei], phaseVolFracOrig[ei], dPhaseMob_dPres[ei], dPres[ei], relTol,
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

        auto dS_dRho = invertLayout( dPhaseMob_dCompDens[ei], NP, NC );
        string var = "compDens[" + components[jc] + "]";

        checkDerivative( phaseMob[ei], phaseVolFracOrig[ei], dS_dRho[jc], dCompDens[ei][jc], relTol,
                         "phaseMob", var, phases );
      }
    }
  } );
}

class CompositionalMultiphaseFlowTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    problemManager = new ProblemManager("Problem", nullptr);
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

    problemManager->InitializePythonInterpreter();
    problemManager->ParseCommandLineInput( argc, argv );
    problemManager->ParseInputFile();

    problemManager->ProblemSetup();

    solver = problemManager->GetPhysicsSolverManager().GetGroup<CompositionalMultiphaseFlow>( "compflow" );
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

TEST_F(CompositionalMultiphaseFlowTest, derivativeNumericalCheck_composition)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-4;

  DomainPartition * domain = problemManager->getDomainPartition();

  testCompositionNumericalDerivatives( solver, domain, eps, tol );
}

TEST_F(CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseVolumeFraction)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition * domain = problemManager->getDomainPartition();

  testPhaseVolumeFractionNumericalDerivatives(solver, domain, eps, tol);
}

TEST_F(CompositionalMultiphaseFlowTest, derivativeNumericalCheck_phaseMobility)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 5e-2; // 5% error margin

  DomainPartition * domain = problemManager->getDomainPartition();

  testPhaseMobilityNumericalDerivatives(solver, domain, eps, tol);
}

TEST_F(CompositionalMultiphaseFlowTest, jacobianNumericalCheck_accumulation)
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
//    targetSolver->AssembleAccumulationTerms( targetDomain, targetJacobian, targetResidual, targetDofManager, time, dt );
//  });
}

TEST_F(CompositionalMultiphaseFlowTest, jacobianNumericalCheck_flux)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  real64 const tol = 1e-1; // 10% error margin

  real64 const time = 0.0;
  real64 const dt = 1e4;

  DomainPartition   * domain = problemManager->getDomainPartition();

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
  });
}


TEST_F(CompositionalMultiphaseFlowTest, jacobianNumericalCheck_volumeBalance)
{
  real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
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
  });
}

int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>( global_argc )];
  for( int i = 0 ; i < argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

  geosx::basicCleanup();

  return result;
}
