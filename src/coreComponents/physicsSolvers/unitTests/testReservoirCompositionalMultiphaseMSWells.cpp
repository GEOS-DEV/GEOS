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
#include "wells/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/CoupledSolvers/ReservoirSolver.hpp"
#include "physicsSolvers/Wells/CompositionalMultiphaseWell.hpp"
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
struct TestReservoirVarContainer
{
  array_slice<real64,DIM> value; // variable value
  array_slice<real64,DIM> dPres; // derivative w.r.t. pressure
  array_slice<real64,DIM> dRate; // derivative w.r.t. rate
};


void testMixtureDensityNumericalDerivatives( CompositionalMultiphaseWell * solver,
                                             DomainPartition * domain,
                                             real64 perturbParameter,
                                             real64 relTol )
{
  localIndex const NC = solver->NumFluidComponents();

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * elemManager = mesh->getElemManager();

  ConstitutiveManager * constitutiveManager = domain->getConstitutiveManager();
  CompositionalMultiphaseFlow * flowSolver = solver->getParent()->GetGroup("compositionalMultiphaseFlow")->group_cast<CompositionalMultiphaseFlow*>();

  MultiFluidBase * fluid = constitutiveManager->GetGroup<MultiFluidBase>( flowSolver->fluidIndex() );
  ASSERT_NE( fluid, nullptr );

  auto const & components = fluid->getReference<string_array>( MultiFluidBase::viewKeyStruct::componentNamesString );

  // bind the stored reservoir views to the current domain
  //solver->ImplicitStepSetup( 0, 0, domain, nullptr );
  
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {

    SCOPED_TRACE( "Well " + subRegion->getName() );

    arrayView1d<real64> & pres =
      subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::pressureString );

    arrayView1d<real64> & dPres =
      subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::deltaPressureString );

    arrayView2d<real64> & compDens =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::globalCompDensityString );

    arrayView2d<real64> & dCompDens =
      subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64> const & wellElemMixtureDensity =
      subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::mixtureDensityString );
  
    arrayView1d<real64> const & dWellElemMixtureDensity_dPres =
      subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dMixtureDensity_dPressureString );

    arrayView2d<real64> const & dWellElemMixtureDensity_dCompDens =
      subRegion->getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

    // reset the solver state to zero out variable updates
    solver->ResetStateToBeginningOfStep( domain );

    // make a copy of unperturbed values of component fractions
    array1d<real64> wellElemMixtureDensityOrig( subRegion->size() );
    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      wellElemMixtureDensityOrig[iwelem] = wellElemMixtureDensity[iwelem];
    }

    // update pressure and check derivatives
    {
      // perturb pressure in each cell
      for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
      {
        real64 const dP = perturbParameter * (pres[iwelem] + perturbParameter);
        dPres[iwelem] = dP;
      }

      // recompute component fractions
      solver->UpdateState( subRegion );

      // check values in each cell
      for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
      {
        SCOPED_TRACE( "Element " + std::to_string(iwelem) );

        checkDerivative( wellElemMixtureDensity[iwelem], wellElemMixtureDensityOrig[iwelem], 
                         dWellElemMixtureDensity_dPres[iwelem], dPres[iwelem], relTol,
                           "wellElemMixtureDensity", "Pres" );
      }

      // update component density and check derivatives
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        // reset the solver state to zero out variable updates (resetting the whole domain is overkill...)
        solver->ResetStateToBeginningOfStep( domain );

        // perturb a single component density in each cell
        for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
        {
          real64 const dRho = perturbParameter * (compDens[iwelem][jc] + perturbParameter);
          dCompDens[iwelem][jc] = dRho;
        }

        // recompute component fractions
        solver->UpdateState( subRegion );

        // check values in each cell
        for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
        {
          SCOPED_TRACE( "Element " + std::to_string(iwelem) );

          string var = "compDens[" + components[jc] + "]";

          checkDerivative( wellElemMixtureDensity[iwelem], wellElemMixtureDensityOrig[iwelem], 
                           dWellElemMixtureDensity_dCompDens[iwelem][jc], dCompDens[iwelem][jc], 
                           relTol, "wellElemMixtureDensity", var );
        }
      }
    }
  });
}


template<typename LAMBDA>
void testNumericalJacobian( ReservoirSolver * solver,
                            DomainPartition * domain,
                            double perturbParameter,
                            double relTol,
                            LAMBDA && assembleFunction )
{
  CompositionalMultiphaseWell * wellSolver = solver->GetWellSolver()->group_cast<CompositionalMultiphaseWell*>();
  CompositionalMultiphaseFlow * flowSolver = solver->GetFlowSolver()->group_cast<CompositionalMultiphaseFlow*>();

  localIndex const NC = flowSolver->numFluidComponents();
  
  ParallelMatrix & jacobian = solver->getSystemMatrix();
  ParallelVector & residual = solver->getSystemRhs();
  DofManager const & dofManager = solver->getDofManager();

  // get a view into local residual vector
  double* localResidual = nullptr;
  residual.extractLocalVector( &localResidual );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  // assemble the analytical residual
  solver->ResetStateToBeginningOfStep( domain );
  residual.zero();
  assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

  // copy the analytical residual
  ParallelVector residualOrig( residual );
  double* localResidualOrig = nullptr;
  residualOrig.extractLocalVector( &localResidualOrig );

  // create the numerical jacobian
  ParallelMatrix jacobianFD( jacobian );
  jacobianFD.zero();

  string const resDofKey  = dofManager.getKey( wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager.getKey( wellSolver->WellElementDofName() );

  // at this point we start assembling the finite-difference block by block

  ////////////////////////////////////////////////
  // Step 1) Compute the terms in J_RR and J_WR //
  ////////////////////////////////////////////////
  
  for (localIndex er = 0; er < elemManager->numRegions(); ++er)
  {
    ElementRegionBase * const elemRegion = elemManager->GetRegion(er);
    elemRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, auto * const subRegion )
    {
      // get the degrees of freedom and ghosting information
      arrayView1d<globalIndex> & dofNumber =
        subRegion-> template getReference<array1d<globalIndex >>( resDofKey );

      arrayView1d<integer> & elemGhostRank =
        subRegion-> template getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

      // get the primary variables on the reservoir elements
      arrayView1d<real64> & pres =
        subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

      arrayView1d<real64> & dPres =
        subRegion-> template getReference<array1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

      arrayView2d<real64> & compDens =
        subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

      arrayView2d<real64> & dCompDens =
        subRegion-> template getReference<array2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );
      
      // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei 
      for (localIndex ei = 0; ei < subRegion->size(); ++ei)
      {
        if (elemGhostRank[ei] >= 0)
        {
          continue;
        }

        globalIndex const eiOffset = dofNumber[ei];

        real64 totalDensity = 0.0;
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          totalDensity += compDens[ei][ic];
        }
        
        {
          solver->ResetStateToBeginningOfStep(domain);

          // here is the perturbation in the pressure of the element
          real64 const dP = perturbParameter * (pres[ei] + perturbParameter);
          dPres[ei] = dP;
          // after perturbing, update the pressure-dependent quantities in the reservoir
          flowSolver->applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion2 )
          {
            flowSolver->UpdateState( subRegion2 );
          });

          residual.zero();
          assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

          globalIndex const dofIndex = eiOffset;

          // consider mass balance eq lid in RESERVOIR elems and WELL elems
          // this is computing J_RR and J_RW
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
 
          flowSolver->applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion2 )
          {
            flowSolver->UpdateState( subRegion2 );
          });

          residual.zero();
          assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

          globalIndex const dofIndex = eiOffset + jc + 1;

          for (localIndex lid = 0; lid < residual.localSize(); ++lid)
          {
            // here is the perturbation in the density of the element
            real64 dRdRho = (localResidual[lid] - localResidualOrig[lid]) / dRho;
            if (std::fabs(dRdRho) > 0.0)
            {
              globalIndex gid = residual.getGlobalRowID( lid );
              jacobianFD.set( gid, dofIndex, dRdRho );
            }
          }
        }
      }
    });
  }

  /////////////////////////////////////////////////
  // Step 2) Compute the terms in J_RW and J_WW //
  /////////////////////////////////////////////////      

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    
    // get the degrees of freedom, ghosting info and next well elem index
    array1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get the primary variables on the well elements
    arrayView1d<real64> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::pressureString );

    arrayView1d<real64> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::deltaPressureString );

    arrayView2d<real64> const & wellElemCompDens =
      subRegion->getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::globalCompDensityString );

    arrayView2d<real64> const & dWellElemCompDens =
      subRegion->getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::deltaGlobalCompDensityString );
    
    arrayView1d<real64> const & connRate  =
      subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::mixtureConnRateString );

    arrayView1d<real64> const & dConnRate =
      subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::deltaMixtureConnRateString );    
    
    // a) compute all the derivatives wrt to the pressure in WELL elem iwelem 
    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] >= 0)
      {
        continue;
      }

      globalIndex const iwelemOffset = wellElemDofNumber[iwelem];

      real64 wellElemTotalDensity = 0.0;
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        wellElemTotalDensity += wellElemCompDens[iwelem][ic];
      }
      
      {
        solver->ResetStateToBeginningOfStep(domain);
        
        // here is the perturbation in the pressure of the well element
        real64 const dP = perturbParameter * (wellElemPressure[iwelem] + perturbParameter);
        dWellElemPressure[iwelem] = dP;
        // after perturbing, update the pressure-dependent quantities in the well
        wellSolver->UpdateState( subRegion );

        residual.zero();
        assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

        globalIndex const dofIndex = iwelemOffset + CompositionalMultiphaseWell::ColOffset::DPRES;

        // consider mass balance eq lid in RESERVOIR elems and WELL elems
        //      this is computing J_RW and J_WW
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

        real64 const dRho = perturbParameter * wellElemTotalDensity;
        dWellElemCompDens[iwelem][jc] = dRho;
        wellSolver->UpdateStateAll( domain );

        residual.zero();
        assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

        globalIndex const dofIndex = iwelemOffset + CompositionalMultiphaseWell::ColOffset::DCOMP + jc;

        for (localIndex lid = 0; lid < residual.localSize(); ++lid)
        {
          // here is the perturbation in the density of the element
          real64 dRdRho = (localResidual[lid] - localResidualOrig[lid]) / dRho;
          if (std::fabs(dRdRho) > 0.0)
          {
            globalIndex gid = residual.getGlobalRowID( lid );
            jacobianFD.set( gid, dofIndex, dRdRho );
          }
        }
      }
    }

    // b) compute all the derivatives wrt to the connection in WELL elem iwelem 
    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] >= 0)
      {
        continue;
      }

      globalIndex const iwelemOffset = wellElemDofNumber[iwelem];

      {
        solver->ResetStateToBeginningOfStep(domain);
        
        // here is the perturbation in the pressure of the well element
        real64 const dRate = perturbParameter * (connRate[iwelem] + perturbParameter);
        dConnRate[iwelem] = dRate;

        residual.zero();
        assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

        globalIndex const dofIndex = iwelemOffset + CompositionalMultiphaseWell::ColOffset::DCOMP + NC;

        // consider mass balance eq lid in RESERVOIR elems and WELL elems
        //      this is computing J_RW and J_WW
        for (localIndex lid = 0; lid < residual.localSize(); ++lid)
        {
          real64 dRdRate = (localResidual[lid] - localResidualOrig[lid]) / dRate;
          if (std::fabs(dRdRate) > 0.0)
          {
            globalIndex gid = residual.getGlobalRowID( lid );
            jacobianFD.set( gid, dofIndex, dRdRate );
          }
        }
      }
    }
  });
        
  jacobianFD.close();

  // assemble the analytical jacobian
  solver->ResetStateToBeginningOfStep( domain );
  jacobian.zero();
  assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

  compareMatrices( jacobian, jacobianFD, relTol );

#if 0
  if (::testing::Test::HasFatalFailure() || ::testing::Test::HasNonfatalFailure())
  {
    jacobian->Print(std::cout);
    jacobianFD->Print(std::cout);
  }
#endif
}

class ReservoirSolverTest : public ::testing::Test
{
protected:

  static void SetUpTestCase()
  {
    problemManager = new ProblemManager("Problem", nullptr);
    char buf[2][1024];

    char const * workdir  = global_argv[1];
    char const * filename = "testReservoirCompositionalMultiphaseMSWells.xml";

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

    solver = problemManager->GetPhysicsSolverManager().GetGroup<ReservoirSolver>( "reservoirSystem" );

    GEOS_ERROR_IF( solver == nullptr, "ReservoirSystem not found" );

  }

  static void TearDownTestCase()
  {
    delete problemManager;
    problemManager = nullptr;
    solver = nullptr;
  }

  static ProblemManager * problemManager;
  static ReservoirSolver * solver;

};

ProblemManager * ReservoirSolverTest::problemManager = nullptr;
ReservoirSolver * ReservoirSolverTest::solver = nullptr;

TEST_F(ReservoirSolverTest, derivativeNumericalCheck_mixtureDensity)
{
  //real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  //real64 const tol = 1e-4;

  //DomainPartition * domain = problemManager->getDomainPartition();

  //GEOS_ERROR_IF( solver == nullptr, "ReservoirSystem not found" );

  //CompositionalMultiphaseWell * wellSolver = solver->GetWellSolver()->group_cast<CompositionalMultiphaseWell*>();
  
  //testMixtureDensityNumericalDerivatives( wellSolver, domain, eps, tol );
}


TEST_F(ReservoirSolverTest, jacobianNumericalCheck_Perforation)
{
  //real64 const eps = sqrt(std::numeric_limits<real64>::epsilon());
  //real64 const tol = 1e-1; // 10% error margin

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

/*
  testNumericalJacobian( solver, domain, system, eps, tol,
                         [&] ( CompositionalMultiphaseWell * const targetSolver,
                               DomainPartition * const targetDomain,
                               Epetra_FECrsMatrix * const targetJacobian,
                               Epetra_FEVector * const targetResidual ) -> void
  {
    targetSolver->AssemblePerforationTerms( targetDomain, targetJacobian, targetResidual, time, dt );
  });
  */
}

TEST_F(ReservoirSolverTest, jacobianNumericalCheck_Flux)
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
                         [&] ( CompositionalMultiphaseWell * const targetSolver,
                               DomainPartition * const targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager ) -> void
  {
    targetSolver->AssembleFluxTerms( time, dt, targetDomain, targetDofManager, targetJacobian, targetResidual );
  });
  
}


TEST_F(ReservoirSolverTest, jacobianNumericalCheck_Control)
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
                         [&] ( CompositionalMultiphaseWell * const targetSolver,
                               DomainPartition * const targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager ) -> void
  {
    targetSolver->FormControlEquation( targetDomain, targetDofManager, targetJacobian, targetResidual );
  });
  
}

TEST_F(ReservoirSolverTest, jacobianNumericalCheck_VolumeBalance)
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
                         [&] ( CompositionalMultiphaseWell * const targetSolver,
                               DomainPartition * const targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager ) -> void
  {
    targetSolver->AssembleVolumeBalanceTerms( time, dt, targetDomain, targetDofManager, targetJacobian, targetResidual );
  });
}

TEST_F(ReservoirSolverTest, jacobianNumericalCheck_PressureRel)
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
                         [&] ( CompositionalMultiphaseWell * const targetSolver,
                               DomainPartition * const targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager ) -> void
  {
    targetSolver->FormPressureRelations( targetDomain, targetDofManager, targetJacobian, targetResidual );
  });
  
}


int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  geosx::basicSetup( argc, argv );

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>(global_argc)];
  for( int i=0 ; i<argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

  geosx::basicCleanup();

  return result;
}

