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
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/CoupledSolvers/ReservoirSolver.hpp"
#include "physicsSolvers/Wells/SinglePhaseWell.hpp"
#include "physicsSolvers/FiniteVolume/SinglePhaseFlow.hpp"

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

template<typename LAMBDA>
void testNumericalJacobian( ReservoirSolver * solver,
                            DomainPartition * domain,
                            double perturbParameter,
                            double relTol,
                            LAMBDA && assembleFunction )
{
  SinglePhaseWell * wellSolver = solver->GetWellSolver()->group_cast<SinglePhaseWell*>();
  SinglePhaseFlow * flowSolver = solver->GetFlowSolver()->group_cast<SinglePhaseFlow*>();

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
    elemRegion->forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                        auto * const subRegion )
    {
      // get the dof numbers and ghosting information
      arrayView1d<globalIndex> & dofNumber =
        subRegion-> template getReference<array1d<globalIndex >>( resDofKey );

      arrayView1d<integer> & elemGhostRank =
        subRegion-> template getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

      // get the primary variables on reservoir elements
      arrayView1d<real64> & pres =
        subRegion-> template getReference<array1d<real64>>( SinglePhaseFlow::viewKeyStruct::pressureString );

      arrayView1d<real64> & dPres =
        subRegion-> template getReference<array1d<real64>>( SinglePhaseFlow::viewKeyStruct::deltaPressureString );

      // a) compute all the derivatives wrt to the pressure in RESERVOIR elem ei 
      for (localIndex ei = 0; ei < subRegion->size(); ++ei)
      {
        if (elemGhostRank[ei] >= 0)
        {
          continue;
        }
        
        globalIndex const eiOffset = dofNumber[ei];

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

          globalIndex const dofIndex = integer_conversion<long long>(eiOffset);

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
      }
    });
  }

  /////////////////////////////////////////////////
  // Step 2) Compute the terms in J_RW and J_WW //
  /////////////////////////////////////////////////      

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    
    // get the degrees of freedom and ghosting information
    array1d<globalIndex const> const & wellElemDofNumber =
      subRegion->getReference<array1d<globalIndex>>( wellDofKey );

    arrayView1d<integer const> const & wellElemGhostRank =
      subRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get the primary variables on well elements
    array1d<real64> const & wellElemPressure =
      subRegion->getReference<array1d<real64>>( SinglePhaseWell::viewKeyStruct::pressureString );

    array1d<real64> const & dWellElemPressure =
      subRegion->getReference<array1d<real64>>( SinglePhaseWell::viewKeyStruct::deltaPressureString );

    array1d<real64> const & connRate  =
      subRegion->getReference<array1d<real64>>( SinglePhaseWell::viewKeyStruct::connRateString );

    array1d<real64> const & dConnRate =
      subRegion->getReference<array1d<real64>>( SinglePhaseWell::viewKeyStruct::deltaConnRateString );    
    
    // a) compute all the derivatives wrt to the pressure in WELL elem iwelem 
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
        real64 const dP = perturbParameter * (wellElemPressure[iwelem] + perturbParameter);
        dWellElemPressure[iwelem] = dP;
        // after perturbing, update the pressure-dependent quantities in the well
        wellSolver->UpdateState( subRegion );

        residual.zero();
        assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

        globalIndex const dofIndex = iwelemOffset + SinglePhaseWell::ColOffset::DPRES;

        // consider mass balance eq lid in RESERVOIR elems and WELL elems
        //      this is computing J_RW and J_WW
        for (int lid = 0; lid < residual.localSize(); ++lid)
        {
          real64 dRdP = (localResidual[lid] - localResidualOrig[lid]) / dP;
          if (std::fabs(dRdP) > 0.0)
          {
            globalIndex gid = residual.getGlobalRowID( lid );
            jacobianFD.set( gid, dofIndex, dRdP );
          }
        }
      }
    }

    // b) compute all the derivatives wrt to the connection in WELL elem iwelem 
    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      globalIndex iwelemOffset = wellElemDofNumber[iwelem];

      {
        solver->ResetStateToBeginningOfStep(domain);
        
        // here is the perturbation in the pressure of the well element
        real64 const dRate = perturbParameter * (connRate[iwelem] + perturbParameter);
        dConnRate[iwelem] = dRate;

        residual.zero();
        assembleFunction( wellSolver, domain, &jacobian, &residual, &dofManager );

        globalIndex const dofIndex = integer_conversion<globalIndex>(iwelemOffset + SinglePhaseWell::ColOffset::DRATE );

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
    char const * filename = "testReservoirSinglePhaseMSWells.xml";

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


TEST_F(ReservoirSolverTest, jacobianNumericalCheck_Perforation)
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
                         [&] ( SinglePhaseWell * const targetSolver,
                               DomainPartition * const targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager ) -> void
  {
    targetSolver->AssemblePerforationTerms( time, dt, targetDomain, targetDofManager, targetJacobian, targetResidual );
  });
  
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
                         [&] ( SinglePhaseWell * const targetSolver,
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
                         [&] ( SinglePhaseWell * const targetSolver,
                               DomainPartition * const targetDomain,
                               ParallelMatrix * targetJacobian,
                               ParallelVector * targetResidual,
                               DofManager const * targetDofManager ) -> void
  {
    targetSolver->FormControlEquation( targetDomain, targetDofManager, targetJacobian, targetResidual );
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
                         [&] ( SinglePhaseWell * const targetSolver,
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
