/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsLagrangianContactSolver.cpp
 *
 */


#include "SinglePhasePoromechanicsLagrangianContactSolver.hpp"

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
//#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactSolver.hpp"
#include "physicsSolvers/multiphysics/LagrangianContactFlowSolver.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "SinglePhasePoromechanicsKernel.hpp"
#include "math/interpolation/Interpolation.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace interpolation;

SinglePhasePoromechanicsLagrangianContactSolver::SinglePhasePoromechanicsLagrangianContactSolver( const string & name,
                                                                Group * const parent ):
  SolverBase( name, parent ),
  m_contactSolverName(),
  m_flowSolverName()

{
  registerWrapper( viewKeyStruct::contactSolverNameString(), &m_contactSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the contact mechanics solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::contactFlowSolverNameString(), &m_contactFlowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the contact mechanics and flow solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::porousMaterialNamesString(), &m_porousMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the material that should be used in the constitutive updates" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhasePoromechanicsLagrangianContactSolver::setupDofs( DomainPartition const & domain,
                                                DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  //m_contactFlowSolver->setupDofs( domain, dofManager );
  //m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addField( keys::TotalDisplacement,
                       DofManager::Location::Node,
                       3,
                       targetRegionNames() );

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connector::Elem );

  // restrict coupling to fracture regions only
  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  dofManager.addField( LagrangianContactSolver::viewKeyStruct::tractionString(),
                       DofManager::Location::Elem,
                       3,
                       fractureRegions );
  dofManager.addCoupling( LagrangianContactSolver::viewKeyStruct::tractionString(),
                          LagrangianContactSolver::viewKeyStruct::tractionString(),
                          DofManager::Connector::Face,
                          fractureRegions );
  dofManager.addCoupling( keys::TotalDisplacement,
                          LagrangianContactSolver::viewKeyStruct::tractionString(),
                          DofManager::Connector::Elem,
                          fractureRegions );

  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          extrinsicMeshData::flow::pressure::key(),
                          DofManager::Connector::Elem );
  dofManager.addCoupling( extrinsicMeshData::flow::pressure::key(),
                          LagrangianContactSolver::viewKeyStruct::tractionString(),
                          DofManager::Connector::None );
//  dofManager.addCoupling( LagrangianContactSolver::viewKeyStruct::tractionString(),
//                          FlowSolverBase::viewKeyStruct::pressureString(),
//                          DofManager::Connector::None );

//  dofManager.addCoupling( FlowSolverBase::viewKeyStruct::pressureString(),
//                          LagrangianContactSolver::viewKeyStruct::tractionString(),
//                          DofManager::Connector::None );
//  dofManager.addCoupling( LagrangianContactSolver::viewKeyStruct::tractionString(),
//                          FlowSolverBase::viewKeyStruct::pressureString(),
//                          DofManager::Connector::None );

  dofManager.printFieldInfo();
}

void SinglePhasePoromechanicsLagrangianContactSolver::setupSystem( DomainPartition & domain,
                                                  DofManager & dofManager,
                                                  CRSMatrix< real64, globalIndex > & localMatrix,
                                                  ParallelVector & rhs,
                                                  ParallelVector & solution,
                                                  bool const setSparsity )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  m_contactFlowSolver->setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  //ParallelMatrix MM;
  //MM.create( localMatrix.toViewConst(), MPI_COMM_GEOSX );
  //MM.print();

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner();
  }
}

void SinglePhasePoromechanicsLagrangianContactSolver::implicitStepSetup( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  m_contactSolver->implicitStepSetup( time_n, dt, domain );
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  //m_contactFlowSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition & domain )
{
  m_contactSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, porousMaterialNames()[targetIndex] );
    porousMaterial.saveConvergedState();
  } );

  // Laura print max pressure and displacement
  //MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  // pressure
  ElementRegionManager & elemManager = mesh.getElemManager();
  elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion & region )
  {
    region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      if( subRegion.hasWrapper( extrinsicMeshData::flow::pressure::key() ) )
      {
        arrayView1d< real64 > pres = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
        double * max_pres = std::max_element(pres.begin(), pres.end());
        GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- max pres {:15.6e}", * max_pres ) );
        double * min_pres = std::min_element(pres.begin(), pres.end());
        GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- min pres {:15.6e}", * min_pres ) );
      }
    } );
  } );
  // displacement
  NodeManager & nodeManager = mesh.getNodeManager();
  arrayView2d< real64 , nodes::TOTAL_DISPLACEMENT_USD > const disp = nodeManager.totalDisplacement();
  double * min_disp = std::min_element(disp.begin(), disp.end());
  GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- min disp {:15.6e}", * min_disp ) );

  real64 const totalFlux = m_flowSolver->computeFluxFaceDirichlet( time_n, dt, domain );
  GEOSX_LOG_RANK_0( GEOSX_FMT( "SinglePhasePoromechanicsLagrangianContactSolver::implicitStepComplete -- total flux through Dirichlet faces {:15.6e}", totalFlux ) );
  // end Laura
}

void SinglePhasePoromechanicsLagrangianContactSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< SinglePhaseBase >( m_flowSolverName );
  //m_contactFlowSolver = &this->getParent().getGroup< LagrangianContactSolver >( m_contactSolverName );
  m_contactSolver = &this->getParent().getGroup< LagrangianContactSolver >( m_contactSolverName );
  m_contactFlowSolver = &this->getParent().getGroup< LagrangianContactFlowSolver >( m_contactFlowSolverName );
}

void SinglePhasePoromechanicsLagrangianContactSolver::initializePostInitialConditionsPreSubGroups()
{
  if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
  }
}

SinglePhasePoromechanicsLagrangianContactSolver::~SinglePhasePoromechanicsLagrangianContactSolver()
{
  // TODO Auto-generated destructor stub
}

void SinglePhasePoromechanicsLagrangianContactSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_contactFlowSolver->resetStateToBeginningOfStep( domain );
}

real64 SinglePhasePoromechanicsLagrangianContactSolver::solverStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   int const cycleNumber,
                                                   DomainPartition & domain )
{
  real64 dt_return = dt;

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_rhs,
               m_solution );

  implicitStepSetup( time_n, dt, domain );

  dt_return = this->nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

real64 SinglePhasePoromechanicsLagrangianContactSolver::nonlinearImplicitStep( real64 const & time_n,
                                                           real64 const & dt,
                                                           integer const cycleNumber,
                                                           DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  // dt may be cut during the course of this step, so we are keeping a local
  // value to track the achieved dt for this step.
  real64 stepDt = dt;

  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

  integer const maxNumberDtCuts = m_nonlinearSolverParameters.m_maxTimeStepCuts;
  real64 const dtCutFactor = m_nonlinearSolverParameters.m_timeStepCutFactor;

  bool const allowNonConverged = m_nonlinearSolverParameters.m_allowNonConverged > 0;

  integer & dtAttempt = m_nonlinearSolverParameters.m_numdtAttempts;

  // a flag to denote whether we have converged
  bool isNewtonConverged = false;

  bool isActiveSetConverged = false;

  //bool useElasticStep = !m_contactSolver->isFractureAllInStickCondition( domain );

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // TMP
    bool useElasticStep = !m_contactSolver->isFractureAllInStickCondition( domain );

    std::cout << "useElasticStep " << useElasticStep << std::endl;

    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      resetStateToBeginningOfStep( domain );
      globalIndex numStick, numSlip, numOpen;
      m_contactSolver->computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
    }

    integer & activeSetIter = m_activeSetIter;
    for( activeSetIter = 0; activeSetIter < m_contactSolver->getActiveSetMaxIter(); ++activeSetIter )
    {
      // *******************************
      // Newton loop: begin
      // *******************************
      isNewtonConverged = false;
      // keep residual from previous iteration in case we need to do a line search
      real64 lastResidual = 1e99;
      integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
      real64 scaleFactor = 1.0;

      // main Newton loop
      bool computeResidual = true;
      for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
      {
        if( getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          char output[55] = {0};
          sprintf( output, "    Attempt: %2d, ActiveSetIter: %2d ; NewtonIter: %2d ; ",
                   dtAttempt, activeSetIter, newtonIter );
          std::cout<<output<<std::endl;
        }

        // zero out matrix/rhs before assembly
        m_localMatrix.zero();
        m_rhs.zero();

        {
          arrayView1d< real64 > const localRhs = m_rhs.open();

          // call assemble to fill the matrix and the rhs
          assembleSystem( time_n,
                          stepDt,
                          domain,
                          m_dofManager,
                          m_localMatrix.toViewConstSizes(),
                          localRhs );

          // apply boundary conditions to system
          applyBoundaryConditions( time_n,
                                   stepDt,
                                   domain,
                                   m_dofManager,
                                   m_localMatrix.toViewConstSizes(),
                                   localRhs );

          m_rhs.close();
        }
        // TODO: maybe add scale function here?
        // Scale()

        real64 residualNorm;
        // get residual norm
        if( computeResidual )
        {
          residualNorm = calculateResidualNorm( domain, m_dofManager, m_rhs.values() );
        }
        else
        {
          residualNorm = lastResidual;
        }

        if( getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          if( newtonIter!=0 )
          {
            char output[46] = {0};
            sprintf( output,
                     "Last LinSolve(iter,tol) = (%4d, %4.2e) ; ",
                     m_linearSolverResult.numIterations,
                     m_linearSolverResult.residualReduction );
            std::cout<<output;
          }
          std::cout<<std::endl;
        }

        // if the residual norm is less than the Newton tolerance we denote that we have
        // converged and break from the Newton loop immediately. 
        if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
        {
          isNewtonConverged = true;
          break;
        }

        // if using adaptive Krylov tolerance scheme, update tolerance.
        LinearSolverParameters::Krylov & krylovParams = m_linearSolverParameters.get().krylov;
        if( krylovParams.useAdaptiveTol )
        {
          krylovParams.relTolerance = eisenstatWalker( residualNorm, lastResidual, krylovParams.weakestTol );
        }

        // Compose parallel LA matrix/rhs out of local LA matrix/rhs
        m_matrix.create( m_localMatrix.toViewConst(), m_dofManager.numLocalDofs(), MPI_COMM_GEOSX );

        // Output the linear system matrix/rhs for debugging purposes
        debugOutputSystem( time_n, cycleNumber, 10*activeSetIter+newtonIter, m_matrix, m_rhs );

        // Solve the linear system
        solveSystem( m_dofManager, m_matrix, m_rhs, m_solution );

        // Output the linear system solution for debugging purposes
        debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );

        scaleFactor = scalingForSystemSolution( domain, m_dofManager, m_solution.values() );

        // do line search in case residual has increased
        if( m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None && newtonIter > 0 )
        {
          bool lineSearchSuccess = lineSearch( time_n,
                                               stepDt,
                                               cycleNumber,
                                               domain,
                                               m_dofManager,
                                               m_localMatrix.toViewConstSizes(),
                                               m_rhs,
                                               m_solution,
                                               scaleFactor,
                                               residualNorm );

          if( !lineSearchSuccess )
          {
            if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Attempt )
            {
              GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search failed to produce reduced residual. Accepting iteration." );
            }
            else if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Require )
            {
              // if line search failed, then break out of the main Newton loop. Timestep will be cut.
              GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search failed to produce reduced residual. Exiting Newton Loop." );
              break;
            }
          }
          // Residual norm already computed in line search and stored in "residualNorm"
          computeResidual = false;
        }
        else
        {
          // apply the system solution to the fields/variables
          applySystemSolution( m_dofManager, m_solution.values(), scaleFactor, domain );
          // Need to compute the residual norm
          computeResidual = true;
        }

        if( !checkSystemSolution( domain, m_dofManager, m_solution.values(), scaleFactor ) )
        {
          // TODO try chopping (similar to line search)
          GEOSX_LOG_RANK_0( "    Solution check failed. Newton loop terminated." );
          break;
        }

        lastResidual = residualNorm;
      }
      // *******************************
      // Newton loop: end
      // *******************************

      // *******************************
      // Active set check: begin
      // *******************************
      bool const isPreviousFractureStateValid = m_contactSolver->updateFractureState( domain );
      GEOSX_LOG_LEVEL_RANK_0( 1, "active set flag: " << std::boolalpha << isPreviousFractureStateValid );

      if( getLogLevel() > 2 )
      {
        globalIndex numStick, numSlip, numOpen;
        m_contactSolver->computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
      }
      // *******************************
      // Active set check: end
      // *******************************

      GEOSX_LOG_LEVEL_RANK_0( 2, "isPreviousFractureStateValid: " << std::boolalpha << isPreviousFractureStateValid <<
                              " | isNewtonConverged: " << isNewtonConverged << " | useElasticStep: " << useElasticStep );
      if( isNewtonConverged )
      {
        isActiveSetConverged = isPreviousFractureStateValid;
        if( isActiveSetConverged )
        {
          break;
        }
      }
      else if( useElasticStep )
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, "Trying with an elastic step" );
        useElasticStep = false;
        resetStateToBeginningOfStep( domain );
        m_contactSolver->setFractureStateForElasticStep( domain );
        m_contactFlowSolver->updateOpeningForFlow( domain );
      }
      else
      {
        GEOSX_LOG_LEVEL_RANK_0( 2, "Newton did not converge in active set loop" );
        break;
      }
    }
    if( !isNewtonConverged )
    {
      // cut timestep, go back to beginning of step and restart the Newton loop
//      if (stepDt >= 0.25)
      stepDt *= dtCutFactor;
      GEOSX_LOG_LEVEL_RANK_0 ( 1, "New dt = " <<  stepDt );
    }
    if( isActiveSetConverged )
    {
      break;
    }
  }

  if( !isNewtonConverged )
  {
    GEOSX_LOG_RANK_0( "Convergence not achieved." );

    if( allowNonConverged )
    {
      GEOSX_LOG_RANK_0( "The accepted solution may be inaccurate." );
    }
    else
    {
      GEOSX_ERROR( "Nonconverged solutions not allowed. Terminating..." );
    }
  }

  if( !isActiveSetConverged )
  {
    //GEOSX_ERROR( "Active set did not reached a solution. Terminating..." );
  }
  else
  {
    GEOSX_LOG_RANK_0( "Number of active set iterations: " << m_activeSetIter );
  }

  // return the achieved timestep
  return stepDt;
}

bool SinglePhasePoromechanicsLagrangianContactSolver::lineSearch( real64 const & time_n,
                                              real64 const & dt,
                                              integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                              DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              ParallelVector & rhs,
                                              ParallelVector & solution,
                                              real64 const scaleFactor,
                                              real64 & lastResidual )
{
  bool lineSearchSuccess = true;

  integer const maxNumberLineSearchCuts = m_nonlinearSolverParameters.m_lineSearchMaxCuts;

  real64 const sigma1 = 0.5;
  real64 const alpha = 1.e-4;

  real64 localScaleFactor = scaleFactor;
  real64 lamm = scaleFactor;
  real64 lamc = localScaleFactor;
  integer lineSearchIteration = 0;

  // get residual norm
  real64 residualNorm0 = lastResidual;

  applySystemSolution( dofManager, solution.values(), scaleFactor, domain );

  // re-assemble system
  localMatrix.zero();
  rhs.zero();

  {
    arrayView1d< real64 > const localRhs = rhs.open();
    assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
    rhs.close();
  }

  // get residual norm
  real64 residualNormT = calculateResidualNorm( domain, dofManager, rhs.values() );

  real64 ff0 = residualNorm0*residualNorm0;
  real64 ffT = residualNormT*residualNormT;
  real64 ffm = ffT;
  real64 cumulativeScale = scaleFactor;

  while( residualNormT >= (1.0 - alpha*localScaleFactor)*residualNorm0 )
  {
    real64 const previousLocalScaleFactor = localScaleFactor;
    // Apply the three point parabolic model
    if( lineSearchIteration == 0 )
    {
      localScaleFactor *= sigma1;
    }
    else
    {
      localScaleFactor = parabolicInterpolationThreePoints( lamc, lamm, ff0, ffT, ffm );
    }

    // Update x; keep the books on lambda
    real64 const deltaLocalScaleFactor = ( localScaleFactor - previousLocalScaleFactor );
    cumulativeScale += deltaLocalScaleFactor;

    if( !checkSystemSolution( domain, dofManager, solution.values(), deltaLocalScaleFactor ) )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search " << lineSearchIteration << ", solution check failed" );
      continue;
    }

    applySystemSolution( dofManager, solution.values(), deltaLocalScaleFactor, domain );
    lamm = lamc;
    lamc = localScaleFactor;

    // Keep the books on the function norms
    // re-assemble system
    // TODO: add a flag to avoid a completely useless Jacobian computation: rhs is enough
    localMatrix.zero();
    rhs.zero();

    {
      arrayView1d< real64 > const localRhs = rhs.open();
      assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );
      applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
      rhs.close();
    }

    if( getLogLevel() >= 1 && logger::internal::rank==0 )
    {
      char output[100];
      sprintf( output, "        Line search @ %0.3f:      ", cumulativeScale );
      std::cout<<output;
    }

    // get residual norm
    residualNormT = calculateResidualNorm( domain, dofManager, rhs.values() );
    ffm = ffT;
    ffT = residualNormT*residualNormT;
    lineSearchIteration += 1;

    if( lineSearchIteration > maxNumberLineSearchCuts )
    {
      lineSearchSuccess = false;
      break;
    }
  }

  lastResidual = residualNormT;

  return lineSearchSuccess;
}

void SinglePhasePoromechanicsLagrangianContactSolver::assembleSystem( real64 const time_n,
                                                     real64 const dt,
                                                     DomainPartition & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{

  GEOSX_MARK_FUNCTION;

  // Need to synchronize the two iteration counters
  m_contactSolver->getNonlinearSolverParameters().m_numNewtonIterations = m_nonlinearSolverParameters.m_numNewtonIterations;

  MeshLevel & mesh = domain.getMeshBodies().getGroup< MeshBody >( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();

  string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

  string const pDofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );

  m_contactSolver->assembleForceResidualDerivativeWrtTraction( domain, dofManager, localMatrix, localRhs );
  m_contactSolver->assembleTractionResidualDerivativeWrtDisplacementAndTraction( domain, dofManager, localMatrix, localRhs );
  m_contactSolver->assembleStabilization( domain, dofManager, localMatrix, localRhs );

  m_contactFlowSolver->updateOpeningForFlow( domain );

  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  m_flowSolver->assembleHydrofracFluxTerms( time_n,
                                            dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs,
                                            m_contactFlowSolver->getDerivativeFluxResidual_dAperture(),
                                            fractureRegions );

  m_contactFlowSolver->assembleForceResidualDerivativeWrtPressure( domain, dofManager, localMatrix, localRhs );
  m_contactFlowSolver->assembleFluidMassResidualDerivativeWrtDisplacement( domain, dofManager, localMatrix, localRhs );
  m_contactFlowSolver->assembleStabilization( domain, dofManager, localMatrix, localRhs );

  m_contactFlowSolver->getRefDerivativeFluxResidual_dAperture()->zero();

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  PoromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                pDofKey,
                                                                dofManager.rankOffset(),
                                                                localMatrix,
                                                                localRhs,
                                                                gravityVectorData,
                                                                m_flowSolver->fluidModelNames() );

  // Cell-based contributions
  finiteElement::
    regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                  constitutive::PorousSolidBase,
                                  CellElementSubRegion >( mesh,
                                                          targetRegionNames(),
                                                          this->getDiscretizationName(),
                                                          porousMaterialNames(),
                                                          kernelFactory );

  m_flowSolver->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              dofManager.getKey( extrinsicMeshData::flow::pressure::key() ) );
}

void SinglePhasePoromechanicsLagrangianContactSolver::applyBoundaryConditions( real64 const time_n,
                                                              real64 const dt,
                                                              DomainPartition & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  m_contactSolver->applyBoundaryConditions( time_n, dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs );

  //ParallelMatrix MM;
  //MM.create( localMatrix.toViewConst(), MPI_COMM_GEOSX );
  //std::cout << MM << std::endl;

  m_flowSolver->applyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64 SinglePhasePoromechanicsLagrangianContactSolver::calculateResidualNorm( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_contactSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output, "    ( Rsolid, Rfluid, Rtotal ) = ( %4.2e, %4.2e, %4.2e )", momementumResidualNorm, massResidualNorm,
             sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm ) );
    std::cout << output << std::endl;
  }

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void SinglePhasePoromechanicsLagrangianContactSolver::createPreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( m_contactFlowSolver->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { keys::TotalDisplacement, { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( m_flowSolver->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { extrinsicMeshData::flow::pressure::key(), { 1, true } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void SinglePhasePoromechanicsLagrangianContactSolver::solveSystem( DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs,
                                                  ParallelVector & solution )
{
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );

  int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  if( rank == 0 )
  {
    string str;
    std::getline( std::cin, str );
    if( str.length() > 0 )
    {
      debugOutputSolution( 99, 99, 99, solution );
      GEOSX_ERROR( "STOP" );
    }
  }
  MpiWrapper::barrier( MPI_COMM_GEOSX );
}

void SinglePhasePoromechanicsLagrangianContactSolver::applySystemSolution( DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localSolution,
                                                          real64 const scalingFactor,
                                                          DomainPartition & domain )
{
  /*
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();
  elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion & region )
  {
    region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      if( subRegion.hasWrapper( SinglePhaseBase::viewKeyStruct::pressureString() ) )
      {
        arrayView1d< real64 > pres = subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString() );
        arrayView1d< real64 > dpres = subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString() );
        std::cout << "pres0 " << pres[0] << " " << dpres[0] << std::endl;
      }
    } );
  } );
  */

  // update displacement field
  m_contactSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

void SinglePhasePoromechanicsLagrangianContactSolver::updateState( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  this->template forTargetSubRegions< CellElementSubRegion >( mesh, [&] ( localIndex const targetIndex,
                                                                          auto & subRegion )
  {
    m_flowSolver->updateFluidState( subRegion, targetIndex );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsLagrangianContactSolver, string const &, Group * const )

} /* namespace geosx */
