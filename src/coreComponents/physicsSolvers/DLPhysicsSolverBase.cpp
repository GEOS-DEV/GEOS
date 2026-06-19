/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "DLPhysicsSolverBase.hpp"
#include "PhysicsSolverManager.hpp"

#include "common/MpiWrapper.hpp"
#include "codingUtilities/RTTypes.hpp"
#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "common/format/LogPart.hpp"
#include "common/TimingMacros.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "mesh/DomainPartition.hpp"
#include "math/interpolation/Interpolation.hpp"
#include "common/Timer.hpp"
#include "common/Units.hpp"

#if defined(GEOS_USE_PYGEOSX)
#include "python/PySolverType.hpp"
#endif

namespace geos
{

using namespace dataRepository;

DLPhysicsSolverBase::DLPhysicsSolverBase( string const & name,
                                          Group * const parent )
  : PhysicsSolverBase( name, parent ),
  m_sharedMemoryManager( "DLSharedMemoryManager", parent )
{}

DLPhysicsSolverBase::~DLPhysicsSolverBase() = default;

void DLPhysicsSolverBase::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution,
                                       bool const setSparsity )
{
  GEOS_MARK_FUNCTION;
  //TODO: Consider cancelling initalization of unused vectors/matrices in PhysicsSolverBase
  PhysicsSolverBase::setupSystem( domain,
                                  dofManager,
                                  localMatrix,
                                  rhs,
                                  solution,
                                  setSparsity );

  // initialize Data Vectors needed for DL simulations
  m_dofXCoords.setName( this->getName() + "/dofXCoords" );
  m_dofXCoords.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  m_dofYCoords.setName( this->getName() + "/dofYCoords" );
  m_dofYCoords.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  m_dofZCoords.setName( this->getName() + "/dofZCoords" );
  m_dofZCoords.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  m_prevSolution.setName( this->getName() + "/prevSolution" );
  m_prevSolution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  m_strainTrace.setName( this->getName() + "/strainTrace" );
  m_strainTrace.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

}

real64 DLPhysicsSolverBase::nonlinearImplicitStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   integer const cycleNumber,
                                                   DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                         GEOS_FMT( "--CustomLog--DLPhysicsSolverBase::nonlinearImplicitStep time:{} dt:{} cycle:{} solver:{} catalogName:{} ",
                                   time_n, dt, cycleNumber, getName(), getCatalogName()));

  // In a DL approach, there is no forseeable need for a dt cut during the course of this step
  real64 stepDt = dt;

  integer const maxNumberDtCuts = m_nonlinearSolverParameters.m_maxTimeStepCuts;
  real64 const dtCutFactor = m_nonlinearSolverParameters.m_timeStepCutFactor;

  bool const allowNonConverged = m_nonlinearSolverParameters.m_allowNonConverged > 0;

  integer & dtAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;

  integer const & maxConfigurationIter = m_nonlinearSolverParameters.m_maxNumConfigurationAttempts;

  integer & configurationLoopIter = m_nonlinearSolverParameters.m_numConfigurationAttempts;

  bool isConfigurationLoopConverged = false;

  GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                         GEOS_FMT(
                           "--CustomLog--DLPhysicsSolverBase::nonlinearImplicitStep  maxNumberDtCuts:{} dtCutFactor:{} allowNonConverged:{} dtAttempt:{} maxConfigurationIter:{} configurationLoopIter:{} ",
                           maxNumberDtCuts, dtCutFactor, allowNonConverged, dtAttempt, maxConfigurationIter, configurationLoopIter ));

  // TODO: check if DL methods need/enable a dt cutting strategy
  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      Timer timer( m_timers.get_inserted( "reset state" ) );

      resetStateToBeginningOfStep( domain );
      resetConfigurationToBeginningOfStep( domain );
    }

    // it's the simplest configuration that can be attempted whenever Newton's fails as a last resource.
    bool attemptedSimplestConfiguration = false;

    // TODO: check if configurationLoop is needed in DL methods, or if a single configuration is sufficient.
    // Configuration loop
    for( configurationLoopIter = 0; configurationLoopIter < maxConfigurationIter; ++configurationLoopIter )
    {
      if( isLogLevelActive< logInfo::NonlinearSolver >( getLogLevel() ) )
      {
        outputConfigurationStatistics( domain );
      }

      bool const isNewtonConverged = solveNonlinearSystemUsingDL( time_n,
                                                                  stepDt,
                                                                  cycleNumber,
                                                                  domain );

      if( isNewtonConverged )
      {
        isConfigurationLoopConverged = updateConfiguration( domain, configurationLoopIter );

        if( isConfigurationLoopConverged )
        {
          break;   // get out of configuration loop coz everything converged.
        }
        else
        {
          // increment the solver statistics for reporting purposes
          getIterationStats().incrementConfigIteration();
          GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                                 "---------- Configuration did not converge. Testing new configuration. ----------" );
        }
      }
      else if( !attemptedSimplestConfiguration )
      {
        resetStateToBeginningOfStep( domain );
        bool const breakLoop = resetConfigurationToDefault( domain );
        attemptedSimplestConfiguration = true;
        if( breakLoop )
        {
          break;
        }
        else
        {
          GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                                 "---------- Restarting Newton loop using default configuration. ----------" );
        }
      }
      else
      {
        break;   // get out of configuration loop and cut the time-step if you can.
      }

    }  // end of configuration loop

    if( isConfigurationLoopConverged )
    {
      // get out of outer loop
      break;
    }
    else
    {
      // TODO: update the logic for dt cutting based on the DL method specifics
      // cut timestep, go back to beginning of step and restart the Newton loop
      const real64 old_stepDt = stepDt;
      stepDt *= dtCutFactor;
      m_numTimestepsSinceLastDtCut = 0;
      GEOS_LOG_LEVEL_RANK_0( logInfo::TimeStep, GEOS_FMT( "\nWarning!\n ----------------------------------- \nTime step cut from {} to {}", old_stepDt, stepDt ));
      GEOS_LOG_LEVEL_RANK_0( logInfo::TimeStep, "Ensure that the DL model can handle the new time step appropriately.\n -----------------------------------" );

      // notify the solver statistics counter that this is a time step cut
      getIterationStats().updateTimeStepCut();
      getIterationStats().writeIterationStatsToTable();
    }
  }   // end of outer loop (dt chopping strategy)

  if( !isConfigurationLoopConverged )
  {
    GEOS_LOG_RANK_0( "Convergence not achieved." );

    if( allowNonConverged )
    {
      GEOS_LOG_RANK_0( "The accepted solution may be inaccurate." );
    }
    else
    {
      GEOS_ERROR( "Nonconverged solutions not allowed. Terminating..." );
    }
  }

  // return the achieved timestep
  return stepDt;
}

bool DLPhysicsSolverBase::solveNonlinearSystemUsingDL( real64 const & time_n,
                                                       real64 const & stepDt,
                                                       integer const cycleNumber,
                                                       DomainPartition & domain )
{
  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer & dtAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;
  integer & configurationLoopIter = m_nonlinearSolverParameters.m_numConfigurationAttempts;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

  // keep residual from previous iteration in case we need to do a line search
  real64 lastResidual = 1e99;
  integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
  real64 scaleFactor = 1.0;

  bool isNewtonConverged = false;

  GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                         GEOS_FMT( "--CustomLog--DLPhysicsSolverBase::solveNonlinearSystemUsingDL time:{} dt:{} cycle:{} solver:{} catalogName:{} ",
                                   time_n, stepDt, cycleNumber, getName(), getCatalogName()));
  GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                         GEOS_FMT( "--CustomLog--DLPhysicsSolverBase::solveNonlinearSystemUsingDL maxNewtonIter:{} dtAttempt:{} configurationLoopIter:{} minNewtonIter:{} newtonTol:{} ",
                                   maxNewtonIter, dtAttempt, configurationLoopIter, minNewtonIter, newtonTol ));

  for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
  {

    GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                           GEOS_FMT( "  Attempt: {:2}, ConfigurationIter: {:2}, NewtonIter: {:2}",
                                     dtAttempt, configurationLoopIter, newtonIter ));

    // assemble system matrix and rhs
    {
      Timer timer( m_timers.get_inserted( "assemble" ) );

      // We sync the nonlinear convergence history. The coupled solver parameters are the one being
      // used. We want to propagate the info to subsolvers. It can be important for solvers that
      // have special treatment for specific iterations.
      synchronizeNonlinearSolverParameters();

      // zero out matrix/rhs before assembly
      m_localMatrix.zero();
      m_rhs.zero();

      arrayView1d< real64 > const localRhs = m_rhs.open();

      // call assemble to fill the matrix and the rhs
      // NOTE: in a DL approach, assembleSystem populates the required inputs for the DL model
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

      if( m_assemblyCallback )
      {
        // Make a copy of LA objects and ship off to the callback
        array1d< real64 > localRhsCopy( m_rhs.localSize());
        localRhsCopy.setValues< parallelDevicePolicy<> >( m_rhs.values());
        m_assemblyCallback( m_localMatrix, std::move( localRhsCopy ));
      }
    }

    // TODO: Consider deleting (residual norm calculation and checks)
    // calculate residual norm
    real64 residualNorm = 0;
    {
      Timer timer( m_timers.get_inserted( "convergence check" ) );

      // get residual norm
      residualNorm = calculateResidualNorm( time_n, stepDt, domain, m_dofManager, m_rhs.values());
      GEOS_LOG_LEVEL_RANK_0( logInfo::ResidualNorm,
                             GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );
      getConvergenceStats().setResidualValue( "R", residualNorm );
      updateAndWriteConvergenceStep( time_n, stepDt, cycleNumber, newtonIter );
    }

    // TODO: Consider deleting (check for Newton convergence)
    // // if the residual norm is less than the Newton tolerance we denote that we have
    // // converged and break from the Newton loop immediately.
    // if (residualNorm < newtonTol && newtonIter >= minNewtonIter)
    // {
    //   isNewtonConverged = true;
    //   break;
    // }

    // TODO: Consider deleting (check for max allowed residual norm)
    // // if the residual norm is above the max allowed residual norm, we break from
    // // the Newton loop to avoid crashes due to Newton divergence
    // if (residualNorm > m_nonlinearSolverParameters.m_maxAllowedResidualNorm)
    // {
    //   string const maxAllowedResidualNormString = NonlinearSolverParameters::viewKeysStruct::maxAllowedResidualNormString();
    //   GEOS_LOG_LEVEL_RANK_0(logInfo::Convergence,
    //                         GEOS_FMT("    The residual norm is above the {} of {}. Newton loop terminated.",
    //                                  maxAllowedResidualNormString,
    //                                  m_nonlinearSolverParameters.m_maxAllowedResidualNorm));
    //   isNewtonConverged = false;
    //   break;
    // }

    // TODO: Consider deleting (line search strategy)
    // // do line search in case residual has increased
    // if (m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None && residualNorm >
    // lastResidual * m_nonlinearSolverParameters.m_lineSearchResidualFactor && newtonIter >=
    // m_nonlinearSolverParameters.m_lineSearchStartingIteration)
    // {
    //   bool lineSearchSuccess = false;
    //   if (m_nonlinearSolverParameters.m_lineSearchInterpType == NonlinearSolverParameters::LineSearchInterpolationType::Linear)
    //   {
    //     residualNorm = lastResidual;
    //     lineSearchSuccess = lineSearch(time_n,
    //                                    stepDt,
    //                                    cycleNumber,
    //                                    domain,
    //                                    m_dofManager,
    //                                    m_localMatrix.toViewConstSizes(),
    //                                    m_rhs,
    //                                    m_solution,
    //                                    scaleFactor,
    //                                    residualNorm);
    //   }
    //   else
    //   {
    //     lineSearchSuccess = lineSearchWithParabolicInterpolation(time_n,
    //                                                              stepDt,
    //                                                              cycleNumber,
    //                                                              domain,
    //                                                              m_dofManager,
    //                                                              m_localMatrix.toViewConstSizes(),
    //                                                              m_rhs,
    //                                                              m_solution,
    //                                                              scaleFactor,
    //                                                              lastResidual,
    //                                                              residualNorm);
    //   }

    //   if (!lineSearchSuccess)
    //   {
    //     if (m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Attempt)
    //     {
    //       GEOS_LOG_LEVEL_RANK_0(logInfo::LineSearch,
    //                             "        Line search failed to produce reduced residual. Accepting iteration.");
    //     }
    //     else if (m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Require)
    //     {
    //       // if line search failed, then break out of the main Newton loop. Timestep will be cut.
    //       GEOS_LOG_LEVEL_RANK_0(logInfo::LineSearch,
    //                             "        Line search failed to produce reduced residual. Exiting Newton Loop.");
    //       break;
    //     }
    //   }
    // }

    // TODO: Applying Linear Solver (to be replaced with DL model calls)
    {
      Timer timer( m_timers.get_inserted( "linear solver total" ) );

      // if using adaptive Krylov tolerance scheme, update tolerance.
      LinearSolverParameters::Krylov & krylovParams = m_linearSolverParameters.get().krylov;
      if( krylovParams.useAdaptiveTol )
      {
        krylovParams.relTolerance = newtonIter > 0 ? eisenstatWalker( residualNorm, lastResidual, krylovParams ) : krylovParams.weakestTol;
      }

      {
        Timer timer_setup( m_timers.get_inserted( "linear solver create" ) );

        // Compose parallel LA matrix/rhs out of local LA matrix/rhs
        //
        m_matrix.create( m_localMatrix.toViewConst(), m_dofManager.numLocalDofs(), MPI_COMM_GEOS );
      }

      // Output the linear system matrix/rhs for debugging purposes
      debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs );

      // Solve the linear system
      solveLinearSystem( m_dofManager, m_matrix, m_rhs, m_solution );

      {
        // TODO: add timer for DL implementation
        // Timer timer_destroy(m_timers["DL imlementaiton"]);
        // Print m_solution size (how many DOFs are being solved)
        GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver, GEOS_FMT( "--CustomLog--DLPhysicsSolverBase::solveNonlinearSystemUsingDL - Solution size: {}", m_solution.globalSize()));
        shareDLModelInputs( time_n, stepDt, cycleNumber, domain );
        readDLModelOutputs( time_n, stepDt, cycleNumber, domain );
      }

      // Increment the solver statistics for reporting purposes
      getIterationStats().updateNonlinearIteration( m_linearSolverResult.numIterations );

      // Output the linear system solution for debugging purposes
      debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );

      // TODO: Consider deleting
      // // Do not allow non converged linear solver - cut time step
      // if (!m_allowNonConvergedLinearSolverSolution && m_linearSolverResult.status == LinearSolverResult::Status::NotConverged)
      // {
      //   return false;
      // }
    }

    // apply solution and update state
    // TODO: check the scale factor logic in DL methods & checkSystemSolution implementation
    {
      Timer timer( m_timers.get_inserted( "apply solution" ) );

      // Compute the scaling factor for the Newton update
      scaleFactor = scalingForSystemSolution( domain, m_dofManager, m_solution.values());

      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Global solution scaling factor = {}", getName(), scaleFactor ));

      // TODO : Consider deleting (check system solution validity)
      // if (!checkSystemSolution(domain, m_dofManager, m_solution.values(), scaleFactor))
      // {
      //   // TODO try chopping (similar to line search)
      //   GEOS_LOG_RANK_0(GEOS_FMT("    {}: Solution check failed. Newton loop terminated.", getName()));
      //   break;
      // }

      // apply the system solution to the fields/variables
      applySystemSolution( m_dofManager, m_solution.values(), scaleFactor, stepDt, domain );
    }

    {
      Timer timer( m_timers.get_inserted( "update state" ) );

      // update non-primary variables (constitutive models)
      updateState( domain );
    }

    lastResidual = residualNorm;

    // Terminate Newton loop if converged
    // TODO : consider future scenarios. For now, DL is assumed to always converge within one iteration.
    {
      isNewtonConverged = true;
      break;
    }
  }

  return isNewtonConverged;
}

void DLPhysicsSolverBase::initializePostInitialConditionsPostSubGroups()
{

  initializeSharedMemories();
}


void DLPhysicsSolverBase::populateDofCoords( DofManager const & /*dofManager*/,
                                             DomainPartition & domain,
                                             string const & elemDofKey )
{
  GEOS_MARK_FUNCTION;

  // TODO : Check implementation of
  // SolidMechanicsLagrangianFEM::registerDataOnMesh
  // PoromechanicsSolver::registerDataOnMesh
  // DLFlowSolverBase::registerDataOnMesh
  // DLPhysicsSolverBase::registerDataOnMesh

  // open device-accessible views for writing
  arrayView1d< real64 > dofX = m_dofXCoords.open();
  arrayView1d< real64 > dofY = m_dofYCoords.open();
  arrayView1d< real64 > dofZ = m_dofZCoords.open();

  // initialize with large values so missing entries are noticeable
  forAll< parallelHostPolicy >( m_dofXCoords.localSize(), [=] ( localIndex const i )
  {
    dofX[i] = 1e99;
    dofY[i] = 1e99;
    dofZ[i] = 1e99;
  } );

  globalIndex const lowerbound = m_solution.ilower();
  globalIndex const upperbound = m_solution.iupper();

  // iterate over mesh bodies & element subregions, collect element-centers and element->dof mapping
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      // element centers
      arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();

      // element to dof mapping
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );

      // for each element, store its center in the corresponding DOF slot
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        globalIndex const eGDof = elemDofNumber[ei];
        // only fill entries owned by this rank
        if( eGDof >= lowerbound && eGDof < upperbound )
        {
          dofX[eGDof - lowerbound] = elemCenter[ei][0];
          dofY[eGDof - lowerbound] = elemCenter[ei][1];
          dofZ[eGDof - lowerbound] = elemCenter[ei][2];
        }
      } );
    } );
  } );

  m_dofXCoords.close();
  m_dofYCoords.close();
  m_dofZCoords.close();


  // // redirect to shared mem stream output
  // static std::ofstream logFile_shared_mem("geosx_sharedMem_rank_" + std::to_string(MpiWrapper::commRank(MPI_COMM_GEOS)) + ".log",
  // std::ios::out | std::ios::app );
  // std::cout.rdbuf(logFile_shared_mem.rdbuf());

  // Temporary debug changes in m_dofXCoords, m_dofYCoords, m_dofZCoords
  GEOS_LOG_LEVEL( logInfo::NonlinearSolver,
                  GEOS_FMT( "--CustomLog--DLSinglePhaseBase::populateDofCoords: Dof X Coords Global Size:{}, Local Size:{} ", m_dofXCoords.globalSize(), m_dofXCoords.localSize()));
  m_dofXCoords.print();

  GEOS_LOG_LEVEL( logInfo::NonlinearSolver,
                  GEOS_FMT( "--CustomLog--DLSinglePhaseBase::populateDofCoords: Dof Y Coords Global Size:{}, Local Size:{} ", m_dofYCoords.globalSize(), m_dofYCoords.localSize()));
  m_dofYCoords.print();

  GEOS_LOG_LEVEL( logInfo::NonlinearSolver,
                  GEOS_FMT( "--CustomLog--DLSinglePhaseBase::populateDofCoords: Dof Z Coords Global Size:{}, Local Size:{} ", m_dofZCoords.globalSize(), m_dofZCoords.localSize()));
  m_dofZCoords.print();


  // // redirect to default stream output
  // static std::ofstream logFile("geosx_output_rank_" + std::to_string(MpiWrapper::commRank(MPI_COMM_GEOS)) + ".log",
  // std::ios::out | std::ios::app );
  // std::cout.rdbuf(logFile.rdbuf());

}



void DLPhysicsSolverBase::populateDofStrainTrace( DofManager const & /*dofManager*/,
                                                  DomainPartition & domain,
                                                  string const & elemDofKey )
{
  GEOS_MARK_FUNCTION;
  // open device-accessible views for writing
  arrayView1d< real64 > dofStrainTrace = m_strainTrace.open();
  // initialize with large values so missing entries are noticeable
  forAll< parallelHostPolicy >( m_strainTrace.localSize(), [=] ( localIndex const i )
  {
    dofStrainTrace[i] = 1e99;
  } );

  globalIndex const lowerbound = m_solution.ilower();
  globalIndex const upperbound = m_solution.iupper();

  // iterate over mesh bodies & element subregions
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );

      fields::solidMechanics::arrayView2dLayoutStrain const strain =
        subRegion.getField< fields::solidMechanics::averageStrain >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex ei )
      {
        // TODO: Check strain components indices
        real64 const trace = strain[ei][0] + strain[ei][1] + strain[ei][2];
        globalIndex const eGDof = elemDofNumber[ei];
        if( eGDof >= lowerbound && eGDof < upperbound )
        {
          dofStrainTrace[eGDof - lowerbound] = trace;
        }
      } );
    } );
  } );

  m_strainTrace.close();


  // // redirect to shared mem stream output
  // static std::ofstream logFile_shared_mem("geosx_sharedMem_rank_" + std::to_string(MpiWrapper::commRank(MPI_COMM_GEOS)) + ".log",
  // std::ios::out | std::ios::app );
  // std::cout.rdbuf(logFile_shared_mem.rdbuf());

  // Temporary debug changes in m_strainTrace
  GEOS_LOG_LEVEL( logInfo::NonlinearSolver,
                  GEOS_FMT( "--CustomLog--DLSinglePhaseBase::populateDofStrainTrace  Dof Strain Trace Global Size:{}, Local Size:{} ", m_strainTrace.globalSize(), m_strainTrace.localSize()));
  m_strainTrace.print();

  // // redirect to default stream output
  // static std::ofstream logFile("geosx_output_rank_" + std::to_string(MpiWrapper::commRank(MPI_COMM_GEOS)) + ".log",
  // std::ios::out | std::ios::app );
  // std::cout.rdbuf(logFile.rdbuf());
}

//TODO: Move to DLPhysicsSolverBase and make it abstract in DLPhysicsSolverBase
void DLPhysicsSolverBase::populateDofPrevSolution( DofManager const & /*dofManager*/,
                                                   DomainPartition & domain,
                                                   string const & elemDofKey )
{
  GEOS_MARK_FUNCTION;

  // Open destination view
  arrayView1d< real64 > prevSolView = m_prevSolution.open();

  // Initialize with a large value for easy spotting of gaps
  forAll< parallelHostPolicy >( m_prevSolution.localSize(), [=] ( localIndex const i )
  {
    prevSolView[i] = 1e99;
  } );

  // Local ownership bounds for DOF indices
  globalIndex const lowerbound = m_solution.ilower();
  globalIndex const upperbound = m_solution.iupper();


  // Iterate mesh targets and subregions
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );

      // Previous pressure time level (pressure_n)
      // TODO: Double check flow::pressure_n vs flow::pressure
      arrayView1d< real64 const > const pres_n =
        subRegion.getField< fields::flow::pressure_n >();

      // Map element scalar to owned DOF slot
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex ei )
      {
        globalIndex const eGDof = elemDofNumber[ei];
        if( eGDof >= lowerbound && eGDof < upperbound )
        {
          prevSolView[eGDof - lowerbound] = pres_n[ei];
        }
      } );
    } );
  } );

  m_prevSolution.close();

  // // redirect to shared mem stream output
  // static std::ofstream logFile_shared_mem("geosx_sharedMem_rank_" + std::to_string(MpiWrapper::commRank(MPI_COMM_GEOS)) + ".log",
  // std::ios::out | std::ios::app );
  // std::cout.rdbuf(logFile_shared_mem.rdbuf());

  // Temporary debug changes in m_prevSolution
  GEOS_LOG_LEVEL( logInfo::NonlinearSolver,
                  GEOS_FMT( "--CustomLog--DLSinglePhaseBase::populateDofPrevSolution  Dof Prev Solution Global Size:{}, Local Size:{} ", m_prevSolution.globalSize(), m_prevSolution.localSize()));
  m_prevSolution.print();

  // // redirect to default stream output
  // static std::ofstream logFile("geosx_output_rank_" + std::to_string(MpiWrapper::commRank(MPI_COMM_GEOS)) + ".log",
  // std::ios::out | std::ios::app );
  // std::cout.rdbuf(logFile.rdbuf());
}

void DLPhysicsSolverBase::initializeSharedMemories()
{
  GEOS_ERROR( "DLPhysicsSolverBase::initializeSharedMemories called!. Should be overridden." );
}

void DLPhysicsSolverBase::shareDLModelInputs( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time_n, dt, cycleNumber, domain );
  GEOS_ERROR( "DLPhysicsSolverBase::shareDLModelInputs called!. Should be overridden." );
}

void DLPhysicsSolverBase::readDLModelOutputs( real64 const & time_n,
                                              real64 const & dt,
                                              integer const cycleNumber,
                                              DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time_n, dt, cycleNumber, domain );
  GEOS_ERROR( "DLPhysicsSolverBase::readDLModelOutputs called!. Should be overridden." );
}

#if defined(GEOS_USE_PYGEOSX)
PyTypeObject *DLPhysicsSolverBase::getPythonType() const
{
  return python::getPySolverType();
}
#endif

} // namespace geos
