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

#include <linearAlgebraInterface/src/LinearSolverParameters.hpp>
#include "SolverBase.hpp"
#include "PhysicsSolverManager.hpp"

#include "common/TimingMacros.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent )
  :
  ExecutableGroup( name, parent ),
  m_verboseLevel( 0 ),
  m_gravityVector( R1Tensor( 0.0 ) ),
  m_systemSolverParameters( groupKeyStruct::systemSolverParametersString, this ),
  m_cflFactor(),
  m_maxStableDt{ 1e99 }
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->RegisterViewWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector, false );

  // This sets a flag to indicate that this object increments time
  this->SetTimestepBehavior( 1 );

  RegisterViewWrapper( viewKeyStruct::verboseLevelString, &m_verboseLevel, false )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Verbosity level for this solver. Higher values will lead to more screen output. For non-debug "
                    " simulations, this should remain at 0." );

  RegisterViewWrapper( viewKeyStruct::cflFactorString, &m_cflFactor, false )->
    setApplyDefaultValue( 0.5 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Factor to apply to the `CFL condition <http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition>`_"
                    " when calculating the maximum allowable time step. Values should be in the interval (0,1] " );

  RegisterViewWrapper( viewKeyStruct::maxStableDtString, &m_maxStableDt, false )->
    setApplyDefaultValue( 0.5 )->
    setInputFlag( InputFlags::FALSE )->
    setDescription( "Value of the Maximum Stable Timestep for this solver." );

  this->RegisterViewWrapper( viewKeyStruct::discretizationString, &m_discretizationName, false )->
    setApplyDefaultValue( "none" )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of discretization object (defined in the :ref:`NumericalMethodsManager`) to use for this "
                    "solver. For instance, if this is a Finite Element Solver, the name of a :ref:`FiniteElement` "
                    "should be specified. If this is a Finite Volume Method, the name of a :ref:`FiniteVolume` "
                    "discretization should be specified." );

  RegisterViewWrapper( viewKeyStruct::targetRegionsString, &m_targetRegions, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Allowable regions that the solver may be applied to. Note that this does not indicate that "
                    "the solver will be applied to these regions, only that allocation will occur such that the "
                    "solver may be applied to these regions. The decision about what regions this solver will be"
                    "applied to rests in the EventManager." );

}

SolverBase::~SolverBase()
{
//  delete m_linearSolverWrapper;
}

SolverBase::CatalogInterface::CatalogType & SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

ManagedGroup * SolverBase::CreateChild( string const & childKey, string const & childName )
{
  ManagedGroup * rval = nullptr;
  if( childKey == SystemSolverParameters::CatalogName() )
  {
    rval = RegisterGroup( childName, &m_systemSolverParameters, 0 );
  }
  else
  {
    GEOS_ERROR( childKey << " is an invalid key to SolverBase child group." );
  }
  return rval;
}

void SolverBase::ExpandObjectCatalogs()
{
  CreateChild( SystemSolverParameters::CatalogName(), SystemSolverParameters::CatalogName() );
}

void SolverBase::PostProcessInput()
{
  if( this->globalGravityVector() != nullptr )
  {
    m_gravityVector = *globalGravityVector();
  }

  SetLinearSolverParameters();
}

void SolverBase::SetLinearSolverParameters()
{
  m_linearSolverParameters.verbosity = m_systemSolverParameters.verbose();

  if ( m_systemSolverParameters.scalingOption() )
  {
    m_linearSolverParameters.scaling.useRowScaling = true;
  }

  if( m_systemSolverParameters.useDirectSolver() )
  {
    m_linearSolverParameters.solverType = "direct";
  }
  else
  {
    m_linearSolverParameters.krylov.maxIterations = m_systemSolverParameters.numKrylovIter();
    m_linearSolverParameters.krylov.tolerance = m_systemSolverParameters.krylovTol();

    if ( m_systemSolverParameters.kspace() > 0 )
    {
      m_linearSolverParameters.krylov.maxRestart = m_systemSolverParameters.kspace();
    }

    if ( m_systemSolverParameters.useBicgstab() )
    {
      m_linearSolverParameters.solverType = "bicgstab";
    }
    else
    {
      m_linearSolverParameters.solverType = "gmres";
    }

    if ( m_systemSolverParameters.useMLPrecond() )
    {
      m_linearSolverParameters.preconditionerType = "amg";

      // TODO hardcoded to match old behavior
      m_linearSolverParameters.amg.cycleType = "W";
      m_linearSolverParameters.amg.smootherType = "ilu";
    }
    else
    {
      m_linearSolverParameters.preconditionerType = "ilut";
      m_linearSolverParameters.ilu.fill = static_cast<int>( m_systemSolverParameters.ilut_fill() );
      m_linearSolverParameters.ilu.threshold = m_systemSolverParameters.ilut_drop();

      // TODO hardcoded to match old behavior
      m_linearSolverParameters.dd.overlap = 1;
    }
  }
}

real64 SolverBase::SolverStep( real64 const & time_n,
                               real64 const & dt,
                               const integer cycleNumber,
                               DomainPartition * const domain )
{
  return 0;
}


void SolverBase::Execute( real64 const time_n,
                          real64 const dt,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const eventProgress,
                          ManagedGroup * const domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtRemaining = dt;

  SystemSolverParameters * const solverParams = getSystemSolverParameters();
  integer const maxSubSteps = solverParams->maxSubSteps();

  for( integer subStep = 0; subStep < maxSubSteps && dtRemaining > 0.0; ++subStep )
  {
    real64 const dtAccepted = SolverStep( time_n + (dt - dtRemaining),
                                          dtRemaining,
                                          cycleNumber,
                                          domain->group_cast<DomainPartition *>() );
    dtRemaining -= dtAccepted;

    if( m_verboseLevel >= 1 && dtRemaining > 0.0 )
    {
      GEOS_LOG_RANK_0( getName() << ": sub-step = " << subStep
                                 << ", accepted dt = " << dtAccepted
                                 << ", remaining dt = " << dtRemaining );
    }
  }

  GEOS_ERROR_IF( dtRemaining > 0.0, "Maximum allowed number of sub-steps reached. Consider increasing maxSubSteps." );
}


real64 SolverBase::LinearImplicitStep( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition * const domain,
                                       DofManager & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution )
{
  // call setup for physics solver. Pre step allocations etc.
  ImplicitStepSetup( time_n, dt, domain, dofManager, matrix, rhs, solution );

  // call assemble to fill the matrix and the rhs
  AssembleSystem( time_n, dt, domain, dofManager, matrix, rhs );

  // apply boundary conditions to system
  ApplyBoundaryConditions( time_n, dt, domain, dofManager, matrix, rhs );

  // call the default linear solver on the system
  SolveSystem( dofManager, matrix, rhs, solution );

  // apply the system solution to the fields/variables
  ApplySystemSolution( dofManager, solution, 1.0, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt, domain );

  // return the achieved timestep
  return dt;
}

bool SolverBase::LineSearch( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition * const domain,
                             DofManager const & dofManager,
                             ParallelMatrix & matrix,
                             ParallelVector & rhs,
                             ParallelVector & solution,
                             real64 & lastResidual )
{
  SystemSolverParameters * const solverParams = getSystemSolverParameters();

  integer const maxNumberLineSearchCuts = solverParams->maxLineSearchCuts();
  real64 const lineSearchCutFactor = solverParams->lineSearchCutFactor();

  // flag to determine if we should solve the system and apply the solution. If the line
  // search fails we just bail.
  bool lineSearchSuccess = false;

  real64 residualNorm = lastResidual;

  // scale factor is value applied to the previous solution. In this case we want to
  // subtract a portion of the previous solution.
  real64 scaleFactor = -1.0;

  // main loop for the line search.
  for( integer lineSearchIteration = 0; lineSearchIteration < maxNumberLineSearchCuts; ++lineSearchIteration )
  {
    // cut the scale factor by half. This means that the scale factors will
    // have values of -0.5, -0.25, -0.125, ...
    scaleFactor *= lineSearchCutFactor;

    if( !CheckSystemSolution( domain, dofManager, solution, scaleFactor ) )
    {
      if( m_verboseLevel >= 1 )
      {
        GEOS_LOG_RANK_0( "Line search: " << lineSearchIteration << ", solution check failed" );
      }
      continue;
    }

    ApplySystemSolution( dofManager, solution, scaleFactor, domain );

    // re-assemble system
    AssembleSystem( time_n, dt, domain, dofManager, matrix, rhs );

    // apply boundary conditions to system
    ApplyBoundaryConditions( time_n, dt, domain, dofManager, matrix, rhs );

    // get residual norm
    residualNorm = CalculateResidualNorm( domain, dofManager, rhs );

    if( m_verboseLevel >= 1 )
    {
      GEOS_LOG_RANK_0( "Line search: " << lineSearchIteration << ", R = " << residualNorm );
    }

    // if the residual norm is less than the last residual, we can proceed to the
    // solution step
    if( residualNorm < lastResidual )
    {
      lineSearchSuccess = true;
      break;
    }
  }

  lastResidual = residualNorm;
  return lineSearchSuccess;
}


real64 SolverBase::NonlinearImplicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          integer const cycleNumber,
                                          DomainPartition * const domain,
                                          DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs,
                                          ParallelVector & solution )
{
  // dt may be cut during the course of this step, so we are keeping a local
  // value to track the achieved dt for this step.
  real64 stepDt = dt;

  SystemSolverParameters * const solverParams = getSystemSolverParameters();

  integer const maxNewtonIter = solverParams->maxIterNewton();
  real64 const newtonTol = solverParams->newtonTol();

  integer const maxNumberDtCuts = solverParams->maxTimeStepCuts();
  real64 const dtCutFactor = solverParams->timeStepCutFactor();

  bool const allowNonConverged = solverParams->allowNonConverged() > 0;

  // a flag to denote whether we have converged
  integer isConverged = 0;

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( int dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      ResetStateToBeginningOfStep( domain );
    }

    // keep residual from previous iteration in case we need to do a line search
    real64 lastResidual = 1e99;
    integer & newtonIter = solverParams->numNewtonIterations();

    // main Newton loop
    for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
    {

      // call assemble to fill the matrix and the rhs
      AssembleSystem( time_n, stepDt, domain, dofManager, matrix, rhs );

      // apply boundary conditions to system
      ApplyBoundaryConditions( time_n, stepDt, domain, dofManager, matrix, rhs );

      // get residual norm
      real64 residualNorm = CalculateResidualNorm( domain, dofManager, rhs );

      if( m_verboseLevel >= 1 )
      {
        GEOS_LOG_RANK_0( "Attempt: " << dtAttempt << ", Newton: " << newtonIter << ", R = " << residualNorm );
      }

      // if the residual norm is less than the Newton tolerance we denote that we have
      // converged and break from the Newton loop immediately.
      if( residualNorm < newtonTol )
      {
        isConverged = 1;
        break;
      }


      // do line search in case residual has increased
      if( residualNorm > lastResidual )
      {

        residualNorm = lastResidual;
        bool lineSearchSuccess = LineSearch( time_n, stepDt, cycleNumber, domain, dofManager,
                                             matrix, rhs, solution, residualNorm );

        // if line search failed, then break out of the main Newton loop. Timestep will be cut.
        if( !lineSearchSuccess )
        {
          break;
        }
      }

      // call the default linear solver on the system
      SolveSystem( dofManager, matrix, rhs, solution );

      if( !CheckSystemSolution( domain, dofManager, solution, 1.0 ) )
      {
        // TODO try chopping (similar to line search)
        GEOS_LOG_RANK_0( "Solution check failed. Newton loop terminated." );
        break;
      }

      // apply the system solution to the fields/variables
      ApplySystemSolution( dofManager, solution, 1.0, domain );

      lastResidual = residualNorm;
    }
    if( isConverged )
    {
      // break out of outer loop
      break;
    }
    else
    {
      // cut timestep, go back to beginning of step and restart the Newton loop
      stepDt *= dtCutFactor;
    }
  }

  if( !isConverged )
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

real64 SolverBase::ExplicitStep( real64 const & time_n,
                                 real64 const & dt,
                                 integer const cycleNumber,
                                 DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::ExplicitStep called!. Should be overridden." );
  return 0;
}

void SolverBase::ImplicitStepSetup( real64 const & time_n,
                                    real64 const & dt,
                                    DomainPartition * const domain,
                                    DofManager & dofManager,
                                    ParallelMatrix & matrix,
                                    ParallelVector & rhs,
                                    ParallelVector & solution )
{
  GEOS_ERROR( "SolverBase::ImplicitStepSetup called!. Should be overridden." );
}

void SolverBase::AssembleSystem( real64 const time,
                                 real64 const dt,
                                 DomainPartition * const domain,
                                 DofManager const & dofManager,
                                 ParallelMatrix & matrix,
                                 ParallelVector & rhs )
{
  GEOS_ERROR( "SolverBase::Assemble called!. Should be overridden." );
}

void SolverBase::ApplyBoundaryConditions( real64 const time,
                                          real64 const dt,
                                          DomainPartition * const domain,
                                          DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs )
{
  GEOS_ERROR( "SolverBase::SolveSystem called!. Should be overridden." );
}

real64
SolverBase::CalculateResidualNorm( DomainPartition const * const domain,
                                   DofManager const & dofManager,
                                   ParallelVector const & rhs )
{
  GEOS_ERROR( "SolverBase::CalculateResidualNorm called!. Should be overridden." );
  return 0;
}

void SolverBase::SolveSystem( DofManager const & dofManager,
                              ParallelMatrix & matrix,
                              ParallelVector & rhs,
                              ParallelVector & solution )
{
  // Create a solver from the parameter list
  LinearSolver solver( m_linearSolverParameters );

  // Solve using the iterative solver and compare norms with true solution
  solver.solve( matrix, solution, rhs );
}

bool SolverBase::CheckSystemSolution( DomainPartition const * const domain,
                                      DofManager const & dofManager,
                                      ParallelVector const & solution,
                                      real64 const scalingFactor )
{
  return true;
}

void SolverBase::ApplySystemSolution( DofManager const & dofManager,
                                      ParallelVector const & solution,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::ApplySystemSolution called!. Should be overridden." );
}

void SolverBase::ResetStateToBeginningOfStep( DomainPartition * const )
{
  GEOS_ERROR( "SolverBase::ResetStateToBeginningOfStep called!. Should be overridden." );
}

void SolverBase::ImplicitStepComplete( real64 const & time,
                                       real64 const & dt,
                                       DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::ImplicitStepComplete called!. Should be overridden." );
}

R1Tensor const * SolverBase::globalGravityVector() const
{
  R1Tensor const * rval = nullptr;
  if( getParent()->getName() == "Solvers" )
  {
    rval = &(getParent()->getReference<R1Tensor>( viewKeyStruct::gravityVectorString ));
  }

  return rval;
}


} /* namespace ANST */
