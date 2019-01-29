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

#include "SolverBase.hpp"
#include "PhysicsSolverManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ):
  ExecutableGroup( name, parent ),
  m_linearSolverWrapper(),
  m_verboseLevel( 0 ),
  m_gravityVector( R1Tensor( 0.0 ) ),
  m_systemSolverParameters( groupKeyStruct::systemSolverParametersString, this ),
  m_cflFactor(),
  m_maxStableDt{1e99}
{

  this->RegisterViewWrapper( viewKeyStruct::verboseLevelString, &m_verboseLevel, 0 );
  this->RegisterViewWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector, 0 );
//  this->RegisterViewWrapper( viewKeyStruct::blockLocalDofNumberString, &m_blockLocalDofNumber, 0 );

  // This sets a flag to indicate that this object increments time
  this->SetTimestepBehavior(1);


  RegisterViewWrapper(viewKeyStruct::verboseLevelString, &m_verboseLevel, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Verbosity level");

  RegisterViewWrapper(viewKeyStruct::cflFactorString, &m_cflFactor, false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Factor to apply to the CFL condition when calculating the maximum allowable time step. "
          "Values should be in the interval (0,1] ");

  RegisterViewWrapper(viewKeyStruct::maxStableDtString, &m_maxStableDt, false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Factor to apply to the CFL condition when calculating the maximum allowable time step. "
          "Values should be in the interval (0,1] ");

  this->RegisterViewWrapper( viewKeyStruct::discretizationString, &m_discretizationName, false )->
    setApplyDefaultValue("none")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of discretization object to use for this solver.");

  RegisterViewWrapper(viewKeyStruct::targetRegionsString, &m_targetRegions, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Allowable regions that the solver may be applied to. Note that this does not indicate that "
                   "the solver will be applied to these regions, only that allocation will occur such that the "
                   "solver may be applied to these regions. The decision about what regions this solver will be"
                   "applied to rests in the EventManager.");

}

SolverBase::~SolverBase()
{
//  delete m_linearSolverWrapper;
}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

ManagedGroup * SolverBase::CreateChild( string const & childKey, string const & childName )
{
  ManagedGroup * rval = nullptr;
  if( childKey==SystemSolverParameters::CatalogName() )
  {
    rval = RegisterGroup( childName, &m_systemSolverParameters, 0 );
  }
  else
  {
    GEOS_ERROR(childKey<<" is an invalid key to SolverBase child group.");
  }
  return rval;
}


void SolverBase::PostProcessInput()
{
  if( this->globalGravityVector() != nullptr )
  {
    m_gravityVector=*globalGravityVector();
  }
}


real64 SolverBase::SolverStep( real64 const& time_n,
                               real64 const& dt,
                               const integer cycleNumber,
                               DomainPartition * domain )
{
  return 0;
}


void SolverBase::Execute( real64 const& time_n,
                          real64 const& dt,
                          integer const cycleNumber,
                          integer const eventCounter,
                          real64 const & eventProgress,
                          ManagedGroup * domain )
{
  if( dt > 0 )
  {
    SolverStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>());
  }
}


real64 SolverBase::LinearImplicitStep( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition * const domain,
                                       systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  // call setup for physics solver. Pre step allocations etc.
  ImplicitStepSetup( time_n, dt, domain, blockSystem );

  // call assemble to fill the matrix and the rhs
  AssembleSystem( domain, blockSystem, time_n+dt, dt );

  // apply boundary conditions to system
  ApplyBoundaryConditions( domain, blockSystem, time_n, dt );

  // call the default linear solver on the system
  SolveSystem( blockSystem,
               getSystemSolverParameters() );

  // apply the system solution to the fields/variables
  ApplySystemSolution( blockSystem, 1.0, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt, domain );

  // return the achieved timestep
  return dt;
}

bool SolverBase::LineSearch( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition * const domain,
                             systemSolverInterface::EpetraBlockSystem * const blockSystem,
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

    if( !CheckSystemSolution( blockSystem, scaleFactor, domain ) )
    {
      if( m_verboseLevel >= 1 )
      {
        GEOS_LOG_RANK_0( "Line search: " << lineSearchIteration << ", solution check failed" );
      }
      continue;
    }

    ApplySystemSolution( blockSystem, scaleFactor, domain );

    // re-assemble system
    AssembleSystem( domain, blockSystem, time_n, dt );

    // apply boundary conditions to system
    ApplyBoundaryConditions( domain, blockSystem, time_n, dt );

    // get residual norm
    residualNorm = CalculateResidualNorm( blockSystem, domain );

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
                                          systemSolverInterface::EpetraBlockSystem * const blockSystem )
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
  for( int dtAttempt = 0 ; dtAttempt<maxNumberDtCuts ; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    ResetStateToBeginningOfStep( domain );

    // main Newton loop
    // keep residual from previous iteration in case we need to do a line search
    real64 lastResidual = 1e99;
    integer & k = solverParams->numNewtonIterations();
    for( k=0 ; k<maxNewtonIter ; ++k )
    {

      // call assemble to fill the matrix and the rhs
      AssembleSystem( domain, blockSystem, time_n, stepDt );

      // apply boundary conditions to system
      ApplyBoundaryConditions( domain, blockSystem, time_n, stepDt );

      // get residual norm
      real64 residualNorm = CalculateResidualNorm( blockSystem, domain );

      if ( m_verboseLevel >= 1 )
      {
        GEOS_LOG_RANK_0( "Attempt: " << dtAttempt  << ", Newton: " << k << ", R = " << residualNorm );
      }

      // if the residual norm is less than the Newton tolerance we denote that we have
      // converged and break from the Newton loop immediately.
      if ( residualNorm < newtonTol )
      {
        isConverged = 1;
        break;
      }


      // do line search in case residual has increased
      if( residualNorm > lastResidual )
      {

        residualNorm = lastResidual;
        bool lineSearchSuccess = LineSearch( time_n, stepDt, cycleNumber, domain, blockSystem, residualNorm );

        // if line search failed, then break out of the main Newton loop. Timestep will be cut.
        if( !lineSearchSuccess )
        {
          break;
        }
      }

      // call the default linear solver on the system
      SolveSystem( blockSystem, getSystemSolverParameters() );

      if ( !CheckSystemSolution( blockSystem, 1.0, domain ) )
      {
        // TODO try chopping (similar to line search)
        GEOS_LOG_RANK_0( "Solution check failed. Newton loop terminated." );
        break;
      }

      // apply the system solution to the fields/variables
      ApplySystemSolution( blockSystem, 1.0, domain );

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

  if ( !isConverged )
  {
    GEOS_LOG_RANK_0( "Convergence not achieved." );

    if ( allowNonConverged )
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

void SolverBase::ImplicitStepSetup( real64 const& time_n,
                                    real64 const& dt,
                                    DomainPartition * const domain,
                                    systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  GEOS_ERROR( "SolverBase::ImplicitStepSetup called!. Should be overridden." );
}

void SolverBase::AssembleSystem( DomainPartition * const domain,
                                 systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                 real64 const time,
                                 real64 const dt )
{
  GEOS_ERROR( "SolverBase::Assemble called!. Should be overridden." );
}

void SolverBase::ApplyBoundaryConditions( DomainPartition * const domain,
                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                          real64 const time,
                                          real64 const dt )
{
  GEOS_ERROR( "SolverBase::SolveSystem called!. Should be overridden." );
}

real64
SolverBase::
CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                       DomainPartition * const domain )
{
  GEOS_ERROR( "SolverBase::CalculateResidualNorm called!. Should be overridden." );
  return 0;
}

void SolverBase::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params )
{
  GEOS_ERROR( "SolverBase::SolveSystem called!. Should be overridden." );
}

void SolverBase::ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const )
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

void SolverBase::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                              SystemSolverParameters const * const params,
                              systemSolverInterface::BlockIDs const blockID )
{
  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( blockID );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( blockID );
  residual->Scale( -1.0 );

  solution->Scale( 0.0 );

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                blockID );

  if( verboseLevel() >= 2 )
  {
    solution->Print( std::cout );
  }

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

systemSolverInterface::EpetraBlockSystem const * SolverBase::getLinearSystemRepository() const
{
  return &( getParent()->
            getReference<systemSolverInterface::EpetraBlockSystem>( PhysicsSolverManager::
                                                                    viewKeyStruct::
                                                                    blockSystemRepositoryString ) );
}

systemSolverInterface::EpetraBlockSystem * SolverBase::getLinearSystemRepository()
{
  return &( getParent()->
            getReference<systemSolverInterface::EpetraBlockSystem>( PhysicsSolverManager::
                                                                    viewKeyStruct::
                                                                    blockSystemRepositoryString ) );
}

bool SolverBase::CheckSystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor, DomainPartition * const domain )
{
  return true;
}


} /* namespace ANST */
