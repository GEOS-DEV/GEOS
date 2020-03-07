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

/**
 * @file ReservoirSolverBase.cpp
 *
 */


#include "ReservoirSolverBase.hpp"

#include "common/TimingMacros.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/WellSolverBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
  
ReservoirSolverBase::ReservoirSolverBase( const std::string& name,
                                          Group * const parent ):
  SolverBase(name,parent),
  m_flowSolverName(),
  m_wellSolverName()
{
  registerWrapper(viewKeyStruct::flowSolverNameString, &m_flowSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the flow solver to use in the reservoir-well system solver");

  registerWrapper(viewKeyStruct::wellSolverNameString, &m_wellSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the well solver to use in the reservoir-well system solver");

}

ReservoirSolverBase::~ReservoirSolverBase()
{
}

void ReservoirSolverBase::PostProcessInput()
{
  SolverBase::PostProcessInput();

  m_flowSolver = this->getParent()->GetGroup<FlowSolverBase>( m_flowSolverName );
  m_wellSolver = this->getParent()->GetGroup<WellSolverBase>( m_wellSolverName );

  GEOSX_ERROR_IF( m_flowSolver == nullptr, "Flow solver not found or invalid type: " << m_flowSolverName );
  GEOSX_ERROR_IF( m_wellSolver == nullptr, "Well solver not found or invalid type: " << m_wellSolverName );

  m_wellSolver->SetFlowSolverName( m_flowSolverName );
  m_flowSolver->setReservoirWellsCoupling();
}

void ReservoirSolverBase::InitializePostInitialConditions_PreSubGroups(Group * const rootGroup)
{
  SolverBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  // bind the stored reservoir views to the current domain
  ResetViews( domain );

}


real64 ReservoirSolverBase::SolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition * domain )
{
  real64 dt_return = dt;

  // setup the coupled linear system
  SetupSystem( domain, m_dofManager, m_matrix, m_rhs, m_solution );
  
  // setup reservoir and well systems
  ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // currently the only method is implicit time integration
  dt_return = this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain,
                                           m_dofManager,
                                           m_matrix,
                                           m_rhs,
                                           m_solution );

  // complete time step in reservoir and well systems
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void ReservoirSolverBase::ImplicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition * const domain,
                                             DofManager & dofManager,
                                             ParallelMatrix & matrix,
                                             ParallelVector & rhs,
                                             ParallelVector & solution )
{
  // setup the individual solvers
  m_flowSolver->ImplicitStepSetup( time_n, dt, domain,
                                   dofManager,
                                   matrix,
                                   rhs,
                                   solution );
  m_wellSolver->ImplicitStepSetup( time_n, dt, domain,
                                   dofManager,
                                   matrix,
                                   rhs,
                                   solution );

}

void ReservoirSolverBase::SetupDofs( DomainPartition const * const GEOSX_UNUSED_PARAM( domain ),
                                     DofManager & GEOSX_UNUSED_PARAM( dofManager ) ) const
{
  GEOSX_ERROR( "ReservoirSolverBase::SetupDofs called!. Should be overridden." );
}

void ReservoirSolverBase::SetupSystem( DomainPartition * const GEOSX_UNUSED_PARAM( domain ),
                                       DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                       ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                       ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                       ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  GEOSX_ERROR( "ReservoirSolverBase::SetupSystem called!. Should be overridden." );
}

void ReservoirSolverBase::AssembleSystem( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition * const domain,
                                          DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs )
{
  // open() and zero() are called from the flow solver

  // assemble J_RR (excluding perforation rates)
  m_flowSolver->AssembleSystem( time_n, dt, domain,
                                dofManager,
                                matrix,
                                rhs );
  // assemble J_WW (excluding perforation rates)
  m_wellSolver->AssembleSystem( time_n, dt, domain,
                                dofManager,
                                matrix,
                                rhs );

  // assemble perforation rates in J_WR, J_RW, J_RR and J_WW
  AssembleCouplingTerms( time_n, dt, domain,
                         &dofManager,
                         &matrix,
                         &rhs );

  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After ReservoirSolverBase::AssembleSystem" );
    GEOSX_LOG_RANK_0("\nJacobian:\n");
    std::cout << matrix;
    GEOSX_LOG_RANK_0("\nResidual:\n");
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOSX_LOG_RANK_0( "After ReservoirSolverBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void ReservoirSolverBase::AssembleCouplingTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                                 DomainPartition * const GEOSX_UNUSED_PARAM( domain ),
                                                 DofManager const * const GEOSX_UNUSED_PARAM( dofManager ),
                                                 ParallelMatrix * const GEOSX_UNUSED_PARAM( matrix ),
                                                 ParallelVector * const GEOSX_UNUSED_PARAM( rhs ) )
{
  GEOSX_ERROR( "ReservoirSolverBase::AssembleCouplingTerms called!. Should be overridden." );
}

void ReservoirSolverBase::ApplyBoundaryConditions( real64 const time_n,
                                                   real64 const dt,
                                                   DomainPartition * const domain,
                                                   DofManager const & dofManager,
                                                   ParallelMatrix & matrix,
                                                   ParallelVector & rhs )
{
  m_flowSolver->ApplyBoundaryConditions( time_n, dt, domain,
                                         dofManager,
                                         matrix,
                                         rhs );
  // no boundary conditions for wells
}

real64 ReservoirSolverBase::CalculateResidualNorm( DomainPartition const * const domain,
                                                   DofManager const & dofManager,
                                                   ParallelVector const & rhs )
{
  // compute norm of reservoir equations residuals
  real64 const reservoirResidualNorm = m_flowSolver->CalculateResidualNorm( domain, dofManager, rhs );
  // compute norm of well equations residuals
  real64 const wellResidualNorm      = m_wellSolver->CalculateResidualNorm( domain, dofManager, rhs );

  return sqrt( reservoirResidualNorm*reservoirResidualNorm
             + wellResidualNorm*wellResidualNorm );
}

void ReservoirSolverBase::SolveSystem( DofManager const & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution )
{
  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0("After SinglePhaseFlow::SolveSystem");
    GEOSX_LOG_RANK_0("\nSolution:\n");
    std::cout << solution;
  }
}

bool ReservoirSolverBase::CheckSystemSolution( DomainPartition const * const domain,
                                               DofManager const & dofManager,
                                               ParallelVector const & solution,
                                               real64 const scalingFactor )
{
  bool const validReservoirSolution = m_flowSolver->CheckSystemSolution( domain, dofManager, solution, scalingFactor );
  bool const validWellSolution      = m_wellSolver->CheckSystemSolution( domain, dofManager, solution, scalingFactor );

  return ( validReservoirSolution && validWellSolution );
}

void ReservoirSolverBase::ApplySystemSolution( DofManager const & dofManager,
                                               ParallelVector const & solution,
                                               real64 const scalingFactor,
                                               DomainPartition * const domain )
{
  // update the reservoir variables
  m_flowSolver->ApplySystemSolution( dofManager, solution, scalingFactor, domain );
  // update the well variables
  m_wellSolver->ApplySystemSolution( dofManager, solution, scalingFactor, domain );
}

void ReservoirSolverBase::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // reset reservoir variables
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  // reset well variables
  m_wellSolver->ResetStateToBeginningOfStep( domain );
}

void ReservoirSolverBase::ImplicitStepComplete( real64 const& time_n,
                                                real64 const& dt,
                                                DomainPartition * const domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_wellSolver->ImplicitStepComplete( time_n, dt, domain );
}

void ReservoirSolverBase::ResetViews( DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{
}


} /* namespace geosx */
