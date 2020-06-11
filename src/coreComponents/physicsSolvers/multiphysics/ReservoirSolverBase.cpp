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

ReservoirSolverBase::ReservoirSolverBase( const std::string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_flowSolverName(),
  m_wellSolverName()
{
  registerWrapper( viewKeyStruct::flowSolverNameString, &m_flowSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the flow solver to use in the reservoir-well system solver" );

  registerWrapper( viewKeyStruct::wellSolverNameString, &m_wellSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the well solver to use in the reservoir-well system solver" );

}

ReservoirSolverBase::~ReservoirSolverBase()
{}

void ReservoirSolverBase::PostProcessInput()
{
  SolverBase::PostProcessInput();

  m_flowSolver = this->getParent()->GetGroup< FlowSolverBase >( m_flowSolverName );
  m_wellSolver = this->getParent()->GetGroup< WellSolverBase >( m_wellSolverName );

  GEOSX_ERROR_IF( m_flowSolver == nullptr, "Flow solver not found or invalid type: " << m_flowSolverName );
  GEOSX_ERROR_IF( m_wellSolver == nullptr, "Well solver not found or invalid type: " << m_wellSolverName );

  m_wellSolver->SetFlowSolverName( m_flowSolverName );
  m_flowSolver->setReservoirWellsCoupling();
}

void ReservoirSolverBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion & subRegion )
  {
    // get the string to access the permeability
    string const permeabilityKey = FlowSolverBase::viewKeyStruct::permeabilityString;

    PerforationData * const perforationData = subRegion.GetPerforationData();

    // compute the Peaceman index (if not read from XML)
    perforationData->ComputeWellTransmissibility( *meshLevel,
                                                  &subRegion,
                                                  permeabilityKey );
  } );

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

void ReservoirSolverBase::SetupDofs( DomainPartition const * const domain,
                                     DofManager & dofManager ) const
{
  m_flowSolver->SetupDofs( domain, dofManager );
  m_wellSolver->SetupDofs( domain, dofManager );
  // TODO: add coupling when dofManager can support perforation connectors
}

void ReservoirSolverBase::SetupSystem( DomainPartition * const domain,
                                       DofManager & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution,
                                       bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalDof = dofManager.numLocalDofs();

  matrix.createWithLocalSize( numLocalDof, numLocalDof, 8, MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );

  if( setSparsity )
  {
    dofManager.setSparsityPattern( matrix, false ); // don't close the matrix

    // by hand, add sparsity pattern induced by well perforations
    AddCouplingSparsityPattern( domain,
                                dofManager,
                                matrix );
  }
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


void ReservoirSolverBase::AssembleSystem( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition * const domain,
                                          DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs )
{

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

  matrix.open();
  rhs.open();

  // assemble perforation rates in J_WR, J_RW, J_RR and J_WW
  AssembleCouplingTerms( time_n, dt, domain,
                         &dofManager,
                         &matrix,
                         &rhs );

  matrix.close();
  rhs.close();

  if( getLogLevel() >= 2 )
  {
    GEOSX_LOG_RANK_0( "After ReservoirSolverBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, LAIOutputFormat::MATRIX_MARKET );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, LAIOutputFormat::MATRIX_MARKET );

    GEOSX_LOG_RANK_0( "After ReservoirSolverBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
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
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();
  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
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

void ReservoirSolverBase::ImplicitStepComplete( real64 const & time_n,
                                                real64 const & dt,
                                                DomainPartition * const domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_wellSolver->ImplicitStepComplete( time_n, dt, domain );
}

void ReservoirSolverBase::ResetViews( DomainPartition * const GEOSX_UNUSED_PARAM( domain ) )
{}


} /* namespace geosx */
