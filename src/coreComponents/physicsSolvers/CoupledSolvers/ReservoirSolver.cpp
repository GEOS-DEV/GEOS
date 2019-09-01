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
 * @file ReservoirSolver.cpp
 *
 */


#include "ReservoirSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "physicsSolvers/FiniteVolume/FlowSolverBase.hpp"
#include "physicsSolvers/Wells/WellSolverBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
  
ReservoirSolver::ReservoirSolver( const std::string& name,
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

ReservoirSolver::~ReservoirSolver()
{
}

void ReservoirSolver::PostProcessInput()
{
  SolverBase::PostProcessInput();

  m_flowSolver = this->getParent()->GetGroup<FlowSolverBase>( m_flowSolverName );
  m_wellSolver = this->getParent()->GetGroup<WellSolverBase>( m_wellSolverName );

  GEOS_ERROR_IF( m_flowSolver == nullptr, "Flow solver not found or invalid type: " << m_flowSolverName );
  GEOS_ERROR_IF( m_wellSolver == nullptr, "Well solver not found or invalid type: " << m_wellSolverName );

  m_wellSolver->SetFlowSolverName( m_flowSolverName );
  m_flowSolver->setReservoirWellsCoupling();
}

real64 ReservoirSolver::SolverStep( real64 const & time_n,
                                    real64 const & dt,
                                    int const cycleNumber,
                                    DomainPartition * domain )
{
  real64 dt_return = dt;

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

void ReservoirSolver::ImplicitStepSetup( real64 const & time_n,
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

  // setup the coupled linear system
  SetupSystem( domain, dofManager, matrix, rhs, solution );
}

void ReservoirSolver::SetupDofs( DomainPartition const * const domain,
                                 DofManager & dofManager ) const
{
  m_flowSolver->SetupDofs( domain, dofManager );
  m_wellSolver->SetupDofs( domain, dofManager );

  // TODO: add coupling when dofManager can support perforation connectors
}

void ReservoirSolver::SetupSystem( DomainPartition * const domain,
                                   DofManager & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
  dofManager.close();

  dofManager.setSparsityPattern( matrix, "", "", false ); // don't close the matrix
  dofManager.setVector( rhs );
  dofManager.setVector( solution );

  // TODO: remove this and just call SolverBase::SetupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  string const resDofKey  = dofManager.getKey( m_wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager.getKey( m_wellSolver->WellElementDofName() );

  localIndex const resNDOF = m_wellSolver->NumDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->NumDofPerWellElement();

  // TODO: overallocating in single phase flow case, find a way around
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( resDofKey );

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {
    PerforationData const * const perforationData = subRegion->GetPerforationData();

    // get the well degrees of freedom and ghosting info
    arrayView1d< globalIndex > const & wellElemDofNumber =
      subRegion->getReference< array1d<globalIndex> >( wellDofKey );

    // get the well element indices corresponding to each perforation
    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::wellElementIndexString );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    stackArray1d< globalIndex, maxNumDof > dofIndexRes( resNDOF );
    stackArray1d< globalIndex, maxNumDof > dofIndexWell( wellNDOF );
    stackArray2d< real64, maxNumDof * maxNumDof > values( resNDOF, wellNDOF );
    values = 1.0;

    // Insert the entries corresponding to reservoir-well perforations
    // This will fill J_WR, and J_RW
    for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];
      localIndex const iwelem = perfWellElemIndex[iperf];

      for( localIndex idof = 0; idof < resNDOF; ++idof )
      {
        dofIndexRes[idof] = resDofNumber[er][esr][ei] + idof;
      }

      for( localIndex idof = 0; idof < wellNDOF; ++idof )
      {
        dofIndexWell[idof] = wellElemDofNumber[iwelem] + idof;
      }

      // fill J_RW
      matrix.insert( dofIndexRes.data(),
                     dofIndexWell.data(),
                     values.data(),
                     resNDOF,
                     wellNDOF );

      // fill J_WR
      matrix.insert( dofIndexWell.data(),
                     dofIndexRes.data(),
                     values.data(),
                     wellNDOF,
                     resNDOF );
    }

   } );

  matrix.close();
}

void ReservoirSolver::AssembleSystem( real64 const time_n,
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
  // assemble J_WW, J_WR, J_RW + perforation rates in J_RR
  m_wellSolver->AssembleSystem( time_n, dt, domain,
                                dofManager,
                                matrix,
                                rhs );

  matrix.close();
  rhs.close();

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After ReservoirSolver::AssembleSystem" );
    GEOS_LOG_RANK_0("\nJacobian:\n");
    std::cout << matrix;
    GEOS_LOG_RANK_0("\nResidual:\n");
    std::cout << rhs;
  }

  if( verboseLevel() >= 3 )
  {
    SystemSolverParameters * const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOS_LOG_RANK_0( "After ReservoirSolver::AssembleSystem" );
    GEOS_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOS_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void ReservoirSolver::ApplyBoundaryConditions( real64 const time_n,
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

real64 ReservoirSolver::CalculateResidualNorm( DomainPartition const * const domain,
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

void ReservoirSolver::SolveSystem( DofManager const & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0("After SinglePhaseFlow::SolveSystem");
    GEOS_LOG_RANK_0("\nSolution:\n");
    std::cout << solution;
  }
}

bool ReservoirSolver::CheckSystemSolution( DomainPartition const * const domain,
                                           DofManager const & dofManager,
                                           ParallelVector const & solution,
                                           real64 const scalingFactor )
{
  bool const validReservoirSolution = m_flowSolver->CheckSystemSolution( domain, dofManager, solution, scalingFactor );
  bool const validWellSolution      = m_wellSolver->CheckSystemSolution( domain, dofManager, solution, scalingFactor );

  return ( validReservoirSolution && validWellSolution );
}

void ReservoirSolver::ApplySystemSolution( DofManager const & dofManager,
                                           ParallelVector const & solution,
                                           real64 const scalingFactor,
                                           DomainPartition * const domain )
{
  // update the reservoir variables
  m_flowSolver->ApplySystemSolution( dofManager, solution, scalingFactor, domain );
  // update the well variables
  m_wellSolver->ApplySystemSolution( dofManager, solution, scalingFactor, domain );
}

void ReservoirSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // reset reservoir variables
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  // reset well variables
  m_wellSolver->ResetStateToBeginningOfStep( domain );
}

void ReservoirSolver::ImplicitStepComplete( real64 const& time_n,
                                            real64 const& dt,
                                            DomainPartition * const domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_wellSolver->ImplicitStepComplete( time_n, dt, domain );
}


REGISTER_CATALOG_ENTRY( SolverBase, ReservoirSolver, std::string const &, Group * const )

} /* namespace geosx */
