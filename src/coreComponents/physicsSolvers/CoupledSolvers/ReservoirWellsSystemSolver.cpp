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
 * @file ReservoirWellsSystemSolver.cpp
 *
 */


#include "ReservoirWellsSystemSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "../FiniteVolume/FlowSolverBase.hpp"
#include "../Wells/WellSolverBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "wells/WellManager.hpp"
#include "wells/Well.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;
  
ReservoirWellsSystemSolver::ReservoirWellsSystemSolver( const std::string& name,
                                                        ManagedGroup * const parent ):
  SolverBase(name,parent),
  m_flowSolverName(),
  m_wellSolverName()
{
  RegisterViewWrapper(viewKeyStruct::flowSolverNameString, &m_flowSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the flow solver to use in the reservoir-well system solver");

  RegisterViewWrapper(viewKeyStruct::wellSolverNameString, &m_wellSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the well solver to use in the reservoir-well system solver");

}


void ReservoirWellsSystemSolver::ImplicitStepSetup( real64 const& time_n,
                                                    real64 const& dt,
                                                    DomainPartition * const domain,
                                                    systemSolverInterface::EpetraBlockSystem * const blockSystem)
{
  std::cout << "ReservoirWellsSystemSolver::ImplicitStepSetup started" << std::endl;
  
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());
  
  // setup the individual solvers
  flowSolver.ImplicitStepSetup( time_n, dt, domain, blockSystem );
  wellSolver.ImplicitStepSetup( time_n, dt, domain, blockSystem );

  // setup the coupled linear system
  SetupSystem( domain, blockSystem );

  std::cout << "ReservoirWellsSystemSolver::ImplicitStep complete" << std::endl;
}

  
void ReservoirWellsSystemSolver::SetupSystem ( DomainPartition * const domain,
                                               systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  std::cout << "ReservoirWellsSystemSolver::SetupSystem started" << std::endl;
  
  // get the reservoir and well solvers
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());
  
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();
  WellManager * const wellManager = domain->getWellManager();
  
  // for this solver, the dof are on the cell center, and the row corrosponds to an element
  localIndex  numResGhostRows   = 0;
  localIndex  numWellGhostRows  = 0;
  localIndex  numResLocalRows   = 0;
  localIndex  numWellLocalRows  = 0;
  globalIndex numResGlobalRows  = 0;
  globalIndex numWellGlobalRows = 0;
  
  // get the number of local reservoir elements, and ghost elements...i.e. local rows and ghost rows
  elementRegionManager->forElementSubRegions( [&]( ObjectManagerBase * const subRegion )
  {
    localIndex subRegionGhosts = subRegion->GetNumberOfGhosts();
    numResGhostRows += subRegionGhosts;
    numResLocalRows += subRegion->size() - subRegionGhosts;
  });
  
  flowSolver.SetNumRowsAndTrilinosIndices( mesh,
                                           numResLocalRows,
                                           numResGlobalRows,
                                           0 );
  /*
  // TODO: uncomment this
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeysSinglePhaseFlow.blockLocalDofNumber.Key());
  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
  */

  std::cout << "SinglePhaseFlow::SetNumRowsAndTrilinosIndices complete" << std::endl;

  // get the number of local well elements, and ghost elements...i.e. local rows and ghost rows
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    // well elements
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
    localIndex subRegionGhosts = wellElementSubRegion->GetNumberOfGhosts();
    numWellGhostRows += subRegionGhosts;
    numWellLocalRows += wellElementSubRegion->size() - subRegionGhosts;
  });

  wellSolver.SetNumRowsAndTrilinosIndices( domain,
                                           numWellLocalRows,
                                           numWellGlobalRows,
					   numResGlobalRows); // offset for well global rows

  // TODO: synchronize the well vars
  
  std::cout << "SinglePhaseWell::SetNumRowsAndTrilinosIndices complete" << std::endl;
  
  // construct row map, and set a pointer to the row map
  globalIndex const nResDOF  = flowSolver.numDofPerCell();
  globalIndex const nWellDOF = wellSolver.numDofPerElement()
                             + wellSolver.numDofPerConnection();
  globalIndex const totalNumberOfDofs = numResGlobalRows  * nResDOF   // dofs in reservoir
                                      + numWellGlobalRows * nWellDOF; // dofs in wells
  Epetra_Map * const
  rowMap = blockSystem->
           SetRowMap( BlockIDs::fluidPressureBlock,
                      std::make_unique<Epetra_Map>( totalNumberOfDofs,
                                                    totalNumberOfDofs,
                                                    0,
                                                    m_linearSolverWrapper.m_epetraComm ) );
  
  // construct sparsity matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const
  sparsity = blockSystem->SetSparsity( BlockIDs::fluidPressureBlock,
                                       BlockIDs::fluidPressureBlock,
                                       std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );

  // setup sparsity pattern for J_RR
  flowSolver.SetSparsityPattern( domain, sparsity );
  // setep sparsity pattern for J_WW, J_WR, and J_RW
  wellSolver.SetSparsityPattern( domain, sparsity, numResGlobalRows, nResDOF ); 
  
  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( BlockIDs::fluidPressureBlock,
                          BlockIDs::fluidPressureBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
								BlockIDs::fluidPressureBlock );
  
  // construct solution vector
  blockSystem->SetSolutionVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  std::cout << "ReservoirWellsSystemSolver::SetupSystem complete" << std::endl;
}
  
void ReservoirWellsSystemSolver::AssembleSystem( DomainPartition * const domain,
                                                 systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                 real64 const time,
                                                 real64 const dt )
{
  std::cout << "ReservoirWellsSystemSolver::AssembleSystem started" << std::endl;
  
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  // set Jacobian and residual to zero
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
								BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  flowSolver.AssembleSystem( domain, blockSystem, time, dt );
  wellSolver.AssembleSystem( domain, blockSystem, time, dt );

  jacobian->GlobalAssemble( true );
  residual->GlobalAssemble();

  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After ReservoirWellsSystemSolver::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }
  
  std::cout << "ReservoirWellsSystemSolver::AssembleSystem complete" << std::endl;
}

void ReservoirWellsSystemSolver::ApplyBoundaryConditions( DomainPartition * const domain,
                                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                          real64 const time,
                                                          real64 const dt )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());

  flowSolver.ApplyBoundaryConditions( domain, blockSystem, time, dt );
  // no boundary conditions for wells
}

real64 ReservoirWellsSystemSolver::CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
                                                          DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  real64 const reservoirResidualNorm = flowSolver.CalculateResidualNorm( blockSystem, domain );
  real64 const wellResidualNorm      = wellSolver.CalculateResidualNorm( blockSystem, domain );
  
  return ( reservoirResidualNorm > wellResidualNorm )
         ? reservoirResidualNorm
         : wellResidualNorm;
}

void ReservoirWellsSystemSolver::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                              SystemSolverParameters const * const params )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());

  // for now, the flow solver is responsible for solving the linear system
  flowSolver.SolveSystem( blockSystem, params );  
}

bool ReservoirWellsSystemSolver::CheckSystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                                      real64 const scalingFactor,
                                                      DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  bool const validReservoirSolution = flowSolver.CheckSystemSolution( blockSystem, scalingFactor, domain );
  bool const validWellSolution      = wellSolver.CheckSystemSolution( blockSystem, scalingFactor, domain );
  
  return ( validReservoirSolution && validWellSolution );
}

void ReservoirWellsSystemSolver::ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                                      real64 const scalingFactor,
                                                      DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  flowSolver.ApplySystemSolution( blockSystem, scalingFactor, domain );
  wellSolver.ApplySystemSolution( blockSystem, scalingFactor, domain );
}

void ReservoirWellsSystemSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  flowSolver.ResetStateToBeginningOfStep( domain );
  wellSolver.ResetStateToBeginningOfStep( domain ); 
}
  
void ReservoirWellsSystemSolver::ImplicitStepComplete( real64 const& time_n,
                                                       real64 const& dt,
                                                       DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  flowSolver.ImplicitStepComplete( time_n, dt, domain );
  wellSolver.ImplicitStepComplete( time_n, dt, domain );
}

void ReservoirWellsSystemSolver::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const problemManager)
{
  this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>()->setReservoirWellsSystemCoupling();
}

ReservoirWellsSystemSolver::~ReservoirWellsSystemSolver()
{
  // TODO Auto-generated destructor stub
}

real64 ReservoirWellsSystemSolver::SolverStep( real64 const & time_n,
                                               real64 const & dt,
                                               int const cycleNumber,
                                               DomainPartition * domain )
{
  real64 dt_return = dt;

  ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

  // currently the only method is implicit time integration
  dt_return= this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          getLinearSystemRepository() );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

  
REGISTER_CATALOG_ENTRY( SolverBase, ReservoirWellsSystemSolver, std::string const &, ManagedGroup * const )

} /* namespace geosx */
