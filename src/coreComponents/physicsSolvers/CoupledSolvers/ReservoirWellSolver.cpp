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
 * @file ReservoirWellSolver.cpp
 *
 */


#include "ReservoirWellSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "../FiniteVolume/FlowSolverBase.hpp"
#include "../../wells/WellSolverBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ReservoirWellSolver::ReservoirWellSolver( const std::string& name,
                                          ManagedGroup * const parent ):
  SolverBase(name,parent),
  m_flowSolverName(),
  m_wellSolverName()
{
  RegisterViewWrapper(viewKeyStruct::flowSolverNameString, &m_flowSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the flow solver to use in the reservoir-well solver");

  RegisterViewWrapper(viewKeyStruct::wellSolverNameString, &m_wellSolverName, 0)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the well solver to use in the reservoir-well solver");

}


void ReservoirWellSolver::ImplicitStepSetup( real64 const& time_n,
                                             real64 const& dt,
                                             DomainPartition * const domain,
                                             systemSolverInterface::EpetraBlockSystem * const blockSystem)
{
  // Todo
}

void ReservoirWellSolver::AssembleSystem( DomainPartition * const domain,
                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                          real64 const time,
                                          real64 const dt )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  flowSolver.AssembleSystem( domain, blockSystem, time, dt );
  wellSolver.AssembleSystem( domain, blockSystem, time, dt );
}
  
void ReservoirWellSolver::ApplyBoundaryConditions( DomainPartition * const domain,
                                                   systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                   real64 const time,
                                                   real64 const dt )
{
  // access the flow and well solvers
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());

  flowSolver.ApplyBoundaryConditions( domain, blockSystem, time, dt );
  // TODO: figure out if we need to apply boundary conditions for the wells
}

real64 ReservoirWellSolver::CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const *const blockSystem,
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

void ReservoirWellSolver::SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                       SystemSolverParameters const * const params )
{
  // TODO: solve the full Jacobian matrix
}

bool ReservoirWellSolver::CheckSystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                               real64 const scalingFactor,
                                               DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  bool const validReservoirSolution = flowSolver.CheckSystemSolution( blockSystem, scalingFactor, domain );
  bool const validWellSolution      = wellSolver.CheckSystemSolution( blockSystem, scalingFactor, domain );
  
  return ( validReservoirSolution && validWellSolution );
}

void ReservoirWellSolver::ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                               real64 const scalingFactor,
                                               DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  flowSolver.ApplySystemSolution( blockSystem, scalingFactor, domain );
  wellSolver.ApplySystemSolution( blockSystem, scalingFactor, domain );
}

void ReservoirWellSolver::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  FlowSolverBase & flowSolver = *(this->getParent()->GetGroup(m_flowSolverName)->group_cast<FlowSolverBase*>());
  WellSolverBase & wellSolver = *(this->getParent()->GetGroup(m_wellSolverName)->group_cast<WellSolverBase*>());

  flowSolver.ResetStateToBeginningOfStep( domain );
  wellSolver.ResetStateToBeginningOfStep( domain ); 
}
  
void ReservoirWellSolver::ImplicitStepComplete( real64 const& time_n,
                                                real64 const& dt,
                                                DomainPartition * const domain)
{
  // Todo
}

void ReservoirWellSolver::PostProcessInput()
{
  // Todo: check that the solver are compatible
}

void ReservoirWellSolver::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const problemManager)
{
  // Todo
}

ReservoirWellSolver::~ReservoirWellSolver()
{
  // TODO Auto-generated destructor stub
}

real64 ReservoirWellSolver::SolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition * domain )
{
  return 0.0;
}

REGISTER_CATALOG_ENTRY( SolverBase, ReservoirWellSolver, std::string const &, ManagedGroup * const )

} /* namespace geosx */
