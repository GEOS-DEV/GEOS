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
 * @file SinglePhaseWell.cpp
 */

#include "SinglePhaseWell.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "wells/WellManager.hpp"
#include "wells/SimpleWell.hpp"
#include "wells/SinglePhaseWell.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

SinglePhaseWell::SinglePhaseWell( const std::string& name,
                                        ManagedGroup * const parent ):
  WellSolverBase(name, parent)
{

}

void SinglePhaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );
}

void SinglePhaseWell::UpdateConstitutiveModels(DomainPartition * const domain)
{

}

void SinglePhaseWell::InitializeWellState( DomainPartition * const domain )
{

}


void SinglePhaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

}


real64 SinglePhaseWell::SolverStep( real64 const& time_n,
                                          real64 const& dt,
                                          const int cycleNumber,
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

void SinglePhaseWell::ImplicitStepSetup( real64 const& time_n,
                                               real64 const& dt,
                                               DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem )
{

}


void SinglePhaseWell::ImplicitStepComplete( real64 const & time_n,
                                                  real64 const & dt,
                                                  DomainPartition * const domain )
{

}

void SinglePhaseWell::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                          localIndex & numLocalRows,
                                                          globalIndex & numGlobalRows,
                                                          localIndex offset )
{
  
}

void SinglePhaseWell::SetupSystem ( real64 const & time_n,
                                          real64 const & dt,
                                          DomainPartition * const domain,
                                          EpetraBlockSystem * const blockSystem )
{

}

void SinglePhaseWell::SetSparsityPattern( real64 const & time_n,
                                                real64 const & dt,
                                                DomainPartition const * const domain,
                                                Epetra_FECrsGraph * const sparsity)
{

}

void SinglePhaseWell::AssembleSystem( DomainPartition * const domain,
                                            EpetraBlockSystem * const blockSystem,
                                            real64 const time_n,
                                            real64 const dt )
{
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock, BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  CheckWellControlSwitch( domain );
  
  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleSourceTerms( domain, blockSystem, time_n, dt );
  
  jacobian->GlobalAssemble(true);
  residual->GlobalAssemble();

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After SinglePhaseWell::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }
}

void
SinglePhaseWell::AssembleAccumulationTerms( DomainPartition const * const domain,
                                                  Epetra_FECrsMatrix * const jacobian,
                                                  Epetra_FEVector * const residual,
                                                  real64 const time_n,
                                                  real64 const dt )
{

}


void SinglePhaseWell::AssembleFluxTerms( DomainPartition const * const domain,
                                               Epetra_FECrsMatrix * const jacobian,
                                               Epetra_FEVector * const residual,
                                               real64 const time_n,
                                               real64 const dt )
{

}


void SinglePhaseWell::AssembleSourceTerms( DomainPartition * const domain,
                                           systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                           real64 const time_n,
                                           real64 const dt )
{

}


void SinglePhaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  WellManager * const wellManager = mesh->getWellManager();

  //  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  //  {
    // check if the well control needs to be switched
  //  });
}

real64
SinglePhaseWell::
CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                       DomainPartition * const domain )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );

  return 0.0;
}

void SinglePhaseWell::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                 real64 const scalingFactor,
                                                 DomainPartition * const domain )
{

}


void SinglePhaseWell::SolveSystem( EpetraBlockSystem * const blockSystem,
                                         SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  residual->Scale(-1.0);

  solution->Scale(0.0);


  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After SinglePhaseWell::SolveSystem");
    GEOS_LOG_RANK("\nsolution\n" << *solution);
  }
}

void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{

}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseWell, std::string const &, ManagedGroup * const )
} /* namespace geosx */
