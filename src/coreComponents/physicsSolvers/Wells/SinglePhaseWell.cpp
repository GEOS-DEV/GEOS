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
#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "wells/WellManager.hpp"
#include "wells/Well.hpp"
#include "wells/PerforationData.hpp"
#include "wells/Perforation.hpp"
#include "wells/ConnectionData.hpp"
#include "wells/Connection.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

SinglePhaseWell::SinglePhaseWell( const string & name,
                                  ManagedGroup * const parent )
  :
  WellSolverBase( name, parent )
{
  m_numDofPerWellElement = 1;
  m_numDofPerConnection  = 1;
}

void SinglePhaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );

  std::cout << "SinglePhaseWell: InitializePreSubGroups" << std::endl;  
}

SingleFluidBase * SinglePhaseWell::GetFluidModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  SingleFluidBase * const fluid = constitutiveModels->GetGroup<SingleFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

SingleFluidBase const * SinglePhaseWell::GetFluidModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  SingleFluidBase const * const fluid = constitutiveModels->GetGroup<SingleFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

void SinglePhaseWell::UpdateFluidModel( ManagedGroup * const dataGroup )
{
  
}

void SinglePhaseWell::UpdateFluidModelAll( DomainPartition * const domain )
{

}
  
void SinglePhaseWell::UpdateState( ManagedGroup * const dataGroup )
{
  UpdateFluidModel( dataGroup );
}

void SinglePhaseWell::UpdateStateAll( DomainPartition * const domain )
{
  UpdateFluidModelAll( domain );
}

void SinglePhaseWell::InitializeWellState( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
     // do something
  });
}

void SinglePhaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  std::cout << "SinglePhaseWell: InitializePostInitialConditions_PreSubGroups" << std::endl;
}

real64 SinglePhaseWell::SolverStep( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer const cycleNumber,
                                                      DomainPartition * const domain )
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

void SinglePhaseWell::BackupFields( DomainPartition * const domain )
{
  
}

void
SinglePhaseWell::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                    DomainPartition * const domain,
                                    EpetraBlockSystem * const blockSystem )
{

}

void SinglePhaseWell::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                    localIndex & numLocalRows,
                                                    globalIndex & numGlobalRows,
                                                    localIndex offset )
{
  
}

void SinglePhaseWell::SetSparsityPattern( DomainPartition const * const domain,
                                          Epetra_FECrsGraph * const sparsity )
{
  
}

void SinglePhaseWell::SetupSystem( DomainPartition * const domain,
                                   EpetraBlockSystem * const blockSystem )
{
  
}

void SinglePhaseWell::AssembleSystem( DomainPartition * const domain,
                                      EpetraBlockSystem * const blockSystem,
                                      real64 const time_n, real64 const dt )
{
  std::cout << "SinglePhaseWell: AssembleSystem" << std::endl;
  
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  CheckWellControlSwitch( domain );
  
  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleSourceTerms( domain, blockSystem, time_n, dt );
  
  jacobian->GlobalAssemble( true );
  residual->GlobalAssemble();

  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After SinglePhaseWell::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

}

void SinglePhaseWell::AssembleAccumulationTerms( DomainPartition * const domain,
                                                 Epetra_FECrsMatrix * const jacobian,
                                                 Epetra_FEVector * const residual,
                                                 real64 const time_n,
                                                 real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    // loop over the segments
  });  
}

void SinglePhaseWell::AssembleFluxTerms( DomainPartition * const domain,
                                         Epetra_FECrsMatrix * const jacobian,
                                         Epetra_FEVector * const residual,
                                         real64 const time_n,
                                         real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  // loop over the wells
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData const * const connectionData = well->getConnections();

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      Connection const * const connection = connectionData->getConnection( iconn );
      std::cout << "SinglePhaseWell: computing flux terms for connection "
	        << connection->getName()
		<< " for well " << well->getName()
	        << std::endl;
    }

  }); 

}


void SinglePhaseWell::AssembleSourceTerms( DomainPartition * const domain,
                                           EpetraBlockSystem * const blockSystem,
                                           real64 const time_n,
                                           real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  // loop over the wells
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData const * const perforationData = well->getPerforations();

    // for each well, loop over the connections
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {      
      Perforation const * const perforation = perforationData->getPerforation( iperf );

      std::cout << "SinglePhaseWell: computing source terms for perforation "
		<< perforation->getName()
		<< " for well " << well->getName()
		<< std::endl;

      //  -- Compute the rates at each perforation
      //  -- Form the control equations
      //  -- Add to residual and Jacobian matrix
      
    }
    
  });    
}

void SinglePhaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    // check if the well control needs to be switched
  });
}


real64
SinglePhaseWell::CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                                        DomainPartition * const domain )
{

  return 0.0;
}


void SinglePhaseWell::SolveSystem( EpetraBlockSystem * const blockSystem,
                                   SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
    solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  Epetra_FEVector * const
    residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  residual->Scale(-1.0);
  solution->Scale(0.0);

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("\nSolution:\n" << *solution);
  }
}

bool
SinglePhaseWell::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  return false;
}

void
SinglePhaseWell::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  
}

void
SinglePhaseWell::ApplyBoundaryConditions( DomainPartition * const domain,
                                          systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                          real64 const time_n,
                                          real64 const dt )
{

}

void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  std::cout << "SinglePhaseWell: ResetStateToBeginningOfStep" << std::endl;
  
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    //   -- set dWellPres = 0;
    //   -- update the other well variables
  });
}

void SinglePhaseWell::ImplicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
  
}


REGISTER_CATALOG_ENTRY(SolverBase, SinglePhaseWell, string const &, ManagedGroup * const)
}// namespace geosx
