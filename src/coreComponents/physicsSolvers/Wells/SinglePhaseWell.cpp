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
#include "mesh/WellElementSubRegion.hpp"
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

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  WellManager * wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::densityString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dDensity_dPresString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::viscosityString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dViscosity_dPresString );

    ConnectionData * connectionData = well->getConnections(); 
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::velocityString );
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaVelocityString );

    PerforationData * perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::flowRateString );
    
  });    
  

  
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
  
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
                                                                BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

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

  // loop over the wells
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      WellElement const * const wellElement = wellElementSubRegion->getWellElement( iwelem );
      std::cout << "SinglePhaseWell: computing flux terms for segment "
	        << wellElement->getName()
		<< " for well " << well->getName()
	        << std::endl;
    }

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
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  WellManager * const wellManager = domain->getWellManager();

  NumericalMethodsManager * const numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager * const fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
                                                                BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

  // ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber = m_dofNumber;

  /*
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resPressure        = m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dResPressure       = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resGravDepth       = m_gravDepth;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resDensity         = m_density;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResDensity_dPres  = m_dDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resViscosity       = m_viscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResViscositydPres = m_dVisc_dPres;
  */

  // temp to make sure that it compiles
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const  resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const  dResPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const  resGravDepth;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const  resDensity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const  dResDensity_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const  resViscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const  dResViscosity_dPres;

  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData * perforationData = well->getPerforations();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get well variables on segments
    array1d<real64> const & wellPressure         = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );
    array1d<real64> const & dWellPressure        = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
    array1d<real64> const & wellDensity          = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::densityString );
    array1d<real64> const & dWellDensity_dPres   = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::dDensity_dPresString );
    array1d<real64> const & wellViscosity        = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::viscosityString );
    array1d<real64> const & dWellViscosity_dPres = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::dViscosity_dPresString );
    // get well variables on perforations
    array1d<real64> const & wellGravDepth = perforationData->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );
    array1d<real64> const & wellFlowRate  = perforationData->getReference<array1d<real64>>( viewKeyStruct::flowRateString );

    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    real64 const densWeight[2] = { 1.0, 0.0 }; // cell / well weights

    // temp hack
    real64 const gravD = 0;
    
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      // local working variables and arrays
      real64 densMean = 0.0;
      stackArray1d<real64, 2> dDensMean_dP( 2 );
      
      stackArray1d<real64, 2> pressure( 2 );
      stackArray1d<real64, 2> gravDepth( 2 );
      stackArray1d<real64, 2> dFlux_dP( 2 );
      
      stackArray1d<real64, 2> density( 2 );
      stackArray1d<real64, 2> dDensity_dP( 2 );

      stackArray1d<real64, 2> viscosity( 2 );
      stackArray1d<real64, 2> dViscosity_dP( 2 );

      stackArray1d<real64, 2> mobility( 2 );
      stackArray1d<real64, 2> dMobility_dP( 2 );

      stackArray1d<localIndex, 2> weightIndex( 2 );
	
      /*
      globalIndex eqnRowIndex = -1;
      stackArray1d<globalIndex, numElems> dofColIndices( numElems );
      dofColIndices = -1;
      */

      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // TODO: put this in a loop
      // TODO: add a tag for Res = 0 and Well = 1
      
      // get reservoir variables
      pressure[0]  = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
      gravDepth[0] = resGravDepth[er][esr][ei]; 
      
      density[0]     = resDensity[er][esr][m_fluidIndex][ei][0];
      dDensity_dP[0] = dResDensity_dPres[er][esr][m_fluidIndex][ei][0];

      viscosity[0]     = resViscosity[er][esr][m_fluidIndex][ei][0];
      dViscosity_dP[0] = dResViscosity_dPres[er][esr][m_fluidIndex][ei][0];

      mobility[0]      = density[0] / viscosity[0];
      dMobility_dP[0]  = dDensity_dP[0] / viscosity[0] - mobility[0] / viscosity[0] * dViscosity_dP[0];
            
      weightIndex[0] = 0;

      densMean += densWeight[weightIndex[0]] * density[0];
      dDensMean_dP[0] = densWeight[weightIndex[0]] * dDensity_dP[0];

      // get well variables
      localIndex const iwelem = 0; // this is a hack
      
      pressure[1]  = wellPressure[iwelem] + dWellPressure[iwelem];
      gravDepth[1] = wellGravDepth[iwelem];
      
      density[1]     = wellDensity[iwelem];
      dDensity_dP[1] = dWellDensity_dPres[iwelem];

      viscosity[1]     = wellViscosity[iwelem];
      dViscosity_dP[1] = dWellViscosity_dPres[iwelem];

      mobility[1]     = density[0] / viscosity[0];
      dMobility_dP[1] = dDensity_dP[0] / viscosity[0] - mobility[0] / viscosity[0] * dViscosity_dP[0];

      weightIndex[1] = 1;

      densMean += densWeight[weightIndex[1]] * density[1];
      dDensMean_dP[1] = densWeight[weightIndex[1]] * dDensity_dP[1];
      
      //***** calculation of flux *****

      // get transmissibility at the interface
      Perforation * perforation = perforationData->getPerforation( iperf ); 
      const real64 trans = perforation->getTransmissibility(); // changed to an array of transmissibilities
      
      // compute potential difference MPFA-style
      real64 potDif = 0.0;
      for (localIndex i = 0; i < 2; ++i)
      {
      
        real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
        real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

	potDif += trans * (pressure[i] + gravTerm);
        dFlux_dP[i] = trans * (1.0 + dGrav_dP);
      }

      // compute the final flux and derivatives
      real64 const flux = mobility[0] * potDif;
      for (localIndex ke = 0; ke < 2; ++ke)
        dFlux_dP[ke] *= mobility[0];
      dFlux_dP[0] += dMobility_dP[0] * potDif;

      //***** end flux terms *****

      // populate local flux vector and derivatives
      // TODO

      // save flux from/into the perforation for calculating well totals
      // wellFlowRate[iperf] = localFlux * sign;

      // Add to global residual/jacobian
      //
    }
  });
    
  // loop over the wells
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    
    // for debugging purposes
    
    PerforationData const * const perforationData = well->getPerforations();

    // for each well, loop over the connections
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {      
      Perforation const * const perforation = perforationData->getPerforation( iperf );

      std::cout << "SinglePhaseWell: computing source terms for perforation "
		<< perforation->getName()
		<< " for well " << well->getName()
		<< std::endl;      
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
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  WellManager * const wellManager = domain->getWellManager();

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  // get the update
  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get a reference to the primary variables on segments
    array1d<real64> const & dWellPressure = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
    // get a reference to the primary variables on connections
    array1d<real64> const & dWellVelocity = connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaVelocityString );
   
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      // TODO: check for ghost segments

      // extract solution and apply to dP
      globalIndex const dummyDofNumber = 0;
      int const lid = rowMap->LID( integer_conversion<int>( dummyDofNumber ) );
      dWellPressure[iwelem] += scalingFactor * local_solution[lid];
      
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {

      // TODO: check for ghost connections if needed
      // TODO: check if there is a primary var defined on this connection

      // extract solution and apply to dP
      globalIndex const dummyDofNumber = 0;
      int const lid = rowMap->LID( integer_conversion<int>( dummyDofNumber ) );
      dWellVelocity[iconn] += scalingFactor * local_solution[lid];
      
    }
  });  

  // TODO: call CommunicationTools::SynchronizeFields

  // update properties
  UpdateStateAll( domain );
  
}

void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get a reference to the primary variables on segments
    array1d<real64> const & dWellPressure = wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
    // get a reference to the primary variables on connections
    array1d<real64> const & dWellVelocity = connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaVelocityString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      // TODO: check for ghost segments

      // extract solution and apply to dP
      dWellPressure[iwelem] = 0;
      
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {

      // TODO: check for ghost connections if needed
      // TODO: check if there is a primary var defined on this connection

      // extract solution and apply to dP
      dWellVelocity[iconn] = 0;
      
    }
  });

  // call constitutive models
  UpdateStateAll( domain );

}

void SinglePhaseWell::ImplicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
  
}


REGISTER_CATALOG_ENTRY(SolverBase, SinglePhaseWell, string const &, ManagedGroup * const)
}// namespace geosx
