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
#include "physicsSolvers/FiniteVolume/SinglePhaseFlow.hpp"
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

void SinglePhaseWell::RegisterDataOnMesh(ManagedGroup * const meshBodies)
{
  std::cout << "SinglePhaseWell::RegisterDataOnMesh started" << std::endl;
  
  WellSolverBase::RegisterDataOnMesh(meshBodies);

  WellManager * wellManager = meshBodies->GetGroup<MeshBody>(0)->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();
    wellElementSubRegion->RegisterViewWrapper<array1d<localIndex>>( viewKeyStruct::dofNumberString ); 
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );

    ConnectionData * connectionData = well->getConnections(); 
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::velocityString );
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaVelocityString );

    PerforationData * perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::flowRateString );
    
  });    

  std::cout << "SinglePhaseWell::RegisterDataOnMesh complete" << std::endl;
}
  
void SinglePhaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  std::cout << "SinglePhaseWell: InitializePreSubGroups started" << std::endl;
  
  WellSolverBase::InitializePreSubGroups( rootGroup );

  std::cout << "SinglePhaseWell: InitializePreSubGroups complete" << std::endl;  
}

SingleFluidBase * SinglePhaseWell::GetFluidModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels =
    dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  SingleFluidBase * const fluid = constitutiveModels->GetGroup<SingleFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

SingleFluidBase const * SinglePhaseWell::GetFluidModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels =
    dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  SingleFluidBase const * const fluid = constitutiveModels->GetGroup<SingleFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

void SinglePhaseWell::UpdateFluidModelAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
    SingleFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

    array1d<real64> const & wellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      real64 const newWellPressure = wellPressure[iwelem] + dWellPressure[iwelem];
      fluid->PointUpdate( newWellPressure, iwelem, 0 ); 
    }

  });    
}
  
void SinglePhaseWell::UpdateStateAll( DomainPartition * const domain )
{
  UpdateFluidModelAll( domain );
}

void SinglePhaseWell::InitializeWellState( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
     // do something
  });
}

void SinglePhaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  std::cout << "SinglePhaseWell::InitializePostInitialConditions_PreSubGroups started" << std::endl;
  
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  std::cout << "SinglePhaseWell::InitializePostInitialConditions_PreSubGroups complete" << std::endl;
}

void SinglePhaseWell::BackupFields( DomainPartition * const domain )
{
}

void
SinglePhaseWell::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                    DomainPartition * const domain,
                                    EpetraBlockSystem * const blockSystem )
{
  std::cout << "SinglePhaseWell::ImplicitStepSetup started" << std::endl;
  
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain );

  // the setup of dof numbers and linear system
  // is done in ReservoirWellSolver

  std::cout << "SinglePhaseWell::ImplicitStepSetup" << std::endl;
}

void SinglePhaseWell::SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
						    localIndex  & numLocalRows,
                                                    globalIndex & numGlobalRows,
                                                    localIndex offset )
{
  std::cout << "SinglePhaseWell::SetNumRowsAndTrilinosIndices started" << std::endl;
  
  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather(numMpiProcesses);

  // communicate the number of local rows to each process
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                &gather.front(),
                                                1 );

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<thisMpiProcess)
      firstLocalRow += gather[p];
  }

  // TODO: double check this for multiple MPI processes
  
  // get the well information
  WellManager const * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

  localIndex localCount = 0;
  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    
    arrayView1d<globalIndex> const & wellDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  
    // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
    for ( localIndex iwelem = 0; iwelem < wellElemGhostRank.size(); ++iwelem )
    {
      wellDofNumber[iwelem] = -1;
    }

    // loop over all well elements and set the dof number if the element is not a ghost
    for ( localIndex iwelem = 0; iwelem < wellElemGhostRank.size(); ++iwelem )
    {
      if (wellElemGhostRank[iwelem] < 0 )
      {
        wellDofNumber[iwelem] = firstLocalRow + localCount + offset;
        localCount += 1;
      }
      else
      {
        wellDofNumber[iwelem] = -1;
      }
    }
  });
  	    

  GEOS_ERROR_IF(localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );

  std::cout << "SinglePhaseWell::SetNumRowsAndTrilinosIndices complete" << std::endl;
}

void SinglePhaseWell::SetSparsityPattern( DomainPartition const * const domain,
                                          Epetra_FECrsGraph * const sparsity )
{
  std::cout << "SinglePhaseWell::SetSparsityPattern started" << std::endl;
  
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();
  WellManager const * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber =
    elementRegionManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( SinglePhaseFlow::viewKeyStruct::blockLocalDofNumberString );

  // TODO: use that
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> const & resElemGhostRank =
    elementRegionManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  localIndex constexpr maxNumDof = 2;      // dofs are pressure (reservoir and well) and velocity (well only)
  localIndex const resNDOF = 1;            // dof is pressure
  localIndex const wellNDOF = resNDOF + 1; // dofs are pressure and velocity

  //**** Loop over all perforations (reservoir-well entries in J_RW and J_WR)
  // Fill in sparsity for all pairs of DOF/elem that are connected by a perforation
  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    PerforationData const * const perforationData = well->getPerforations();
    ConnectionData const * const connectionData   = well->getConnections();

    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( resNDOF + wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( resNDOF + wellNDOF );

      globalIndex const resOffset  = resNDOF * resDofNumber[er][esr][ei];
      elementLocalDofIndexRow[ElemTag::RES * wellNDOF] = resOffset;
      elementLocalDofIndexCol[ElemTag::RES * wellNDOF] = resOffset;
      
      // TODO: globalIndex const wellOffset = firstWellOffset + wellNDOF * iwell * wellElemLocalDofNumber[iwelem]
      globalIndex const wellOffset = 0;
      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
	elementLocalDofIndexRow[ElemTag::WELL * resNDOF + idof] = wellOffset + idof;
	elementLocalDofIndexCol[ElemTag::WELL * resNDOF + idof] = wellOffset + idof;
      }      

      sparsity->InsertGlobalIndices( integer_conversion<int>( resNDOF + wellNDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( resNDOF + wellNDOF ),
                                     elementLocalDofIndexCol.data() );
    }

    //**** Loop over all connections between segments in the wellbore (off-diagonal entries in J_WW)
    //     Fill in sparsity for all pairs of DOF/elem that are connected by a connection
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // this is not needed for now
      
      // TODO: check if the connection has a primary var
      
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( 2 * wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( 2 * wellNDOF );

      // TODO: get previous segment

      // TODO: get second segment

      // TODO: globalIndex const wellOffsetPrev = firstWellOffset + wellNDOF * iwell * wellElemLocalDofNumber[iwelemPrev]
      globalIndex const wellOffsetPrev = 0;
      // TODO: globalIndex const wellOffsetNext = firstWellOffset + wellNDOF * iwell * wellElemLocalDofNumber[iwelemNext]
      globalIndex const wellOffsetNext = 0;
      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
	elementLocalDofIndexRow[idof]            = wellOffsetPrev + idof;
	elementLocalDofIndexRow[wellNDOF + idof] = wellOffsetNext + idof;
	elementLocalDofIndexCol[idof]            = wellOffsetPrev + idof;
	elementLocalDofIndexCol[wellNDOF + idof] = wellOffsetNext + idof;
      }      

      sparsity->InsertGlobalIndices( integer_conversion<int>( 2 * wellNDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( 2 * wellNDOF ),
                                     elementLocalDofIndexCol.data() );
    }
    
  });
  
  std::cout << "SinglePhaseWell::SetSparsityPattern complete" << std::endl;
}

  
void SinglePhaseWell::AssembleSystem( DomainPartition * const domain,
                                      EpetraBlockSystem * const blockSystem,
                                      real64 const time_n, real64 const dt )
{
  std::cout << "SinglePhaseWell: AssembleSystem started" << std::endl;

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
                                                                BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

  CheckWellControlSwitch( domain );

  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleSourceTerms( domain, jacobian, residual, time_n, dt );

  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After SinglePhaseWell::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

  std::cout << "SinglePhaseWell: AssembleSystem complete" << std::endl;
}

void SinglePhaseWell::AssembleAccumulationTerms( DomainPartition * const domain,
                                                 Epetra_FECrsMatrix * const jacobian,
                                                 Epetra_FEVector * const residual,
                                                 real64 const time_n,
                                                 real64 const dt )
{
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

  // NOT NEEDED FOR NOW

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
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

  // NOT NEEDED FOR NOW

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
                                           Epetra_FECrsMatrix * const jacobian,
                                           Epetra_FEVector * const residual,
                                           real64 const time_n,
                                           real64 const dt )
{
  std::cout << "SinglePhaseWell::AssembleSourceTerms started" << std::endl;
  
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();


  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber = m_resDofNumber;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resPressure         = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dResPressure        = m_deltaResPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resGravDepth        = m_resGravDepth;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resDensity          = m_resDensity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResDensity_dPres   = m_dResDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resViscosity        = m_resViscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResViscosity_dPres = m_dResVisc_dPres;

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData * perforationData = well->getPerforations();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom 
    array1d<globalIndex> const & wellDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get well primary variables on segments
    array1d<real64> const & wellPressure         =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64> const & dWellPressure        =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get well secondary variables on segments
    SingleFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

    array2d<real64> const & wellDensity          =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

    array2d<real64> const & dWellDensity_dPres   =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

    array2d<real64> const & wellViscosity        =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::viscosityString );

    array2d<real64> const & dWellViscosity_dPres =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dVisc_dPresString );

    // get well variables on perforations
    array1d<real64> const & wellFlowRate  =
      perforationData->getReference<array1d<real64>>( viewKeyStruct::flowRateString );

    array1d<real64> const & wellGravDepth =
      perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::gravityDepthString );

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
      stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
      stackArray1d<globalIndex, 2> dofColIndices( 2 );

      stackArray1d<real64, 2> localFlux( 2 );
      stackArray2d<real64, 4> localFluxJacobian(2, 2);

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
      stackArray1d<localIndex, 2> multiplier( 2 );

      // clear working arrays
      eqnRowIndices = -1;
      dDensMean_dP = 0.0;

      // TODO: add localFlux / localFluxJacobian

      // 1) Reservoir side

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // row index on reservoir side
      eqnRowIndices[ElemTag::RES] = resDofNumber[er][esr][ei];

      // column index on reservoir side
      dofColIndices[ElemTag::RES] = resDofNumber[er][esr][ei];

      // get reservoir variables
      pressure[ElemTag::RES]  = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
      gravDepth[ElemTag::RES] = resGravDepth[er][esr][ei]; 

      density[ElemTag::RES]     = resDensity[er][esr][m_fluidIndex][ei][0];
      dDensity_dP[ElemTag::RES] = dResDensity_dPres[er][esr][m_fluidIndex][ei][0];

      viscosity[ElemTag::RES]     = resViscosity[er][esr][m_fluidIndex][ei][0];
      dViscosity_dP[ElemTag::RES] = dResViscosity_dPres[er][esr][m_fluidIndex][ei][0];

      // TODO: check whether we need this weight
      weightIndex[ElemTag::RES] = 0;
      multiplier[ElemTag::RES]  = 1;

      densMean += densWeight[weightIndex[ElemTag::RES]] * density[ElemTag::RES];
      dDensMean_dP[ElemTag::RES] = densWeight[weightIndex[ElemTag::RES]] * dDensity_dP[ElemTag::RES];

      // 2) Well side

      localIndex const iwelem = 0; // this is a hack

      // row index on well side
      eqnRowIndices[ElemTag::WELL] = wellDofNumber[iwelem];

      // column index on well side
      dofColIndices[ElemTag::WELL] = wellDofNumber[iwelem];

      // get well variables
      pressure[ElemTag::WELL]  = wellPressure[iwelem] + dWellPressure[iwelem];
      gravDepth[ElemTag::WELL] = wellGravDepth[iwelem];

      density[ElemTag::WELL]     = wellDensity[iwelem][0];
      dDensity_dP[ElemTag::WELL] = dWellDensity_dPres[iwelem][0];

      viscosity[ElemTag::WELL]     = wellViscosity[iwelem][0];
      dViscosity_dP[ElemTag::WELL] = dWellViscosity_dPres[iwelem][0];

      weightIndex[ElemTag::WELL] = 1;
      multiplier[ElemTag::WELL]  = -1;

      densMean += densWeight[weightIndex[ElemTag::WELL]] * density[ElemTag::WELL];
      dDensMean_dP[ElemTag::WELL] = densWeight[weightIndex[ElemTag::WELL]] * dDensity_dP[ElemTag::WELL];

      //***** calculation of flux *****

      // TODO: use distinct treatments for injector and producer

      // get transmissibility at the interface
      Perforation * perforation = perforationData->getPerforation( iperf ); 
      const real64 trans = perforation->getTransmissibility(); // TODO: change to an array of transmissibilities in PerforationData

      real64 potDif = 0.0;
      for (localIndex i = 0; i < 2; ++i)
      {
        // compute mobility
        mobility[i]     = density[i] / viscosity[i];
        dMobility_dP[i] = dDensity_dP[i] / viscosity[i]
   	                - mobility[i] / viscosity[i] * dViscosity_dP[i];

	// compute gravity terms
	real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
        real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

	// compute potential difference
	potDif += multiplier[i] * trans * (pressure[i] + gravTerm);
        dFlux_dP[i] = multiplier[i] * trans * (1.0 + dGrav_dP);
      }

      // compute the final flux and derivatives

      // choose upstream cell
      localIndex const k_up = (potDif >= 0) ? ElemTag::RES : ElemTag::WELL;

      real64 const flux = mobility[k_up] * potDif;
      for (localIndex ke = 0; ke < 2; ++ke)
        dFlux_dP[ke] *= mobility[k_up];
      dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

      //***** end flux terms *****

      // populate local flux vector and derivatives
      localFlux[ElemTag::RES]  =  dt * flux;
      localFlux[ElemTag::WELL] = -localFlux[ElemTag::RES];

      for (localIndex ke = 0; ke < 2; ++ke)
      {
        localFluxJacobian[ElemTag::RES][ke]  =   dt * dFlux_dP[ke];
        localFluxJacobian[ElemTag::WELL][ke] = - dt * dFlux_dP[ke];
      }

      // Add to global residual/jacobian
      jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 ),
                                     eqnRowIndices.data(),
                                     integer_conversion<int>( 2 ),
                                     dofColIndices.data(),
                                     localFluxJacobian.data() );

      residual->SumIntoGlobalValues( integer_conversion<int>( 2 ),
				     eqnRowIndices.data(),
				     localFlux.data() );

      // TODO:
      // save flux from/into the perforation for constraint and calculating well totals
    }
  });

  std::cout << "SinglePhaseWell::AssembleSourceTerms complete" << std::endl;
}

void SinglePhaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

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
  std::cout << "SinglePhaseWell::ApplySystemSolution started" << std::endl;
  
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

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
    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get a reference to the primary variables on connections
    array1d<real64> const & dWellVelocity =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaVelocityString );

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

  std::cout << "SinglePhaseWell::ApplySystemSolution complete" << std::endl;
}

void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  std::cout << "SinglePhaseWell::ResetStateToBeginningOfStep started" << std::endl;
  
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get a reference to the primary variables on segments
    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get a reference to the primary variables on connections
    array1d<real64> const & dWellVelocity =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaVelocityString );

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

  std::cout << "SinglePhaseWell::ResetStateToBeginningOfStep complete" << std::endl;
}

void SinglePhaseWell::ResetViews(DomainPartition * const domain)
{
  std::cout << "SinglePhaseWell::ResetViews started" << std::endl;
  
  WellSolverBase::ResetViews(domain);

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_resDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( SinglePhaseFlow::viewKeyStruct::blockLocalDofNumberString );
  
  m_resPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( SinglePhaseFlow::viewKeyStruct::pressureString );

  m_deltaResPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( SinglePhaseFlow::viewKeyStruct::deltaPressureString );

  m_resDensity =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::densityString,
                                                                                      constitutiveManager );
  m_dResDens_dPres =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString,
                                                                                      constitutiveManager );
  m_resViscosity =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::viscosityString,
                                                                                      constitutiveManager );
  m_dResVisc_dPres =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::dVisc_dPresString,
                                                                                      constitutiveManager );

  std::cout << "SinglePhaseWell::ResetViews complete" << std::endl;
}
  
void SinglePhaseWell::ImplicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
  std::cout << "SinglePhaseWell::ImplicitStepComplete started" << std::endl;
  
  WellManager * const wellManager = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();
    
    // get a reference to the primary variables on segments
    array1d<real64> const & wellPressure  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get a reference to the primary variables on connections
    array1d<real64> const & wellVelocity  =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::velocityString );

    array1d<real64> const & dWellVelocity =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaVelocityString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      wellPressure[iwelem] += dWellPressure[iwelem];
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // TODO: check if there is a variable on this connection
      wellVelocity[iconn] += dWellVelocity[iconn];
    }    
  });

  std::cout << "SinglePhaseWell::ImplicitStepComplete complete" << std::endl;
}


REGISTER_CATALOG_ENTRY(SolverBase, SinglePhaseWell, string const &, ManagedGroup * const)
}// namespace geosx
