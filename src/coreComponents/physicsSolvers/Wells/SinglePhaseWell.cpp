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
  m_numDofPerElement    = 1;
  m_numDofPerConnection = 1;

}

void SinglePhaseWell::RegisterDataOnMesh(ManagedGroup * const meshBodies)
{
  std::cout << "SinglePhaseWell::RegisterDataOnMesh started" << std::endl;
  
  WellSolverBase::RegisterDataOnMesh(meshBodies);

  WellManager * wellManager = meshBodies->getParent()->group_cast<DomainPartition *>()->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();
    wellElementSubRegion->RegisterViewWrapper<array1d<globalIndex>>( viewKeyStruct::dofNumberString ); 
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );

    ConnectionData * connectionData = well->getConnections(); 
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::rateString );
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaRateString );

    PerforationData * perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::perforationRateString );
    perforationData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
    
  });    

  std::cout << "SinglePhaseWell::RegisterDataOnMesh complete" << std::endl;
}
  
void SinglePhaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  std::cout << "SinglePhaseWell: InitializePreSubGroups started" << std::endl;
  
  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  WellManager * wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData * perforationData = well->getPerforations();
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString ).resizeDimension<1>(2);
  });
				   
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

void SinglePhaseWell::UpdateFluidModel( Well * const well )
{
  WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
  SingleFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

  array1d<real64 const> const & wellPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  array1d<real64 const> const & dWellPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
  {
    real64 const newWellPressure = wellPressure[iwelem] + dWellPressure[iwelem];
    fluid->PointUpdate( newWellPressure, iwelem, 0 ); 
  }
}
  
void SinglePhaseWell::UpdateStateAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    UpdateState( well );
  });    
}

void SinglePhaseWell::UpdateState( Well * const well )
{
  UpdateFluidModel( well );
}

void SinglePhaseWell::InitializeWellState( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resDensity = m_resDensity;
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    ConnectionData const * const connectionData = well->getConnections();
    PerforationData const * const perforationData = well->getPerforations();

    // get well primary variables on segments
    array1d<real64> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    // get the info stored on segments
    array1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );

    array1d<real64> const & wellElemSumRates =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::sumRatesString );

    // get a reference to the primary variables on connections
    array1d<real64> const & connRate  =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::rateString );
    
    // get the well element indices corresponding to each connect
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    arrayView1d<localIndex const> const & prevWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::prevWellElementIndexString );    
    
    // get well data on perforations
    array1d<real64> const & perfRate =
      perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );    
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // get constitutive data
    SingleFluidBase const * const fluid = GetFluidModel( wellElementSubRegion );
    
    array2d<real64> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );
    
    // TODO: fix this in the case of multiple MPI ranks working on the same well
    // In particular, steps 1), 3) and 7) will need to be rewritten
    
    // 1) Loop over all perforations to compute an average density
    real64 avgDensity = 0;
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];
      
      avgDensity += resDensity[er][esr][m_fluidIndex][ei][0];
    }
    globalIndex const numPerforationsGlobal = perforationData->numPerforationsGlobal();
    avgDensity /= numPerforationsGlobal;

    // get the reference data for this well
    localIndex const iwelemRef = 0; // here, we assume that the first segment is on this MPI rank, which will NOT be the case in general
    real64 const refGravDepth = well->getReferenceDepth();

    // 2) Initialize the reference pressure
    real64 const targetBHP = well->getTargetBHP();
    if (well->getControl() == Well::Control::BHP)
    {
      // if pressure constraint, set the ref pressure at the constraint
      wellElemPressure[iwelemRef] = targetBHP;
    }
    else // rate control
    {  
      // if rate constraint, set the ref pressure slightly above/below the target pressure
      wellElemPressure[iwelemRef] = (well->getType() == Well::Type::PRODUCER)
                                  ? 0.9 * targetBHP
		    	          : 1.1 * targetBHP;
    }
      
    // 3) Estimate the pressures in the segments using this avgDensity
    // TODO: implement this in parallel with the communication of the pressures
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (iwelem != iwelemRef)
        wellElemPressure[iwelem] = wellElemPressure[iwelemRef]
 	  + ( m_gravityFlag ? avgDensity * wellElemGravDepth[iwelem] - refGravDepth : 0 );
    }

    // 4) Recompute the pressure-dependent properties
    UpdateState( well );

    // 5) Compute the perforation rates
    ComputeAllPerforationRates( well );
    
    // 6) Collect all the perforation rates
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
      wellElemSumRates[iwelem] = 0;
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {  
      localIndex const iwelem = perfWellElemIndex[iperf];
      wellElemSumRates[iwelem] += perfRate[iperf];
    }

    // 7) Estimate the connection rates
    // TODO: implement this in parallel with the communication of wellElemSumRates
    // TODO: improve this later
    localIndex const lastConnectionIndex = connectionData->numConnectionsLocal()-1;
    for (localIndex iconn = lastConnectionIndex; iconn >= 0; --iconn)
    {
      localIndex const iwelemNext = nextWellElemIndex[iconn];
      connRate[iconn] = wellElemSumRates[iwelemNext] / wellElemDensity[iwelemNext][0];
                      
      localIndex const iwelemPrev = prevWellElemIndex[iconn];
      if (iwelemPrev >= 0)
	wellElemSumRates[iwelemPrev] += connRate[iconn];
    }
    
  });
}

void SinglePhaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  std::cout << "SinglePhaseWell::InitializePostInitialConditions_PreSubGroups started" << std::endl;
  
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>( keys::domain );

  // update properties
  UpdateStateAll( domain );
  
  std::cout << "SinglePhaseWell::InitializePostInitialConditions_PreSubGroups complete" << std::endl;
}

void SinglePhaseWell::BackupFields( DomainPartition * const domain )
{
  // will be needed if we decide to compute the accumulation term in each segment
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
  WellManager const * const wellManager = domain->getWellManager();

  localIndex localCount = 0;
  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    arrayView1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  
    // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
    for ( localIndex iwelem = 0; iwelem < wellElemGhostRank.size(); ++iwelem )
    {
      wellElemDofNumber[iwelem] = -1;
    }

    // loop over all well elements and set the dof number if the element is not a ghost
    for ( localIndex iwelem = 0; iwelem < wellElemGhostRank.size(); ++iwelem )
    {
      if (wellElemGhostRank[iwelem] < 0 )
      {
        wellElemDofNumber[iwelem] = firstLocalRow + localCount + offset;
	std::cout << "wells: currentDofNumber = " << firstLocalRow + localCount + offset << std::endl;
        localCount += 1;
      }
      else
      {
        wellElemDofNumber[iwelem] = -1;
      }
    }
  });
  	    

  GEOS_ERROR_IF(localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );

  std::cout << "SinglePhaseWell::SetNumRowsAndTrilinosIndices complete" << std::endl;
}

void SinglePhaseWell::SetSparsityPattern( DomainPartition const * const domain,
                                          Epetra_FECrsGraph * const sparsity,
					  globalIndex firstWellElemDofNumber,
					  localIndex numDofPerResElement)
{
  std::cout << "SinglePhaseWell::SetSparsityPattern started" << std::endl;
  
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();
  WellManager const * const wellManager = domain->getWellManager();

  // get the reservoir degrees of freedom
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber =
    elementRegionManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( SinglePhaseFlow::viewKeyStruct::blockLocalDofNumberString );

  // save these numbers (reused to compute the well elem offsets in multiple functions)
  m_firstWellElemDofNumber = firstWellElemDofNumber;
  m_numDofPerResElement    = numDofPerResElement;
  
  // dofs are pressure (reservoir and well) and rate (well only) 
  localIndex constexpr maxNumDof = 2; 
  localIndex const resNDOF  = numDofPerResElement; 
  localIndex const wellNDOF = numDofPerElement()
                            + numDofPerConnection();

  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    PerforationData const * const perforationData = well->getPerforations();
    ConnectionData const * const connectionData   = well->getConnections();

    // get the well degrees of freedom 
    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get the well element indices corresponding to each connect
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    arrayView1d<localIndex const> const & prevWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::prevWellElementIndexString );    
    
    // get the well element indices corresponding to each perforation
    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );    
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // 1) Insert the entries corresponding to reservoir-well perforations
    //    This will fill J_WW, J_WR, and J_RW
    //    That is all we need for single-segmented wells
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( resNDOF + wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( resNDOF + wellNDOF );

      std::cout << "er = "   << er
		<< " esr = " << esr
		<< " ei = "  << ei
		<< " resDofNumber = " << resDofNumber[er][esr][ei]
		<< std::endl;
      
      // get the offset of the reservoir element equation
      globalIndex const resOffset  = resNDOF * resDofNumber[er][esr][ei];
      // specify the reservoir equation number
      elementLocalDofIndexRow[ElemTag::RES * wellNDOF] = resOffset;
      // specify the reservoir variable number
      elementLocalDofIndexCol[ElemTag::RES * wellNDOF] = resOffset;
      std::cout << "iperf = " << iperf << " resOffset = " << resOffset << std::endl;
      
      /*
       * firstWellElemDofNumber denotes the first DOF number of the well segments, for all the wells (i.e., first segment of first well)
       * currentElemDofNumber denotes the DOF number of the current segment
       *
       * The coordinates of this element's 2x2 block in J_WW can be accessed using:
       *
       * IndexRow = firstWellElemDofNumber * resNDOF ( = all the equations in J_RR)
       *          + (currentElemDofNumber - firstWellElemDofNumber ) * wellNDOF ( = offset of current segment in J_WW)
       *          + idof ( = local dofs for this segment, pressure and rate)
       *	   
       * This is needed because resNDOF is not equal to wellNDOF
       */

      localIndex const iwelem = perfWellElemIndex[iperf]; 
      globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );
      
      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
	// specify the well equation number
	elementLocalDofIndexRow[ElemTag::WELL * resNDOF + idof] = elemOffset + idof;
	// specify the reservoir variable number
	elementLocalDofIndexCol[ElemTag::WELL * resNDOF + idof] = elemOffset + idof;
	std::cout << "iperf = " << iperf << " idof = " << idof << " wellOffset = "
		  << elementLocalDofIndexRow[ElemTag::WELL * resNDOF + idof] << std::endl;
      }      

      sparsity->InsertGlobalIndices( integer_conversion<int>( resNDOF + wellNDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( resNDOF + wellNDOF ),
                                     elementLocalDofIndexCol.data() );
    }

    // 2) Insert the entries corresponding to well-well connection between segments
    //    This will fill J_WW only
    //    That is not needed for single-segmented wells
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // TODO: figure out how to do that with multiple MPI ranks
      
      // get previous segment index
      localIndex iwelemPrev =  prevWellElemIndex[iconn];

      // get second segment index
      localIndex iwelemNext = nextWellElemIndex[iconn];

      // check if this is not an entry or exit
      if (iwelemPrev < 0 || iwelemNext < 0)
	continue;

      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( 2 * wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( 2 * wellNDOF );
      
      // get the offset of the well element equations
      globalIndex const prevElemOffset = getElementOffset( wellElemDofNumber[iwelemPrev] );
      globalIndex const nextElemOffset = getElementOffset( wellElemDofNumber[iwelemNext] );

      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
	elementLocalDofIndexRow[idof]            = prevElemOffset + idof;
	elementLocalDofIndexRow[wellNDOF + idof] = nextElemOffset + idof;
	elementLocalDofIndexCol[idof]            = prevElemOffset + idof;
	elementLocalDofIndexCol[wellNDOF + idof] = nextElemOffset + idof;
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

  FormControlEquation( domain, jacobian, residual );
  
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
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get a reference to the degree of freedom numbers
    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // get a reference to the next/prev well element index for each connection
    arrayView1d<localIndex const> const & prevWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::prevWellElementIndexString );
    
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    // get a reference to the primary variables on connections
    array1d<real64 const> const & connRate  =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::rateString );

    array1d<real64 const> const & dConnRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaRateString );
      
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      localIndex const iwelemNext = nextWellElemIndex[iconn];
      localIndex const iwelemPrev = prevWellElemIndex[iconn];
      
      if (iwelemNext < 0) // nothing to do
	continue;
      else if (iwelemPrev < 0) // exit connection
      {	
        // get well secondary variables on segments
        SingleFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

        array2d<real64> const & wellElemDensity =
          fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

        array2d<real64> const & dWellElemDensity_dPres =
          fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

	// TODO: make sure that naming is consistent (flux contains dt here?)
        real64 const currentConnRate = connRate[iconn] + dConnRate[iconn];
        real64 const flux = dt * wellElemDensity[iwelemNext][0] * currentConnRate;	
        real64 const dFlux_dRate = dt * wellElemDensity[iwelemNext][0];
	real64 const dFlux_dPres = dt * dWellElemDensity_dPres[iwelemNext][0] * currentConnRate;
	
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelemNext] );
	globalIndex const eqnRowIndex     = elemOffset + RowOffset::MASSBAL;
        globalIndex const dofColIndexPres = elemOffset + ColOffset::DPRES;
	globalIndex const dofColIndexRate = elemOffset + ColOffset::DRATE;
	
        // add contribution to global residual and jacobian
        residual->SumIntoGlobalValues( 1, &eqnRowIndex, &flux );
        jacobian->SumIntoGlobalValues( 1, &eqnRowIndex, 1, &dofColIndexPres, &dFlux_dPres );
        jacobian->SumIntoGlobalValues( 1, &eqnRowIndex, 1, &dofColIndexRate, &dFlux_dRate );
		
      }
      else // connection between two segments
      {
        // not needed for now
      }
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
  
  WellManager * const wellManager = domain->getWellManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber = m_resDofNumber;
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData * perforationData = well->getPerforations();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // compute the local rates for this well
    ComputeAllPerforationRates( well );
    
    // get the degrees of freedom 
    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get well variables on perforations
    array1d<real64> const & perfRate =
      perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

    array2d<real64> const & dPerfRate_dPres =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );

    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // loop over the perforations and add the rates to the residual and jacobian
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      // local working variables and arrays
      stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
      stackArray1d<globalIndex, 2> dofColIndices( 2 );

      stackArray1d<real64, 2> localFlux( 2 );
      stackArray2d<real64, 4> localFluxJacobian(2, 2);

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      
      // row index on reservoir side
      eqnRowIndices[ElemTag::RES] = resDofNumber[er][esr][ei];

      // column index on reservoir side
      dofColIndices[ElemTag::RES] = resDofNumber[er][esr][ei];

      // row index on well side
      eqnRowIndices[ElemTag::WELL] = wellElemDofNumber[iwelem] + RowOffset::MASSBAL;

      // column index on well side
      dofColIndices[ElemTag::WELL] = wellElemDofNumber[iwelem] + ColOffset::DPRES;

      // populate local flux vector and derivatives
      localFlux[ElemTag::RES]  =  dt * perfRate[iperf];
      localFlux[ElemTag::WELL] = -localFlux[ElemTag::RES];

      for (localIndex ke = 0; ke < 2; ++ke)
      {
        localFluxJacobian[ElemTag::RES][ke]  = dt * dPerfRate_dPres[iperf][ke];
        localFluxJacobian[ElemTag::WELL][ke] = - localFluxJacobian[ElemTag::RES][ke];
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

    }
    
  });

  std::cout << "SinglePhaseWell::AssembleSourceTerms complete" << std::endl;
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
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);
 
  WellManager * const wellManager = domain->getWellManager();

  real64 localResidualNorm = 0;
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers
    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

	localIndex const wellNDOF = numDofPerElement() + numDofPerConnection(); 
        for (localIndex idof = 0; idof < wellNDOF; ++idof)
        {
          int const lid = rowMap->LID(integer_conversion<int>( elemOffset + idof ) );
          real64 const val = localResidual[lid];
          localResidualNorm += val * val;
        }
      }
    }
  });

  // compute global residual norm
  real64 globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
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

    // get the degree of freedom numbers on segments
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // get a reference to the primary variables on segments
    array1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
    
    // get a reference to the primary variables on connections
    array1d<real64> const & dConnRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaRateString );

    // get a reference to the next well element index for each connection
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    
    
    // get the ghosting information on segments and on connections
    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract solution and apply to dP
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );
        int const lid = rowMap->LID( integer_conversion<int>( elemOffset + ColOffset::DPRES) );
        dWellElemPressure[iwelem] += scalingFactor * local_solution[lid];
      }
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // as a temporary solution, the connection variable is looked up with the next element's DOF number
      localIndex iwelemNext = nextWellElemIndex[iconn];

      // check if there is a variable defined here
      if (iwelemNext < 0)
        continue;

      if (wellElemGhostRank[iwelemNext] < 0)
      {      
        // extract solution and apply to dQ
	globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelemNext] );
        int const lid = rowMap->LID( integer_conversion<int>( elemOffset + ColOffset::DRATE ) );
        dConnRate[iconn] += scalingFactor * local_solution[lid];
      }
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
  
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the ghosting information
    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
    
    // get a reference to the primary variables on segments
    array1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get a reference to the primary variables on connections
    array1d<real64> const & dConnRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaRateString );

    // get a reference to the well element index for each connection
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    
    
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
        dWellElemPressure[iwelem] = 0;
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // as a temporary solution, the connection variable is looked up with the next element's DOF number
      localIndex iwelemNext = nextWellElemIndex[iconn];

      // check if there is a variable defined here
      if (iwelemNext < 0)
        continue;

      dConnRate[iconn] = 0;
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


void SinglePhaseWell::ComputeAllPerforationRates( Well * well )
{

  // get the reservoir data
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber = m_resDofNumber;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resPressure         = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dResPressure        = m_deltaResPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resGravDepth        = m_resGravDepth;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resDensity          = m_resDensity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResDensity_dPres   = m_dResDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & resViscosity        = m_resViscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dResViscosity_dPres = m_dResVisc_dPres;

  // get the well data
  PerforationData * perforationData = well->getPerforations();
  WellElementSubRegion * wellElementSubRegion = well->getWellElements();

  // get the degrees of freedom 
  array1d<globalIndex> const & wellElemDofNumber =
    wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

  // get well primary variables on segments
  array1d<real64> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  array1d<real64> const & dWellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  array1d<real64> const & wellElemGravDepth =
    wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );
  
  // get well secondary variables on segments
  SingleFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

  array2d<real64> const & wellElemDensity =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

  array2d<real64> const & dWellElemDensity_dPres =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

  array2d<real64> const & wellElemViscosity =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::viscosityString );

  array2d<real64> const & dWellElemViscosity_dPres =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dVisc_dPresString );

  // get well variables on perforations
  array1d<real64> const & perfGravDepth =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::gravityDepthString );

  arrayView1d<localIndex const> const & perfWellElemIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d<real64 const> const & perfTransmissibility =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::transmissibilityString );

  arrayView1d<real64> const & perfRate =
    perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

  arrayView2d<real64> const & dPerfRate_dPres =
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
  
  // get the element region, subregion, index
  arrayView1d<localIndex const> const & resElementRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

  arrayView1d<localIndex const> const & resElementSubRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

  arrayView1d<localIndex const> const & resElementIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

  for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
  {
    // local working variables and arrays
    stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
    stackArray1d<globalIndex, 2> dofColIndices( 2 );

    stackArray1d<real64, 2> pressure( 2 );
    stackArray1d<real64, 2> dPressure_dP( 2 );
    
    stackArray1d<real64, 2> dFlux_dP( 2 );

    stackArray1d<real64, 2> density( 2 );
    stackArray1d<real64, 2> dDensity_dP( 2 );

    stackArray1d<real64, 2> viscosity( 2 );
    stackArray1d<real64, 2> dViscosity_dP( 2 );

    stackArray1d<real64, 2> mobility( 2 );
    stackArray1d<real64, 2> dMobility_dP( 2 );

    stackArray1d<localIndex, 2> multiplier( 2 );

    // clear working arrays
    eqnRowIndices = -1;

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
    density[ElemTag::RES] = resDensity[er][esr][m_fluidIndex][ei][0];
    dDensity_dP[ElemTag::RES] = dResDensity_dPres[er][esr][m_fluidIndex][ei][0];

    viscosity[ElemTag::RES] = resViscosity[er][esr][m_fluidIndex][ei][0];
    dViscosity_dP[ElemTag::RES] = dResViscosity_dPres[er][esr][m_fluidIndex][ei][0];

    // TODO: ask about Ecl/Intersect treatment here
    // In AD-GPRS, depth difference between the reservoir elem center and perf is not accounted for
    pressure[ElemTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[ElemTag::RES] = 1;
    
    multiplier[ElemTag::RES] = 1;

    // 2) Well side

    localIndex const iwelem = perfWellElemIndex[iperf]; 

    // row index on well side
    eqnRowIndices[ElemTag::WELL] = wellElemDofNumber[iwelem] + RowOffset::MASSBAL;

    // column index on well side
    dofColIndices[ElemTag::WELL] = wellElemDofNumber[iwelem] + ColOffset::DPRES;

    // get well variables
    density[ElemTag::WELL] = wellElemDensity[iwelem][0];
    dDensity_dP[ElemTag::WELL] = dWellElemDensity_dPres[iwelem][0];

    viscosity[ElemTag::WELL] = wellElemViscosity[iwelem][0];
    dViscosity_dP[ElemTag::WELL] = dWellElemViscosity_dPres[iwelem][0];

    pressure[ElemTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    pressure[ElemTag::WELL] += m_gravityFlag
			     ? density[ElemTag::WELL] * ( perfGravDepth[iperf] - wellElemGravDepth[iwelem] )
			     : 0; 		
    dPressure_dP[ElemTag::WELL] = 1;
    dPressure_dP[ElemTag::WELL] += m_gravityFlag
                                 ? dDensity_dP[ElemTag::WELL] * ( perfGravDepth[iperf] - wellElemGravDepth[iwelem] )
                                 : 0;
   
    multiplier[ElemTag::WELL] = -1;

    //***** calculation of flux *****
      
    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf]; 

    real64 potDif = 0.0;
    for (localIndex i = 0; i < 2; ++i)
    {
      // compute mobility
      mobility[i]     = density[i] / viscosity[i];
      dMobility_dP[i] = dDensity_dP[i] / viscosity[i]
   	              - mobility[i] / viscosity[i] * dViscosity_dP[i];

      // compute potential difference
      potDif += multiplier[i] * trans * pressure[i];
      dPerfRate_dPres[iperf][i] = multiplier[i] * trans * dPressure_dP[i];
    }

    // compute the final flux and derivatives

    // choose upstream cell
    localIndex const k_up = (potDif >= 0) ? ElemTag::RES : ElemTag::WELL;

    perfRate[iperf] = mobility[k_up] * potDif;
    for (localIndex ke = 0; ke < 2; ++ke)
      dPerfRate_dPres[iperf][ke] *= mobility[k_up];
    dPerfRate_dPres[iperf][k_up] += dMobility_dP[k_up] * potDif;

    //***** end flux terms *****

  }
}
  
void SinglePhaseWell::FormControlEquation( DomainPartition * const domain,
                                           Epetra_FECrsMatrix * const jacobian,
                                           Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom 
    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // TODO: check that the first connection is on my rank
    
    // the well control equation is defined on the first connection
    localIndex iconnControl = 0;
    // again we assume here that the first segment is on this MPI rank
    localIndex iwelemControl = 0;

    // get well control
    Well::Control control = well->getControl();

    // BHP control
    if (control == Well::Control::BHP)
    {
      // get pressure data
      array1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      array1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // form the control equation now
      real64 const currentBHP = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];
      real64 const targetBHP  = well->getTargetBHP();
      
      real64 const controlEqn = currentBHP - targetBHP;
      real64 const dControlEqn_dPres = 1;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DPRES;
      
      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex, &controlEqn );
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex, 1, &dofColIndex, &dControlEqn_dPres );
    }
    // rate control
    else
    {

      // get a reference to the primary variables on connections
      array1d<real64 const> const & connRate  =
        connectionData->getReference<array1d<real64>>( viewKeyStruct::rateString );

      array1d<real64 const> const & dConnRate =
        connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaRateString );
      
      // for the control equation here
      real64 const currentConnRate = connRate[iconnControl] + dConnRate[iconnControl];
      real64 const targetConnRate  = well->getTargetRate();

      real64 const controlEqn = currentConnRate - targetConnRate;
      real64 const dControlEqn_dRate = 1;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DRATE;

      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex, &controlEqn );
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex, 1, &dofColIndex, &dControlEqn_dRate );
    }
  });
}

  
void SinglePhaseWell::ImplicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
  std::cout << "SinglePhaseWell::ImplicitStepComplete started" << std::endl;
  
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();
    
    // get a reference to the primary variables on segments
    array1d<real64> const & wellElemPressure  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get a reference to the next well element index
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    
    
    // get a reference to the primary variables on connections
    array1d<real64> const & connRate  =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::rateString );

    array1d<real64> const & dConnRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaRateString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // as a temporary solution, the connection variable is looked up with the next element's DOF number
      localIndex iwelemNext = nextWellElemIndex[iconn];

      // check if there is a variable defined here
      if (iwelemNext < 0)
        continue;

      connRate[iconn] += dConnRate[iconn];
    }    
  });

  std::cout << "SinglePhaseWell::ImplicitStepComplete complete" << std::endl;
}


REGISTER_CATALOG_ENTRY(SolverBase, SinglePhaseWell, string const &, ManagedGroup * const)
}// namespace geosx
