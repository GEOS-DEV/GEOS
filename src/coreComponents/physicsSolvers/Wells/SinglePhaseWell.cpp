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
  WellSolverBase::RegisterDataOnMesh(meshBodies);

  WellManager * const wellManager = meshBodies->getParent()->group_cast<DomainPartition *>()->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
    wellElementSubRegion->RegisterViewWrapper<array1d<globalIndex>>( viewKeyStruct::dofNumberString ); 
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::massBalanceNormalizerString );
    
    ConnectionData * const connectionData = well->getConnections(); 
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::rateString );
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaRateString );

    PerforationData * const perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::perforationRateString );
    perforationData->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
    
  });    
}
  
void SinglePhaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData * const perforationData = well->getPerforations();
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString ).resizeDimension<1>(2);
  });
}

SingleFluidBase * SinglePhaseWell::GetFluidModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels =
    dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  SingleFluidBase * const fluid = constitutiveModels->GetGroup<SingleFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" )
;
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

void SinglePhaseWell::InitializeWells( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & resPressure = m_resPressure;
  
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

    array1d<real64> const & wellElemSumPerfRates =
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

    // TODO: communicate avgDens
    // TODO: collect individual dens and then communicate
    
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
      real64 resPres = well->getTargetBHP();
      if (perforationData->numPerforationsLocal() > 0)
      {
        // TODO: again, this will need to be improved
        // we have to make sure that iperf = 0 is the first perforation
        localIndex const er  = resElementRegion[0];
        localIndex const esr = resElementSubRegion[0];
        localIndex const ei  = resElementIndex[0];
        resPres = resPressure[er][esr][ei];
      }
      // if rate constraint, set the ref pressure slightly above/below the target pressure
      wellElemPressure[iwelemRef] = (well->getType() == Well::Type::PRODUCER)
	                          ? 0.9 * resPres
		    	          : 1.1 * resPres;
    }

    // TODO: communicate ref pressure
    // TODO: communicate avgDens
    
    // 3) Estimate the pressures in the segments using this avgDensity
    // TODO: implement this in parallel with the communication of the pressures
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (iwelem != iwelemRef)
        wellElemPressure[iwelem] = wellElemPressure[iwelemRef]
 	  + ( m_gravityFlag ? avgDensity * ( wellElemGravDepth[iwelem] - refGravDepth ) : 0 );
    }
    
    // 4) Recompute the pressure-dependent properties
    UpdateState( well );

    // 5) Compute the perforation rates
    ComputeAllPerforationRates( well );
    
    // 6) Collect all the perforation rates
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
      wellElemSumPerfRates[iwelem] = 0.0; // TODO: ask for a better way to do that!
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {  
      localIndex const iwelem = perfWellElemIndex[iperf];
      wellElemSumPerfRates[iwelem] += perfRate[iperf];
    }

    // 7) Estimate the connection rates
    // TODO: implement this in parallel with the communication of wellElemSumPerfRates
    real64 prevConnRate = 0;
    localIndex const lastConnectionIndex = connectionData->numConnectionsLocal()-1;
    for (localIndex iconn = lastConnectionIndex; iconn >= 0; --iconn)
    {
      localIndex const iwelemNext = nextWellElemIndex[iconn];
      connRate[iconn] = prevConnRate
	              - wellElemSumPerfRates[iwelemNext] / wellElemDensity[iwelemNext][0];
      prevConnRate = connRate[iconn];
    }
  });
}

void SinglePhaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );

  // Initialize the primary and secondary variables
  InitializeWells( domain );
}

void SinglePhaseWell::BackupFields( DomainPartition * const domain,
				    real64 const & dt )
{
  // get the well information
  WellManager const * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    ConnectionData const * const connectionData = well->getConnections();
    
    // get the normalizer for the mass balance equations on well elements
    array1d<real64> const & massBalanceNormalizer =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::massBalanceNormalizerString );

    // get a reference to the primary variables on connections
    array1d<real64 const> const & connRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::rateString );

    // get the well element indices corresponding to each connection
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    arrayView1d<localIndex const> const & prevWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::prevWellElementIndexString );
    
    // get constitutive data on well elements to compute the normalizer
    SingleFluidBase const * const fluid = GetFluidModel( wellElementSubRegion );
    
    array2d<real64 const> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );
    
    // set the normalizer to zero
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem) 
      massBalanceNormalizer[iwelem] = 0.0; // TODO: ask for a better way to do that!

    // loop over all connections and assign a normalizer to the well elements
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // get previous segment index
      localIndex const iwelemPrev =  prevWellElemIndex[iconn];

      // get second segment index
      localIndex const iwelemNext = nextWellElemIndex[iconn];
      
      if (iwelemPrev >= 0)
      {
	if (massBalanceNormalizer[iwelemPrev] < std::numeric_limits<real64>::epsilon())
	{
	  // get the normalizer
          real64 const absConnRate = fabs( connRate[iconn] );
	  real64 const normalizer = dt * wellElemDensity[iwelemPrev][0] * absConnRate;
	  if (normalizer > std::numeric_limits<real64>::epsilon())
	    massBalanceNormalizer[iwelemPrev] = normalizer;
	}
      }
      if (iwelemNext >= 0)
      {
	if (massBalanceNormalizer[iwelemNext] < std::numeric_limits<real64>::epsilon())
	{
	  // get the normalizer
	  real64 const absConnRate = fabs( connRate[iconn] );
	  real64 const normalizer = dt * wellElemDensity[iwelemNext][0] * absConnRate;
	  if (normalizer > std::numeric_limits<real64>::epsilon())
	    massBalanceNormalizer[iwelemNext] = normalizer;
	}
      }
    }
  });
}

void
SinglePhaseWell::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                    DomainPartition * const domain,
                                    EpetraBlockSystem * const blockSystem )
{
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain, dt );

  // the setup of dof numbers and linear system
  // is done in ReservoirWellSolver

}

void SinglePhaseWell::SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
						    localIndex  & numLocalRows,
                                                    globalIndex & numGlobalRows,
                                                    localIndex offset )
{
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

    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    array1d<integer> const & wellElemGhostRank =
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
        localCount += 1;
      }
      else
      {
        wellElemDofNumber[iwelem] = -1;
      }
    }
  });
  	    

  GEOS_ERROR_IF(localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}

void SinglePhaseWell::SetSparsityPattern( DomainPartition const * const domain,
                                          Epetra_FECrsGraph * const sparsity,
					  globalIndex firstWellElemDofNumber,
					  localIndex numDofPerResElement)
{
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
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get the well element indices corresponding to each connection
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

      // get the offset of the reservoir element equation
      globalIndex const resOffset  = resNDOF * resDofNumber[er][esr][ei];
      // specify the reservoir equation number
      elementLocalDofIndexRow[SubRegionTag::RES * wellNDOF] = resOffset;
      // specify the reservoir variable number
      elementLocalDofIndexCol[SubRegionTag::RES * wellNDOF] = resOffset;
      
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
	elementLocalDofIndexRow[SubRegionTag::WELL * resNDOF + idof] = elemOffset + idof;
	// specify the reservoir variable number
	elementLocalDofIndexCol[SubRegionTag::WELL * resNDOF + idof] = elemOffset + idof;
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
      
      // get previous segment index
      localIndex const iwelemPrev =  prevWellElemIndex[iconn];

      // get second segment index
      localIndex const iwelemNext = nextWellElemIndex[iconn];

      // check if this is not an entry or exit
      if (iwelemPrev < 0 || iwelemNext < 0)
	continue;

      // TODO: check ghost well element here
       
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
}

  
void SinglePhaseWell::AssembleSystem( DomainPartition * const domain,
                                      EpetraBlockSystem * const blockSystem,
                                      real64 const time_n, real64 const dt )
{
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
                                                                BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

  // first deal with the well control and switch if necessary
  CheckWellControlSwitch( domain );

  // then assemble the mass balance equations
  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleSourceTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );

  // then assemble the connection equations
  FormMomentumEquations( domain, jacobian, residual );

  // finally assemble the well control equation
  FormControlEquation( domain, jacobian, residual );
  
  //  if( verboseLevel() >= 3 )
  //{
    GEOS_LOG_RANK("After SinglePhaseWell::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  //}
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
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get a reference to the degree-of-freedom numbers and to ghosting info
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
    
    // get the normalizer for the mass balance equations on well elements
    array1d<real64 const> const & massBalanceNormalizer =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::massBalanceNormalizerString );
    
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

    /* @Francois try to remember
     *  Q > 0 flow from prev to next; Q < 0 flow from next to prev
     *  We assume that prev is the top segment / next is to bottom segment
     *  With this assumption, production implies Q < 0 and injection Q > 0
     *  Remember that in the residual, +Q in prev and -Q in next
     */
    
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      localIndex const iwelemNext = nextWellElemIndex[iconn];
      localIndex const iwelemPrev = prevWellElemIndex[iconn];
      
      if (iwelemNext < 0) // nothing to do
	continue;
      else 
      {
	// get well secondary variables on segments
        SingleFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

        array2d<real64> const & wellElemDensity =
          fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

        array2d<real64> const & dWellElemDensity_dPres =
          fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

        // 1) Compute the flux and its derivatives
	
	real64 const currentConnRate = connRate[iconn] + dConnRate[iconn];
        real64 const volFlux = dt * currentConnRate;
        real64 const dVolFlux_dRate = dt;

	real64 densityUp = 0.0;
	real64 dDensityUp_dPresPrev = 0.0;
	real64 dDensityUp_dPresNext = 0.0;
	
        if (iwelemPrev < 0 || // exit connection
	    currentConnRate < 0) // iwelemNext is upstream
	{	
	  densityUp = wellElemDensity[iwelemNext][0];
	  dDensityUp_dPresNext = dWellElemDensity_dPres[iwelemNext][0];
	}
	else // iwelemNext is downstream
	{
	  densityUp = wellElemDensity[iwelemPrev][0];
	  dDensityUp_dPresPrev = dWellElemDensity_dPres[iwelemPrev][0];
	}

        real64 const flux = densityUp * volFlux;
	real64 const dFlux_dRate = densityUp * dVolFlux_dRate;
	real64 const dFlux_dPresPrev = dDensityUp_dPresPrev * volFlux;
	real64 const dFlux_dPresNext = dDensityUp_dPresNext * volFlux;

	// 2) Assemble the flux into residual and Jacobian

        // TODO: check for ghost well elements
			 
	if ( iwelemPrev < 0)
	{
          real64 normalizer = 1;
	  if (m_normalizeMassBalanceEqnsFlag
              && massBalanceNormalizer[iwelemNext] >= std::numeric_limits<real64>::epsilon())
            normalizer = 1.0 / massBalanceNormalizer[iwelemNext]; 

	  // flux terms
          real64 const localFlux = - flux * normalizer;
	  real64 const localFluxJacobian_dRate = - dFlux_dRate * normalizer;
	  real64 const localFluxJacobian_dPres = - dFlux_dPresNext * normalizer; 

          // jacobian indices
	  globalIndex const offset = getElementOffset( wellElemDofNumber[iwelemNext] );
          globalIndex const eqnRowIndex = offset + RowOffset::MASSBAL;
	  globalIndex const dofColIndex_dPres = offset + ColOffset::DPRES;
	  globalIndex const dofColIndex_dRate = offset + ColOffset::DRATE;

          // add contribution to global residual and jacobian
          residual->SumIntoGlobalValues( 1, &eqnRowIndex, &localFlux );
          jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
					 1, &dofColIndex_dPres,
					 &localFluxJacobian_dPres );
          jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
					 1, &dofColIndex_dRate,
					 &localFluxJacobian_dRate );
	  
	}
	else
	{
	  
          // local working variables and arrays
          stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
          globalIndex dofColIndex_dRate;
          stackArray1d<globalIndex, 2> dofColIndices_dPres( 2 );
	  
          stackArray1d<real64, 2> localFlux( 2 );
          stackArray1d<real64, 2> localFluxJacobian_dRate( 2 );
	  stackArray2d<real64, 4> localFluxJacobian_dPres( 2, 2);

	  stackArray1d<real64, 2> normalizer( 2 );
          normalizer = 1.0;
	  
	  if (m_normalizeMassBalanceEqnsFlag)
	  {
	    if (massBalanceNormalizer[iwelemPrev] >= std::numeric_limits<real64>::epsilon())
              normalizer[ElemTag::PREV] = 1.0 / massBalanceNormalizer[iwelemPrev]; 
	    if (massBalanceNormalizer[iwelemNext] >= std::numeric_limits<real64>::epsilon())
              normalizer[ElemTag::NEXT] = 1.0 / massBalanceNormalizer[iwelemNext]; 
	  }
	  
	  // flux terms	  
          localFlux[ElemTag::PREV] =   flux * normalizer[ElemTag::PREV];
	  localFlux[ElemTag::NEXT] = - flux * normalizer[ElemTag::NEXT];

          localFluxJacobian_dRate[ElemTag::PREV] =   dFlux_dRate * normalizer[ElemTag::PREV];
	  localFluxJacobian_dRate[ElemTag::NEXT] = - dFlux_dRate * normalizer[ElemTag::NEXT];

          localFluxJacobian_dPres[ElemTag::PREV][ElemTag::PREV] =   dFlux_dPresPrev * normalizer[ElemTag::PREV];
          localFluxJacobian_dPres[ElemTag::PREV][ElemTag::NEXT] =   dFlux_dPresNext * normalizer[ElemTag::PREV];
          localFluxJacobian_dPres[ElemTag::NEXT][ElemTag::PREV] = - dFlux_dPresPrev * normalizer[ElemTag::NEXT];
          localFluxJacobian_dPres[ElemTag::NEXT][ElemTag::NEXT] = - dFlux_dPresNext * normalizer[ElemTag::NEXT];

          // indices
	  globalIndex const offsetNext = getElementOffset( wellElemDofNumber[iwelemNext] );
	  globalIndex const offsetPrev = getElementOffset( wellElemDofNumber[iwelemPrev] );
	  eqnRowIndices[ElemTag::PREV] = offsetPrev + RowOffset::MASSBAL;
	  eqnRowIndices[ElemTag::NEXT] = offsetNext + RowOffset::MASSBAL;
	  dofColIndex_dRate = offsetNext + ColOffset::DRATE; // this follows the convention used elsewhere
	  dofColIndices_dPres[ElemTag::PREV] = offsetPrev + ColOffset::DPRES;
	  dofColIndices_dPres[ElemTag::NEXT] = offsetNext + ColOffset::DPRES;

	  // Add to global residual/jacobian
          residual->SumIntoGlobalValues( integer_conversion<int>( 2 ),
				         eqnRowIndices.data(),
				         localFlux.data() );
	  // Add rate derivatives
	  jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 ),
					 eqnRowIndices.data(),
					 integer_conversion<int>( 1 ),
					 &dofColIndex_dRate,
                                         localFluxJacobian_dRate.data(),
					 Epetra_FECrsMatrix::ROW_MAJOR);					 
          jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 ),
                                         eqnRowIndices.data(),
                                         integer_conversion<int>( 2 ),
                                         dofColIndices_dPres.data(),
                                         localFluxJacobian_dPres.data(),
                                         Epetra_FECrsMatrix::ROW_MAJOR);
	}
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
  WellManager * const wellManager = domain->getWellManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber = m_resDofNumber;
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData const * const perforationData = well->getPerforations();
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // compute the local rates for this well
    ComputeAllPerforationRates( well );
    
    // get the degrees of freedom 
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    array1d<real64 const> const & massBalanceNormalizer =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::massBalanceNormalizerString );
    
    // get well variables on perforations
    array1d<real64 const> const & perfRate =
      perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

    array2d<real64 const> const & dPerfRate_dPres =
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
      globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );
      
      // row index on reservoir side
      eqnRowIndices[SubRegionTag::RES] = resDofNumber[er][esr][ei];

      // column index on reservoir side
      dofColIndices[SubRegionTag::RES] = resDofNumber[er][esr][ei];

      // row index on well side
      eqnRowIndices[SubRegionTag::WELL] = elemOffset + RowOffset::MASSBAL;
      
      // column index on well side
      dofColIndices[SubRegionTag::WELL] = elemOffset + ColOffset::DPRES;

      // get the normalizer for the mass balance equation
      real64 normalizer = 1.0;
      if (m_normalizeMassBalanceEqnsFlag
          && massBalanceNormalizer[iwelem] >= std::numeric_limits<real64>::epsilon())
        normalizer = 1.0 / massBalanceNormalizer[iwelem]; 
      
      // populate local flux vector and derivatives
      localFlux[SubRegionTag::RES]  =  dt * perfRate[iperf];
      localFlux[SubRegionTag::WELL] = -localFlux[SubRegionTag::RES] * normalizer;

      for (localIndex ke = 0; ke < 2; ++ke)
      {
        localFluxJacobian[SubRegionTag::RES][ke]  = dt * dPerfRate_dPres[iperf][ke];
        localFluxJacobian[SubRegionTag::WELL][ke] = - localFluxJacobian[SubRegionTag::RES][ke] * normalizer;
      }

      // Add to global residual/jacobian
      residual->SumIntoGlobalValues( integer_conversion<int>( 2 ),
				     eqnRowIndices.data(),
				     localFlux.data() );

      jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 ),
                                     eqnRowIndices.data(),
                                     integer_conversion<int>( 2 ),
                                     dofColIndices.data(),
                                     localFluxJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR);
    }
    
  });
}


void SinglePhaseWell::FormMomentumEquations( DomainPartition * const domain,
                                             Epetra_FECrsMatrix * const jacobian,
                                             Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData const * const connectionData = well->getConnections();
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom 
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get pressure data
    array1d<real64 const> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );
    
    // get the well element indices corresponding to each connect
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    arrayView1d<localIndex const> const & prevWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::prevWellElementIndexString );    

    // get well secondary variables on segments
    SingleFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

    array2d<real64 const> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

    array2d<real64 const> const & dWellElemDensity_dPres =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

    
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      localIndex const iwelemNext = nextWellElemIndex[iconn];
      localIndex const iwelemPrev = prevWellElemIndex[iconn];

      // TODO: check for ghost well elements
      
      if ( iwelemNext < 0 || // exit: nothing to do
	   iwelemPrev < 0 )  // first connection: form control equation, not momentum
	continue;
      else 
      {
        // local working variables and arrays
        stackArray1d<globalIndex, 2> dofColIndices( 2 );
        stackArray1d<real64, 2> localFluxJacobian( 2 );
 
        // compute avg density
        real64 const avgDensity = 0.5 * ( wellElemDensity[iwelemPrev][0] + wellElemDensity[iwelemNext][0] );
        real64 const dAvgDensity_dPresPrev = 0.5 * dWellElemDensity_dPres[iwelemPrev][0];
        real64 const dAvgDensity_dPresNext = 0.5 * dWellElemDensity_dPres[iwelemNext][0];

        // compute depth diff times acceleration
	real64 const gravD = ( wellElemGravDepth[iwelemPrev] - wellElemGravDepth[iwelemNext] );

	// compute the current pressure in the two well elements
	real64 const pressurePrev = wellElemPressure[iwelemPrev] + dWellElemPressure[iwelemPrev];
	real64 const pressureNext = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];

        // compute a coefficient to normalize the momentum equation
	real64 normalizer = (pressurePrev > pressureNext)
			  ? pressurePrev 
			  : pressureNext;
        normalizer = 1. / normalizer;
	
	// compute momentum flux and derivatives
        real64 const localFlux = ( pressurePrev - pressureNext - avgDensity * gravD ) * normalizer;
        localFluxJacobian[ElemTag::PREV] = ( 1 - dAvgDensity_dPresPrev * gravD ) * normalizer;
	localFluxJacobian[ElemTag::NEXT] = (-1 - dAvgDensity_dPresNext * gravD ) * normalizer;

        // TODO: add friction and acceleration terms
	
	// jacobian indices
	globalIndex const offsetPrev = getElementOffset( wellElemDofNumber[iwelemPrev] ); 
	globalIndex const offsetNext = getElementOffset( wellElemDofNumber[iwelemNext] ); 
        globalIndex const eqnRowIndex = offsetNext + RowOffset::CONTROL; // according to convention used in this class
	dofColIndices[ElemTag::PREV]  = offsetPrev + ColOffset::DPRES;
	dofColIndices[ElemTag::NEXT]  = offsetNext + ColOffset::DPRES;

        // add contribution to global residual and jacobian
        residual->SumIntoGlobalValues( 1, &eqnRowIndex, &localFlux );
        jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
				       integer_conversion<int>( 2 ), dofColIndices.data(),
				       localFluxJacobian.data() );
      }
    }
  });
}

void SinglePhaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {    
    // if isViable is true at the end of the following checks, no need to switch
    bool controlIsViable = false;

    // get well control and type
    Well::Control const currentControl = well->getControl();
    Well::Type const type = well->getType();
    
    // BHP control
    if (currentControl == Well::Control::BHP)
    {
      PerforationData * const perforationData = well->getPerforations();

      arrayView1d<real64 const> const & perfRate =
        perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

      // compute the local rates for this well
      // TODO: this is a waste of computations, we only need to know the sign of the potential difference
      ComputeAllPerforationRates( well );

      // the control is viable if at least one of the  perforations can inject / produce without cross flow
      // TODO: this loop will require communication
      for (localIndex iperf = 0; !controlIsViable && iperf < perforationData->numPerforationsLocal(); ++iperf)
      {	
	// producer should have positive difference; injector should have negative difference
	controlIsViable = ( ( type == Well::Type::PRODUCER && perfRate[iperf] > 0 ) ||
                            ( type == Well::Type::INJECTOR && perfRate[iperf] < 0 ) );		
      }
    }
    else // rate control
    {
      WellElementSubRegion * wellElementSubRegion = well->getWellElements();
      
      // get pressure data
      array1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      array1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // again we assume here that the first segment is on this MPI rank
      localIndex const iwelemRef = 0;
      
      // the control is viable if the reference pressure is below/above the max/min pressure
      real64 const refPressure = wellElemPressure[iwelemRef] + dWellElemPressure[iwelemRef];

      if ( type == Well::Type::PRODUCER )
      {
        real64 const minPressure = well->getTargetBHP(); // targetBHP specifies a min pressure here
	controlIsViable = ( refPressure >= minPressure );
      }
      else
      {
        real64 const maxPressure = well->getTargetBHP(); // targetBHP specifies a max pressure here
	controlIsViable = ( refPressure <= maxPressure );
      }
    }

    if (!controlIsViable)
    {
      if ( currentControl == Well::Control::BHP )
      {
	well->setControl( Well::Control::OILRATE );
        if ( m_verboseLevel >= 1 )
        {
          GEOS_LOG_RANK_0( "Control switch for well " << well->getName()
			   << " from BHP constraint to rate constraint" );
        }
      }
      else // rate control
      {
	well->setControl( Well::Control::BHP );
	if ( m_verboseLevel >= 1 )
        {
          GEOS_LOG_RANK_0( "Control switch for well " << well->getName()
			   << " from rate constraint to BHP constraint" );
        }
      }
    }
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
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
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

  return true;
}

void
SinglePhaseWell::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
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
	std::cout << "pressure: local_solution " << local_solution[lid] << std::endl;
        dWellElemPressure[iwelem] += scalingFactor * local_solution[lid];
      }
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // as a temporary solution, the connection variable is looked up with the next element's DOF number
      localIndex const iwelemNext = nextWellElemIndex[iconn];

      // check if there is a variable defined here
      if (iwelemNext < 0)
        continue;

      if (wellElemGhostRank[iwelemNext] < 0)
      {      
        // extract solution and apply to dQ
	globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelemNext] );
        int const lid = rowMap->LID( integer_conversion<int>( elemOffset + ColOffset::DRATE ) );
	std::cout << "rate: local_solution " << local_solution[lid] << std::endl;
        dConnRate[iconn] += scalingFactor * local_solution[lid];
      }
    }
  });  

  // TODO: call CommunicationTools::SynchronizeFields

  // update properties
  UpdateStateAll( domain );

}

void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  
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

}

void SinglePhaseWell::ResetViews(DomainPartition * const domain)
{
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
  array1d<globalIndex const> const & wellElemDofNumber =
    wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

  // get well primary variables on segments
  array1d<real64> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  array1d<real64> const & dWellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  array1d<real64 const> const & wellElemGravDepth =
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
    eqnRowIndices[SubRegionTag::RES] = resDofNumber[er][esr][ei];

    // column index on reservoir side
    dofColIndices[SubRegionTag::RES] = resDofNumber[er][esr][ei];

    // get reservoir variables
    density[SubRegionTag::RES] = resDensity[er][esr][m_fluidIndex][ei][0];
    dDensity_dP[SubRegionTag::RES] = dResDensity_dPres[er][esr][m_fluidIndex][ei][0];

    viscosity[SubRegionTag::RES] = resViscosity[er][esr][m_fluidIndex][ei][0];
    dViscosity_dP[SubRegionTag::RES] = dResViscosity_dPres[er][esr][m_fluidIndex][ei][0];

    // TODO: ask about Ecl/Intersect treatment here
    // In AD-GPRS, depth difference between the reservoir elem center and perf is not accounted for
    pressure[SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[SubRegionTag::RES] = 1;
    
    multiplier[SubRegionTag::RES] = 1;

    // 2) Well side

    localIndex const iwelem = perfWellElemIndex[iperf]; 

    // row index on well side
    eqnRowIndices[SubRegionTag::WELL] = wellElemDofNumber[iwelem] + RowOffset::MASSBAL;

    // column index on well side
    dofColIndices[SubRegionTag::WELL] = wellElemDofNumber[iwelem] + ColOffset::DPRES;

    // get well variables
    density[SubRegionTag::WELL] = wellElemDensity[iwelem][0];
    dDensity_dP[SubRegionTag::WELL] = dWellElemDensity_dPres[iwelem][0];

    viscosity[SubRegionTag::WELL] = wellElemViscosity[iwelem][0];
    dViscosity_dP[SubRegionTag::WELL] = dWellElemViscosity_dPres[iwelem][0];

    pressure[SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    pressure[SubRegionTag::WELL] += m_gravityFlag
                                  ? density[SubRegionTag::WELL] * ( perfGravDepth[iperf] - wellElemGravDepth[iwelem] )
                                  : 0; 		
    dPressure_dP[SubRegionTag::WELL] = 1;
    dPressure_dP[SubRegionTag::WELL] += m_gravityFlag
                                      ? dDensity_dP[SubRegionTag::WELL] * ( perfGravDepth[iperf] - wellElemGravDepth[iwelem] )
                                      : 0;
   
    multiplier[SubRegionTag::WELL] = -1;

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
    localIndex const k_up = (potDif >= 0) ? SubRegionTag::RES : SubRegionTag::WELL;

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
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // TODO: check that the first connection is on my rank
    
    // the well control equation is defined on the first connection
    localIndex const iconnControl = 0;
    // again we assume here that the first segment is on this MPI rank
    localIndex const iwelemControl = 0;

    // get well control
    Well::Control const control = well->getControl();

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
      real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
	                      ? 1.0 / targetBHP
	                      : 1.0;
      
      real64 const controlEqn = ( currentBHP - targetBHP ) * normalizer;
      real64 const dControlEqn_dPres = normalizer;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DPRES;
      
      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     &controlEqn );
      
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     1, &dofColIndex,
                                     &dControlEqn_dPres );
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
      real64 const normalizer      = targetConnRate > std::numeric_limits<real64>::min()
	                           ? 1.0 / targetConnRate
	                           : 1.0;

      real64 const controlEqn = ( currentConnRate - targetConnRate ) * normalizer;
      real64 const dControlEqn_dRate = normalizer;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DRATE;

      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     &controlEqn );
      
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     1, &dofColIndex,
                                     &dControlEqn_dRate );
    }
  });

}


void SinglePhaseWell::ImplicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
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

    RecordWellData( well );
  });

}


void SinglePhaseWell::RecordWellData( Well * well )
{
  ConnectionData * connectionData = well->getConnections();
  WellElementSubRegion * wellElementSubRegion = well->getWellElements();
  PerforationData * perforationData = well->getPerforations();

  array1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  array1d<real64 const> const & connRate =
    connectionData->getReference<array1d<real64>>( viewKeyStruct::rateString );

  array1d<real64 const> const & perfRate =
    perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

  // here, we will save the well info
  // for now, brute force: output to terminal

  std::cout << "Well : " << well->getName() << std::endl;
  if (well->getType() == Well::Type::PRODUCER)
    std::cout << "Type : PRODUCER" << std::endl;
  else 
    std::cout << "Type : INJECTOR" << std::endl;
  if (well->getControl() == Well::Control::BHP)
    std::cout << "Control : BHP" << std::endl;
  else
    std::cout << "Control : RATE" << std::endl;

  std::cout << "Below, positive perforation rate means flow from reservoir to well" << std::endl;
  std::cout << "Negative perforation rate means flow from well to reservoir" << std::endl;
  
  // output perforation rates
  for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
  {
    std::cout << "Rate at perforation #" << iperf << ": " << perfRate[iperf] << std::endl;
  }

  // output the reference pressure
  localIndex const iwelemRef = 0; // be careful here for the parallel case
  real64 const pressure = wellElemPressure[iwelemRef];
  real64 const targetPressure = well->getTargetBHP();

  if (well->getControl() == Well::Control::BHP)
  {
    std::cout << "Current reference pressure = " << pressure
  	      << ", targetPressure = "           << targetPressure
	      << std::endl;
  }
  else
  {
    if (well->getType() == Well::Type::PRODUCER)
    {
      std::cout << "Current reference pressure = " << pressure
		<< ", min pressure = " << targetPressure
		<< std::endl;
    }
    else
    {
      std::cout << "Current reference pressure = " << pressure
		<< ", max pressure = " << targetPressure
		<< std::endl;
    }
  }

  std::cout << "Below, negative connection rate means production" << std::endl;
  std::cout << "Positive connection rate means injection" << std::endl;
  
  for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
  {
    if (well->getControl() == Well::Control::BHP || iconn > 0)
      std::cout << "Rate at connection #" << iconn << ": " << connRate[iconn] << std::endl;
    else
      std::cout << "Rate at connection #" << iconn << ": " << connRate[iconn]
		<< ", target rate : " << well->getTargetRate() << std::endl;
  }
}

REGISTER_CATALOG_ENTRY(SolverBase, SinglePhaseWell, string const &, ManagedGroup * const)
}// namespace geosx
