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
  m_numDofPerElement = 2;
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
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::connRateString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaConnRateString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::sumPerforationRateString );
    
    PerforationData * const perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::perforationRateString );
    perforationData->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString );
    
  });    
}
  
void SinglePhaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain  = rootGroup->GetGroup<DomainPartition>( keys::domain );
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData * const perforationData = well->getPerforations();
    perforationData->getReference<array2d<real64>>( viewKeyStruct::dPerforationRate_dPresString ).resizeDimension<1>(2);
  });
}

void SinglePhaseWell::UpdateFluidModel( Well * const well )
{
  WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
  SingleFluidBase * const fluid = GetConstitutiveModel<SingleFluidBase>( wellElementSubRegion, m_fluidName );

  arrayView1d<real64 const> const & wellPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & dWellPressure =
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
    PerforationData const * const perforationData = well->getPerforations();

    // get the info stored on well elements
    arrayView1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    
    arrayView1d<real64> const & wellElemSumPerfRates =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::sumPerforationRateString );

    // get well primary variables on well elements
    arrayView1d<real64> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64> const & connRate  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );
    
    // get well data on perforations
    arrayView1d<real64> const & perfRate =
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
    SingleFluidBase const * const fluid = GetConstitutiveModel<SingleFluidBase>( wellElementSubRegion, m_fluidName );
    
    arrayView2d<real64> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );
    
    // TODO: fix this in the case of multiple MPI ranks working on the same well
    // In particular, steps 1), 3) and 7) will need to be rewritten
 
   /*
    * The methodology used here is based on Yifan Zhou's PhD dissertation, Stanford University
    */
    
    // 1) Loop over all perforations to compute an average density
    real64 avgDensity = 0;
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];
      
      avgDensity += resDensity[er][esr][m_resFluidIndex][ei][0];
    }

    // TODO: communicate avgDens
    // TODO: collect individual dens and then communicate
    
    globalIndex const numPerforationsGlobal = perforationData->numPerforationsGlobal();
    avgDensity /= numPerforationsGlobal;

    // get the reference data for this well
    localIndex const iwelemRef = 0; // here, we assume that the first well element is on this MPI rank, which will NOT be the case in general
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
                                  ? 0.9 * resPres // hard-coded values come from AD-GPRS
                                  : 1.1 * resPres;
    }

    // TODO: communicate ref pressure
    // TODO: communicate avgDens
    
    // 3) Estimate the pressures in the well elements using this avgDensity
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

    // 7) Estimate the connection rates between well elements
    // Currently this assumes that the well elements are ordered
    // TODO: implement this in parallel with the communication of wellElemSumPerfRates
    real64 prevConnRate = 0;
    for (localIndex iwelem = wellElementSubRegion->numWellElementsLocal()-1; 
         iwelem >= 0; --iwelem)
    {
      connRate[iwelem] = prevConnRate
                       - wellElemSumPerfRates[iwelem] / wellElemDensity[iwelem][0];
      prevConnRate = connRate[iwelem];
    }
  });
}

void SinglePhaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );
}

void
SinglePhaseWell::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                    DomainPartition * const domain,
                                    EpetraBlockSystem * const blockSystem )
{
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // Initialize the primary and secondary variables
  InitializeWells( domain );
  
  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );
  
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
  localIndex const resNDOF       = numDofPerResElement; 
  localIndex const wellNDOF      = numDofPerElement();

  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    PerforationData const * const perforationData = well->getPerforations();

    // get the well degrees of freedom and next well elem index
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );        
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
      
      localIndex const iwelem      = perfWellElemIndex[iperf]; 
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

    // 2) Insert the entries corresponding to fluxes between well elements
    //    This will fill J_WW only
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      
      // get next well element index
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      // check if this is not an exit
      if (iwelemNext < 0)
      {
        continue;
      }

      // TODO: check ghost well element here
       
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( 2 * wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( 2 * wellNDOF );
      
      // get the offset of the well element equations
      globalIndex const currentElemOffset = getElementOffset( wellElemDofNumber[iwelem] );
      globalIndex const nextElemOffset    = getElementOffset( wellElemDofNumber[iwelemNext] );

      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
        elementLocalDofIndexRow[idof]            = currentElemOffset + idof;
        elementLocalDofIndexRow[wellNDOF + idof] = nextElemOffset    + idof;
        elementLocalDofIndexCol[idof]            = currentElemOffset + idof;
        elementLocalDofIndexCol[wellNDOF + idof] = nextElemOffset    + idof;
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

  // then assemble the pressure relations
  FormPressureRelations( domain, jacobian, residual );

  // finally assemble the well control equation
  FormControlEquation( domain, jacobian, residual );
  
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
  // will not be implemented
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
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get a reference to the degree-of-freedom numbers and to ghosting info
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
        
    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get a reference to the primary variables on well elements
    arrayView1d<real64 const> const & connRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

    arrayView1d<real64 const> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    // get well constitutive data
    SingleFluidBase * const fluid = GetConstitutiveModel<SingleFluidBase>( wellElementSubRegion, m_fluidName );

    arrayView2d<real64> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

    arrayView2d<real64> const & dWellElemDensity_dPres =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

    // local working variables and arrays
    stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
    stackArray1d<globalIndex, 2> dofColIndices_dPres( 2 );
          
    stackArray1d<real64, 2> localFlux( 2 );
    stackArray1d<real64, 2> localFluxJacobian_dRate( 2 );
    stackArray2d<real64, 4> localFluxJacobian_dPres( 2, 2);

    // loop over the well elements to compute the fluxes between elements
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      localIndex const iwelemNext = nextWellElemIndex[iwelem];
      
      // 1) Compute the flux and its derivatives

      real64 const currentConnRate = connRate[iwelem] + dConnRate[iwelem];
      real64 const volFlux = dt * currentConnRate;
      real64 const dVolFlux_dRate = dt;

      real64 densityUp = 0.0;
      real64 dDensityUp_dPresNext    = 0.0;
      real64 dDensityUp_dPresCurrent = 0.0;

      /*  currentConnRate < 0 flow from iwelem to iwelemNext
       *  currentConnRate > 0 flow from iwelemNext to iwelem
       *  With this convention, currentConnRate < 0 at the last connection for a producer
       *                        currentConnRate > 0 at the last connection for a injector
       */
        
      if (iwelemNext < 0 ||    // exit connection
          currentConnRate < 0) // iwelem is upstream
      {       
        densityUp = wellElemDensity[iwelem][0];
        dDensityUp_dPresCurrent = dWellElemDensity_dPres[iwelem][0];
      }
      else // iwelemNext is upstream
      {
        densityUp = wellElemDensity[iwelemNext][0];
        dDensityUp_dPresNext = dWellElemDensity_dPres[iwelemNext][0];
      }

      real64 const flux = densityUp * volFlux;
      real64 const dFlux_dRate = densityUp * dVolFlux_dRate;
      real64 const dFlux_dPresCurrent = dDensityUp_dPresCurrent * volFlux;
      real64 const dFlux_dPresNext = dDensityUp_dPresNext * volFlux;

      // 2) Assemble the flux into residual and Jacobian

      // TODO: check for ghost well elements
                         
      if ( iwelemNext < 0) // exit connection
      {
        // flux terms
        real64 const oneSidedLocalFlux = - flux;
        real64 const oneSidedLocalFluxJacobian_dRate = - dFlux_dRate;
        real64 const oneSidedLocalFluxJacobian_dPres = - dFlux_dPresCurrent;

        // jacobian indices
        globalIndex const offset = getElementOffset( wellElemDofNumber[iwelem] );
        globalIndex const oneSidedEqnRowIndex = offset + RowOffset::MASSBAL;
        globalIndex const oneSidedDofColIndex_dPres = offset + ColOffset::DPRES;
        globalIndex const oneSidedDofColIndex_dRate = offset + ColOffset::DRATE;

        // add contribution to global residual and jacobian
        residual->SumIntoGlobalValues( 1, &oneSidedEqnRowIndex, &oneSidedLocalFlux );
        jacobian->SumIntoGlobalValues( 1, &oneSidedEqnRowIndex,
                                       1, &oneSidedDofColIndex_dPres,
                                       &oneSidedLocalFluxJacobian_dPres );
        jacobian->SumIntoGlobalValues( 1, &oneSidedEqnRowIndex,
                                       1, &oneSidedDofColIndex_dRate,
                                       &oneSidedLocalFluxJacobian_dRate );
          
      }
      else // not an exit connection
      {
        // local working variables and arrays
        eqnRowIndices = -1;
        dofColIndices_dPres = -1;
          
        localFlux = 0.0;
        localFluxJacobian_dRate = 0.0;
        localFluxJacobian_dPres = 0.0;

        // flux terms   
        localFlux[ElemTag::NEXT]    =   flux;
        localFlux[ElemTag::CURRENT] = - flux;

        localFluxJacobian_dRate[ElemTag::NEXT]    =   dFlux_dRate;
        localFluxJacobian_dRate[ElemTag::CURRENT] = - dFlux_dRate;

        localFluxJacobian_dPres[ElemTag::NEXT][ElemTag::NEXT]       =   dFlux_dPresNext;
        localFluxJacobian_dPres[ElemTag::NEXT][ElemTag::CURRENT]    =   dFlux_dPresCurrent;
        localFluxJacobian_dPres[ElemTag::CURRENT][ElemTag::NEXT]    = - dFlux_dPresNext;
        localFluxJacobian_dPres[ElemTag::CURRENT][ElemTag::CURRENT] = - dFlux_dPresCurrent;

        // indices
        globalIndex const offsetCurrent = getElementOffset( wellElemDofNumber[iwelem] );
        globalIndex const offsetNext    = getElementOffset( wellElemDofNumber[iwelemNext] );
        eqnRowIndices[ElemTag::NEXT]    = offsetNext + RowOffset::MASSBAL;
        eqnRowIndices[ElemTag::CURRENT] = offsetCurrent + RowOffset::MASSBAL;
        dofColIndices_dPres[ElemTag::NEXT]    = offsetNext + ColOffset::DPRES;
        dofColIndices_dPres[ElemTag::CURRENT] = offsetCurrent + ColOffset::DPRES;
        globalIndex const dofColIndex_dRate   = offsetCurrent + ColOffset::DRATE; 

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
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get well variables on perforations
    arrayView1d<real64 const> const & perfRate =
      perforationData->getReference<array1d<real64>>( viewKeyStruct::perforationRateString );

    arrayView2d<real64 const> const & dPerfRate_dPres =
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

    // local working variables and arrays
    stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
    stackArray1d<globalIndex, 2> dofColIndices( 2 );

    stackArray1d<real64, 2> localFlux( 2 );
    stackArray2d<real64, 4> localFluxJacobian(2, 2);

    // loop over the perforations and add the rates to the residual and jacobian
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      eqnRowIndices = -1;
      dofColIndices = -1;
      
      localFlux = 0;
      localFluxJacobian = 0;

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
      
      // populate local flux vector and derivatives
      localFlux[SubRegionTag::RES]  =  dt * perfRate[iperf];
      localFlux[SubRegionTag::WELL] = -localFlux[SubRegionTag::RES];

      for (localIndex ke = 0; ke < 2; ++ke)
      {
        localFluxJacobian[SubRegionTag::RES][ke]  = dt * dPerfRate_dPres[iperf][ke];
        localFluxJacobian[SubRegionTag::WELL][ke] = - localFluxJacobian[SubRegionTag::RES][ke];
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


void SinglePhaseWell::FormPressureRelations( DomainPartition * const domain,
                                             Epetra_FECrsMatrix * const jacobian,
                                             Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom numbers, depth, next well elem index
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );
    
    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    
    // get primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get well constitutive data
    SingleFluidBase * const fluid = GetConstitutiveModel<SingleFluidBase>( wellElementSubRegion, m_fluidName );

    arrayView2d<real64 const> const & wellElemDensity =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

    arrayView2d<real64 const> const & dWellElemDensity_dPres =
      fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );
    
    // local working variables and arrays
    stackArray1d<globalIndex, 2> dofColIndices( 2 );
    stackArray1d<real64, 2> localFluxJacobian( 2 );

    // loop over the well elements to compute the pressure relations between well elements
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      // TODO: check for ghost well elements
      
      if ( iwelemNext < 0 )  // if iwelemNext < 0, form control equation, not momentum
      {
        continue;
      }

      dofColIndices     = -1;
      localFluxJacobian = 0;

      // compute avg density
      real64 const avgDensity = 0.5 * ( wellElemDensity[iwelem][0] + wellElemDensity[iwelemNext][0] );
      real64 const dAvgDensity_dPresNext    = 0.5 * dWellElemDensity_dPres[iwelemNext][0];
      real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemDensity_dPres[iwelem][0];

      // compute depth diff times acceleration
      real64 const gravD = ( wellElemGravDepth[iwelemNext] - wellElemGravDepth[iwelem] );

      // compute the current pressure in the two well elements
      real64 const pressureCurrent = wellElemPressure[iwelem]     + dWellElemPressure[iwelem];
      real64 const pressureNext    = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];

      // compute a coefficient to normalize the momentum equation
      real64 normalizer = (pressureNext > pressureCurrent)
                        ? pressureNext 
                        : pressureCurrent;
      normalizer = 1. / normalizer;
        
      // compute momentum flux and derivatives
      real64 const localFlux = ( pressureNext - pressureCurrent - avgDensity * gravD ) * normalizer;
      localFluxJacobian[ElemTag::NEXT]    = ( 1 - dAvgDensity_dPresNext * gravD )    * normalizer;
      localFluxJacobian[ElemTag::CURRENT] = (-1 - dAvgDensity_dPresCurrent * gravD ) * normalizer;

      // TODO: add friction and acceleration terms
        
      // jacobian indices
      globalIndex const offsetNext     = getElementOffset( wellElemDofNumber[iwelemNext] ); 
      globalIndex const offsetCurrent  = getElementOffset( wellElemDofNumber[iwelem] ); 
      globalIndex const eqnRowIndex    = offsetCurrent + RowOffset::CONTROL; 
      dofColIndices[ElemTag::NEXT]     = offsetNext + ColOffset::DPRES;
      dofColIndices[ElemTag::CURRENT]  = offsetCurrent + ColOffset::DPRES;

      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex, &localFlux );
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     integer_conversion<int>( 2 ), dofColIndices.data(),
                                     localFluxJacobian.data() );
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
      arrayView1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      arrayView1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // again we assume here that the first well element is on this MPI rank
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
        well->setControl( Well::Control::LIQUIDRATE );
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

    // get the degree of freedom numbers and ghosting info
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

        localIndex const wellNDOF = numDofPerElement(); 
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
  // not implemented
}

bool
SinglePhaseWell::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
                                      real64 const scalingFactor,
                                      DomainPartition * const domain )
{
  // not implemented yet

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
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers on well elements and ghosting info
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
    
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
    
    arrayView1d<real64> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );
    
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract solution and apply to dP
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

        // pressure 
        int lid = rowMap->LID( integer_conversion<int>( elemOffset + ColOffset::DPRES) );
        std::cout << "pressure: local_solution " << local_solution[lid] << std::endl;
        dWellElemPressure[iwelem] += scalingFactor * local_solution[lid];

        // rate
        lid = rowMap->LID( integer_conversion<int>( elemOffset + ColOffset::DRATE ) );
        std::cout << "rate: local_solution " << local_solution[lid] << std::endl;
        dConnRate[iwelem] += scalingFactor * local_solution[lid];
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
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the ghosting information
    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
    
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView1d<real64> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        dWellElemPressure[iwelem] = 0;
        dConnRate[iwelem] = 0;
      }
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
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::densityString,
                                                                                          constitutiveManager );
  m_dResDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString,
                                                                                          constitutiveManager );
  m_resViscosity =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::viscosityString,
                                                                                          constitutiveManager );
  m_dResVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( SingleFluidBase::viewKeyStruct::dVisc_dPresString,
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

  // get the degrees of freedom and depth
  arrayView1d<globalIndex const> const & wellElemDofNumber =
    wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

  arrayView1d<real64 const> const & wellElemGravDepth =
    wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );

  // get well primary variables on well elements
  arrayView1d<real64> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64> const & dWellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
  
  // get well constitutive data
  SingleFluidBase * const fluid = GetConstitutiveModel<SingleFluidBase>( wellElementSubRegion, m_fluidName );

  arrayView2d<real64> const & wellElemDensity =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::densityString );

  arrayView2d<real64> const & dWellElemDensity_dPres =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dDens_dPresString );

  arrayView2d<real64> const & wellElemViscosity =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::viscosityString );

  arrayView2d<real64> const & dWellElemViscosity_dPres =
    fluid->getReference<array2d<real64>>( SingleFluidBase::viewKeyStruct::dVisc_dPresString );

  // get well variables on perforations
  arrayView1d<real64> const & perfGravDepth =
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

  // local working variables and arrays
  stackArray1d<globalIndex, 2> eqnRowIndices( 2 );
  stackArray1d<globalIndex, 2> dofColIndices( 2 );

  stackArray1d<real64, 2> pressure( 2 );
  stackArray1d<real64, 2> dPressure_dP( 2 );

  stackArray1d<localIndex, 2> multiplier( 2 );

  // loop over the perforations to compute the perforation rates 
  for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
  {
    eqnRowIndices = -1;
    dofColIndices = -1;
  
    pressure = 0;
    dPressure_dP = 0;

    multiplier = 0;

    // 1) Reservoir side

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    // get reservoir variables
    pressure[SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[SubRegionTag::RES] = 1;
    
    // multiplier for reservoir side in the flux 
    multiplier[SubRegionTag::RES] = 1;

    // 2) Well side

    // get the local index of the well element
    localIndex const iwelem = perfWellElemIndex[iperf]; 

    // get well variables
    pressure[SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    dPressure_dP[SubRegionTag::WELL] = 1.0;

    if (m_gravityFlag)
    {
      real64 const gravD = ( perfGravDepth[iperf] - wellElemGravDepth[iwelem] );
      pressure[SubRegionTag::WELL]     += wellElemDensity[iwelem][0] * gravD;
      dPressure_dP[SubRegionTag::WELL] += dWellElemDensity_dPres[iwelem][0] * gravD;
    }
    
    // multiplier for well side in the flux
    multiplier[SubRegionTag::WELL] = -1;

    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf]; 

    // compute potential difference
    real64 potDif = 0.0;
    for (localIndex i = 0; i < 2; ++i)
    {
      potDif += multiplier[i] * trans * pressure[i];
      dPerfRate_dPres[iperf][i] = multiplier[i] * trans * dPressure_dP[i];
    }

    // choose upstream cell based on potential difference
    localIndex const k_up = (potDif >= 0) ? SubRegionTag::RES : SubRegionTag::WELL;

    // compute upstream density, viscosity, and mobility
    real64 densityUp       = 0.0;
    real64 dDensityUp_dP   = 0.0;  
    real64 viscosityUp     = 0.0;
    real64 dViscosityUp_dP = 0.0;

    // upwinding the variables
    if (k_up == SubRegionTag::RES) // use reservoir vars
    {
      densityUp     = resDensity[er][esr][m_resFluidIndex][ei][0];
      dDensityUp_dP = dResDensity_dPres[er][esr][m_resFluidIndex][ei][0];

      viscosityUp     = resViscosity[er][esr][m_resFluidIndex][ei][0];
      dViscosityUp_dP = dResViscosity_dPres[er][esr][m_resFluidIndex][ei][0];
    }
    else // use well vars
    {
      densityUp = wellElemDensity[iwelem][0];
      dDensityUp_dP = dWellElemDensity_dPres[iwelem][0];

      viscosityUp = wellElemViscosity[iwelem][0];
      dViscosityUp_dP = dWellElemViscosity_dPres[iwelem][0];
    }

    // compute mobility
    real64 const mobilityUp     = densityUp / viscosityUp;
    real64 const dMobilityUp_dP = dDensityUp_dP / viscosityUp
                                - mobilityUp / viscosityUp * dViscosityUp_dP;

    perfRate[iperf] = mobilityUp * potDif;
    for (localIndex ke = 0; ke < 2; ++ke)
    {
      dPerfRate_dPres[iperf][ke] *= mobilityUp;
    }
    dPerfRate_dPres[iperf][k_up] += dMobilityUp_dP * potDif;
  }
}
  
void SinglePhaseWell::FormControlEquation( DomainPartition * const domain,
                                           Epetra_FECrsMatrix * const jacobian,
                                           Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom 
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // TODO: check that the first well element is on my rank
    // if not, there is nothing to do here

    // for now, we assume here that the first well element is on this MPI rank
    localIndex const iwelemControl = 0;

    // get well control
    Well::Control const control = well->getControl();

    // BHP control
    if (control == Well::Control::BHP)
    {
      // get primary variables on well elements
      arrayView1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      arrayView1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // get pressures and compute normalizer
      real64 const currentBHP = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];
      real64 const targetBHP  = well->getTargetBHP();
      real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
                              ? 1.0 / targetBHP
                              : 1.0;
      
      // control equation is a normalized difference between current pressure and target pressure
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

      // get a reference to the primary variables on well element
      arrayView1d<real64 const> const & connRate  =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

      arrayView1d<real64 const> const & dConnRate =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );
      
      // get rates and compute normalizer
      real64 const currentConnRate = connRate[iwelemControl] + dConnRate[iwelemControl];
      real64 const targetConnRate  = well->getTargetRate();
      real64 const normalizer      = targetConnRate > std::numeric_limits<real64>::min()
                                   ? 1.0 / ( 1e-2 * targetConnRate ) // hard-coded value comes from AD-GPRS
                                   : 1.0;

      // control equation is a normalized difference between current rate and target rate
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
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();
    
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & wellElemPressure  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView1d<real64> const & connRate  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

    arrayView1d<real64 const> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaConnRateString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      connRate[iwelem]         += dConnRate[iwelem];
    }

    // TODO: improve well data output 
    RecordWellData( well );
  });
}


void SinglePhaseWell::RecordWellData( Well * well )
{
  // This function is for debug and will go away

  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
  PerforationData const * perforationData = well->getPerforations();

  arrayView1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & connRate =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::connRateString );

  arrayView1d<real64 const> const & perfRate =
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
    std::cout << "Mass rate at perforation #" << iperf << ": " << perfRate[iperf] << std::endl;
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
  
  for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
  {
    if (well->getControl() == Well::Control::BHP)
      std::cout << "Volumetric rate at connection #" << iwelem << ": " << connRate[iwelem] << std::endl;
    else
      std::cout << "Volumetric rate at connection #" << iwelem << ": " << connRate[iwelem]
                << ", target rate : " << well->getTargetRate() << std::endl;
  }

}

REGISTER_CATALOG_ENTRY(SolverBase, SinglePhaseWell, string const &, ManagedGroup * const)
}// namespace geosx
