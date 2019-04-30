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
 * @file WellSolverBase.cpp
 */

#include "WellSolverBase.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "wells/WellManager.hpp"
#include "wells/Well.hpp"
#include "wells/PerforationData.hpp"
#include "physicsSolvers/FiniteVolume/FlowSolverBase.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;
  
WellSolverBase::WellSolverBase( std::string const & name,
                                ManagedGroup * const parent )
  : SolverBase( name, parent ),
    m_gravityFlag(1),
    m_fluidName(),
    m_resFluidIndex(),
    m_numDofPerElement(0),
    m_numDofPerResElement(0),
    m_firstWellElemDofNumber(-1)
{
  RegisterViewWrapper( viewKeyStruct::gravityFlagString, &m_gravityFlag, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Flag that enables/disables gravity");

  this->RegisterViewWrapper( viewKeyStruct::fluidNameString,  &m_fluidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of fluid constitutive object to use for this solver.");

  this->RegisterViewWrapper( viewKeyStruct::resFluidIndexString, &m_resFluidIndex, false );

}

void WellSolverBase::RegisterDataOnMesh( ManagedGroup * const meshBodies )
{
  SolverBase::RegisterDataOnMesh( meshBodies );
}

void WellSolverBase::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                        DomainPartition * const domain,
                                        EpetraBlockSystem * const blockSystem )
{
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // Initialize the primary and secondary variables for the first time step
  if (time_n <= 0.0)
  {
    InitializeWells( domain );
  }

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );
  
  // assumes that the setup of dof numbers and linear system
  // is done in ReservoirWellSolver
}

void WellSolverBase::AssembleSystem( DomainPartition * const domain,
                                     EpetraBlockSystem * const blockSystem,
                                     real64 const time_n, real64 const dt )
{ 
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  // TODO: implement parallelization here
  // We want to have the possibility to solve different wells on different MPI ranks,
  // but we will first assume that all the well elements of a given well belong to the same rank

  // TODO here: gather all the property info from the perforated cells 

  // first deal with the well control and switch if necessary
  CheckWellControlSwitch( domain );

  // then assemble the mass balance equations
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssemblePerforationTerms( domain, jacobian, residual, time_n, dt );

  // then assemble the volume balance equations
  AssembleVolumeBalanceTerms( domain, jacobian, residual, time_n, dt );

  // then assemble the pressure relations between well elements
  FormPressureRelations( domain, jacobian, residual );

  // finally assemble the well control equation
  FormControlEquation( domain, jacobian, residual );
  
  // TODO here: scatter all the wells terms (perforation derivatives) to other ranks before linear solve

  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After CompositionalMultiphaseWell::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

}

void WellSolverBase::ImplicitStepComplete( real64 const & time,
                                           real64 const & dt,
                                           DomainPartition * const domain )
{
  // implemented in derived classes
}

void WellSolverBase::UpdateStateAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * const well ) -> void
  {
    UpdateState( well );
  });
}
 
void WellSolverBase::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * fluid  = cm->GetConstitituveRelation<ConstitutiveBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Fluid model " + m_fluidName + " not found" );

  m_resFluidIndex = fluid->getIndexInParent(); 
}
  
void WellSolverBase::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  
  // bind the stored reservoir views to the current domain
  ResetViews( domain );
  
  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);
}

void WellSolverBase::PrecomputeData(DomainPartition * const domain)
{
  R1Tensor const & gravityVector = getGravityVector();

  WellManager * const wellManager = domain->getWellManager();
   
  wellManager->forSubGroups<Well>( [&] ( Well * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    PerforationData const * const perforationData = well->getPerforations();

    arrayView1d<R1Tensor const> const & wellElemLocation = 
      wellElementSubRegion->getReference<array1d<R1Tensor>>( ElementSubRegionBase::viewKeyStruct::elementCenterString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    

    arrayView1d<real64> const & wellElemGravDepth = 
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );

    arrayView1d<R1Tensor const> const & perfLocation = 
      perforationData->getReference<array1d<R1Tensor>>( PerforationData::viewKeyStruct::locationString );

    arrayView1d<real64> const & perfGravDepth = 
      perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::gravityDepthString );

    localIndex iwelemControl = -1;
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      // precompute the depth of the well elements
      wellElemGravDepth[iwelem] = Dot( wellElemLocation[iwelem], gravityVector );

      // set the first well element of the well
      localIndex const iwelemNext = nextWellElemIndex[iwelem];
      if (iwelemNext < 0)
      {
        GEOS_ERROR_IF( iwelemControl >= 0, 
                      "Invalid well definition: well " << well->getName() 
                   << " has multiple well heads");
        iwelemControl = iwelem;
      }
    }
    
    // save the index of reference well element (used to enforce constraints) 
    GEOS_ERROR_IF( iwelemControl < 0, 
                  "Invalid well definition: well " << well->getName() 
               << " has no well head");
    well->setReferenceWellElementIndex( iwelemControl );

    forall_in_range( 0, perforationData->numPerforationsLocal(), GEOSX_LAMBDA ( localIndex const iperf )
    {
      // precompute the depth of the perforations
      perfGravDepth[iperf] = Dot( perfLocation[iperf], gravityVector );
    });
  });
}

void WellSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_resGravDepth =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FlowSolverBase::viewKeyStruct::gravityDepthString );
}

  
WellSolverBase::~WellSolverBase() = default;

globalIndex WellSolverBase::getElementOffset( globalIndex welemDofNumber ) const
{
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

  localIndex const resNDOF  = numDofPerResElement(); // dof is pressure
  localIndex const wellNDOF = numDofPerElement(); // dofs are pressure and rate
  
  globalIndex const firstElemDofNumber = getFirstWellElementDofNumber();
  globalIndex const currentElemOffset  = firstElemDofNumber * resNDOF // number of eqns in J_RR
                                       + (welemDofNumber - firstElemDofNumber) * wellNDOF; // number of eqns in J_WW, before this element's equations

  return currentElemOffset;
}
  
} // namespace geosx
