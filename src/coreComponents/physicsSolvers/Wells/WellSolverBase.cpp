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
#include "wells/WellControls.hpp"
#include "wells/WellElementRegion.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "wells/PerforationData.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
  
WellSolverBase::WellSolverBase( std::string const & name,
                                ManagedGroup * const parent )
  : SolverBase( name, parent ),
    m_flowSolverName(""),
    m_gravityFlag(1),
    m_fluidName(),
    m_resFluidIndex(),
    m_numDofPerWellElement(0),
    m_numDofPerResElement(0)
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

ManagedGroup * WellSolverBase::CreateChild( string const & childKey, string const & childName )
{
  ManagedGroup * rval = nullptr;

  if( childKey == keys::wellControls )
  {
    rval = RegisterGroup<WellControls>( childName );
  }
  else
  {
    GEOS_ERROR(childKey<<" is an invalid key to SolverBase child group.");
  }
  return rval;
}

WellSolverBase::~WellSolverBase() = default;

void WellSolverBase::RegisterDataOnMesh( ManagedGroup * const meshBodies )
{
  SolverBase::RegisterDataOnMesh( meshBodies );

  MeshLevel * const meshLevel = meshBodies->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    subRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::gravityDepthString );

    PerforationData * const perforationData = subRegion->GetPerforationData();
    perforationData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::gravityDepthString );
  });
}

void WellSolverBase::ImplicitStepSetup( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition * const domain,
                                        DofManager & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs,
                                        ParallelVector & solution )
{
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  FlowSolverBase const * const flowSolver = getParent()->GetGroup<FlowSolverBase>( GetFlowSolverName() );
  m_numDofPerResElement = flowSolver->numDofPerCell();

  // Initialize the primary and secondary variables for the first time step
  if (time_n <= 0.0)
  {
    InitializeWells( domain );
  }

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );
}

void WellSolverBase::AssembleSystem( real64 const time,
                                     real64 const dt,
                                     DomainPartition * const domain,
                                     DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs )
{ 

  // first deal with the well control and switch if necessary
  CheckWellControlSwitch( domain );

  // then assemble the mass balance equations
  AssembleFluxTerms( time, dt, domain, &dofManager, &matrix, &rhs );
  AssemblePerforationTerms( time, dt, domain, &dofManager, &matrix, &rhs );

  // then assemble the volume balance equations
  AssembleVolumeBalanceTerms( time, dt, domain, &dofManager, &matrix, &rhs );

  // then assemble the pressure relations between well elements
  FormPressureRelations( domain, &dofManager, &matrix, &rhs );

  // finally assemble the well control equation
  FormControlEquation( domain, &dofManager, &matrix, &rhs );

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After WellSolverBase::AssembleSystem" );
    GEOS_LOG_RANK_0("\nJacobian:\n" << matrix);
    GEOS_LOG_RANK_0("\nResidual:\n" << rhs);
  }

}

void WellSolverBase::UpdateStateAll( DomainPartition * const domain )
{

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {
    UpdateState( subRegion );
  });

}
 
void WellSolverBase::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * fluid  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Fluid model " + m_fluidName + " not found" );

  m_resFluidIndex = fluid->getIndexInParent(); 
}
  
void WellSolverBase::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion * const subRegion )
  {
    subRegion->ReconstructLocalConnectivity();
  });

  // bind the stored reservoir views to the current domain
  ResetViews( domain );
  
  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);
}

void WellSolverBase::PrecomputeData(DomainPartition * const domain)
{
  R1Tensor const & gravityVector = getGravityVector();

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // loop over the wells
  elemManager->forElementSubRegions<WellElementSubRegion>([&]( WellElementSubRegion * const subRegion )
  {
    WellControls * const wellControls = GetWellControls( subRegion );

    PerforationData const * const perforationData = subRegion->GetPerforationData();

    arrayView1d<R1Tensor const> const & wellElemLocation = 
      subRegion->getReference<array1d<R1Tensor>>( ElementSubRegionBase::viewKeyStruct::elementCenterString );

    arrayView1d<real64> const & wellElemGravDepth = 
      subRegion->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

    arrayView1d<R1Tensor const> const & perfLocation = 
      perforationData->getReference<array1d<R1Tensor>>( PerforationData::viewKeyStruct::locationString );

    arrayView1d<real64> const & perfGravDepth = 
      perforationData->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

    for (localIndex iwelem = 0; iwelem < subRegion->size(); ++iwelem)
    {
      // precompute the depth of the well elements
      wellElemGravDepth[iwelem] = Dot( wellElemLocation[iwelem], gravityVector );

    }

    forall_in_range( 0, perforationData->size(), GEOSX_LAMBDA ( localIndex const iperf )
    {
      // precompute the depth of the perforations
      perfGravDepth[iperf] = Dot( perfLocation[iperf], gravityVector );
    });


    // set the first well element of the well    
    if (subRegion->IsLocallyOwned())
    {
 
      localIndex const iwelemControl = subRegion->GetTopWellElementIndex();

      GEOS_ERROR_IF( iwelemControl < 0, 
                     "Invalid well definition: well " << subRegion->getName() 
                  << " has no well head");
  
      // save the index of reference well element (used to enforce constraints) 
      wellControls->SetReferenceWellElementIndex( iwelemControl );
    }
  });
}

void WellSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_resGravDepth =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FlowSolverBase::viewKeyStruct::gravityDepthString );
}
  
WellControls * WellSolverBase::GetWellControls( WellElementSubRegion const * const subRegion )
{
  string const & name = subRegion->GetWellControlsName();

  WellControls * wellControls = this->GetGroup<WellControls>( name );
  GEOS_ERROR_IF( wellControls == nullptr, "Well constraint " + name + " not found" ); 

  return wellControls;
}

WellControls const * WellSolverBase::GetWellControls( WellElementSubRegion const * const subRegion ) const
{
  string const & name = subRegion->GetWellControlsName();

  WellControls const * wellControls = this->GetGroup<WellControls>( name );
  GEOS_ERROR_IF( wellControls == nullptr, "Well constraint " + name + " not found" ); 

  return wellControls;
}

  
} // namespace geosx
