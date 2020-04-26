/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
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
                                Group * const parent )
  : SolverBase( name, parent ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 )
{
  this->registerWrapper( viewKeyStruct::fluidNamesString, &m_fluidModelNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of fluid constitutive object to use for this solver." );
}

Group * WellSolverBase::CreateChild( string const & childKey, string const & childName )
{
  Group * rval = nullptr;

  if( childKey == keys::wellControls )
  {
    rval = RegisterGroup< WellControls >( childName );
  }
  else
  {
    SolverBase::CreateChild( childKey, childName );
  }
  return rval;
}

void WellSolverBase::ExpandObjectCatalogs()
{
  CreateChild( keys::wellControls, keys::wellControls );
}


WellSolverBase::~WellSolverBase() = default;

void WellSolverBase::PostProcessInput()
{
  SolverBase::PostProcessInput();
  CheckModelNames( m_fluidModelNames, viewKeyStruct::fluidNamesString );
}

void WellSolverBase::RegisterDataOnMesh( Group * const meshBodies )
{
  SolverBase::RegisterDataOnMesh( meshBodies );

  MeshLevel & meshLevel = *meshBodies->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    PerforationData * const perforationData = subRegion.GetPerforationData();
    perforationData->registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString );
  } );
}

void WellSolverBase::SetupDofs( DomainPartition const * const domain,
                                DofManager & dofManager ) const
{
  MeshLevel const & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  array1d< string > regions;
  forTargetRegions< WellElementRegion >( meshLevel, [&]( localIndex const,
                                                         WellElementRegion const & region )
  {
    regions.push_back( region.getName() );
  } );

  dofManager.addField( WellElementDofName(),
                       DofManager::Location::Elem,
                       NumDofPerWellElement(),
                       regions );

  dofManager.addCoupling( WellElementDofName(),
                          WellElementDofName(),
                          DofManager::Connector::Node );
}

void WellSolverBase::ImplicitStepSetup( real64 const & time_n,
                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                        DomainPartition * const domain,
                                        DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                        ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                        ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                        ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // Initialize the primary and secondary variables for the first time step
  if( time_n <= 0.0 )
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

  // then assemble the mass balance equations (exception perforation terms)
  matrix.open();
  rhs.open();

  // then assemble the mass balance equations
  AssembleFluxTerms( time, dt, domain, &dofManager, &matrix, &rhs );

  // then assemble the volume balance equations
  AssembleVolumeBalanceTerms( time, dt, domain, &dofManager, &matrix, &rhs );

  // then assemble the pressure relations between well elements
  FormPressureRelations( domain, &dofManager, &matrix, &rhs );

  // finally assemble the well control equation
  FormControlEquation( domain, &dofManager, &matrix, &rhs );

  matrix.close();
  rhs.close();

  if( getLogLevel() >= 2 )
  {
    GEOSX_LOG_RANK_0( "After WellSolverBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }
}

void WellSolverBase::UpdateStateAll( DomainPartition & domain )
{

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    UpdateState( subRegion, targetIndex );
  } );

}

void WellSolverBase::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping( *meshLevel.getElemManager(), m_fluidModelNames );
  }

  FlowSolverBase const * const flowSolver = getParent()->GetGroup< FlowSolverBase >( GetFlowSolverName() );
  m_numDofPerResElement = flowSolver->numDofPerCell();
}

void WellSolverBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  forTargetSubRegions< WellElementSubRegion >( mesh, [&]( localIndex const,
                                                          WellElementSubRegion & subRegion )
  {
    subRegion.ReconstructLocalConnectivity();
  } );

  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // Precompute solver-specific constant data (e.g. gravity-coefficient)
  PrecomputeData( domain );
}

void WellSolverBase::PrecomputeData( DomainPartition * const domain )
{
  R1Tensor const gravVector = gravityVector();
  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    PerforationData * const perforationData = subRegion.GetPerforationData();

    arrayView1d< R1Tensor const > const & wellElemLocation =
      subRegion.getReference< array1d< R1Tensor > >( ElementSubRegionBase::viewKeyStruct::elementCenterString );

    arrayView1d< real64 > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    arrayView1d< R1Tensor const > const & perfLocation =
      perforationData->getReference< array1d< R1Tensor > >( PerforationData::viewKeyStruct::locationString );

    arrayView1d< real64 > const & perfGravCoef =
      perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      // precompute the depth of the well elements
      wellElemGravCoef[iwelem] = Dot( wellElemLocation[iwelem], gravVector );
    }

    forAll< serialPolicy >( perforationData->size(), [=]( localIndex const iperf )
    {
      // precompute the depth of the perforations
      perfGravCoef[iperf] = Dot( perfLocation[iperf], gravVector );
    } );


    // set the first well element of the well
    if( subRegion.IsLocallyOwned() )
    {

      localIndex const iwelemControl = subRegion.GetTopWellElementIndex();

      GEOSX_ERROR_IF( iwelemControl < 0,
                      "Invalid well definition: well " << subRegion.getName()
                                                       << " has no well head" );

      // save the index of reference well element (used to enforce constraints)
      WellControls & wellControls = GetWellControls( subRegion );
      wellControls.SetReferenceWellElementIndex( iwelemControl );
    }
  } );
}

void WellSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_resGravCoef =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( FlowSolverBase::viewKeyStruct::gravityCoefString );
}

WellControls & WellSolverBase::GetWellControls( WellElementSubRegion const & subRegion )
{
  string const & name = subRegion.GetWellControlsName();

  WellControls * wellControls = this->GetGroup< WellControls >( name );
  GEOSX_ERROR_IF( wellControls == nullptr, "Well constraint " + name + " not found" );

  return *wellControls;
}

WellControls const & WellSolverBase::GetWellControls( WellElementSubRegion const & subRegion ) const
{
  string const & name = subRegion.GetWellControlsName();

  WellControls const * wellControls = this->GetGroup< WellControls >( name );
  GEOSX_ERROR_IF( wellControls == nullptr, "Well constraint " + name + " not found" );

  return *wellControls;
}


} // namespace geosx
