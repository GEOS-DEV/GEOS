/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file WellSolverBase.cpp
 */

#include "WellSolverBase.hpp"

#include "managers/DomainPartition.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "meshUtilities/PerforationData.hpp"

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

  this->getWrapper< string >( viewKeyStruct::discretizationString )->
    setInputFlag( InputFlags::FALSE );
}

Group * WellSolverBase::createChild( string const & childKey, string const & childName )
{
  Group * rval = nullptr;

  if( childKey == keys::wellControls )
  {
    rval = registerGroup< WellControls >( childName );
  }
  else
  {
    SolverBase::createChild( childKey, childName );
  }
  return rval;
}

void WellSolverBase::expandObjectCatalogs()
{
  createChild( keys::wellControls, keys::wellControls );
}


WellSolverBase::~WellSolverBase() = default;

void WellSolverBase::postProcessInput()
{
  SolverBase::postProcessInput();
  checkModelNames( m_fluidModelNames, viewKeyStruct::fluidNamesString );
}

void WellSolverBase::registerDataOnMesh( Group * const meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  MeshLevel & meshLevel = *meshBodies->getGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    PerforationData * const perforationData = subRegion.getPerforationData();
    perforationData->registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString );
  } );
}

void WellSolverBase::setupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  array1d< string > regions;
  forTargetRegions< WellElementRegion >( meshLevel, [&]( localIndex const,
                                                         WellElementRegion const & region )
  {
    regions.emplace_back( region.getName() );
  } );

  dofManager.addField( wellElementDofName(),
                       DofManager::Location::Elem,
                       numDofPerWellElement(),
                       regions );

  dofManager.addCoupling( wellElementDofName(),
                          wellElementDofName(),
                          DofManager::Connector::Node );
}

void WellSolverBase::implicitStepSetup( real64 const & time_n,
                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                        DomainPartition & domain )
{
  // bind the stored reservoir views to the current domain
  resetViews( domain );

  // Initialize the primary and secondary variables for the first time step
  if( time_n <= 0.0 )
  {
    initializeWells( domain );
  }

  // set deltas to zero and recompute dependent quantities
  resetStateToBeginningOfStep( domain );
}

void WellSolverBase::assembleSystem( real64 const time,
                                     real64 const dt,
                                     DomainPartition & domain,
                                     DofManager const & dofManager,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
{
  // then assemble the mass balance equations
  assembleFluxTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the volume balance equations
  assembleVolumeBalanceTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the pressure relations between well elements
  formPressureRelations( domain, dofManager, localMatrix, localRhs );
}

void WellSolverBase::updateStateAll( DomainPartition & domain )
{

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    updateState( subRegion, targetIndex );
  } );

}

void WellSolverBase::initializePreSubGroups( Group * const rootGroup )
{
  SolverBase::initializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->getGroup< DomainPartition >( keys::domain );

  for( auto & mesh : domain->getMeshBodies()->getSubGroups() )
  {
    MeshLevel & meshLevel = *Group::groupCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping( *meshLevel.getElemManager(), m_fluidModelNames );
  }

  FlowSolverBase const * const flowSolver = getParent()->getGroup< FlowSolverBase >( getFlowSolverName() );
  m_numDofPerResElement = flowSolver->numDofPerCell();
}

void WellSolverBase::initializePostInitialConditionsPreSubGroups( Group * const rootGroup )
{
  SolverBase::initializePostInitialConditionsPreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->getGroup< DomainPartition >( keys::domain );

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  forTargetSubRegions< WellElementSubRegion >( mesh, [&]( localIndex const,
                                                          WellElementSubRegion & subRegion )
  {
    subRegion.reconstructLocalConnectivity();
  } );

  // bind the stored reservoir views to the current domain
  resetViews( domain );

  // Precompute solver-specific constant data (e.g. gravity-coefficient)
  precomputeData( domain );
}

void WellSolverBase::precomputeData( DomainPartition & domain )
{
  R1Tensor const gravVector = gravityVector();
  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    PerforationData * const perforationData = subRegion.getPerforationData();

    arrayView2d< real64 const > const wellElemLocation = subRegion.getElementCenter();

    arrayView1d< real64 > const wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    arrayView2d< real64 const > const perfLocation =
      perforationData->getReference< array2d< real64 > >( PerforationData::viewKeyStruct::locationString );

    arrayView1d< real64 > const perfGravCoef =
      perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      // precompute the depth of the well elements
      wellElemGravCoef[iwelem] = wellElemLocation( iwelem, 0 ) * gravVector[ 0 ] + wellElemLocation( iwelem, 1 ) * gravVector[ 1 ] + wellElemLocation( iwelem,
                                                                                                                                                       2 ) *
                                 gravVector[ 2 ];
    }

    forAll< serialPolicy >( perforationData->size(), [=]( localIndex const iperf )
    {
      // precompute the depth of the perforations
      perfGravCoef[iperf] = LvArray::tensorOps::AiBi< 3 >( perfLocation[iperf], gravVector );
    } );


    // set the first well element of the well
    if( subRegion.isLocallyOwned() )
    {

      localIndex const iwelemControl = subRegion.getTopWellElementIndex();

      GEOSX_ERROR_IF( iwelemControl < 0,
                      "Invalid well definition: well " << subRegion.getName()
                                                       << " has no well head" );

      // save the index of reference well element (used to enforce constraints)
      WellControls & wellControls = getWellControls( subRegion );
      wellControls.setReferenceWellElementIndex( iwelemControl );
    }
  } );
}

void WellSolverBase::resetViews( DomainPartition & domain )
{
  MeshLevel * const mesh = domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh->getElemManager();

  m_resGravCoef.clear();
  m_resGravCoef = elemManager.ConstructArrayViewAccessor< real64, 1 >( FlowSolverBase::viewKeyStruct::gravityCoefString );
  m_resGravCoef.setName( getName() + "/accessors/" + FlowSolverBase::viewKeyStruct::gravityCoefString );

}

WellControls & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion )
{
  string const & name = subRegion.getWellControlsName();

  WellControls * wellControls = this->getGroup< WellControls >( name );
  GEOSX_ERROR_IF( wellControls == nullptr, "Well constraint " + name + " not found" );

  return *wellControls;
}

WellControls const & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion ) const
{
  string const & name = subRegion.getWellControlsName();

  WellControls const * wellControls = this->getGroup< WellControls >( name );
  GEOSX_ERROR_IF( wellControls == nullptr, "Well constraint " + name + " not found" );

  return *wellControls;
}

std::vector< string > WellSolverBase::getConstitutiveRelations( string const & regionName ) const
{

  localIndex const regionIndex = this->targetRegionIndex( regionName );

  std::vector< string > rval{  m_fluidModelNames[regionIndex] };

  return rval;
}



} // namespace geosx
