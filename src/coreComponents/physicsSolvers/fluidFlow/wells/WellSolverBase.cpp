/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
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

#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationData.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseExtrinsicData.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

WellSolverBase::WellSolverBase( string const & name,
                                Group * const parent )
  : SolverBase( name, parent ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_currentTime( 0 ),
  m_currentDt( 0 )
{
  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );
}

Group * WellSolverBase::createChild( string const & childKey, string const & childName )
{
  Group * rval = nullptr;

  if( childKey == keys::wellControls )
  {
    rval = &registerGroup< WellControls >( childName );
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
}

void WellSolverBase::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );


  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & meshLevel,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            WellElementSubRegion & subRegion )
    {

    subRegion.registerExtrinsicData< extrinsicMeshData::well::gravityCoefficient >( getName() );

      subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      PerforationData * const perforationData = subRegion.getPerforationData();
      perforationData->registerExtrinsicData< extrinsicMeshData::well::gravityCoefficient >( getName() );
    } );
  } );
}

void WellSolverBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );
}

void WellSolverBase::setupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  map< string, array1d< string > > meshTargets;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                MeshLevel const & meshLevel,
                                                arrayView1d< string const > const & regionNames )
  {
    array1d< string > regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[meshBodyName] = std::move( regions );
  } );

  dofManager.addField( wellElementDofName(),
                       DofManager::Location::Elem,
                       numDofPerWellElement(),
                       meshTargets );

  dofManager.addCoupling( wellElementDofName(),
                          wellElementDofName(),
                          DofManager::Connector::Node );
}

void WellSolverBase::implicitStepSetup( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition & domain )
{
  // saved time and current dt for residual normalization and time-dependent tables
  m_currentDt = dt;
  m_currentTime = time_n;

  // Initialize the primary and secondary variables for the first time step
  if( time_n <= 0.0 )
  {
    initializeWells( domain );
  }

  // set deltas to zero and recompute dependent quantities
  resetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    backupFields( mesh, regionNames );
  } );
}

void WellSolverBase::assembleSystem( real64 const time,
                                     real64 const dt,
                                     DomainPartition & domain,
                                     DofManager const & dofManager,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
{
  // assemble the accumulation term in the mass balance equations
  assembleAccumulationTerms( domain, dofManager, localMatrix, localRhs );

  // then assemble the flux terms in the mass balance equations
  assembleFluxTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the volume balance equations
  assembleVolumeBalanceTerms( domain, dofManager, localMatrix, localRhs );

  // then assemble the pressure relations between well elements
  assemblePressureRelations( domain, dofManager, localMatrix, localRhs );
}

void WellSolverBase::updateState( DomainPartition & domain )
{

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      updateSubRegionState( mesh, subRegion );
    } );
  } );


}

void WellSolverBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  FlowSolverBase const & flowSolver = getParent().getGroup< FlowSolverBase >( getFlowSolverName() );
  m_numDofPerResElement = flowSolver.numDofPerCell();
}

void WellSolverBase::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      subRegion.reconstructLocalConnectivity();
    } );
  } );

  // Precompute solver-specific constant data (e.g. gravity-coefficient)
  precomputeData( domain );
}

void WellSolverBase::precomputeData( DomainPartition & domain )
{
  R1Tensor const gravVector = gravityVector();
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      PerforationData & perforationData = *subRegion.getPerforationData();
      WellControls & wellControls = getWellControls( subRegion );
      real64 const refElev = wellControls.getReferenceElevation();

    arrayView2d< real64 const > const wellElemLocation = subRegion.getElementCenter();
    arrayView1d< real64 > const wellElemGravCoef =
      subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

    arrayView2d< real64 const > const perfLocation =
      perforationData.getReference< array2d< real64 > >( PerforationData::viewKeyStruct::locationString() );

    arrayView1d< real64 > const perfGravCoef =
      perforationData.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

      forAll< serialPolicy >( perforationData.size(), [=]( localIndex const iperf )
      {
        // precompute the depth of the perforations
        perfGravCoef[iperf] = LvArray::tensorOps::AiBi< 3 >( perfLocation[iperf], gravVector );
      } );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
      {
        // precompute the depth of the well elements
        wellElemGravCoef[iwelem] = LvArray::tensorOps::AiBi< 3 >( wellElemLocation[iwelem], gravVector );
      } );

      // set the reference well element where the BHP control is applied
      wellControls.setReferenceGravityCoef( refElev * gravVector[ 2 ] );

    } );
  } );
}

WellControls & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion )
{ return this->getGroup< WellControls >( subRegion.getWellControlsName() ); }

WellControls const & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion ) const
{ return this->getGroup< WellControls >( subRegion.getWellControlsName() ); }

} // namespace geosx
