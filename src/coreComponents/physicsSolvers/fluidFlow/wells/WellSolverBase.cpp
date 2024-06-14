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
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/ThermalCompositionalMultiphaseWellKernels.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

WellSolverBase::WellSolverBase( string const & name,
                                Group * const parent )
  : SolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_isThermal( 0 ),
  m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), name + "_rates" ))
{
  registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not." );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  this->registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Write rates into a CSV file" );
}

Group *WellSolverBase::createChild( string const & childKey, string const & childName )
{
  Group *rval = nullptr;

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

  // 1. Set key dimensions of the problem
  m_numDofPerWellElement = m_isThermal ?    m_numComponents + 2 : m_numComponents + 1; // 1 pressure  connectionRate + temp if thermal
  m_numDofPerResElement = m_isThermal ? m_numComponents  + 1: m_numComponents;   // 1 pressure   + temp if thermal


  // create dir for rates output
  if( m_writeCSV > 0 )
  {
    if( MpiWrapper::commRank() == 0 )
    {
      makeDirsForPath( m_ratesOutputDir );
    }
    // wait till the dir is created by rank 0
    MPI_Barrier( MPI_COMM_WORLD );
  }
}

void WellSolverBase::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  //this->registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
  //setApplyDefaultValue( 0 ).
  //setInputFlag( InputFlags::OPTIONAL ).
  //  setDescription( "Flag indicating whether the problem is thermal or not. Set isThermal=\"1\" to enable the thermal coupling" );

  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & meshLevel,
                                                   arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            WellElementSubRegion & subRegion )
    {

      subRegion.registerField< fields::well::gravityCoefficient >( getName() );

      subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      PerforationData * const perforationData = subRegion.getPerforationData();
      perforationData->registerField< fields::well::gravityCoefficient >( getName() );
    } );
  } );
}

void WellSolverBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  SolverBase::setConstitutiveNamesCallSuper( subRegion );
  subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString()).setPlotLevel( PlotLevel::NOPLOT ).setRestartFlags( RestartFlags::NO_WRITE ).setSizedFromParent( 0 );
}

void WellSolverBase::setupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  map< std::pair< string, string >, array1d< string > > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
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
    auto const key = std::make_pair( meshBodyName, meshLevel.getName());
    meshTargets[key] = std::move( regions );
  } );

  dofManager.addField( wellElementDofName(),
                       FieldLocation::Elem,
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
  // Initialize the primary and secondary variables for the first time step

  initializeWells( domain, time_n, dt );
}


void WellSolverBase::shutInWell( real64 const time_n,
                                 real64 const dt,
                                 DomainPartition const & domain,
                                 DofManager const & dofManager,
                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      globalIndex const rankOffset = dofManager.rankOffset();

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;


        real64 unity = 1.0;
        for( integer i=0; i < m_numDofPerWellElement; i++ )
        {
          globalIndex const cindex =  dofNumber[ei] + i;
          globalIndex const rindex = localRow+i;
          localMatrix.template addToRow< serialAtomic >( rindex,
                                                         &cindex,
                                                         &unity,
                                                         1 );
          localRhs[cindex] = 0.0;
        }

      } );
    } );
  } );
}


void WellSolverBase::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    { updateSubRegionState( subRegion ); } );
  } );
}

void WellSolverBase::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    { subRegion.reconstructLocalConnectivity(); } );
  } );

  // Precompute solver-specific constant data (e.g. gravity-coefficient)
  precomputeData( domain );
}

void WellSolverBase::precomputeData( DomainPartition & domain )
{
  R1Tensor const gravVector = gravityVector();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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
      arrayView1d< real64 > const wellElemGravCoef = subRegion.getField< fields::well::gravityCoefficient >();

      arrayView2d< real64 const > const perfLocation = perforationData.getField< fields::perforation::location >();
      arrayView1d< real64 > const perfGravCoef = perforationData.getField< fields::well::gravityCoefficient >();

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
      wellControls.setReferenceGravityCoef( refElev * gravVector[2] );
    } );
  } );
}

WellControls & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

WellControls const & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion ) const
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

} // namespace geos
