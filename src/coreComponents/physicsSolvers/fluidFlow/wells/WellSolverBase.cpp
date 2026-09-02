/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file WellSolverBase.cpp
 */

#include "WellSolverBase.hpp"

#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "functions/FunctionManager.hpp"
namespace geos
{

using namespace dataRepository;
using namespace fields;

WellSolverBase::WellSolverBase( string const & name,
                                Group * const parent )
  : WellControls( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_isThermal( 0 ),
  m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), name + "_rates" ) ),
  m_keepVariablesConstantDuringInitStep( false )

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
    setDescription( "When set to 1, write the rates into a CSV file." );

  this->registerWrapper( viewKeyStruct::timeStepFromTablesFlagString(), &m_timeStepFromTables ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Choose time step to honor rates/bhp tables time intervals" );


  addLogLevel< logInfo::WellControl >();
}

Group * WellSolverBase::createChild( string const & childKey, string const & childName )
{
  Group * baseChild = WellControls::createChild( childKey, childName );
  if( baseChild != nullptr )
  {
    return baseChild;
  }
  static std::set< string > const childTypes = {
    //keys::wellControls,
    PhysicsSolverBase::groupKeyStruct::linearSolverParametersString(),
    PhysicsSolverBase::groupKeyStruct::nonlinearSolverParametersString(),
  };
  GEOS_ERROR_IF( childTypes.count( childKey ) == 0,
                 CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ),
                 getDataContext() );

  PhysicsSolverBase::createChild( childKey, childName );
  return nullptr;

}

void WellSolverBase::expandObjectCatalogs()
{
  //createChild( keys::wellControls, keys::wellControls );
}

WellSolverBase::~WellSolverBase() = default;

void WellSolverBase::postInputInitialization()
{
  WellControls::postInputInitialization();

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
    MpiWrapper::barrier( MPI_COMM_WORLD );
  }
}

void WellSolverBase::registerDataOnMesh( Group & meshBodies )
{
  PhysicsSolverBase::registerDataOnMesh( meshBodies );

  // loop over the wells
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            WellElementSubRegion & subRegion )
    {
      subRegion.registerField< well::pressure >( getName() );
      subRegion.registerField< well::pressure_n >( getName() );

      subRegion.registerField< well::temperature >( getName() );
      if( isThermal() )
      {
        subRegion.registerField< well::temperature_n >( getName() );
      }

      subRegion.registerField< well::gravityCoefficient >( getName() );

      PerforationData * const perforationData = subRegion.getPerforationData();
      perforationData->registerField< well::gravityCoefficient >( getName() );
    } );
  } );
}

void WellSolverBase::initializePostSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  FunctionManager & functionManager = FunctionManager::getInstance();
  Group & meshBodies = domain.getMeshBodies();
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      validateWellConstraints( 0, 0, subRegion );

      // validate perforation status table
      PerforationData & perforationData = *subRegion.getPerforationData();
      string_array const & perfStatusTableName = perforationData.getPerfStatusTableName();
      for( integer i=0; i<perforationData.size(); i++ )
      {
        TableFunction * tableFunction =  functionManager.getGroupPointer< TableFunction >( perfStatusTableName[i] );
        GEOS_THROW_IF( tableFunction->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                       GEOS_FMT( "The interpolation method for the perforation status table {} "
                                 "should be TableFunction::InterpolationType::Lower",
                                 tableFunction->getName() ),
                       InputError, getDataContext() );
      }
    } );
  } );
}

void WellSolverBase::setupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  map< std::pair< string, string >, string_array > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                               MeshLevel const & meshLevel,
                                                               string_array const & regionNames )
  {
    string_array regions;
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



void WellSolverBase::selectWellConstraint( real64 const & time_n,
                                           real64 const & dt,
                                           const integer coupledIterationNumber,
                                           DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( dt );
  GEOS_UNUSED_VAR( coupledIterationNumber );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                               MeshLevel & meshLevel,
                                                               string_array const & regionNames )
  {
    GEOS_UNUSED_VAR( meshBodyName );
    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion & region )
    {
      WellElementSubRegion & subRegion = region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                                           .getGroup< WellElementSubRegion >( region.getSubRegionName() );
      WellControls & wellControls = getWellControls( subRegion );
      // Intiialize well if it is open
      // Well state estimated from reservoir conditions
      if( wellControls.isWellOpen() )
      {
        if( !wellControls.getWellState() )
        {
          wellControls.setWellState( 1 );

          initializeWell( domain, meshLevel, subRegion, time_n );
        }
      }
      else
      {
        wellControls.setWellState( 0 );
      }

      if( wellControls.getWellState())
      {
        evaluateConstraints( time_n,
                             subRegion );

        // If a well is opened and then timestep is cut resulting in the well being shut, if the well is opened
        // the well initialization code requires control type to by synced
        integer owner = -1;
        // Only subregion owner evaluates well control and control changes need to be broadcast to all ranks
        if( subRegion.isLocallyOwned() )
        {
          owner = MpiWrapper::commRank( MPI_COMM_GEOS );
        }
        owner = MpiWrapper::max( owner );
        WellControls::Control wellControl = wellControls.getControl();
        MpiWrapper::broadcast( wellControl, owner );
        wellControls.setControl( wellControl );
      }
    } );
  } );

}

void WellSolverBase::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion & subRegion )
    { updateSubRegionState( subRegion ); } );
  } );
}

void WellSolverBase::assembleSystem( real64 const time,
                                     real64 const dt,
                                     DomainPartition & domain,
                                     DofManager const & dofManager,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
{


  // selects constraints one of 2 ways
  //  wellEstimator flag set to 0 => orginal logic rates are computed during update state and constraints are selected every newton
  // iteration
  //  wellEstimator flag > 0 =>   well esitmator solved for each constraint and then selects the constraint
  //                         =>   estimator solve only performed first "wellEstimator" iterations
  NonlinearSolverParameters const & nonlinearParams =  getNonlinearSolverParameters();
  selectWellConstraint( time, dt, nonlinearParams.m_numNewtonIterations, domain );



  // assemble the accumulation term in the mass balance equations
  assembleAccumulationTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the pressure relations between well elements
  //assemblePressureRelations( time, dt, domain, dofManager, localMatrix, localRhs );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion & region )
    {
      WellElementSubRegion & subRegion = region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                                           .getGroup< WellElementSubRegion >( region.getSubRegionName() );
      assembleWellConstraintTerms( time, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
    } );
  } );

  // then compute the perforation rates (later assembled by the coupled solver)
  computePerforationRates( time, dt, domain );

  // then assemble the flux terms in the mass balance equations
  // get a reference to the degree-of-freedom numbers
  // then assemble the flux terms in the mass balance equations
  assembleFluxTerms( time, dt, domain, dofManager, localMatrix, localRhs );
}

void WellSolverBase::initializePostInitialConditionsPreSubGroups()
{
  PhysicsSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
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
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      PerforationData & perforationData = *subRegion.getPerforationData();
      WellControls & wellControls = getWellControls( subRegion );
      real64 const refElev = wellControls.getReferenceElevation();

      arrayView2d< real64 const > const wellElemLocation = subRegion.getElementCenter();
      arrayView1d< real64 > const wellElemGravCoef = subRegion.getField< well::gravityCoefficient >();

      arrayView2d< real64 const > const perfLocation = perforationData.getField< perforation::location >();
      arrayView1d< real64 > const perfGravCoef = perforationData.getField< well::gravityCoefficient >();

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

      wellControls.forSubGroups< MinimumBHPConstraint, MaximumBHPConstraint >( [&]( auto & constraint )
      {
        // set the reference well element where the BHP control is applied
        real64 const refElev1 = constraint.getReferenceElevation();
        constraint.setReferenceGravityCoef( refElev1 * gravVector[2] );
      } );

      // set the reference well element where the BHP control is applied
      wellControls.setReferenceGravityCoef( refElev * gravVector[2] );  // tjb remove
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

real64 WellSolverBase::setNextDt( real64 const & currentTime, const real64 & currentDt, geos::DomainPartition & domain )
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  real64 nextDt = PhysicsSolverBase::setNextDt( currentTime, currentDt, domain );

  if( m_timeStepFromTables )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            WellElementSubRegion & subRegion )
      {
        real64 nextDt_perf=nextDt;
        WellControls & wellControls = getWellControls( subRegion );
        // Find min dt from perf status tables
        PerforationData & perforationData = *subRegion.getPerforationData();
        string_array const & perfStatusTableName = perforationData.getPerfStatusTableName();

        // Get dt for local perforations
        for( integer i=0; i<perforationData.size(); i++ )
        {
          TableFunction * tableFunction =  functionManager.getGroupPointer< TableFunction >( perfStatusTableName[i] );
          WellControls::setNextDtFromTable( tableFunction, currentTime, nextDt_perf );
        }
        nextDt = MpiWrapper::min< real64 >( nextDt_perf );
        // Find min dt including rate and status tables
        real64 const nextDt_orig = nextDt;
        wellControls.setNextDtFromTables( currentTime, nextDt );
        if( m_nonlinearSolverParameters.getLogLevel() > 0 && nextDt < nextDt_orig )
          GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on tables coordinates = {}", getName(), nextDt ));
      } );
    } );
  }

  return nextDt;
}

} // namespace geos
