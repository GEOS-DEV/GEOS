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
 * @file WellManager.cpp
 */

#include "WellManager.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/SolutionCheckHelpers.hpp"

#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "functions/FunctionManager.hpp"
namespace geos
{

using namespace dataRepository;
using namespace fields;

WellManager::WellManager( string const & name,
                          Group * const parent )
  : PhysicsSolverBase( name, parent ),
  m_useMass( false ),
  m_useTotalMassEquation( 1 ),
  m_isThermal( 0 ),
  m_isCompositional( true ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 )
{
  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not." );


  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::useTotalMassEquationString(), &m_useTotalMassEquation ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use total mass equation" );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );

  this->registerWrapper( viewKeyStruct::timeStepFromTablesFlagString(), &m_timeStepFromTables ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription ( "Choose time step to honor rates/bhp tables time intervals" );

}
Group * WellManager::createChild( string const & childKey, string const & childName )
{
  static std::set< string > const childTypes = {
    keys::compositionalMultiphaseWell,
    keys::singlePhaseWell,
    PhysicsSolverBase::groupKeyStruct::linearSolverParametersString(),
    PhysicsSolverBase::groupKeyStruct::nonlinearSolverParametersString(),
  };
  GEOS_ERROR_IF( childTypes.count( childKey ) == 0,
                 CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ),
                 getDataContext() );
  if( childKey == keys::compositionalMultiphaseWell )
  {
    setCompositional( true );
    return &registerGroup< CompositionalMultiphaseWell >( childName );
  }
  else if( childKey == keys::singlePhaseWell )
  {
    setCompositional( false );
    return &registerGroup< SinglePhaseWell >( childName );
  }
  else
  {
    PhysicsSolverBase::createChild( childKey, childName );
    return nullptr;
  }
}

void WellManager::expandObjectCatalogs()
{
  createChild( keys::compositionalMultiphaseWell, keys::compositionalMultiphaseWell );
  createChild( keys::singlePhaseWell, keys::singlePhaseWell );
}

void WellManager::registerDataOnMesh( Group & meshBodies )
{

  //std::string const & flowSolverName = getParent().getName();//getGroup< CompositionalMultiphaseBase >().getName();
  // CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      WellControls & well =  getWellControls( subRegion );
      //well.setFlowSolverName( flowSolver.getName() );
      well.setThermal( isThermal() );
      well.setUseMass( m_useMass );
      well.registerWellDataOnMesh( subRegion );
      m_numFluidPhases = well.numFluidPhases();
      m_numFluidComponents = well.numFluidComponents();

    } );
  } );
  // 1. Set key dimensions of the problem
  // Empty check needed to avoid errors when running in schema generation mode.

  // 1 pressure + NC compositions + 1 connectionRate + temp if thermal
  m_numDofPerWellElement = isThermal() ? m_numFluidComponents + 3 : m_numFluidComponents + 2;
  // 1 pressure + NC compositions + temp if thermal
  m_numDofPerResElement = isThermal() ? m_numFluidComponents + 2 : m_numFluidComponents + 1;

}

WellControls & WellManager::getWell( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}
WellControls & WellManager::getWell( std::string const & wellControlsName )
{
  return this->getGroup< WellControls >( wellControlsName );
}

WellControls const & WellManager::getWell( std::string const & wellControlsName ) const
{
  return this->getGroup< WellControls >( wellControlsName );
}
void WellManager::implicitStepSetup( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain )
{

  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementRegions< WellElementRegion >( regionNames,
                                                        [&]( localIndex const,
                                                             WellElementRegion & region )
    {
      // TODO expose a getWellControlsName() on WellElementRegion because Wrapper look-up is not useful here.
      WellControls & well = getWell(
        region.getReference< string >( WellElementRegion::viewKeyStruct::wellControlsString() ) );
      if( well.estimateSolution() )
      {
        well.setupWellDofs( domain, region, meshBodyName, mesh );
      }

    } )
    ;
  } );

  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls & well = getWell( subRegion );
      well.implicitStepSetup( time_n, dt, domain, meshBodyName, elemManager, subRegion );
    } );
  } );
}
real64
WellManager::setNextDt( real64 const & currentTime, const real64 & currentDt, geos::DomainPartition & domain )
{

  real64 nextDt = PhysicsSolverBase::setNextDt( currentTime, currentDt, domain );

  if( m_timeStepFromTables )
  {
    real64 nextDt_orig = nextDt;
    real64 nextDtLocal = nextDt;
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            WellElementSubRegion & subRegion )
      {

        WellControls & wellControls = getWellControls( subRegion );
        real64 nextDtWell = wellControls.setNextDt( currentTime, nextDt, subRegion );
        if( nextDtWell < nextDtLocal )
        {
          nextDtLocal = nextDtWell;
        }
      } );
    } );
    // get the minimum across all ranks
    nextDt = MpiWrapper::min< real64 >( nextDtLocal );
    if( getLogLevel() > 0 && nextDt < nextDt_orig )
      GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on tables coordinates = {}", getName(), nextDt ));
  }

  return nextDt;
}
localIndex WellManager::numDofPerWellElement() const
{
  return m_numDofPerWellElement;
}

localIndex WellManager::numDofPerResElement() const
{
  return m_numDofPerResElement;
}
integer WellManager::isThermal() const
{
  return m_isThermal;
}

string WellManager::wellElementDofName() const
{
  return viewKeyStruct::dofFieldString();
}

string WellManager::resElementDofName() const
{
  if( isCompositional() )
    return CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString();
  else
    return SinglePhaseBase::viewKeyStruct::elemDofFieldString();
}

localIndex WellManager::numFluidComponents() const
{
  return m_numFluidComponents;
}

localIndex WellManager::numFluidPhases() const
{
  return m_numFluidPhases;
}
WellControls & WellManager::getWellControls( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

WellControls const & WellManager::getWellControls( WellElementSubRegion const & subRegion ) const
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

CompositionalMultiphaseWell & WellManager::getCompositionalMultiphaseWell( WellElementSubRegion const & subRegion )
{
  return this->getGroup< CompositionalMultiphaseWell >( subRegion.getWellControlsName());
}

CompositionalMultiphaseWell const & WellManager::getCompositionalMultiphaseWell( WellElementSubRegion const & subRegion ) const
{
  return this->getGroup< CompositionalMultiphaseWell >( subRegion.getWellControlsName());
}
void WellManager::initializePostSubGroups()
{
  GEOS_MARK_FUNCTION;
  // Validate constitutive models
  if( isCompositional() )
  {
    DomainPartition & domain = ProblemRepository::getManager< DomainPartition >( *this );
    constitutive::ConstitutiveManager const & cm = domain.getConstitutiveManager();
    CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
    string const referenceFluidName = flowSolver.referenceFluidModelName();
    constitutive::MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< constitutive::MultiFluidBase >( referenceFluidName );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel const & mesh,
                                                                 string_array const & regionNames )
    {

      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            WellElementSubRegion const & subRegion )
      {
        string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseWell::viewKeyStruct::fluidNamesString() );
        constitutive::MultiFluidBase const & fluid = getConstitutiveModel< constitutive::MultiFluidBase >( subRegion, fluidName );
        CompositionalMultiphaseWell * wellControls = dynamic_cast< CompositionalMultiphaseWell * >(&getWellControls ( subRegion ));
        wellControls->validateFluidModel( fluid, referenceFluid );
      } );

    } );
  }
  else
  {
    // Single phase validation can be added here in the future
  }

}

void WellManager::setupDofs( DomainPartition const & domain,
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

void WellManager::assembleSystem( real64 const time,
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
  IterationsStatistics const & iterationsStatistics =  getIterationStats();
  //selectWellConstraint( time, dt, solverStatistics.m_numNewtonIterations, domain );

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
      wellControls.selectWellConstraint( time,
                                         dt,
                                         iterationsStatistics.getNumTimeSteps(),
                                         nonlinearParams.m_numNewtonIterations,
                                         domain,
                                         meshBodyName,
                                         meshLevel,
                                         elementRegionManager,
                                         subRegion,
                                         dofManager );

      // assemble the accumulation term in the mass balance equations
      wellControls.assembleWellAccumulationTerms( time, dt, subRegion, dofManager, localMatrix, localRhs );
      if( wellControls.isWellOpen() )
      {
        // assemble the pressure relations between well elements
        wellControls.assembleWellPressureRelations( time, dt, subRegion, dofManager, localMatrix, localRhs );
        // assemble well constraint terms
        wellControls.assembleWellConstraintTerms( time, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
        // compute the perforation rates (later assembled by the coupled solver)
        wellControls.computeWellPerforationRates( time, dt, elementRegionManager, subRegion );
        // assemble the flux terms in the mass balance equations
        wellControls.assembleWellFluxTerms( time, dt, subRegion, dofManager, localMatrix, localRhs );
      }
    } );
  } );

}


void WellManager::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls & wellControls = getWellControls( subRegion );
      wellControls.resetStateToBeginningOfStep( domain, meshBodyName, elemManager, subRegion );


    } );
  } );
}

void WellManager::implicitStepComplete( real64 const & time,
                                        real64 const & dt,
                                        DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls & wellControls = getWellControls( subRegion );
      wellControls.implicitStepComplete( time, dt, subRegion );
    } );
  } );
}

void WellManager::postRestartInitialization()
{}
void WellManager::initializePostInitialConditionsPreSubGroups()
{
  PhysicsSolverBase::initializePostInitialConditionsPreSubGroups();
  DomainPartition & domain = ProblemRepository::getManager< DomainPartition >( *this );
  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {

    // loop over the wells
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      // reconstruct local connectivity needed for flux calculations
      subRegion.reconstructLocalConnectivity();
      WellControls & wellControls = getWellControls( subRegion );
      wellControls.initializeWellPostInitialConditionsPreSubGroups( subRegion );


    } );
  } );
}
void WellManager::setKeepVariablesConstantDuringInitStep( bool const keepVariablesConstantDuringInitStep )
{
  DomainPartition & domain = ProblemRepository::getManager< DomainPartition >( *this );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    // loop over the wells
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )

    {
      WellControls & wellControls = getWellControls( subRegion );
      wellControls.setKeepVariablesConstantDuringInitStep( keepVariablesConstantDuringInitStep );

    } );
  } );
}
void WellManager::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  real64 maxPhaseVolFrac = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion & subRegion )
    {
      WellControls & wellControls = getWellControls( subRegion );
      if( wellControls.getWellState())
      {

        real64 const maxRegionPhaseVolFrac = wellControls.updateWellState( domain.getMeshBody( meshBodyName ), elemManager, subRegion );

        maxPhaseVolFrac = LvArray::math::max( maxRegionPhaseVolFrac, maxPhaseVolFrac );
      }
    } );
  } );
  maxPhaseVolFrac = MpiWrapper::max( maxPhaseVolFrac );

  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                         GEOS_FMT( "        {}: Max well phase volume fraction change = {}",
                                   getName(), fmt::format( "{:.{}f}", maxPhaseVolFrac, 4 ) ) );

}

real64
WellManager::calculateResidualNorm( real64 const & time_n,
                                    real64 const & dt,
                                    DomainPartition const & domain,
                                    DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer numNorm = 1;       // mass balance
  array1d< real64 > localResidualNorm, wellResidualNorm;
  array1d< real64 > localResidualNormalizer;

  if( isThermal() )
  {
    numNorm = 2;       // mass balance and energy balance
  }
  localResidualNorm.resize( numNorm );

  localResidualNormalizer.resize( numNorm );

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {


    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      WellControls & wellControls = getWellControls( subRegion );

      // step 1: compute the norm in the subRegion
      if( wellControls.isWellOpen( ) )
      {
        wellResidualNorm = wellControls.calculateLocalWellResidualNorm( time_n,
                                                                        dt,
                                                                        m_nonlinearSolverParameters,
                                                                        subRegion,
                                                                        dofManager,
                                                                        localRhs );
        for( integer i=0; i<numNorm; i++ )
        {
          if( wellResidualNorm[i] > localResidualNorm[i] )
          {
            localResidualNorm[i] = wellResidualNorm[i];
          }
        }
      }
      else
      {
        for( integer i=0; i<numNorm; i++ )
        {
          localResidualNorm[i] = 0.0;
        }

      }
    } );
  } );

  // step 3: second reduction across MPI ranks
  real64 resNorm=localResidualNorm[0];
  if( isThermal() )
  {
    real64 globalResidualNorm[2]{};
    globalResidualNorm[0] = MpiWrapper::max( localResidualNorm[0] );
    globalResidualNorm[1] = MpiWrapper::max( localResidualNorm[1] );
    resNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                                                coupledSolverAttributePrefix(), globalResidualNorm[0], globalResidualNorm[1] ));

    getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), globalResidualNorm[0] );
    getConvergenceStats().setResidualValue( "Renergy", globalResidualNorm[1] );
  }
  else
  {
    resNorm = MpiWrapper::max( resNorm );

    GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )",
                                                                coupledSolverAttributePrefix(), resNorm ));
    getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), resNorm );
  }
  return resNorm;
}
real64
WellManager::scalingForSystemSolution( DomainPartition & domain,
                                       DofManager const & dofManager,
                                       arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  real64 scalingFactor = 1.0;
  real64 localScalingFactor = 1.0;
  if( isCompositional() )
  {


    real64 maxDeltaPres = 0.0, maxDeltaCompDens = 0.0, maxDeltaTemp = 0.0;
    real64 minPresScalingFactor = 1.0, minCompDensScalingFactor = 1.0, minTempScalingFactor = 1.0;
    real64 localMaxDeltaPres = 0.0, localMaxDeltaCompDens = 0.0, localMaxDeltaTemp = 0.0;
    real64 localMinPresScalingFactor = 1.0, localMinCompDensScalingFactor = 1.0, localMinTempScalingFactor = 1.0;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {

      ElementRegionManager & elemManager = mesh.getElemManager();
      elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     WellElementSubRegion & subRegion )

      {
        CompositionalMultiphaseWell * wellControls = dynamic_cast< CompositionalMultiphaseWell * >(&getWellControls ( subRegion ));
        localScalingFactor = wellControls->scalingForLocalSystemSolution( subRegion,
                                                                          dofManager,
                                                                          localMaxDeltaPres,
                                                                          localMaxDeltaCompDens,
                                                                          localMaxDeltaTemp,
                                                                          localMinPresScalingFactor,
                                                                          localMinCompDensScalingFactor,
                                                                          localMinTempScalingFactor,
                                                                          localSolution );
        maxDeltaPres = LvArray::math::max( localMaxDeltaPres, maxDeltaPres );
        maxDeltaCompDens = LvArray::math::max( localMaxDeltaCompDens, maxDeltaCompDens );
        maxDeltaTemp = LvArray::math::max( localMaxDeltaTemp, maxDeltaTemp );
        minPresScalingFactor = LvArray::math::min( localMinPresScalingFactor, minPresScalingFactor );
        minCompDensScalingFactor = LvArray::math::min( localMinCompDensScalingFactor, minCompDensScalingFactor );
        minTempScalingFactor = LvArray::math::min( localMinTempScalingFactor, minTempScalingFactor );
        scalingFactor = LvArray::math::min( localScalingFactor, scalingFactor );

      } );
    } );

    scalingFactor = MpiWrapper::min( scalingFactor );
    maxDeltaPres  = MpiWrapper::max( maxDeltaPres );
    maxDeltaCompDens = MpiWrapper::max( maxDeltaCompDens );
    minPresScalingFactor = MpiWrapper::min( minPresScalingFactor );
    minCompDensScalingFactor = MpiWrapper::min( minCompDensScalingFactor );

    string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Max well pressure change: {} Pa (before scaling)",
                                     getName(), GEOS_FMT( "{:.{}f}", maxDeltaPres, 3 ) ) );
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Max well component density change: {} {} (before scaling)",
                                     getName(), GEOS_FMT( "{:.{}f}", maxDeltaCompDens, 3 ), massUnit ) );

    if( m_isThermal )
    {
      maxDeltaTemp = MpiWrapper::max( maxDeltaTemp );
      minTempScalingFactor = MpiWrapper::min( minTempScalingFactor );
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Max well temperature change: {} K (before scaling)",
                                       getName(), GEOS_FMT( "{:.{}f}", maxDeltaTemp, 3 ) ) );
    }


    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Min well pressure scaling factor: {}",
                                     getName(), minPresScalingFactor ) );
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Min well component density scaling factor: {}",
                                     getName(), minCompDensScalingFactor ) );
    if( m_isThermal )
    {
      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Min well temperature scaling factor: {}",
                                       getName(), minTempScalingFactor ) );
    }

  }
  else
  {
    // Single phase well scaling- not implemented yet
    scalingFactor=1.0;
  }
  return LvArray::math::max( scalingFactor, m_minScalingFactor );

}

bool
WellManager::checkSystemSolution( DomainPartition & domain,
                                  DofManager const & dofManager,
                                  arrayView1d< real64 const > const & localSolution,
                                  real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  integer globalCheck = 1;

  real64 minPressure = 0.0, minDensity = 0.0, minTotalDensity = 0.0;

  bool const solutionLogActive = isLogLevelActive< logInfo::Solution >( getLogLevel() );
  bool const solutionDetailsLogActive = isLogLevelActive< logInfo::SolutionDetails >( getLogLevel() );
  ElementsReporterBuffer rankNegPressureIds{ solutionLogActive, solutionDetailsLogActive ? 16 : 0 };
  ElementsReporterBuffer rankNegDensityIds{ solutionLogActive, solutionDetailsLogActive ? 16 : 0 };
  // output only total density sum, not cell details
  ElementsReporterBuffer rankTotalNegDensityIds{ solutionLogActive, 0 };

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )

    {
      WellControls & wellControls = getWellControls( subRegion );
      integer localCheck = wellControls.checkWellSystemSolution( subRegion,
                                                                 dofManager,
                                                                 localSolution,
                                                                 scalingFactor,
                                                                 minPressure,
                                                                 minDensity,
                                                                 minTotalDensity,
                                                                 rankNegPressureIds,
                                                                 rankNegDensityIds,
                                                                 rankTotalNegDensityIds );
      globalCheck = std::min( localCheck, globalCheck );
    } );
  } );
  globalCheck  = MpiWrapper::min( globalCheck );
  minPressure  = MpiWrapper::min( minPressure );
  minDensity = MpiWrapper::min( minDensity );
  minTotalDensity = MpiWrapper::min( minTotalDensity );

  units::Unit const massUnit = m_useMass ? units::Unit::Density : units::Unit::MolarDensity;
  rankNegPressureIds.createOutput()
    .outputTooLowValues( GEOS_FMT( "        {}: ", getName() ),
                         "negative pressure", minPressure, units::Unit::Pressure );
  rankNegDensityIds.createOutput()
    .outputTooLowValues( GEOS_FMT( "        {}: ", getName() ),
                         "negative component density", minDensity, massUnit );
  rankTotalNegDensityIds.createOutput()
    .outputTooLowValues( GEOS_FMT( "        {}: ", getName() ),
                         "negative components total density", minTotalDensity, massUnit );

  return globalCheck;
}

void
WellManager::applySystemSolution( DofManager const & dofManager,
                                  arrayView1d< real64 const > const & localSolution,
                                  real64 const scalingFactor,
                                  real64 const dt,
                                  DomainPartition & domain )
{


  DofManager::CompMask pressureMask( m_numDofPerWellElement, 0, 1 );

  DofManager::CompMask connRateMask( m_numDofPerWellElement, numFluidComponents()+1, numFluidComponents()+2 );
  GEOS_UNUSED_VAR( dt );
  // update all the fields using the global damping coefficients
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               well::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               well::connectionRate::key(),
                               scalingFactor,
                               connRateMask );
  if( isCompositional())
  {
    DofManager::CompMask componentMask( m_numDofPerWellElement, 1, numFluidComponents()+1 );
    dofManager.addVectorToField( localSolution,
                                 wellElementDofName(),
                                 well::globalCompDensity::key(),
                                 scalingFactor,
                                 componentMask );

  }
  if( isThermal() )
  {
    DofManager::CompMask temperatureMask( m_numDofPerWellElement, numFluidComponents()+2, numFluidComponents()+3 );

    dofManager.addVectorToField( localSolution,
                                 wellElementDofName(),
                                 well::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );

  }
  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( isCompositional() && m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    stdVector< string > propNames;
    propNames.emplace_back( well::pressure::key() );

    propNames.emplace_back( well::connectionRate::key() );
    if( isCompositional())
    {
      propNames.emplace_back( well::globalCompDensity::key() );
    }
    if( isThermal() )
    {
      propNames.emplace_back( well::temperature::key() );
    }
    // synchronize
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( propNames,
                                     regionNames );


    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );

}

void WellManager::chopNegativeDensities( DomainPartition & domain )
{
  integer const numComp = m_numFluidComponents;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getField< well::globalCompDensity >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
      {
        if( wellElemGhostRank[iwelem] < 0 )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            // we allowed for some densities to be slightly negative in CheckSystemSolution
            // if the new density is negative, chop back to zero
            if( wellElemCompDens[iwelem][ic] < 0 )
            {
              wellElemCompDens[iwelem][ic] = 0;
            }
          }
        }
      } );
    } );

  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, WellManager, string const &, Group * const )
}     // namespace geos
