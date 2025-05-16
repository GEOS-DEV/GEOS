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

#include "physicsSolvers/LogLevelsInfo.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/ThermalCompositionalMultiphaseWellKernels.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

WellSolverBase::WellSolverBase( string const & name,
                                Group * const parent )
  : PhysicsSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_isThermal( 0 ),
  m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), name + "_rates" ) ),
  m_keepVariablesConstantDuringInitStep( 0 ),
  m_estimateSolution( 0 )
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
    setDescription( "When set to 1, write the rates into a CSV file" );

  this->registerWrapper( viewKeyStruct::timeStepFromTablesFlagString(), &m_timeStepFromTables ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Choose time step to honor rates/bhp tables time intervals" );

  this->registerWrapper( viewKeyStruct::estimateWellSolutionString(), &m_estimateSolution ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to esitmate well solution prior to coupled reservoir and well solve." );

}

Group *WellSolverBase::createChild( string const & childKey, string const & childName )
{
  const auto childTypes = { keys::wellControls };
  GEOS_ERROR_IF( childKey != keys::wellControls,
                 PhysicsSolverBase::CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ) );
  return &registerGroup< WellControls >( childName );
}

void WellSolverBase::expandObjectCatalogs()
{
  createChild( keys::wellControls, keys::wellControls );
}

WellSolverBase::~WellSolverBase() = default;

void WellSolverBase::postInputInitialization()
{
  PhysicsSolverBase::postInputInitialization();

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
      subRegion.registerField< fields::well::pressure >( getName() );
      subRegion.registerField< fields::well::pressure_n >( getName() );

      subRegion.registerField< fields::well::temperature >( getName() );
      if( isThermal() )
      {
        subRegion.registerField< fields::well::temperature_n >( getName() );
      }

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

void WellSolverBase::initializePostSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      validateWellConstraints( 0, 0, subRegion );
    } );
  } );
}

void WellSolverBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  PhysicsSolverBase::setConstitutiveNamesCallSuper( subRegion );
  subRegion.registerWrapper< string >( viewKeyStruct::fluidNamesString() ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setSizedFromParent( 0 );
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

void WellSolverBase::implicitStepSetup( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition & domain )
{
  // Initialize the primary and secondary variables for the first time step

  initializeWells( domain, time_n );

}

void WellSolverBase::setupWellDofs( DomainPartition & domain )
{
  if( m_estimatorDoFManager.empty() )
  {

    map< std::pair< string, string >, string_array > meshTargets;
    string_array regions;
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                                 MeshLevel & meshLevel,
                                                                 string_array const & regionNames )
    {
      ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
      elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                   [&]( localIndex const,
                                                                        WellElementRegion & region )
      {
        meshTargets.clear();
        regions.clear();
        regions.emplace_back( region.getName() );
        auto const key = std::make_pair( meshBodyName, meshLevel.getName());
        meshTargets[key] = std::move( regions );

        DofManager regionDoFManager( region.getName());
        regionDoFManager.setDomain( domain );
        regionDoFManager.addField( wellElementDofName(),
                                   FieldLocation::Elem,
                                   numDofPerWellElement(),
                                   meshTargets );

        regionDoFManager.addCoupling( wellElementDofName(),
                                      wellElementDofName(),
                                      DofManager::Connector::Node );

        regionDoFManager.reorderByRank();
        m_estimatorDoFManager.emplace( region.getName(), std::move( regionDoFManager ));
      } );
    } );
  }
}
void WellSolverBase::estimateWellSolution( real64 const & time_n,
                                           real64 const & dt,
                                           const integer cycleNumber,
                                           DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  if( !estimateSolution() )
    return;
  setupWellDofs( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                               MeshLevel & meshLevel,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion & region )
    {
      auto it = m_estimatorDoFManager.find( region.getName());
      if( it == m_estimatorDoFManager.end())
      {
        throw std::runtime_error( "DofManager for region " + region.getName() + " not found." );
      }
      DofManager & dofManager = it->second;
      WellElementSubRegion & subRegion = region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                                           .getGroup< WellElementSubRegion >( region.getSubRegionName() );
// Only build the sparsity pattern if the mesh has changed
      Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

      if( meshModificationTimestamp > getSystemSetupTimestamp() )
      {
        setupWellSystem( domain, dofManager, m_localMatrix, m_rhs, m_solution );
        //setSystemSetupTimestamp( meshModificationTimestamp );

        //std::ostringstream oss;
        //m_dofManager.printFieldInfo( oss );
        //GEOS_LOG_LEVEL( logInfo::Fields, oss.str())
      }

      //implicitStepSetup( time_n, dt, domain );

// currently the only method is implicit time integration
      //real64 const dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

// final step for completion of timestep. typically secondary variable updates and cleanup.
      //implicitStepComplete( time_n, dt_return, domain );

      solveNonlinearSystem( time_n,
                            dt,
                            cycleNumber,
                            domain,
                            elementRegionManager,
                            subRegion,
                            dofManager );

    } );
  } );
}

void WellSolverBase::setupWellSystem( DomainPartition & domain,
                                      DofManager & dofManager,
                                      CRSMatrix< real64, globalIndex > & localMatrix,
                                      ParallelVector & rhs,
                                      ParallelVector & solution,
                                      bool const setSparsity )
{
  GEOS_MARK_FUNCTION;

  setupWellDofs( domain );

  if( setSparsity )
  {
    SparsityPattern< globalIndex > pattern;
    dofManager.setSparsityPattern( pattern );
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  }
  localMatrix.setName( this->getName() + "/matrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );
}
void WellSolverBase::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    { updateSubRegionState( subRegion ); } );
  } );
}

void WellSolverBase::assembleWellSystem( real64 const time_n,
                                         real64 const dt,
                                         ElementRegionManager const & elementRegionManager,
                                         WellElementSubRegion & subRegion,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs )
{
  assembleWellAccumulationTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  assembleWellPressureRelations( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  computeWellPerforationRates( time_n, dt, elementRegionManager, subRegion );
  assembleWellFluxTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
}

void WellSolverBase::assembleSystem( real64 const time,
                                     real64 const dt,
                                     DomainPartition & domain,
                                     DofManager const & dofManager,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
{
  string const wellDofKey = dofManager.getKey( wellElementDofName());


  // assemble the accumulation term in the mass balance equations
  assembleAccumulationTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the pressure relations between well elements
  assemblePressureRelations( time, dt, domain, dofManager, localMatrix, localRhs );
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

real64 WellSolverBase::setNextDt( real64 const & currentTime, const real64 & currentDt, geos::DomainPartition & domain )
{
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
        WellControls & wellControls = getWellControls( subRegion );
        real64 const nextDt_orig = nextDt;
        wellControls.setNextDtFromTables( currentTime, nextDt );
        if( m_nonlinearSolverParameters.getLogLevel() > 0 && nextDt < nextDt_orig )
          GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on tables coordinates = {}", getName(), nextDt ));
      } );
    } );
  }

  return nextDt;
}

bool WellSolverBase::solveNonlinearSystem( real64 const & time_n,
                                           real64 const & stepDt,
                                           integer const cycleNumber,
                                           DomainPartition & domain,
                                           ElementRegionManager & elemManager,
                                           WellElementSubRegion & subRegion,
                                           DofManager const & dofManager )
{
  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer & dtAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;
  integer & configurationLoopIter = m_nonlinearSolverParameters.m_numConfigurationAttempts;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

// keep residual from previous iteration in case we need to do a line search
  real64 lastResidual = 1e99;
  integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
  real64 scaleFactor = 1.0;

  bool isNewtonConverged = false;

  for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
  {

    GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                           GEOS_FMT( "    Attempt: {:2}, ConfigurationIter: {:2}, NewtonIter: {:2}", dtAttempt, configurationLoopIter, newtonIter ));

    {
      Timer timer( m_timers["assemble"] );

// We sync the nonlinear convergence history. The coupled solver parameters are the one being
// used. We want to propagate the info to subsolvers. It can be important for solvers that
// have special treatment for specific iterations.
      synchronizeNonlinearSolverParameters();

// zero out matrix/rhs before assembly
      m_localMatrix.zero();
      m_rhs.zero();

      arrayView1d< real64 > const localRhs = m_rhs.open();

// call assemble to fill the matrix and the rhs
      assembleWellSystem( time_n,
                          stepDt,
                          elemManager,
                          subRegion,
                          dofManager,
                          m_localMatrix.toViewConstSizes(),
                          localRhs );

// apply boundary conditions to system
      applyWellBoundaryConditions( time_n,
                                   stepDt,
                                   elemManager,
                                   subRegion,
                                   dofManager,
                                   localRhs,
                                   m_localMatrix.toViewConstSizes() );

      m_rhs.close();

      if( m_assemblyCallback )
      {
// Make a copy of LA objects and ship off to the callback
        array1d< real64 > localRhsCopy( m_rhs.localSize() );
        localRhsCopy.setValues< parallelDevicePolicy<> >( m_rhs.values() );
        m_assemblyCallback( m_localMatrix, std::move( localRhsCopy ) );
      }
    }

    real64 residualNorm = 0;
    {
      Timer timer( m_timers["convergence check"] );

// get residual norm
      residualNorm = calculateWellResidualNorm( time_n, stepDt, subRegion, dofManager, m_rhs.values() );
      GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                             GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );
    }

// if the residual norm is less than the Newton tolerance we denote that we have
// converged and break from the Newton loop immediately.
    if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
    {
      isNewtonConverged = true;
      break;
    }

// if the residual norm is above the max allowed residual norm, we break from
// the Newton loop to avoid crashes due to Newton divergence
    if( residualNorm > m_nonlinearSolverParameters.m_maxAllowedResidualNorm )
    {
      string const maxAllowedResidualNormString = NonlinearSolverParameters::viewKeysStruct::maxAllowedResidualNormString();
      GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                             GEOS_FMT( "    The residual norm is above the {} of {}. Newton loop terminated.",
                                       maxAllowedResidualNormString,
                                       m_nonlinearSolverParameters.m_maxAllowedResidualNorm )  );
      isNewtonConverged = false;
      break;
    }

// do line search in case residual has increased
    if( m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None
        && residualNorm > lastResidual * m_nonlinearSolverParameters.m_lineSearchResidualFactor
        && newtonIter >= m_nonlinearSolverParameters.m_lineSearchStartingIteration )
    {
      bool lineSearchSuccess = false;
      if( m_nonlinearSolverParameters.m_lineSearchInterpType == NonlinearSolverParameters::LineSearchInterpolationType::Linear )
      {
        residualNorm = lastResidual;
        lineSearchSuccess = lineSearch( time_n,
                                        stepDt,
                                        cycleNumber,
                                        domain,
                                        dofManager,
                                        m_localMatrix.toViewConstSizes(),
                                        m_rhs,
                                        m_solution,
                                        scaleFactor,
                                        residualNorm );
      }
      else
      {
        lineSearchSuccess = lineSearchWithParabolicInterpolation( time_n,
                                                                  stepDt,
                                                                  cycleNumber,
                                                                  domain,
                                                                  dofManager,
                                                                  m_localMatrix.toViewConstSizes(),
                                                                  m_rhs,
                                                                  m_solution,
                                                                  scaleFactor,
                                                                  lastResidual,
                                                                  residualNorm );
      }

      if( !lineSearchSuccess )
      {
        if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Attempt )
        {
          GEOS_LOG_LEVEL_RANK_0( logInfo::LineSearch,
                                 "        Line search failed to produce reduced residual. Accepting iteration." );
        }
        else if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Require )
        {
// if line search failed, then break out of the main Newton loop. Timestep will be cut.
          GEOS_LOG_LEVEL_RANK_0( logInfo::LineSearch,
                                 "        Line search failed to produce reduced residual. Exiting Newton Loop." );
          break;
        }
      }
    }

    {
      Timer timer( m_timers["linear solver total"] );

// if using adaptive Krylov tolerance scheme, update tolerance.
      LinearSolverParameters::Krylov & krylovParams = m_linearSolverParameters.get().krylov;
      if( krylovParams.useAdaptiveTol )
      {
        krylovParams.relTolerance = newtonIter > 0 ? eisenstatWalker( residualNorm, lastResidual, krylovParams ) : krylovParams.weakestTol;
      }

// TODO: Trilinos currently requires this, re-evaluate after moving to Tpetra-based solvers
      if( m_precond )
      {
        m_precond->clear();
      }

      {
        Timer timer_setup( m_timers["linear solver create"] );

// Compose parallel LA matrix/rhs out of local LA matrix/rhs
//
        m_matrix.create( m_localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
      }

// Output the linear system matrix/rhs for debugging purposes
      debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs );

// Solve the linear system
      solveLinearSystem( dofManager, m_matrix, m_rhs, m_solution );

// Increment the solver statistics for reporting purposes
      m_solverStatistics.logNonlinearIteration( m_linearSolverResult.numIterations );

// Output the linear system solution for debugging purposes
      debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );
    }

    {
      Timer timer( m_timers["apply solution"] );

// Compute the scaling factor for the Newton update
      scaleFactor = scalingForWellSystemSolution( subRegion, dofManager, m_solution.values() );

      GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                             GEOS_FMT( "        {}: Global solution scaling factor = {}", getName(), scaleFactor ) );

      if( !checkWellSystemSolution( subRegion, dofManager, m_solution.values(), scaleFactor ) )
      {
// TODO try chopping (similar to line search)
        GEOS_LOG_RANK_0( GEOS_FMT( "    {}: Solution check failed. Newton loop terminated.", getName()) );
        break;
      }

// apply the system solution to the fields/variables
      applyWellSystemSolution( dofManager, m_solution.values(), scaleFactor, stepDt, domain );
    }

    {
      Timer timer( m_timers["update state"] );

// update non-primary variables (constitutive models)
      updateWellState( subRegion );
    }

    lastResidual = residualNorm;
  }

  return isNewtonConverged;
}


} // namespace geos
