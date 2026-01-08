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
 * @file CompositionalMultiphaseWell.cpp
 */

#include "CompositionalMultiphaseWell.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/ThermalCompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/PerforationFluxKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalSolutionScalingKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/SolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalSolutionCheckKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/GlobalComponentFractionKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/PhaseVolumeFractionKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/ThermalPhaseVolumeFractionKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/FluidUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseStatistics.hpp"

#include "physicsSolvers/fluidFlow/wells/WellBHPConstraints.hpp"

#include "physicsSolvers/fluidFlow/wells/WellVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPhaseVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellMassRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellLiquidRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/PipeFlowTableFunction.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellConstraintKernels.hpp"

// tjb wrong place
#include "physicsSolvers/multiphysics/CoupledReservoirAndWellKernels.hpp"
#include "common/MpiWrapper.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "events/EventManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;

CompositionalMultiphaseWell::CompositionalMultiphaseWell( const string & name,
                                                          Group * const parent )
  :
  WellSolverBase( name, parent ),
  m_useMass( false ),
  m_useTotalMassEquation( 1 ),
  m_maxCompFracChange( 1.0 ),
  m_maxRelativePresChange( 0.2 ),
  m_maxAbsolutePresChange( -1 ), // disabled by default
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 ),
  m_targetPhaseIndex( -1 ),
  m_wellDebugInit( false )
{
  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::useTotalMassEquationString(), &m_useTotalMassEquation ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use total mass equation" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (absolute) allowed change in the composition fraction of any component within a single time step, in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::maxRelativeCompDensChangeString(), &m_maxRelativeCompDensChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( LvArray::NumericLimits< real64 >::max/1.0e100 ). // disabled by default
    setDescription( "Maximum (relative) change in a component density between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString(), &m_maxRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (relative) change in pressure between two Newton iterations (recommended with rate control)" );

  this->registerWrapper( viewKeyStruct::maxAbsolutePresChangeString(), &m_maxAbsolutePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).       // disabled by default
    setDescription( "Maximum (absolute) pressure change in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::maxRelativeTempChangeString(), &m_maxRelativeTempChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (relative) change in temperature between two Newton iterations  " );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );
}

void CompositionalMultiphaseWell::postInputInitialization()
{
  WellSolverBase::postInputInitialization();

  GEOS_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) <<
                        ": The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOS_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) <<
                        ": The maximum absolute change in component fraction must larger or equal to 0.0" );

  GEOS_ERROR_IF_LE_MSG( m_maxRelativeCompDensChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::maxRelativeCompDensChangeString() ) <<
                        ": The maximum relative change in component density must be larger than 0.0" );
}

void CompositionalMultiphaseWell::registerDataOnMesh( Group & meshBodies )
{
  WellSolverBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_referenceFluidModelName.empty() )
      {
        m_referenceFluidModelName = getConstitutiveName< MultiFluidBase >( subRegion );
      }
    } );
  } );

  // 1. Set key dimensions of the problem
  // Empty check needed to avoid errors when running in schema generation mode.
  if( !m_referenceFluidModelName.empty() )
  {
    MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );
    m_numPhases = fluid0.numFluidPhases();
    m_numComponents = fluid0.numFluidComponents();
  }
  // 1 pressure + NC compositions + 1 connectionRate + temp if thermal
  m_numDofPerWellElement = isThermal() ? m_numComponents + 3 : m_numComponents + 2;
  // 1 pressure + NC compositions + temp if thermal
  m_numDofPerResElement = isThermal() ? m_numComponents + 2 : m_numComponents + 1;

  // loop over the wells
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      string const & fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.

      subRegion.registerField< well::globalCompDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< well::globalCompDensity_n >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      subRegion.registerField< well::mixtureConnectionRate >( getName() );
      subRegion.registerField< well::mixtureConnectionRate_n >( getName() );

      subRegion.registerField< well::globalCompFraction >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerField< well::dGlobalCompFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

      subRegion.registerField< well::phaseVolumeFraction >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< well::dPhaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents + 2 ); // dP, dT, dC

      subRegion.registerField< well::totalMassDensity >( getName() );
      subRegion.registerField< well::dTotalMassDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents +2 ); // dP, dT, dC

      subRegion.registerField< well::phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< well::pressureScalingFactor >( getName() );
      subRegion.registerField< well::temperatureScalingFactor >( getName() );
      subRegion.registerField< well::globalCompDensityScalingFactor >( getName() );

      PerforationData & perforationData = *subRegion.getPerforationData();
      perforationData.registerField< well::compPerforationRate >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      perforationData.registerField< well::dCompPerforationRate >( getName() ).
        reference().resizeDimension< 1, 2, 3 >( 2, m_numComponents, m_numComponents+ 2 );
      if( fluid.isThermal() )
      {
        perforationData.registerField< well::energyPerforationFlux >( getName() );
        perforationData.registerField< well::dEnergyPerforationFlux >( getName() ).
          reference().resizeDimension< 1, 2 >( 2, m_numComponents+2 );
      }

      WellControls & wellControls = getWellControls( subRegion );
      wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString() );
      if( wellControls.hasMinimumWHPConstraint() )
        wellControls.registerWrapper< real64 >( viewKeyStruct::currentWHPString() );
      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numPhases );

      wellControls.registerWrapper< real64 >( viewKeyStruct::massDensityString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::massDensityString() );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentMassRateString() );

      // write rates output header
      // the rank that owns the reference well element is responsible
      if( m_writeCSV > 0 && subRegion.isLocallyOwned() )
      {
        string const fileName = GEOS_FMT( "{}/{}.csv", m_ratesOutputDir, wellControls.getName() );
        string const massUnit = m_useMass ? "kg" : "mol";
        integer const useSurfaceConditions = wellControls.useSurfaceConditions();
        string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
        string const unitKey = useSurfaceConditions ? "s" : "r";
        integer const numPhase = m_numPhases;
        integer const numComp = m_numComponents;
        // format: time,bhp,total_rate,total_vol_rate,phase0_vol_rate,phase1_vol_rate,...
        makeDirsForPath( m_ratesOutputDir );
        GEOS_LOG( GEOS_FMT( "{}: Rates CSV generated at {}", getName(), fileName ) );
        std::ofstream outputFile( fileName );
        outputFile << "Time [s],dt[s],BHP [Pa]";
        if( wellControls.hasMinimumWHPConstraint() )
          outputFile << ",WHP [Pa]";
        outputFile << ",Total rate [" << massUnit << "/s],Total " << conditionKey << " volumetric rate [" << unitKey << "m3/s]";
        for( integer ip = 0; ip < numPhase; ++ip )
        {
          outputFile << ",Phase" << ip << " " << conditionKey << " volumetric rate [" << unitKey << "m3/s]";
        }
        for( integer ic = 0; ic < numComp; ++ic )
        {
          outputFile << ",Component" << ic << " rate [" << massUnit << "/s]";
        }
        outputFile << std::endl;
        outputFile.close();
      }
    } );
  } );

}

void CompositionalMultiphaseWell::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  setConstitutiveName< MultiFluidBase >( subRegion, viewKeyStruct::fluidNamesString(), "multiphase fluid" );
}

namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                        GEOS_FMT( "Mismatch in number of phases between constitutive models {} and {}",
                                  lhs.getDataContext(), rhs.getDataContext() ),
                        InputError );

  for( integer ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOS_THROW_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                          GEOS_FMT( "Mismatch in phase names between constitutive models {} and {}",
                                    lhs.getDataContext(), rhs.getDataContext() ),
                          InputError );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                        GEOS_FMT( "Mismatch in number of components between constitutive models {} and {}",
                                  lhs.getDataContext(), rhs.getDataContext() ),
                        InputError );

  for( integer ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOS_THROW_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                          GEOS_FMT( "Mismatch in component names between constitutive models {} and {}",
                                    lhs.getDataContext(), rhs.getDataContext() ),
                          InputError );
  }
}

}

/**
 * @brief Checks if the WellControls parameters are within the fluid tables ranges
 * @param fluid the fluid to check
 */
void CompositionalMultiphaseWell::validateWellControlsForFluid( WellControls const & wellControls,
                                                                MultiFluidBase const & fluid ) const
{
  if( wellControls.useSurfaceConditions() )
  {
    try
    {
      real64 const & surfaceTemp = wellControls.getSurfaceTemperature();
      real64 const & surfacePres = wellControls.getSurfacePressure();
      fluid.checkTablesParameters( surfacePres, surfaceTemp );
    } catch( SimulationError const & ex )
    {
      string const errorMsg = GEOS_FMT( "{}: wrong surface pressure / temperature.\n", getDataContext() );
      ErrorLogger::global().currentErrorMsg()
        .addToMsg( errorMsg )
        .addContextInfo( getDataContext().getContextInfo().setPriority( 1 ) );
      throw SimulationError( ex, errorMsg );
    }
  }
}

void CompositionalMultiphaseWell::validateConstitutiveModels( DomainPartition const & domain ) const
{
  GEOS_MARK_FUNCTION;

  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
  string const referenceFluidName = flowSolver.referenceFluidModelName();
  MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion const & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      compareMultiphaseModels( fluid, referenceFluid );
      compareMulticomponentModels( fluid, referenceFluid );

      WellControls const & wellControls = getWellControls( subRegion );
      validateWellControlsForFluid( wellControls, fluid );
    } );

  } );
}



void CompositionalMultiphaseWell::validateWellConstraints( real64 const & time_n,
                                                           real64 const & GEOS_UNUSED_PARAM( dt ),
                                                           WellElementSubRegion const & subRegion )
{
  GEOS_UNUSED_VAR( time_n );
  WellControls & wellControls = getWellControls( subRegion );
  if( !wellControls.useSurfaceConditions() )
  {
    bool const useSeg =wellControls.referenceReservoirRegion().empty();
    GEOS_WARNING_IF( useSeg,
                     "WellControls " <<WellControls::viewKeyStruct::referenceReservoirRegionString() <<
                     " not set and well constraint fluid property calculations will use top segement pressure and temp " );
    if( useSeg )
    {
      wellControls.setRegionAveragePressure( -1 );
      wellControls.setRegionAverageTemperature( -1 );
    }
    else
    {
      // Check if region name exists in list of Reservoir's target regions
      string const regionName = wellControls.referenceReservoirRegion();
      CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
      string_array const & targetRegionsNames = flowSolver.getTargetRegionNames();
      auto const pos = std::find( targetRegionsNames.begin(), targetRegionsNames.end(), regionName );
      GEOS_ERROR_IF( pos == targetRegionsNames.end(),
                     GEOS_FMT( "{}: Region {} is not a target of the reservoir solver and cannot be used for referenceReservoirRegion in WellControl {}.",
                               getDataContext(), regionName, wellControls.getName() ) );


    }
  }
  // Disable BHP constraint if well head pressure constraint is active
  if( false &&wellControls.hasMinimumWHPConstraint() )
  {
    wellControls.forSubGroups< MinimumBHPConstraint >( [&]( auto & constraint )
    {
      std::cout << "Disable BHP Constraint " << std::endl;
      constraint.setConstraintActive( false );
    } );
  }

  // Check constraint phase types
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  // tjb
  //wellControls.forSubGroups< InjectionConstraint<LiquidRateConstraint>,InjectionConstraint< PhaseVolumeRateConstraint >,
  // ProductionConstraint< PhaseVolumeRateConstraint >,ProductionConstraint<LiquidRateConstraint> >( [&]( auto & constraint )

  wellControls.forSubGroups< InjectionConstraint< PhaseVolumeRateConstraint >, ProductionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
  {
    constraint.validatePhaseType( fluid );
  } );



}

void CompositionalMultiphaseWell::initializePostSubGroups()
{
  WellSolverBase::initializePostSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  validateConstitutiveModels( domain );

}

void CompositionalMultiphaseWell::outputWellDebug( real64 const time,
                                                   real64 const dt,
                                                   integer num_timesteps,
                                                   integer num_timestep_cuts,
                                                   integer current_newton_iteration,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  // tjbreturn;
  GEOS_UNUSED_VAR( time );
  GEOS_UNUSED_VAR( dofManager );
  GEOS_UNUSED_VAR( localMatrix );
  GEOS_UNUSED_VAR( localRhs );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {

      if( m_writeSegDebug > 1 )
      {
        CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
        auto solver_names = getParent().getSubGroupsNames();
        //integer n = solver_names.size();
        // Bit of a hack, cases with > 3 solvers we need to find the base solver for wells
        // Assume that solver definition order follows coupledreswell, res, and then well
        //std::string coupled_solver_name = solver_names[n-3];

        //GeosxState & gs = getGlobalState();

        //CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > * solver =
        //  &(gs.getProblemManager().getPhysicsSolverManager().getGroup< geos::CompositionalMultiphaseReservoirAndWells<
        // geos::CompositionalMultiphaseBase > >( coupled_solver_name ));

        EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
        real64 const & ctime = event.getReference< real64 >( EventManager::viewKeyStruct::timeString() );
        //real64 const  dt = event.getReference< real64 >( EventManager::viewKeyStruct::dtString() );
        integer const & cycle = event.getReference< integer >( EventManager::viewKeyStruct::cycleString() );
        integer const & subevent = event.getReference< integer >( EventManager::viewKeyStruct::currentSubEventString() );


        // std::cout << "tjbtime1 " << ctime <<  " " << m_globalNumTimeSteps <<  " " << dt << " " << cycle << " " << subevent
        // << " "  << m_numTimeStepCuts << " " << m_currentNewtonIteration << std::endl;
        if( true ) // need to fix for restarts cycle >= m_writeSegDebug   )
        {
          //SolverStatistics & solver_stat = solver->getSolverStatistics();
          //integer num_timesteps = solver_stat.getReference< integer >( SolverStatistics::viewKeyStruct::numTimeStepsString());
          //integer current_newton_iteration = solver_stat.getReference< integer >(
          // SolverStatistics::viewKeyStruct::numCurrentNonlinearIterationsString());
          //integer num_timestep_cuts = solver_stat.getReference< integer >( SolverStatistics::viewKeyStruct::numTimeStepCutsString());
          //std::cout << "tjbtime2 " << ctime <<  " " << m_globalNumTimeSteps <<  " " << dt << " " << cycle << " " << subevent
          //<< " "  << m_numTimeStepCuts << " " << m_currentNewtonIteration << std::endl;
          string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
          fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
          WellControls & wellControls = getWellControls( subRegion );


          MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
          PerforationData & perforationData = *subRegion.getPerforationData();
          using CompFlowAccessors =
            StencilAccessors< fields::flow::pressure,
                              fields::flow::temperature,
                              fields::flow::phaseVolumeFraction,



                              fields::flow::dPhaseVolumeFraction,
                              fields::flow::globalCompDensity,
                              fields::flow::dGlobalCompFraction_dGlobalCompDensity >;
          string const flowSolverName = flowSolver.getName();
          CompFlowAccessors compFlowAccessors( mesh.getElemManager(), flowSolver.getName() );

          using MultiFluidAccessors =
            StencilMaterialAccessors< MultiFluidBase,
                                      fields::multifluid::phaseEnthalpy,
                                      fields::multifluid::phaseDensity,
                                      fields::multifluid::phaseViscosity,
                                      fields::multifluid::dPhaseDensity,
                                      fields::multifluid::phaseViscosity,
                                      fields::multifluid::phaseInternalEnergy,
                                      fields::multifluid::dPhaseViscosity,
                                      fields::multifluid::phaseCompFraction,
                                      fields::multifluid::dPhaseCompFraction >;
          MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), flowSolver.getName() );

          using RelPermAccessors =
            StencilMaterialAccessors< RelativePermeabilityBase,
                                      fields::relperm::phaseRelPerm,
                                      fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;

          RelPermAccessors relPermAccessors( mesh.getElemManager(), flowSolver.getName() );


          string const srn = subRegion.getName();

          std::vector< string > cp_der {"dP", "dT"};
          for( integer i=0; i<m_numComponents; i++ )
          {
            cp_der.push_back( "dRho"+std::to_string( i+1 ));
          }
          if( !m_wellPropWriter[srn].initialized() )
          {
            integer my_rank = MpiWrapper::commRank( MPI_COMM_GEOS );
            m_wellPropWriter[srn].initialize_perf( my_rank, m_ratesOutputDir, wellControls.getName(), perforationData );
            m_wellPropWriter[srn].initialize_seg( my_rank, m_ratesOutputDir, wellControls.getName(), fluid.phaseNames(), fluid.componentNames(), subRegion );
          }
          m_wellPropWriter[srn].registerSeg2dProp( {"X", "Y", "Z"}, subRegion.getElementCenter());
          m_wellPropWriter[srn].registerSegProp( "Pressure", subRegion.getField< fields::well::pressure >());
          if( isThermal() )
          {
            m_wellPropWriter[srn].registerSegProp( "Temperature", subRegion.getField< fields::well::temperature >());
          }
          m_wellPropWriter[srn].registerSegComponentProp( "ComponentDensity", subRegion.getField< fields::well::globalCompDensity >());

          m_wellPropWriter[srn].registerSegProp( "TotalRate", subRegion.getField< fields::well::mixtureConnectionRate >());

          m_wellPropWriter[srn].registerSegProp( "MassDensity", subRegion.getField< fields::well::totalMassDensity >());

          m_wellPropWriter[srn].registerSegComponentProp( "CompFraction", subRegion.getField< fields::well::globalCompFraction >());

          m_wellPropWriter[srn].registerSegPhasePropf( "PhaseDensity", fluid.phaseMassDensity());
          m_wellPropWriter[srn].registerSegPhasePropf( "PhaseViscosity", fluid.phaseViscosity());
          m_wellPropWriter[srn].registerSegPhasePropDerf( "dPhaseDensity", cp_der, fluid.dPhaseMassDensity());
          m_wellPropWriter[srn].registerSegPhaseProp( "PhaseVolumeFraction", subRegion.getField< fields::well::phaseVolumeFraction >());
          m_wellPropWriter[srn].registerSegPhasePropDer( "dPhaseVolume", cp_der, subRegion.getField< fields::well::dPhaseVolumeFraction >());
          if( isThermal() )
          {
            m_wellPropWriter[srn].registerSegPhasePropf( "InternalEnergy", fluid.phaseInternalEnergy());
            m_wellPropWriter[srn].registerSegPhasePropDerf( "dPhaseEnthalpy", cp_der, fluid.dPhaseEnthalpy());

            m_wellPropWriter[srn].registerSegPhasePropf( "PhaseEnthalpy", fluid.phaseEnthalpy());
            m_wellPropWriter[srn].registerSegPhasePropDerf( "dPhaseInternalEnergy", cp_der, fluid.dPhaseInternalEnergy());
          }
          m_wellPropWriter[srn].registerSegPhaseComponentPropf( "PhaseCompFrac", fluid.phaseCompFraction());

          // Perforation properties
          m_wellPropWriter[srn].registerPerf2dProp( {"X", "Y", "Z"}, perforationData.getLocation());
          m_wellPropWriter[srn].registerPerfResProp( "Pressure", compFlowAccessors.get( fields::flow::pressure {} ));
          if( isThermal() )
          {
            m_wellPropWriter[srn].registerPerfResProp( "Temperature", compFlowAccessors.get( fields::flow::temperature{} ));
            m_wellPropWriter[srn].registerPerfResPhasePropf( "PhaseEnthalpy", multiFluidAccessors.get( fields::multifluid::phaseEnthalpy{} ));
            m_wellPropWriter[srn].registerPerfResPhasePropf( "PhaseInternalEnergy", multiFluidAccessors.get( fields::multifluid::phaseInternalEnergy{} ));
          }
          m_wellPropWriter[srn].registerPerfComponentProp( "CompPerfRate", perforationData.getField< fields::well::compPerforationRate >());
          // m_wellPropWriter[srn].registerPerfResComponentProp( "ComponentDensity", compFlowAccessors.get(
          // fields::flow::globalCompDensity{} ));
          m_wellPropWriter[srn].registerPerfResPhaseComponentProp( "PhaseCompFrac", multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ));
          m_wellPropWriter[srn].registerPerfResPhaseProp( "PhaseVolFrac", compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ));
          m_wellPropWriter[srn].registerPerfResPhasePropf( "Viscosity", multiFluidAccessors.get( fields::multifluid::phaseViscosity{} ));
          m_wellPropWriter[srn].registerPerfResPhasePropf( "RelPerm", relPermAccessors.get( fields::relperm::phaseRelPerm{} ));

          m_wellPropWriter[srn].write( ctime, dt, cycle, subevent, num_timesteps,
                                       current_newton_iteration,
                                       num_timestep_cuts );
        }
      }
    } );
  } );

}

void CompositionalMultiphaseWell::outputSingleWellDebug( real64 const time,
                                                         real64 const dt,
                                                         integer num_timesteps,
                                                         integer current_newton_iteration,
                                                         integer num_timestep_cuts,
                                                         MeshLevel & mesh,
                                                         WellElementSubRegion & subRegion,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< const real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time );
  GEOS_UNUSED_VAR( dofManager );
  GEOS_UNUSED_VAR( localMatrix );
  GEOS_UNUSED_VAR( localRhs );


  if( m_writeSegDebug > 1 )
  {
    CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
    auto solver_names = getParent().getSubGroupsNames();
//integer n = solver_names.size();
// Bit of a hack, cases with > 3 solvers we need to find the base solver for wells
// Assume that solver definition order follows coupledreswell, res, and then well
//std::string coupled_solver_name = solver_names[n-3];

//GeosxState & gs = getGlobalState();

//CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > * solver =
//  &(gs.getProblemManager().getPhysicsSolverManager().getGroup< geos::CompositionalMultiphaseReservoirAndWells<
// geos::CompositionalMultiphaseBase > >( coupled_solver_name ));

    EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
    //  real64 const & ctime = event.getReference< real64 >( EventManager::viewKeyStruct::timeString() );
//real64 const  dt = event.getReference< real64 >( EventManager::viewKeyStruct::dtString() );
    integer const & cycle = event.getReference< integer >( EventManager::viewKeyStruct::cycleString() );
    integer const & subevent = event.getReference< integer >( EventManager::viewKeyStruct::currentSubEventString() );


// std::cout << "tjbtime1 " << ctime <<  " " << m_globalNumTimeSteps <<  " " << dt << " " << cycle << " " << subevent
// << " "  << m_numTimeStepCuts << " " << m_currentNewtonIteration << std::endl;
    if( true ) // need to fix for restarts cycle >= m_writeSegDebug   )
    {
//SolverStatistics & solver_stat = solver->getSolverStatistics();
//integer num_timesteps = solver_stat.getReference< integer >( SolverStatistics::viewKeyStruct::numTimeStepsString());
//integer current_newton_iteration = solver_stat.getReference< integer >(
// SolverStatistics::viewKeyStruct::numCurrentNonlinearIterationsString());
//integer num_timestep_cuts = solver_stat.getReference< integer >( SolverStatistics::viewKeyStruct::numTimeStepCutsString());
//std::cout << "tjbtime2 " << ctime <<  " " << m_globalNumTimeSteps <<  " " << dt << " " << cycle << " " << subevent
//<< " "  << m_numTimeStepCuts << " " << m_currentNewtonIteration << std::endl;
      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
      WellControls & wellControls = getWellControls( subRegion );


      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      PerforationData & perforationData = *subRegion.getPerforationData();
      using CompFlowAccessors =
        StencilAccessors< fields::flow::pressure,
                          fields::flow::temperature,
                          fields::flow::phaseVolumeFraction,



                          fields::flow::dPhaseVolumeFraction,
                          fields::flow::globalCompDensity,
                          fields::flow::dGlobalCompFraction_dGlobalCompDensity >;
      string const flowSolverName = flowSolver.getName();
      CompFlowAccessors compFlowAccessors( mesh.getElemManager(), flowSolver.getName() );

      using MultiFluidAccessors =
        StencilMaterialAccessors< MultiFluidBase,
                                  fields::multifluid::phaseEnthalpy,
                                  fields::multifluid::phaseDensity,
                                  fields::multifluid::phaseViscosity,
                                  fields::multifluid::dPhaseDensity,
                                  fields::multifluid::phaseViscosity,
                                  fields::multifluid::phaseInternalEnergy,
                                  fields::multifluid::dPhaseViscosity,
                                  fields::multifluid::phaseCompFraction,
                                  fields::multifluid::dPhaseCompFraction >;
      MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), flowSolver.getName() );

      using RelPermAccessors =
        StencilMaterialAccessors< RelativePermeabilityBase,
                                  fields::relperm::phaseRelPerm,
                                  fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;

      RelPermAccessors relPermAccessors( mesh.getElemManager(), flowSolver.getName() );


      string const srn = subRegion.getName();

      std::vector< string > cp_der {"dP", "dT"};
      for( integer i=0; i<m_numComponents; i++ )
      {
        cp_der.push_back( "dRho"+std::to_string( i+1 ));
      }
      if( !m_wellPropWriter[srn].initialized() )
      {
        integer my_rank = MpiWrapper::commRank( MPI_COMM_GEOS );
        m_wellPropWriter[srn].initialize_perf( my_rank, m_ratesOutputDir, wellControls.getName(), perforationData );
        m_wellPropWriter[srn].initialize_seg( my_rank, m_ratesOutputDir, wellControls.getName(), fluid.phaseNames(), fluid.componentNames(), subRegion );
        m_wellDebugInit=true;
      }
      m_wellPropWriter[srn].registerSeg2dProp( {"X", "Y", "Z"}, subRegion.getElementCenter());
      m_wellPropWriter[srn].registerSegProp( "Pressure", subRegion.getField< fields::well::pressure >());
      if( isThermal() )
      {
        m_wellPropWriter[srn].registerSegProp( "Temperature", subRegion.getField< fields::well::temperature >());
      }
      m_wellPropWriter[srn].registerSegComponentProp( "ComponentDensity", subRegion.getField< fields::well::globalCompDensity >());

      m_wellPropWriter[srn].registerSegProp( "TotalRate", subRegion.getField< fields::well::mixtureConnectionRate >());

      m_wellPropWriter[srn].registerSegProp( "MassDensity", subRegion.getField< fields::well::totalMassDensity >());

      m_wellPropWriter[srn].registerSegComponentProp( "CompFraction", subRegion.getField< fields::well::globalCompFraction >());

      m_wellPropWriter[srn].registerSegPhasePropf( "PhaseDensity", fluid.phaseMassDensity());
      m_wellPropWriter[srn].registerSegPhasePropf( "PhaseViscosity", fluid.phaseViscosity());
      m_wellPropWriter[srn].registerSegPhasePropDerf( "dPhaseDensity", cp_der, fluid.dPhaseMassDensity());
      m_wellPropWriter[srn].registerSegPhaseProp( "PhaseVolumeFraction", subRegion.getField< fields::well::phaseVolumeFraction >());
      m_wellPropWriter[srn].registerSegPhasePropDer( "dPhaseVolume", cp_der, subRegion.getField< fields::well::dPhaseVolumeFraction >());
      if( isThermal() )
      {
        m_wellPropWriter[srn].registerSegPhasePropf( "InternalEnergy", fluid.phaseInternalEnergy());
        m_wellPropWriter[srn].registerSegPhasePropDerf( "dPhaseEnthalpy", cp_der, fluid.dPhaseEnthalpy());

        m_wellPropWriter[srn].registerSegPhasePropf( "PhaseEnthalpy", fluid.phaseEnthalpy());
        m_wellPropWriter[srn].registerSegPhasePropDerf( "dPhaseInternalEnergy", cp_der, fluid.dPhaseInternalEnergy());
      }
      m_wellPropWriter[srn].registerSegPhaseComponentPropf( "PhaseCompFrac", fluid.phaseCompFraction());

// Perforation properties
      m_wellPropWriter[srn].registerPerf2dProp( {"X", "Y", "Z"}, perforationData.getLocation());
      m_wellPropWriter[srn].registerPerf1dProp( {"Trans"}, perforationData.getWellTransmissibility());
      m_wellPropWriter[srn].registerPerfResProp( "Pressure", compFlowAccessors.get( fields::flow::pressure {} ));
      if( isThermal() )
      {
        m_wellPropWriter[srn].registerPerfResProp( "Temperature", compFlowAccessors.get( fields::flow::temperature{} ));
        m_wellPropWriter[srn].registerPerfResPhasePropf( "PhaseEnthalpy", multiFluidAccessors.get( fields::multifluid::phaseEnthalpy{} ));
        m_wellPropWriter[srn].registerPerfResPhasePropf( "PhaseInternalEnergy", multiFluidAccessors.get( fields::multifluid::phaseInternalEnergy{} ));
      }
      m_wellPropWriter[srn].registerPerfComponentProp( "CompPerfRate", perforationData.getField< fields::well::compPerforationRate >());
      //m_wellPropWriter[srn].registerPerfResComponentProp( "ComponentDensity", compFlowAccessors.get( fields::flow::globalCompDensity{} ));
      m_wellPropWriter[srn].registerPerfResPhaseComponentProp( "PhaseCompFrac", multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ));
      m_wellPropWriter[srn].registerPerfResPhaseProp( "PhaseVolFrac", compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ));
      m_wellPropWriter[srn].registerPerfResPhasePropf( "Viscosity", multiFluidAccessors.get( fields::multifluid::phaseViscosity{} ));
      m_wellPropWriter[srn].registerPerfResPhasePropf( "RelPerm", relPermAccessors.get( fields::relperm::phaseRelPerm{} ));

      m_wellPropWriter[srn].write( my_ctime, dt, cycle, subevent, num_timesteps,
                                   current_newton_iteration,
                                   num_timestep_cuts );
    }

  }


}


void CompositionalMultiphaseWell::printSegRates( real64 const & time,
                                                 real64 const & dt,
                                                 integer num_timesteps,
                                                 integer num_timestep_cuts,
                                                 integer current_newton_iteration,
                                                 DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time );


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {

      if( m_writeSegDebug > 0 )
      {
        CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
        auto solver_names = getParent().getSubGroupsNames();
        //integer n = solver_names.size();
        // Bit of a hack, cases with > 3 solvers we need to find the base solver for wells
        // Assume that solver definition order follows coupledreswell, res, and then well
        //std::string coupled_solver_name = solver_names[n-3];

        //GeosxState & gs = getGlobalState();

        //CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > * solver =
        //  &(gs.getProblemManager().getPhysicsSolverManager().getGroup< geos::CompositionalMultiphaseReservoirAndWells<
        // geos::CompositionalMultiphaseBase > >( coupled_solver_name ));

        EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
        real64 const & ctime = event.getReference< real64 >( EventManager::viewKeyStruct::timeString() );
        //real64 const & dt = event.getReference< real64 >( EventManager::viewKeyStruct::dtString() );
        integer const & cycle = event.getReference< integer >( EventManager::viewKeyStruct::cycleString() );
        integer const & subevent = event.getReference< integer >( EventManager::viewKeyStruct::currentSubEventString() );


        // std::cout << "tjbtime1 " << ctime <<  " " << m_globalNumTimeSteps <<  " " << dt << " " << cycle << " " << subevent
        // << " "  << m_numTimeStepCuts << " " << m_currentNewtonIteration << std::endl;
        if( true ) // need to fix for restarts cycle >= m_writeSegDebug   )
        {
          //SolverStatistics & solver_stat = solver->getSolverStatistics();
          //integer num_timesteps = solver_stat.getReference< integer >( SolverStatistics::viewKeyStruct::numTimeStepsString());
          //integer current_newton_iteration = solver_stat.getReference< integer >(
          // SolverStatistics::viewKeyStruct::numCurrentNonlinearIterationsString());
          //integer num_timestep_cuts = solver_stat.getReference< integer >( SolverStatistics::viewKeyStruct::numTimeStepCutsString());
          //std::cout << "tjbtime2 " << ctime <<  " " << m_globalNumTimeSteps <<  " " << dt << " " << cycle << " " << subevent
          //<< " "  << m_numTimeStepCuts << " " << m_currentNewtonIteration << std::endl;
          string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
          fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
          WellControls & wellControls = getWellControls( subRegion );


          MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
          PerforationData & perforationData = *subRegion.getPerforationData();
          using CompFlowAccessors =
            StencilAccessors< fields::flow::pressure,
                              fields::flow::temperature,
                              fields::flow::phaseVolumeFraction,
                              fields::flow::dPhaseVolumeFraction,
                              fields::flow::globalCompDensity,
                              fields::flow::dGlobalCompFraction_dGlobalCompDensity >;
          string const flowSolverName = flowSolver.getName();
          CompFlowAccessors compFlowAccessors( mesh.getElemManager(), flowSolver.getName() );

          using MultiFluidAccessors =
            StencilMaterialAccessors< MultiFluidBase,
                                      fields::multifluid::phaseDensity,
                                      fields::multifluid::dPhaseDensity,
                                      fields::multifluid::phaseInternalEnergy,
                                      fields::multifluid::phaseEnthalpy,
                                      fields::multifluid::dPhaseDensity,
                                      fields::multifluid::phaseViscosity,
                                      fields::multifluid::dPhaseViscosity,
                                      fields::multifluid::phaseCompFraction,
                                      fields::multifluid::dPhaseCompFraction >;
          MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), flowSolver.getName() );

          using RelPermAccessors =
            StencilMaterialAccessors< RelativePermeabilityBase,
                                      fields::relperm::phaseRelPerm,
                                      fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;

          RelPermAccessors relPermAccessors( mesh.getElemManager(), flowSolver.getName() );


          string const srn = subRegion.getName();

          std::vector< string > cp_der {"dP", "dT"};
          for( integer i=0; i<m_numComponents; i++ )
          {
            cp_der.push_back( "dRho"+std::to_string( i+1 ));
          }
          if( num_timesteps == 0 && num_timestep_cuts == 0 )
          {
            integer my_rank = MpiWrapper::commRank( MPI_COMM_GEOS );
            m_wellPropWriter_eot[srn].initialize_perf( my_rank, m_ratesOutputDir, wellControls.getName()+"_eot", perforationData );
            m_wellPropWriter_eot[srn].initialize_seg( my_rank, m_ratesOutputDir, wellControls.getName()+"_eot", fluid.phaseNames(), fluid.componentNames(), subRegion );
          }
          m_wellPropWriter_eot[srn].registerSeg2dProp( {"X", "Y", "Z"}, subRegion.getElementCenter());
          m_wellPropWriter_eot[srn].registerSegProp( "Pressure", subRegion.getField< fields::well::pressure >());
          if( isThermal() )
          {
            m_wellPropWriter_eot[srn].registerSegProp( "Temperature", subRegion.getField< fields::well::temperature >());
          }
          m_wellPropWriter_eot[srn].registerSegComponentProp( "ComponentDensity", subRegion.getField< fields::well::globalCompDensity >());

          m_wellPropWriter_eot[srn].registerSegProp( "TotalRate", subRegion.getField< fields::well::mixtureConnectionRate >());

          m_wellPropWriter_eot[srn].registerSegProp( "MassDensity", subRegion.getField< fields::well::totalMassDensity >());

          m_wellPropWriter_eot[srn].registerSegComponentProp( "CompFraction", subRegion.getField< fields::well::globalCompFraction >());

          m_wellPropWriter_eot[srn].registerSegPhasePropf( "PhaseDensity", fluid.phaseMassDensity());
          m_wellPropWriter_eot[srn].registerSegPhasePropf( "PhaseViscosity", fluid.phaseViscosity());
          m_wellPropWriter_eot[srn].registerSegPhaseProp( "PhaseVolumeFraction", subRegion.getField< fields::well::phaseVolumeFraction >());
          if( isThermal() )
          {
            m_wellPropWriter_eot[srn].registerSegPhasePropf( "InternalEnergy", fluid.phaseInternalEnergy());
            m_wellPropWriter_eot[srn].registerSegPhasePropf( "PhaseEnthalpy", fluid.phaseEnthalpy());
          }
          m_wellPropWriter_eot[srn].registerSegPhaseComponentPropf( "PhaseCompFrac", fluid.phaseCompFraction());

          // Perforation properties
          m_wellPropWriter_eot[srn].registerPerf2dProp( {"X", "Y", "Z"}, perforationData.getLocation());


          m_wellPropWriter_eot[srn].registerPerfResProp( "Pressure", compFlowAccessors.get( fields::flow::pressure {} ));
          if( isThermal() )
          {
            m_wellPropWriter_eot[srn].registerPerfResProp( "Temperature", compFlowAccessors.get( fields::flow::temperature{} ));
            m_wellPropWriter_eot[srn].registerPerfResPhasePropf( "PhaseEnthalpy", multiFluidAccessors.get( fields::multifluid::phaseEnthalpy{} ));
            m_wellPropWriter_eot[srn].registerPerfResPhasePropf( "PhaseInternalEnergy", multiFluidAccessors.get( fields::multifluid::phaseInternalEnergy{} ));
          }
          m_wellPropWriter_eot[srn].registerPerfComponentProp( "CompPerfRate", perforationData.getField< fields::well::compPerforationRate >());
          //m_wellPropWriter_eot[srn].registerPerfResComponentProp( "ComponentDensity",
          // compFlowAccessors.get(fields::flow::globalCompDensity{} ));
          m_wellPropWriter_eot[srn].registerPerfResPhaseComponentProp( "PhaseCompFrac", multiFluidAccessors.get( fields::multifluid::phaseCompFraction{} ));
          m_wellPropWriter_eot[srn].registerPerfResPhaseProp( "PhaseVolFrac", compFlowAccessors.get( fields::flow::phaseVolumeFraction{} ));
          m_wellPropWriter_eot[srn].registerPerfResPhasePropf( "Viscosity", multiFluidAccessors.get( fields::multifluid::phaseViscosity{} ));
          m_wellPropWriter_eot[srn].registerPerfResPhasePropf( "RelPerm", relPermAccessors.get( fields::relperm::phaseRelPerm{} ));

          m_wellPropWriter_eot[srn].write( ctime, dt, cycle, subevent, num_timesteps,
                                           current_newton_iteration,
                                           num_timestep_cuts );
        }
      }
    } );
  } );

}
void CompositionalMultiphaseWell::initializePostInitialConditionsPreSubGroups()
{
  WellSolverBase::initializePostInitialConditionsPreSubGroups();
  createSeparator();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    // loop over the wells
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      WellControls & wellControls = getWellControls( subRegion );
      // setup internal constraints if needed
      if( wellControls.hasMinimumWHPConstraint())
      {
        wellControls.createMinBHPConstraintForWHP();
        wellControls.createMaxLiquidConstraintForWHP();
      }
    } );
  } );
}

void CompositionalMultiphaseWell::postRestartInitialization()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    // loop over the wells
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )

    {
      // setup fluid separator
      WellControls & wellControls = getWellControls( subRegion );
      constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();
      fluidSeparator.allocateConstitutiveData( wellControls, 1 );
      fluidSeparator.resize( 1 );
    } );
  } );
}

void CompositionalMultiphaseWell::createSeparator()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    // loop over the wells
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      fluid.setMassFlag( m_useMass );
      // setup fluid separator
      WellControls & wellControls = getWellControls( subRegion );
      string const fluidSeparatorName = wellControls.getName() + "Separator";
      std::unique_ptr< constitutive::ConstitutiveBase >  fluidSeparatorPtr  = fluid.deliverClone( fluidSeparatorName, &wellControls );
      fluidSeparatorPtr->allocateConstitutiveData( wellControls, 1 );
      fluidSeparatorPtr->resize( 1 );
      wellControls.setFluidSeparator( std::move( fluidSeparatorPtr ));
    } );
  } );
}
void CompositionalMultiphaseWell::updateGlobalComponentFraction( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  isothermalCompositionalMultiphaseBaseKernels::
    GlobalComponentFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               subRegion );

}

void CompositionalMultiphaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();
  //localIndex const iwelemRef = subRegion.size()-1;
  // subRegion data
  arrayView1d< real64 const > const & pres = subRegion.getField< well::pressure >();
  arrayView1d< real64 > const & totalMassDens = subRegion.getField< well::totalMassDensity >();
  arrayView1d< real64 const > const wellElemGravCoef = subRegion.getField< well::gravityCoefficient >();

  // control data
  WellControls & wellControls = getWellControls( subRegion );
  string const wellControlsName = wellControls.getName();
  real64 const & refGravCoef = wellControls.getReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [ pres,
                               totalMassDens,
                               wellElemGravCoef,
                               &currentBHP,
                               &iwelemRef,
                               &refGravCoef] ( localIndex const )
  {
    real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
    currentBHP = pres[iwelemRef] + totalMassDens[iwelemRef] * diffGravCoef;
    std::cout << "tjbwellbhp " << pres[iwelemRef] << " " << totalMassDens[iwelemRef] << " "
              << wellElemGravCoef[iwelemRef] << " " << refGravCoef << " " << diffGravCoef << " " << currentBHP << std::endl;
  } );


  GEOS_LOG_LEVEL_BY_RANK( logInfo::BoundaryConditions,
                          GEOS_FMT( "{}: BHP (at the specified reference elevation) = {} Pa",
                                    wellControlsName, currentBHP ) );

}

void CompositionalMultiphaseWell::updateVolRatesForConstraint( WellElementSubRegion const & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }


  integer const numPhase = m_numPhases;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  WellControls & wellControls = getWellControls( subRegion );

  // subRegion data

  arrayView1d< real64 const > const & connRate = subRegion.getField< well::mixtureConnectionRate >();

  // fluid data
  constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();

  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac = fluidSeparator.phaseFraction();
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & totalDens = fluidSeparator.totalDensity();
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens = fluidSeparator.phaseDensity();

  // control data

  string const wellControlsName = wellControls.getName();

  arrayView1d< real64 > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );

  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

  real64 & currentMassRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() );

  real64 & massDensity =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::massDensityString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [&numPhase,
                              connRate,
                              totalDens,
                              phaseDens,
                              phaseFrac,
                              &currentTotalVolRate,
                              currentPhaseVolRate,
                              &currentMassRate,
                              &iwelemRef,
                              &massDensity] ( localIndex const )
  {
    // Step 1: update the total volume rate

    real64 const currentTotalRate = connRate[iwelemRef];
    // Assumes useMass is true
    currentMassRate = currentTotalRate;
    // Step 1.1: compute the inverse of the total density and derivatives
    massDensity = totalDens[iwelemRef][0];
    real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];

    // Step 1.2: divide the total mass/molar rate by the total density to get the total volumetric rate
    currentTotalVolRate = currentTotalRate * totalDensInv;

    // Step 2: update the phase volume rate
    for( integer ip = 0; ip < numPhase; ++ip )
    {
      // Step 2.1: compute the inverse of the (phase density * phase fraction) and derivatives

      // skip the rest of this function if phase ip is absent
      bool const phaseExists = (phaseFrac[iwelemRef][0][ip] > 0);
      if( !phaseExists )
      {
        continue;
      }

      real64 const phaseDensInv =  1.0 / phaseDens[iwelemRef][0][ip];
      real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;

      // Step 2.2: divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
      currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
    }
  } );

}

void CompositionalMultiphaseWell::calculateReferenceElementRates( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  integer const numPhase = m_numPhases;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  WellControls & wellControls = getWellControls( subRegion );

  // subRegion data
  arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::mixtureConnectionRate >();

  // fluid data
  constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac = fluidSeparator.phaseFraction();
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & totalDens = fluidSeparator.totalDensity();
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens = fluidSeparator.phaseDensity();

  // control data
  string const wellControlsName = wellControls.getName();
  bool const logSurfaceCondition = isLogLevelActive< logInfo::BoundaryConditions >( wellControls.getLogLevel());
  string const massUnit = m_useMass ? "kg" : "mol";

  integer const useSurfaceConditions = wellControls.useSurfaceConditions();

  arrayView1d< real64 > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );

  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

  real64 & currentMassRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() );

  real64 & massDensity =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::massDensityString() );

  constitutive::constitutiveUpdatePassThru( fluidSeparator, [&] ( auto & castedFluidSeparator )
  {
    // typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    typename TYPEOFREF( castedFluidSeparator ) ::KernelWrapper fluidSeparatorWrapper = castedFluidSeparator.createKernelWrapper();

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [fluidSeparatorWrapper,
                                &numPhase,
                                connRate,
                                totalDens,
                                phaseDens,
                                phaseFrac,
                                logSurfaceCondition,
                                &useSurfaceConditions,
                                &currentTotalVolRate,
                                currentPhaseVolRate,
                                &currentMassRate,
                                &iwelemRef,
                                &wellControlsName,
                                &massUnit,
                                &massDensity] ( localIndex const )
    {
      GEOS_UNUSED_VAR( massUnit );


      // Step 2: update the total volume rate

      real64 const currentTotalRate = connRate[iwelemRef];
      // Assumes useMass is true
      currentMassRate = currentTotalRate;
      // Step 2.1: compute the inverse of the total density
      massDensity = totalDens[iwelemRef][0];
      real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];


      // Step 2.2: divide the total mass/molar rate by the total density to get the total volumetric rate
      currentTotalVolRate = currentTotalRate * totalDensInv;


      if( logSurfaceCondition && useSurfaceConditions )
      {
        GEOS_LOG_RANK( GEOS_FMT( "{}: total fluid density at surface conditions = {} {}/sm3, total rate = {} {}/s, total surface volumetric rate = {} sm3/s",
                                 wellControlsName, totalDens[iwelemRef][0], massUnit, connRate[iwelemRef], massUnit, currentTotalVolRate ) );
      }

      // Step 3: update the phase volume rate
      for( integer ip = 0; ip < numPhase; ++ip )
      {

        // Step 3.1: compute the inverse of the (phase density * phase fraction)

        // skip the rest of this function if phase ip is absent
        bool const phaseExists = (phaseFrac[iwelemRef][0][ip] > 0);
        if( !phaseExists )
        {
          continue;
        }

        real64 const phaseDensInv =  1.0 / phaseDens[iwelemRef][0][ip];
        real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;

        // Step 3.2: divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
        currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
        if( currentPhaseVolRate[ip] > 0 && wellControlsName != "I1" && wellControlsName != "I2"  && wellControlsName != "I3" )
        {
          std::cout << "tjbwellphaserate " << wellControlsName << " " << ip << " "
                    << phaseDens[iwelemRef][0][ip] << " "
                    << phaseFrac[iwelemRef][0][ip] << " "
                    << currentPhaseVolRate[ip] << std::endl;
        }
        if( logSurfaceCondition && useSurfaceConditions )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: density of phase {} at surface conditions = {} {}/sm3, phase surface volumetric rate = {} sm3/s",
                                   wellControlsName, ip, phaseDens[iwelemRef][0][ip], massUnit, currentPhaseVolRate[ip] )  );
        }
      }
    } );
  } );
}


void CompositionalMultiphaseWell::updateFluidModel( WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 const > const & pres = subRegion.getField< well::pressure >();
  arrayView1d< real64 const > const & temp = subRegion.getField< well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< well::globalCompFraction >();

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  std::cout << "tjbwellupdatefluid " << subRegion.getName() << std::endl;
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    thermalCompositionalMultiphaseBaseKernels::
      FluidUpdateKernel::
      launch< ExecPolicy >( subRegion.size(),
                            fluidWrapper,
                            pres,
                            temp,
                            compFrac );
  } );

}

void CompositionalMultiphaseWell::updateSeparator( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  WellControls & wellControls = getWellControls( subRegion );

  // subRegion data
  arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const & temp = subRegion.getField< fields::well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();


  // fluid data
  constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();

  // control data

  string const wellControlsName = wellControls.getName();
  bool const logSurfaceCondition = isLogLevelActive< logInfo::BoundaryConditions >( wellControls.getLogLevel());
  string const massUnit = m_useMass ? "kg" : "mol";

  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 flashPressure;
  real64 flashTemperature;
  if( useSurfaceConditions )
  {
    // use surface conditions
    flashPressure = wellControls.getSurfacePressure();
    flashTemperature = wellControls.getSurfaceTemperature();
  }
  else
  {
    if( !wellControls.referenceReservoirRegion().empty() )
    {
      ElementRegionBase const & region = elemManager.getRegion( wellControls.referenceReservoirRegion() );
      GEOS_ERROR_IF ( !region.hasWrapper( CompositionalMultiphaseStatistics::regionStatisticsName()),
                      GEOS_FMT( "{}: WellControl {} referenceReservoirRegion field requires CompositionalMultiphaseStatistics to be configured for region {} ",
                                getDataContext(), wellControls.getName(), wellControls.referenceReservoirRegion() ) );

      CompositionalMultiphaseStatistics::RegionStatistics const & stats = region.getReference< CompositionalMultiphaseStatistics::RegionStatistics >(
        CompositionalMultiphaseStatistics::regionStatisticsName() );
      GEOS_ERROR_IF( stats.averagePressure <= 0.0,
                     GEOS_FMT(
                       "{}: No region average quantities computed.  WellControl {} referenceReservoirRegion field requires CompositionalMultiphaseStatistics to be configured for region {} ",
                       getDataContext(), wellControls.getName(), wellControls.referenceReservoirRegion() ));
      wellControls.setRegionAveragePressure( stats.averagePressure );
      wellControls.setRegionAverageTemperature( stats.averageTemperature );
    }
    // If flashPressure is not set by region the value is defaulted to -1 and indicates to use top segment conditions
    flashPressure = wellControls.getRegionAveragePressure();
    if( flashPressure < 0.0 )
    {
      // region name not set, use segment conditions
      flashPressure   = pres[iwelemRef];
      flashTemperature = temp[iwelemRef];
    }
    else
    {
      // use reservoir region averages
      flashTemperature = wellControls.getRegionAverageTemperature();
    }
  }

  constitutive::constitutiveUpdatePassThru( fluidSeparator, [&] ( auto & castedFluidSeparator )
  {
    // typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    typename TYPEOFREF( castedFluidSeparator ) ::KernelWrapper fluidSeparatorWrapper = castedFluidSeparator.createKernelWrapper();
    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [fluidSeparatorWrapper,
                                wellControlsName,
                                useSurfaceConditions,
                                flashPressure,
                                flashTemperature,
                                logSurfaceCondition,
                                iwelemRef,
                                pres,
                                temp,
                                compFrac] ( localIndex const )
    {
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element
      if( useSurfaceConditions )
      {
        // we need to compute the surface density
        //fluidWrapper.update( iwelemRef, 0, surfacePres, surfaceTemp, compFrac[iwelemRef] );
        fluidSeparatorWrapper.update( iwelemRef, 0, flashPressure, flashTemperature, compFrac[iwelemRef] );
        std::cout << "flash " << flashPressure << " " << flashTemperature << " " << compFrac[iwelemRef] << std::endl;
        if( logSurfaceCondition )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: surface density computed with P_surface = {} Pa and T_surface = {} K",
                                   wellControlsName, flashPressure, flashTemperature ) );
        }
#ifdef GEOS_USE_HIP
        GEOS_UNUSED_VAR( wellControlsName );
#endif
      }
      else
      {
        fluidSeparatorWrapper.update( iwelemRef, 0, flashPressure, flashTemperature, compFrac[iwelemRef] );
      }
    } );
  } );
}


real64 CompositionalMultiphaseWell::updatePhaseVolumeFraction( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  real64 maxDeltaPhaseVolFrac  =
    m_isThermal ?
    thermalCompositionalMultiphaseBaseKernels::
      PhaseVolumeFractionKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 subRegion,
                                                 fluid )
:    isothermalCompositionalMultiphaseBaseKernels::
      PhaseVolumeFractionKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 subRegion,
                                                 fluid );

  return maxDeltaPhaseVolFrac;
}

void CompositionalMultiphaseWell::updateTotalMassDensity( WellElementSubRegion & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  fluid.isThermal() ?
  thermalCompositionalMultiphaseWellKernels::
    TotalMassDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid )
  :
  compositionalMultiphaseWellKernels::
    TotalMassDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid );

}

real64 CompositionalMultiphaseWell::updateWellState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;

  // Update state and derived quantities for the well subregion
  real64 maxPhaseVolFrac =  updateSubRegionState( elemManager, subRegion );

  WellControls & wellControls = getWellControls( subRegion );
  WellConstraintBase * currentConstraint = wellControls.getCurrentConstraint();
  if( currentConstraint != nullptr )
  {
    // Store computed well quantities for this constraint
    currentConstraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
    currentConstraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                               CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
    currentConstraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                              CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
    currentConstraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
    std::cout << "updateWellState " << wellControls.getName() << " currentConstraint " << currentConstraint->getName() << " "
              << currentConstraint->bottomHolePressure() << " "
              << currentConstraint->phaseVolumeRates() << " "
              << currentConstraint->totalVolumeRate() << " "
              << currentConstraint->massRate() << std::endl;
  }
  else
  {
    std::cout << " no current constraint set for well " << wellControls.getName() << std::endl;
  }

  return maxPhaseVolFrac;

}

void CompositionalMultiphaseWell::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  //tjb
  real64 maxPhaseVolFrac = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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

        real64 const maxRegionPhaseVolFrac = updateWellState( elemManager, subRegion );

        maxPhaseVolFrac = LvArray::math::max( maxRegionPhaseVolFrac, maxPhaseVolFrac );
      }
    } );
  } );
  maxPhaseVolFrac = MpiWrapper::max( maxPhaseVolFrac );

  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                         GEOS_FMT( "        {}: Max well phase volume fraction change = {}",
                                   getName(), fmt::format( "{:.{}f}", maxPhaseVolFrac, 4 ) ) );

}

real64 CompositionalMultiphaseWell::updateSubRegionState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{
  real64 maxPhaseVolChange=0.0;
  WellControls & wellControls = getWellControls( subRegion );
  if( wellControls.getWellState())
  {

    // update properties
    updateGlobalComponentFraction( subRegion );

    // update densities, phase fractions, phase volume fractions

    updateFluidModel( subRegion );   //  Calculate fluid properties

    updateSeparator( elemManager, subRegion );   //  Calculate fluid properties at control conditions

    updateVolRatesForConstraint( subRegion );    // remove tjb ??

    maxPhaseVolChange = updatePhaseVolumeFraction( subRegion );
    updateTotalMassDensity( subRegion );

    // Calculate the reference element rates
    calculateReferenceElementRates( subRegion );

    // update the current BHP
    updateBHPForConstraint( subRegion );

    // Broad case the updated well state to other ranks
    real64 & currentBHP =
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
    array1d< real64 >   currentPhaseVolRate =
      wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
    real64 & currentTotalVolRate =
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
    real64 & currentMassRate =
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() );
    integer topRank = subRegion.getTopRank();
    MpiWrapper::broadcast( currentBHP, topRank );
    MpiWrapper::bcast( currentPhaseVolRate.data(), LvArray::integerConversion< int >( currentPhaseVolRate.size() ), topRank );
    MpiWrapper::broadcast( currentTotalVolRate, topRank );
    MpiWrapper::broadcast( currentMassRate, topRank );
    if( !subRegion.isLocallyOwned() )
    {
      wellControls.getReference< array1d< real64 > >(
        CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) =currentPhaseVolRate;


    }


  }
  return maxPhaseVolChange;
}

void CompositionalMultiphaseWell::initializeWell( DomainPartition & domain, MeshLevel & mesh, WellElementSubRegion & subRegion, real64 const & time_n )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( domain );
  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  m_nextDt = -1;
  // TODO: change the way we access the flowSolver here
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  compositionalMultiphaseWellKernels::PresTempCompFracInitializationKernel::CompFlowAccessors
    resCompFlowAccessors( elemManager, flowSolver.getName() );
  compositionalMultiphaseWellKernels::PresTempCompFracInitializationKernel::MultiFluidAccessors
    resMultiFluidAccessors( elemManager, flowSolver.getName() );

  WellControls & wellControls = getWellControls( subRegion );
  PerforationData const & perforationData = *subRegion.getPerforationData();
  arrayView2d< real64 const > const compPerfRate = perforationData.getField< fields::well::compPerforationRate >();

  bool const hasNonZeroRate = MpiWrapper::max< integer >( hasNonZero( compPerfRate ));

  if( time_n <= 0.0  || ( wellControls.isWellOpen(  ) && !hasNonZeroRate ) )
  {
    wellControls.setWellState( true );
    if( wellControls.getCurrentConstraint() == nullptr )
    {
      // tjb needed for backward compatibility. and these 2 lists must be consistent
      ConstraintTypeId inputControl = ConstraintTypeId( wellControls.getInputControl());
      if( wellControls.isProducer() )
      {
        wellControls.forSubGroups< MinimumBHPConstraint, ProductionConstraint< VolumeRateConstraint >, ProductionConstraint< MassRateConstraint >,
                                   ProductionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
        {
          if( constraint.getControl() == inputControl && constraint.isConstraintActive() )
          {
            wellControls.setCurrentConstraint( &constraint );
            wellControls.setControl( static_cast< WellControls::Control >(inputControl) );  // tjb old
          }

        } );
      }
      else
      {

        wellControls.forSubGroups< MaximumBHPConstraint, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint >,
                                   InjectionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
        {
          if( constraint.getControl() == inputControl && constraint.isConstraintActive() )
          {
            wellControls.setCurrentConstraint( &constraint );
            wellControls.setControl( static_cast< WellControls::Control >(inputControl) );   // tjb old
          }
        } );
      }
    }
    // get well primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure = subRegion.getField< well::pressure >();
    arrayView1d< real64 > const & wellElemTemp = subRegion.getField< well::temperature >();
    arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens = subRegion.getField< well::globalCompDensity >();
    arrayView1d< real64 > const & connRate = subRegion.getField< well::mixtureConnectionRate >();

    // get the info stored on well elements
    arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac = subRegion.getField< well::globalCompFraction >();
    arrayView1d< real64 const > const & wellElemGravCoef = subRegion.getField< well::gravityCoefficient >();

    // get the element region, subregion, index
    arrayView1d< localIndex const > const resElementRegion = perforationData.getField< perforation::reservoirElementRegion >();
    arrayView1d< localIndex const > const resElementSubRegion = perforationData.getField< perforation::reservoirElementSubRegion >();
    arrayView1d< localIndex const > const resElementIndex = perforationData.getField< perforation::reservoirElementIndex >();

    arrayView1d< real64 const > const & perfGravCoef = perforationData.getField< fields::well::gravityCoefficient >();
    arrayView1d< integer const > const & perfStatus = perforationData.getField< fields::perforation::perforationStatus >();
    arrayView1d< localIndex const > const & perfWellElemIndex =perforationData.getField< fields::perforation::wellElementIndex >();
    // 1) Loop over all perforations to compute an average mixture density and component fraction
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    compositionalMultiphaseWellKernels::
      PresTempCompFracInitializationKernel::
      launch( perforationData.size(),
              subRegion.size(),
              numComp,
              numPhase,
              wellControls,
              0.0,     // initialization done at t = 0
              resCompFlowAccessors.get( flow::pressure{} ),
              resCompFlowAccessors.get( flow::temperature{} ),
              resCompFlowAccessors.get( flow::globalCompDensity{} ),
              resCompFlowAccessors.get( flow::phaseVolumeFraction{} ),
              resMultiFluidAccessors.get( fields::multifluid::phaseMassDensity{} ),
              resElementRegion,
              resElementSubRegion,
              resElementIndex,
              perfGravCoef,
              perfStatus,
              perfWellElemIndex,
              wellElemGravCoef,
              wellElemPressure,
              wellElemTemp,
              wellElemCompFrac );

    // get well secondary variables on well elements
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & wellElemPhaseDens = fluid.phaseDensity();
    arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & wellElemTotalDens = fluid.totalDensity();

    // 4) Back calculate component densities
    constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      thermalCompositionalMultiphaseBaseKernels::
        FluidUpdateKernel::
        launch< serialPolicy >( subRegion.size(),
                                fluidWrapper,
                                wellElemPressure,
                                wellElemTemp,
                                wellElemCompFrac );
    } );

    compositionalMultiphaseWellKernels::
      CompDensInitializationKernel::launch( subRegion.size(),
                                            numComp,
                                            wellElemCompFrac,
                                            wellElemTotalDens,
                                            wellElemCompDens );

    // 5) Recompute the pressure-dependent properties

    updateSubRegionState( elemManager, subRegion );

    // 6) Estimate the well rates
    // TODO: initialize rates using perforation rates
    compositionalMultiphaseWellKernels::
      RateInitializationKernel::
      launch( subRegion.size(),
              wellControls,
              time_n,       // initialization done at time_n
              wellElemPhaseDens,
              wellElemTotalDens,
              connRate );

    updateVolRatesForConstraint( subRegion );
    //  Since this is a well manager class the rates need to be pushed into the WellControls class, which represnets the well
    WellConstraintBase * constraint =  wellControls.getCurrentConstraint();
    constraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
    constraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                        CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
    constraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                       CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
    constraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
    // 7) Copy well / fluid dofs to "prop"_n variables
    saveState( subRegion );
  }
  else if( !hasNonZeroRate )
  {
    wellControls.setWellState( false );
    GEOS_LOG_RANK_0( "tjb shut wells "<< subRegion.getName());
  }
  else
  {
    wellControls.setWellState( true );
    // setup if restart
    if( wellControls.getCurrentConstraint() == nullptr )
    {
      updateSubRegionState( elemManager, subRegion );
      if( wellControls.isProducer() )
      {
        wellControls.forSubGroups< MinimumBHPConstraint, ProductionConstraint< VolumeRateConstraint >, ProductionConstraint< MassRateConstraint >,
                                   ProductionConstraint< PhaseVolumeRateConstraint > >( [&](
                                                                                          auto
                                                                                          & constraint )
        {
          if( ConstraintTypeId( wellControls.getControl()) == constraint.getControl() && constraint.isConstraintActive() )
          {
            wellControls.setCurrentConstraint( &constraint );
          }
        } );
      }
      else
      {
        wellControls.forSubGroups< MaximumBHPConstraint, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint >, InjectionConstraint< PhaseVolumeRateConstraint > >( [&](
                                                                                                                                                                                                       auto
                                                                                                                                                                                                       &
                                                                                                                                                                                                       constraint )
        {
          if( ConstraintTypeId( wellControls.getControl()) == constraint.getControl() && constraint.isConstraintActive() )
          {
            wellControls.setCurrentConstraint( &constraint );
          }
        } );
      }
    }

  }
}


void CompositionalMultiphaseWell::initializeWells( DomainPartition & domain, real64 const & time_n )
{
  GEOS_MARK_FUNCTION;
  m_nextDt = -1;
  // TODO: change the way we access the flowSolver here
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  // loop over the wells
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();
    compositionalMultiphaseWellKernels::PresTempCompFracInitializationKernel::CompFlowAccessors
    resCompFlowAccessors( mesh.getElemManager(), flowSolver.getName() );
    compositionalMultiphaseWellKernels::PresTempCompFracInitializationKernel::MultiFluidAccessors
    resMultiFluidAccessors( mesh.getElemManager(), flowSolver.getName() );

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      initializeWell( domain, mesh, subRegion, time_n );
    } );

  } );
}

void CompositionalMultiphaseWell::assembleWellFluxTerms( real64 const & time,
                                                         real64 const & dt,
                                                         WellElementSubRegion const & subRegion,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time );

  BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
  if( m_useTotalMassEquation )
    kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );

  string const wellDofKey = dofManager.getKey( wellElementDofName());

  WellControls const & well_controls = getWellControls( subRegion );
  if( well_controls.isWellOpen( ) && !m_keepVariablesConstantDuringInitStep )
  {
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    int numComponents = fluid.numFluidComponents();

    if( isThermal() )
    {
      thermalCompositionalMultiphaseWellKernels::
        FaceBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                   dt,
                                                   dofManager.rankOffset(),
                                                   kernelFlags,
                                                   wellDofKey,
                                                   well_controls,
                                                   subRegion,
                                                   fluid,
                                                   localMatrix,
                                                   localRhs );
    }
    else
    {
      compositionalMultiphaseWellKernels::
        FaceBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                   dt,
                                                   dofManager.rankOffset(),
                                                   kernelFlags,
                                                   wellDofKey,
                                                   well_controls,
                                                   subRegion,
                                                   localMatrix,
                                                   localRhs );
    }
  }


}

void CompositionalMultiphaseWell::assembleFluxTerms( real64 const & time,
                                                     real64 const & dt,
                                                     DomainPartition & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time );

  BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
  if( m_useTotalMassEquation )
    kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );

  string const wellDofKey = dofManager.getKey( wellElementDofName());
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {

      assembleWellFluxTerms( time, dt, subRegion, dofManager, localMatrix, localRhs );

    } );
  } );

}

void CompositionalMultiphaseWell::assembleAccumulationTerms( real64 const & time,
                                                             real64 const & dt,
                                                             DomainPartition & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time );
  GEOS_UNUSED_VAR( dt );

  BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
  if( m_useTotalMassEquation )
    kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );

  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             WellElementSubRegion & subRegion )
    {

      assembleWellAccumulationTerms( time, dt, subRegion, dofManager, localMatrix, localRhs );
    } );
  } );


}

void CompositionalMultiphaseWell::assembleWellAccumulationTerms( real64 const & time,
                                                                 real64 const & dt,
                                                                 WellElementSubRegion & subRegion,
                                                                 DofManager const & dofManager,
                                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time );
  GEOS_UNUSED_VAR( dt );

  BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
  if( m_useTotalMassEquation )
    kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
  integer const numPhases = fluid.numFluidPhases();
  integer const numComponents = fluid.numFluidComponents();
  WellControls & wellControls = getWellControls( subRegion );
  if( wellControls.getWellStatus() == WellControls::Status::OPEN && !m_keepVariablesConstantDuringInitStep )
  {
    if( isThermal() )
    {

      thermalCompositionalMultiphaseWellKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                   numPhases,
                                                   wellControls,
                                                   wellControls.isProducer(),
                                                   dofManager.rankOffset(),
                                                   kernelFlags,
                                                   wellDofKey,
                                                   subRegion,
                                                   fluid,
                                                   localMatrix,
                                                   localRhs );
    }
    else
    {
      compositionalMultiphaseWellKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                   numPhases,
                                                   wellControls,
                                                   wellControls.isProducer(),
                                                   dofManager.rankOffset(),
                                                   kernelFlags,
                                                   wellDofKey,
                                                   subRegion,
                                                   fluid,
                                                   localMatrix,
                                                   localRhs );
    }
    // get the degrees of freedom and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const wellElemGhostRank = subRegion.ghostRank();
    arrayView1d< integer const > const elemStatus = subRegion.getLocalWellElementStatus();
    arrayView1d< real64 > const mixConnRate = subRegion.getField< fields::well::mixtureConnectionRate >();
    localIndex rank_offset = dofManager.rankOffset();
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( wellElemGhostRank[ei] < 0 )
      {
        if( elemStatus[ei]==WellElementSubRegion::WellElemStatus::CLOSED )
        {
          mixConnRate[ei] = 0.0;
          globalIndex const dofIndex = wellElemDofNumber[ei];
          localIndex const localRow = dofIndex - rank_offset;

          real64 const unity = 1.0;
          for( integer i=0; i < m_numDofPerWellElement; i++ )
          {
            globalIndex const rindex = localRow+i;
            globalIndex const cindex =dofIndex + i;
            localMatrix.template addToRow< serialAtomic >( rindex,
                                                           &cindex,
                                                           &unity,
                                                           1 );
            localRhs[rindex] = 0.0;
          }
        }
      }
    } );
  }
  else
  {
//wellControls.setWellOpen(false);
// get the degrees of freedom and ghosting info
    // get the degrees of freedom and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const wellElemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 >  mixConnRate = subRegion.getField< fields::well::mixtureConnectionRate >();
    localIndex rank_offset = dofManager.rankOffset();
    //bool wellClosed =  wellControls.getWellStatus() != WellControls::Status::OPEN;
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( wellElemGhostRank[ei] < 0 )
      {
        mixConnRate[ei] = 0.0;
        globalIndex const dofIndex = wellElemDofNumber[ei];
        localIndex const localRow = dofIndex - rank_offset;

        real64 const unity = 1.0;
        for( integer i=0; i < m_numDofPerWellElement; i++ )
        {
          globalIndex const rindex = localRow+i;
          globalIndex const cindex =dofIndex + i;
          localMatrix.template addToRow< serialAtomic >( rindex,
                                                         &cindex,
                                                         &unity,
                                                         1 );
          localRhs[rindex] = 0.0;
        }
      }
    } );
  }

}

void
CompositionalMultiphaseWell::applyWellBoundaryConditions( real64 const time_n,
                                                          real64 const dt,
                                                          ElementRegionManager & elemManager,
                                                          WellElementSubRegion & subRegion,
                                                          DofManager const & dofManager,
                                                          arrayView1d< real64 > const & localRhs,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix )
{
  GEOS_UNUSED_VAR( elemManager );
  GEOS_UNUSED_VAR( time_n );

  using namespace compositionalMultiphaseUtilities;

  BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
  if( useTotalMassEquation() )
    kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );

  integer const numComps = numFluidComponents();

  globalIndex const rankOffset = dofManager.rankOffset();


  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  // if the well is shut, we neglect reservoir-well flow that may occur despite the zero rate
  // therefore, we do not want to compute perforation rates and we simply assume they are zero
  WellControls const & wellControls = getWellControls( subRegion );
  //bool const detectCrossflow =
  //  ( wellControls.isInjector() ) && wellControls.isCrossflowEnabled() &&
  //  getLogLevel() >= 1;     // since detect crossflow requires communication, we detect it only if the logLevel is sufficiently high

  if( !wellControls.isWellOpen( ) )
  {
    return;
  }

  PerforationData const * const perforationData = subRegion.getPerforationData();

  // get the degrees of freedom
  string const wellDofKey = dofManager.getKey( wellElementDofName() );
#if 1
  if( isThermal ( )  )
  {
    coupledReservoirAndWellKernels::
      ThermalCompositionalMultiPhaseWellFluxKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( numComps,
                                                 wellControls,
                                                 wellControls.isProducer(),
                                                 dt,
                                                 rankOffset,
                                                 wellDofKey,
                                                 subRegion,
                                                 perforationData,
                                                 fluid,
                                                 kernelFlags,
                                                 localRhs,
                                                 localMatrix );

  }
  else
  {
    coupledReservoirAndWellKernels::
      IsothermalCompositionalMultiPhaseWellFluxKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( numComps,
                                                 dt,
                                                 rankOffset,
                                                 wellDofKey,
                                                 subRegion,
                                                 perforationData,
                                                 fluid,
                                                 localRhs,
                                                 localMatrix,
                                                 kernelFlags );
  }
#endif

}

void
CompositionalMultiphaseWell::applyBoundaryConditions( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                      real64 const GEOS_UNUSED_PARAM( dt ),
                                                      DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                      DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                      CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                      arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{}

real64
CompositionalMultiphaseWell::calculateWellResidualNorm( real64 const & time_n,
                                                        real64 const & dt,
                                                        WellElementSubRegion const & subRegion,
                                                        DofManager const & dofManager,
                                                        arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  integer numNorm = 1;     // mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;

  if( isThermal() )
  {
    numNorm = 2;     // mass balance and energy balance
  }
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );


  globalIndex const rankOffset = dofManager.rankOffset();
  string const wellDofKey = dofManager.getKey( wellElementDofName() );



  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  WellControls const & wellControls = getWellControls( subRegion );

  if( wellControls.isWellOpen( ) )
  {
    // step 1: compute the norm in the subRegion
    if( isThermal() )
    {
      real64 subRegionResidualNorm[2]{};

      thermalCompositionalMultiphaseWellKernels::ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   rankOffset,
                                                   wellDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   fluid,
                                                   wellControls,
                                                   time_n,
                                                   dt,
                                                   m_nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm );
      // step 2: reduction across meshBodies/regions/subRegions

      for( integer i=0; i<numNorm; i++ )
      {
        if( subRegionResidualNorm[i] > localResidualNorm[i] )
        {
          localResidualNorm[i] = subRegionResidualNorm[i];
        }
      }

    }
    else
    {
      real64 subRegionResidualNorm[1]{};
      compositionalMultiphaseWellKernels::ResidualNormKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   numDofPerWellElement(),
                                                   rankOffset,
                                                   wellDofKey,
                                                   localRhs,
                                                   subRegion,
                                                   fluid,
                                                   wellControls,
                                                   time_n,
                                                   dt,
                                                   m_nonlinearSolverParameters.m_minNormalizer,
                                                   subRegionResidualNorm );



      // step 2: reduction across meshBodies/regions/subRegions

      if( subRegionResidualNorm[0] > localResidualNorm[0] )
      {
        localResidualNorm[0] = subRegionResidualNorm[0];
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
  // step 3: second reduction across MPI ranks
  real64 resNorm=localResidualNorm[0];
  if( isThermal() )
  {
    real64 globalResidualNorm[2]{};
    globalResidualNorm[0] = MpiWrapper::max( localResidualNorm[0] );
    globalResidualNorm[1] = MpiWrapper::max( localResidualNorm[1] );
    resNorm=sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    GEOS_LOG_LEVEL_RANK_0( logInfo::ResidualNorm, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                                            coupledSolverAttributePrefix(), globalResidualNorm[0], globalResidualNorm[1] ));

  }
  else
  {
    resNorm= MpiWrapper::max( resNorm );

    GEOS_LOG_LEVEL_RANK_0( logInfo::ResidualNorm, GEOS_FMT( "        ( R{} ) = ( {:4.2e} )",
                                                            coupledSolverAttributePrefix(), resNorm ));
  }
  return resNorm;
}


real64
CompositionalMultiphaseWell::calculateResidualNorm( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer numNorm = 1;     // mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;

  if( isThermal() )
  {
    numNorm = 2;     // mass balance and energy balance
  }
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );


  globalIndex const rankOffset = dofManager.rankOffset();
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


      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      WellControls const & wellControls = getWellControls( subRegion );

      // step 1: compute the norm in the subRegion
      if( true )  // tjb wellControls.isWellOpen( ) )
      {
        if( isThermal() )
        {
          real64 subRegionResidualNorm[2]{};

          thermalCompositionalMultiphaseWellKernels::ResidualNormKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       rankOffset,
                                                       wellDofKey,
                                                       localRhs,
                                                       subRegion,
                                                       fluid,
                                                       wellControls,
                                                       time_n,
                                                       dt,
                                                       m_nonlinearSolverParameters.m_minNormalizer,
                                                       subRegionResidualNorm );
          // step 2: reduction across meshBodies/regions/subRegions

          for( integer i=0; i<numNorm; i++ )
          {
            if( subRegionResidualNorm[i] > localResidualNorm[i] )
            {
              localResidualNorm[i] = subRegionResidualNorm[i];
            }
          }

        }
        else
        {
          real64 subRegionResidualNorm[1]{};
          compositionalMultiphaseWellKernels::ResidualNormKernelFactory::
            createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                       numDofPerWellElement(),
                                                       rankOffset,
                                                       wellDofKey,
                                                       localRhs,
                                                       subRegion,
                                                       fluid,
                                                       wellControls,
                                                       time_n,
                                                       dt,
                                                       m_nonlinearSolverParameters.m_minNormalizer,
                                                       subRegionResidualNorm );



          // step 2: reduction across meshBodies/regions/subRegions

          if( subRegionResidualNorm[0] > localResidualNorm[0] )
          {
            localResidualNorm[0] = subRegionResidualNorm[0];
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
CompositionalMultiphaseWell::scalingForWellSystemSolution( ElementSubRegionBase & subRegion,
                                                           DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  real64 scalingFactor = 1.0;
  real64 maxDeltaPres = 0.0, maxDeltaCompDens = 0.0, maxDeltaTemp = 0.0;
  real64 minPresScalingFactor = 1.0, minCompDensScalingFactor = 1.0, minTempScalingFactor = 1.0;

  arrayView1d< real64 const > const pressure = subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const temperature = subRegion.getField< fields::well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::well::globalCompDensity >();
  arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::well::pressureScalingFactor >();
  arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::well::temperatureScalingFactor >();
  arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::well::globalCompDensityScalingFactor >();
  const integer temperatureOffset = m_numComponents+2;
  auto const subRegionData =
    m_isThermal
  ? thermalCompositionalMultiphaseBaseKernels::
      SolutionScalingKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                 m_maxAbsolutePresChange,
                                                 m_maxRelativeTempChange,
                                                 m_maxCompFracChange,
                                                 m_maxRelativeCompDensChange,
                                                 pressure,
                                                 temperature,
                                                 compDens,
                                                 pressureScalingFactor,
                                                 compDensScalingFactor,
                                                 temperatureScalingFactor,
                                                 dofManager.rankOffset(),
                                                 m_numComponents,
                                                 wellDofKey,
                                                 subRegion,
                                                 localSolution,
                                                 temperatureOffset )
  : isothermalCompositionalMultiphaseBaseKernels::
      SolutionScalingKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                 m_maxAbsolutePresChange,
                                                 m_maxCompFracChange,
                                                 m_maxRelativeCompDensChange,
                                                 pressure,
                                                 compDens,
                                                 pressureScalingFactor,
                                                 compDensScalingFactor,
                                                 dofManager.rankOffset(),
                                                 m_numComponents,
                                                 wellDofKey,
                                                 subRegion,
                                                 localSolution );


  scalingFactor = std::min( subRegionData.localMinVal, scalingFactor );

  maxDeltaPres  = std::max( maxDeltaPres, subRegionData.localMaxDeltaPres );
  maxDeltaCompDens = std::max( maxDeltaCompDens, subRegionData.localMaxDeltaCompDens );
  maxDeltaTemp = std::max( maxDeltaTemp, subRegionData.localMaxDeltaTemp );
  minPresScalingFactor = std::min( minPresScalingFactor, subRegionData.localMinPresScalingFactor );
  minCompDensScalingFactor = std::min( minCompDensScalingFactor, subRegionData.localMinCompDensScalingFactor );
  minTempScalingFactor = std::min( minTempScalingFactor, subRegionData.localMinTempScalingFactor );


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


  return LvArray::math::max( scalingFactor, m_minScalingFactor );

}

real64
CompositionalMultiphaseWell::scalingForSystemSolution( DomainPartition & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  real64 scalingFactor = 1.0;
  real64 maxDeltaPres = 0.0, maxDeltaCompDens = 0.0, maxDeltaTemp = 0.0;
  real64 minPresScalingFactor = 1.0, minCompDensScalingFactor = 1.0, minTempScalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure = subRegion.getField< well::pressure >();
      arrayView1d< real64 const > const temperature = subRegion.getField< well::temperature >();
      arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< well::globalCompDensity >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< well::pressureScalingFactor >();
      arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< well::temperatureScalingFactor >();
      arrayView1d< real64 > compDensScalingFactor = subRegion.getField< well::globalCompDensityScalingFactor >();
      const integer temperatureOffset = m_numComponents+2;
      auto const subRegionData =
        m_isThermal
  ? thermalCompositionalMultiphaseBaseKernels::
          SolutionScalingKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxAbsolutePresChange,
                                                     m_maxRelativeTempChange,
                                                     m_maxCompFracChange,
                                                     m_maxRelativeCompDensChange,
                                                     pressure,
                                                     temperature,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     temperatureScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     wellDofKey,
                                                     subRegion,
                                                     localSolution,
                                                     temperatureOffset )
  : isothermalCompositionalMultiphaseBaseKernels::
          SolutionScalingKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                     m_maxAbsolutePresChange,
                                                     m_maxCompFracChange,
                                                     m_maxRelativeCompDensChange,
                                                     pressure,
                                                     compDens,
                                                     pressureScalingFactor,
                                                     compDensScalingFactor,
                                                     dofManager.rankOffset(),
                                                     m_numComponents,
                                                     wellDofKey,
                                                     subRegion,
                                                     localSolution );


      scalingFactor = std::min( subRegionData.localMinVal, scalingFactor );

      maxDeltaPres  = std::max( maxDeltaPres, subRegionData.localMaxDeltaPres );
      maxDeltaCompDens = std::max( maxDeltaCompDens, subRegionData.localMaxDeltaCompDens );
      maxDeltaTemp = std::max( maxDeltaTemp, subRegionData.localMaxDeltaTemp );
      minPresScalingFactor = std::min( minPresScalingFactor, subRegionData.localMinPresScalingFactor );
      minCompDensScalingFactor = std::min( minCompDensScalingFactor, subRegionData.localMinCompDensScalingFactor );
      minTempScalingFactor = std::min( minTempScalingFactor, subRegionData.localMinTempScalingFactor );
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


  return LvArray::math::max( scalingFactor, m_minScalingFactor );

}

bool
CompositionalMultiphaseWell::checkWellSystemSolution( ElementSubRegionBase & subRegion,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  integer localCheck = 1;


  real64 minPres = 0.0, minDens = 0.0, minTotalDens = 0.0;
  integer numNegPres = 0, numNegDens = 0, numNegTotalDens = 0;

  const std::string wellName = subRegion.getName();

  //integer const m_allowCompDensChopping(true);
  integer const m_allowNegativePressure( false );
  compositionalMultiphaseUtilities::ScalingType const m_scalingType( compositionalMultiphaseUtilities::ScalingType::Global );
  arrayView1d< real64 const > const pressure =
    subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const temperature =
    subRegion.getField< fields::well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const compDens =
    subRegion.getField< fields::well::globalCompDensity >();
  arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::well::pressureScalingFactor >();
  arrayView1d< real64 > temperatureScalingFactor = subRegion.getField< fields::well::temperatureScalingFactor >();
  arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::well::globalCompDensityScalingFactor >();

  // check that pressure and component densities are non-negative
  // for thermal, check that temperature is above 273.15 K
  const integer temperatureOffset = m_numComponents+2;
  auto const subRegionData =
    m_isThermal
  ? thermalCompositionalMultiphaseBaseKernels::
      SolutionCheckKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                 m_allowNegativePressure,
                                                 m_scalingType,
                                                 scalingFactor,
                                                 pressure,
                                                 temperature,
                                                 compDens,
                                                 pressureScalingFactor,
                                                 temperatureScalingFactor,
                                                 compDensScalingFactor,
                                                 dofManager.rankOffset(),
                                                 m_numComponents,
                                                 wellDofKey,
                                                 subRegion,
                                                 localSolution,
                                                 temperatureOffset )
  : isothermalCompositionalMultiphaseBaseKernels::
      SolutionCheckKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_allowCompDensChopping,
                                                 m_allowNegativePressure,
                                                 m_scalingType,
                                                 scalingFactor,
                                                 pressure,
                                                 compDens,
                                                 pressureScalingFactor,
                                                 compDensScalingFactor,
                                                 dofManager.rankOffset(),
                                                 m_numComponents,
                                                 wellDofKey,
                                                 subRegion,
                                                 localSolution );

  localCheck = std::min( localCheck, subRegionData.localMinVal );

  minPres  = std::min( minPres, subRegionData.localMinPres );
  minDens = std::min( minDens, subRegionData.localMinDens );
  minTotalDens = std::min( minTotalDens, subRegionData.localMinTotalDens );
  numNegPres += subRegionData.localNumNegPressures;
  numNegDens += subRegionData.localNumNegDens;
  numNegTotalDens += subRegionData.localNumNegTotalDens;


  minPres  = MpiWrapper::min( minPres );
  minDens = MpiWrapper::min( minDens );
  minTotalDens = MpiWrapper::min( minTotalDens );
  numNegPres = MpiWrapper::sum( numNegPres );
  numNegDens = MpiWrapper::sum( numNegDens );
  numNegTotalDens = MpiWrapper::sum( numNegTotalDens );

  if( numNegPres > 0 )
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Number of negative well pressure values: {}, minimum value: {} Pa",
                                     subRegion.getName(), numNegPres, fmt::format( "{:.{}f}", minPres, 3 ) ) );
  string const massUnit = m_useMass ? "kg/m3" : "mol/m3";
  if( numNegDens > 0 )
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Number of negative well component density values: {}, minimum value: {} {} ",
                                     subRegion.getName(), numNegDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );
  if( minTotalDens > 0 )
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                           GEOS_FMT( "        {}: Number of negative total well density values: {}, minimum value: {} {} ",
                                     subRegion.getName(), minTotalDens, fmt::format( "{:.{}f}", minDens, 3 ), massUnit ) );



  return MpiWrapper::min( localCheck );
}

bool
CompositionalMultiphaseWell::checkSystemSolution( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );
  integer globalCheck = 1;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      integer localCheck = checkWellSystemSolution( subRegion, dofManager, localSolution, scalingFactor );
      globalCheck =  MpiWrapper::min( localCheck );
    } );
  } );
  return globalCheck;
}

void CompositionalMultiphaseWell::computeWellPerforationRates( real64 const & time_n,
                                                               real64 const & GEOS_UNUSED_PARAM( dt ),
                                                               ElementRegionManager const & elemManager,
                                                               WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );

  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  PerforationData * const perforationData = subRegion.getPerforationData();
  WellControls const & wellControls = getWellControls( subRegion );
  if( wellControls.isWellOpen() && !m_keepVariablesConstantDuringInitStep )
  {

    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    bool isThermal = fluid.isThermal();

    if( isThermal )
    {
      thermalPerforationFluxKernels::
        PerforationFluxKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   flowSolver.getName(),
                                                   perforationData,
                                                   subRegion,
                                                   fluid,
                                                   elemManager,
                                                   wellControls.isInjector(),
                                                   wellControls.isCrossflowEnabled());
    }
    else
    {
      isothermalPerforationFluxKernels::
        PerforationFluxKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   flowSolver.getName(),
                                                   perforationData,
                                                   subRegion,
                                                   elemManager,
                                                   wellControls.isInjector(),
                                                   wellControls.isCrossflowEnabled() );
    }
  }
  else
  {
// Zero completion flow rate
    arrayView2d< real64 > const compPerfRate = perforationData->getField< fields::well::compPerforationRate >();
    for( integer iperf=0; iperf<perforationData->size(); iperf++ )
    {
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        compPerfRate[iperf][ic] = 0.0;
      }
    }
  }


}


void CompositionalMultiphaseWell::computePerforationRates( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( dt );
  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {

    // TODO: change the way we access the flowSolver here

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion & subRegion )
    {
      computeWellPerforationRates( time_n, dt, elemManager, subRegion );

    } );

  } );

}

void
CompositionalMultiphaseWell::applyWellSystemSolution( DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution,
                                                      real64 const scalingFactor,
                                                      real64 const dt,
                                                      DomainPartition & domain,
                                                      MeshLevel & mesh,
                                                      WellElementSubRegion & subRegion )
{

  GEOS_UNUSED_VAR( domain );
  DofManager::CompMask pressureMask( m_numDofPerWellElement, 0, 1 );
  DofManager::CompMask componentMask( m_numDofPerWellElement, 1, numFluidComponents()+1 );
  DofManager::CompMask connRateMask( m_numDofPerWellElement, numFluidComponents()+1, numFluidComponents()+2 );
  GEOS_UNUSED_VAR( dt );
  // update all the fields using the global damping coefficients
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::pressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::globalCompDensity::key(),
                               scalingFactor,
                               componentMask );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               fields::well::mixtureConnectionRate::key(),
                               scalingFactor,
                               connRateMask );

  if( isThermal() )
  {
    DofManager::CompMask temperatureMask( m_numDofPerWellElement, numFluidComponents()+2, numFluidComponents()+3 );

    dofManager.addVectorToField( localSolution,
                                 wellElementDofName(),
                                 fields::well::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );

  }

#if 1
  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( subRegion );
  }
#endif

  // synchronize
  FieldIdentifiers fieldsToBeSync;
  if( isThermal() )
  {
    fieldsToBeSync.addElementFields( { fields::well::pressure::key(),
                                       fields::well::globalCompDensity::key(),
                                       fields::well::mixtureConnectionRate::key(),
                                       fields::well::temperature::key() },
                                     getTargetRegionNames() );
  }
  else
  {
    fieldsToBeSync.addElementFields( { fields::well::pressure::key(),
                                       fields::well::globalCompDensity::key(),
                                       fields::well::mixtureConnectionRate::key() },
                                     getTargetRegionNames() );
  }
  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                       mesh,
                                                       domain.getNeighbors(),
                                                       true );
}

void
CompositionalMultiphaseWell::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  real64 const dt,
                                                  DomainPartition & domain )
{


  DofManager::CompMask pressureMask( m_numDofPerWellElement, 0, 1 );
  DofManager::CompMask componentMask( m_numDofPerWellElement, 1, numFluidComponents()+1 );
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
                               well::globalCompDensity::key(),
                               scalingFactor,
                               componentMask );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               well::mixtureConnectionRate::key(),
                               scalingFactor,
                               connRateMask );

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
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    // synchronize
    FieldIdentifiers fieldsToBeSync;
    if( isThermal() )
    {
      fieldsToBeSync.addElementFields( { well::pressure::key(),
                                         well::globalCompDensity::key(),
                                         well::mixtureConnectionRate::key(),
                                         well::temperature::key() },
                                       regionNames );
    }
    else
    {
      fieldsToBeSync.addElementFields( { well::pressure::key(),
                                         well::globalCompDensity::key(),
                                         well::mixtureConnectionRate::key() },
                                       regionNames );
    }
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );


}

void CompositionalMultiphaseWell::chopNegativeDensities( WellElementSubRegion & subRegion )
{
  integer const numComp = m_numComponents;


  arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

  arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
    subRegion.getField< fields::well::globalCompDensity >();

  //arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens_n =
  // subRegion.getField< fields::well::globalCompDensity_n >();
  WellControls const & wellControls = getWellControls( subRegion );
  real64 minDensity = 1.0e-10;
  if( wellControls.isInjector() )
  {
    minDensity = 0.0;
  }
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
          wellElemCompDens[iwelem][ic] =minDensity;
        }
      }
    }
  } );


  if( isThermal() && !wellControls.isInjector()  )
  {
    real64 minTemp =290.0;
    arrayView1d< real64 > const & wellElemTemperature =
      subRegion.getField< well::temperature >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        if( wellElemTemperature[iwelem] < minTemp )
        {
          GEOS_LOG_LEVEL_BY_RANK( logInfo::Solution,
                                  GEOS_FMT( "{}: Bad temperature update, value set to injection temp = {} K",
                                            wellElemTemperature[iwelem], minTemp ) );
          wellElemTemperature[iwelem] = minTemp;
        }
      }
    } );
  }

}

void CompositionalMultiphaseWell::chopNegativeDensities( DomainPartition & domain )
{
  integer const numComp = m_numComponents;

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
      WellControls const & wellControls = getWellControls( subRegion );
      if( isThermal() && !wellControls.isInjector() )
      {

        real64 minTemp =290.0;
        arrayView1d< real64 > const & wellElemTemperature =
          subRegion.getField< well::temperature >();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
        {
          if( wellElemGhostRank[iwelem] < 0 )
          {
            if( wellElemTemperature[iwelem] < minTemp )
            {
              // GEOS_LOG_LEVEL_BY_RANK( logInfo::Solution,
              //                         GEOS_FMT( "{}: Bad temperature update, value set to injection temp = {} K",
              //                                   wellElemTemperature[iwelem], minTemp ) );
              wellElemTemperature[iwelem] = minTemp;
            }
          }
        } );
      }
    } );

  } );
}


void CompositionalMultiphaseWell::resetStateToBeginningOfStep( DomainPartition & domain )
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
      // get a reference to the primary variables on well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getField< well::pressure >();
      arrayView1d< real64 const > const & wellElemPressure_n =
        subRegion.getField< well::pressure_n >();
      wellElemPressure.setValues< parallelDevicePolicy<> >( wellElemPressure_n );

      if( isThermal() )
      {
        // get a reference to the primary variables on well elements
        arrayView1d< real64 > const & wellElemTemperature =
          subRegion.getField< well::temperature >();
        arrayView1d< real64 const > const & wellElemTemperature_n =
          subRegion.getField< well::temperature_n >();
        wellElemTemperature.setValues< parallelDevicePolicy<> >( wellElemTemperature_n );
      }
      arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity =
        subRegion.getField< well::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemGlobalCompDensity_n =
        subRegion.getField< well::globalCompDensity_n >();
      wellElemGlobalCompDensity.setValues< parallelDevicePolicy<> >( wellElemGlobalCompDensity_n );

      arrayView1d< real64 > const & connRate =
        subRegion.getField< well::mixtureConnectionRate >();
      arrayView1d< real64 const > const & connRate_n =
        subRegion.getField< well::mixtureConnectionRate_n >();
      connRate.setValues< parallelDevicePolicy<> >( connRate_n );
      WellControls & wellControls = getWellControls( subRegion );

      if( wellControls.isWellOpen( )  )
      {
        updateSubRegionState( elemManager, subRegion );
      }
    } );
  } );
}

void CompositionalMultiphaseWell::assembleWellConstraintTerms( real64 const & time_n,
                                                               real64 const & GEOS_UNUSED_PARAM( dt ),
                                                               WellElementSubRegion const & subRegion,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;


  // the rank that owns the reference well element is responsible for the calculations below.
  WellControls & wellControls = getWellControls( subRegion );
  if( !subRegion.isLocallyOwned() || !( wellControls.getWellStatus() == WellControls::Status::OPEN ))
  {
    return;
  }

  if( wellControls.isProducer() )
  {
    wellControls.forSubGroups< MinimumBHPConstraint, ProductionConstraint< PhaseVolumeRateConstraint >, ProductionConstraint< MassRateConstraint >, ProductionConstraint< VolumeRateConstraint >,
                               ProductionConstraint< LiquidRateConstraint >
                               >( [&]( auto & constraint )
    {
      if( constraint.getName() == wellControls.getCurrentConstraint()->getName())
      {
        // found limiting constraint
        std::cout << "Assembling constraint: " << constraint.getName() << std::endl;

        // fluid data
        constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();
        integer isThermal = fluidSeparator.isThermal();
        integer const numComp = fluidSeparator.numFluidComponents();
        geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
        {
          integer constexpr NUM_COMP = NC();
          integer constexpr IS_THERMAL = ISTHERMAL();

          wellConstraintKernels::ConstraintHelper< NUM_COMP, IS_THERMAL >::assembleConstraintEquation( time_n,
                                                                                                       wellControls,
                                                                                                       constraint,
                                                                                                       subRegion,
                                                                                                       dofManager.getKey( wellElementDofName() ),
                                                                                                       dofManager.rankOffset(),
                                                                                                       localMatrix,
                                                                                                       localRhs );
        } );
      }
    } );
  }
  else
  {
    wellControls.forSubGroups< MaximumBHPConstraint, InjectionConstraint< PhaseVolumeRateConstraint >, InjectionConstraint< MassRateConstraint >,
                               InjectionConstraint< VolumeRateConstraint >,
                               InjectionConstraint< LiquidRateConstraint >
                               >( [&]( auto & constraint )
    {
      if( constraint.getName() == wellControls.getCurrentConstraint()->getName())
      {
        // found limiting constraint

        // fluid data
        constitutive::MultiFluidBase & fluidSeparator =  wellControls.getMultiFluidSeparator();
        integer isThermal = fluidSeparator.isThermal();
        integer const numComp = fluidSeparator.numFluidComponents();
        geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
        {
          integer constexpr NUM_COMP = NC();
          integer constexpr IS_THERMAL = ISTHERMAL();

          wellConstraintKernels::ConstraintHelper< NUM_COMP, IS_THERMAL >::assembleConstraintEquation( time_n,
                                                                                                       wellControls,
                                                                                                       constraint,
                                                                                                       subRegion,
                                                                                                       dofManager.getKey( wellElementDofName() ),
                                                                                                       dofManager.rankOffset(),
                                                                                                       localMatrix,
                                                                                                       localRhs );
        } );
      }
    } );
  }
}

void CompositionalMultiphaseWell::assembleWellPressureRelations( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                                 real64 const & GEOS_UNUSED_PARAM( dt ),
                                                                 WellElementSubRegion const & subRegion,
                                                                 DofManager const & dofManager,
                                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                 arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;


  WellControls & wellControls = getWellControls( subRegion );

  if( wellControls.isWellOpen( ) && !m_keepVariablesConstantDuringInitStep )
  {
    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
    bool const isThermal = fluid.isThermal();
    // get the degrees of freedom, depth info, next welem index
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getField< well::gravityCoefficient >();
    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

    // get primary variables on well elements
    arrayView1d< real64 const > const & wellElemPres =
      subRegion.getField< well::pressure >();

    // get total mass density on well elements (for potential calculations)
    arrayView1d< real64 const > const & wellElemTotalMassDens =
      subRegion.getField< well::totalMassDensity >();
    arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens =
      subRegion.getField< well::dTotalMassDensity >();

    // segment status
    arrayView1d< integer const > const elemStatus =subRegion.getLocalWellElementStatus();

    bool controlHasSwitched = false;
    isothermalCompositionalMultiphaseBaseKernels::
      KernelLaunchSelectorCompTherm< compositionalMultiphaseWellKernels::PressureRelationKernel >
      ( numFluidComponents(),
      isThermal,
      subRegion.size(),
      dofManager.rankOffset(),
      elemStatus,
      wellElemDofNumber,
      wellElemGravCoef,
      nextWellElemIndex,
      wellElemPres,
      wellElemTotalMassDens,
      dWellElemTotalMassDens,
      controlHasSwitched,
      localMatrix,
      localRhs );

  }

}

void CompositionalMultiphaseWell::assemblePressureRelations( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_PARAM( dt )
  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel const & mesh,
                                                                 string_array const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      assembleWellPressureRelations( time_n, dt, subRegion, dofManager, localMatrix, localRhs );
    } );
  } );

}

void CompositionalMultiphaseWell::saveState( WellElementSubRegion & subRegion )
{


  // get a reference to the primary variables on well elements
  arrayView1d< real64 const > const & wellElemPressure =
    subRegion.getField< fields::well::pressure >();
  arrayView2d< real64 const, compflow::USD_COMP > const & wellElemGlobalCompDensity =
    subRegion.getField< fields::well::globalCompDensity >();
  arrayView1d< real64 const > const & wellElemTemperature =
    subRegion.getField< fields::well::temperature >();

  arrayView1d< real64 > const & wellElemPressure_n =
    subRegion.getField< fields::well::pressure_n >();
  wellElemPressure_n.setValues< parallelDevicePolicy<> >( wellElemPressure );

  if( isThermal() )
  {

    arrayView1d< real64 > const & wellElemTemperature_n =
      subRegion.getField< fields::well::temperature_n >();
    wellElemTemperature_n.setValues< parallelDevicePolicy<> >( wellElemTemperature );
  }

  arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity_n =
    subRegion.getField< fields::well::globalCompDensity_n >();
  wellElemGlobalCompDensity_n.setValues< parallelDevicePolicy<> >( wellElemGlobalCompDensity );

  arrayView1d< real64 const > const & connRate =
    subRegion.getField< fields::well::mixtureConnectionRate >();
  arrayView1d< real64 > const & connRate_n =
    subRegion.getField< fields::well::mixtureConnectionRate_n >();
  connRate_n.setValues< parallelDevicePolicy<> >( connRate );

  arrayView2d< real64 const, compflow::USD_PHASE > const wellElemPhaseVolFrac =
    subRegion.getField< fields::well::phaseVolumeFraction >();
  arrayView2d< real64, compflow::USD_PHASE > const wellElemPhaseVolFrac_n =
    subRegion.getField< fields::well::phaseVolumeFraction_n >();
  wellElemPhaseVolFrac_n.setValues< parallelDevicePolicy<> >( wellElemPhaseVolFrac );

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
  fluid.saveConvergedState();


}

void CompositionalMultiphaseWell::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  WellSolverBase::implicitStepSetup( time_n, dt, domain );

  forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&] ( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {

      WellControls & wellControls = getWellControls( subRegion );
      if( wellControls.isWellOpen() )
      {
        // get a reference to the primary variables on well elements
        arrayView1d< real64 const > const & wellElemPressure =
          subRegion.getField< fields::well::pressure >();
        arrayView2d< real64 const, compflow::USD_COMP > const & wellElemGlobalCompDensity =
          subRegion.getField< fields::well::globalCompDensity >();
        arrayView1d< real64 const > const & wellElemTemperature =
          subRegion.getField< fields::well::temperature >();

        arrayView1d< real64 > const & wellElemPressure_n =
          subRegion.getField< fields::well::pressure_n >();
        wellElemPressure_n.setValues< parallelDevicePolicy<> >( wellElemPressure );

        if( isThermal() )
        {
          arrayView1d< real64 > const & wellElemTemperature_n =
            subRegion.getField< fields::well::temperature_n >();
          wellElemTemperature_n.setValues< parallelDevicePolicy<> >( wellElemTemperature );
        }

        arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity_n =
          subRegion.getField< fields::well::globalCompDensity_n >();
        wellElemGlobalCompDensity_n.setValues< parallelDevicePolicy<> >( wellElemGlobalCompDensity );

        arrayView1d< real64 const > const & connRate =
          subRegion.getField< fields::well::mixtureConnectionRate >();
        arrayView1d< real64 > const & connRate_n =
          subRegion.getField< fields::well::mixtureConnectionRate_n >();
        connRate_n.setValues< parallelDevicePolicy<> >( connRate );

        arrayView2d< real64 const, compflow::USD_PHASE > const wellElemPhaseVolFrac =
          subRegion.getField< fields::well::phaseVolumeFraction >();
        arrayView2d< real64, compflow::USD_PHASE > const wellElemPhaseVolFrac_n =
          subRegion.getField< fields::well::phaseVolumeFraction_n >();
        wellElemPhaseVolFrac_n.setValues< parallelDevicePolicy<> >( wellElemPhaseVolFrac );

        string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
        MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
        fluid.saveConvergedState();

        validateWellConstraints( time_n, dt, subRegion );

        updateSubRegionState( elemManager, subRegion );
      }
    } )
    ;
  } );
}

void CompositionalMultiphaseWell::implicitStepComplete( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  WellSolverBase::implicitStepComplete( time_n, dt, domain );

  //if( getLogLevel() > 0 )
  //{
  printRates( time_n, dt, domain );
  //auto iterInfo = currentIter( time_n, dt );
  //printSegRates( time_n, dt, std::get< 0 >( iterInfo ), std::get< 1 >( iterInfo ), std::get< 2 >( iterInfo ), domain );
  //}
}

void CompositionalMultiphaseWell::printRates( real64 const & time_n,
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
      integer const numPhase = m_numPhases;
      integer const numComp = m_numComponents;
      integer const numPerf = subRegion.getPerforationData()->size();

      // control data
      WellControls const & wellControls = getWellControls( subRegion );

      stdVector< double > compRate( numComp, 0.0 );
      if( m_writeCSV > 0 && wellControls.isWellOpen( ) )
      {
        arrayView2d< real64 > const compPerfRate = subRegion.getPerforationData()->getField< fields::well::compPerforationRate >();

        // bring everything back to host, capture the scalars by reference
        forAll< serialPolicy >( 1, [&numComp,
                                    &numPerf,
                                    compPerfRate,
                                    &compRate] ( localIndex const )
        {
          for( integer ic = 0; ic < numComp; ++ic )
          {
            for( integer iperf = 0; iperf < numPerf; iperf++ )
            {
              compRate[ic] += compPerfRate[iperf][ic];
            }
          }
        } );
        for( integer ic = 0; ic < numComp; ++ic )
        {
          compRate[ic] = MpiWrapper::sum( compRate[ic] );
        }
      }

      // the rank that owns the reference well element is responsible for the calculations below.
      if( !subRegion.isLocallyOwned() )
      {
        return;
      }

      string const wellControlsName = wellControls.getName();

      // format: time,total_rate,total_vol_rate,phase0_vol_rate,phase1_vol_rate,...
      std::ofstream outputFile;
      if( m_writeCSV > 0 )
      {
        outputFile.open( m_ratesOutputDir + "/" + wellControlsName + ".csv", std::ios_base::app );
        outputFile << time_n << "," << dt;
      }

      if( wellControls.getWellStatus() == WellControls::Status::CLOSED )
      {
        GEOS_LOG( GEOS_FMT( "{}: well is shut", wellControlsName ) );
        if( outputFile.is_open())
        {
          // print all zeros in the rates file
          outputFile << ",0.0,0.0,0.0";
          if( wellControls.hasMinimumWHPConstraint() )
            outputFile << ",0.0";
          for( integer ip = 0; ip < numPhase; ++ip )
          {
            outputFile << ",0.0";
          }
          for( integer ic = 0; ic < numComp; ++ic )
          {
            outputFile << ",0.0";
          }
          outputFile << std::endl;
          outputFile.close();
        }
        return;
      }

      localIndex const iwelemRef = subRegion.getTopWellElementIndex();
      string const massUnit = m_useMass ? "kg" : "mol";

      // subRegion data

      arrayView1d< real64 const > const & connRate =
        subRegion.getField< well::mixtureConnectionRate >();

      integer const useSurfaceConditions = wellControls.useSurfaceConditions();
      real64 currentWHP =0.0;

      real64 const & currentBHP =
        wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
      if( wellControls.hasMinimumWHPConstraint() )
        currentWHP =
          wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentWHPString() );
      arrayView1d< real64 const > const & currentPhaseVolRate =
        wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
      real64 const & currentTotalVolRate =
        wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
      integer const hasWHP = wellControls.hasMinimumWHPConstraint();
      // bring everything back to host, capture the scalars by reference
      forAll< serialPolicy >( 1, [&numPhase,
                                  &numComp,
                                  &useSurfaceConditions,
                                  &currentBHP,
                                  &currentWHP,
                                  connRate,
                                  &currentTotalVolRate,
                                  currentPhaseVolRate,
                                  &compRate,
                                  &iwelemRef,
                                  &wellControlsName,
                                  &hasWHP,
                                  &massUnit,
                                  &outputFile] ( localIndex const )
      {
        string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
        string const unitKey = useSurfaceConditions ? "s" : "r";

        real64 const currentTotalRate = connRate[iwelemRef];
        GEOS_LOG( GEOS_FMT( "{}: BHP (at the specified reference elevation): {} Pa",
                            wellControlsName, currentBHP ) );
        GEOS_LOG( GEOS_FMT( "{}: Total rate: {} {}/s; total {} volumetric rate: {} {}m3/s",
                            wellControlsName, currentTotalRate, massUnit, conditionKey, currentTotalVolRate, unitKey ) );
        for( integer ip = 0; ip < numPhase; ++ip )
          GEOS_LOG( GEOS_FMT( "{}: Phase {} {} volumetric rate: {} {}m3/s",
                              wellControlsName, ip, conditionKey, currentPhaseVolRate[ip], unitKey ) );
        if( outputFile.is_open())
        {
          outputFile << "," << currentBHP;
          if( hasWHP )
            outputFile << "," << currentWHP;
          outputFile << "," << currentTotalRate << "," << currentTotalVolRate;
          for( integer ip = 0; ip < numPhase; ++ip )
          {
            outputFile << "," << currentPhaseVolRate[ip];
          }
          for( integer ic = 0; ic < numComp; ++ic )
          {
            outputFile << "," << compRate[ic];
          }
          outputFile << std::endl;
          outputFile.close();
        }
      } );
    } );
  } );
}


bool CompositionalMultiphaseWell::evaluateConstraints( real64 const & time_n,
                                                       real64 const & dt,
                                                       integer const cycleNumber,
                                                       integer const coupledIterationNumber,
                                                       DomainPartition & domain,
                                                       MeshLevel & mesh,
                                                       ElementRegionManager & elemManager,
                                                       WellElementSubRegion & subRegion,
                                                       DofManager const & dofManager )
{
  WellControls & wellControls = getWellControls( subRegion );
  bool useEstimator =   coupledIterationNumber <  wellControls.estimateSolution();


  if( useEstimator )
  {
    // create list of all constraints to solve
    // note that initializeWells sets the initial constraint
    // tjb reactive control schema to allow use to set if needed
    std::vector< WellConstraintBase * >  constraintList;
    WellConstraintBase *  limitingConstraint = wellControls.getCurrentConstraint();
    if( wellControls.isProducer() )
    {
      constraintList = wellControls.getProdRateConstraints();
      if( limitingConstraint->getControl() != ConstraintTypeId::BHP )
      {
        { // remove from list and add BHP constraint
          auto it = std::find( constraintList.begin(), constraintList.end(), limitingConstraint );
          if( it != constraintList.end() )
          {
            constraintList.erase( it );
          }
          if( wellControls.getMinBHPConstraint()->isConstraintActive() )
          {
            constraintList.push_back( wellControls.getMinBHPConstraint() );
          }
          constraintList.insert( constraintList.begin(), limitingConstraint );
        }
      }
      // Solve minimum bhp constraint first
      if( false && wellControls.getMinBHPConstraint()->isConstraintActive() )
      {
        // this is related to WHP option which introduces a new BHP constraint
        limitingConstraint = wellControls.getMinBHPConstraint();
      }
      else if( limitingConstraint == nullptr )
      {
        limitingConstraint = constraintList[0];
      }
    }
    else
    {
      constraintList = wellControls.getInjRateConstraints();
      // remove the limiting constraint from the list if present
      {
        if( limitingConstraint->getControl() != ConstraintTypeId::BHP )
        { // remove from list and add BHP constraint
          //auto it = std::find( constraintList.begin(), constraintList.end(), limitingConstraint );
          //if( it != constraintList.end() )
          //{
          //  constraintList.erase( it );
          //}
          constraintList.push_back( wellControls.getMaxBHPConstraint() );
        }
      }
    }
    if( wellControls.isoThermalEstimatorEnabled() )
    {
      wellControls.enableThermalEffects( false );
      solveConstraint ( limitingConstraint, time_n,
                        dt,
                        cycleNumber,
                        coupledIterationNumber,
                        domain,
                        mesh,
                        elemManager,
                        subRegion,
                        dofManager );
      wellControls.enableThermalEffects( true );
    }
    solveConstraint ( limitingConstraint, time_n,
                      dt,
                      cycleNumber,
                      coupledIterationNumber,
                      domain,
                      mesh,
                      elemManager,
                      subRegion,
                      dofManager );

    std::vector< int > constraintChecked( constraintList.size(), 0 );
    for( int i = 0; i < static_cast< int >(constraintList.size()); ++i )
    {
      auto & constraint = constraintList[i];
      GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                         " Well " << subRegion.getName() << " Constraint " << constraint->getName() << " active " << constraint->isConstraintActive() <<
                         " value " << constraint->getConstraintValue( time_n ) << " limiting " << limitingConstraint->getName() );
      if( limitingConstraint == nullptr )
      {
        std::cout << "error" << std::endl;
      }
      if( !constraintChecked[i] &&  constraint->isConstraintActive()  &&  constraint->checkViolation( *limitingConstraint, time_n ))
      {
        limitingConstraint=constraint;
        wellControls.setControl( static_cast< WellControls::Control >(constraint->getControl()) );                     // tjb old
        wellControls.setCurrentConstraint( limitingConstraint );
        GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                           " Well " << subRegion.getName() << " New Limiting Constraint " << constraint->getName() << " active " << constraint->isConstraintActive() <<
                           " value " << constraint->getConstraintValue( time_n ) );
        solveConstraint ( constraint, time_n,
                          dt,
                          cycleNumber,
                          coupledIterationNumber,
                          domain,
                          mesh,
                          elemManager,
                          subRegion,
                          dofManager );
        // tjb. this is likely not needed. set in update state
        constraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
        constraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
        constraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                           CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
        constraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));

      }
      constraintChecked[i]=1;
    }

    bool whpLimiting= solveWHPConstraint ( time_n,
                                           dt,
                                           cycleNumber,
                                           coupledIterationNumber,
                                           domain,
                                           mesh,
                                           elemManager,
                                           subRegion,
                                           dofManager );

    if( whpLimiting )
      limitingConstraint= wellControls.getCurrentConstraint();
    /*
       solveConstraint ( limitingConstraint, time_n,

                    dt,
                    cycleNumber,
                    coupledIterationNumber,
                    domain,
                    mesh,
                    elemManager,
                    subRegion,
                    dofManager );*/
    limitingConstraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
    limitingConstraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                                CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
    limitingConstraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                               CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
    limitingConstraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
    GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                       " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " << limitingConstraint->phaseVolumeRates() << " " <<
                       limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());
  }
  else
  {
    // create list of all constraints to process
    std::vector< WellConstraintBase * > constraintList;
    if( wellControls.isProducer() )
    {
      constraintList = wellControls.getProdRateConstraints();

      // Esitmator is not applied at this Newton iteration
      //if constraints arent updated with estimator and WHP is binding dont check allow constraint to be switch during remainder of timestep
      if( wellControls.hasMinimumWHPConstraint()  )
      {
        MinimumBHPConstraint * minBHPForWHP = wellControls.getMinimumBHPConstraintForWHP();
        if( minBHPForWHP != nullptr && minBHPForWHP->isConstraintActive())
        {
          std::cout << "we not active " << subRegion.getName() << " Constraint " << minBHPForWHP->getName() << " active " << minBHPForWHP->isConstraintActive() <<
            " value " << minBHPForWHP->getConstraintValue( time_n ) << std::endl;
          constraintList.insert( constraintList.begin(), minBHPForWHP );

        }
        else
        {
          ProductionConstraint< LiquidRateConstraint > * maxLiqForWHP = wellControls.getMaxLiquidConstraintForWHP();
          if( maxLiqForWHP != nullptr && maxLiqForWHP->isConstraintActive())
          {
            std::cout << "we  not active " << subRegion.getName() << " Constraint " << maxLiqForWHP->getName() << " active " << maxLiqForWHP->isConstraintActive() <<
              " value " << maxLiqForWHP->getConstraintValue( time_n ) << std::endl;
            constraintList.insert( constraintList.begin(), maxLiqForWHP );

          }
          else
          {
            // Solve minimum bhp constraint first
            if( wellControls.getMinBHPConstraint()->isConstraintActive() )
            {
              std::cout << "we  not active " << subRegion.getName() << " Constraint add minbp " << std::endl;
              constraintList.insert( constraintList.begin(), wellControls.getMinBHPConstraint() );
            }
          }
        }
      }

    }
    else
    {
      constraintList = wellControls.getInjRateConstraints();
      // Solve maximum bhp constraint first;
      constraintList.insert( constraintList.begin(), wellControls.getMaxBHPConstraint() );
    }
    // remove limiting constraint from list
    for( auto it = constraintList.begin(); it != constraintList.end(); ++it )
    {
      if( ( *it )->getName() == wellControls.getCurrentConstraint()->getName() )
      {
        constraintList.erase( it );
        break;
      }
    }
    // Get current constraint
    WellConstraintBase *  limitingConstraint = wellControls.getCurrentConstraint();
    for( auto & constraint : constraintList )
    {
      std::cout << " Well we=0 " << subRegion.getName() << " Constraint " << constraint->getName() << " active " << constraint->isConstraintActive() <<
        " value " << constraint->getConstraintValue( time_n ) << std::endl;
      if( constraint->getName() == wellControls.getCurrentConstraint()->getName())
      {
        limitingConstraint =  constraint;
        // tjb. this is likely not needed. set in update state

        GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                           " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " << limitingConstraint->phaseVolumeRates() << " " <<
                           limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());
      }
    }
     
    // Check current against other constraints
    std::vector< int > constraintChecked( constraintList.size(), 0 );
    for( int i = 0; i < static_cast< int >(constraintList.size()); ++i )
    {
      auto & constraint = constraintList[i];
      std::cout << " Well we=0 " << subRegion.getName() << " Checking Constraint " << constraint->getName() << " active " << constraint->isConstraintActive() <<
        " value " << constraint->getConstraintValue( time_n ) << " "  << limitingConstraint->getName() << std::endl;
      if( !constraintChecked[i] && limitingConstraint->getName() != constraint->getName())
      {
        if( constraint->checkViolation( *limitingConstraint, time_n ) )
        {
          limitingConstraint=constraint;
          wellControls.setControl( static_cast< WellControls::Control >(constraint->getControl()) );   // tjb old
          wellControls.setCurrentConstraint( constraint );
        constraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
        constraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
        constraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                           CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
        constraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));

          GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                             " Well " << subRegion.getName() << " New Limiting Constraint " << constraint->getName() << " "  << constraint->getConstraintValue( time_n )  );
          std::cout << wellControls.getName() <<  " New Limiting Constraint " << constraint->getName() << " "  << constraint->getConstraintValue( time_n )  << std::endl;
        }

      }
      constraintChecked[i]=1;
    }
    GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                       " Well " << subRegion.getName() << " Limiting Constraint " << limitingConstraint->getName() << " "  << limitingConstraint->bottomHolePressure() << " " << limitingConstraint->phaseVolumeRates() << " " <<
                       limitingConstraint->totalVolumeRate() << " " << limitingConstraint->massRate());

  }
  WellConstraintBase *  limitingConstraint = wellControls.getCurrentConstraint();
  std::cout << time_n << " " << dt <<" " << coupledIterationNumber << " " << wellControls.getName() <<  " limitingConstraint   " << wellControls.getCurrentConstraint()->getName() << " " <<
    limitingConstraint->getConstraintValue( time_n ) << " " << limitingConstraint->bottomHolePressure() << " " << limitingConstraint->phaseVolumeRates() << " " <<
    limitingConstraint->totalVolumeRate() <<   std::endl;
  return true;
}
void CompositionalMultiphaseWell::solveConstraint( WellConstraintBase *constraint,
                                                   real64 const & time_n,
                                                   real64 const & dt,
                                                   integer const cycleNumber,
                                                   integer const coupledIterationNumber,
                                                   DomainPartition & domain,
                                                   MeshLevel & mesh,
                                                   ElementRegionManager & elemManager,
                                                   WellElementSubRegion & subRegion,
                                                   DofManager const & dofManager )
{

  WellControls & wellControls = getWellControls( subRegion );
  bool useEstimator =   coupledIterationNumber <  wellControls.estimateSolution();
  if( useEstimator )
  {

    //if( getLogLevel() > 4 )
    {
      GEOS_LOG_RANK_0( "Well " << wellControls.getName() << " Evaluating constraint " << constraint->getName() << " value " << constraint->getConstraintValue( time_n ) << " active " <<
                       constraint->isConstraintActive() );
    }
    if( constraint->isConstraintActive() )
    {
      wellControls.setControl( static_cast< WellControls::Control >(constraint->getControl()) );     // tjb old
      wellControls.setCurrentConstraint( constraint );
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
      solveNonlinearSystem( time_n,
                            dt,
                            cycleNumber,
                            domain,
                            mesh,
                            elemManager,
                            subRegion,
                            dofManager );

      // Store computed well quantities for this constraint
      constraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      constraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                          CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      constraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                         CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      constraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
      if( getLogLevel() > 4 )
      {
        GEOS_LOG_RANK_0( "Well " << wellControls.getName() << " aft solve Constraint rates " << constraint->getName() << " bhp " << constraint->bottomHolePressure() << " phaseVolRate " <<
                         constraint->phaseVolumeRates() << " totalVolRate " << constraint->totalVolumeRate() << " massRate " << constraint->massRate());
      }
    }
  }
}

bool CompositionalMultiphaseWell::solveWHPConstraint( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer const cycleNumber,
                                                      integer const coupledIterationNumber,
                                                      DomainPartition & domain,
                                                      MeshLevel & mesh,
                                                      ElementRegionManager & elemManager,
                                                      WellElementSubRegion & subRegion,
                                                      DofManager const & dofManager )
{
  bool whpLimiting = false;
  WellControls & wellControls = getWellControls( subRegion );

  MinimumWHPConstraint * whpConstraint = wellControls.getMinWHPConstraint();
  if( whpConstraint == nullptr || !whpConstraint->isConstraintActive() )
    return whpLimiting;

  real64 & currentBHP =   wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  array1d< real64 > & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

  real64 currentBHP_local = currentBHP;
  array1d< real64 > currentPhaseVolRate_local = currentPhaseVolRate;
  real64 currentTotalVolRate_local = currentTotalVolRate;
  // Turn off BHP for WHP constraint if active, will be reset if WHP is limiting
  MinimumBHPConstraint * bhpConstraint=  wellControls.getMinimumBHPConstraintForWHP();
  bhpConstraint->setConstraintActive( false );
  real64 constraintWHP = whpConstraint->getConstraintValue( time_n );
  real64 currentWHP = constraintWHP;
  integer owner = -1;

  // Get the flow table function
  FunctionManager & functionManager = FunctionManager::getInstance();
  const PipeFlowTableFunction & m_flowTable =  functionManager.getGroup< PipeFlowTableFunction const >( whpConstraint->getFlowTableName());
  //m_flowTable.writeTable();
  integer flowTableSolveState;

  // this will be deleted with next merge
  if( subRegion.isLocallyOwned() )
  {
    owner = MpiWrapper::commRank( MPI_COMM_GEOS );
  }
  owner = MpiWrapper::max( owner );

  MpiWrapper::broadcast( currentWHP,
                         owner );

  // get current WHP from flow table
  m_flowTable.calculateWHP( wellControls.getName(), currentBHP, currentPhaseVolRate, currentWHP, flowTableSolveState );
  wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentWHPString() ) = currentWHP;
  //wellControls.getMinBHPConstraint()->setConstraintActive( true );

  std::cout << wellControls.getName() << " " <<coupledIterationNumber <<  " WHP constraint check " << whpLimiting << " " << currentWHP << " < " << constraintWHP << " bhp " << currentBHP <<
    " phase rates " << currentPhaseVolRate <<
    " total vol " << currentTotalVolRate <<
    " whpLimiting " << whpLimiting << std::endl;
  //whpLimiting=false;
  // check stability
  bool stabCheck = false;
  if( stabCheck )
  {

    real64 ql0, ql1, bhp0, bhp1;
    real64 dP_dQ_table = m_flowTable.calculatedPdQ( currentPhaseVolRate, currentWHP, ql0, ql1, bhp0, bhp1 );

    if( dP_dQ_table < 0.0 )
    {
      std::cout << time_n << " "  << wellControls.getName()<<  " VFPBracketSolve  negative dP/dQ table " << dP_dQ_table << " " << currentWHP << " " <<
        " " << ql0 <<  " " << ql1 << " " << bhp0 << " " << bhp1  << std::endl;
      ProductionConstraint< LiquidRateConstraint > *  liqConstraint=  wellControls.getMaxLiquidConstraintForWHP();
      wellControls.setCurrentConstraint( liqConstraint );
      liqConstraint->setConstraintActive( true );

      wellControls.setControl( static_cast< WellControls::Control >(liqConstraint->getControl()) );        // tjb old
      WellControls::Control wellControl = wellControls.getControl();
      MpiWrapper::broadcast( wellControl, owner );
      wellControls.setControl( wellControl );

      // lower bracker IPR solve
      liqConstraint->setConstraintValue( -ql0 );
      solveNonlinearSystem( time_n,
                            dt,
                            cycleNumber,
                            domain,
                            mesh,
                            elemManager,
                            subRegion,
                            dofManager );
      real64 iprBHP0 = currentBHP;
      std::cout << time_n << " IPRBracketSolve0 " << wellControls.getName() << " whp " << currentWHP << " BHP " << " " << currentBHP << " " << liqConstraint->bottomHolePressure() << " " << -ql0 <<
        " " << liqConstraint->phaseVolumeRates() << " " << liqConstraint->totalVolumeRate() << " " <<
        liqConstraint->massRate() << std::endl;
      // upper bracket IPR solve
      liqConstraint->setConstraintValue( -ql1 );
      solveNonlinearSystem( time_n,
                            dt,
                            cycleNumber,
                            domain,
                            mesh,
                            elemManager,
                            subRegion,
                            dofManager );
      real64 iprBHP1 = currentBHP;
      std::cout << time_n << " IPRBracketSolve1 " << wellControls.getName() << " whp " << currentWHP << " BHP  "  << currentBHP << " " << liqConstraint->bottomHolePressure() << " " << -ql1 << " " <<
        liqConstraint->phaseVolumeRates() << " " << liqConstraint->totalVolumeRate() << " " <<
        liqConstraint->massRate() << std::endl;
      liqConstraint->setConstraintActive( false );
      real64 dP_dQ_ipr = ( iprBHP1 - iprBHP0 ) / ( -ql1 - (-ql0) );
      if( dP_dQ_ipr > dP_dQ_table )
      {
        std::cout << time_n << " " << wellControls.getName() << " WHP 0 constraint stability dP/dQ table " << dP_dQ_table << " dP/dQ ipr " << dP_dQ_ipr  << " " << currentWHP << " " <<
          (currentWHP < constraintWHP) << std::endl;
        dP_dQ_table = dP_dQ_ipr;
        whpLimiting = currentWHP < constraintWHP;
      }
      else
      {
        // set so well operates at minwhp
        currentWHP = constraintWHP;
        std::cout << time_n << " " << wellControls.getName() << " WHP 1 constraint stability dP/dQ table " << dP_dQ_table << " dP/dQ ipr " << dP_dQ_ipr  << " " << currentWHP << " " <<
          (currentWHP < constraintWHP)  << std::endl;

        whpLimiting = true;
      }
      currentBHP = currentBHP_local;
      currentPhaseVolRate = currentPhaseVolRate_local;
      currentTotalVolRate = currentTotalVolRate_local;
    }
    else
    {
      // currentWHP is stable value
      whpLimiting = currentWHP < constraintWHP;
      std::cout << time_n << " " << wellControls.getName() << " WHP 2 constraint stability dP/dQ table " << dP_dQ_table << " " << currentWHP << " " << (currentWHP < constraintWHP)  << std::endl;


    }

  }
  else
  {
    // no stab check
    whpLimiting = currentWHP < constraintWHP;
  }

  if( whpLimiting )
  {

    std::cout << wellControls.getName() << " WHP constraint violated " << currentWHP << " < " << constraintWHP << " bhp " << currentBHP << " phase rates " << currentPhaseVolRate << " total vol " <<
      currentTotalVolRate <<
      std::endl;
    // WHP is limiting  set WHP to constraint value
    currentWHP = constraintWHP;



    // sets. tjb cleanup
    ProductionConstraint< LiquidRateConstraint > *  liqConstraint=  wellControls.getMaxLiquidConstraintForWHP();
    wellControls.setCurrentConstraint( liqConstraint );
    liqConstraint->setConstraintActive( true );
    wellControls.setControl( static_cast< WellControls::Control >(liqConstraint->getControl()) );         // tjb old
    WellControls::Control wellControl = wellControls.getControl();
    MpiWrapper::broadcast( wellControl, owner );
    wellControls.setControl( wellControl );
    std::ofstream of;
    of.open( "fl.csv" );
    of << "liq ,bhp ,tablebhp"<< std::endl;
    // Liquid constraint is used to find intersection of IPR and VLP
    const array1d< real64 > & liquidRates = m_flowTable.getRates();
    std::cout << liquidRates << std::endl;
    integer numRates = liquidRates.size();
    real64 liqRate0;
    real64 bhp0;
    real64 tableBHP0;

    liqRate0= liquidRates[0];
    for( integer i=0; i < 1; ++i )    // numRates; ++i )
    {
      of << liquidRates[i] << ",";
      liqRate0 = liquidRates[i];
      liqConstraint->setConstraintValue( liqRate0 );

      solveNonlinearSystem( time_n,
                            dt,
                            cycleNumber,
                            domain,
                            mesh,
                            elemManager,
                            subRegion,
                            dofManager );
      bhp0 = currentBHP;
      flowTableSolveState=1;
      m_flowTable.calculateBHP( currentPhaseVolRate, currentWHP, tableBHP0,
                                flowTableSolveState );
      std::cout << wellControls.getName() << " Solve at liquid rate [0] " << liqRate0 << " whp " << constraintWHP << " bhp " << bhp0 << " table bhp " << tableBHP0<< " phase rates " <<
        currentPhaseVolRate <<
        " total vol " << currentTotalVolRate << std::endl;
      of << bhp0 << "," << tableBHP0 << std::endl;
    }
    of.close();
    bool cSolve=false;
    integer currentRateIndex=numRates;
    while( !cSolve && currentRateIndex > 0 )
    {
      currentRateIndex--;
      liqConstraint->setConstraintValue( liquidRates[currentRateIndex] );
      cSolve = solveNonlinearSystem( time_n,
                                     dt,
                                     cycleNumber,
                                     domain,
                                     mesh,
                                     elemManager,
                                     subRegion,
                                     dofManager );

    }
    if( !cSolve )
    {
      throw("ft solve ");
    }
    real64 bhp1 = currentBHP;
    real64 liqRate1 = liquidRates[currentRateIndex];
    real64 tableBHP1;
    m_flowTable.calculateBHP( currentPhaseVolRate, currentWHP, tableBHP1,
                              flowTableSolveState );

    std::cout << time_n <<wellControls.getName() << "IPRSolve at liquid rate [index] " << currentRateIndex << " liqr " << liqRate1 << " whp  " << constraintWHP << " bhp " << bhp1 << " table bhp " <<
      tableBHP1<< " phase rates " <<
      currentPhaseVolRate <<
      " total vol " << currentTotalVolRate << " bhp res " << bhp1-tableBHP1 << std::endl;

    std::cout << wellControls.getName() << " Solve at liq brackets found " <<  std::endl;
    wellControls.setCurrentConstraint( bhpConstraint );
    wellControls.setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );

    integer const maxIters=100;
    real64 const tol = 1;
    integer iter = 0;

    of.open( "fls_"+wellControls.getName() + "_"+std::to_string( time_n )+"_"+std::to_string( dt )+"_"+std::to_string( coupledIterationNumber )+".csv" );
    of << " error,wellbhp,tablebhp,whp,orate,grate,wrate "  << std::endl;
    bhpConstraint->setConstraintActive( true );
    while( iter < maxIters && std::abs( tableBHP1 - bhp1 )  > tol )
    {
      // update whp
      bhp1=bhp1+0.50*(tableBHP1-bhp1);

      bhpConstraint->setConstraintValue( bhp1 );
      std::cout << "SolveLinear system " << iter << std::endl;
      solveNonlinearSystem( time_n,
                            dt,
                            cycleNumber,
                            domain,
                            mesh,
                            elemManager,
                            subRegion,
                            dofManager );
      //bhp1 = currentBHP;   // current(s) were updated in solveNonlinearSystem

      m_flowTable.calculateBHP( currentPhaseVolRate, currentWHP, tableBHP1,
                                flowTableSolveState );
      of << std::abs( bhp1-tableBHP1 ) << ","<<bhp1 << "," << tableBHP1<< ","<< currentWHP <<  "," << currentPhaseVolRate[0] << "," << currentPhaseVolRate[1] << " ," <<
        currentPhaseVolRate[2] << std::endl;
      std::cout <<  time_n << " " << wellControls.getName() << "IPRVFPSolve " << iter<< " whp  " <<  " whp  " << constraintWHP << " bhp " << bhp1 << " table bhp " << tableBHP1 << " phase rates " <<
        currentPhaseVolRate <<
        " total vol " << currentTotalVolRate << " bhp res " << bhp1-tableBHP1 << " " << (tableBHP1 - bhp1 )  << std::endl;

      bhpConstraint->setConstraintValue( bhp1 );

      ++iter;
    }
    of.close();
    if( false )
    {
      wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentWHPString() ) = currentWHP;
      // Use liquid constraint to find intersection of IPR and VLP
      wellControls.setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );                       // tjb old
      wellControls.setCurrentConstraint( bhpConstraint );
      wellControls.setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );     // tjb old
      wellControls.setCurrentConstraint( bhpConstraint );
      wellControls.getMinBHPConstraint()->setConstraintActive( false );
      bhpConstraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      bhpConstraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                             CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      bhpConstraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      bhpConstraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));

      std::cout << bhpConstraint->getName() << " " << bhpConstraint->bottomHolePressure() << " " << bhpConstraint->phaseVolumeRates() << " " << bhpConstraint->totalVolumeRate() << " " <<
        bhpConstraint->massRate() <<
        std::endl;
      GEOS_LOG_RANK_IF ( getLogLevel() > 4 && subRegion.isLocallyOwned(),
                         " Well " << subRegion.getName() << " New Limiting Constraint " << bhpConstraint->getName() << " active " << bhpConstraint->isConstraintActive() <<
                         " value " << bhpConstraint->getConstraintValue( time_n ) );
      std::cout << " WHP limiting solved in " << iter << " iters " << std::endl;
      std::cout << " Final WHP " << currentWHP << " BHP " << currentBHP << " phase rates " << currentPhaseVolRate << " total vol " << currentTotalVolRate << std::endl;
    }
    else
    {
      bhpConstraint->setConstraintActive( false );
      liqConstraint->setConstraintValue( -currentPhaseVolRate[0]- currentPhaseVolRate[2] );
      liqConstraint->setConstraintActive( true );
      wellControls.setCurrentConstraint( liqConstraint );
      wellControls.setControl( static_cast< WellControls::Control >(liqConstraint->getControl()) );       // tjb old
      liqConstraint->setBHP ( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      liqConstraint->setPhaseVolumeRates ( wellControls.getReference< array1d< real64 > >(
                                             CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      liqConstraint->setTotalVolumeRate ( wellControls.getReference< real64 >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      liqConstraint->setMassRate( wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
      std::cout << wellControls.getName() << " WHPConstraint " << liqConstraint->getName() << " limiting " << whpLimiting << " " << liqConstraint->bottomHolePressure() << " " <<
        liqConstraint->phaseVolumeRates() << " " << liqConstraint->totalVolumeRate() << " " <<
        liqConstraint->massRate() <<
        std::endl;

    }

  }
  return whpLimiting;
}
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, CompositionalMultiphaseWell, string const &, Group * const )
}   // namespace geos
