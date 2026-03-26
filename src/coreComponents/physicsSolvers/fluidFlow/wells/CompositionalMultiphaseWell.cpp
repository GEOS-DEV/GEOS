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
#include "physicsSolvers/fluidFlow/wells/WellWHPConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/ProdPipeFlowTableFunction.hpp"
#include "physicsSolvers/fluidFlow/wells/InjPipeFlowTableFunction.hpp"
#include "physicsSolvers/fluidFlow/wells/WellVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellPhaseVolumeRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellMassRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/WellLiquidRateConstraint.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellConstraintKernels.hpp"
#include "physicsSolvers/multiphysics/CoupledReservoirAndWellKernels.hpp"

#include "events/EventManager.hpp"
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
  WellControls( name, parent ),
  m_useTotalMassEquation( 1 ),
  m_maxCompFracChange( 1.0 ),
  m_maxRelativePresChange( 0.2 ),
  m_maxAbsolutePresChange( -1 ), // disabled by default
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 )
{

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
  WellControls::postInputInitialization();

  GEOS_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                        "The maximum absolute change in component fraction must smaller or equal to 1.0",
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) );
  GEOS_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                        "The maximum absolute change in component fraction must larger or equal to 0.0",
                        getWrapperDataContext( viewKeyStruct::maxCompFracChangeString() ) );

  GEOS_ERROR_IF_LE_MSG( m_maxRelativeCompDensChange, 0.0,
                        "The maximum relative change in component density must be larger than 0.0",
                        getWrapperDataContext( viewKeyStruct::maxRelativeCompDensChangeString() ) );

}
void CompositionalMultiphaseWell::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  setConstitutiveName< MultiFluidBase >( subRegion, viewKeyStruct::fluidNamesString(), "multiphase fluid" );
}

void CompositionalMultiphaseWell::registerWellDataOnMesh( WellElementSubRegion & subRegion )
{


  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  setConstitutiveNames ( subRegion );
  if( m_referenceFluidModelName.empty() )
  {
    m_referenceFluidModelName = getConstitutiveName< MultiFluidBase >( subRegion );
  }

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


  WellControls::registerWellDataOnMesh( subRegion );



  string const & fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
  // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.
  subRegion.registerField< well::pressure >( getName() );
  subRegion.registerField< well::pressure_n >( getName() );

  subRegion.registerField< well::temperature >( getName() );
  if( isThermal() )
  {
    subRegion.registerField< well::temperature_n >( getName() );
  }

  subRegion.registerField< well::gravityCoefficient >( getName() );


  subRegion.registerField< well::globalCompDensity >( getName() ).
    reference().resizeDimension< 1 >( m_numComponents );
  subRegion.registerField< well::globalCompDensity_n >( getName() ).
    reference().resizeDimension< 1 >( m_numComponents );

  subRegion.registerField< well::connectionRate >( getName() );
  subRegion.registerField< well::connectionRate_n >( getName() );

  subRegion.registerField< well::globalCompFraction >( getName() ).
    setDimLabels( 1, fluid.componentNames() ).
    reference().resizeDimension< 1 >( m_numComponents );
  subRegion.registerField< well::dGlobalCompFraction_dGlobalCompDensity >( getName() ).
    reference().resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

  subRegion.registerField< well::phaseVolumeFraction >( getName() ).
    setDimLabels( 1, fluid.phaseNames() ).
    reference().resizeDimension< 1 >( m_numPhases );
  subRegion.registerField< well::dPhaseVolumeFraction >( getName() ).
    reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents + 2 );     // dP, dT, dC

  subRegion.registerField< well::totalMassDensity >( getName() );
  subRegion.registerField< well::dTotalMassDensity >( getName() ).
    reference().resizeDimension< 1 >( m_numComponents +2 );     // dP, dT, dC

  subRegion.registerField< well::phaseVolumeFraction_n >( getName() ).
    reference().resizeDimension< 1 >( m_numPhases );

  subRegion.registerField< well::pressureScalingFactor >( getName() );
  subRegion.registerField< well::temperatureScalingFactor >( getName() );
  subRegion.registerField< well::globalCompDensityScalingFactor >( getName() );

  PerforationData & perforationData = *subRegion.getPerforationData();

  perforationData.registerField< well::gravityCoefficient >( getName() );
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

  registerWrapper< real64 >( viewKeyStruct::currentBHPString() );

  registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
    setSizedFromParent( 0 ).
    reference().resizeDimension< 0 >( m_numPhases );

  registerWrapper< real64 >( viewKeyStruct::massDensityString() );

  registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );

  registerWrapper< real64 >( viewKeyStruct::massDensityString() );

  registerWrapper< real64 >( viewKeyStruct::currentMassRateString() );

  // write rates output header
  // the rank that owns the reference well element is responsible
  if( m_writeCSV > 0 && subRegion.isLocallyOwned() )
  {
    string const fileName = GEOS_FMT( "{}/{}.csv", m_ratesOutputDir, getName() );
    string const massUnit = m_useMass ? "kg" : "mol";
    integer const useSurfaceConditions =  this->useSurfaceConditions();
    string const conditionKey = useSurfaceConditions ? "surface" : "reservoir";
    string const unitKey = useSurfaceConditions ? "s" : "r";
    integer const numPhase = m_numPhases;
    integer const numComp = m_numComponents;
    // format: time,bhp,total_rate,total_vol_rate,phase0_vol_rate,phase1_vol_rate,...
    makeDirsForPath( m_ratesOutputDir );
    GEOS_LOG( GEOS_FMT( "{}: Rates CSV generated at {}", getName(), fileName ) );
    std::ofstream outputFile( fileName );
    outputFile << "Time [s],dt[s],BHP [Pa]";

    if( hasMinimumWHPConstraint() || hasMaximumWHPConstraint() )
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


}


namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                        GEOS_FMT( "Mismatch in number of phases between constitutive models {} and {}",
                                  lhs.getName(), rhs.getName() ),
                        InputError, lhs.getDataContext(), rhs.getDataContext() );

  for( integer ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOS_THROW_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                          GEOS_FMT( "Mismatch in phase names between constitutive models {} and {}",
                                    lhs.getName(), rhs.getName() ),
                          InputError, lhs.getDataContext(), rhs.getDataContext() );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOS_THROW_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                        GEOS_FMT( "Mismatch in number of components between constitutive models {} and {}",
                                  lhs.getName(), rhs.getName() ),
                        InputError, lhs.getDataContext(), rhs.getDataContext() );

  for( integer ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOS_THROW_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                          GEOS_FMT( "Mismatch in component names between constitutive models {} and {}",
                                    lhs.getName(), rhs.getName() ),
                          InputError, lhs.getDataContext(), rhs.getDataContext() );
  }
}

}

/**
 * @brief Checks if the WellControls parameters are within the fluid tables ranges
 * @param fluid the fluid to check
 */
void CompositionalMultiphaseWell::validateFluidModel(
  constitutive::MultiFluidBase const & fluid, constitutive::MultiFluidBase const & referenceFluid ) const
{
  compareMultiphaseModels( fluid, referenceFluid );
  compareMulticomponentModels( fluid, referenceFluid );
  if( useSurfaceConditions() )
  {
    try
    {
      real64 const & surfaceTemp = getSurfaceTemperature();
      real64 const & surfacePres = getSurfacePressure();
      fluid.checkTablesParameters( surfacePres, surfaceTemp );
    } catch( SimulationError const & ex )
    {
      string const errorMsg = GEOS_FMT( "{}: wrong surface pressure / temperature.\n", getDataContext() );
      ErrorLogger::global().modifyCurrentExceptionMessage()
        .addToMsg( errorMsg )
        .addContextInfo( getDataContext().getContextInfo().setPriority( 1 ) );
      throw SimulationError( ex, errorMsg );
    }
  }
}

void CompositionalMultiphaseWell::validateWellConstraints( real64 const & time_n,
                                                           real64 const & GEOS_UNUSED_PARAM( dt ),
                                                           WellElementSubRegion const & subRegion )
{
  GEOS_UNUSED_VAR( time_n );

  if( !useSurfaceConditions() )
  {
    bool const useSeg = getReferenceReservoirRegion().empty();
    GEOS_WARNING_IF( useSeg,
                     "WellControls " <<WellControls::viewKeyStruct::referenceReservoirRegionString() <<
                     " not set and well constraint fluid property calculations will use top segement pressure and temp " );
    if( useSeg )
    {
      setRegionAveragePressure( -1 );
      setRegionAverageTemperature( -1 );
    }
    else
    {
      // Check if region name exists in list of Reservoir's target regions
      string const regionName =  getReferenceReservoirRegion();
      CompositionalMultiphaseBase const & flowSolver = getParent().getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
      string_array const & targetRegionsNames = flowSolver.getTargetRegionNames();
      auto const pos = std::find( targetRegionsNames.begin(), targetRegionsNames.end(), regionName );
      GEOS_ERROR_IF( pos == targetRegionsNames.end(),
                     GEOS_FMT( "{}: Region {} is not a target of the reservoir solver and cannot be used for referenceReservoirRegion in WellControl {}.",
                               getDataContext(), regionName, getName() ) );


    }
  }
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  forSubGroups< InjectionConstraint< PhaseVolumeRateConstraint >, ProductionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
  {
    constraint.validatePhaseType( fluid );
  } );



}

void CompositionalMultiphaseWell::initializePostSubGroups()
{
  WellControls::initializePostSubGroups();

}

void CompositionalMultiphaseWell::initializePostInitialConditionsPreSubGroups()
{
  WellControls::initializePostInitialConditionsPreSubGroups();

}

void CompositionalMultiphaseWell::initializeWellPostInitialConditionsPreSubGroups( WellElementSubRegion & subRegion )
{
  // set gravity coefficient
  setGravCoef( subRegion, getParent().getParent().getReference< R1Tensor >( PhysicsSolverManager::viewKeyStruct::gravityVectorString() ));

  // setup fluid model
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  constitutive::MultiFluidBase & fluid = subRegion.getConstitutiveModel< constitutive::MultiFluidBase >( fluidName );
  fluid.setMassFlag( m_useMass );
  createSeparator( subRegion );

  // Wellhead pressure constraints
  if( hasMinimumWHPConstraint())
  {
    createMinBHPConstraintForWHP();
    createMaxLiquidConstraintForWHP();
  }
  if( hasMaximumWHPConstraint())
  {
    createMaxBHPConstraintForWHP();
    createMaxVolumeInjConstraintForWHP();
  }

}
void CompositionalMultiphaseWell::postRestartInitialization( )
{
  // setup fluid separator
  constitutive::MultiFluidBase & fluidSeparator =   getMultiFluidSeparator();
  fluidSeparator.allocateConstitutiveData( *this, 1 );
  fluidSeparator.resize( 1 );
  // Wellhead pressure constraints
  if( hasMinimumWHPConstraint())
  {
    createMinBHPConstraintForWHP();
    createMaxLiquidConstraintForWHP();
  }
}

void CompositionalMultiphaseWell::createSeparator( WellElementSubRegion & subRegion )
{

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  fluid.setMassFlag( m_useMass );
  // setup fluid separator
  string const fluidSeparatorName = getName() + "Separator";
  std::unique_ptr< constitutive::ConstitutiveBase >  fluidSeparatorPtr  = fluid.deliverClone( fluidSeparatorName, this );
  fluidSeparatorPtr->allocateConstitutiveData( *this, 1 );
  fluidSeparatorPtr->resize( 1 );
  setFluidSeparator( std::move( fluidSeparatorPtr ));

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

  // subRegion data
  arrayView1d< real64 const > const & pres = subRegion.getField< well::pressure >();
  arrayView1d< real64 > const & totalMassDens = subRegion.getField< well::totalMassDensity >();
  arrayView1d< real64 const > const wellElemGravCoef = subRegion.getField< well::gravityCoefficient >();

  // control data

  string const wellControlsName =  getName();
  real64 const & refGravCoef = getReferenceGravityCoef();

  real64 & currentBHP =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
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
    //   std::cout << "tjbwellbhp " << pres[iwelemRef] << " " << totalMassDens[iwelemRef] << " "
    //         << wellElemGravCoef[iwelemRef] << " " << refGravCoef << " " << diffGravCoef << " " << currentBHP << std::endl;
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

  // subRegion data

  arrayView1d< real64 const > const & connRate = subRegion.getField< well::connectionRate >();

  // fluid data
  constitutive::MultiFluidBase & fluidSeparator =  getMultiFluidSeparator();

  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac = fluidSeparator.phaseFraction();
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & totalDens = fluidSeparator.totalDensity();
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens = fluidSeparator.phaseDensity();

  // control data

  string const wellControlsName = getName();

  arrayView1d< real64 > const & currentPhaseVolRate =
    getReference< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() );

  real64 & currentTotalVolRate =
    getReference< real64 >( viewKeyStruct::currentTotalVolRateString() );

  real64 & currentMassRate =
    getReference< real64 >( viewKeyStruct::currentMassRateString() );

  real64 & massDensity =
    getReference< real64 >( viewKeyStruct::massDensityString() );

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



  // subRegion data
  arrayView1d< real64 const > const & connRate = subRegion.getField< fields::well::connectionRate >();

  // fluid data
  constitutive::MultiFluidBase & fluidSeparator =      getMultiFluidSeparator();
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac = fluidSeparator.phaseFraction();
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & totalDens = fluidSeparator.totalDensity();
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens = fluidSeparator.phaseDensity();

  // control data
  string const wellControlsName =  getName();
  bool const logSurfaceCondition = isLogLevelActive< logInfo::BoundaryConditions >( getLogLevel());
  string const massUnit = m_useMass ? "kg" : "mol";

  integer const useSurfCond =  useSurfaceConditions();

  arrayView1d< real64 > const & currentPhaseVolRate =
    getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );

  real64 & currentTotalVolRate =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

  real64 & currentMassRate =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() );

  real64 & massDensity =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::massDensityString() );

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
                                &useSurfCond,
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


      if( logSurfaceCondition && useSurfCond )
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

        if( logSurfaceCondition && useSurfCond )
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
#if 0
  array2d< real64 > xcomp( 2, m_numComponents );
  arrayView2d< real64, compflow::USD_COMP >    compFrac1 = subRegion.getField< well::globalCompFraction >();

  compFrac1[0][0] = 0.565253; compFrac1[0][1] = 0.19773; compFrac1[0][2] = 0.225733; compFrac1[0][3] = 0.0107744; compFrac1[0][4] = 0.000509843;
  compFrac1[1][0] = 0.590029; compFrac1[1][1] = 0.205601; compFrac1[1][2] = 0.194712; compFrac1[1][3] = 0.0092061; compFrac1[1][4]   = 0.000452149;


  arrayView1d< real64 >   p = subRegion.getField< well::pressure >();
  arrayView1d< real64 >   T = subRegion.getField< well::temperature >();

  p[0] = 31209727.218176923692226; p[1] = 31311724.592292115092278;
  T[0] = 391.779; T[1] = 391.779;
#endif
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  //fluid.initializeState();  // tjb
  //std::cout << "Well " << getName() << " updating fluid model with pressure = " << pres << ", temperature = " << temp << ", compFrac = "
  // << compFrac    << std::endl;
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
#if 0
    thermalCompositionalMultiphaseBaseKernels::
      FluidUpdateKernel::
      launch< ExecPolicy >( subRegion.size(),
                            fluidWrapper,
                            p,
                            T,
                            compFrac1 );
#endif
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

  // subRegion data
  arrayView1d< real64 const > const & pres = subRegion.getField< fields::well::pressure >();
  arrayView1d< real64 const > const & temp = subRegion.getField< fields::well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac = subRegion.getField< fields::well::globalCompFraction >();


  // fluid data
  constitutive::MultiFluidBase & fluidSeparator =   getMultiFluidSeparator();

  // control data

  string const wellControlsName =  getName();
  bool const logSurfaceCondition = isLogLevelActive< logInfo::BoundaryConditions >( getLogLevel());
  string const massUnit = m_useMass ? "kg" : "mol";

  integer const useSurfaceCond =  useSurfaceConditions();
  real64 flashPressure;
  real64 flashTemperature;
  if( useSurfaceCond )
  {
    // use surface conditions
    flashPressure = getSurfacePressure();
    flashTemperature = getSurfaceTemperature();
  }
  else
  {
    if( !getReferenceReservoirRegion().empty() )
    {
      ElementRegionBase const & region = elemManager.getRegion( getReferenceReservoirRegion() );
      GEOS_ERROR_IF ( !region.hasWrapper( CompositionalMultiphaseStatistics::regionStatisticsName()),
                      GEOS_FMT( "{}: WellControl {} referenceReservoirRegion field requires CompositionalMultiphaseStatistics to be configured for region {} ",
                                getDataContext(), getName(), getReferenceReservoirRegion() ) );

      CompositionalMultiphaseStatistics::RegionStatistics const & stats = region.getReference< CompositionalMultiphaseStatistics::RegionStatistics >(
        CompositionalMultiphaseStatistics::regionStatisticsName() );
      GEOS_ERROR_IF( stats.averagePressure <= 0.0,
                     GEOS_FMT(
                       "{}: No region average quantities computed.  WellControl {} referenceReservoirRegion field requires CompositionalMultiphaseStatistics to be configured for region {} ",
                       getDataContext(), getName(), getReferenceReservoirRegion() ));
      setRegionAveragePressure( stats.averagePressure );
      setRegionAverageTemperature( stats.averageTemperature );
    }
    // If flashPressure is not set by region the value is defaulted to -1 and indicates to use top segment conditions
    flashPressure = getRegionAveragePressure();
    if( flashPressure < 0.0 )
    {
      // region name not set, use segment conditions
      flashPressure   = pres[iwelemRef];
      flashTemperature = temp[iwelemRef];
    }
    else
    {
      // use reservoir region averages
      flashTemperature = getRegionAverageTemperature();
    }
  }
  fluidSeparator.initializeState();  // tjb
  //std::cout << "Well " << getName() << " updating separator fluid model with pressure = " << pres << ", temperature = " << temp << ",
  // compFrac = " << compFrac    << std::endl;

  constitutive::constitutiveUpdatePassThru( fluidSeparator, [&] ( auto & castedFluidSeparator )
  {
    // typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    typename TYPEOFREF( castedFluidSeparator ) ::KernelWrapper fluidSeparatorWrapper = castedFluidSeparator.createKernelWrapper();
    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [fluidSeparatorWrapper,
                                wellControlsName,
                                useSurfaceCond,
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
      if( useSurfaceCond )
      {
        // we need to compute the surface density
        //fluidWrapper.update( iwelemRef, 0, surfacePres, surfaceTemp, compFrac[iwelemRef] );
        fluidSeparatorWrapper.update( iwelemRef, 0, flashPressure, flashTemperature, compFrac[iwelemRef] );
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

  real64 maxPhaseVolFrac =  updateSubRegionState( elemManager, subRegion );
  return maxPhaseVolFrac;

}

real64 CompositionalMultiphaseWell::updateSubRegionState( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
{
  real64 maxPhaseVolChange=0.0;

  if( getWellState())
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
      getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
    array1d< real64 >   currentPhaseVolRate =
      getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
    real64 & currentTotalVolRate =
      getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
    real64 & currentMassRate =
      getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() );
    integer topRank = subRegion.getTopRank();
    MpiWrapper::broadcast( currentBHP, topRank );
    MpiWrapper::bcast( currentPhaseVolRate.data(), LvArray::integerConversion< int >( currentPhaseVolRate.size() ), topRank );
    MpiWrapper::broadcast( currentTotalVolRate, topRank );
    MpiWrapper::broadcast( currentMassRate, topRank );
    if( !subRegion.isLocallyOwned() )
    {
      getReference< array1d< real64 > >(
        CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) =currentPhaseVolRate;

    }

    WellConstraintBase * constraint = getCurrentConstraint();
    std::cout << "Constraint ptr in well control is " << constraint->getName() << std::endl;
    if( constraint != nullptr )
    {
      constraint->setBHP ( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      constraint->setPhaseVolumeRates ( getReference< array1d< real64 > >(
                                          CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      constraint->setTotalVolumeRate ( getReference< real64 >(
                                         CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      constraint->setMassRate( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
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


  // TODO: change the way we access the flowSolver here
  ElementRegionManager const & elemManager = mesh.getElemManager();

  compositionalMultiphaseWellKernels::PresTempCompFracInitializationKernel::CompFlowAccessors
    resCompFlowAccessors( elemManager, getFlowSolverName() );
  compositionalMultiphaseWellKernels::PresTempCompFracInitializationKernel::MultiFluidAccessors
    resMultiFluidAccessors( elemManager, getFlowSolverName() );

  PerforationData const & perforationData = *subRegion.getPerforationData();
  arrayView2d< real64 const > const compPerfRate = perforationData.getField< fields::well::compPerforationRate >();

  bool const hasNonZeroRate = MpiWrapper::max< integer >( hasNonZero( compPerfRate ));
  std::cout << "Initializing well " << getName() << " at time " << time_n << " " <<  isWellOpen(  ) << " " << hasNonZeroRate << std::endl;
  if( time_n <= 0.0  || (  isWellOpen(  ) && !hasNonZeroRate ) )
  {
    setWellState( true );
    if( getCurrentConstraint() == nullptr )
    {
      // tjb needed for backward compatibility. and these 2 lists must be consistent
      ConstraintTypeId inputControl = ConstraintTypeId( getInputControl());
      if( isProducer() )
      {
        forSubGroups< MinimumBHPConstraint, ProductionConstraint< VolumeRateConstraint >, ProductionConstraint< MassRateConstraint >,
                      ProductionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
        {
          if( constraint.getControl() == inputControl && constraint.isConstraintActive())
          {
            setCurrentConstraint( &constraint );
            setControl( static_cast< WellControls::Control >(inputControl) );  // tjb old
          }
        } );
      }
      else
      {

        forSubGroups< MaximumBHPConstraint, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint >,
                      InjectionConstraint< PhaseVolumeRateConstraint > >( [&]( auto & constraint )
        {
          if( constraint.getControl() == inputControl && constraint.isConstraintActive())
          {
            setCurrentConstraint( &constraint );
            setControl( static_cast< WellControls::Control >(inputControl) );   // tjb old
          }
        } );
      }
    }

    GEOS_ERROR_IF( getCurrentConstraint() == nullptr, GEOS_FMT( "No active constraint found for well {} with input control {}", getName(), getInputControl() ) );
    // get well primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure = subRegion.getField< well::pressure >();
    arrayView1d< real64 > const & wellElemTemp = subRegion.getField< well::temperature >();
    arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens = subRegion.getField< well::globalCompDensity >();
    arrayView1d< real64 > const & connRate = subRegion.getField< well::connectionRate >();

    // get the info stored on well elements
    arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac = subRegion.getField< well::globalCompFraction >();
    arrayView1d< real64 const > const & wellElemGravCoef = subRegion.getField< well::gravityCoefficient >();

    // get the element region, subregion, index
    arrayView1d< localIndex const > const resElementRegion = perforationData.getField< perforation::reservoirElementRegion >();
    arrayView1d< localIndex const > const resElementSubRegion = perforationData.getField< perforation::reservoirElementSubRegion >();
    arrayView1d< localIndex const > const resElementIndex = perforationData.getField< perforation::reservoirElementIndex >();

    arrayView1d< real64 const > const & perfGravCoef = perforationData.getField< fields::well::gravityCoefficient >();
    arrayView1d< integer const > const & perfStatus = perforationData.getField< fields::perforation::perforationStatus >();

    // 1) Loop over all perforations to compute an average mixture density and component fraction
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    compositionalMultiphaseWellKernels::
      PresTempCompFracInitializationKernel::
      launch( perforationData.size(),
              subRegion.size(),
              numComp,
              numPhase,
              *this,
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
              *this,
              time_n,       // initialization done at time_n
              wellElemPhaseDens,
              wellElemTotalDens,
              connRate );

    updateVolRatesForConstraint( subRegion );
    //  Since this is a well manager class the rates need to be pushed into the WellControls class, which represnets the well
    WellConstraintBase * constraint =  getCurrentConstraint();
    constraint->setBHP ( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
    constraint->setPhaseVolumeRates ( getReference< array1d< real64 > >(
                                        CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
    constraint->setTotalVolumeRate ( getReference< real64 >(
                                       CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
    constraint->setMassRate( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
    // 7) Copy well / fluid dofs to "prop"_n variables
    saveState( subRegion );
  }
  else if( !hasNonZeroRate )
  {
    setWellState( false );
  }
  else
  {
    setWellState( true );
    // setup if restart
    if( getCurrentConstraint() == nullptr )
    {
      if( isProducer() )
      {
        forSubGroups< MinimumBHPConstraint, ProductionConstraint< VolumeRateConstraint >, ProductionConstraint< MassRateConstraint >,
                      ProductionConstraint< PhaseVolumeRateConstraint > >( [&](
                                                                             auto
                                                                             & constraint )
        {
          if( ConstraintTypeId( getControl()) == constraint.getControl() && constraint.isConstraintActive() )
          {
            setCurrentConstraint( &constraint );
          }
        } );
      }
      else
      {
        forSubGroups< MaximumBHPConstraint, InjectionConstraint< VolumeRateConstraint >, InjectionConstraint< MassRateConstraint >, InjectionConstraint< PhaseVolumeRateConstraint > >( [&](
                                                                                                                                                                                          auto
                                                                                                                                                                                          &
                                                                                                                                                                                          constraint )
        {
          if( ConstraintTypeId( getControl()) == constraint.getControl() && constraint.isConstraintActive() )
          {
            setCurrentConstraint( &constraint );
          }
        } );
      }
      updateSubRegionState( elemManager, subRegion );
    }

  }
  std::cout << "Finished initializing well " << getName() << " at time " << time_n << " " <<  int(getWellStatus()) << std::endl;
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


  if( isWellOpen( ) && !m_keepVariablesConstantDuringInitStep )
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
                                                   *this,
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
                                                   *this,
                                                   subRegion,
                                                   localMatrix,
                                                   localRhs );
    }
  }


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

  if( getWellStatus() == WellControls::Status::OPEN && !m_keepVariablesConstantDuringInitStep )
  {
    if( isThermal() )
    {

      thermalCompositionalMultiphaseWellKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( numComponents,
                                                   numPhases,
                                                   thermalEffectsEnabled(),
                                                   isProducer(),
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
                                                   thermalEffectsEnabled(),
                                                   isProducer(),
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
    arrayView1d< real64 > const mixConnRate = subRegion.getField< fields::well::connectionRate >();
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
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const wellElemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 >  mixConnRate = subRegion.getField< fields::well::connectionRate >();
    localIndex rank_offset = dofManager.rankOffset();
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
    // zero out current state constraint quantities
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() )=0.0;
    getReference< array1d< real64 > >(
      CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ).zero();
    getReference< real64 >(
      CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() )=0.0;
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() )=0.0;
  }
}

array1d< real64 >
CompositionalMultiphaseWell::calculateLocalWellResidualNorm( real64 const & time_n,
                                                             real64 const & dt,
                                                             NonlinearSolverParameters const & nonlinearSolverParameters,
                                                             WellElementSubRegion const & subRegion,
                                                             DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer const numNorm = isThermal() ? 2 : 1;     // mass balance and energy balance for thermal
  array1d< real64 > localResidualNorm( numNorm );
  array1d< real64 > localResidualNormalizer( numNorm );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const wellDofKey = dofManager.getKey( wellElementDofName() );



  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );


  // step 1: compute the norm in the subRegion

  if( !isWellOpen() )
  {
    for( integer i = 0; i < numNorm; ++i )
    {
      localResidualNorm[i] = 0.0;
      localResidualNormalizer[i] =  nonlinearSolverParameters.m_minNormalizer;
    }
  }
  else if( isThermal() )
  {
    real64 subRegionResidualNorm[2]{};

    thermalCompositionalMultiphaseWellKernels::ResidualNormKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 rankOffset,
                                                 wellDofKey,
                                                 localRhs,
                                                 subRegion,
                                                 fluid,
                                                 *this,
                                                 time_n,
                                                 dt,
                                                 nonlinearSolverParameters.m_minNormalizer,
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
                                                 m_numDofPerWellElement,
                                                 rankOffset,
                                                 wellDofKey,
                                                 localRhs,
                                                 subRegion,
                                                 fluid,
                                                 *this,
                                                 time_n,
                                                 dt,
                                                 nonlinearSolverParameters.m_minNormalizer,
                                                 subRegionResidualNorm );



    // step 2: reduction across meshBodies/regions/subRegions

    if( subRegionResidualNorm[0] > localResidualNorm[0] )
    {
      localResidualNorm[0] = subRegionResidualNorm[0];
    }
  }

  return localResidualNorm;

}

real64
CompositionalMultiphaseWell::calculateWellResidualNorm( real64 const & time_n,
                                                        real64 const & dt,
                                                        NonlinearSolverParameters const & nonlinearSolverParameters,
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


  //globalIndex const rankOffset = dofManager.rankOffset();
  string const wellDofKey = dofManager.getKey( wellElementDofName() );



  //string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  //MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  if( isWellOpen( ) )
  {
    localResidualNorm = calculateLocalWellResidualNorm( time_n,
                                                        dt,
                                                        nonlinearSolverParameters,
                                                        subRegion,
                                                        dofManager,
                                                        localRhs );

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
CompositionalMultiphaseWell::scalingForLocalSystemSolution( WellElementSubRegion & subRegion,
                                                            DofManager const & dofManager,
                                                            real64 & maxDeltaPres,
                                                            real64 & maxDeltaCompDens,
                                                            real64 & maxDeltaTemp,
                                                            real64 & minPresScalingFactor,
                                                            real64 & minCompDensScalingFactor,
                                                            real64 & minTempScalingFactor,
                                                            arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  real64 scalingFactor = 1.0;
  maxDeltaPres = 0.0;
  maxDeltaCompDens = 0.0;
  maxDeltaTemp = 0.0;
  minPresScalingFactor = 1.0;
  minCompDensScalingFactor = 1.0;
  minTempScalingFactor = 1.0;

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
CompositionalMultiphaseWell::scalingForWellSystemSolution( WellElementSubRegion & subRegion,
                                                           DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  real64 maxDeltaPres = 0.0, maxDeltaCompDens = 0.0, maxDeltaTemp = 0.0;
  real64 minPresScalingFactor = 1.0, minCompDensScalingFactor = 1.0, minTempScalingFactor = 1.0;

  real64 scalingFactor =scalingForLocalSystemSolution( subRegion,
                                                       dofManager,
                                                       maxDeltaPres,
                                                       maxDeltaCompDens,
                                                       maxDeltaTemp,
                                                       minPresScalingFactor,
                                                       minCompDensScalingFactor,
                                                       minTempScalingFactor,
                                                       localSolution );
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
CompositionalMultiphaseWell::checkWellSystemSolution( WellElementSubRegion & subRegion,
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

void CompositionalMultiphaseWell::computeWellPerforationRates( real64 const & time_n,
                                                               real64 const & GEOS_UNUSED_PARAM( dt ),
                                                               ElementRegionManager & elemManager,
                                                               WellElementSubRegion & subRegion )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n );

  //CompositionalMultiphaseBase const & flowSolver = getParent().getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  PerforationData * const perforationData = subRegion.getPerforationData();

  if( isWellOpen() && !m_keepVariablesConstantDuringInitStep )
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
                                                   getFlowSolverName(),
                                                   perforationData,
                                                   subRegion,
                                                   fluid,
                                                   elemManager,
                                                   isInjector(),
                                                   isCrossflowEnabled());
    }
    else
    {
      isothermalPerforationFluxKernels::
        PerforationFluxKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                   m_numPhases,
                                                   getFlowSolverName(),
                                                   perforationData,
                                                   subRegion,
                                                   elemManager,
                                                   isInjector(),
                                                   isCrossflowEnabled() );
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

  //bool const detectCrossflow =
  //  ( isInjector() ) && isCrossflowEnabled() &&
  //  getLogLevel() >= 1;     // since detect crossflow requires communication, we detect it only if the logLevel is sufficiently high

  if( !isWellOpen( ) )
  {
    return;
  }

  PerforationData const * const perforationData = subRegion.getPerforationData();

  // get the degrees of freedom
  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  if( isThermal ( )  )
  {
    coupledReservoirAndWellKernels::
      ThermalCompositionalMultiPhaseWellFluxKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( numComps,
                                                 *this,
                                                 isProducer(),
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
                               fields::well::connectionRate::key(),
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

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  //if( m_allowCompDensChopping )
  {
    chopNegativeDensities( subRegion );
  }

  // synchronize
  FieldIdentifiers fieldsToBeSync;
  if( isThermal() )
  {
    fieldsToBeSync.addElementFields( { fields::well::pressure::key(),
                                       fields::well::globalCompDensity::key(),
                                       fields::well::connectionRate::key(),
                                       fields::well::temperature::key() },
                                     getTargetRegionNames() );
  }
  else
  {
    fieldsToBeSync.addElementFields( { fields::well::pressure::key(),
                                       fields::well::globalCompDensity::key(),
                                       fields::well::connectionRate::key() },
                                     getTargetRegionNames() );
  }
  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                       mesh,
                                                       domain.getNeighbors(),
                                                       true );

}



void CompositionalMultiphaseWell::chopNegativeDensities( WellElementSubRegion & subRegion )
{
  integer const numComp = m_numComponents;


  arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

  arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
    subRegion.getField< fields::well::globalCompDensity >();
  // arrayView1d< real64 > const & wellElemTemperature =
  //  subRegion.getField< well::temperature >();
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    if( wellElemGhostRank[iwelem] < 0 )
    {
      //std::cout << " Temperature: " << iwelem << " " << wellElemTemperature[iwelem] << std::endl;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        // we allowed for some densities to be slightly negative in CheckSystemSolution
        // if the new density is negative, chop back to zero
        if( wellElemCompDens[iwelem][ic] < 0 )
        {
          std::cout << "Chopping negative density to zero for well element " << iwelem << " component " << ic << " value before chopping: " << wellElemCompDens[iwelem][ic] << std::endl;
          wellElemCompDens[iwelem][ic] = 0.0;
        }
      }
    }
  } );

}


void CompositionalMultiphaseWell::resetStateToBeginningOfStep( ElementRegionManager const & elemManager, WellElementSubRegion & subRegion )
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
    subRegion.getField< well::connectionRate >();
  arrayView1d< real64 const > const & connRate_n =
    subRegion.getField< well::connectionRate_n >();
  connRate.setValues< parallelDevicePolicy<> >( connRate_n );


  if( isWellOpen( )  )
  {
    updateSubRegionState( elemManager, subRegion );
  }

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

  if( !subRegion.isLocallyOwned() || !(  getWellStatus() == WellControls::Status::OPEN ))
  {
    return;
  }

  if( isProducer() )
  {
    forSubGroups< MinimumBHPConstraint, ProductionConstraint< PhaseVolumeRateConstraint >, ProductionConstraint< MassRateConstraint >, ProductionConstraint< VolumeRateConstraint >,
                  ProductionConstraint< LiquidRateConstraint >
                  >( [&]( auto & constraint )
    {
      // Need to use name since there could be multiple constraints of the same type
      if( constraint.getName() ==  getCurrentConstraint()->getName())
      {
        //std::cout << "Assembling constraint: " << constraint.getName() << std::endl;
        // found limiting constraint
        constitutive::MultiFluidBase & fluidSeparator =   getMultiFluidSeparator();
        integer isThermal = fluidSeparator.isThermal();
        integer const numComp = fluidSeparator.numFluidComponents();
        geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
        {
          integer constexpr NUM_COMP = NC();
          integer constexpr IS_THERMAL = ISTHERMAL();

          wellConstraintKernels::ConstraintHelper< NUM_COMP, IS_THERMAL >::assembleConstraintEquation( time_n,
                                                                                                       *this,
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
    forSubGroups< MaximumBHPConstraint, InjectionConstraint< PhaseVolumeRateConstraint >, InjectionConstraint< MassRateConstraint >,
                  InjectionConstraint< VolumeRateConstraint >,
                  InjectionConstraint< LiquidRateConstraint >
                  >( [&]( auto & constraint )
    {
      if( constraint.getName() ==  getCurrentConstraint()->getName())
      {
        // found limiting constraint
        constitutive::MultiFluidBase & fluidSeparator =   getMultiFluidSeparator();
        integer isThermal = fluidSeparator.isThermal();
        integer const numComp = fluidSeparator.numFluidComponents();
        geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
        {
          integer constexpr NUM_COMP = NC();
          integer constexpr IS_THERMAL = ISTHERMAL();

          wellConstraintKernels::ConstraintHelper< NUM_COMP, IS_THERMAL >::assembleConstraintEquation( time_n,
                                                                                                       *this,
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

  if( isWellOpen( ) && !m_keepVariablesConstantDuringInitStep )
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
    subRegion.getField< fields::well::connectionRate >();
  arrayView1d< real64 > const & connRate_n =
    subRegion.getField< fields::well::connectionRate_n >();
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
                                                     ElementRegionManager & elemManager,
                                                     WellElementSubRegion & subRegion )
{
  WellControls::implicitStepSetup( time_n, dt, elemManager, subRegion );

  if( isWellOpen() )
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
      subRegion.getField< fields::well::connectionRate >();
    arrayView1d< real64 > const & connRate_n =
      subRegion.getField< fields::well::connectionRate_n >();
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

}

void CompositionalMultiphaseWell::implicitStepComplete( real64 const & time_n,
                                                        real64 const & dt,
                                                        WellElementSubRegion const & subRegion )
{
  printRates( time_n, dt, subRegion );
}

void CompositionalMultiphaseWell::printRates( real64 const & time_n,
                                              real64 const & dt,
                                              WellElementSubRegion const & subRegion )
{

  integer const numPhase = m_numPhases;
  integer const numComp = m_numComponents;
  integer const numPerf = subRegion.getPerforationData()->size();



  stdVector< double > compRate( numComp, 0.0 );
  if( m_writeCSV > 0 &&  isWellOpen( ) )
  {
    arrayView2d< real64 const > const & compPerfRate = subRegion.getPerforationData()->getField< fields::well::compPerforationRate >();

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

  string const wellControlsName =  getName();

  // format: time,total_rate,total_vol_rate,phase0_vol_rate,phase1_vol_rate,...
  std::ofstream outputFile;
  if( m_writeCSV > 0 )
  {
    outputFile.open( m_ratesOutputDir + "/" + wellControlsName + ".csv", std::ios_base::app );
    outputFile << time_n << "," << dt;
  }

  if( getWellStatus() == WellControls::Status::CLOSED )
  {
    GEOS_LOG( GEOS_FMT( "{}: well is shut", wellControlsName ) );
    if( outputFile.is_open())
    {
      // print all zeros in the rates file
      if( hasMinimumWHPConstraint() || hasMaximumWHPConstraint() )
        outputFile << ",0.0";
      outputFile << ",0.0,0.0,0.0";
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
    subRegion.getField< well::connectionRate >();

  integer const useSurfaceCond =  useSurfaceConditions();

  real64 const & currentBHP =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  real64 const & currentWHP =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentWHPString() );
  arrayView1d< real64 const > const & currentPhaseVolRate =
    getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  real64 const & currentTotalVolRate =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

  integer const hasWHP =  hasMinimumWHPConstraint()||hasMaximumWHPConstraint();
  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [&numPhase,
                              &numComp,
                              &useSurfaceCond,
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
    string const conditionKey = useSurfaceCond ? "surface" : "reservoir";
    string const unitKey = useSurfaceCond ? "s" : "r";

    real64 const currentTotalRate = connRate[iwelemRef];
    GEOS_LOG( GEOS_FMT( "{}: BHP (at the specified reference elevation): {} Pa",
                        wellControlsName, currentBHP ) );
    if( hasWHP )
      GEOS_LOG( GEOS_FMT( "{}: WHP (at the well head): {} Pa",
                          wellControlsName, currentWHP ) );

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

}

bool CompositionalMultiphaseWell::solveMinWHPConstraint( real64 const & time_n,
                                                         real64 const & dt,
                                                         integer const cycleNumber,
                                                         integer const coupledIterationNumber,
                                                         DomainPartition & domain,
                                                         MeshLevel & mesh,
                                                         ElementRegionManager & elemManager,
                                                         WellElementSubRegion & subRegion )
{
  bool whpLimiting = false;

  MinimumWHPConstraint * whpConstraint = getMinWHPConstraint();
  if( whpConstraint == nullptr || !whpConstraint->isConstraintActive() )
    return whpLimiting;

  real64 & currentBHP =   getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  array1d< real64 > & currentPhaseVolRate =
    getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  real64 & currentTotalVolRate =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

  real64 currentBHP_local = currentBHP;
  array1d< real64 > currentPhaseVolRate_local = currentPhaseVolRate;
  real64 currentTotalVolRate_local = currentTotalVolRate;
  // Turn off BHP for WHP constraint if active, will be reset if WHP is limiting
  MinimumBHPConstraint * bhpConstraint=  getMinimumBHPConstraintForWHP();
  bhpConstraint->setConstraintActive( false );
  real64 constraintWHP = whpConstraint->getConstraintValue( time_n );
  real64 currentWHP = constraintWHP;
  integer owner = -1;

  // Get the flow table function
  FunctionManager & functionManager = FunctionManager::getInstance();
  const ProdPipeFlowTableFunction & m_flowTable =  functionManager.getGroup< ProdPipeFlowTableFunction const >( whpConstraint->getFlowTableName());
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
  m_flowTable.calculateWHP( getName(), currentBHP, currentPhaseVolRate, currentWHP, flowTableSolveState );
  getReference< real64 >( viewKeyStruct::currentWHPString() ) = currentWHP;
  //getMinBHPConstraint()->setConstraintActive( true );

  std::cout << getName() << " " <<coupledIterationNumber <<  " WHP constraint check " << whpLimiting << " " << currentWHP << " < " << constraintWHP << " bhp " << currentBHP <<
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
      std::cout << time_n << " "  << getName()<<  " VFPBracketSolve  negative dP/dQ table " << dP_dQ_table << " " << currentWHP << " " <<
        " " << ql0 <<  " " << ql1 << " " << bhp0 << " " << bhp1  << std::endl;
      ProductionConstraint< LiquidRateConstraint > *  liqConstraint=  getMaxLiquidConstraintForWHP();
      setCurrentConstraint( liqConstraint );
      liqConstraint->setConstraintActive( true );

      setControl( static_cast< WellControls::Control >(liqConstraint->getControl()) );        // tjb old
      WellControls::Control wellControl = getControl();
      MpiWrapper::broadcast( wellControl, owner );
      setControl( wellControl );

      // lower bracker IPR solve
      liqConstraint->setConstraintValue( -ql0 );
      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      real64 iprBHP0 = currentBHP;
      std::cout << time_n << " IPRBracketSolve0 " << getName() << " whp " << currentWHP << " BHP " << " " << currentBHP << " " << liqConstraint->bottomHolePressure() << " " << -ql0 <<
        " " << liqConstraint->phaseVolumeRates() << " " << liqConstraint->totalVolumeRate() << " " <<
        liqConstraint->massRate() << std::endl;
      // upper bracket IPR solve
      liqConstraint->setConstraintValue( -ql1 );
      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      real64 iprBHP1 = currentBHP;
      std::cout << time_n << " IPRBracketSolve1 " << getName() << " whp " << currentWHP << " BHP  "  << currentBHP << " " << liqConstraint->bottomHolePressure() << " " << -ql1 << " " <<
        liqConstraint->phaseVolumeRates() << " " << liqConstraint->totalVolumeRate() << " " <<
        liqConstraint->massRate() << std::endl;
      liqConstraint->setConstraintActive( false );
      real64 dP_dQ_ipr = ( iprBHP1 - iprBHP0 ) / ( -ql1 - (-ql0) );
      if( dP_dQ_ipr > dP_dQ_table )
      {
        std::cout << time_n << " " << getName() << " WHP 0 constraint stability dP/dQ table " << dP_dQ_table << " dP/dQ ipr " << dP_dQ_ipr  << " " << currentWHP << " " <<
          (currentWHP < constraintWHP) << std::endl;
        dP_dQ_table = dP_dQ_ipr;
        whpLimiting = currentWHP < constraintWHP;
      }
      else
      {
        // set so well operates at minwhp
        currentWHP = constraintWHP;
        std::cout << time_n << " " << getName() << " WHP 1 constraint stability dP/dQ table " << dP_dQ_table << " dP/dQ ipr " << dP_dQ_ipr  << " " << currentWHP << " " <<
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
      std::cout << time_n << " " << getName() << " WHP 2 constraint stability dP/dQ table " << dP_dQ_table << " " << currentWHP << " " << (currentWHP < constraintWHP)  << std::endl;


    }

  }
  else
  {
    // no stab check
    whpLimiting = currentWHP < constraintWHP;
  }

  if( whpLimiting )
  {

    std::cout << getName() << " WHP constraint violated " << currentWHP << " < " << constraintWHP << " bhp " << currentBHP << " phase rates " << currentPhaseVolRate << " total vol " <<
      currentTotalVolRate <<
      std::endl;
    // WHP is limiting  set WHP to constraint value
    currentWHP = constraintWHP;



    // sets. tjb cleanup
    ProductionConstraint< LiquidRateConstraint > *  liqConstraint=  getMaxLiquidConstraintForWHP();
    setCurrentConstraint( liqConstraint );
    liqConstraint->setConstraintActive( true );
    setControl( static_cast< WellControls::Control >(liqConstraint->getControl()) );         // tjb old
    WellControls::Control wellControl = getControl();
    MpiWrapper::broadcast( wellControl, owner );
    setControl( wellControl );
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

      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      bhp0 = currentBHP;
      flowTableSolveState=1;
      m_flowTable.calculateBHP( currentPhaseVolRate, currentWHP, tableBHP0,
                                flowTableSolveState );
      std::cout << getName() << " Solve at liquid rate [0] " << liqRate0 << " whp " << constraintWHP << " bhp " << bhp0 << " table bhp " << tableBHP0<< " phase rates " <<
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
      cSolve = m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                                        dt,
                                                        cycleNumber,
                                                        domain,
                                                        mesh,
                                                        elemManager,
                                                        subRegion );

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

    std::cout << time_n <<getName() << "IPRSolve at liquid rate [index] " << currentRateIndex << " liqr " << liqRate1 << " whp  " << constraintWHP << " bhp " << bhp1 << " table bhp " <<
      tableBHP1<< " phase rates " <<
      currentPhaseVolRate <<
      " total vol " << currentTotalVolRate << " bhp res " << bhp1-tableBHP1 << std::endl;

    std::cout << getName() << " Solve at liq brackets found " <<  std::endl;
    setCurrentConstraint( bhpConstraint );
    setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );

    integer const maxIters=100;
    real64 const tol = 1;
    integer iter = 0;

    of.open( "fls_"+getName() + "_"+std::to_string( time_n )+"_"+std::to_string( dt )+"_"+std::to_string( coupledIterationNumber )+".csv" );
    of << " error,wellbhp,tablebhp,whp,orate,grate,wrate "  << std::endl;
    bhpConstraint->setConstraintActive( true );
    while( iter < maxIters && std::abs( tableBHP1 - bhp1 )  > tol )
    {
      // update whp
      bhp1=bhp1+0.50*(tableBHP1-bhp1);

      bhpConstraint->setConstraintValue( bhp1 );
      std::cout << "SolveLinear system " << iter << std::endl;
      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      //bhp1 = currentBHP;   // current(s) were updated in solveNonlinearSystem

      m_flowTable.calculateBHP( currentPhaseVolRate, currentWHP, tableBHP1,
                                flowTableSolveState );
      of << std::abs( bhp1-tableBHP1 ) << ","<<bhp1 << "," << tableBHP1<< ","<< currentWHP <<  "," << currentPhaseVolRate[0] << "," << currentPhaseVolRate[1] << " ," <<
        currentPhaseVolRate[2] << std::endl;
      std::cout <<  time_n << " " << getName() << "IPRVFPSolve " << iter<< " whp  " <<  " whp  " << constraintWHP << " bhp " << bhp1 << " table bhp " << tableBHP1 << " phase rates " <<
        currentPhaseVolRate <<
        " total vol " << currentTotalVolRate << " bhp res " << bhp1-tableBHP1 << " " << (tableBHP1 - bhp1 )  << std::endl;

      bhpConstraint->setConstraintValue( bhp1 );

      ++iter;
    }
    of.close();
    if( false )
    {
      getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentWHPString() ) = currentWHP;
      // Use liquid constraint to find intersection of IPR and VLP
      setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );                       // tjb old
      setCurrentConstraint( bhpConstraint );
      setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );     // tjb old
      setCurrentConstraint( bhpConstraint );
      getMinBHPConstraint()->setConstraintActive( false );
      bhpConstraint->setBHP ( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      bhpConstraint->setPhaseVolumeRates ( getReference< array1d< real64 > >(
                                             CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      bhpConstraint->setTotalVolumeRate ( getReference< real64 >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      bhpConstraint->setMassRate( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));

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
      setCurrentConstraint( liqConstraint );
      setControl( static_cast< WellControls::Control >(liqConstraint->getControl()) );       // tjb old
      liqConstraint->setBHP ( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      liqConstraint->setPhaseVolumeRates ( getReference< array1d< real64 > >(
                                             CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      liqConstraint->setTotalVolumeRate ( getReference< real64 >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      liqConstraint->setMassRate( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
      std::cout << getName() << " WHPConstraint " << liqConstraint->getName() << " limiting " << whpLimiting << " " << liqConstraint->bottomHolePressure() << " " <<
        liqConstraint->phaseVolumeRates() << " " << liqConstraint->totalVolumeRate() << " " <<
        liqConstraint->massRate() <<
        std::endl;

    }

  }
  return whpLimiting;
}
bool CompositionalMultiphaseWell::solveMaxWHPConstraint( real64 const & time_n,
                                                         real64 const & dt,
                                                         integer const cycleNumber,
                                                         integer const coupledIterationNumber,
                                                         DomainPartition & domain,
                                                         MeshLevel & mesh,
                                                         ElementRegionManager & elemManager,
                                                         WellElementSubRegion & subRegion )
{
  bool whpLimiting = false;

  MaximumWHPConstraint * whpConstraint = getMaxWHPConstraint();
  if( whpConstraint == nullptr || !whpConstraint->isConstraintActive() )
    return whpLimiting;

  real64 & currentBHP =   getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  array1d< real64 > & currentPhaseVolRate =
    getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  real64 & currentTotalVolRate =
    getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );

  real64 currentBHP_local = currentBHP;
  array1d< real64 > currentPhaseVolRate_local = currentPhaseVolRate;
  real64 currentTotalVolRate_local = currentTotalVolRate;
  std::cout << getName() << " Max WHP constraint check " << whpLimiting << " " << currentBHP << " phase rates " << currentPhaseVolRate <<
    " total vol " << currentTotalVolRate <<
    std::endl;
  // Turn off BHP for WHP constraint if active, will be reset if WHP is limiting
  MaximumBHPConstraint * bhpConstraint=  getMaximumBHPConstraintForWHP();
  bhpConstraint->setConstraintActive( false );
  real64 constraintWHP = whpConstraint->getConstraintValue( time_n );
  real64 currentWHP = constraintWHP;
  integer owner = -1;

  // Get the flow table function
  FunctionManager & functionManager = FunctionManager::getInstance();
  const InjPipeFlowTableFunction & m_flowTable =  functionManager.getGroup< InjPipeFlowTableFunction const >( whpConstraint->getFlowTableName());
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
  m_flowTable.calculateWHP( getName(), currentBHP, currentTotalVolRate, currentWHP, flowTableSolveState );
  getReference< real64 >( viewKeyStruct::currentWHPString() ) = currentWHP;
  //getMinBHPConstraint()->setConstraintActive( true );

  std::cout << getName() << " " <<coupledIterationNumber <<  " WHP constraint check " << whpLimiting << " " << currentWHP << " < " << constraintWHP << " bhp " << currentBHP <<
    " phase rates " << currentTotalVolRate <<
    " total vol " << currentTotalVolRate <<
    " whpLimiting " << whpLimiting << std::endl;
  //whpLimiting=false;
  // check stability
  bool stabCheck = false;
  if( stabCheck )
  {

    real64 ql0, ql1, bhp0, bhp1;
    real64 dP_dQ_table = m_flowTable.calculatedPdQ( currentTotalVolRate, currentWHP, ql0, ql1, bhp0, bhp1 );

    if( dP_dQ_table < 0.0 )
    {
      std::cout << time_n << " "  << getName()<<  " VFPBracketSolve  negative dP/dQ table " << dP_dQ_table << " " << currentWHP << " " <<
        " " << ql0 <<  " " << ql1 << " " << bhp0 << " " << bhp1  << std::endl;
      InjectionConstraint< PhaseVolumeRateConstraint > *  volConstraint=  getMaxPhaseVolumeConstraintForWHP();
      setCurrentConstraint( volConstraint );
      volConstraint->setConstraintActive( true );

      setControl( static_cast< WellControls::Control >(volConstraint->getControl()) );        // tjb old
      WellControls::Control wellControl = getControl();
      MpiWrapper::broadcast( wellControl, owner );
      setControl( wellControl );

      // lower bracker IPR solve
      volConstraint->setConstraintValue( -ql0 );
      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      real64 iprBHP0 = currentBHP;
      std::cout << time_n << " IPRBracketSolve0 " << getName() << " whp " << currentWHP << " BHP " << " " << currentBHP << " " << volConstraint->bottomHolePressure() << " " << -ql0 <<
        " " << volConstraint->phaseVolumeRates() << " " << volConstraint->totalVolumeRate() << " " <<
        volConstraint->massRate() << std::endl;
      // upper bracket IPR solve
      volConstraint->setConstraintValue( -ql1 );
      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      real64 iprBHP1 = currentBHP;
      std::cout << time_n << " IPRBracketSolve1 " << getName() << " whp " << currentWHP << " BHP  "  << currentBHP << " " << volConstraint->bottomHolePressure() << " " << -ql1 << " " <<
        volConstraint->phaseVolumeRates() << " " << volConstraint->totalVolumeRate() << " " <<
        volConstraint->massRate() << std::endl;
      volConstraint->setConstraintActive( false );
      real64 dP_dQ_ipr = ( iprBHP1 - iprBHP0 ) / ( -ql1 - (-ql0) );
      if( dP_dQ_ipr > dP_dQ_table )
      {
        std::cout << time_n << " " << getName() << " WHP 0 constraint stability dP/dQ table " << dP_dQ_table << " dP/dQ ipr " << dP_dQ_ipr  << " " << currentWHP << " " <<
          (currentWHP < constraintWHP) << std::endl;
        dP_dQ_table = dP_dQ_ipr;
        whpLimiting = currentWHP < constraintWHP;
      }
      else
      {
        // set so well operates at minwhp
        currentWHP = constraintWHP;
        std::cout << time_n << " " << getName() << " WHP 1 constraint stability dP/dQ table " << dP_dQ_table << " dP/dQ ipr " << dP_dQ_ipr  << " " << currentWHP << " " <<
          (currentWHP > constraintWHP)  << std::endl;

        whpLimiting = true;
      }
      currentBHP = currentBHP_local;
      currentPhaseVolRate = currentPhaseVolRate_local;
      currentTotalVolRate = currentTotalVolRate_local;
    }
    else
    {
      // currentWHP is stable value
      whpLimiting = currentWHP > constraintWHP;
      std::cout << time_n << " " << getName() << " WHP 2 constraint stability dP/dQ table " << dP_dQ_table << " " << currentWHP << " " << (currentWHP < constraintWHP)  << std::endl;


    }

  }
  else
  {
    // no stab check
    whpLimiting = currentWHP > constraintWHP;
  }

  if( whpLimiting )
  {

    std::cout << getName() << " WHP constraint violated " << currentWHP << " < " << constraintWHP << " bhp " << currentBHP << " phase rates " << currentPhaseVolRate << " total vol " <<
      currentTotalVolRate <<
      std::endl;
    // WHP is limiting  set WHP to constraint value
    currentWHP = constraintWHP;



    // sets. tjb cleanup
    InjectionConstraint< PhaseVolumeRateConstraint > *  volConstraint=  getMaxPhaseVolumeConstraintForWHP();
    setCurrentConstraint( volConstraint );
    volConstraint->setConstraintActive( true );
    setControl( static_cast< WellControls::Control >(volConstraint->getControl()) );         // tjb old
    WellControls::Control wellControl = getControl();
    MpiWrapper::broadcast( wellControl, owner );
    setControl( wellControl );
    std::ofstream of;
    of.open( "fl.csv" );
    of << "liq ,bhp ,tablebhp"<< std::endl;
    // Liquid constraint is used to find intersection of IPR and VLP
    const array1d< real64 > & tableRates = m_flowTable.getRates();
    std::cout << tableRates << std::endl;
    integer numRates = tableRates.size();
    real64 liqRate0;
    real64 bhp0;
    real64 tableBHP0;

    liqRate0= tableRates[0];
    for( integer i=0; i < 1; ++i )    // numRates; ++i )
    {
      of << tableRates[i] << ",";
      liqRate0 = tableRates[i];
      volConstraint->setConstraintValue( liqRate0 );

      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      bhp0 = currentBHP;
      flowTableSolveState=1;
      m_flowTable.calculateBHP( currentTotalVolRate, currentWHP, tableBHP0,
                                flowTableSolveState );
      std::cout << getName() << " Solve at liquid rate [0] " << liqRate0 << " whp " << constraintWHP << " bhp " << bhp0 << " table bhp " << tableBHP0<< " phase rates " <<
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
      volConstraint->setConstraintValue( tableRates[currentRateIndex] );
      cSolve = m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                                        dt,
                                                        cycleNumber,
                                                        domain,
                                                        mesh,
                                                        elemManager,
                                                        subRegion );

    }
    if( !cSolve )
    {
      throw("ft solve ");
    }
    real64 bhp1 = currentBHP;
    real64 liqRate1 = tableRates[currentRateIndex];
    real64 tableBHP1;
    m_flowTable.calculateBHP( currentTotalVolRate, currentWHP, tableBHP1,
                              flowTableSolveState );

    std::cout << time_n <<getName() << " IPRSolve at liquid rate [index] " << currentRateIndex << " liqr " << liqRate1 << " whp  " << constraintWHP << " bhp " << bhp1 << " table bhp " <<
      tableBHP1<< " phase rates " <<
      currentPhaseVolRate <<
      " total vol " << currentTotalVolRate << " bhp res " << bhp1-tableBHP1 << std::endl;

    std::cout << getName() << " Solve at liq brackets found " <<  std::endl;
    setCurrentConstraint( bhpConstraint );
    setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );

    integer const maxIters=100;
    real64 const tol = 1;
    integer iter = 0;

    of.open( "fls_"+getName() + "_"+std::to_string( time_n )+"_"+std::to_string( dt )+"_"+std::to_string( coupledIterationNumber )+".csv" );
    of << " error,wellbhp,tablebhp,whp,orate,grate,wrate "  << std::endl;
    bhpConstraint->setConstraintActive( true );
    while( iter < maxIters && std::abs( tableBHP1 - bhp1 )  > tol )
    {
      // update whp
      bhp1=bhp1+0.50*(tableBHP1-bhp1);

      bhpConstraint->setConstraintValue( bhp1 );
      std::cout << "SolveLinear system " << iter << std::endl;
      m_wellNewtonSolver.solveNonlinearSystem( *this, time_n,
                                               dt,
                                               cycleNumber,
                                               domain,
                                               mesh,
                                               elemManager,
                                               subRegion );
      //bhp1 = currentBHP;   // current(s) were updated in solveNonlinearSystem

      m_flowTable.calculateBHP( currentTotalVolRate, currentWHP, tableBHP1,
                                flowTableSolveState );
      of << std::abs( bhp1-tableBHP1 ) << ","<<bhp1 << "," << tableBHP1<< ","<< currentWHP <<  "," <<currentTotalVolRate << std::endl;
      std::cout <<  time_n << " " << getName() << "IPRVFPSolve " << iter<< " whp  " <<  " whp  " << constraintWHP << " bhp " << bhp1 << " table bhp " << tableBHP1 << " phase rates " <<

        " total vol " << currentTotalVolRate << " bhp res " << bhp1-tableBHP1 << " " << (tableBHP1 - bhp1 )  << std::endl;

      bhpConstraint->setConstraintValue( bhp1 );

      ++iter;
    }
    of.close();
    if( false )
    {
      getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentWHPString() ) = currentWHP;
      // Use liquid constraint to find intersection of IPR and VLP
      setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );                       // tjb old
      setCurrentConstraint( bhpConstraint );
      setControl( static_cast< WellControls::Control >(bhpConstraint->getControl()) );     // tjb old
      setCurrentConstraint( bhpConstraint );
      getMinBHPConstraint()->setConstraintActive( false );
      bhpConstraint->setBHP ( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      bhpConstraint->setPhaseVolumeRates ( getReference< array1d< real64 > >(
                                             CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      bhpConstraint->setTotalVolumeRate ( getReference< real64 >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      bhpConstraint->setMassRate( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));

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
      volConstraint->setConstraintValue( currentTotalVolRate );
      volConstraint->setConstraintActive( true );
      setCurrentConstraint( volConstraint );
      setControl( static_cast< WellControls::Control >(volConstraint->getControl()) );       // tjb old
      volConstraint->setBHP ( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() ));
      volConstraint->setPhaseVolumeRates ( getReference< array1d< real64 > >(
                                             CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() ) );
      volConstraint->setTotalVolumeRate ( getReference< real64 >(
                                            CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() ));
      volConstraint->setMassRate( getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentMassRateString() ));
      std::cout << getName() << " WHPConstraint " << volConstraint->getName() << " limiting " << whpLimiting << " " << volConstraint->bottomHolePressure() << " " <<
        volConstraint->phaseVolumeRates() << " " << volConstraint->totalVolumeRate() << " " <<
        volConstraint->massRate() <<
        std::endl;

    }

  }
  return whpLimiting;
}

void CompositionalMultiphaseWell::outputSingleWellDebug( real64 const time,
                                                         real64 const dt,
                                                         integer current_newton_iteration,
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

  integer num_timestep_cuts  =0;
  if( m_writeSegDebug > 1 )
  {
    // CompositionalMultiphaseBase const & flowSolver = .getParent().getParent().getGroup< CompositionalMultiphaseBase >(
    // getFlowSolverName() );
    // auto solver_names = getParent().getSubGroupsNames();
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



      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      PerforationData & perforationData = *subRegion.getPerforationData();
      using CompFlowAccessors =
        StencilAccessors< fields::flow::pressure,
                          fields::flow::temperature,
                          fields::flow::phaseVolumeFraction,



                          fields::flow::dPhaseVolumeFraction,
                          fields::flow::globalCompDensity,
                          fields::flow::dGlobalCompFraction_dGlobalCompDensity >;

      CompFlowAccessors compFlowAccessors( mesh.getElemManager(), getFlowSolverName() );

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
      MultiFluidAccessors multiFluidAccessors( mesh.getElemManager(), getFlowSolverName() );

      using RelPermAccessors =
        StencilMaterialAccessors< RelativePermeabilityBase,
                                  fields::relperm::phaseRelPerm,
                                  fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;

      RelPermAccessors relPermAccessors( mesh.getElemManager(), getFlowSolverName() );


      string const srn = subRegion.getName();

      std::vector< string > cp_der {"dP", "dT"};
      for( integer i=0; i<m_numComponents; i++ )
      {
        cp_der.push_back( "dRho"+std::to_string( i+1 ));
      }
      if( !m_wellPropWriter[srn].initialized() )
      {
        integer my_rank = MpiWrapper::commRank( MPI_COMM_GEOS );
        m_wellPropWriter[srn].initialize_perf( my_rank, m_ratesOutputDir, getName(), perforationData );
        m_wellPropWriter[srn].initialize_seg( my_rank, m_ratesOutputDir, getName(), fluid.phaseNames(), fluid.componentNames(), subRegion );
        m_wellDebugInit=true;
      }
      m_wellPropWriter[srn].registerSeg2dProp( {"X", "Y", "Z"}, subRegion.getElementCenter());
      m_wellPropWriter[srn].registerSegProp( "Pressure", subRegion.getField< fields::well::pressure >());
      if( isThermal() )
      {
        m_wellPropWriter[srn].registerSegProp( "Temperature", subRegion.getField< fields::well::temperature >());
      }
      m_wellPropWriter[srn].registerSegComponentProp( "ComponentDensity", subRegion.getField< fields::well::globalCompDensity >());

      m_wellPropWriter[srn].registerSegProp( "TotalRate", subRegion.getField< fields::well::connectionRate >());

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

      m_wellPropWriter[srn].write( m_numTimesteps, dt, cycle, subevent, m_numTimesteps,
                                   current_newton_iteration,
                                   num_timestep_cuts );
    }

  }


}



}   // namespace geos
