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
 * @file ImmiscibleMultiphaseFlow.cpp
 */

#include "ImmiscibleMultiphaseFlow.hpp"

#include "FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleTrustRegionKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/RelativePermeabilityUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CapillaryPressureUpdateKernel.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureFields.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"

#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "physicsSolvers/fluidFlow/SourceFluxStatistics.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"

#include <cmath>

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields::immiscibleMultiphaseFlow;
using namespace immiscibleMultiphaseKernels;
using namespace immiscibleFlowUtilities;


ImmiscibleMultiphaseFlow::ImmiscibleMultiphaseFlow( const string & name,
                                                    Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 2 ),
  m_hasCapPressure( 0 ),
  m_fluxInflection ( 0 )  
{
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Temperature" );

  this->registerWrapper( viewKeyStruct::useTotalMassEquationString(), &m_useTotalMassEquation ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag indicating whether total mass equation is used" );

  this->registerWrapper( viewKeyStruct::gravityDensitySchemeString(), &m_gravityDensityScheme ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( GravityDensityScheme::ArithmeticAverage ).
    setDescription( "Scheme for density treatment in gravity" );

  this->registerWrapper( viewKeyStruct::solutionChangeScalingFactorString(), &m_solutionChangeScalingFactor ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Damping factor for solution change targets" );

  this->registerWrapper( viewKeyStruct::targetRelativePresChangeString(), &m_targetRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.2 ).
    setDescription( "Target (relative) change in pressure in a time step (expected value between 0 and 1)" );
    
  this->registerWrapper( viewKeyStruct::targetPhaseVolFracChangeString(), &m_targetPhaseVolFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.2 ).
    setDescription( "Target (absolute) change in phase volume fraction in a time step" );

  this->registerWrapper( viewKeyStruct::scalingTypeString(), &m_scalingType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( ScalingType::Local ).
    setDescription( "Solution scaling type."
                    "Valid options:\n* " + EnumStrings< ScalingType >::concat( "\n* " ) );

  this->registerWrapper( viewKeyStruct::scalingFactorTypeString(), &m_scalingFactorType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( ScalingFactorType::MaxVariation ).
    setDescription( "Scaling factor type."
                    "Valid options:\n* " + EnumStrings< ScalingFactorType >::concat( "\n* " ) );

  this->registerWrapper( viewKeyStruct::minScalingFactorString(), &m_minScalingFactor ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "Minimum value for solution scaling factor" );

  this->registerWrapper( viewKeyStruct::maxAbsolutePresChangeString(), &m_maxAbsolutePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).       // disabled by default
    setDescription( "Maximum (absolute) pressure change in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString(), &m_maxRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).
    setDescription( "Maximum (relative) change in pressure in a Newton iteration" );  
  
  this->registerWrapper( viewKeyStruct::maxAbsoluteSatChangeString(), &m_maxAbsoluteSatChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.2 ).       // << Appleyard chopping limit >>
    setDescription( "Maximum (absolute) saturation change in a Newton iteration" );

  this->registerWrapper( viewKeyStruct::maxRelativeSatChangeString(), &m_maxRelativeSatChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1.0 ).     // disabled by default
    setDescription( "Maximum (relative) change in saturation in a Newton iteration" );  

  this->registerWrapper( viewKeyStruct::allowOutOfBoundPressureString(), &m_allowOutOfBoundPressure ).    
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ). // negative pressure is not allowed by default
    setDescription( "Flag indicating if negative pressure is allowed" );

  this->registerWrapper( viewKeyStruct::allowLocalPresChoppingString(), &m_allowPresChopping ).    
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ). 
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative pressure is allowed" );

  this->registerWrapper( viewKeyStruct::allowOutOfBoundSatString(), &m_allowOutOfBoundSaturation ).    
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ). // out of bound saturation is allowed by default
    setDescription( "Flag indicating if out of bound saturation is allowed" );

  this->registerWrapper( viewKeyStruct::allowLocalSatChoppingString(), &m_allowSatChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ). // local saturation chopping is allowed by default
    setDescription( "Flag indicating whether local (cell-wise) chopping of out of bound saturation is allowed" );

  this->registerWrapper( viewKeyStruct::trustRegionMinNewtonIterString(), &m_trustRegionParams.minNewtonIter ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Minimum number of Newton iterations taken before applying trust region" );

  this->registerWrapper( viewKeyStruct::trustRegionMinGradientString(), &m_trustRegionParams.minGradient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Minimum gradient necessary for applying trust region" );

  this->registerWrapper( viewKeyStruct::trustRegionMaxIterString(), &m_trustRegionParams.maxIter ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 5 ).
    setDescription( "Maximum number of trust region iterations" );

  this->registerWrapper( viewKeyStruct::trustRegionMinPotentialDiffString(), &m_trustRegionParams.dPhiMin ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Minimum potential difference to apply a damping factor" );

  this->registerWrapper( viewKeyStruct::trustRegionMinDerivativeString(), &m_trustRegionParams.d2RMin ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "Minimum directional second derivative to apply a damping factor" );

  this->registerWrapper( viewKeyStruct::trustRegionMinKinkFactorString(), &m_trustRegionParams.minKinkFactor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "Minimum discontinuity damping factor" );

  this->registerWrapper( viewKeyStruct::trustRegionMinInfFactorString(), &m_trustRegionParams.minInfFactor ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "Minimum inflection damping factor" );

  this->registerWrapper( viewKeyStruct::trustRegionKinkFactorDeltaString(), &m_trustRegionParams.kinkFactorDelta ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.05 ).
    setDescription( "Stretching factor applied to damping factor to allow for a small crossing of discontinuities" );

  this->registerWrapper( viewKeyStruct::trustRegionRelResThresString(), &m_trustRegionParams.relResThres ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Minimum relative residual threshold for applying damping factor" );

  this->registerWrapper( viewKeyStruct::trustRegionAbsResThresString(), &m_trustRegionParams.absResThres ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e2 ). // 1e-1
    setDescription( "Minimum absolute residual threshold for applying damping factor" );

  this->registerWrapper( viewKeyStruct::trustRegionUseAccumString(), &m_trustRegionParams.useAccum ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag for inclusion of accumulation term" );
  
}

void ImmiscibleMultiphaseFlow::postInputInitialization()
{
  FlowSolverBase::postInputInitialization();

  GEOS_ERROR_IF_GT_MSG( m_maxAbsoluteSatChange, 1.0,
                        getWrapperDataContext( viewKeyStruct::maxAbsoluteSatChangeString() ) <<
                        ": The maximum absolute change in component saturation in a Newton iteration must be smaller than or equal to 1.0" ); 
  GEOS_ERROR_IF_LE_MSG( m_targetRelativePresChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::targetRelativePresChangeString() ) <<
                        ": The target relative change in pressure in a time step must be larger than 0.0" );  
  GEOS_ERROR_IF_LE_MSG( m_targetPhaseVolFracChange, 0.0,
                        getWrapperDataContext( viewKeyStruct::targetPhaseVolFracChangeString() ) <<
                        ": The target change in phase volume fraction in a time step must be larger than or equal to 0.0" );
  GEOS_ERROR_IF_LT_MSG( m_solutionChangeScalingFactor, 0.0,
                        getWrapperDataContext( viewKeyStruct::solutionChangeScalingFactorString() ) <<
                        ": The solution change scaling factor must be larger or equal to 0.0" );
  GEOS_ERROR_IF_GT_MSG( m_solutionChangeScalingFactor, 1.0,
                        getWrapperDataContext( viewKeyStruct::solutionChangeScalingFactorString() ) <<
                        ": The solution change scaling factor must be smaller or equal to 1.0" );
  GEOS_ERROR_IF_LE_MSG( m_minScalingFactor, 0.0,
                        getWrapperDataContext( viewKeyStruct::minScalingFactorString() ) <<
                        ": The minumum scaling factor must be larger than 0.0" );
  GEOS_ERROR_IF_GT_MSG( m_minScalingFactor, 1.0,
                        getWrapperDataContext( viewKeyStruct::minScalingFactorString() ) <<
                        ": The minumum scaling factor must be smaller than or equal to 1.0" );
}

void ImmiscibleMultiphaseFlow::registerDataOnMesh( Group & meshBodies )
{  
  FlowSolverBase::registerDataOnMesh( meshBodies );

  // 0. Find a "reference" fluid model name (at this point, models are already attached to subregions)
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // If at least one region has a capillary pressure model, consider it enabled for all
      string const capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
      if( !capPresName.empty() )
      {
        m_hasCapPressure = true;
      }
    } );
  } );

  m_numDofPerCell = m_numPhases;

  // 2. Register and resize all fields as necessary
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_hasCapPressure )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::capPressureNamesString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRestartFlags( RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 ).
          setDescription( "Name of the capillary pressure constitutive model to use" ).
          reference();

        string & capPresName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
        GEOS_THROW_IF( capPresName.empty(),
                       GEOS_FMT( "{}: Capillary pressure model not found on subregion {}",
                                 getDataContext(), subRegion.getDataContext() ),
                       InputError );
      }

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.
      subRegion.registerField< phaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< bcPhaseVolumeFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMass >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMass_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< phaseMobility >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< dPhaseMobility >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numPhases ); // dP, dS

      subRegion.registerField< solutionUpdate >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );

      subRegion.registerField< fields::flow::pressureScalingFactor >( getName() );
      subRegion.registerField< phaseVolumeFractionScalingFactor >( getName() );

    } );

  } );
}

void ImmiscibleMultiphaseFlow::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{

  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< TwoPhaseImmiscibleFluid >( subRegion );
  GEOS_ERROR_IF( fluidName.empty(), GEOS_FMT( "{}: Fluid model not found on subregion {}",
                                              getDataContext(), subRegion.getName() ) );



  string & relPermName = subRegion.registerWrapper< string >( viewKeyStruct::relPermNamesString() ).
                           setPlotLevel( PlotLevel::NOPLOT ).
                           setRestartFlags( RestartFlags::NO_WRITE ).
                           setSizedFromParent( 0 ).
                           setDescription( "Name of the relative permeability constitutive model to use" ).
                           reference();

  relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );

  GEOS_THROW_IF( relPermName.empty(),
                 GEOS_FMT( "{}: Relative permeability model not found on subregion {}",
                           getDataContext(), subRegion.getDataContext() ),
                 InputError );
}

void ImmiscibleMultiphaseFlow::initializePreSubGroups()
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::immiscibleMultiphaseFVM;

  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getField< fields::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );
}


void ImmiscibleMultiphaseFlow::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();

  TwoPhaseImmiscibleFluid & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    FluidUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(), fluidWrapper, pres );
  } );
}


void ImmiscibleMultiphaseFlow::updateRelPermModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;


  GEOS_UNUSED_VAR( dataGroup );

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
    dataGroup.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

  string const & relPermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase & relPerm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relPermName );

  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

    isothermalCompositionalMultiphaseBaseKernels::
      RelativePermeabilityUpdateKernel::
      launch< parallelDevicePolicy<> >( dataGroup.size(),
                                        relPermWrapper,
                                        phaseVolFrac );
  } );
}

void ImmiscibleMultiphaseFlow::updateCapPressureModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  if( m_hasCapPressure )
  {
    arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
      dataGroup.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

    string const & cappresName = dataGroup.getReference< string >( viewKeyStruct::capPressureNamesString() );
    CapillaryPressureBase & capPressure = getConstitutiveModel< CapillaryPressureBase >( dataGroup, cappresName );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      isothermalCompositionalMultiphaseBaseKernels::
        CapillaryPressureUpdateKernel::
        launch< parallelDevicePolicy<> >( dataGroup.size(),
                                          capPresWrapper,
                                          phaseVolFrac );
    } );
  }
}


void ImmiscibleMultiphaseFlow::updateFluidState( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  updateFluidModel( subRegion );
  updateVolumeConstraint( subRegion );
  updatePhaseMass( subRegion );
  updateRelPermModel( subRegion );
  updatePhaseMobility( subRegion );
  updateCapPressureModel( subRegion );
}


void ImmiscibleMultiphaseFlow::updatePhaseMass( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );

  TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );

  CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView2d< real64 const > const porosity = solid.getPorosity();

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac= subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView3d< real64 const, multifluid::USD_PHASE > phaseDens = fluid.phaseDensity();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseMass = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass >();

  // Might be needed for geomechanics????? if so, need to change the accumulation as well?
  //arrayView1d< real64 > const deltaVolume = subRegion.getField< fields::flow::deltaVolume >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const poreVolume = volume[ei] * porosity[ei][0];
    for( integer ip = 0; ip < 2; ++ip )
    {
      phaseMass[ei][ip] = poreVolume * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip];
    }
  } );
}


void ImmiscibleMultiphaseFlow::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  // note that the phase mobility computed here also includes phase density
  string const & fluidName = dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() );
  TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( dataGroup, fluidName );

  string const & relpermName = dataGroup.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, relpermName );

  immiscibleMultiphaseKernels::
    PhaseMobilityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                               dataGroup,
                                               fluid,
                                               relperm );
}

void ImmiscibleMultiphaseFlow::initializeFluidState( MeshLevel & mesh,
                                                     string_array const & regionNames )
{
  GEOS_MARK_FUNCTION;

  mesh.getElemManager().forElementSubRegions( regionNames,
                                              [&]( localIndex const,
                                                   ElementSubRegionBase & subRegion )
  {
    // 2. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    updateFluidModel( subRegion );

  } );

  // for some reason CUDA does not want the host_device lambda to be defined inside the generic lambda
  // I need the exact type of the subRegion for updateSolidflowProperties to work well.
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                              SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                           auto & subRegion )
  {
    // 4. Initialize/update dependent state quantities

    // 4.1 Update the constitutive models that only depend on
    //      - the primary variables
    //      - the fluid constitutive quantities (as they have already been updated)
    // We postpone the other constitutive models for now
    // In addition, to avoid multiplying permeability/porosity bay netToGross in the assembly kernel, we do it once and for all here
    arrayView1d< real64 const > const netToGross = subRegion.template getField< fields::flow::netToGross >();
    CoupledSolidBase const & porousSolid =
      getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
    PermeabilityBase const & permeabilityModel =
      getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() ) );
    permeabilityModel.scaleHorizontalPermeability( netToGross );
    porousSolid.scaleReferencePorosity( netToGross );
    saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
    updatePorosityAndPermeability( subRegion );

    // Now, we initialize and update each constitutive model one by one

    // 4.2 Save the computed porosity into the old porosity
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity
    string const & solidName = subRegion.template getReference< string >( viewKeyStruct::solidNamesString() );
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
    porousMaterial.initializeState();

    // 4.3 Initialize/update the relative permeability model using the initial phase volume fraction
    //     This is needed to handle relative permeability hysteresis
    //     Also, initialize the fluid model
    //
    // Note:
    // - This must be called after updateVolumeConstraint
    // - This step depends on phaseVolFraction

    // initialized phase volume fraction
    arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
      subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

    string const & relpermName = subRegion.template getReference< string >( viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase & relPermMaterial =
      getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
    relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac ); // this needs to happen before calling updateRelPermModel
    updateRelPermModel( subRegion );
    relPermMaterial.saveConvergedState(); // this needs to happen after calling updateRelPermModel

    // 4.4 Then, we initialize/update the capillary pressure model
    //
    // Note:
    // - This must be called after updatePorosityAndPermeability
    // - This step depends on porosity and permeability
    if( m_hasCapPressure )
    {
      // initialized porosity
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      string const & permName = subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeabilityMaterial = getConstitutiveModel< PermeabilityBase >( subRegion, permName );
      // initialized permeability
      arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

      string const & capPressureName = subRegion.template getReference< string >( viewKeyStruct::capPressureNamesString() );
      CapillaryPressureBase const & capPressureMaterial =
        getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressureName );
      capPressureMaterial.initializeRockState( porosity, permeability ); // this needs to happen before calling updateCapPressureModel
      updateCapPressureModel( subRegion );
    }

    // 4.5 Update the phase mobility
    //
    // Note:
    // - This must be called after updateRelPermModel
    // - This step depends phaseRelPerm
    updatePhaseMobility( subRegion );

  } );

  // 5. Save initial pressure
  mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 > const initPres = subRegion.getField< fields::flow::initialPressure >();
    arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();
    arrayView1d< real64 > const initTemp = subRegion.template getField< fields::flow::initialTemperature >();
    initPres.setValues< parallelDevicePolicy<> >( pres );
    initTemp.setValues< parallelDevicePolicy<> >( temp );

    // TODO: Missing updatePhaseMass?
  } );
}


void ImmiscibleMultiphaseFlow::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::pressure::key(),
                                       fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );
  } );

  initializeState( domain );
}


void
ImmiscibleMultiphaseFlow::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                             real64 const & GEOS_UNUSED_PARAM( dt ),
                                             DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      saveConvergedState( subRegion );

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateVolumeConstraint( subRegion );
      updateFluidState( subRegion );

      // after the update, save the new saturation
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
      phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );

      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const & phaseMass =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const & phaseMass_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();
      phaseMass_n.setValues< parallelDevicePolicy<> >( phaseMass );

    } );
  } );
}

void ImmiscibleMultiphaseFlow::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                               real64 const dt,
                                               DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationTerm( domain,
                            dofManager,
                            localMatrix,
                            localRhs );


  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}



void ImmiscibleMultiphaseFlow::assembleAccumulationTerm( DomainPartition & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );

      TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
        
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      immiscibleMultiphaseKernels::
        AccumulationKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                   dofManager.rankOffset(),
                                                   m_useTotalMassEquation,
                                                   dofKey,
                                                   subRegion,
                                                   fluid,
                                                   solid,
                                                   localMatrix,
                                                   localRhs );

    } );
  } );
}

void ImmiscibleMultiphaseFlow::assembleFluxTerms( real64 const dt,
                                                  DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & )
  {
    fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
      immiscibleMultiphaseKernels::
        FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                                             dofManager.rankOffset(),
                                                                             dofKey,
                                                                             m_hasCapPressure,
                                                                             m_useTotalMassEquation,
                                                                             m_gravityDensityScheme == GravityDensityScheme::PhasePresence,
                                                                             getName(),
                                                                             mesh.getElemManager(),
                                                                             stencilWrapper,
                                                                             dt,
                                                                             localMatrix.toViewConstSizes(),
                                                                             localRhs.toView() );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::setupDofs( DomainPartition const & domain,
                                          DofManager & dofManager ) const
{
  GEOS_UNUSED_VAR( domain, dofManager );
  // add a field for the cell-centered degrees of freedom
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  //// this call with instruct GEOS to reorder the dof numbers
  //dofManager.setLocalReorderingType( viewKeyStruct::elemDofFieldString(),
  //                                   DofManager::LocalReorderingType::ReverseCutHillMcKee );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}

void ImmiscibleMultiphaseFlow::applyBoundaryConditions( real64 const time_n,
                                                        real64 const dt,
                                                        DomainPartition & domain,
                                                        DofManager const & dofManager,
                                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                        arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // apply pressure boundary conditions.
  applyDirichletBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

  // apply flux boundary conditions
  applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
}


namespace
{
char const bcLogMessage[] =
  "ImmiscibleMultiphaseFlow {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";
}

bool ImmiscibleMultiphaseFlow::validateDirichletBC( DomainPartition & domain,
                                                    real64 const time ) const
{
  constexpr integer MAX_NP = 2;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  bool bcConsistent = true;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    // map: regionName -> subRegionName -> setName -> numPhases to check pressure/phase are present consistent
    map< string, map< string, map< string, ComponentMask< MAX_NP > > > > bcPresCompStatusMap;

    // 1. Check pressure Dirichlet BCs
    fsManager.apply< ElementSubRegionBase >( time,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&]( FieldSpecificationBase const &,
                                                  string const & setName,
                                                  SortedArrayView< localIndex const > const &,
                                                  ElementSubRegionBase & subRegion,
                                                  string const & )
    {
      // Check whether pressure has already been applied to this set
      string const & subRegionName = subRegion.getName();
      string const & regionName = subRegion.getParent().getParent().getName();

      auto & subRegionSetMap = bcPresCompStatusMap[regionName][subRegionName];
      if( subRegionSetMap.count( setName ) > 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( BCMessage::pressureConflict( regionName, subRegionName, setName,
                                                   fields::flow::pressure::key() ) );
      }
      subRegionSetMap[setName].setNumComp( m_numPhases );
    } );
    // 2. Check saturation Dirichlet BCs
    fsManager.apply< ElementSubRegionBase >( time,
                                             mesh,
                                             fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key(),
                                             [&] ( FieldSpecificationBase const & fs,
                                                   string const & setName,
                                                   SortedArrayView< localIndex const > const &,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      string const & subRegionName = subRegion.getName(   );
      string const & regionName = subRegion.getParent().getParent().getName();
      integer const comp = fs.getComponent();

      auto & subRegionSetMap = bcPresCompStatusMap[regionName][subRegionName];
      if( subRegionSetMap.count( setName ) == 0 )
      {
        bcConsistent = false;
        GEOS_WARNING( BCMessage::missingPressure( regionName, subRegionName, setName,
                                                  fields::flow::pressure::key() ) );
      }
      if( comp < 0 || comp >= m_numPhases )
      {
        bcConsistent = false;
        GEOS_WARNING( BCMessage::invalidComponentIndex( comp, fs.getName(),
                                                        fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() ) );
        return; // can't check next part with invalid component id
      }

      ComponentMask< MAX_NP > & compMask = subRegionSetMap[setName];
      if( compMask[comp] )
      {
        bcConsistent = false;
        fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & bc )
        {
          string_array const & componentNames = bc.getComponentNames();
          GEOS_WARNING( BCMessage::conflictingComposition( comp, componentNames[comp],
                                                           regionName, subRegionName, setName,
                                                           fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() ) );
        } );
      }
      compMask.set( comp );
    } );

    // 3.2 Check consistency between composition BC applied to sets
    // Note: for a temperature-only boundary condition, this loop does not do anything
    for( auto const & regionEntry : bcPresCompStatusMap )
    {
      for( auto const & subRegionEntry : regionEntry.second )
      {
        for( auto const & setEntry : subRegionEntry.second )
        {
          ComponentMask< MAX_NP > const & compMask = setEntry.second;

          fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & fs )
          {
            string_array const & componentNames = fs.getComponentNames();
            for( size_t ic = 0; ic < componentNames.size(); ic++ )
            {
              if( !compMask[ic] )
              {
                bcConsistent = false;
                GEOS_WARNING( BCMessage::notAppliedOnRegion( ic, componentNames[ic],
                                                             regionEntry.first, subRegionEntry.first, setEntry.first,
                                                             fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() ) );
              }
            }
          } );
        }
      }
    }
  } );

  return bcConsistent;
}

void ImmiscibleMultiphaseFlow::applyDirichletBC( real64 const time_n,
                                                 real64 const dt,
                                                 DofManager const & dofManager,
                                                 DomainPartition & domain,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateDirichletBC( domain, time_n + dt );
    GEOS_ERROR_IF( !bcConsistent, GEOS_FMT( "ImmiscibleMultiphaseFlow {}: inconsistent boundary conditions", getDataContext() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {

    // 1. Apply pressure Dirichlet BCs, store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::pressure::key(), fields::flow::bcPressure::key() );
    // 2. Apply saturation BC (phase volume fraction) and store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key(),
                                             fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() ); // fields::immiscibleMultiphaseFlow::bcPhaseVolumeFraction::key()

    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    // 3. Call constitutive update
    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {

      arrayView1d< real64 const > const bcPres =
        subRegion.getReference< array1d< real64 > >( fields::flow::bcPressure::key() );
      // arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const bcPhaseVolFraction =
      //   subRegion.getReference< array2d< real64, immiscibleFlow::LAYOUT_PHASE > >(
      //     fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() ); // fields::immiscibleMultiphaseFlow::bcPhaseVolumeFraction::key()

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< real64 const > const pres =
        subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFraction =
        subRegion.getReference< array2d< real64, immiscibleFlow::LAYOUT_PHASE > >(
          fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() );

      integer const numPhase = m_numPhases;


      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 3.1. Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    bcPres[ei],
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 3.2. For each phase, apply target saturation value
        for( integer ip = 0; ip < numPhase-1; ++ip )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ip + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      phaseVolFraction[ei][ip],   // bcPhaseVolFraction[ei][ip]
                                                      phaseVolFraction[ei][ip] );
          localRhs[localRow + ip + 1] = rhsValue;
        }
      } );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::applySourceFluxBC( real64 const time,
                                                  real64 const dt,
                                                  DofManager const & dofManager,
                                                  DomainPartition & domain,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Step 1: count individual source flux boundary conditions

  std::map< string, localIndex > bcNameToBcId;
  localIndex bcCounter = 0;

  fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&] ( SourceFluxBoundaryCondition const & bc )
  {
    // collect all the bc names to idx
    bcNameToBcId[bc.getName()] = bcCounter;
    bcCounter++;
  } );

  if( bcCounter == 0 )
  {
    return;
  }

  // Step 2: count the set size for each source flux (each source flux may have multiple target sets)

  array1d< globalIndex > bcAllSetsSize( bcNameToBcId.size() );

  computeSourceFluxSizeScalingFactor( time,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {

    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&]( SourceFluxBoundaryCondition const & fs,
                                                         string const & setName,
                                                         SortedArrayView< localIndex const > const & targetSet,
                                                         ElementSubRegionBase & subRegion,
                                                         string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
                                   getName(), time+dt, fs.getCatalogName(), fs.getName(),
                                   setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      if( targetSet.size() == 0 )
      {
        return;
      }
      if( !subRegion.hasWrapper( dofKey ) )
      {
        if( fs.getLogLevel() >= 1 )
        {
          GEOS_LOG_RANK( GEOS_FMT( "{}: trying to apply SourceFlux, but its targetSet named '{}' intersects with non-simulated region named '{}'.",
                                   getDataContext(), setName, subRegion.getName() ) );
        }
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      // Step 3.1: get the values of the source boundary condition that need to be added to the rhs

      array1d< globalIndex > dofArray( targetSet.size() );
      array1d< real64 > rhsContributionArray( targetSet.size() );
      arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
      localIndex const rankOffset = dofManager.rankOffset();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd( 0.0 );

      // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
      fs.computeRhsContribution< FieldSpecificationAdd,
                                 parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                           time + dt,
                                                           dt,
                                                           subRegion,
                                                           dofNumber,
                                                           rankOffset,
                                                           localMatrix,
                                                           dofArray.toView(),
                                                           rhsContributionArrayView,
                                                           [] GEOS_HOST_DEVICE ( localIndex const )
      {
        return 0.0;
      } );

      // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];
      integer const fluidPhaseId = fs.getComponent();
      integer const numFluidPhases = m_numPhases;
      integer useTotalMassEquation = m_useTotalMassEquation;
      forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                           targetSet,
                                                           rankOffset,
                                                           ghostRank,
                                                           fluidPhaseId,
                                                           numFluidPhases,
                                                           useTotalMassEquation,
                                                           dofNumber,
                                                           rhsContributionArrayView,
                                                           localRhs,
                                                           massProd] GEOS_HOST_DEVICE ( localIndex const a )
      {
        // we need to filter out ghosts here, because targetSet may contain them
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor here!
        massProd += rhsValue;
        if( useTotalMassEquation > 0 )
        {
          // for all "fluid components", we add the value to the total mass balance equation
          globalIndex const totalMassBalanceRow = dofNumber[ei] - rankOffset;
          localRhs[totalMassBalanceRow] += rhsValue;
          if( fluidPhaseId < numFluidPhases - 1 )
          {
            globalIndex const compMassBalanceRow = totalMassBalanceRow + fluidPhaseId + 1; // component mass bal equations are shifted
            localRhs[compMassBalanceRow] += rhsValue;
          }
        }
        else
        {
          globalIndex const compMassBalanceRow = dofNumber[ei] - rankOffset + fluidPhaseId;
          localRhs[compMassBalanceRow] += rhsValue;
        }
      } );

      SourceFluxStatsAggregator::forAllFluxStatWrappers( subRegion, fs.getName(),
                                                         [&]( SourceFluxStatsAggregator::WrappedStats & wrapper )
      {
        // set the new sub-region statistics for this timestep
        array1d< real64 > massProdArr{ m_numPhases };
        massProdArr[fluidPhaseId] = massProd.get();
        wrapper.gatherTimeStepStats( time, dt, massProdArr.toViewConst(), targetSet.size() );
      } );
    } );
  } );
}

real64 ImmiscibleMultiphaseFlow::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                        real64 const & GEOS_UNUSED_PARAM( dt ),
                                                        DomainPartition const & domain,
                                                        DofManager const & dofManager,
                                                        arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[numNorm]{};
      real64 subRegionResidualNormalizer[numNorm]{};

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // step 1: compute the norm in the subRegion

      real64 subRegionFlowResidualNorm[1]{};
      real64 subRegionFlowResidualNormalizer[1]{};

      immiscibleMultiphaseKernels::
        ResidualNormKernelFactory::createAndLaunch< parallelDevicePolicy<> >( normType,
                                                                              2,
                                                                              rankOffset,
                                                                              dofKey,
                                                                              localRhs,
                                                                              subRegion,
                                                                              solid,
                                                                              m_nonlinearSolverParameters.m_minNormalizer,
                                                                              subRegionFlowResidualNorm,
                                                                              subRegionFlowResidualNormalizer );
      subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
      subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];

      // step 2: first reduction across meshBodies/regions/subRegions
      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        physicsSolverBaseKernels::LinfResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, localResidualNorm );
      }
      else
      {
        physicsSolverBaseKernels::L2ResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, subRegionResidualNormalizer, localResidualNorm, localResidualNormalizer );
      }
    } );
  } );

  real64 residualNorm = 0.0;
  residualNorm = localResidualNorm[0];
  if( normType == physicsSolverBaseKernels::NormType::Linf )
  {
    physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localResidualNorm[0], residualNorm );
  }
  else
  {
    physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localResidualNorm[0], localResidualNormalizer[0], residualNorm );
  }

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residualNorm );
  }

  return residualNorm;
}

// real64
// ImmiscibleMultiphaseFlow::scalingForSystemSolution( DomainPartition & domain,
//                                                     DofManager const & dofManager,
//                                                     arrayView1d< real64 const > const & localSolution,
//                                                     arrayView1d< real64 const > const & localResidual )
// {
//   GEOS_MARK_FUNCTION; // Trust Region Solver
//   if ( m_trustRegion == 0 )
//   {
//     return 1.0;
//   }

//   // Update solution field
//   updateSolutionField( dofManager, localSolution, domain );

//   // Compute residual norm
//   real64 resNorm = calculateResidualNorm( 0, 0, domain, dofManager, localResidual );

//   // Compute kink scaling factor 
//   real64 localKinkFactor = 1.0;

//   NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
//   FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
//   FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  
//   string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

//   forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
//                                                                MeshLevel const & mesh,
//                                                                string_array const & )
//   {
//     fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
//     {
//       real64 stencilKinkFactor;            

//       // step 1.1: compute the kink damping factor in the stencil

//       typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
//       immiscibleMultiphaseKernels::
//         KinkFactorKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
//                                                                             dofManager.rankOffset(),                                                                            
//                                                                             dofKey,
//                                                                             localSolution,
//                                                                             localResidual,
//                                                                             getName(),
//                                                                             mesh.getElemManager(),
//                                                                             stencilWrapper,
//                                                                             m_hasCapPressure,
//                                                                             resNorm,
//                                                                             stencilKinkFactor );           

//       // step 1.2: local reduction across stencils

//       localKinkFactor = fmin( localKinkFactor, stencilKinkFactor );
//     } );
//   } );

//   // step 1.3: global reduction across mpi ranks

//   real64 globalKinkFactor = MpiWrapper::min( localKinkFactor );

//   // Compute inflection scaling factor
//   real64 localInflectionFactor = 1.0;  

//   forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
//                                                                MeshLevel const & mesh,
//                                                                string_array const & regionNames )
//   {    
//     mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
//                                                                                           CellElementSubRegion const & subRegion )
//     {      
//       // Build wrappers to the fluid, relative permeability and capillary pressure model objects
//       string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
//       TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
//       TwoPhaseImmiscibleFluid::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

//       string const & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
//       BrooksCoreyRelativePermeability const & relPerm = getConstitutiveModel< BrooksCoreyRelativePermeability >( subRegion, relPermName );
//       BrooksCoreyRelativePermeability::KernelWrapper relPermWrapper = relPerm.createKernelWrapper();

//       string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
//       CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

//       BrooksCoreyCapillaryPressure::KernelWrapper* capPressureWrapper = nullptr;
//       if ( m_hasCapPressure )
//       {
//         string const & cappresName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
//         BrooksCoreyCapillaryPressure const & capPressure = getConstitutiveModel< BrooksCoreyCapillaryPressure >( subRegion, cappresName );
//         BrooksCoreyCapillaryPressure::KernelWrapper wrapper = capPressure.createKernelWrapper();
//         capPressureWrapper = &wrapper;
//       }

//       fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
//       {
//         real64 stencilInflectionFactor;

//         // Use flux inflection analysis
//         if ( m_fluxInflection == 1 && stencil.size() > 0 && subRegion.size() > 0 )
//         {
//           // step 2.1a: compute the inflection damping factor in the subRegion

//           typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
//           immiscibleMultiphaseKernels::
//             FluxInflectionFactorKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
//                                                                                           dofManager.rankOffset(),                                                                            
//                                                                                           dofKey,
//                                                                                           localSolution,
//                                                                                           globalKinkFactor,
//                                                                                           getName(),
//                                                                                           mesh.getElemManager(),
//                                                                                           stencilWrapper,
//                                                                                           fluidWrapper,
//                                                                                           relPermWrapper,
//                                                                                           capPressureWrapper,
//                                                                                           m_hasCapPressure,
//                                                                                           stencilInflectionFactor );           

//           // step 2.2a: local reduction across meshBodies/regions/subRegions/stencils

//           localInflectionFactor = fmin( localInflectionFactor, stencilInflectionFactor );
//         }
//         // Use residual inflection analysis
//         else if ( stencil.size() > 0 && subRegion.size() > 0 )
//         {
//           // step 2.1b: compute the inflection damping factor in the subRegion

//           typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
//           immiscibleMultiphaseKernels::
//             ResidualInflectionFactorKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
//                                                                                               dofManager.rankOffset(),                                                                            
//                                                                                               dofKey,
//                                                                                               localSolution,
//                                                                                               localResidual,
//                                                                                               globalKinkFactor,
//                                                                                               getName(),                                                                                              
//                                                                                               mesh.getElemManager(),
//                                                                                               subRegion,
//                                                                                               solid,
//                                                                                               stencilWrapper,
//                                                                                               fluidWrapper,
//                                                                                               relPermWrapper,
//                                                                                               capPressureWrapper,
//                                                                                               m_hasCapPressure,
//                                                                                               resNorm,
//                                                                                               stencilInflectionFactor );           

//           // step 2.2b: local reduction across meshBodies/regions/subRegions/stencils

//           localInflectionFactor = fmin( localInflectionFactor, stencilInflectionFactor );
//         }      
//       } );    
//     } );
//   } );

//   // step 2.3: global reduction across mpi ranks

//   real64 globalInflectionFactor = MpiWrapper::min( localInflectionFactor );

//   globalKinkFactor = fmax( globalKinkFactor, 0.1 ); // 0.8
//   globalInflectionFactor = fmax( globalInflectionFactor, 0.1 ); // 0.4

//   // step 2.4: global combined damping factor
//   real64 scalingFactor = fmax( globalKinkFactor * globalInflectionFactor, 0.1 ); /////// 0.01, 0.1, 0.2, 0.4, 1.0

//   return scalingFactor;
// }

real64
ImmiscibleMultiphaseFlow::scalingForSystemSolution( DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 > const & localSolution,
                                                    arrayView1d< real64 const > const & localResidual,
                                                    real64 const dt,
                                                    real64 const residualNorm,
                                                    integer const newtonIter )
{
  GEOS_MARK_FUNCTION;

  m_currentScaling = m_scalingType;
  
  // Trust Region Solver
  if( m_scalingFactorType == ScalingFactorType::TrustRegion || m_scalingFactorType == ScalingFactorType::TrustRegionFlux )
  {
    if( newtonIter > m_trustRegionParams.minNewtonIter )
    {
      real64 maxGrad = checkMaxGradient( domain, dofManager, localSolution );
      if ( maxGrad > m_trustRegionParams.minGradient )
      {
        return scalingForSystemSolutionTrustRegion( domain, dofManager, localSolution, localResidual, dt, residualNorm, newtonIter );
      }    
    }
    m_currentScaling = ScalingType::Local; // if trust region is not applied, use local scaling
  }

  // Relative and maximum variation scaling
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 scalingFactor = 1.0;
  real64 minPresScalingFactor = 1.0, minSatScalingFactor = 1.0;

  stdVector< MpiWrapper::PairType< real64, globalIndex > > regionDeltaPresMaxLoc;
  stdVector< MpiWrapper::PairType< real64, globalIndex > > regionDeltaSatMaxLoc;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< globalIndex const > const localToGlobalMap = subRegion.localToGlobalMap();
      arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > saturationScalingFactor = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFractionScalingFactor >();

       auto const subRegionData =
       immiscibleMultiphaseKernels::
         SolutionScalingKernelFactory::
         createAndLaunch< parallelDevicePolicy<> >( m_maxRelativePresChange,
                                                    m_maxAbsolutePresChange,                                                    
                                                    m_maxRelativeSatChange,
                                                    m_maxAbsoluteSatChange,
                                                    pressure,
                                                    phaseVolFrac,
                                                    pressureScalingFactor,
                                                    saturationScalingFactor,
                                                    dofManager.rankOffset(),
                                                    m_numPhases,
                                                    dofKey,
                                                    subRegion,
                                                    localSolution );

      if( subRegion.size() > 0 || subRegion.size() !=  subRegion.getNumberOfGhosts() )
      {
        if( m_scalingType == ScalingType::Global )
        {
          scalingFactor = std::min( scalingFactor, subRegionData.localMinVal );
        }

        regionDeltaPresMaxLoc.push_back( { subRegionData.localMaxDeltaPres,
                                           subRegionData.localMaxDeltaPresLoc >= 0 ? localToGlobalMap[subRegionData.localMaxDeltaPresLoc] : -1 } );
        minPresScalingFactor = std::min( minPresScalingFactor, subRegionData.localMinPresScalingFactor );

        regionDeltaSatMaxLoc.push_back( { subRegionData.localMaxDeltaSat,
                                          subRegionData.localMaxDeltaSatLoc >= 0 ? localToGlobalMap[subRegionData.localMaxDeltaSatLoc] : -1 } );
        minSatScalingFactor = std::min( minSatScalingFactor, subRegionData.localMinSatScalingFactor );
      }
    } );
  } );

  auto globalDeltaPresMax = MpiWrapper::max< real64, globalIndex >( regionDeltaPresMaxLoc );
  auto globalDeltaSatMax = MpiWrapper::max< real64, globalIndex >( regionDeltaSatMaxLoc );

  scalingFactor = MpiWrapper::min( scalingFactor );
  minPresScalingFactor = MpiWrapper::min( minPresScalingFactor );
  minSatScalingFactor = MpiWrapper::min( minSatScalingFactor );

  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                         GEOS_FMT( "        {}: Max pressure change = {:.3f} Pa (before scaling) at cell {}",
                                   getName(), globalDeltaPresMax.first, globalDeltaPresMax.second ) );

  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                         GEOS_FMT( "        {}: Max component saturation change = {:.3f} (before scaling) at cell {}",
                                   getName(), globalDeltaSatMax.first, globalDeltaSatMax.second ) );

  if( m_scalingType == ScalingType::Local )
  {
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min pressure scaling factor = {}", getName(), minPresScalingFactor ) );
    GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Min component saturation scaling factor = {}", getName(), minSatScalingFactor ) );
  }

  return LvArray::math::max( scalingFactor, m_minScalingFactor );
}

real64
ImmiscibleMultiphaseFlow::scalingForSystemSolutionTrustRegion( DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               arrayView1d< real64 > const & localSolution,
                                                               arrayView1d< real64 const > const & localResidual,
                                                               real64 const dt,
                                                               real64 const resNorm,
                                                               integer const GEOS_UNUSED_PARAM( newtonIter ) )
{
  GEOS_MARK_FUNCTION;

  if( m_scalingFactorType == ScalingFactorType::TrustRegionFlux && m_scalingType == ScalingType::Local )
  {
    GEOS_ERROR( "Local trust region with flux function not implemented. Terminating..." );   
  }

  // Keep primary variables within physical bounds  
  if ( m_allowSatChopping )
  {
    avoidOutOfBoundPhaseVolFrac( domain, dofManager, localSolution ); 
  }
  if ( m_allowPresChopping )
  {
    avoidOutOfBoundPressure( domain, dofManager, localSolution );
  }

  // Update solution field
  updateSolutionField( dofManager, localSolution, domain );  

  // Reset local restriction factors
  if ( m_scalingType == ScalingType::Local )
  {
    resetLocalScalingFactors( domain );
  }  

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Compute kink scaling factor 
  real64 localKinkFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & subRegion )
    {
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > saturationScalingFactor = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFractionScalingFactor >(); 
    
      fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
      {
        real64 stencilKinkFactor;            

        if ( stencil.size() > 0 && subRegion.size() > 0 )
        {
          // step 1.1: compute the kink damping factor in the stencil

          typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
          immiscibleMultiphaseKernels::
            KinkFactorKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                                                dofManager.rankOffset(),                                                                            
                                                                                dofKey,
                                                                                localSolution,
                                                                                localResidual,
                                                                                getName(),
                                                                                mesh.getElemManager(),
                                                                                stencilWrapper,
                                                                                m_hasCapPressure,
                                                                                resNorm,                                                                       
                                                                                m_trustRegionParams,
                                                                                m_scalingType,
                                                                                pressureScalingFactor,
                                                                                saturationScalingFactor,
                                                                                stencilKinkFactor );           

          // step 1.2: local reduction across meshBodies/regions/subRegions/stencils

          localKinkFactor = fmin( localKinkFactor, stencilKinkFactor );
        }
      } );
    } );
  } );

  // step 1.3: global reduction across mpi ranks

  real64 globalKinkFactor = MpiWrapper::min( localKinkFactor );

  // Compute inflection scaling factor
  real64 localInflectionFactor = 1.0;  

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {    
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & subRegion )
    {      
      // Build wrappers to the fluid, relative permeability and capillary pressure model objects
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
      TwoPhaseImmiscibleFluid::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

      string const & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      BrooksCoreyRelativePermeability const & relPerm = getConstitutiveModel< BrooksCoreyRelativePermeability >( subRegion, relPermName );
      BrooksCoreyRelativePermeability::KernelWrapper relPermWrapper = relPerm.createKernelWrapper();

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      // string const & poroName = subRegion.getReference< string >( viewKeyStruct::porosityModelNameString() );
      // PressurePorosity const & porosity = getConstitutiveModel< PressurePorosity >( subRegion, poroName );
      // PressurePorosity::KernelWrapper poroWrapper = porosity.createKernelUpdates();
      // PressurePorosity::KernelWrapper* poroWrapper = nullptr;

      BrooksCoreyCapillaryPressure::KernelWrapper* capPressureWrapper = nullptr;
      if ( m_hasCapPressure )
      {
        string const & cappresName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        BrooksCoreyCapillaryPressure const & capPressure = getConstitutiveModel< BrooksCoreyCapillaryPressure >( subRegion, cappresName );
        BrooksCoreyCapillaryPressure::KernelWrapper wrapper = capPressure.createKernelWrapper();
        capPressureWrapper = &wrapper;
      }

      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > saturationScalingFactor = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFractionScalingFactor >();      

      fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
      {
        real64 stencilInflectionFactor;

        // Use flux inflection analysis
        if ( m_scalingFactorType == ScalingFactorType::TrustRegionFlux && stencil.size() > 0 && subRegion.size() > 0 )
        {
          // step 2.1a: compute the inflection damping factor in the subRegion

          typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
          immiscibleMultiphaseKernels::
            FluxInflectionFactorKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                                                          dofManager.rankOffset(),                                                                            
                                                                                          dofKey,
                                                                                          localSolution,
                                                                                          localResidual,
                                                                                          globalKinkFactor,
                                                                                          getName(),
                                                                                          mesh.getElemManager(),
                                                                                          stencilWrapper,
                                                                                          fluidWrapper,
                                                                                          relPermWrapper,
                                                                                          capPressureWrapper,
                                                                                          m_hasCapPressure,
                                                                                          resNorm,
                                                                                          m_trustRegionParams,
                                                                                          m_scalingType,
                                                                                          stencilInflectionFactor );           

          // step 2.2a: local reduction across meshBodies/regions/subRegions/stencils
          localInflectionFactor = fmin( localInflectionFactor, stencilInflectionFactor );
        }
        // Use residual inflection analysis
        else if ( stencil.size() > 0 && subRegion.size() > 0 )
        {
          // step 2.1b: compute the inflection damping factor in the subRegion

          typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();
          immiscibleMultiphaseKernels::
            ResidualInflectionFactorKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                                                              dofManager.rankOffset(),                                                                            
                                                                                              dofKey,
                                                                                              dt,
                                                                                              localSolution,
                                                                                              localResidual,
                                                                                              globalKinkFactor,
                                                                                              getName(),                                                                                              
                                                                                              mesh.getElemManager(),
                                                                                              subRegion,
                                                                                              fluid,
                                                                                              solid,
                                                                                              stencilWrapper,
                                                                                              fluidWrapper,
                                                                                              relPermWrapper,
                                                                                              capPressureWrapper,                                                                                              
                                                                                              m_hasCapPressure,
                                                                                              resNorm,
                                                                                              m_trustRegionParams,
                                                                                              m_scalingType,
                                                                                              pressureScalingFactor,
                                                                                              saturationScalingFactor,
                                                                                              stencilInflectionFactor );           

          // step 2.2b: local reduction across meshBodies/regions/subRegions/stencils

          localInflectionFactor = fmin( localInflectionFactor, stencilInflectionFactor );
        }      
      } );    
    } );
  } );

  // step 2.3: global reduction across mpi ranks

  real64 globalInflectionFactor = MpiWrapper::min( localInflectionFactor );
  
  // step 2.4: global combined damping factor
  real64 scalingFactor = fmax( globalKinkFactor * globalInflectionFactor, m_minScalingFactor );

  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Global discontinuity scaling factor = {}", getName(), fmt::format( "{:.{}f}", globalKinkFactor, 3 ) ) ); 
  GEOS_LOG_LEVEL_RANK_0( logInfo::Solution, GEOS_FMT( "        {}: Global inflection scaling factor = {}", getName(), fmt::format( "{:.{}f}", globalInflectionFactor, 3 ) ) );

  return scalingFactor;
}

bool ImmiscibleMultiphaseFlow::checkSystemSolution( DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localSolution,
                                                    real64 const scalingFactor )
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer localCheck = 1;
  real64 minPres = 0.0, minSat = 0.0, maxSat = 1.0;
  integer numNegPres = 0, numNegSat = 0, numUnitySat = 0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pressure =
        subRegion.getField< fields::flow::pressure >();
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > saturationScalingFactor = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFractionScalingFactor >();

      auto subRegionData =
      immiscibleMultiphaseKernels::
        SolutionCheckKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                   subRegion,
                                                   dofManager.rankOffset(),
                                                   dofKey,
                                                   localSolution,                                                   
                                                   pressure,
                                                   phaseVolFrac,
                                                   m_currentScaling,
                                                   scalingFactor,
                                                   pressureScalingFactor,
                                                   saturationScalingFactor,                                                   
                                                   m_allowOutOfBoundPressure,
                                                   m_allowOutOfBoundSaturation );

    localCheck = fmin( localCheck, subRegionData.localMinVal );

    minPres  = fmin( minPres, subRegionData.localMinPres );
    minSat = fmin( minSat, subRegionData.localMinSat );
    maxSat = fmax( maxSat, subRegionData.localMaxSat );
    numNegPres += subRegionData.localNumNegPressures;
    numNegSat += subRegionData.localNumNegSat;
    numUnitySat += subRegionData.localNumUnitySat;
    } );
  } );

  minPres  = MpiWrapper::min( minPres );
  minSat = MpiWrapper::min( minSat );
  maxSat = MpiWrapper::max( maxSat );
  numNegPres = MpiWrapper::sum( numNegPres );
  numNegSat = MpiWrapper::sum( numNegSat );
  numUnitySat = MpiWrapper::sum( numUnitySat );  

  if( numNegPres > 0 )
    {GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                            GEOS_FMT( "        {}: Number of negative pressure values: {}, minimum value: {} Pa",
                                      getName(), numNegPres, fmt::format( "{:.{}f}", minPres, 3 ) ) );}

  if( numNegSat > 0 )
    {GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                            GEOS_FMT( "        {}: Number of negative component saturation values: {}, minimum value: {}",
                                      getName(), numNegSat, fmt::format( "{:.{}f}", minSat, 3 ) ) );}

  if( numUnitySat > 0 )
    {GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                            GEOS_FMT( "        {}: Number of saturation values above unity: {}, maximum value: {}",
                                      getName(), numUnitySat, fmt::format( "{:.{}f}", maxSat, 3 ) ) );}  
  
  return MpiWrapper::min( localCheck );
}                                                      

void ImmiscibleMultiphaseFlow::applySystemSolution( DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localSolution,
                                                    real64 const scalingFactor,
                                                    real64 const dt,
                                                    DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( dt );

  bool const globalScaling = m_currentScaling == ScalingType::Global;

  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );

  // 1. apply the pressure update
  if ( globalScaling )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 fields::flow::pressureScalingFactor::key(),
                                 pressureMask );
  }  

  // 2. apply the phaseVolumeFraction update
  if( globalScaling )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key(),
                                 scalingFactor,
                                 ~pressureMask );
  }
  else
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key(),
                                 fields::immiscibleMultiphaseFlow::phaseVolumeFractionScalingFactor::key(),
                                 ~pressureMask );
  }

  // 3. ensure primary variables are within physical bounds
  if ( m_allowSatChopping )
  {
    chopOutOfBoundPhaseVolFrac( domain );
  }
  if ( m_allowPresChopping )
  {
    chopOutOfBoundPressure( domain );
  }

  // 4. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    std::vector< string > fields{ fields::flow::pressure::key(), fields::immiscibleMultiphaseFlow::phaseVolumeFraction::key() };

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}


void ImmiscibleMultiphaseFlow::updateSolutionField( DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localSolution,
                                                    DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  DofManager::CompMask allDofMask( m_numDofPerCell, 0, 2 );

  // copy updates to field
  dofManager.copyVectorToField( localSolution,
                                viewKeyStruct::elemDofFieldString(),
                                fields::immiscibleMultiphaseFlow::solutionUpdate::key(),
                                1.0,
                                allDofMask );

  // synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    std::vector< string > fields{ fields::immiscibleMultiphaseFlow::solutionUpdate::key() };

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}


void ImmiscibleMultiphaseFlow::updateVolumeConstraint( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolumeFraction = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    phaseVolumeFraction[ei][1] = 1.0 - phaseVolumeFraction[ei][0];
  } );
}


void ImmiscibleMultiphaseFlow::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      arrayView1d< real64 > const & pres =
        subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const & pres_n =
        subRegion.template getField< fields::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      // after the update, save the new saturation
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      phaseVolFrac.setValues< parallelDevicePolicy<> >( phaseVolFrac_n );

      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const & phaseMass_n =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const & phaseMass =
        subRegion.template getField< fields::immiscibleMultiphaseFlow::phaseMass >();
      phaseMass.setValues< parallelDevicePolicy<> >( phaseMass_n );

      if( m_isThermal )
      {
        arrayView1d< real64 > const & temp =
          subRegion.template getField< fields::flow::temperature >();
        arrayView1d< real64 const > const & temp_n =
          subRegion.template getField< fields::flow::temperature_n >();
        temp.setValues< parallelDevicePolicy<> >( temp_n );
      }

      // update porosity, permeability
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateFluidState( subRegion );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::implicitStepComplete( real64 const & time,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  // Step 1: save the converged aquifer state
  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // Step 3: save the converged solid state
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      if( m_keepVariablesConstantDuringInitStep )
      {
        porousMaterial.ignoreConvergedState(); // newPorosity <- porosity_n
      }
      else
      {
        porousMaterial.saveConvergedState(); // porosity_n <- porosity
      }

      // Step 4: save converged state for the relperm model to handle hysteresis
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      string const & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPermMaterial =
        getConstitutiveModel< RelativePermeabilityBase >( subRegion, relPermName );
      relPermMaterial.saveConvergedPhaseVolFractionState( phaseVolFrac );

      // Step 5: if capillary pressure is supported, send the converged porosity and permeability to the capillary pressure model
      // note: this is needed when the capillary pressure depends on porosity and permeability (Leverett J-function for instance)
      if( m_hasCapPressure )
      {
        arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

        string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
        PermeabilityBase const & permeabilityMaterial =
          getConstitutiveModel< PermeabilityBase >( subRegion, permName );
        arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

        string const & capPressName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        CapillaryPressureBase const & capPressureMaterial =
          getConstitutiveModel< CapillaryPressureBase >( subRegion, capPressName );
        capPressureMaterial.saveConvergedRockState( porosity, permeability );
      }
    } );
  } );
}

void ImmiscibleMultiphaseFlow::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::saveConvergedState( subRegion );

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
  phaseVolFrac_n.setValues< parallelDevicePolicy<> >( phaseVolFrac );

  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseMass =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass >();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseMass_n =
    subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass_n >();
  phaseMass_n.setValues< parallelDevicePolicy<> >( phaseMass );

}

void ImmiscibleMultiphaseFlow::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             auto & subRegion )
    {
      // update porosity, permeability, and solid internal energy
      updatePorosityAndPermeability( subRegion );
      // update all fluid properties
      updateVolumeConstraint( subRegion );
      updateFluidState( subRegion );
    } );
  } );
}

real64 ImmiscibleMultiphaseFlow::setNextDtBasedOnStateChange( real64 const & currentDt,
                                                              DomainPartition & domain )
{
  if( m_targetRelativePresChange >= 1.0 &&
      m_targetPhaseVolFracChange >= 1.0 )
  {
    return LvArray::NumericLimits< real64 >::max;
  }

  real64 maxRelativePresChange = 0.0;
  real64 maxAbsolutePhaseVolFracChange = 0.0;

  integer const numPhase = m_numPhases;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n = subRegion.getField< fields::flow::pressure_n >();
      arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac =
        subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac_n =
        subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction_n >();

      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPresChange( 0.0 );
      RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPhaseVolFracChange( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          // switch from relative to absolute when values less than 1
          subRegionMaxPresChange.max( LvArray::math::abs( pres[ei] - pres_n[ei] ) / LvArray::math::max( LvArray::math::abs( pres_n[ei] ), 1.0 ) );
          for( integer ip = 0; ip < numPhase; ++ip )
          {
            subRegionMaxPhaseVolFracChange.max( LvArray::math::abs( phaseVolFrac[ei][ip] - phaseVolFrac_n[ei][ip] ) );
          }
        }
      } );

      maxRelativePresChange = LvArray::math::max( maxRelativePresChange, subRegionMaxPresChange.get() );
      maxAbsolutePhaseVolFracChange = LvArray::math::max( maxAbsolutePhaseVolFracChange, subRegionMaxPhaseVolFracChange.get() );

    } );
  } );

  maxRelativePresChange = MpiWrapper::max( maxRelativePresChange );
  maxAbsolutePhaseVolFracChange = MpiWrapper::max( maxAbsolutePhaseVolFracChange );

  GEOS_LOG_LEVEL_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max relative pressure change during time step = {} %",
                                                      getName(), GEOS_FMT( "{:.{}f}", 100*maxRelativePresChange, 3 ) ) );
  GEOS_LOG_LEVEL_RANK_0( logInfo::TimeStep, GEOS_FMT( "{}: max absolute phase volume fraction change during time step = {}",
                                                      getName(), GEOS_FMT( "{:.{}f}", maxAbsolutePhaseVolFracChange, 3 ) ) );

  real64 const eps = LvArray::NumericLimits< real64 >::epsilon;

  real64 const nextDtPressure = currentDt * ( 1.0 + m_solutionChangeScalingFactor ) * m_targetRelativePresChange
                                / std::max( eps, maxRelativePresChange + m_solutionChangeScalingFactor * m_targetRelativePresChange );
  if( m_nonlinearSolverParameters.getLogLevel() > 0 )
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on pressure change = {}", getName(), nextDtPressure ));
  real64 const nextDtPhaseVolFrac = currentDt * ( 1.0 + m_solutionChangeScalingFactor ) * m_targetPhaseVolFracChange
                                    / std::max( eps, maxAbsolutePhaseVolFracChange + m_solutionChangeScalingFactor * m_targetPhaseVolFracChange );
  if( m_nonlinearSolverParameters.getLogLevel() > 0 )
    GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on phase volume fraction change = {}", getName(), nextDtPhaseVolFrac ));

  return std::min( nextDtPressure, nextDtPhaseVolFrac );

}

void ImmiscibleMultiphaseFlow::chopOutOfBoundPhaseVolFrac( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  integer const numPhase = m_numPhases;
  real64 const minSat = 0.0; // technically should be the irreducible saturation
  real64 const maxSat = 1.0; // technically should be the residual saturation

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac =
      subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          for( integer ip = 0; ip < numPhase; ++ip )
          {
            if( phaseVolFrac[ei][ip] < minSat )
            {
              phaseVolFrac[ei][ip] = minSat;
            }
            else if( phaseVolFrac[ei][ip] > maxSat )
            {
              phaseVolFrac[ei][ip] = maxSat;
            }
          }
        }
      } );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::chopOutOfBoundPressure( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  real64 const minPres = 0.0; // could be the minimum presure for which the fluid model is valid
  real64 const maxPres = LvArray::NumericLimits< real64 >::max; // could be the maximum presure for which the fluid model is valid

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView1d< real64 > const pressure = subRegion.getField< fields::flow::pressure >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          if( pressure[ei] < minPres )
          {
            pressure[ei] = minPres;
          }
          else if( pressure[ei] > maxPres )
          {
            pressure[ei] = maxPres;
          }
        }
      } );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::avoidOutOfBoundPressure( DomainPartition & domain,
                                                        DofManager const & dofManager,
                                                        arrayView1d< real64 > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  real64 const minPres = 0.0; // could be the minimum presure for which the fluid model is valid  

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );  

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();      
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 > const pressure = subRegion.getField< fields::flow::pressure >();      

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          globalIndex const globalRow = dofNumber[ei] - rankOffset;
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );
          
          if( pressure[ei] + localSolution[localRow] < minPres )
          {
            localSolution[localRow] = minPres - pressure[ei];
          }          
        }
      } );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::avoidOutOfBoundPhaseVolFrac( DomainPartition & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 > const & localSolution )
{
  GEOS_MARK_FUNCTION;
  
  real64 const minSat = 0.0; // technically should be the irreducible saturation
  real64 const maxSat = 1.0; // technically should be the residual saturation  

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVolFrac = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();      

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          globalIndex const globalRow = dofNumber[ei] - rankOffset;
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );

          if( phaseVolFrac[ei][0] + localSolution[localRow + 1] < minSat )
          {
            localSolution[localRow + 1] = -phaseVolFrac[ei][0];
          }
          else if( phaseVolFrac[ei][0] + localSolution[localRow + 1] > maxSat )     
          {
            localSolution[localRow + 1] = maxSat - phaseVolFrac[ei][0];
          }  
        }
      } );
    } );
  } );
}

void ImmiscibleMultiphaseFlow::resetLocalScalingFactors( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

      arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::flow::pressureScalingFactor >();
      arrayView1d< real64 > saturationScalingFactor = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFractionScalingFactor >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          pressureScalingFactor[ei] = 1.0;
          saturationScalingFactor[ei] = 1.0;
        }
      } );
    } );
  } );
}

real64 ImmiscibleMultiphaseFlow::checkMaxGradient( DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   arrayView1d< real64 const > const & localSolution )
{
  GEOS_MARK_FUNCTION;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  real64 maxGrad = 0.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      RAJA::ReduceMax< parallelDeviceReduce, real64 > maxSubRegionGrad( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] < 0 )
        {
          globalIndex const globalRow = dofNumber[ei] - rankOffset;
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOS_ASSERT_GE( localRow, 0 );

          maxSubRegionGrad.max( LvArray::math::abs( localSolution[localRow + 1] ) );     
        }
      } );

      maxGrad = LvArray::math::max( maxGrad, maxSubRegionGrad.get() );
    } );
  } );

  maxGrad = MpiWrapper::max( maxGrad );

  return maxGrad;
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ImmiscibleMultiphaseFlow, string const &, Group * const )

} // namespace geos
