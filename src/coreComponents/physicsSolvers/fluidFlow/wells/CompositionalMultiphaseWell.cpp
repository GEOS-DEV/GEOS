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
 * @file CompositionalMultiphaseWell.cpp
 */

#include "CompositionalMultiphaseWell.hpp"


#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/MultiFluidExtrinsicData.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityExtrinsicData.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationData.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace compositionalMultiphaseWellKernels;

CompositionalMultiphaseWell::CompositionalMultiphaseWell( const string & name,
                                                          Group * const parent )
  :
  WellSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_useMass( false ),
  m_maxCompFracChange( 1.0 ),
  m_maxRelativePresChange( 0.2 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 ),
  m_targetPhaseIndex( -1 )
{
  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (absolute) change in a component fraction between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString(), &m_maxRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (relative) change in pressure between two Newton iterations (recommended with rate control)" );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );
}

void CompositionalMultiphaseWell::postProcessInput()
{
  WellSolverBase::postProcessInput();

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         "CompositionalMultiphaseWell named " << getName() <<
                         ": The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         "CompositionalMultiphaseWell named " << getName() <<
                         ": The maximum absolute change in component fraction must larger or equal to 0.0" );

}

void CompositionalMultiphaseWell::registerDataOnMesh( Group & meshBodies )
{
  WellSolverBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  forMeshTargets( meshBodies, [&]( string const &,
                                   MeshLevel & mesh,
                                   arrayView1d< string const > const & regionNames )
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
  m_numDofPerWellElement = m_numComponents + 2; // 1 pressure + NC compositions + 1 connectionRate

  // loop over the wells
  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      string const & fluidName = getConstitutiveName< MultiFluidBase >( subRegion );

      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::pressure >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::pressure_n >( getName() );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::temperature >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::temperature_n >( getName() );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.

      subRegion.registerExtrinsicData< extrinsicMeshData::well::globalCompDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::globalCompDensity_n >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate_n >( getName() );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::globalCompFraction >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dPressure >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::totalMassDensity >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dPressure >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      subRegion.registerExtrinsicData< extrinsicMeshData::well::totalDensity_n >( getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::phaseDensity_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< extrinsicMeshData::well::phaseComponentFraction_n >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

      PerforationData & perforationData = *subRegion.getPerforationData();
      perforationData.registerExtrinsicData< extrinsicMeshData::well::compPerforationRate >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );
      perforationData.registerExtrinsicData< extrinsicMeshData::well::dCompPerforationRate_dPres >( getName() ).
        reference().resizeDimension< 1, 2 >( 2, m_numComponents );
      perforationData.registerExtrinsicData< extrinsicMeshData::well::dCompPerforationRate_dComp >( getName() ).
        reference().resizeDimension< 1, 2, 3 >( 2, m_numComponents, m_numComponents );

      WellControls & wellControls = getWellControls( subRegion );
      wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentBHP_dPresString() ).
        setRestartFlags( RestartFlags::NO_WRITE );
      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentBHP_dCompDensString() ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numComponents );

      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numPhases );
      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dPresString() ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numPhases );
      wellControls.registerWrapper< array2d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0, 1 >( m_numPhases, m_numComponents );
      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dRateString() ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numPhases );

      wellControls.registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentTotalVolRate_dPresString() ).
        setRestartFlags( RestartFlags::NO_WRITE );
      wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentTotalVolRate_dCompDensString() ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 ).
        reference().resizeDimension< 0 >( m_numComponents );
      wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentTotalVolRate_dRateString() ).
        setRestartFlags( RestartFlags::NO_WRITE );
    } );
  } );

}

void CompositionalMultiphaseWell::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{

  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
  GEOSX_THROW_IF( fluidName.empty(),
                  GEOSX_FMT( "Fluid model not found on subregion {}", subRegion.getName() ),
                  InputError );

  string & relPermName = subRegion.registerWrapper< string >( viewKeyStruct::relPermNamesString() ).
                           setPlotLevel( PlotLevel::NOPLOT ).
                           setRestartFlags( RestartFlags::NO_WRITE ).
                           setSizedFromParent( 0 ).
                           setDescription( "Name of the relative permeability constitutive model to use" ).
                           reference();

  relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );

  GEOSX_THROW_IF( relPermName.empty(),
                  GEOSX_FMT( "Relative permeability model not found on subregion {}", subRegion.getName() ),
                  InputError );

}
namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_THROW_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                         GEOSX_FMT( "Mismatch in number of phases between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                         InputError );

  for( integer ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOSX_THROW_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                           GEOSX_FMT( "Mismatch in phase names between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                           InputError );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_THROW_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                         GEOSX_FMT( "Mismatch in number of components between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                         InputError );

  for( integer ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOSX_THROW_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                           GEOSX_FMT( "Mismatch in component names between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                           InputError );
  }
}

}

void CompositionalMultiphaseWell::validateConstitutiveModels( DomainPartition const & domain ) const
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
  string const referenceFluidName = flowSolver.referenceFluidModelName();
  MultiFluidBase const & referenceFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_referenceFluidModelName );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion const & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      compareMultiphaseModels( fluid, referenceFluid );
      compareMulticomponentModels( fluid, referenceFluid );

      string const & relpermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      RelativePermeabilityBase const & relPerm = getConstitutiveModel< RelativePermeabilityBase >( subRegion, relpermName );
      compareMultiphaseModels( relPerm, referenceFluid );
    } );

  } );
}

void CompositionalMultiphaseWell::validateInjectionStreams( WellElementSubRegion const & subRegion ) const
{
  WellControls const & wellControls = getWellControls( subRegion );

  // check well injection stream for injectors
  if( wellControls.isInjector() )
  {
    arrayView1d< real64 const > const & injection = wellControls.getInjectionStream();
    real64 compFracSum = 0;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      real64 const compFrac = injection[ic];
      GEOSX_THROW_IF( ( compFrac < 0.0 ) || ( compFrac > 1.0 ),
                      "WellControls '" << wellControls.getName() << "'" <<
                      ": Invalid injection stream for well " << subRegion.getName(),
                      InputError );
      compFracSum += compFrac;
    }
    GEOSX_THROW_IF( ( compFracSum < 1.0 - std::numeric_limits< real64 >::epsilon() ) ||
                    ( compFracSum > 1.0 + std::numeric_limits< real64 >::epsilon() ),
                    "WellControls '" << wellControls.getName() << "'" <<
                    ": Invalid injection stream for well " << subRegion.getName(),
                    InputError );
  }
}

void CompositionalMultiphaseWell::validateWellConstraints( WellElementSubRegion const & subRegion )
{
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());
  MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  // now that we know we are single-phase, we can check a few things in the constraints
  WellControls const & wellControls = getWellControls( subRegion );
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const & targetTotalRate = wellControls.getTargetTotalRate( m_currentTime + m_currentDt );
  real64 const & targetPhaseRate = wellControls.getTargetPhaseRate( m_currentTime + m_currentDt );
  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 const & surfaceTemp = wellControls.getSurfaceTemperature();

  GEOSX_THROW_IF( wellControls.isInjector() && currentControl == WellControls::Control::PHASEVOLRATE,
                  "WellControls named " << wellControls.getName() <<
                  ": Phase rate control is not available for injectors",
                  InputError );
  GEOSX_THROW_IF( wellControls.isProducer() && currentControl == WellControls::Control::TOTALVOLRATE,
                  "WellControls named " << wellControls.getName() <<
                  ": Phase rate control is not available for producers",
                  InputError );

  GEOSX_THROW_IF( wellControls.isInjector() && targetTotalRate < 0.0,
                  "WellControls named " << wellControls.getName() <<
                  ": Target total rate cannot be negative for injectors",
                  InputError );
  GEOSX_THROW_IF( wellControls.isInjector() && !isZero( targetPhaseRate ),
                  "WellControls named " << wellControls.getName() <<
                  ": Target phase rate cannot be used for injectors",
                  InputError );

  // The user always provides positive rates, but these rates are later multiplied by -1 internally for producers
  GEOSX_THROW_IF( wellControls.isProducer() && targetPhaseRate > 0.0,
                  "WellControls named " << wellControls.getName() <<
                  ": Target phase rate cannot be negative for producers",
                  InputError );
  GEOSX_THROW_IF( wellControls.isProducer() && !isZero( targetTotalRate ),
                  "WellControls named " << wellControls.getName() <<
                  ": Target total rate cannot be used for producers",
                  InputError );

  GEOSX_THROW_IF( useSurfaceConditions && surfaceTemp <= 0,
                  "WellControls named " << wellControls.getName() <<
                  ": Surface temperature must be set to a strictly positive value",
                  InputError );

  // Find target phase index for phase rate constraint
  for( integer ip = 0; ip < fluid.numFluidPhases(); ++ip )
  {
    if( fluid.phaseNames()[ip] == wellControls.getTargetPhaseName() )
    {
      m_targetPhaseIndex = ip;
    }
  }
  GEOSX_THROW_IF( wellControls.isProducer() && m_targetPhaseIndex == -1,
                  "WellControls '" << wellControls.getName() <<
                  "': Phase " << wellControls.getTargetPhaseName() << " not found for well control " << wellControls.getName(),
                  InputError );
}

void CompositionalMultiphaseWell::initializePostSubGroups()
{
  WellSolverBase::initializePostSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  validateConstitutiveModels( domain );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< MultiFluidBase >( subRegion );
      GEOSX_THROW_IF( fluidName.empty(),
                      GEOSX_FMT( "Fluid model not found on subregion {}", subRegion.getName() ),
                      InputError );

      string & relPermName = subRegion.getReference< string >( viewKeyStruct::relPermNamesString() );
      relPermName = getConstitutiveName< RelativePermeabilityBase >( subRegion );
      GEOSX_THROW_IF( relPermName.empty(),
                      GEOSX_FMT( "Fluid model not found on subregion {}", subRegion.getName() ),
                      InputError );

      validateInjectionStreams( subRegion );
      validateWellConstraints( subRegion );

    } );
  } );
}

void CompositionalMultiphaseWell::initializePostInitialConditionsPreSubGroups()
{
  WellSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    // loop over the wells
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      fluid.setMassFlag( m_useMass );
    } );
  } );
}

void CompositionalMultiphaseWell::updateComponentFraction( WellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  isothermalCompositionalMultiphaseBaseKernels::
    ComponentFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               subRegion );

}

void CompositionalMultiphaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  integer const numComp = m_numComponents;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & pres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();

  arrayView1d< real64 > const & totalMassDens =
    subRegion.getExtrinsicData< extrinsicMeshData::well::totalMassDensity >();
  arrayView1d< real64 > const & dTotalMassDens_dPres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dPressure >();
  arrayView2d< real64, compflow::USD_FLUID_DC > const & dTotalMassDens_dCompDens =
    subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dGlobalCompDensity >();

  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

  // control data

  WellControls & wellControls = getWellControls( subRegion );

  real64 const & refGravCoef = wellControls.getReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  real64 & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dPresString() );
  arrayView1d< real64 > const & dCurrentBHP_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dCompDensString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [&numComp,
                              pres,
                              totalMassDens,
                              dTotalMassDens_dPres,
                              dTotalMassDens_dCompDens,
                              wellElemGravCoef,
                              &currentBHP,
                              &dCurrentBHP_dPres,
                              dCurrentBHP_dCompDens,
                              &iwelemRef,
                              &refGravCoef] ( localIndex const )
  {
    real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
    currentBHP = pres[iwelemRef] + totalMassDens[iwelemRef] * diffGravCoef;
    dCurrentBHP_dPres = 1 + dTotalMassDens_dPres[iwelemRef] * diffGravCoef;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      dCurrentBHP_dCompDens[ic] = dTotalMassDens_dCompDens[iwelemRef][ic] * diffGravCoef;
    }
  } );
}

void CompositionalMultiphaseWell::updateVolRatesForConstraint( WellElementSubRegion & subRegion )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  integer constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & pres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const & temp =
    subRegion.getExtrinsicData< extrinsicMeshData::well::temperature >();
  arrayView1d< real64 const > const & connRate =
    subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >();

  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac =
    subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();
  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens =
    subRegion.getExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >();

  // fluid data

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac = fluid.phaseFraction();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac = fluid.dPhaseFraction();

  arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens = fluid.totalDensity();
  arrayView3d< real64 const, multifluid::USD_FLUID_DC > const & dTotalDens = fluid.dTotalDensity();

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens = fluid.phaseDensity();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens = fluid.dPhaseDensity();

  // control data

  WellControls & wellControls = getWellControls( subRegion );

  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 const & surfacePres = wellControls.getSurfacePressure();
  real64 const & surfaceTemp = wellControls.getSurfaceTemperature();

  arrayView1d< real64 > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dPres =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dPresString() );
  arrayView2d< real64 > const & dCurrentPhaseVolRate_dCompDens =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dRateString() );

  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
  real64 & dCurrentTotalVolRate_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dPresString() );
  arrayView1d< real64 > const & dCurrentTotalVolRate_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dCompDensString() );
  real64 & dCurrentTotalVolRate_dRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dRateString() );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [&numComp,
                                &numPhase,
                                fluidWrapper,
                                pres,
                                temp,
                                compFrac,
                                dCompFrac_dCompDens,
                                connRate,
                                totalDens,
                                dTotalDens,
                                phaseDens,
                                dPhaseDens,
                                phaseFrac,
                                dPhaseFrac,
                                &useSurfaceConditions,
                                &surfacePres,
                                &surfaceTemp,
                                &currentTotalVolRate,
                                &dCurrentTotalVolRate_dPres,
                                dCurrentTotalVolRate_dCompDens,
                                &dCurrentTotalVolRate_dRate,
                                currentPhaseVolRate,
                                dCurrentPhaseVolRate_dPres,
                                dCurrentPhaseVolRate_dCompDens,
                                dCurrentPhaseVolRate_dRate,
                                &iwelemRef] ( localIndex const )
    {
      using Deriv = multifluid::DerivativeOffset;

      stackArray1d< real64, maxNumComp > work( numComp );

      // Step 1: evaluate the phase and total density in the reference element

      //    We need to evaluate the density as follows:
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element
      if( useSurfaceConditions )
      {
        // we need to compute the surface density
        fluidWrapper.update( iwelemRef, 0, surfacePres, surfaceTemp, compFrac[iwelemRef] );
      }
      else
      {
        real64 const refPres = pres[iwelemRef];
        fluidWrapper.update( iwelemRef, 0, refPres, temp[iwelemRef], compFrac[iwelemRef] );
      }

      // Step 2: update the total volume rate

      real64 const currentTotalRate = connRate[iwelemRef];

      // Step 2.1: compute the inverse of the total density and derivatives

      real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];
      real64 const dTotalDensInv_dPres = -dTotalDens[iwelemRef][0][Deriv::dP] * totalDensInv * totalDensInv;
      stackArray1d< real64, maxNumComp > dTotalDensInv_dCompDens( numComp );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalDensInv_dCompDens[ic] = -dTotalDens[iwelemRef][0][Deriv::dC+ic] * totalDensInv * totalDensInv;
      }
      applyChainRuleInPlace( numComp, dCompFrac_dCompDens[iwelemRef], dTotalDensInv_dCompDens, work.data() );

      // Step 2.2: divide the total mass/molar rate by the total density to get the total volumetric rate
      currentTotalVolRate = currentTotalRate * totalDensInv;
      dCurrentTotalVolRate_dPres = ( useSurfaceConditions ==  0 ) * currentTotalRate * dTotalDensInv_dPres;
      dCurrentTotalVolRate_dRate = totalDensInv;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dCurrentTotalVolRate_dCompDens[ic] = currentTotalRate * dTotalDensInv_dCompDens[ic];
      }

      // Step 3: update the phase volume rate
      for( integer ip = 0; ip < numPhase; ++ip )
      {

        // Step 3.1: compute the inverse of the (phase density * phase fraction) and derivatives

        // skip the rest of this function if phase ip is absent
        bool const phaseExists = (phaseFrac[iwelemRef][0][ip] > 0);
        if( !phaseExists )
        {
          continue;
        }

        real64 const phaseDensInv =  1.0 / phaseDens[iwelemRef][0][ip];
        real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;
        real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac[iwelemRef][0][ip][Deriv::dP] * phaseDensInv
                                                         - dPhaseDens[iwelemRef][0][ip][Deriv::dP] * phaseFracTimesPhaseDensInv * phaseDensInv;


        // Step 3.2: divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
        currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
        dCurrentPhaseVolRate_dPres[ip] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
        dCurrentPhaseVolRate_dRate[ip] = phaseFracTimesPhaseDensInv;
        for( integer ic = 0; ic < numComp; ++ic )
        {
          dCurrentPhaseVolRate_dCompDens[ip][ic] = -phaseFracTimesPhaseDensInv * dPhaseDens[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] += dPhaseFrac[iwelemRef][0][ip][Deriv::dC+ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] *= currentTotalRate;
        }
        applyChainRuleInPlace( numComp, dCompFrac_dCompDens[iwelemRef], dCurrentPhaseVolRate_dCompDens[ip], work.data() );
      }
    } );
  } );
}


void CompositionalMultiphaseWell::updateFluidModel( WellElementSubRegion & subRegion )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres = subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const & temp = subRegion.getExtrinsicData< extrinsicMeshData::well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac =
    subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

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

void CompositionalMultiphaseWell::updatePhaseVolumeFraction( WellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
  isothermalCompositionalMultiphaseBaseKernels::
    PhaseVolumeFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid );

}

void CompositionalMultiphaseWell::updateTotalMassDensity( WellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  MultiFluidBase & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );

  TotalMassDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid );
}

void CompositionalMultiphaseWell::updateSubRegionState( WellElementSubRegion & subRegion )
{
  // update properties
  updateComponentFraction( subRegion );

  // update volumetric rates for the well constraints
  // note: this must be called before updateFluidModel
  updateVolRatesForConstraint( subRegion );

  // update densities, phase fractions, phase volume fractions
  updateFluidModel( subRegion );
  updatePhaseVolumeFraction( subRegion );
  updateTotalMassDensity( subRegion );

  // update the current BHP pressure
  updateBHPForConstraint( subRegion );

  // note: the perforation rates are updated separately
}

void CompositionalMultiphaseWell::initializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // TODO: change the way we access the flowSolver here
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  // loop over the wells
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();
    PresTempCompFracInitializationKernel::CompFlowAccessors resCompFlowAccessors( mesh.getElemManager(), flowSolver.getName() );
    PresTempCompFracInitializationKernel::MultiFluidAccessors resMultiFluidAccessors( mesh.getElemManager(), flowSolver.getName() );

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls const & wellControls = getWellControls( subRegion );
      PerforationData const & perforationData = *subRegion.getPerforationData();

      // get well primary variables on well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView1d< real64 > const & wellElemTemp =
        subRegion.getExtrinsicData< extrinsicMeshData::well::temperature >();
      arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
      arrayView1d< real64 > const & connRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >();

      // get the info stored on well elements
      arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();
      arrayView1d< real64 const > const & wellElemGravCoef =
        subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

      arrayView1d< real64 const > const & perfGravCoef =
        perforationData.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();


      // 1) Loop over all perforations to compute an average mixture density and component fraction
      // 2) Initialize the reference pressure
      // 3) Estimate the pressures in the well elements using the average density
      PresTempCompFracInitializationKernel::
        launch( perforationData.size(),
                subRegion.size(),
                numComp,
                numPhase,
                perforationData.getNumPerforationsGlobal(),
                wellControls,
                0.0, // initialization done at t = 0
                resCompFlowAccessors.get( extrinsicMeshData::flow::pressure{} ),
                resCompFlowAccessors.get( extrinsicMeshData::flow::temperature{} ),
                resCompFlowAccessors.get( extrinsicMeshData::flow::globalCompDensity{} ),
                resCompFlowAccessors.get( extrinsicMeshData::flow::phaseVolumeFraction{} ),
                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseMassDensity{} ),
                resElementRegion,
                resElementSubRegion,
                resElementIndex,
                perfGravCoef,
                wellElemGravCoef,
                wellElemPressure,
                wellElemTemp,
                wellElemCompFrac );

      // get well secondary variables on well elements
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );
      arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens = fluid.phaseDensity();
      arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens = fluid.totalDensity();

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

      CompDensInitializationKernel::launch( subRegion.size(),
                                            numComp,
                                            wellElemCompFrac,
                                            wellElemTotalDens,
                                            wellElemCompDens );

      // 5) Recompute the pressure-dependent properties
      updateSubRegionState( subRegion );

      // 6) Estimate the well rates
      // TODO: initialize rates using perforation rates
      compositionalMultiphaseWellKernels::
        RateInitializationKernel::
        launch( subRegion.size(),
                m_targetPhaseIndex,
                wellControls,
                0.0, // initialization done at t = 0
                wellElemPhaseDens,
                wellElemTotalDens,
                connRate );
    } );

  } );
}

void CompositionalMultiphaseWell::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                     real64 const dt,
                                                     DomainPartition const & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // loop over the wells
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      WellControls const & wellControls = getWellControls( subRegion );

      // get a reference to the degree-of-freedom numbers
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< localIndex const > const & nextWellElemIndex =
        subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & connRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >();

      // get the info stored on well elements
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();
      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >();

      isothermalCompositionalMultiphaseBaseKernels::
        KernelLaunchSelector1< FluxKernel >( numFluidComponents(),
                                             subRegion.size(),
                                             dofManager.rankOffset(),
                                             wellControls,
                                             wellElemDofNumber,
                                             nextWellElemIndex,
                                             connRate,
                                             wellElemCompFrac,
                                             dWellElemCompFrac_dCompDens,
                                             dt,
                                             localMatrix,
                                             localRhs );
    } );
  } );
}

void CompositionalMultiphaseWell::assembleAccumulationTerms( DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      // get the degrees of freedom and ghosting info
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      // get the properties on the well element
      arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction >();
      arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dPressure >();
      arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dComp =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dGlobalCompDensity >();

      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >();

      arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction_n >();
      arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseDens_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::phaseDensity_n >();
      arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & wellElemPhaseCompFrac_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::phaseComponentFraction_n >();

      arrayView1d< real64 const > const & wellElemVolume = subRegion.getElementVolume();

      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
      arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens = fluid.phaseDensity();
      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseDens = fluid.dPhaseDensity();
      arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac = fluid.phaseCompFraction();
      arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dWellElemPhaseCompFrac = fluid.dPhaseCompFraction();

      isothermalCompositionalMultiphaseBaseKernels::
        KernelLaunchSelector1< AccumulationKernel >( numFluidComponents(),
                                                     subRegion.size(),
                                                     numFluidPhases(),
                                                     dofManager.rankOffset(),
                                                     wellElemDofNumber,
                                                     wellElemGhostRank,
                                                     wellElemVolume,
                                                     wellElemPhaseVolFrac,
                                                     dWellElemPhaseVolFrac_dPres,
                                                     dWellElemPhaseVolFrac_dComp,
                                                     dWellElemCompFrac_dCompDens,
                                                     wellElemPhaseDens,
                                                     dWellElemPhaseDens,
                                                     wellElemPhaseCompFrac,
                                                     dWellElemPhaseCompFrac,
                                                     wellElemPhaseVolFrac_n,
                                                     wellElemPhaseDens_n,
                                                     wellElemPhaseCompFrac_n,
                                                     localMatrix,
                                                     localRhs );
    } );
  } );
}


void CompositionalMultiphaseWell::assembleVolumeBalanceTerms( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      // get the degrees of freedom and ghosting info
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      // get the properties on the well element
      arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction >();
      arrayView2d< real64 const, compflow::USD_PHASE > const & dWellElemPhaseVolFrac_dPres =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dPressure >();
      arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac_dComp =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dPhaseVolumeFraction_dGlobalCompDensity >();


      arrayView1d< real64 const > const & wellElemVolume =
        subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

      isothermalCompositionalMultiphaseBaseKernels::
        KernelLaunchSelector1< VolumeBalanceKernel >( numFluidComponents(),
                                                      subRegion.size(),
                                                      numFluidPhases(),
                                                      dofManager.rankOffset(),
                                                      wellElemDofNumber,
                                                      wellElemGhostRank,
                                                      wellElemPhaseVolFrac,
                                                      dWellElemPhaseVolFrac_dPres,
                                                      dWellElemPhaseVolFrac_dComp,
                                                      wellElemVolume,
                                                      localMatrix,
                                                      localRhs );
    } );
  } );
}


real64
CompositionalMultiphaseWell::calculateResidualNorm( DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  real64 localResidualNorm = 0;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {


    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      // get the degree of freedom numbers
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & wellElemVolume =
        subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );
      arrayView2d< real64 const, compflow::USD_PHASE > const wellElemPhaseDens_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::phaseDensity_n >();
      arrayView1d< real64 const > const wellElemTotalDens_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::totalDensity_n >();

      WellControls const & wellControls = getWellControls( subRegion );

      ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                            dofManager.rankOffset(),
                                                            subRegion.isLocallyOwned(),
                                                            subRegion.getTopWellElementIndex(),
                                                            m_numComponents,
                                                            numDofPerWellElement(),
                                                            m_targetPhaseIndex,
                                                            wellControls,
                                                            wellElemDofNumber,
                                                            wellElemGhostRank,
                                                            wellElemVolume,
                                                            wellElemPhaseDens_n,
                                                            wellElemTotalDens_n,
                                                            m_currentTime + m_currentDt, // residual normalized with rate of the end of the
                                                            // time
                                                            // interval
                                                            m_currentDt,
                                                            &localResidualNorm );
    } );
  } );
  return sqrt( MpiWrapper::sum( localResidualNorm, MPI_COMM_GEOSX ) );
}

real64
CompositionalMultiphaseWell::scalingForSystemSolution( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  // check if we want to rescale the Newton update
  if( m_maxCompFracChange >= 1.0 && m_maxRelativePresChange >= 1.0 )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  real64 scalingFactor = 1.0;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      // get the degree of freedom numbers on well elements and ghosting info
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();


      real64 const subRegionScalingFactor =
        SolutionScalingKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                                 dofManager.rankOffset(),
                                                                 numFluidComponents(),
                                                                 wellElemDofNumber,
                                                                 wellElemGhostRank,
                                                                 wellElemPressure,
                                                                 wellElemCompDens,
                                                                 m_maxRelativePresChange,
                                                                 m_maxCompFracChange );


      if( subRegionScalingFactor < scalingFactor )
      {
        scalingFactor = subRegionScalingFactor;
      }
    } );

  } );

  return LvArray::math::max( MpiWrapper::min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
}

bool
CompositionalMultiphaseWell::checkSystemSolution( DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  int localCheck = 1;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {
      // get the degree of freedom numbers on well elements and ghosting info
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();

      localIndex const subRegionSolutionCheck =
        SolutionCheckKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                               dofManager.rankOffset(),
                                                               numFluidComponents(),
                                                               wellElemDofNumber,
                                                               wellElemGhostRank,
                                                               wellElemPressure,
                                                               wellElemCompDens,
                                                               m_allowCompDensChopping,
                                                               scalingFactor );

      if( subRegionSolutionCheck == 0 )
      {
        localCheck = 0;
      }
    } );
  } );
  return MpiWrapper::min( localCheck );
}

void CompositionalMultiphaseWell::computePerforationRates( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    // TODO: change the way we access the flowSolver here
    CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
    PerforationKernel::CompFlowAccessors resCompFlowAccessors( mesh.getElemManager(), flowSolver.getName() );
    PerforationKernel::MultiFluidAccessors resMultiFluidAccessors( mesh.getElemManager(), flowSolver.getName() );
    PerforationKernel::RelPermAccessors resRelPermAccessors( mesh.getElemManager(), flowSolver.getName() );

    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {

      WellControls const & wellControls = getWellControls( subRegion );
      bool const disableReservoirToWellFlow = wellControls.isInjector() and !wellControls.isCrossflowEnabled();

      PerforationData * const perforationData = subRegion.getPerforationData();

      // get depth
      arrayView1d< real64 const > const & wellElemGravCoef =
        subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

      // get well primary variables on well elements
      arrayView1d< real64 const > const & wellElemPres =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();

      arrayView1d< real64 const > const & wellElemTotalMassDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::totalMassDensity >();
      arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dPressure >();
      arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dGlobalCompDensity >();

      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();
      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >();

      // get well variables on perforations
      arrayView1d< real64 const > const & perfGravCoef =
        perforationData->getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();
      arrayView1d< localIndex const > const & perfWellElemIndex =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );
      arrayView1d< real64 const > const & perfTrans =
        perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::wellTransmissibilityString() );

      arrayView2d< real64 > const & compPerfRate =
        perforationData->getExtrinsicData< extrinsicMeshData::well::compPerforationRate >();
      arrayView3d< real64 > const & dCompPerfRate_dPres =
        perforationData->getExtrinsicData< extrinsicMeshData::well::dCompPerforationRate_dPres >();
      arrayView4d< real64 > const & dCompPerfRate_dComp =
        perforationData->getExtrinsicData< extrinsicMeshData::well::dCompPerforationRate_dComp >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

      isothermalCompositionalMultiphaseBaseKernels::
        KernelLaunchSelector2< PerforationKernel >( numFluidComponents(),
                                                    numFluidPhases(),
                                                    perforationData->size(),
                                                    disableReservoirToWellFlow,
                                                    resCompFlowAccessors.get( extrinsicMeshData::flow::pressure{} ),
                                                    resCompFlowAccessors.get( extrinsicMeshData::flow::phaseVolumeFraction{} ),
                                                    resCompFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure{} ),
                                                    resCompFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity{} ),
                                                    resCompFlowAccessors.get( extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                                    resMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseDensity{} ),
                                                    resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseDensity{} ),
                                                    resMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseViscosity{} ),
                                                    resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseViscosity{} ),
                                                    resMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseCompFraction{} ),
                                                    resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseCompFraction{} ),
                                                    resRelPermAccessors.get( extrinsicMeshData::relperm::phaseRelPerm{} ),
                                                    resRelPermAccessors.get( extrinsicMeshData::relperm::dPhaseRelPerm_dPhaseVolFraction{} ),
                                                    wellElemGravCoef,
                                                    wellElemPres,
                                                    wellElemCompDens,
                                                    wellElemTotalMassDens,
                                                    dWellElemTotalMassDens_dPres,
                                                    dWellElemTotalMassDens_dCompDens,
                                                    wellElemCompFrac,
                                                    dWellElemCompFrac_dCompDens,
                                                    perfGravCoef,
                                                    perfWellElemIndex,
                                                    perfTrans,
                                                    resElementRegion,
                                                    resElementSubRegion,
                                                    resElementIndex,
                                                    compPerfRate,
                                                    dCompPerfRate_dPres,
                                                    dCompPerfRate_dComp );

    } );
  } );

}


void
CompositionalMultiphaseWell::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  // update all the fields using the global damping coefficients
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               extrinsicMeshData::well::pressure::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 0, 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               extrinsicMeshData::well::globalCompDensity::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 1, m_numDofPerWellElement - 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               extrinsicMeshData::well::mixtureConnectionRate::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, m_numDofPerWellElement - 1, m_numDofPerWellElement } );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  // synchronize
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::well::pressure::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::well::globalCompDensity::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::well::mixtureConnectionRate::key() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       domain.getMeshBody( 0 ).getMeshLevel( MeshLevel::groupStructKeys::baseDiscretizationString() ),
                                                       domain.getNeighbors(),
                                                       true );
}

void CompositionalMultiphaseWell::chopNegativeDensities( DomainPartition & domain )
{
  integer const numComp = m_numComponents;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
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


void CompositionalMultiphaseWell::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      // get a reference to the primary variables on well elements
      arrayView1d< real64 > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView1d< real64 const > const & wellElemPressure_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure_n >();
      wellElemPressure.setValues< parallelDevicePolicy<> >( wellElemPressure_n );

      arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemGlobalCompDensity_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity_n >();
      wellElemGlobalCompDensity.setValues< parallelDevicePolicy<> >( wellElemGlobalCompDensity_n );

      arrayView1d< real64 > const & connRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >();
      arrayView1d< real64 const > const & connRate_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate_n >();
      connRate.setValues< parallelDevicePolicy<> >( connRate_n );

      updateSubRegionState( subRegion );
    } );
  } );
}

void CompositionalMultiphaseWell::backupFields( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) const
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // backup some fields used in time derivative approximation

  ElementRegionManager & elemManager = mesh.getElemManager();

  elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 WellElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const wellElemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 const, compflow::USD_PHASE > const wellElemPhaseVolFrac =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction >();

    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
    MultiFluidBase const & fluid = subRegion.getConstitutiveModel< MultiFluidBase >( fluidName );
    arrayView2d< real64 const, multifluid::USD_FLUID > const wellElemTotalDens = fluid.totalDensity();
    arrayView3d< real64 const, multifluid::USD_PHASE > const wellElemPhaseDens = fluid.phaseDensity();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const wellElemPhaseCompFrac = fluid.phaseCompFraction();

    arrayView1d< real64 > const wellElemTotalDens_n =
      subRegion.getExtrinsicData< extrinsicMeshData::well::totalDensity_n >();
    arrayView2d< real64, compflow::USD_PHASE > const wellElemPhaseDens_n =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseDensity_n >();
    arrayView2d< real64, compflow::USD_PHASE > const wellElemPhaseVolFrac_n =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction_n >();
    arrayView3d< real64, compflow::USD_PHASE_COMP > const wellElemPhaseCompFrac_n =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseComponentFraction_n >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( wellElemGhostRank[ei] >= 0 )
      {
        return;
      }

      wellElemTotalDens_n[ei] = wellElemTotalDens[ei][0];

      for( integer ip = 0; ip < numPhase; ++ip )
      {
        wellElemPhaseDens_n[ei][ip] = wellElemPhaseDens[ei][0][ip];
        wellElemPhaseVolFrac_n[ei][ip] = wellElemPhaseVolFrac[ei][ip];

        for( integer ic = 0; ic < numComp; ++ic )
        {
          wellElemPhaseCompFrac_n[ei][ip][ic] = wellElemPhaseCompFrac[ei][0][ip][ic];
        }
      }
    } );
  } );
}

void CompositionalMultiphaseWell::assemblePressureRelations( DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion const & subRegion )
    {

      WellControls & wellControls = getWellControls( subRegion );

      // get the degrees of freedom, depth info, next welem index
      string const wellDofKey = dofManager.getKey( wellElementDofName() );
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );
      arrayView1d< real64 const > const & wellElemGravCoef =
        subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();
      arrayView1d< localIndex const > const & nextWellElemIndex =
        subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

      // get primary variables on well elements
      arrayView1d< real64 const > const & wellElemPres =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();

      // get total mass density on well elements (for potential calculations)
      arrayView1d< real64 const > const & wellElemTotalMassDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::totalMassDensity >();
      arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dPressure >();
      arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens =
        subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dGlobalCompDensity >();


      bool controlHasSwitched = false;
      isothermalCompositionalMultiphaseBaseKernels::
        KernelLaunchSelector1< PressureRelationKernel >( numFluidComponents(),
                                                         subRegion.size(),
                                                         dofManager.rankOffset(),
                                                         subRegion.isLocallyOwned(),
                                                         subRegion.getTopWellElementIndex(),
                                                         m_targetPhaseIndex,
                                                         wellControls,
                                                         m_currentTime + m_currentDt,     // controls evaluated with BHP/rate of the end of
                                                                                          // the
                                                                                          // time interval
                                                         wellElemDofNumber,
                                                         wellElemGravCoef,
                                                         nextWellElemIndex,
                                                         wellElemPres,
                                                         wellElemTotalMassDens,
                                                         dWellElemTotalMassDens_dPres,
                                                         dWellElemTotalMassDens_dCompDens,
                                                         controlHasSwitched,
                                                         localMatrix,
                                                         localRhs );

      if( controlHasSwitched )
      {
        // TODO: move the switch logic into wellControls
        // TODO: implement a more general switch when more then two constraints per well type are allowed

        real64 const timeAtEndOfStep = m_currentTime + m_currentDt;

        if( wellControls.getControl() == WellControls::Control::BHP )
        {
          if( wellControls.isProducer() )
          {
            wellControls.switchToPhaseRateControl( wellControls.getTargetPhaseRate( timeAtEndOfStep ) );
            GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                  << " from BHP constraint to phase volumetric rate constraint" );
          }
          else
          {
            wellControls.switchToTotalRateControl( wellControls.getTargetTotalRate( timeAtEndOfStep ) );
            GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                  << " from BHP constraint to total volumetric rate constraint" );
          }
        }
        else
        {
          wellControls.switchToBHPControl( wellControls.getTargetBHP( timeAtEndOfStep ) );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                << " from rate constraint to BHP constraint" );
        }
      }
    } );
  } );
}

void CompositionalMultiphaseWell::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  WellSolverBase::implicitStepSetup( time_n, dt, domain );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
      arrayView1d< real64 > const & wellElemPressure_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::pressure_n >();
      wellElemPressure_n.setValues< parallelDevicePolicy<> >( wellElemPressure );

      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemGlobalCompDensity =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
      arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity_n >();
      wellElemGlobalCompDensity_n.setValues< parallelDevicePolicy<> >( wellElemGlobalCompDensity );

      arrayView1d< real64 const > const & connRate =
        subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >();
      arrayView1d< real64 > const & connRate_n =
        subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate_n >();
      connRate_n.setValues< parallelDevicePolicy<> >( connRate );

      validateWellConstraints( subRegion );

      updateSubRegionState( subRegion );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseWell, string const &, Group * const )
} // namespace geosx
