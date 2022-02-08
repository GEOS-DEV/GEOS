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
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/PerforationData.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"
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
using namespace CompositionalMultiphaseWellKernels;

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

  this->registerWrapper( viewKeyStruct::relPermNamesString(), &m_relPermModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Names of relative permeability constitutive models to use" );

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
  checkModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString() );

  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
  GEOSX_UNUSED_VAR( flowSolver );

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
  MeshLevel & meshLevel = meshBodies.getGroup< MeshBody >( 0 ).getMeshLevel( 0 );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 1. Set key dimensions of the problem
  // Empty check needed to avoid accessing m_fluidModelNames when running in schema generation mode.
  if( !m_fluidModelNames.empty() )
  {
    MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
    m_numPhases = fluid0.numFluidPhases();
    m_numComponents = fluid0.numFluidComponents();
  }
  m_numDofPerWellElement = m_numComponents + 2; // 1 pressure + NC compositions + 1 connectionRate

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

    subRegion.registerExtrinsicData< extrinsicMeshData::well::pressure >( getName() );
    subRegion.registerExtrinsicData< extrinsicMeshData::well::deltaPressure >( getName() );

    subRegion.registerExtrinsicData< extrinsicMeshData::well::temperature >( getName() );

    // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
    // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.

    subRegion.registerExtrinsicData< extrinsicMeshData::well::globalCompDensity >( getName() ).
      reference().resizeDimension< 1 >( m_numComponents );
    subRegion.registerExtrinsicData< extrinsicMeshData::well::deltaGlobalCompDensity >( getName() ).
      reference().resizeDimension< 1 >( m_numComponents );

    subRegion.registerExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >( getName() );
    subRegion.registerExtrinsicData< extrinsicMeshData::well::deltaMixtureConnectionRate >( getName() );

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

    subRegion.registerExtrinsicData< extrinsicMeshData::well::totalDensityOld >( getName() );
    subRegion.registerExtrinsicData< extrinsicMeshData::well::phaseVolumeFractionOld >( getName() ).
      reference().resizeDimension< 1 >( m_numPhases );
    subRegion.registerExtrinsicData< extrinsicMeshData::well::phaseDensityOld >( getName() ).
      reference().resizeDimension< 1 >( m_numPhases );
    subRegion.registerExtrinsicData< extrinsicMeshData::well::phaseComponentFractionOld >( getName() ).
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

void CompositionalMultiphaseWell::validateConstitutiveModels( MeshLevel const & meshLevel, ConstitutiveManager const & cm ) const
{
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  arrayView1d< string const > const & flowTargetRegionNames = flowSolver.targetRegionNames();
  arrayView1d< string const > const & flowFluidModels = flowSolver.fluidModelNames();
  arrayView1d< string const > const & flowRelPermModels = flowSolver.relPermModelNames();

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {
    // Make a set of unique reservoir region indices the well is connected to
    PerforationData const & perforationData = *subRegion.getPerforationData();
    arrayView1d< localIndex const > const & resElementRegion = perforationData.getMeshElements().m_toElementRegion;
    std::set< string > reservoirRegionNames;
    for( localIndex const ei : resElementRegion )
    {
      reservoirRegionNames.insert( meshLevel.getElemManager().getRegion( ei ).getName() );
    }

    // Check that each well model is compatible with all models in perforated reservoir regions
    MultiFluidBase const & wellFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[targetIndex] );
    RelativePermeabilityBase const & wellRelPerm = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[targetIndex] );
    for( localIndex resTargetIndex = 0; resTargetIndex < flowTargetRegionNames.size(); ++resTargetIndex )
    {
      if( reservoirRegionNames.count( flowTargetRegionNames[resTargetIndex] ) > 0 )
      {
        MultiFluidBase const & resFluid = cm.getConstitutiveRelation< MultiFluidBase >( flowFluidModels[resTargetIndex] );
        compareMultiphaseModels( wellFluid, resFluid );
        compareMulticomponentModels( wellFluid, resFluid );

        RelativePermeabilityBase const & resRelPerm = cm.getConstitutiveRelation< RelativePermeabilityBase >( flowRelPermModels[resTargetIndex] );
        compareMultiphaseModels( wellRelPerm, resRelPerm );
      }
    }
  } );
}

void CompositionalMultiphaseWell::validateInjectionStreams( MeshLevel const & meshLevel ) const
{
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
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
                        "WellControls named " << wellControls.getName() <<
                        ": Invalid injection stream for well " << subRegion.getName(),
                        InputError );
        compFracSum += compFrac;
      }
      GEOSX_THROW_IF( ( compFracSum < 1.0 - std::numeric_limits< real64 >::epsilon() ) ||
                      ( compFracSum > 1.0 + std::numeric_limits< real64 >::epsilon() ),
                      "WellControls named " << wellControls.getName() <<
                      ": Invalid injection stream for well " << subRegion.getName(),
                      InputError );
    }
  } );
}

void CompositionalMultiphaseWell::validateWellConstraints( MeshLevel const & meshLevel, MultiFluidBase const & fluid )
{
  // now that we know we are single-phase, we can check a few things in the constraints
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
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
                    "WellControls named " << wellControls.getName() <<
                    ": Phase " << wellControls.getTargetPhaseName() << " not found for well control " << wellControls.getName(),
                    InputError );

  } );

}

void CompositionalMultiphaseWell::initializePreSubGroups()
{
  WellSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    validateConstitutiveModels( meshLevel, cm );

    MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
    validateModelMapping< MultiFluidBase >( meshLevel.getElemManager(), m_fluidModelNames );
    validateModelMapping< RelativePermeabilityBase >( meshLevel.getElemManager(), m_relPermModelNames );
    validateInjectionStreams( meshLevel );
    validateWellConstraints( meshLevel, fluid0 );
  }
}

void CompositionalMultiphaseWell::initializePostInitialConditionsPreSubGroups()
{
  WellSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    fluid.setMassFlag( m_useMass );
  } );
}

void CompositionalMultiphaseWell::updateComponentFraction( WellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  CompositionalMultiphaseBaseKernels::
    ComponentFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               subRegion );

}

void CompositionalMultiphaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion,
                                                          localIndex const targetIndex )
{
  GEOSX_UNUSED_VAR( targetIndex );

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
  arrayView1d< real64 const > const & dPres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

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
                              dPres,
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
    currentBHP = pres[iwelemRef] + dPres[iwelemRef] + totalMassDens[iwelemRef] * diffGravCoef;
    dCurrentBHP_dPres = 1 + dTotalMassDens_dPres[iwelemRef] * diffGravCoef;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      dCurrentBHP_dCompDens[ic] = dTotalMassDens_dCompDens[iwelemRef][ic] * diffGravCoef;
    }
  } );
}

void CompositionalMultiphaseWell::updateVolRatesForConstraint( WellElementSubRegion & subRegion,
                                                               localIndex const targetIndex )
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
  arrayView1d< real64 const > const & dPres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

  arrayView1d< real64 const > const & temp =
    subRegion.getExtrinsicData< extrinsicMeshData::well::temperature >();

  arrayView1d< real64 const > const & connRate =
    subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >();
  arrayView1d< real64 const > const & dConnRate =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaMixtureConnectionRate >();

  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac =
    subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();
  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens =
    subRegion.getExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >();

  // fluid data

  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens = fluid.totalDensity();
  arrayView2d< real64 const, multifluid::USD_FLUID > const & dTotalDens_dPres = fluid.dTotalDensity_dPressure();
  arrayView3d< real64 const, multifluid::USD_FLUID_DC > const & dTotalDens_dComp = fluid.dTotalDensity_dGlobalCompFraction();

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

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
                                dPres,
                                temp,
                                compFrac,
                                dCompFrac_dCompDens,
                                connRate,
                                dConnRate,
                                totalDens,
                                dTotalDens_dPres,
                                dTotalDens_dComp,
                                phaseDens,
                                dPhaseDens_dPres,
                                dPhaseDens_dComp,
                                phaseFrac,
                                dPhaseFrac_dPres,
                                dPhaseFrac_dComp,
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
        real64 const refPres = pres[iwelemRef] + dPres[iwelemRef];
        fluidWrapper.update( iwelemRef, 0, refPres, temp[iwelemRef], compFrac[iwelemRef] );
      }

      // Step 2: update the total volume rate

      real64 const currentTotalRate = connRate[iwelemRef] + dConnRate[iwelemRef];

      // Step 2.1: compute the inverse of the total density and derivatives

      real64 const totalDensInv = 1.0 / totalDens[iwelemRef][0];
      real64 const dTotalDensInv_dPres = -dTotalDens_dPres[iwelemRef][0] * totalDensInv * totalDensInv;
      stackArray1d< real64, maxNumComp > dTotalDensInv_dCompDens( numComp );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalDensInv_dCompDens[ic] = -dTotalDens_dComp[iwelemRef][0][ic] * totalDensInv * totalDensInv;
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
        real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac_dPres[iwelemRef][0][ip] * phaseDensInv
                                                         - dPhaseDens_dPres[iwelemRef][0][ip] * phaseFracTimesPhaseDensInv * phaseDensInv;


        // Step 3.2: divide the total mass/molar rate by the (phase density * phase fraction) to get the phase volumetric rate
        currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
        dCurrentPhaseVolRate_dPres[ip] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
        dCurrentPhaseVolRate_dRate[ip] = phaseFracTimesPhaseDensInv;
        for( integer ic = 0; ic < numComp; ++ic )
        {
          dCurrentPhaseVolRate_dCompDens[ip][ic] = -phaseFracTimesPhaseDensInv * dPhaseDens_dComp[iwelemRef][0][ip][ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] += dPhaseFrac_dComp[iwelemRef][0][ip][ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] *= currentTotalRate;
        }
        applyChainRuleInPlace( numComp, dCompFrac_dCompDens[iwelemRef], dCurrentPhaseVolRate_dCompDens[ip], work.data() );
      }
    } );
  } );
}


void CompositionalMultiphaseWell::updateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres = subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const & dPres = subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();
  arrayView1d< real64 const > const & temp = subRegion.getExtrinsicData< extrinsicMeshData::well::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const & compFrac =
    subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();

  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    CompositionalMultiphaseBaseKernels::FluidUpdateKernel::launch< ExecPolicy >( subRegion.size(),
                                                                                 fluidWrapper,
                                                                                 pres,
                                                                                 dPres,
                                                                                 temp,
                                                                                 compFrac );
  } );
}

void CompositionalMultiphaseWell::updatePhaseVolumeFraction( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  CompositionalMultiphaseBaseKernels::
    PhaseVolumeFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid );

}

void CompositionalMultiphaseWell::updateTotalMassDensity( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  TotalMassDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               subRegion,
                                               fluid );

}

void CompositionalMultiphaseWell::updateSubRegionState( MeshLevel const & meshLevel, WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  // update properties
  updateComponentFraction( subRegion );

  // update volumetric rates for the well constraints
  // note: this must be called before updateFluidModel
  updateVolRatesForConstraint( subRegion, targetIndex );

  // update densities, phase fractions, phase volume fractions
  updateFluidModel( subRegion, targetIndex );
  updatePhaseVolumeFraction( subRegion, targetIndex );
  updateTotalMassDensity( subRegion, targetIndex );

  // update the current BHP pressure
  updateBHPForConstraint( subRegion, targetIndex );

  // update perforation rates
  computePerforationRates( meshLevel, subRegion, targetIndex );
}

void CompositionalMultiphaseWell::initializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
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

    // TODO: change the way we access the flowSolver here
    CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
    PresTempCompFracInitializationKernel::CompFlowAccessors resCompFlowAccessors( meshLevel.getElemManager(), flowSolver.getName() );
    PresTempCompFracInitializationKernel::MultiFluidAccessors resMultiFluidAccessors( meshLevel.getElemManager(), flowSolver.getName(), flowSolver.targetRegionNames(), flowSolver.fluidModelNames() );

    // 1) Loop over all perforations to compute an average mixture density and component fraction
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    PresTempCompFracInitializationKernel::launch( perforationData.size(),
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
                                                  wellElemGravCoef,
                                                  wellElemPressure,
                                                  wellElemTemp,
                                                  wellElemCompFrac );

    // get well secondary variables on well elements
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens = fluid.phaseDensity();
    arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens = fluid.totalDensity();

    // 4) Back calculate component densities
    constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      CompositionalMultiphaseBaseKernels::FluidUpdateKernel::launch< serialPolicy >( subRegion.size(),
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
    // Note: I am leaving that here because I would like to use the perforationRates (computed in UpdateState)
    //       to better initialize the rates
    updateSubRegionState( meshLevel, subRegion, targetIndex );

    // 6) Estimate the well rates
    // TODO: initialize rates using perforation rates
    CompositionalMultiphaseWellKernels::RateInitializationKernel::launch( subRegion.size(),
                                                                          m_targetPhaseIndex,
                                                                          wellControls,
                                                                          0.0, // initialization done at t = 0
                                                                          wellElemPhaseDens,
                                                                          wellElemTotalDens,
                                                                          connRate );
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

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
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
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaMixtureConnectionRate >();

    // get the info stored on well elements
    arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac =
      subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompFraction >();
    arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::dGlobalCompFraction_dGlobalCompDensity >();

    CompositionalMultiphaseBaseKernels::
      KernelLaunchSelector1< FluxKernel >( numFluidComponents(),
                                           subRegion.size(),
                                           dofManager.rankOffset(),
                                           wellControls,
                                           wellElemDofNumber,
                                           nextWellElemIndex,
                                           connRate,
                                           dConnRate,
                                           wellElemCompFrac,
                                           dWellElemCompFrac_dCompDens,
                                           dt,
                                           localMatrix,
                                           localRhs );
  } );
}

void CompositionalMultiphaseWell::assembleAccumulationTerms( DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
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

    arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFracOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFractionOld >();
    arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseDensOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseDensityOld >();
    arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & wellElemPhaseCompFracOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseComponentFractionOld >();

    arrayView1d< real64 const > const & wellElemVolume = subRegion.getElementVolume();

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens = fluid.phaseDensity();
    arrayView3d< real64 const, multifluid::USD_PHASE > const & dWellElemPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
    arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac = fluid.phaseCompFraction();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dWellElemPhaseCompFrac_dPres = fluid.dPhaseCompFraction_dPressure();
    arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dWellElemPhaseCompFrac_dComp = fluid.dPhaseCompFraction_dGlobalCompFraction();

    CompositionalMultiphaseBaseKernels::
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
                                                   dWellElemPhaseDens_dPres,
                                                   dWellElemPhaseDens_dComp,
                                                   wellElemPhaseCompFrac,
                                                   dWellElemPhaseCompFrac_dPres,
                                                   dWellElemPhaseCompFrac_dComp,
                                                   wellElemPhaseVolFracOld,
                                                   wellElemPhaseDensOld,
                                                   wellElemPhaseCompFracOld,
                                                   localMatrix,
                                                   localRhs );
  } );
}


void CompositionalMultiphaseWell::assembleVolumeBalanceTerms( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const wellDofKey = dofManager.getKey( wellElementDofName() );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
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

    CompositionalMultiphaseBaseKernels::
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
}


real64
CompositionalMultiphaseWell::calculateResidualNorm( DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  real64 localResidualNorm = 0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & wellElemVolume =
      subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const wellElemPhaseDensOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseDensityOld >();
    arrayView1d< real64 const > const wellElemTotalDensOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::totalDensityOld >();

    WellControls const & wellControls = getWellControls( subRegion );

    ResidualNormKernel::launch< parallelDevicePolicy<>,
                                parallelDeviceReduce >( localRhs,
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
                                                        wellElemPhaseDensOld,
                                                        wellElemTotalDensOld,
                                                        m_currentTime + m_currentDt, // residual normalized with rate of the end of the time
                                                                                     // interval
                                                        m_currentDt,
                                                        &localResidualNorm );
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

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  real64 scalingFactor = 1.0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
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
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

    arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
    arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaGlobalCompDensity >();

    real64 const subRegionScalingFactor =
      SolutionScalingKernel::launch< parallelDevicePolicy<>,
                                     parallelDeviceReduce >( localSolution,
                                                             dofManager.rankOffset(),
                                                             numFluidComponents(),
                                                             wellElemDofNumber,
                                                             wellElemGhostRank,
                                                             wellElemPressure,
                                                             dWellElemPressure,
                                                             wellElemCompDens,
                                                             dWellElemCompDens,
                                                             m_maxRelativePresChange,
                                                             m_maxCompFracChange );


    if( subRegionScalingFactor < scalingFactor )
    {
      scalingFactor = subRegionScalingFactor;
    }
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

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  int localCheck = 1;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
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
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

    arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
    arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaGlobalCompDensity >();

    localIndex const subRegionSolutionCheck =
      SolutionCheckKernel::launch< parallelDevicePolicy<>,
                                   parallelDeviceReduce >( localSolution,
                                                           dofManager.rankOffset(),
                                                           numFluidComponents(),
                                                           wellElemDofNumber,
                                                           wellElemGhostRank,
                                                           wellElemPressure,
                                                           dWellElemPressure,
                                                           wellElemCompDens,
                                                           dWellElemCompDens,
                                                           m_allowCompDensChopping,
                                                           scalingFactor );

    if( subRegionSolutionCheck == 0 )
    {
      localCheck = 0;
    }
  } );
  return MpiWrapper::min( localCheck );
}

void CompositionalMultiphaseWell::computePerforationRates( MeshLevel const & meshLevel,
                                                           WellElementSubRegion & subRegion,
                                                           localIndex const GEOSX_UNUSED_PARAM( targetIndex ) )
{
  GEOSX_MARK_FUNCTION;

  // if the well is shut, we neglect reservoir-well flow that may occur despite the zero rate
  // therefore, we do not want to compute perforation rates and we simply assume they are zero
  WellControls const & wellControls = getWellControls( subRegion );
  bool const disableReservoirToWellFlow = wellControls.isInjector() and wellControls.isCrossflowDisabled();
  if( !wellControls.isWellOpen( m_currentTime + m_currentDt ) )
  {
    return;
  }

  PerforationData * const perforationData = subRegion.getPerforationData();

  // get depth
  arrayView1d< real64 const > const & wellElemGravCoef =
    subRegion.getExtrinsicData< extrinsicMeshData::well::gravityCoefficient >();

  // get well primary variables on well elements
  arrayView1d< real64 const > const & wellElemPres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
  arrayView1d< real64 const > const & dWellElemPres =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

  arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens =
    subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
  arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemCompDens =
    subRegion.getExtrinsicData< extrinsicMeshData::well::deltaGlobalCompDensity >();

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

  // TODO: change the way we access the flowSolver here
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
  PerforationKernel::CompFlowAccessors resCompFlowAccessors( meshLevel.getElemManager(), flowSolver.getName() );
  PerforationKernel::MultiFluidAccessors resMultiFluidAccessors( meshLevel.getElemManager(), flowSolver.getName(), flowSolver.targetRegionNames(), flowSolver.fluidModelNames() );
  PerforationKernel::RelPermAccessors resRelPermAccessors( meshLevel.getElemManager(), flowSolver.getName(), flowSolver.targetRegionNames(), flowSolver.relPermModelNames() );

  CompositionalMultiphaseBaseKernels::
    KernelLaunchSelector2< PerforationKernel >( numFluidComponents(),
                                                numFluidPhases(),
                                                perforationData->size(),
                                                disableReservoirToWellFlow,
                                                resCompFlowAccessors.get( extrinsicMeshData::flow::pressure{} ),
                                                resCompFlowAccessors.get( extrinsicMeshData::flow::deltaPressure{} ),
                                                resCompFlowAccessors.get( extrinsicMeshData::flow::phaseVolumeFraction{} ),
                                                resCompFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction_dPressure{} ),
                                                resCompFlowAccessors.get( extrinsicMeshData::flow::dPhaseVolumeFraction_dGlobalCompDensity{} ),
                                                resCompFlowAccessors.get( extrinsicMeshData::flow::dGlobalCompFraction_dGlobalCompDensity{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseDensity{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseDensity_dPressure{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseDensity_dGlobalCompFraction{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseViscosity{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseViscosity_dPressure{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseViscosity_dGlobalCompFraction{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::phaseCompFraction{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseCompFraction_dPressure{} ),
                                                resMultiFluidAccessors.get( extrinsicMeshData::multifluid::dPhaseCompFraction_dGlobalCompFraction{} ),
                                                resRelPermAccessors.get( extrinsicMeshData::relperm::phaseRelPerm{} ),
                                                resRelPermAccessors.get( extrinsicMeshData::relperm::dPhaseRelPerm_dPhaseVolFraction{} ),
                                                wellElemGravCoef,
                                                wellElemPres,
                                                dWellElemPres,
                                                wellElemCompDens,
                                                dWellElemCompDens,
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
                               extrinsicMeshData::well::deltaPressure::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 0, 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               extrinsicMeshData::well::deltaGlobalCompDensity::key(),
                               scalingFactor,
                               { m_numDofPerWellElement, 1, m_numDofPerWellElement - 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               extrinsicMeshData::well::deltaMixtureConnectionRate::key(),
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
  fieldNames["elems"].emplace_back( extrinsicMeshData::well::deltaPressure::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::well::deltaGlobalCompDensity::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::well::deltaMixtureConnectionRate::key() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       true );
}

void CompositionalMultiphaseWell::chopNegativeDensities( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  integer const numComp = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
    arrayView2d< real64, compflow::USD_COMP > const & dWellElemCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaGlobalCompDensity >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        for( integer ic = 0; ic < numComp; ++ic )
        {
          // get the latest component density (i.e., after the update of dWellElemCompDens)
          real64 const newDens = wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic];
          // we allowed for some densities to be slightly negative in CheckSystemSolution
          // if the new density is negative, chop back to zero
          if( newDens < 0 )
          {
            dWellElemCompDens[iwelem][ic] = -wellElemCompDens[iwelem][ic];
          }
        }
      }
    } );
  } );
}


void CompositionalMultiphaseWell::resetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  integer const numComp = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const & dWellElemPressure =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();
    arrayView2d< real64, compflow::USD_COMP > const & dWellElemGlobalCompDensity =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaGlobalCompDensity >();
    arrayView1d< real64 > const & dConnRate =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaMixtureConnectionRate >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      // extract solution and apply to dP
      dWellElemPressure[iwelem] = 0;
      dConnRate[iwelem] = 0;
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dWellElemGlobalCompDensity[iwelem][ic] = 0;
      }
    } );
  } );

  // call constitutive models
  updateState( domain );
}

void CompositionalMultiphaseWell::backupFields( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // backup some fields used in time derivative approximation
  forTargetSubRegions< WellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                          WellElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const wellElemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 const, compflow::USD_PHASE > const wellElemPhaseVolFrac =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFraction >();

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView2d< real64 const, multifluid::USD_FLUID > const wellElemTotalDens = fluid.totalDensity();
    arrayView3d< real64 const, multifluid::USD_PHASE > const wellElemPhaseDens = fluid.phaseDensity();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const wellElemPhaseCompFrac = fluid.phaseCompFraction();

    arrayView1d< real64 > const wellElemTotalDensOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::totalDensityOld >();
    arrayView2d< real64, compflow::USD_PHASE > const wellElemPhaseDensOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseDensityOld >();
    arrayView2d< real64, compflow::USD_PHASE > const wellElemPhaseVolFracOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseVolumeFractionOld >();
    arrayView3d< real64, compflow::USD_PHASE_COMP > const wellElemPhaseCompFracOld =
      subRegion.getExtrinsicData< extrinsicMeshData::well::phaseComponentFractionOld >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( wellElemGhostRank[ei] >= 0 )
      {
        return;
      }

      wellElemTotalDensOld[ei] = wellElemTotalDens[ei][0];

      for( integer ip = 0; ip < numPhase; ++ip )
      {
        wellElemPhaseDensOld[ei][ip] = wellElemPhaseDens[ei][0][ip];
        wellElemPhaseVolFracOld[ei][ip] = wellElemPhaseVolFrac[ei][ip];

        for( integer ic = 0; ic < numComp; ++ic )
        {
          wellElemPhaseCompFracOld[ei][ip][ic] = wellElemPhaseCompFrac[ei][0][ip][ic];
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

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
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
    arrayView1d< real64 const > const & dWellElemPres =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

    // get total mass density on well elements (for potential calculations)
    arrayView1d< real64 const > const & wellElemTotalMassDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::totalMassDensity >();
    arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres =
      subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dPressure >();
    arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::well::dTotalMassDensity_dGlobalCompDensity >();


    bool controlHasSwitched = false;
    CompositionalMultiphaseBaseKernels::
      KernelLaunchSelector1< PressureRelationKernel >( numFluidComponents(),
                                                       subRegion.size(),
                                                       dofManager.rankOffset(),
                                                       subRegion.isLocallyOwned(),
                                                       subRegion.getTopWellElementIndex(),
                                                       m_targetPhaseIndex,
                                                       wellControls,
                                                       m_currentTime + m_currentDt, // controls evaluated with BHP/rate of the end of the
                                                                                    // time interval
                                                       wellElemDofNumber,
                                                       wellElemGravCoef,
                                                       nextWellElemIndex,
                                                       wellElemPres,
                                                       dWellElemPres,
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
}

void CompositionalMultiphaseWell::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  WellSolverBase::implicitStepSetup( time_n, dt, domain );

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
  validateWellConstraints( meshLevel, fluid0 );
}


void CompositionalMultiphaseWell::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  integer const numComp = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure =
      subRegion.getExtrinsicData< extrinsicMeshData::well::pressure >();
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaPressure >();

    arrayView2d< real64, compflow::USD_COMP > const & wellElemGlobalCompDensity =
      subRegion.getExtrinsicData< extrinsicMeshData::well::globalCompDensity >();
    arrayView2d< real64 const, compflow::USD_COMP > const & dWellElemGlobalCompDensity =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaGlobalCompDensity >();

    arrayView1d< real64 > const & connRate =
      subRegion.getExtrinsicData< extrinsicMeshData::well::mixtureConnectionRate >();
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getExtrinsicData< extrinsicMeshData::well::deltaMixtureConnectionRate >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      for( integer ic = 0; ic < numComp; ++ic )
      {
        wellElemGlobalCompDensity[iwelem][ic] += dWellElemGlobalCompDensity[iwelem][ic];
      }
      connRate[iwelem] += dConnRate[iwelem];
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseWell, string const &, Group * const )
}// namespace geosx
