/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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

#include "mpiCommunications/CommunicationTools.hpp"
#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/PerforationData.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
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
  m_temperature( 0.0 ),
  m_useMass( false ),
  m_maxCompFracChange( 1.0 ),
  m_maxRelativePresChange( 0.2 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 ),
  m_oilPhaseIndex( -1 )
{
  this->registerWrapper( viewKeyStruct::temperatureString, &m_temperature )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Temperature" );

  this->registerWrapper( viewKeyStruct::useMassFlagString, &m_useMass )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::relPermNamesString, &m_relPermModelNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Names of relative permeability constitutive models to use" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString, &m_maxCompFracChange )->
    setSizedFromParent( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setApplyDefaultValue( 1.0 )->
    setDescription( "Maximum (absolute) change in a component fraction between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString, &m_maxRelativePresChange )->
    setSizedFromParent( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setApplyDefaultValue( 1.0 )->
    setDescription( "Maximum (relative) change in pressure between two Newton iterations (recommended with rate control)" );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString, &m_allowCompDensChopping )->
    setSizedFromParent( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setApplyDefaultValue( 1 )->
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );

}

void CompositionalMultiphaseWell::PostProcessInput()
{
  WellSolverBase::PostProcessInput();
  CheckModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString );

  CompositionalMultiphaseFlow const * const flowSolver = getParent()->GetGroup< CompositionalMultiphaseFlow >( GetFlowSolverName() );
  GEOSX_ERROR_IF( flowSolver == nullptr,
                  "Flow solver " << GetFlowSolverName() << " not found or incompatible type "
                                                           "(referenced from well solver " << getName() << ")" );

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         "The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         "The maximum absolute change in component fraction must larger or equal to 0.0" );

}

void CompositionalMultiphaseWell::RegisterDataOnMesh( Group * const meshBodies )
{
  WellSolverBase::RegisterDataOnMesh( meshBodies );

  MeshLevel & meshLevel = *meshBodies->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompDensityString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::mixtureConnRateString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompFractionString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString )->
      setRestartFlags( RestartFlags::NO_WRITE );
    subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::totalMassDensityString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString )->
      setRestartFlags( RestartFlags::NO_WRITE );
    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    PerforationData & perforationData = *subRegion.GetPerforationData();
    perforationData.registerWrapper< array2d< real64 > >( viewKeyStruct::compPerforationRateString );
    perforationData.registerWrapper< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString )->
      setRestartFlags( RestartFlags::NO_WRITE );
    perforationData.registerWrapper< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    WellControls & wellControls = GetWellControls( subRegion );
    wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentBHP_dPresString )->
      setRestartFlags( RestartFlags::NO_WRITE );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentBHP_dCompDensString )->
      setRestartFlags( RestartFlags::NO_WRITE )->
      setSizedFromParent( 0 );

    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString )->
      setSizedFromParent( 0 );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dPresString )->
      setRestartFlags( RestartFlags::NO_WRITE )->
      setSizedFromParent( 0 );
    wellControls.registerWrapper< array2d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dCompDensString )->
      setRestartFlags( RestartFlags::NO_WRITE )->
      setSizedFromParent( 0 );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dRateString )->
      setRestartFlags( RestartFlags::NO_WRITE )->
      setSizedFromParent( 0 );

    wellControls.registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentTotalVolRate_dPresString )->
      setRestartFlags( RestartFlags::NO_WRITE );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentTotalVolRate_dCompDensString )->
      setRestartFlags( RestartFlags::NO_WRITE )->
      setSizedFromParent( 0 );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentTotalVolRate_dRateString )->
      setRestartFlags( RestartFlags::NO_WRITE );
  } );

}

namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void CompareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_ERROR_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                         "Mismatch in number of phases between constitutive models "
                         << lhs.getName() << " and " << rhs.getName() );

  for( localIndex ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOSX_ERROR_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                           "Mismatch in phase names between constitutive models "
                           << lhs.getName() << " and " << rhs.getName() );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void CompareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_ERROR_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                         "Mismatch in number of components between constitutive models "
                         << lhs.getName() << " and " << rhs.getName() );

  for( localIndex ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOSX_ERROR_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                           "Mismatch in component names between constitutive models "
                           << lhs.getName() << " and " << rhs.getName() );
  }
}

}

void CompositionalMultiphaseWell::ValidateConstitutiveModels( MeshLevel const & meshLevel, ConstitutiveManager const & cm ) const
{
  CompositionalMultiphaseFlow const & flowSolver = *getParent()->GetGroup< CompositionalMultiphaseFlow >( GetFlowSolverName() );
  arrayView1d< string const > const & flowTargetRegionNames = flowSolver.targetRegionNames();
  arrayView1d< string const > const & flowFluidModels = flowSolver.fluidModelNames();
  arrayView1d< string const > const & flowRelPermModels = flowSolver.relPermModelNames();

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {
    // Make a set of unique reservoir region indices the well is connected to
    PerforationData const & perforationData = *subRegion.GetPerforationData();
    arrayView1d< localIndex const > const & resElementRegion = perforationData.GetMeshElements().m_toElementRegion;
    std::set< string > reservoirRegionNames;
    for( localIndex const ei : resElementRegion )
    {
      reservoirRegionNames.insert( meshLevel.getElemManager()->GetRegion( ei )->getName() );
    }

    // Check that each well model is compatible with all models in perforated reservoir regions
    MultiFluidBase const & wellFluid = *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[targetIndex] );
    RelativePermeabilityBase const & wellRelPerm = *cm.GetConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[targetIndex] );
    for( localIndex resTargetIndex = 0; resTargetIndex < flowTargetRegionNames.size(); ++resTargetIndex )
    {
      if( reservoirRegionNames.count( flowTargetRegionNames[resTargetIndex] ) > 0 )
      {
        MultiFluidBase const & resFluid = *cm.GetConstitutiveRelation< MultiFluidBase >( flowFluidModels[resTargetIndex] );
        CompareMultiphaseModels( wellFluid, resFluid );
        CompareMulticomponentModels( wellFluid, resFluid );

        RelativePermeabilityBase const & resRelPerm = *cm.GetConstitutiveRelation< RelativePermeabilityBase >( flowRelPermModels[resTargetIndex] );
        CompareMultiphaseModels( wellRelPerm, resRelPerm );
      }
    }
  } );
}

void CompositionalMultiphaseWell::ValidateInjectionStreams( MeshLevel const & meshLevel ) const
{
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    WellControls const & wellControls = GetWellControls( subRegion );

    // check well injection stream for injectors
    if( wellControls.GetType() == WellControls::Type::INJECTOR )
    {
      arrayView1d< real64 const > const & injection = wellControls.GetInjectionStream();
      real64 compFracSum = 0;
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const compFrac = injection[ic];
        GEOSX_ERROR_IF( ( compFrac < 0.0 ) || ( compFrac > 1.0 ),
                        "Invalid injection stream for well " << subRegion.getName() );
        compFracSum += compFrac;
      }
      GEOSX_ERROR_IF( ( compFracSum < 1.0 - std::numeric_limits< real64 >::epsilon() ) ||
                      ( compFracSum > 1.0 + std::numeric_limits< real64 >::epsilon() ),
                      "Invalid injection stream for well " << subRegion.getName() );
    }
  } );
}

void CompositionalMultiphaseWell::InitializePreSubGroups( Group * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ConstitutiveManager const & cm = *domain->getConstitutiveManager();

  ValidateConstitutiveModels( meshLevel, cm );

  // Set key dimensions (phases, components) from one of the fluids - they should all be compatible
  MultiFluidBase const & fluid0 = *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
  m_numPhases     = fluid0.numFluidPhases();
  m_numComponents = fluid0.numFluidComponents();
  m_numDofPerWellElement = m_numComponents + 2; // 1 pressure + NC compositions + 1 connectionRate

  ValidateModelMapping< MultiFluidBase >( *meshLevel.getElemManager(), m_fluidModelNames );
  ValidateModelMapping< RelativePermeabilityBase >( *meshLevel.getElemManager(), m_relPermModelNames );
  ValidateInjectionStreams( meshLevel );

  // Find oil phase index for oil rate constraint
  for( localIndex ip = 0; ip < fluid0.numFluidPhases(); ++ip )
  {
    if( fluid0.phaseNames()[ip] == "Oil" || fluid0.phaseNames()[ip] == "oil" )
    {
      m_oilPhaseIndex = ip;
    }
  }

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    ResizeFields( subRegion );
  } );
}

void CompositionalMultiphaseWell::ResizeFields( WellElementSubRegion & subRegion )
{
  PerforationData & perforationData = *subRegion.GetPerforationData();
  WellControls & wellControls = GetWellControls( subRegion );

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString ).resizeDimension< 1 >( NC );
  subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString ).resizeDimension< 1 >( NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString ).resizeDimension< 1 >( NC );
  subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString ).resizeDimension< 1, 2 >( NC, NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString ).resizeDimension< 1 >( NP );
  subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString ).resizeDimension< 1 >( NP );
  subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString ).resizeDimension< 1, 2 >( NP, NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString ).resizeDimension< 1 >( NC );

  perforationData.getReference< array2d< real64 > >( viewKeyStruct::compPerforationRateString ).resizeDimension< 1 >( NC );
  perforationData.getReference< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString ).resizeDimension< 1, 2 >( 2, NC );
  perforationData.getReference< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString ).resizeDimension< 1, 2, 3 >( 2, NC, NC );

  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentBHP_dCompDensString ).resizeDimension< 0 >( NC );

  wellControls.getReference< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString ).resizeDimension< 0 >( NP );
  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dPresString ).resizeDimension< 0 >( NP );
  wellControls.getReference< array2d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dCompDensString ).resizeDimension< 0, 1 >( NP, NC );
  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dRateString ).resizeDimension< 0 >( NP );

  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentTotalVolRate_dCompDensString ).resizeDimension< 0 >( NC );
}

void CompositionalMultiphaseWell::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    fluid.setMassFlag( m_useMass );
  } );
}

void CompositionalMultiphaseWell::UpdateComponentFraction( WellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs
  arrayView2d< real64 > const & compFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
  arrayView3d< real64 > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // inputs
  arrayView2d< real64 const > const & compDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
  arrayView2d< real64 const > const & dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  CompositionalMultiphaseFlowKernels::KernelLaunchSelector1< CompositionalMultiphaseFlowKernels::ComponentFractionKernel
                                                             >( NumFluidComponents(),
                                                                subRegion.size(),
                                                                compDens,
                                                                dCompDens,
                                                                compFrac,
                                                                dCompFrac_dCompDens );

}


void CompositionalMultiphaseWell::UpdateBHPForConstraint( WellElementSubRegion & subRegion,
                                                          localIndex const targetIndex )
{
  GEOSX_UNUSED_VAR( targetIndex );

  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.IsLocallyOwned() )
  {
    return;
  }

  localIndex const NC = m_numComponents;
  localIndex const iwelemRef = subRegion.GetTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & pres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView1d< real64 > const & totalMassDens =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString );
  arrayView1d< real64 > const & dTotalMassDens_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString );
  arrayView2d< real64 > const & dTotalMassDens_dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString );

  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // control data

  WellControls & wellControls = GetWellControls( subRegion );

  real64 const & refGravCoef = wellControls.GetReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString );
  real64 & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dPresString );
  arrayView1d< real64 > const & dCurrentBHP_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dCompDensString );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [&NC,
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
                              &refGravCoef] GEOSX_HOST_DEVICE ( localIndex const )
  {
    real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
    currentBHP = pres[iwelemRef] + dPres[iwelemRef] + totalMassDens[iwelemRef] * diffGravCoef;
    dCurrentBHP_dPres = 1 + dTotalMassDens_dPres[iwelemRef] * diffGravCoef;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dCurrentBHP_dCompDens[ic] = dTotalMassDens_dCompDens[iwelemRef][ic] * diffGravCoef;
    }
  } );
}

void CompositionalMultiphaseWell::UpdateVolRatesForConstraint( WellElementSubRegion & subRegion,
                                                               localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.IsLocallyOwned() )
  {
    return;
  }

  localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;
  localIndex const iwelemRef = subRegion.GetTopWellElementIndex();
  real64 const temp = m_temperature;

  // subRegion data

  arrayView1d< real64 const > const & pres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView2d< real64 const > const & compDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
  arrayView2d< real64 const > const & dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  arrayView1d< real64 const > const & connRate =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );
  arrayView1d< real64 const > const & dConnRate =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString );

  arrayView2d< real64 const > const & compFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // fluid data

  MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  //arrayView2d< real64 const > const & totalDens = fluid.totalDensity();
  //arrayView2d< real64 const > const & dTotalDens_dPres = fluid.dTotalDensity_dPressure();

  arrayView3d< real64 const > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  // control data

  WellControls & wellControls = GetWellControls( subRegion );

  integer const useSurfaceConditions = wellControls.UseSurfaceConditions();

  arrayView1d< real64 > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dPres =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dPresString );
  arrayView2d< real64 > const & dCurrentPhaseVolRate_dCompDens =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dCompDensString );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dRateString );

  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString );
  real64 & dCurrentTotalVolRate_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dPresString );
  arrayView1d< real64 > const & dCurrentTotalVolRate_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dCompDensString );
  real64 & dCurrentTotalVolRate_dRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dRateString );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [&NC,
                                &NP,
                                fluidWrapper,
                                pres,
                                dPres,
                                compDens,
                                dCompDens,
                                compFrac,
                                dCompFrac_dCompDens,
                                connRate,
                                dConnRate,
                                phaseDens,
                                dPhaseDens_dPres,
                                dPhaseDens_dComp,
                                phaseFrac,
                                dPhaseFrac_dPres,
                                dPhaseFrac_dComp,
                                &useSurfaceConditions,
                                &temp,
                                &currentTotalVolRate,
                                &dCurrentTotalVolRate_dPres,
                                &dCurrentTotalVolRate_dCompDens,
                                &dCurrentTotalVolRate_dRate,
                                currentPhaseVolRate,
                                dCurrentPhaseVolRate_dPres,
                                dCurrentPhaseVolRate_dCompDens,
                                dCurrentPhaseVolRate_dRate,
                                &iwelemRef] GEOSX_HOST_DEVICE ( localIndex const )
    {
      //    We need to evaluate the density as follows:
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element

      if( useSurfaceConditions )
      {
        real64 const surfacePres = 101325.0;
        real64 const surfaceTemp = 293.15;
        // we need to compute the surface density
        fluidWrapper.Update( iwelemRef, 0, surfacePres, surfaceTemp, compFrac[iwelemRef] );
      }
      else
      {
        real64 const refPres = pres[iwelemRef] + dPres[iwelemRef];
        fluidWrapper.Update( iwelemRef, 0, refPres, temp, compFrac[iwelemRef] );
      }

      // total volume rate
      real64 const currentTotalRate = connRate[iwelemRef] + dConnRate[iwelemRef];
      real64 totalDens = 0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        totalDens += compDens[iwelemRef][ic] + dCompDens[iwelemRef][ic];
      }
      real64 const totalDensInv = 1.0 / totalDens;
      currentTotalVolRate = currentTotalRate * totalDensInv;
      dCurrentTotalVolRate_dPres = 0; // using the fact that totalDens = \sum_c compDens (independent of pressure)
      dCurrentTotalVolRate_dRate = totalDensInv;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dCurrentTotalVolRate_dCompDens[ic] = -currentTotalVolRate * totalDensInv; // using the fact that totalDens = \sum_c compDens
      }

      // phase volume rate
      stackArray1d< real64, maxNumComp > work( NC );
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        real64 const phaseDensInv = 1.0 / phaseDens[iwelemRef][0][ip];
        real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;
        real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac_dPres[iwelemRef][0][ip] * phaseDensInv
                                                         - dPhaseDens_dPres[iwelemRef][0][ip] * phaseFracTimesPhaseDensInv * phaseDensInv;

        currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
        dCurrentPhaseVolRate_dPres[ip] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
        dCurrentPhaseVolRate_dRate[ip] = phaseFracTimesPhaseDensInv;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dCurrentPhaseVolRate_dCompDens[ip][ic] = -phaseFracTimesPhaseDensInv * dPhaseDens_dComp[iwelemRef][0][ip][ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] += dPhaseFrac_dComp[iwelemRef][0][ip][ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] *= currentTotalRate;
        }
        applyChainRuleInPlace( NC, dCompFrac_dCompDens[iwelemRef], dCurrentPhaseVolRate_dCompDens[ip], work.data() );
      }
    } );
  } );
}


void CompositionalMultiphaseWell::UpdateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
  arrayView2d< real64 const > const & compFrac = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );

  MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    CompositionalMultiphaseFlowKernels::FluidUpdateKernel::Launch< serialPolicy >( subRegion.size(),
                                                                                   fluidWrapper,
                                                                                   pres,
                                                                                   dPres,
                                                                                   m_temperature,
                                                                                   compFrac );
  } );
}

void CompositionalMultiphaseWell::UpdatePhaseVolumeFraction( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseVolFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );
  arrayView2d< real64 > const & dPhaseVolFrac_dPres =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  arrayView3d< real64 > const & dPhaseVolFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  // inputs

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
  arrayView2d< real64 const > const & compDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
  arrayView2d< real64 const > const & dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  CompositionalMultiphaseFlowKernels::KernelLaunchSelector2< CompositionalMultiphaseFlowKernels::PhaseVolumeFractionKernel
                                                             >( NumFluidComponents(), NumFluidPhases(),
                                                                subRegion.size(),
                                                                compDens,
                                                                dCompDens,
                                                                dCompFrac_dCompDens,
                                                                phaseDens,
                                                                dPhaseDens_dPres,
                                                                dPhaseDens_dComp,
                                                                phaseFrac,
                                                                dPhaseFrac_dPres,
                                                                dPhaseFrac_dComp,
                                                                phaseVolFrac,
                                                                dPhaseVolFrac_dPres,
                                                                dPhaseVolFrac_dCompDens );
}

void CompositionalMultiphaseWell::UpdateTotalMassDensity( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  // outputs

  arrayView1d< real64 > const & totalMassDens =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString );
  arrayView1d< real64 > const & dTotalMassDens_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString );
  arrayView2d< real64 > const & dTotalMassDens_dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const & phaseVolFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );
  arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseMassDens = fluid.phaseMassDensity();
  arrayView3d< real64 const > const & dPhaseMassDens_dPres = fluid.dPhaseMassDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseMassDens_dComp = fluid.dPhaseMassDensity_dGlobalCompFraction();

  TotalMassDensityKernel::Launch< parallelDevicePolicy<> >( subRegion.size(),
                                                            NumFluidComponents(),
                                                            NumFluidPhases(),
                                                            phaseVolFrac,
                                                            dPhaseVolFrac_dPres,
                                                            dPhaseVolFrac_dCompDens,
                                                            dCompFrac_dCompDens,
                                                            phaseMassDens,
                                                            dPhaseMassDens_dPres,
                                                            dPhaseMassDens_dComp,
                                                            totalMassDens,
                                                            dTotalMassDens_dPres,
                                                            dTotalMassDens_dCompDens );

}

void CompositionalMultiphaseWell::UpdateState( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  // update properties
  UpdateComponentFraction( subRegion );

  // update volumetric rates for the well constraints
  // note: this must be called before UpdateFluidModel
  UpdateVolRatesForConstraint( subRegion, targetIndex );

  // update densities, phase fractions, phase volume fractions
  UpdateFluidModel( subRegion, targetIndex );
  UpdatePhaseVolumeFraction( subRegion, targetIndex );
  UpdateTotalMassDensity( subRegion, targetIndex );

  // update the current BHP pressure
  UpdateBHPForConstraint( subRegion, targetIndex );

  // update perforation rates
  ComputePerforationRates( subRegion, targetIndex );
}

void CompositionalMultiphaseWell::InitializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {

    WellControls const & wellControls = GetWellControls( subRegion );
    PerforationData const & perforationData = *subRegion.GetPerforationData();

    // get well primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView2d< real64 > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView1d< real64 > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );

    // get the info stored on well elements
    arrayView2d< real64 > const & wellElemCompFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // 1) Loop over all perforations to compute an average mixture density and component fraction
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    PresCompFracInitializationKernel::Launch< parallelDevicePolicy<> >( perforationData.size(),
                                                                        subRegion.size(),
                                                                        NC,
                                                                        NP,
                                                                        perforationData.GetNumPerforationsGlobal(),
                                                                        wellControls,
                                                                        m_resPres.toNestedViewConst(),
                                                                        m_resCompDens.toNestedViewConst(),
                                                                        m_resPhaseVolFrac.toNestedViewConst(),
                                                                        m_resPhaseMassDens.toNestedViewConst(),
                                                                        resElementRegion,
                                                                        resElementSubRegion,
                                                                        resElementIndex,
                                                                        wellElemGravCoef,
                                                                        wellElemPressure,
                                                                        wellElemCompFrac );

    // get well secondary variables on well elements
    MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView3d< real64 const > const & wellElemPhaseDens = fluid.phaseDensity();
    arrayView2d< real64 const > const & wellElemTotalDens = fluid.totalDensity();


    // 4) Back calculate component densities
    constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      CompositionalMultiphaseFlowKernels::FluidUpdateKernel::Launch< serialPolicy >( subRegion.size(),
                                                                                     fluidWrapper,
                                                                                     wellElemPressure,
                                                                                     m_temperature,
                                                                                     wellElemCompFrac );
    } );

    CompDensInitializationKernel::Launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                    NC,
                                                                    wellElemCompFrac,
                                                                    wellElemTotalDens,
                                                                    wellElemCompDens );

    // 5) Recompute the pressure-dependent properties
    // Note: I am leaving that here because I would like to use the perforationRates (computed in UpdateState)
    //       to better initialize the rates
    UpdateState( subRegion, targetIndex );

    // 6) Estimate the well rates
    // TODO: initialize rates using perforation rates
    CompositionalMultiphaseWellKernels::RateInitializationKernel::Launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                                                    m_oilPhaseIndex,
                                                                                                    wellControls,
                                                                                                    wellElemPhaseDens,
                                                                                                    wellElemTotalDens,
                                                                                                    connRate );


  } );
}

void CompositionalMultiphaseWell::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                     real64 const dt,
                                                     DomainPartition const & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // saved current dt for residual normalization
  m_currentDt = dt;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    WellControls const & wellControls = GetWellControls( subRegion );

    // get a reference to the degree-of-freedom numbers
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString );

    // get the info stored on well elements
    arrayView2d< real64 const > const & wellElemCompFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
    arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    FluxKernel::Launch< parallelDevicePolicy<> >( subRegion.size(),
                                                  dofManager.rankOffset(),
                                                  NumFluidComponents(),
                                                  NumDofPerResElement(),
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

void CompositionalMultiphaseWell::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const GEOSX_UNUSED_PARAM( dt ),
                                                              DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degrees of freedom and ghosting info
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get the properties on the well element
    arrayView2d< real64 const > const & wellElemPhaseVolFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );
    arrayView2d< real64 const > const & dWellElemPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
    arrayView3d< real64 const > const & dWellElemPhaseVolFrac_dComp =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    arrayView1d< real64 const > const & wellElemVolume =
      subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

    VolumeBalanceKernel::Launch< parallelDevicePolicy<> >( subRegion.size(),
                                                           NumFluidComponents(),
                                                           NumFluidPhases(),
                                                           NumDofPerWellElement(),
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
CompositionalMultiphaseWell::CalculateResidualNorm( DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  real64 localResidualNorm = 0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    WellControls const & wellControls = GetWellControls( subRegion );

    ResidualNormKernel::Launch< parallelDevicePolicy<>,
                                parallelDeviceReduce >( localRhs,
                                                        dofManager.rankOffset(),
                                                        subRegion.IsLocallyOwned(),
                                                        subRegion.GetTopWellElementIndex(),
                                                        NumDofPerWellElement(),
                                                        m_oilPhaseIndex,
                                                        wellControls,
                                                        wellElemDofNumber,
                                                        wellElemGhostRank,
                                                        phaseDens,
                                                        totalDens,
                                                        m_currentDt,
                                                        &localResidualNorm );
  } );
  return sqrt( MpiWrapper::Sum( localResidualNorm, MPI_COMM_GEOSX ) );
}

real64
CompositionalMultiphaseWell::ScalingForSystemSolution( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  // check if we want to rescale the Newton update
  if( m_maxCompFracChange >= 1.0 )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  real64 scalingFactor = 1.0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers on well elements and ghosting info
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    arrayView2d< real64 const > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dWellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

    real64 const subRegionScalingFactor =
      SolutionScalingKernel::Launch< parallelDevicePolicy<>,
                                     parallelDeviceReduce >( localSolution,
                                                             dofManager.rankOffset(),
                                                             NumFluidComponents(),
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

  return LvArray::math::max( MpiWrapper::Min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
}

bool
CompositionalMultiphaseWell::CheckSystemSolution( DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  int localCheck = 1;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers on well elements and ghosting info
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    arrayView2d< real64 const > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dWellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

    localIndex const subRegionSolutionCheck =
      SolutionCheckKernel::Launch< parallelDevicePolicy<>,
                                   parallelDeviceReduce >( localSolution,
                                                           dofManager.rankOffset(),
                                                           NumFluidComponents(),
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
  return MpiWrapper::Min( localCheck );
}

void CompositionalMultiphaseWell::ComputePerforationRates( WellElementSubRegion & subRegion,
                                                           localIndex const GEOSX_UNUSED_PARAM( targetIndex ) )
{
  GEOSX_MARK_FUNCTION;

  PerforationData * const perforationData = subRegion.GetPerforationData();

  // get depth
  arrayView1d< real64 const > const & wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // get well primary variables on well elements
  arrayView1d< real64 const > const & wellElemPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dWellElemPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView2d< real64 const > const & wellElemCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
  arrayView2d< real64 const > const & dWellElemCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  arrayView1d< real64 const > const & wellElemTotalMassDens =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString );
  arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString );
  arrayView2d< real64 const > const & dWellElemTotalMassDens_dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString );

  arrayView2d< real64 const > const & wellElemCompFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
  arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // get well variables on perforations
  arrayView1d< real64 const > const & perfGravCoef =
    perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
  arrayView1d< localIndex const > const & perfWellElemIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );
  arrayView1d< real64 const > const & perfTrans =
    perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::wellTransmissibilityString );

  arrayView2d< real64 > const & compPerfRate =
    perforationData->getReference< array2d< real64 > >( viewKeyStruct::compPerforationRateString );
  arrayView3d< real64 > const & dCompPerfRate_dPres =
    perforationData->getReference< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString );
  arrayView4d< real64 > const & dCompPerfRate_dComp =
    perforationData->getReference< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString );

  // get the element region, subregion, index
  arrayView1d< localIndex const > const & resElementRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
  arrayView1d< localIndex const > const & resElementSubRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
  arrayView1d< localIndex const > const & resElementIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

  PerforationKernel::Launch< parallelDevicePolicy<> >( perforationData->size(),
                                                       NumFluidComponents(),
                                                       NumFluidPhases(),
                                                       m_resPres.toNestedViewConst(),
                                                       m_deltaResPres.toNestedViewConst(),
                                                       m_resPhaseMob.toNestedViewConst(),
                                                       m_dResPhaseMob_dPres.toNestedViewConst(),
                                                       m_dResPhaseMob_dCompDens.toNestedViewConst(),
                                                       m_dResPhaseVolFrac_dPres.toNestedViewConst(),
                                                       m_dResPhaseVolFrac_dCompDens.toNestedViewConst(),
                                                       m_dResCompFrac_dCompDens.toNestedViewConst(),
                                                       m_resPhaseVisc.toNestedViewConst(),
                                                       m_dResPhaseVisc_dPres.toNestedViewConst(),
                                                       m_dResPhaseVisc_dComp.toNestedViewConst(),
                                                       m_resPhaseCompFrac.toNestedViewConst(),
                                                       m_dResPhaseCompFrac_dPres.toNestedViewConst(),
                                                       m_dResPhaseCompFrac_dComp.toNestedViewConst(),
                                                       m_resPhaseRelPerm.toNestedViewConst(),
                                                       m_dResPhaseRelPerm_dPhaseVolFrac.toNestedViewConst(),
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
CompositionalMultiphaseWell::ApplySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  // update all the fields using the global damping coefficients
  dofManager.addVectorToField( localSolution,
                               WellElementDofName(),
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( localSolution,
                               WellElementDofName(),
                               viewKeyStruct::deltaGlobalCompDensityString,
                               scalingFactor,
                               1, m_numDofPerWellElement - 1 );

  dofManager.addVectorToField( localSolution,
                               WellElementDofName(),
                               viewKeyStruct::deltaMixtureConnRateString,
                               scalingFactor,
                               m_numDofPerWellElement - 1, m_numDofPerWellElement );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    ChopNegativeDensities( domain );
  }

  // synchronize
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaGlobalCompDensityString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaMixtureConnRateString ) );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );

  // update properties
  UpdateStateAll( domain );

}

void CompositionalMultiphaseWell::ChopNegativeDensities( DomainPartition & domain )
{
  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 const > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 > const & dWellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
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


void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {

    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView2d< real64 > const & dWellElemGlobalCompDensity =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );
    arrayView1d< real64 > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      // extract solution and apply to dP
      dWellElemPressure[iwelem] = 0;
      dConnRate[iwelem] = 0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dWellElemGlobalCompDensity[iwelem][ic] = 0;
      }
    } );
  } );

  // call constitutive models
  UpdateStateAll( domain );
}


void CompositionalMultiphaseWell::ResetViews( DomainPartition & domain )
{
  WellSolverBase::ResetViews( domain );

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager & elemManager = *mesh.getElemManager();

  CompositionalMultiphaseFlow & flowSolver = *getParent()->GetGroup< CompositionalMultiphaseFlow >( GetFlowSolverName() );

  {
    using keys = CompositionalMultiphaseFlow::viewKeyStruct;

    m_resPres.clear();
    m_resPres = elemManager.ConstructArrayViewAccessor< real64, 1 >( keys::pressureString );
    m_resPres.setName( getName() + "/accessors/" + keys::pressureString );

    m_deltaResPres.clear();
    m_deltaResPres = elemManager.ConstructArrayViewAccessor< real64, 1 >( keys::deltaPressureString );
    m_deltaResPres.setName( getName() + "/accessors/" + keys::deltaPressureString );

    m_resCompDens.clear();
    m_resCompDens = elemManager.ConstructArrayViewAccessor< real64, 2 >( keys::globalCompDensityString );
    m_resCompDens.setName( getName() + "/accessors/" + keys::globalCompDensityString );

    m_dResCompFrac_dCompDens.clear();
    m_dResCompFrac_dCompDens = elemManager.ConstructArrayViewAccessor< real64, 3 >( keys::dGlobalCompFraction_dGlobalCompDensityString );
    m_dResCompFrac_dCompDens.setName( getName() + "/accessors/" + keys::dGlobalCompFraction_dGlobalCompDensityString );

    m_resPhaseMob.clear();
    m_resPhaseMob = elemManager.ConstructArrayViewAccessor< real64, 2 >( keys::phaseMobilityString );
    m_resPhaseMob.setName( getName() + "/accessors/" + keys::phaseMobilityString );

    m_dResPhaseMob_dPres.clear();
    m_dResPhaseMob_dPres = elemManager.ConstructArrayViewAccessor< real64, 2 >( keys::dPhaseMobility_dPressureString );
    m_dResPhaseMob_dPres.setName( getName() + "/accessors/" + keys::dPhaseMobility_dPressureString );

    m_dResPhaseMob_dCompDens.clear();
    m_dResPhaseMob_dCompDens = elemManager.ConstructArrayViewAccessor< real64, 3 >( keys::dPhaseMobility_dGlobalCompDensityString );
    m_dResPhaseMob_dCompDens.setName( getName() + "/accessors/" + keys::dPhaseMobility_dGlobalCompDensityString );

    m_resPhaseVolFrac.clear();
    m_resPhaseVolFrac = elemManager.ConstructArrayViewAccessor< real64, 2 >( keys::phaseVolumeFractionString );
    m_resPhaseVolFrac.setName( getName() + "/accessors/" + keys::phaseVolumeFractionString );

    m_dResPhaseVolFrac_dPres.clear();
    m_dResPhaseVolFrac_dPres = elemManager.ConstructArrayViewAccessor< real64, 2 >( keys::dPhaseVolumeFraction_dPressureString );
    m_dResPhaseVolFrac_dPres.setName( getName() + "/accessors/" + keys::dPhaseVolumeFraction_dPressureString );

    m_dResPhaseVolFrac_dCompDens.clear();
    m_dResPhaseVolFrac_dCompDens = elemManager.ConstructArrayViewAccessor< real64, 3 >( keys::dPhaseVolumeFraction_dGlobalCompDensityString );
    m_dResPhaseVolFrac_dCompDens.setName( getName() + "/accessors/" + keys::dPhaseVolumeFraction_dGlobalCompDensityString );

  }
  {
    using keys = MultiFluidBase::viewKeyStruct;

    m_resPhaseMassDens.clear();
    m_resPhaseMassDens = elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >( keys::phaseMassDensityString,
                                                                                      flowSolver.targetRegionNames(),
                                                                                      flowSolver.fluidModelNames() );
    m_resPhaseMassDens.setName( getName() + "/accessors/" + keys::phaseMassDensityString );

    m_resPhaseVisc.clear();
    m_resPhaseVisc = elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >( keys::phaseViscosityString,
                                                                                  flowSolver.targetRegionNames(),
                                                                                  flowSolver.fluidModelNames() );
    m_resPhaseVisc.setName( getName() + "/accessors/" + keys::phaseViscosityString );

    m_dResPhaseVisc_dPres.clear();
    m_dResPhaseVisc_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >( keys::dPhaseViscosity_dPressureString,
                                                                                         flowSolver.targetRegionNames(),
                                                                                         flowSolver.fluidModelNames() );
    m_dResPhaseVisc_dPres.setName( getName() + "/accessors/" + keys::dPhaseViscosity_dPressureString );

    m_dResPhaseVisc_dComp.clear();
    m_dResPhaseVisc_dComp = elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseViscosity_dGlobalCompFractionString,
                                                                                         flowSolver.targetRegionNames(),
                                                                                         flowSolver.fluidModelNames() );
    m_dResPhaseVisc_dComp.setName( getName() + "/accessors/" + keys::dPhaseViscosity_dGlobalCompFractionString );

    m_resPhaseCompFrac.clear();
    m_resPhaseCompFrac = elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >( keys::phaseCompFractionString,
                                                                                      flowSolver.targetRegionNames(),
                                                                                      flowSolver.fluidModelNames() );
    m_resPhaseCompFrac.setName( getName() + "/accessors/" + keys::phaseCompFractionString );

    m_dResPhaseCompFrac_dPres.clear();
    m_dResPhaseCompFrac_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseCompFraction_dPressureString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseCompFrac_dPres.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dPressureString );

    m_dResPhaseCompFrac_dComp.clear();
    m_dResPhaseCompFrac_dComp = elemManager.ConstructMaterialArrayViewAccessor< real64, 5 >( keys::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseCompFrac_dComp.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dGlobalCompFractionString );

  }
  {
    using keys = RelativePermeabilityBase::viewKeyStruct;

    m_resPhaseRelPerm.clear();
    m_resPhaseRelPerm = elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >( keys::phaseRelPermString,
                                                                                     flowSolver.targetRegionNames(),
                                                                                     flowSolver.relPermModelNames() );
    m_resPhaseRelPerm.setName( getName() + "/accessors/" + keys::phaseRelPermString );

    m_dResPhaseRelPerm_dPhaseVolFrac.clear();
    m_dResPhaseRelPerm_dPhaseVolFrac = elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                                    flowSolver.targetRegionNames(),
                                                                                                    flowSolver.relPermModelNames() );
    m_dResPhaseRelPerm_dPhaseVolFrac.setName( getName() + "/accessors/" + keys::dPhaseRelPerm_dPhaseVolFractionString );

  }
}


void CompositionalMultiphaseWell::FormPressureRelations( DomainPartition const & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {

    WellControls & wellControls = GetWellControls( subRegion );

    // get the degrees of freedom, depth info, next welem index
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get primary variables on well elements
    arrayView1d< real64 const > const & wellElemPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    // get total mass density on well elements (for potential calculations)
    arrayView1d< real64 const > const & wellElemTotalMassDens =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString );
    arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString );
    arrayView2d< real64 const > const & dWellElemTotalMassDens_dCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString );


    localIndex const controlHasSwitched =
      PressureRelationKernel::Launch< parallelDevicePolicy<>,
                                      parallelDeviceReduce >( subRegion.size(),
                                                              dofManager.rankOffset(),
                                                              subRegion.IsLocallyOwned(),
                                                              subRegion.GetTopWellElementIndex(),
                                                              NumFluidComponents(),
                                                              m_oilPhaseIndex,
                                                              NumDofPerResElement(),
                                                              wellControls,
                                                              wellElemDofNumber,
                                                              wellElemGravCoef,
                                                              nextWellElemIndex,
                                                              wellElemPres,
                                                              dWellElemPres,
                                                              wellElemTotalMassDens,
                                                              dWellElemTotalMassDens_dPres,
                                                              dWellElemTotalMassDens_dCompDens,
                                                              localMatrix,
                                                              localRhs );

    if( controlHasSwitched == 1 )
    {
      // TODO: move the switch logic into wellControls

      // TODO: implement a more general switch when more then two constraints per well type are allowed

      if( wellControls.GetControl() == WellControls::Control::BHP )
      {
        WellControls::Type const wellType = wellControls.GetType();
        if( wellType == WellControls::Type::PRODUCER )
        {
          wellControls.SetControl( WellControls::Control::OILVOLRATE,
                                   wellControls.GetTargetRate() );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                << " from BHP constraint to oil volumetric rate constraint" );
        }
        else
        {
          wellControls.SetControl( WellControls::Control::TOTALVOLRATE,
                                   wellControls.GetTargetRate() );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                << " from BHP constraint to total volumetric rate constraint" );
        }
      }
      else
      {
        wellControls.SetControl( WellControls::Control::BHP,
                                 wellControls.GetTargetBHP() );

        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from rate constraint to BHP constraint" );
      }
    }
  } );
}


void CompositionalMultiphaseWell::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition & domain )
{
  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {

    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    arrayView2d< real64 > const & wellElemGlobalCompDensity =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dWellElemGlobalCompDensity =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d< real64 > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        wellElemGlobalCompDensity[iwelem][ic] += dWellElemGlobalCompDensity[iwelem][ic];
      }
      connRate[iwelem] += dConnRate[iwelem];
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseWell, string const &, Group * const )
}// namespace geosx
