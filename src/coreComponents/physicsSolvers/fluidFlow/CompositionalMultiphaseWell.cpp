/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "managers/DomainPartition.hpp"
#include "wells/PerforationData.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "wells/WellControls.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseFlowKernels;

CompositionalMultiphaseWell::CompositionalMultiphaseWell( const string & name,
                                                          Group * const parent )
  :
  WellSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_temperature( 0.0 ),
  m_useMass( false )
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
}

void CompositionalMultiphaseWell::PostProcessInput()
{
  WellSolverBase::PostProcessInput();
  CheckModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString );

  CompositionalMultiphaseFlow const * const flowSolver = getParent()->GetGroup< CompositionalMultiphaseFlow >( GetFlowSolverName() );
  GEOSX_ERROR_IF( flowSolver == nullptr,
                  "Flow solver " << GetFlowSolverName() << " not found or incompatible type "
                                                           "(referenced from well solver " << getName() << ")" );
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

    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::mixtureDensityString )->setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::dMixtureDensity_dPressureString )->
      setRestartFlags( RestartFlags::NO_WRITE );
    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString )->
      setRestartFlags( RestartFlags::NO_WRITE );

    PerforationData * const perforationData = subRegion.GetPerforationData();
    perforationData->registerWrapper< array2d< real64 > >( viewKeyStruct::compPerforationRateString );
    perforationData->registerWrapper< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString )->
      setRestartFlags( RestartFlags::NO_WRITE );
    perforationData->registerWrapper< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString )->
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

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    ResizeFields( subRegion );
  } );
}

void CompositionalMultiphaseWell::ResizeFields( WellElementSubRegion & subRegion )
{
  PerforationData * const perforationData = subRegion.GetPerforationData();

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString ).resizeDimension< 1 >( NC );
  subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString ).resizeDimension< 1 >( NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString ).resizeDimension< 1 >( NC );
  subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString ).resizeDimension< 1, 2 >( NC, NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString ).resizeDimension< 1 >( NP );
  subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString ).resizeDimension< 1 >( NP );
  subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString ).resizeDimension< 1, 2 >( NP, NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString ).resizeDimension< 1 >( NC );

  perforationData->getReference< array2d< real64 > >( viewKeyStruct::compPerforationRateString ).resizeDimension< 1 >( NC );
  perforationData->getReference< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString ).resizeDimension< 1, 2 >( 2, NC );
  perforationData->getReference< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString ).resizeDimension< 1, 2, 3 >( 2, NC, NC );

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


void CompositionalMultiphaseWell::UpdateMixtureDensity( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  // get well secondary variables on well elements
  arrayView1d< real64 > const & wellElemMixtureDensity =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureDensityString );
  arrayView1d< real64 > const & dWellElemMixtureDensity_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dMixtureDensity_dPressureString );
  arrayView2d< real64 > const & dWellElemMixtureDensity_dComp =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

  arrayView2d< real64 const > const & wellElemPhaseVolFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );
  arrayView2d< real64 const > const & dWellElemPhaseVolFrac_dPres =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  arrayView3d< real64 const > const & dWellElemPhaseVolFrac_dComp =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // get constitutive data
  MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & wellElemPhaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dWellElemPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dWellElemPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  stackArray1d< real64, maxNumComp > dDens_dC( NC );

  for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
  {

    // reset to zero
    wellElemMixtureDensity[iwelem]        = 0.0;
    dWellElemMixtureDensity_dPres[iwelem] = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dWellElemMixtureDensity_dComp[iwelem][ic] = 0.0;
    }

    // increment mixture velocity
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      wellElemMixtureDensity[iwelem] += wellElemPhaseVolFrac[iwelem][ip] * wellElemPhaseDens[iwelem][0][ip];
      dWellElemMixtureDensity_dPres[iwelem] += dWellElemPhaseVolFrac_dPres[iwelem][ip] * wellElemPhaseDens[iwelem][0][ip]
                                               + wellElemPhaseVolFrac[iwelem][ip] * dWellElemPhaseDens_dPres[iwelem][0][ip];

      dDens_dC = 0.0;

      applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseDens_dComp[iwelem][0][ip], dDens_dC );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dWellElemMixtureDensity_dComp[iwelem][ic] += dWellElemPhaseVolFrac_dComp[iwelem][ip][ic] * wellElemPhaseDens[iwelem][0][ip]
                                                     + wellElemPhaseVolFrac[iwelem][ip] * dDens_dC[ic];
      }
    }
  }
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

  KernelLaunchSelector1< ComponentFractionKernel >( m_numComponents,
                                                    subRegion.size(),
                                                    compDens,
                                                    dCompDens,
                                                    compFrac,
                                                    dCompFrac_dCompDens );
}

void CompositionalMultiphaseWell::UpdateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView1d< real64 const > const & pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
  arrayView2d< real64 const > const & compFrac = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );

  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forAll< serialPolicy >( subRegion.size(), [&] ( localIndex const a )
  {
    fluid.PointUpdate( pres[a] + dPres[a], m_temperature, compFrac[a], a, 0 );
  } );
  //fluid->BatchUpdate( pres, temp, compFrac );
}

void CompositionalMultiphaseWell::UpdatePhaseVolumeFraction( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseVolFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );

  arrayView2d< real64 > const & dPhaseVolFrac_dPres =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 > const & dPhaseVolFrac_dComp =
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

  KernelLaunchSelector2< PhaseVolumeFractionKernel >( m_numComponents, m_numPhases,
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
                                                      dPhaseVolFrac_dComp );
}

void CompositionalMultiphaseWell::UpdateState( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  // update properties
  UpdateComponentFraction( subRegion );
  UpdateFluidModel( subRegion, targetIndex );
  UpdatePhaseVolumeFraction( subRegion, targetIndex );
  UpdateMixtureDensity( subRegion, targetIndex );

  // update perforation rates
  ComputePerforationRates( subRegion, targetIndex );
}

void CompositionalMultiphaseWell::InitializeWells( DomainPartition * const domain )
{

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & resPressure = m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resCompDens = m_resGlobalCompDensity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resPhaseVolFrac = m_resPhaseVolFrac;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & resPhaseDens    = m_resPhaseDens;

  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {

    WellControls const & wellControls = GetWellControls( subRegion );
    PerforationData const * const perforationData = subRegion.GetPerforationData();

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
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // get well secondary variables on well elements
    MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    // 1) Loop over all perforations to compute an average mixture density
    //    and component fraction
    real64 avgTotalDensity = 0.0;
    real64 avgMixtureDensity = 0.0;
    stackArray1d< real64, maxNumComp > avgCompFrac( NC );

    avgCompFrac = 0.0;

    // define a reservoir pressure used for initialization
    real64 resPres = ( wellControls.GetType() == WellControls::Type::PRODUCER ) ? 1e20 : 0;

    for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      // save min pressure for producer
      if( wellControls.GetType() == WellControls::Type::PRODUCER && resPres > resPressure[er][esr][ei] )
      {
        resPres = resPressure[er][esr][ei];
      }
      // save max pressure for injector
      else if( wellControls.GetType() == WellControls::Type::INJECTOR && resPres < resPressure[er][esr][ei] )
      {
        resPres = resPressure[er][esr][ei];
      }

      // increment the average mixture density
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        real64 const resDensity = resPhaseDens[er][esr][ei][0][ip];
        real64 const resVolFrac = resPhaseVolFrac[er][esr][ei][ip];
        avgMixtureDensity += resVolFrac * resDensity;
      }

      // increment the average total density
      real64 perfTotalDensity = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        perfTotalDensity += resCompDens[er][esr][ei][ic];
      }
      avgTotalDensity += perfTotalDensity;

      // increment the average component fraction
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        avgCompFrac[ic] += resCompDens[er][esr][ei][ic] / perfTotalDensity;

      }
    }

    // communicate the pressures to the ranks without perforations
    // this will be used to initialize the pressure, starting by the owner rank
    if( wellControls.GetType() == WellControls::Type::PRODUCER )
    {
      resPres = MpiWrapper::Min( resPres );
    }
    else if( wellControls.GetType() == WellControls::Type::INJECTOR )
    {
      resPres = MpiWrapper::Max( resPres );
    }

    // compute average densities
    globalIndex const numPerforationsGlobal = perforationData->GetNumPerforationsGlobal();

    avgMixtureDensity = MpiWrapper::Sum( avgMixtureDensity );
    avgTotalDensity = MpiWrapper::Sum( avgTotalDensity );

    avgMixtureDensity /= numPerforationsGlobal;
    avgTotalDensity /= numPerforationsGlobal;

    // compute average component fraction
    if( wellControls.GetType() == WellControls::Type::PRODUCER )
    {
      // use average comp frac from reservoir
      real64 compFracSum = 0;
      real64 const tol = 1e-7;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        avgCompFrac[ic] = MpiWrapper::Sum( avgCompFrac[ic] );
        avgCompFrac[ic] /= numPerforationsGlobal;
        compFracSum += avgCompFrac[ic];
      }
      GEOSX_ERROR_IF( compFracSum < 1 - tol || compFracSum > 1 + tol,
                      "Invalid well initialization: negative pressure was found" );
    }
    else // injector
    {
      arrayView1d< real64 const > const & injection = wellControls.GetInjectionStream();

      // use average comp frac from XML file
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        avgCompFrac[ic] = injection[ic];
      }
    }

    // set the global component fractions to avgCompFrac
    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        wellElemCompFrac[iwelem][ic] = avgCompFrac[ic];
      }
    }

    real64 pressureControl = 0.0;
    real64 gravCoefControl = 0.0;

    if( subRegion.IsLocallyOwned() )
    {

      localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();
      gravCoefControl = wellElemGravCoef[iwelemControl];

      // 2) Initialize the reference pressure
      real64 const & targetBHP = wellControls.GetTargetBHP();
      if( wellControls.GetControl() == WellControls::Control::BHP )
      {
        // if pressure constraint, set the ref pressure at the constraint
        pressureControl = targetBHP;
      }
      else // rate control
      {
        // if rate constraint, set the ref pressure slightly
        // above/below the target pressure depending on well type
        pressureControl = ( wellControls.GetType() == WellControls::Type::PRODUCER )
                          ? 0.5 * resPres // hard-coded values come from personal communication with Hui
                          : 2.0 * resPres;
      }

      wellElemPressure[iwelemControl] = pressureControl;
    }

    // TODO optimize
    MpiWrapper::Broadcast( pressureControl, subRegion.GetTopRank() );
    MpiWrapper::Broadcast( gravCoefControl, subRegion.GetTopRank() );

    GEOSX_ERROR_IF( pressureControl <= 0, "Invalid well initialization: negative pressure was found" );

    // 3) Estimate the pressures in the well elements using this avgDensity
    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {
      wellElemPressure[iwelem] = pressureControl
                                 + avgMixtureDensity * ( wellElemGravCoef[iwelem] - gravCoefControl );

    } );


    // 4) Back calculate component densities
    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      fluid.PointUpdate( wellElemPressure[iwelem], m_temperature,
                         wellElemCompFrac[iwelem], iwelem, 0 );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        wellElemCompDens[iwelem][ic] = avgCompFrac[ic] * totalDens[iwelem][0];
      }
    }

    // 5) Recompute all the pressure-dependent properties
    UpdateState( subRegion, targetIndex );


    // 6) Estimate the connection rates based on the min/max pressure
    real64 const targetRate = wellControls.GetTargetRate();
    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      if( wellControls.GetControl() == WellControls::Control::BHP )
      {
        // if BHP constraint set rate below the absolute max rate
        // with the appropriate sign (negative for prod, positive for inj)
        connRate[iwelem] = ( wellControls.GetType() == WellControls::Type::PRODUCER )
                           ? std::max( 0.1 * targetRate, -1e3 )  // hard-coded values come from personal communication
                                                                 // with Hui
                           : std::min( 0.1 * targetRate, 1e3 );
      }
      else
      {
        connRate[iwelem] = targetRate;
      }
    }
  } );
}

void CompositionalMultiphaseWell::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                     real64 const dt,
                                                     DomainPartition const * const domain,
                                                     DofManager const * const dofManager,
                                                     ParallelMatrix * const matrix,
                                                     ParallelVector * const rhs )
{
  MeshLevel const & meshLevel = *domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC      = m_numComponents;
  localIndex const resNDOF = NumDofPerResElement();

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    WellControls const & wellControls = GetWellControls( subRegion );

    // get a reference to the degree-of-freedom numbers
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

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

    WellControls::Type const wellType = wellControls.GetType();
    arrayView1d< real64 const > const & injection = wellControls.GetInjectionStream();

    // loop over the well elements to compute the fluxes between elements
    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {

      // rank that owns the current well element assembles the flux between current and next
      if( wellElemGhostRank[iwelem] >= 0 )
      {
        return;
      }

      // create local work arrays
      stackArray1d< real64, maxNumComp > compFracUp( NC );
      stackArray1d< real64, maxNumComp > dCompFrac_dPresUp( NC );
      stackArray2d< real64, maxNumComp * maxNumComp > dCompFrac_dCompDensUp( NC, NC );

      stackArray1d< real64, maxNumComp > compFlux( NC );
      stackArray1d< real64, maxNumComp > dCompFlux_dRate( NC );
      stackArray1d< real64, maxNumComp > dCompFlux_dPresUp( NC );
      stackArray2d< real64, maxNumComp * maxNumComp > dCompFlux_dCompDensUp( NC, NC );

      stackArray1d< real64, 2 * maxNumComp > localFlux( 2 * NC );
      stackArray2d< real64, 2 * maxNumComp > localFluxJacobian_dRate( 2 * NC, 1 );
      stackArray2d< real64, 2 * maxNumComp * maxNumDof > localFluxJacobian_dPresCompUp( 2 * NC, resNDOF );

      stackArray1d< globalIndex, 2 * maxNumComp > eqnRowIndices( 2 * NC );
      stackArray1d< globalIndex, maxNumDof > dofColIndices_dPresCompUp( resNDOF );
      stackArray1d< globalIndex, 1 > dofColIndices_dRate( 1 );

      // Step 1) decide the upwind well element

      /*  currentConnRate < 0 flow from iwelem to iwelemNext
       *  currentConnRate > 0 flow from iwelemNext to iwelem
       *  With this convention, currentConnRate < 0 at the last connection for a producer
       *                        currentConnRate > 0 at the last connection for a injector
       */

      localIndex const iwelemNext = nextWellElemIndex[iwelem];
      real64 const currentConnRate = connRate[iwelem] + dConnRate[iwelem];
      localIndex iwelemUp = -1;

      if( iwelemNext < 0 && wellType == WellControls::Type::INJECTOR ) // exit connection, injector
      {
        // we still need to define iwelemUp for Jacobian assembly
        iwelemUp = iwelem;

        // just copy the injection stream into compFrac
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          compFracUp[ic] = injection[ic];
        }
        // zero out the derivatives wrt composition
        dCompFrac_dCompDensUp = 0.0;

      }
      else
      {
        // first set iwelemUp to the upstream cell
        if( ( iwelemNext < 0 && wellType == WellControls::Type::PRODUCER )  // exit connection, producer
            || currentConnRate < 0 ) // not an exit connection, iwelem is upstream
        {
          iwelemUp = iwelem;
        }
        else // not an exit connection, iwelemNext is upstream
        {
          iwelemUp = iwelemNext;
        }

        // copy the vars of iwelemUp into compFrac
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          compFracUp[ic] = wellElemCompFrac[iwelemUp][ic];
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCompFrac_dCompDensUp[ic][jc] = dWellElemCompFrac_dCompDens[iwelemUp][ic][jc];
          }
        }
      }

      // Step 2) compute upstream transport coefficient

      compFlux = 0.0;
      dCompFlux_dRate = 0.0;
      dCompFlux_dPresUp = 0.0;
      dCompFlux_dCompDensUp = 0.0;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compFlux[ic] = compFracUp[ic] * currentConnRate;
        dCompFlux_dRate[ic] = compFracUp[ic];
        dCompFlux_dPresUp[ic] = 0.0; // none of these quantities depend on pressure
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFlux_dCompDensUp[ic][jc] = dCompFrac_dCompDensUp[ic][jc] * currentConnRate;
        }
      }

      globalIndex const offsetUp = wellElemDofNumber[iwelemUp];
      globalIndex const offsetCurrent = wellElemDofNumber[iwelem];

      if( iwelemNext < 0 )  // exit connection
      {

        // for this case, we only need NC mass conservation equations
        // so we do not use the arrays initialized before the loop
        stackArray1d< real64, maxNumComp > oneSidedFlux( NC );
        stackArray2d< real64, maxNumComp > oneSidedFluxJacobian_dRate( NC, 1 );
        stackArray2d< real64, maxNumComp * maxNumDof > oneSidedFluxJacobian_dPresCompUp( NC, resNDOF );

        stackArray1d< globalIndex, maxNumComp > oneSidedEqnRowIndices( NC );
        stackArray1d< globalIndex, maxNumDof > oneSidedDofColIndices_dPresCompUp( resNDOF );
        stackArray1d< globalIndex, 1 > oneSidedDofColIndices_dRate( 1 );

        // flux terms
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          oneSidedFlux[ic] = -dt * compFlux[ic];

          // derivative with respect to rate
          oneSidedFluxJacobian_dRate( ic, 0 ) = -dt * dCompFlux_dRate[ic];

          // derivative with respect to upstream pressure
          oneSidedFluxJacobian_dPresCompUp[ic][0] = -dt * dCompFlux_dPresUp[ic];

          // derivatives with respect to upstream component densities
          for( localIndex jdof = 0; jdof < NC; ++jdof )
          {
            oneSidedFluxJacobian_dPresCompUp[ic][jdof + 1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
          }

        }

        // jacobian indices
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          // mass balance equations for all components
          oneSidedEqnRowIndices[ic] = offsetUp + RowOffset::MASSBAL + ic;
        }

        // in the dof ordering used in this class, there are 1 pressure dofs
        // and NC compDens dofs before the rate dof in this block
        localIndex const dRateColOffset = ColOffset::DCOMP + NC;
        oneSidedDofColIndices_dRate[0] = offsetCurrent + dRateColOffset;

        for( localIndex jdof = 0; jdof < resNDOF; ++jdof )
        {
          // dofs are the **upstream** pressure and component densities
          oneSidedDofColIndices_dPresCompUp[jdof] = offsetUp + ColOffset::DPRES + jdof;
        }

        rhs->add( oneSidedEqnRowIndices,
                  oneSidedFlux );

        matrix->add( oneSidedEqnRowIndices,
                     oneSidedDofColIndices_dRate,
                     oneSidedFluxJacobian_dRate );

        matrix->add( oneSidedEqnRowIndices,
                     oneSidedDofColIndices_dPresCompUp,
                     oneSidedFluxJacobian_dPresCompUp );

      }
      else // not an exit connection
      {
        globalIndex const offsetNext = wellElemDofNumber[iwelemNext];

        // flux terms
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          localFlux[ElemTag::NEXT * NC + ic] = dt * compFlux[ic];
          localFlux[ElemTag::CURRENT * NC + ic] = -dt * compFlux[ic];

          // derivative with respect to rate
          localFluxJacobian_dRate( ElemTag::NEXT * NC + ic, 0 ) = dt * dCompFlux_dRate[ic];
          localFluxJacobian_dRate( ElemTag::CURRENT * NC + ic, 0 ) = -dt * dCompFlux_dRate[ic];

          // derivative with respect to upstream pressure
          localFluxJacobian_dPresCompUp[ElemTag::NEXT * NC + ic][0] = dt * dCompFlux_dPresUp[ic];
          localFluxJacobian_dPresCompUp[ElemTag::CURRENT * NC + ic][0] = -dt * dCompFlux_dPresUp[ic];

          // derivatives with respect to upstream component densities
          for( localIndex jdof = 0; jdof < NC; ++jdof )
          {
            localFluxJacobian_dPresCompUp[ElemTag::NEXT * NC + ic][jdof + 1] =
              dt * dCompFlux_dCompDensUp[ic][jdof];
            localFluxJacobian_dPresCompUp[ElemTag::CURRENT * NC + ic][jdof + 1] =
              -dt * dCompFlux_dCompDensUp[ic][jdof];
          }

        }

        // jacobian indices
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          // mass balance equations for all components
          eqnRowIndices[ElemTag::NEXT * NC + ic] = offsetNext + RowOffset::MASSBAL + ic;
          eqnRowIndices[ElemTag::CURRENT * NC + ic] = offsetCurrent + RowOffset::MASSBAL + ic;
        }

        // in the dof ordering used in this class, there are 1 pressure dofs
        // and NC compDens dofs before the rate dof in this block
        localIndex const dRateColOffset = ColOffset::DCOMP + NC;
        dofColIndices_dRate[0] = offsetCurrent + dRateColOffset;

        for( localIndex jdof = 0; jdof < resNDOF; ++jdof )
        {
          // dofs are the **upstream** pressure and component densities
          dofColIndices_dPresCompUp[jdof] = offsetUp + ColOffset::DPRES + jdof;
        }

        rhs->add( eqnRowIndices,
                  localFlux );

        matrix->add( eqnRowIndices,
                     dofColIndices_dRate,
                     localFluxJacobian_dRate );

        matrix->add( eqnRowIndices,
                     dofColIndices_dPresCompUp,
                     localFluxJacobian_dPresCompUp );
      }
    } );
  } );
}

void CompositionalMultiphaseWell::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const GEOSX_UNUSED_PARAM( dt ),
                                                              DomainPartition const * const domain,
                                                              DofManager const * const dofManager,
                                                              ParallelMatrix * const matrix,
                                                              ParallelVector * const rhs )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC        = m_numComponents;
  localIndex const NP        = m_numPhases;
  localIndex const welemNDOF = m_numDofPerWellElement;

  MeshLevel const & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degrees of freedom and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get the properties on the well element
    arrayView2d< real64 const > const & wellElemPhaseVolFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );
    arrayView2d< real64 const > const & dWellElemPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
    arrayView3d< real64 const > const & dWellElemPhaseVolFrac_dComp =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    arrayView1d< real64 const > const & wellElemVolume =
      subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {

      if( wellElemGhostRank[iwelem] >= 0 )
      {
        return;
      }

      stackArray1d< globalIndex, maxNumDof > localVolBalanceDOF( welemNDOF );
      stackArray1d< real64, maxNumDof > localVolBalanceJacobian( welemNDOF );

      localVolBalanceDOF = 0;
      localVolBalanceJacobian = 0;

      // get equation/dof indices
      globalIndex const offset = wellElemDofNumber[iwelem];
      localIndex const volBalRowOffset = RowOffset::MASSBAL + NC;
      globalIndex const localVolBalanceEqnIndex = offset + volBalRowOffset;
      for( localIndex jdof = 0; jdof < welemNDOF; ++jdof )
      {
        localVolBalanceDOF[jdof] = offset + ColOffset::DPRES + jdof;
      }

      real64 localVolBalance = 1.0;
      localVolBalanceJacobian = 0.0;

      // sum contributions to component accumulation from each phase
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        localVolBalance -= wellElemPhaseVolFrac[iwelem][ip];
        localVolBalanceJacobian[0] -= dWellElemPhaseVolFrac_dPres[iwelem][ip];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          localVolBalanceJacobian[jc + 1] -= dWellElemPhaseVolFrac_dComp[iwelem][ip][jc];
        }
      }

      // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
      for( localIndex idof = 0; idof < welemNDOF; ++idof )
      {
        localVolBalanceJacobian[idof] *= wellElemVolume[iwelem];
      }
      localVolBalance *= wellElemVolume[iwelem];

      rhs->add( localVolBalanceEqnIndex,
                localVolBalance );

      matrix->add( localVolBalanceEqnIndex,
                   localVolBalanceDOF,
                   localVolBalanceJacobian );
    } );
  } );
}


real64
CompositionalMultiphaseWell::CalculateResidualNorm( DomainPartition const * const domain,
                                                    DofManager const & dofManager,
                                                    ParallelVector const & rhs )
{
  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  MeshLevel const & meshLevel = *domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  string const wellDofKey = dofManager.getKey( WellElementDofName() );

  real64 residualNorm = 0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d< real64 const > const & wellElemVolume =
      subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

    MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        for( localIndex idof = 0; idof < NumDofPerWellElement(); ++idof )
        {
          real64 const normalizer = ( idof >= RowOffset::MASSBAL && idof < RowOffset::MASSBAL + NC )
                                    ? totalDens[iwelem][0] * wellElemVolume[iwelem]
                                    : 1;
          localIndex const lid = rhs.getLocalRowID( wellElemDofNumber[iwelem] + idof );
          real64 const val = localResidual[lid] / normalizer;
          residualNorm += val * val;
        }
      }
    }
  } );

  real64 const globalResidualNorm = MpiWrapper::Sum( residualNorm, MPI_COMM_GEOSX );
  return sqrt( globalResidualNorm );
}

bool
CompositionalMultiphaseWell::CheckSystemSolution( DomainPartition const * const domain,
                                                  DofManager const & dofManager,
                                                  ParallelVector const & solution,
                                                  real64 const scalingFactor )
{
  // get the update
  real64 const * localSolution = solution.extractLocalVector();

  MeshLevel const & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  string const wellDofKey = dofManager.getKey( WellElementDofName() );

  int isInvalidLocal = 0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {

    // get the degree of freedom numbers on well elements and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    arrayView2d< real64 const > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dWellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

    for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        // pressure
        localIndex lid = solution.getLocalRowID( wellElemDofNumber[iwelem] + ColOffset::DPRES );
        real64 const newPres = wellElemPressure[iwelem] + dWellElemPressure[iwelem]
                               + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          isInvalidLocal = 1;
        }

        // comp densities
        for( localIndex ic = 0; ic < m_numComponents; ++ic )
        {
          lid = solution.getLocalRowID( wellElemDofNumber[iwelem] + ic + 1 );
          real64 const newDens = wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic]
                                 + scalingFactor * localSolution[lid];
          if( newDens < 0.0 )
          {
            isInvalidLocal = 1;
          }
        }
      }
    }
  } );

  return !MpiWrapper::Max( isInvalidLocal, MPI_COMM_GEOSX );
}

void CompositionalMultiphaseWell::ComputePerforationRates( WellElementSubRegion & subRegion,
                                                           localIndex const GEOSX_UNUSED_PARAM( targetIndex ) )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;

  localIndex const NC = NumFluidComponents();
  localIndex const NP = NumFluidPhases();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & resPressure  = m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & dResPressure = m_deltaResPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resPhaseMob = m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dResPhaseMob_dPres = m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResPhaseMob_dComp = m_dResPhaseMob_dCompDens;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dResPhaseVolFrac_dPres = m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResPhaseVolFrac_dComp = m_dResPhaseVolFrac_dCompDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResCompFrac_dCompDens = m_dResCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & resPhaseVisc        = m_resPhaseVisc;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResPhaseVisc_dPres = m_dResPhaseVisc_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & dResPhaseVisc_dComp = m_dResPhaseVisc_dComp;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & resPhaseCompFrac    = m_resPhaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & dResPhaseCompFrac_dPres = m_dResPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 > > const & dResPhaseCompFrac_dComp = m_dResPhaseCompFrac_dComp;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & resPhaseRelPerm                = m_resPhaseRelPerm;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & dResPhaseRelPerm_dPhaseVolFrac = m_dResPhaseRelPerm_dPhaseVolFrac;

  PerforationData * const perforationData = subRegion.GetPerforationData();

  // get depth
  arrayView1d< real64 const > const &
  wellElemGravCoef = subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // get well primary variables on well elements
  arrayView1d< real64 const > const &
  wellElemPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const &
  dWellElemPressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  // get secondary well data on well elements
  arrayView1d< real64 const > const &
  wellElemMixtureDensity = subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureDensityString );
  arrayView1d< real64 const > const &
  dWellElemMixtureDensity_dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::dMixtureDensity_dPressureString );
  arrayView2d< real64 const > const &
  dWellElemMixtureDensity_dComp = subRegion.getReference< array2d< real64 > >( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

  arrayView2d< real64 const > const &
  wellElemCompFrac = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
  arrayView3d< real64 const > const &
  dWellElemCompFrac_dCompDens = subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // get well variables on perforations
  arrayView1d< real64 const > const &
  perfGravCoef = perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  arrayView1d< localIndex const > const &
  perfWellElemIndex = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d< real64 const > const &
  perfTransmissibility = perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::wellTransmissibilityString );

  arrayView2d< real64 > const &
  compPerfRate = perforationData->getReference< array2d< real64 > >( viewKeyStruct::compPerforationRateString );
  arrayView3d< real64 > const &
  dCompPerfRate_dPres = perforationData->getReference< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString );
  arrayView4d< real64 > const &
  dCompPerfRate_dComp = perforationData->getReference< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString );

  // get the element region, subregion, index
  arrayView1d< localIndex const > const &
  resElementRegion = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
  arrayView1d< localIndex const > const &
  resElementSubRegion = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
  arrayView1d< localIndex const > const &
  resElementIndex = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

  // local working variables and arrays
  stackArray1d< real64, maxNumComp > dPhaseCompFrac_dCompDens( NC );

  stackArray1d< real64, 2 >              pressure( 2 );
  stackArray1d< real64, 2 >              dPressure_dP( 2 );
  stackArray2d< real64, 2 * maxNumComp > dPressure_dC( 2, NC );

  stackArray1d< real64, 2 >              dFlux_dP( 2 );
  stackArray2d< real64, 2 * maxNumComp > dFlux_dC( 2, NC );

  stackArray1d< real64, 2 >              dMult_dP( 2 );
  stackArray2d< real64, 2 * maxNumComp > dMult_dC( 2, NC );

  stackArray1d< real64, maxNumComp > dMixtureDensity_dC( NC );
  stackArray1d< real64, maxNumComp > dResTotalMobility_dC( NC );

  stackArray2d< real64, 2 * maxNumComp >              phaseCompFrac( 2, NC );
  stackArray2d< real64, 2 * maxNumComp >              dPhaseCompFrac_dP( 2, NC );
  stackArray3d< real64, 2 * maxNumComp * maxNumComp > dPhaseCompFrac_dC( 2, NC, NC );

  stackArray1d< real64, maxNumComp > dVisc_dC( NC );
  stackArray1d< real64, maxNumComp > dRelPerm_dC( NC );

  stackArray1d< real64, 2 >              dPotDiff_dP( 2 );
  stackArray2d< real64, 2 * maxNumComp > dPotDiff_dC( 2, NC );

  stackArray1d< real64, 2 > multiplier( 2 );


  // loop over the perforations to compute the perforation rates
  for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
  {

    // reset the perforation rates
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      compPerfRate[iperf][ic] = 0.0;
      for( localIndex ke = 0; ke < 2; ++ke )
      {
        dCompPerfRate_dPres[iperf][ke][ic] = 0.0;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompPerfRate_dComp[iperf][ke][ic][jc] = 0.0;
        }
      }
    }

    // clear working arrays
    pressure = 0.0;
    dPressure_dP = 0.0;
    dPressure_dC = 0.0;

    // 1) copy the variables from the reservoir and well element

    // a) get reservoir variables

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];
    // get the index of the well elem
    localIndex const iwelem = perfWellElemIndex[iperf];

    pressure[SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[SubRegionTag::RES] = 1.0;

    // TODO: add a buoyancy term for the reservoir side here

    multiplier[SubRegionTag::RES] = 1.0;

    // b) get well variables

    pressure[SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    dPressure_dP[SubRegionTag::WELL] = 1.0;

    multiplier[SubRegionTag::WELL] = -1.0;

    real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );

    pressure[SubRegionTag::WELL] += wellElemMixtureDensity[iwelem] * gravD;
    dPressure_dP[SubRegionTag::WELL] += dWellElemMixtureDensity_dPres[iwelem] * gravD;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dPressure_dC[SubRegionTag::WELL][ic] +=
        dWellElemMixtureDensity_dComp[iwelem][ic] * gravD;
    }

    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf];

    // 2) compute potential difference

    real64 potDiff = 0.0;
    dPotDiff_dP = 0.0;
    dPotDiff_dC = 0.0;

    for( localIndex i = 0; i < 2; ++i )
    {
      potDiff += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
      dPotDiff_dP[i] += multiplier[i] * trans * dPressure_dP[i];

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPotDiff_dC[i][ic] += multiplier[i] * trans * dPressure_dC[i][ic];
      }
    }

    real64 flux = 0.0;
    dFlux_dP = 0.0;
    dFlux_dC = 0.0;

    // 3) upwinding

    if( potDiff >= 0 )  // ** reservoir cell is upstream **
    {

      phaseCompFrac = 0.0;
      dPhaseCompFrac_dP = 0.0;
      dPhaseCompFrac_dC = 0.0;

      dPhaseCompFrac_dCompDens = 0;

      // loop over phases, compute and upwind phase flux
      // and sum contributions to each component's perforation rate
      for( localIndex ip = 0; ip < NP; ++ip )
      {

        // compute the phase flux and derivatives using upstream cell mobility
        flux = resPhaseMob[er][esr][ei][ip] * potDiff;

        dFlux_dP[SubRegionTag::RES] =
          dResPhaseMob_dPres[er][esr][ei][ip] * potDiff
          + resPhaseMob[er][esr][ei][ip] * dPotDiff_dP[SubRegionTag::RES];

        dFlux_dP[SubRegionTag::WELL] =
          resPhaseMob[er][esr][ei][ip] * dPotDiff_dP[SubRegionTag::WELL];

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dFlux_dC[SubRegionTag::RES][ic] =
            dResPhaseMob_dComp[er][esr][ei][ip][ic] * potDiff
            + resPhaseMob[er][esr][ei][ip] * dPotDiff_dC[SubRegionTag::RES][ic];

          dFlux_dC[SubRegionTag::WELL][ic] =
            resPhaseMob[er][esr][ei][ip] * dPotDiff_dC[SubRegionTag::RES][ic];
        }

        // increment component fluxes
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          compPerfRate[iperf][ic] += flux * resPhaseCompFrac[er][esr][ei][0][ip][ic];

          dCompPerfRate_dPres[iperf][SubRegionTag::RES][ic] +=
            resPhaseCompFrac[er][esr][ei][0][ip][ic] * dFlux_dP[SubRegionTag::RES];

          dCompPerfRate_dPres[iperf][SubRegionTag::RES][ic] +=
            dResPhaseCompFrac_dPres[er][esr][ei][0][ip][ic] * flux;

          dCompPerfRate_dPres[iperf][SubRegionTag::WELL][ic] +=
            resPhaseCompFrac[er][esr][ei][0][ip][ic] * dFlux_dP[SubRegionTag::WELL];

          applyChainRule( NC,
                          dResCompFrac_dCompDens[er][esr][ei],
                          dResPhaseCompFrac_dComp[er][esr][ei][0][ip][ic],
                          dPhaseCompFrac_dCompDens );

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCompPerfRate_dComp[iperf][SubRegionTag::RES][ic][jc] +=
              dFlux_dC[SubRegionTag::RES][jc]
              * resPhaseCompFrac[er][esr][ei][0][ip][ic];

            dCompPerfRate_dComp[iperf][SubRegionTag::RES][ic][jc] +=
              flux * dPhaseCompFrac_dCompDens[jc];

            dCompPerfRate_dComp[iperf][SubRegionTag::WELL][ic][jc] +=
              dFlux_dC[SubRegionTag::WELL][jc]
              * resPhaseCompFrac[er][esr][ei][0][ip][ic];
          }
        }
      }
    }
    else // ** well is upstream **
    {

      real64 resTotalMobility     = 0.0;
      real64 dResTotalMobility_dP = 0.0;
      dResTotalMobility_dC        = 0.0;

      // first, compute the reservoir total mobitity (excluding phase density)
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        // viscosity
        real64 const resViscosity = resPhaseVisc[er][esr][ei][0][ip];
        real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][ei][0][ip];
        applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei],
                        dResPhaseVisc_dComp[er][esr][ei][0][ip],
                        dVisc_dC );

        // relative permeability
        real64 const resRelPerm = resPhaseRelPerm[er][esr][ei][0][ip];
        real64 dResRelPerm_dP = 0.0;
        dRelPerm_dC = 0.0;
        for( localIndex jp = 0; jp < NP; ++jp )
        {
          real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[er][esr][ei][jp];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[er][esr][ei][jp][jc];
          }
        }

        // increment total mobility
        resTotalMobility     += resRelPerm / resViscosity;
        dResTotalMobility_dP += ( dResRelPerm_dP * resViscosity - resRelPerm * dResVisc_dP )
                                / ( resViscosity * resViscosity );
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dResTotalMobility_dC[ic] += ( dRelPerm_dC[ic] * resViscosity - resRelPerm * dVisc_dC[ic] )
                                      / ( resViscosity * resViscosity );
        }
      }

      // compute a potdiff multiplier = wellElemMixtureDensity * resTotalMobility
      real64 const mult = wellElemMixtureDensity[iwelem] * resTotalMobility;
      dMult_dP[SubRegionTag::RES]  =
        wellElemMixtureDensity[iwelem] * dResTotalMobility_dP;
      dMult_dP[SubRegionTag::WELL] =
        dWellElemMixtureDensity_dPres[iwelem] * resTotalMobility;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dMult_dC[SubRegionTag::RES][ic]  =
          wellElemMixtureDensity[iwelem] * dResTotalMobility_dC[ic];

        dMult_dC[SubRegionTag::WELL][ic] =
          dWellElemMixtureDensity_dComp[iwelem][ic] * resTotalMobility;
      }

      // compute the volumetric flux and derivatives using upstream cell mobility
      flux = mult * potDiff;

      dFlux_dP[SubRegionTag::RES] =
        dMult_dP[SubRegionTag::RES] * potDiff
        + mult * dPotDiff_dP[SubRegionTag::RES];

      dFlux_dP[SubRegionTag::WELL] =
        dMult_dP[SubRegionTag::WELL] * potDiff
        + mult * dPotDiff_dP[SubRegionTag::WELL];

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dFlux_dC[SubRegionTag::RES][ic] =
          dMult_dC[SubRegionTag::RES][ic] * potDiff
          + mult * dPotDiff_dC[SubRegionTag::RES][ic];

        dFlux_dC[SubRegionTag::WELL][ic] =
          dMult_dC[SubRegionTag::WELL][ic] * potDiff
          + mult * dPotDiff_dC[SubRegionTag::RES][ic];
      }

      // compute component fluxes
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compPerfRate[iperf][ic] += wellElemCompFrac[iwelem][ic] * flux;

        dCompPerfRate_dPres[iperf][SubRegionTag::RES][ic] =
          wellElemCompFrac[iwelem][ic] * dFlux_dP[SubRegionTag::RES];

        dCompPerfRate_dPres[iperf][SubRegionTag::WELL][ic] =
          wellElemCompFrac[iwelem][ic] * dFlux_dP[SubRegionTag::WELL];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompPerfRate_dComp[iperf][SubRegionTag::RES][ic][jc]  +=
            wellElemCompFrac[iwelem][ic] * dFlux_dC[SubRegionTag::RES][jc];

          dCompPerfRate_dComp[iperf][SubRegionTag::WELL][ic][jc] +=
            wellElemCompFrac[iwelem][ic] * dFlux_dC[SubRegionTag::WELL][jc];

          dCompPerfRate_dComp[iperf][SubRegionTag::WELL][ic][jc] +=
            dWellElemCompFrac_dCompDens[iwelem][ic][jc] * flux;
        }
      }
    }
  }
}


void
CompositionalMultiphaseWell::ApplySystemSolution( DofManager const & dofManager,
                                                  ParallelVector const & solution,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  dofManager.addVectorToField( solution,
                               WellElementDofName(),
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( solution,
                               WellElementDofName(),
                               viewKeyStruct::deltaGlobalCompDensityString,
                               scalingFactor,
                               1, m_numDofPerWellElement - 1 );

  dofManager.addVectorToField( solution,
                               WellElementDofName(),
                               viewKeyStruct::deltaMixtureConnRateString,
                               scalingFactor,
                               m_numDofPerWellElement - 1, m_numDofPerWellElement );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaGlobalCompDensityString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaMixtureConnRateString );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain->getNeighbors() );

  // update properties
  UpdateStateAll( *domain );

}

void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{

  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

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

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {
      // extract solution and apply to dP
      dWellElemPressure[iwelem] = 0;
      dConnRate[iwelem] = 0;
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        dWellElemGlobalCompDensity[iwelem][ic] = 0;
      }
    } );
  } );

  // call constitutive models
  UpdateStateAll( *domain );
}


void CompositionalMultiphaseWell::ResetViews( DomainPartition * const domain )
{
  WellSolverBase::ResetViews( domain );

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager & elemManager = *mesh.getElemManager();

  CompositionalMultiphaseFlow & flowSolver = *getParent()->GetGroup< CompositionalMultiphaseFlow >( GetFlowSolverName() );

  {
    using keys = CompositionalMultiphaseFlow::viewKeyStruct;

    m_resPressure =
      elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::pressureString );

    m_deltaResPressure =
      elemManager.ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( keys::deltaPressureString );

    m_resGlobalCompDensity =
      elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::globalCompDensityString );

    m_dResCompFrac_dCompDens =
      elemManager.ConstructViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dGlobalCompFraction_dGlobalCompDensityString );

    m_resPhaseMob =
      elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::phaseMobilityString );

    m_dResPhaseMob_dPres =
      elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dPhaseMobility_dPressureString );

    m_dResPhaseMob_dCompDens =
      elemManager.ConstructViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dPhaseMobility_dGlobalCompDensityString );

    m_resPhaseVolFrac =
      elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::phaseVolumeFractionString );

    m_dResPhaseVolFrac_dPres =
      elemManager.ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dPhaseVolumeFraction_dPressureString );

    m_dResPhaseVolFrac_dCompDens =
      elemManager.ConstructViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dPhaseVolumeFraction_dGlobalCompDensityString );
  }
  {
    using keys = MultiFluidBase::viewKeyStruct;

    m_resPhaseDens =
      elemManager.ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseDensityString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_resPhaseVisc =
      elemManager.ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseViscosityString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseVisc_dPres =
      elemManager.ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dPhaseViscosity_dPressureString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseVisc_dComp =
      elemManager.ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseViscosity_dGlobalCompFractionString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_resPhaseCompFrac =
      elemManager.ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::phaseCompFractionString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseCompFrac_dPres =
      elemManager.ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseCompFraction_dPressureString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseCompFrac_dComp =
      elemManager.ConstructMaterialViewAccessor< array5d< real64 >, arrayView5d< real64 > >( keys::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
  }
  {
    using keys = RelativePermeabilityBase::viewKeyStruct;

    m_resPhaseRelPerm =
      elemManager.ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseRelPermString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.relPermModelNames() );
    m_dResPhaseRelPerm_dPhaseVolFrac =
      elemManager.ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.relPermModelNames() );
  }
}


void CompositionalMultiphaseWell::FormPressureRelations( DomainPartition const * const domain,
                                                         DofManager const * const dofManager,
                                                         ParallelMatrix * const matrix,
                                                         ParallelVector * const rhs )
{

  MeshLevel const & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1; // here, dofs are 1 pressure, NC comp densities (reservoir and well)

  localIndex const resNDOF  = NumDofPerResElement();
  localIndex const NC = m_numComponents;

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {

    WellControls const & wellControls = GetWellControls( subRegion );

    // get the degrees of freedom, depth info, next welem index
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    // get secondary data on well elements
    arrayView1d< real64 const > const & wellElemMixtureDensity =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureDensityString );
    arrayView1d< real64 const > const & dWellElemMixtureDensity_dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::dMixtureDensity_dPressureString );
    arrayView2d< real64 const > const & dWellElemMixtureDensity_dComp =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

    real64 const targetBHP = wellControls.GetTargetBHP();

    // compute a coefficient to normalize the momentum equation
    real64 const normalizer = targetBHP > std::numeric_limits< real64 >::epsilon()
                              ? 1.0 / targetBHP
                              : 1.0;

    // loop over the well elements to compute the pressure relations between well elements
    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {

      if( wellElemGhostRank[iwelem] >= 0 )
      {
        return;
      }

      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      if( iwelemNext >= 0 )   // if iwelemNext < 0, form control equation, not momentum
      {

        // local working variables and arrays
        stackArray1d< globalIndex, 2 * maxNumDof > dofColIndices( 2 * resNDOF );
        stackArray1d< real64, 2 * maxNumDof > localPresRelJacobian( 2 * resNDOF );

        stackArray1d< real64, maxNumComp > dAvgDensity_dCompCurrent( NC );
        stackArray1d< real64, maxNumComp > dAvgDensity_dCompNext( NC );

        // reset local working variables and arrays
        dofColIndices = -1;
        localPresRelJacobian = 0.0;

        // compute the average density at the interface between well elements
        real64 const avgDensity = 0.5 * ( wellElemMixtureDensity[iwelemNext]
                                          + wellElemMixtureDensity[iwelem] );
        real64 const dAvgDensity_dPresNext = 0.5 * dWellElemMixtureDensity_dPres[iwelemNext];
        real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemMixtureDensity_dPres[iwelem];
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dAvgDensity_dCompNext[ic] = 0.5 * dWellElemMixtureDensity_dComp[iwelemNext][ic];
          dAvgDensity_dCompCurrent[ic] = 0.5 * dWellElemMixtureDensity_dComp[iwelem][ic];
        }

        // compute depth diff times acceleration
        real64 const gravD = wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem];

        // compute the current pressure in the two well elements
        real64 const pressureNext = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];
        real64 const pressureCurrent = wellElemPressure[iwelem] + dWellElemPressure[iwelem];

        // compute momentum flux and derivatives
        localIndex const localDofIndexPresNext = ElemTag::NEXT * resNDOF;
        localIndex const localDofIndexPresCurrent = ElemTag::CURRENT * resNDOF;

        globalIndex const offsetNext = wellElemDofNumber[iwelemNext];
        globalIndex const offsetCurrent = wellElemDofNumber[iwelem];

        globalIndex const eqnRowIndex = offsetCurrent + RowOffset::CONTROL;

        dofColIndices[localDofIndexPresNext] = offsetNext + ColOffset::DPRES;
        dofColIndices[localDofIndexPresCurrent] = offsetCurrent + ColOffset::DPRES;

        real64 const localPresRel = ( pressureNext - pressureCurrent - avgDensity * gravD ) * normalizer;

        localPresRelJacobian[localDofIndexPresNext] = ( 1 - dAvgDensity_dPresNext * gravD ) * normalizer;
        localPresRelJacobian[localDofIndexPresCurrent] = ( -1 - dAvgDensity_dPresCurrent * gravD ) * normalizer;

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          localIndex const localDofIndexCompNext = localDofIndexPresNext + ic + 1;
          localIndex const localDofIndexCompCurrent = localDofIndexPresCurrent + ic + 1;

          dofColIndices[localDofIndexCompNext] = offsetNext + ColOffset::DCOMP + ic;
          dofColIndices[localDofIndexCompCurrent] = offsetCurrent + ColOffset::DCOMP + ic;

          localPresRelJacobian[localDofIndexCompNext] = -dAvgDensity_dCompNext[ic] * gravD * normalizer;
          localPresRelJacobian[localDofIndexCompCurrent] = -dAvgDensity_dCompCurrent[ic] * gravD * normalizer;
        }

        // TODO: add friction and acceleration terms

        rhs->add( eqnRowIndex,
                  localPresRel );

        matrix->add( eqnRowIndex,
                     dofColIndices,
                     localPresRelJacobian );
      }
    } );
  } );
}

void CompositionalMultiphaseWell::FormControlEquation( DomainPartition const * const domain,
                                                       DofManager const * const dofManager,
                                                       ParallelMatrix * const matrix,
                                                       ParallelVector * const rhs )
{
  MeshLevel const & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  string const wellDofKey = dofManager->getKey( WellElementDofName() );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {

    if( !subRegion.IsLocallyOwned() )
    {
      return;
    }

    WellControls const & wellControls = GetWellControls( subRegion );

    // get the degrees of freedom
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get the index of the well element where the control is enforced
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    // get well control
    WellControls::Control const control = wellControls.GetControl();

    // BHP control
    if( control == WellControls::Control::BHP )
    {
      // get primary variables on well elements
      arrayView1d< real64 const > const & wellElemPressure =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
      arrayView1d< real64 const > const & dWellElemPressure =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

      // get the pressure and compute normalizer
      real64 const currentBHP = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];
      real64 const targetBHP = wellControls.GetTargetBHP();
      real64 const normalizer = targetBHP > std::numeric_limits< real64 >::epsilon()
                                ? 1.0 / targetBHP
                                : 1.0;

      // control equation is a normalized difference
      // between current pressure and target pressure
      real64 const controlEqn = ( currentBHP - targetBHP ) * normalizer;
      real64 const dControlEqn_dPres = normalizer;

      globalIndex const elemOffset = wellElemDofNumber[iwelemControl];
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DPRES;

      rhs->add( eqnRowIndex,
                controlEqn );

      matrix->add( eqnRowIndex,
                   dofColIndex,
                   dControlEqn_dPres );

    }
    else if( control == WellControls::Control::LIQUIDRATE ) // liquid rate control
    {
      localIndex const NC = m_numComponents;

      // get a reference to the primary variables on well elements
      arrayView1d< real64 const > const & connRate =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );
      arrayView1d< real64 const > const & dConnRate =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString );

      // get rates and compute normalizer
      real64 const currentConnRate = connRate[iwelemControl] + dConnRate[iwelemControl];
      real64 const targetConnRate = wellControls.GetTargetRate();
      real64 const normalizer = targetConnRate > std::numeric_limits< real64 >::min()
                                ? 1.0 / ( 1e-2 * targetConnRate ) // hard-coded value comes from AD-GPRS
                                : 1.0;

      // control equation is a normalized difference
      // between current rate and target rate
      real64 const controlEqn = ( currentConnRate - targetConnRate ) * normalizer;
      real64 const dControlEqn_dRate = normalizer;

      globalIndex const elemOffset = wellElemDofNumber[iwelemControl];
      localIndex const dRateColOffset = ColOffset::DCOMP + NC;
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + dRateColOffset;

      rhs->add( eqnRowIndex,
                controlEqn );

      matrix->add( eqnRowIndex,
                   dofColIndex,
                   dControlEqn_dRate );

    }
    else
    {
      GEOSX_ERROR_IF( ( control != WellControls::Control::BHP ) && ( control != WellControls::Control::LIQUIDRATE ),
                      "Phase rate constraints for CompositionalMultiphaseWell will be implemented later" );
    }
  } );
}


void CompositionalMultiphaseWell::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition * const domain )
{

  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

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

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        wellElemGlobalCompDensity[iwelem][ic] += dWellElemGlobalCompDensity[iwelem][ic];
      }
      connRate[iwelem] += dConnRate[iwelem];
    } );
  } );
}


void CompositionalMultiphaseWell::RecordWellData( WellElementSubRegion const & subRegion )
{
  // note: this function is for debug and will go away

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;

  WellControls const & wellControls = GetWellControls( subRegion );
  PerforationData const * const perforationData = subRegion.GetPerforationData();

  arrayView1d< real64 const > const & wellElemPressure =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );

  arrayView1d< real64 const > const & connRate =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );

  arrayView2d< real64 const > const & wellElemCompFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );

  arrayView2d< real64 const > const & compPerfRate =
    perforationData->getReference< array2d< real64 > >( viewKeyStruct::compPerforationRateString );

  // here, we will save the well info
  // for now, brute force: output to terminal

  std::cout << "Well : " << wellControls.getName() << std::endl;
  if( wellControls.GetType() == WellControls::Type::PRODUCER )
    std::cout << "Type : PRODUCER" << std::endl;
  else
    std::cout << "Type : INJECTOR" << std::endl;
  if( wellControls.GetControl() == WellControls::Control::BHP )
    std::cout << "Control : BHP" << std::endl;
  else
    std::cout << "Control : RATE" << std::endl;

  std::cout << "Below, positive perforation rate means flow from reservoir to well" << std::endl;
  std::cout << "Negative perforation rate means flow from well to reservoir" << std::endl;

  // output perforation rates
  for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
  {
    for( localIndex ic = 0; ic < NC; ++ic )
      std::cout << "Mass rate at perforation #" << iperf << " for component #" << ic << ": " << compPerfRate[iperf][ic] << std::endl;
  }

  // output the reference pressure
  if( subRegion.IsLocallyOwned())
  {
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();
    real64 const pressure = wellElemPressure[iwelemControl];
    real64 const targetPressure = wellControls.GetTargetBHP();

    if( wellControls.GetControl() == WellControls::Control::BHP )
    {
      std::cout << "Current reference pressure = " << pressure
                << ", targetPressure = "           << targetPressure
                << std::endl;
    }
    else
    {
      if( wellControls.GetType() == WellControls::Type::PRODUCER )
      {
        std::cout << "Current reference pressure = " << pressure
                  << ", min pressure = " << targetPressure
                  << std::endl;
      }
      else
      {
        std::cout << "Current reference pressure = " << pressure
                  << ", max pressure = " << targetPressure
                  << std::endl;
      }
    }
  }

  std::cout << "Below, negative connection rate means production" << std::endl;
  std::cout << "Positive connection rate means injection" << std::endl;

  for( localIndex iwelem = 0; iwelem < subRegion.size(); ++iwelem )
  {
    if( iwelem > 0 || wellControls.GetControl() == WellControls::Control::BHP )
      std::cout << "Mass rate at connection #" << iwelem << ": " << connRate[iwelem] << std::endl;
    else
    {
      std::cout << "Mass rate at connection #" << iwelem << ": " << connRate[iwelem]
                << ", target rate : " << wellControls.GetTargetRate() << std::endl;
    }

    if( iwelem > 0 )
    {
      continue;
    }

    stackArray1d< real64, maxNumComp > compStream( NC );

    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const compFrac = wellElemCompFrac[iwelem][ic];
      compStream[ic] = compFrac * connRate[iwelem];
    }

    if( wellControls.GetType() == WellControls::Type::INJECTOR )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
        std::cout << "Mass injection rate for component #" << ic
                  << " at connection #" << iwelem
                  << ": " << compStream[ic] << std::endl;
    }
    else
    {
      for( localIndex ic = 0; ic < NC; ++ic )
        std::cout << "Mass production rate for component #" << ic
                  << " at connection #" << iwelem
                  << ": " << compStream[ic] << std::endl;
    }
  }
}


void CompositionalMultiphaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {

    if( !subRegion.IsLocallyOwned() )
    {
      return;
    }

    WellControls & wellControls = GetWellControls( subRegion );

    // get the primary variables
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString );

    // if isViable is true at the end of the following checks, no need to switch
    bool controlIsViable = false;

    // get well control and type
    WellControls::Control const currentControl = wellControls.GetControl();
    WellControls::Type const type = wellControls.GetType();

    // get the index of the well element where the control is enforced
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    real64 const refRate = connRate[iwelemControl] + dConnRate[iwelemControl];
    real64 const refPressure = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];

    // TODO: check all inactive constraints (possibly more than one) and switch the one which is most violated
    // TODO: for the rate, use surface conditions (flash for compositional, easier for BO)

    // BHP control
    if( currentControl == WellControls::Control::BHP )
    {
      // the control is viable if the reference rate is below/above the max/min rate
      // targetRate specifies a max rate here
      real64 const & maxRate = wellControls.GetTargetRate();
      controlIsViable = ( fabs( refRate ) <= fabs( maxRate ) );
    }
    else // rate control
    {

      // the control is viable if the reference pressure is below/above the max/min pressure
      if( type == WellControls::Type::PRODUCER )
      {
        // targetBHP specifies a min pressure here
        real64 const & minPressure = wellControls.GetTargetBHP();
        controlIsViable = ( refPressure >= minPressure );
      }
      else
      {
        // targetBHP specifies a max pressure here
        real64 const & maxPressure = wellControls.GetTargetBHP();
        controlIsViable = ( refPressure <= maxPressure );
      }
    }

    if( !controlIsViable )
    {
      if( currentControl == WellControls::Control::BHP )
      {
        wellControls.SetControl( WellControls::Control::LIQUIDRATE,
                                 wellControls.GetTargetRate() );

        // Debug information for logLevel >= 1
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from BHP constraint to rate constraint" );

      }
      else // rate control
      {
        wellControls.SetControl( WellControls::Control::BHP,
                                 wellControls.GetTargetBHP() );
        // Debug information for logLevel >= 1
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from rate constraint to BHP constraint" );
      }
    }
  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseWell, string const &, Group * const )
}// namespace geosx
