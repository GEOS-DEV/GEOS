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
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "managers/DomainPartition.hpp"
#include "wells/PerforationData.hpp"
#include "wells/WellElementSubRegion.hpp"
#include "wells/WellControls.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseWellKernels.hpp"

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
  GEOSX_MARK_FUNCTION;

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

  MixtureDensityKernel::Launch< serialPolicy >( subRegion.size(),
                                                NumFluidComponents(),
                                                NumFluidPhases(),
                                                wellElemPhaseVolFrac,
                                                dWellElemPhaseVolFrac_dPres,
                                                dWellElemPhaseVolFrac_dComp,
                                                dWellElemCompFrac_dCompDens,
                                                wellElemPhaseDens,
                                                dWellElemPhaseDens_dPres,
                                                dWellElemPhaseDens_dComp,
                                                wellElemMixtureDensity,
                                                dWellElemMixtureDensity_dPres,
                                                dWellElemMixtureDensity_dComp );

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

void CompositionalMultiphaseWell::InitializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & resPressure = m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resCompDens = m_resGlobalCompDensity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resPhaseVolFrac = m_resPhaseVolFrac;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & resPhaseDens    = m_resPhaseDens;

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

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
    constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      CompositionalMultiphaseFlowKernels::FluidUpdateKernel::Launch< serialPolicy >( subRegion.size(),
                                                                                     fluidWrapper,
                                                                                     wellElemPressure,
                                                                                     m_temperature,
                                                                                     wellElemCompFrac );
    } );

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        wellElemCompDens[iwelem][ic] = avgCompFrac[ic] * totalDens[iwelem][0];
      }
    } );

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
                                                     DomainPartition const & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

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

    FluxKernel::Launch< serialPolicy >( subRegion.size(),
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

    VolumeBalanceKernel::Launch< serialPolicy >( subRegion.size(),
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
    arrayView1d< integer const > const & wellElemGhostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
    arrayView1d< real64 const > const & wellElemVolume =
      subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

    MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    ResidualNormKernel::Launch< serialPolicy, serialReduce >( localRhs,
                                                              dofManager.rankOffset(),
                                                              NumFluidComponents(),
                                                              NumDofPerWellElement(),
                                                              wellElemDofNumber,
                                                              wellElemGhostRank,
                                                              wellElemVolume,
                                                              totalDens,
                                                              &localResidualNorm );
  } );
  return sqrt( MpiWrapper::Sum( localResidualNorm, MPI_COMM_GEOSX ) );
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

    localIndex const subRegionSolutionCheck =
      SolutionCheckKernel::Launch< serialPolicy, serialReduce >( localSolution,
                                                                 dofManager.rankOffset(),
                                                                 NumFluidComponents(),
                                                                 wellElemDofNumber,
                                                                 wellElemGhostRank,
                                                                 wellElemPressure,
                                                                 dWellElemPressure,
                                                                 wellElemCompDens,
                                                                 dWellElemCompDens,
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

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & resPressure = m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const & dResPressure = m_deltaResPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & resPhaseMob = m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dResPhaseMob_dPres = m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResPhaseMob_dComp = m_dResPhaseMob_dCompDens;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & dResPhaseVolFrac_dPres = m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResPhaseVolFrac_dComp = m_dResPhaseVolFrac_dCompDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResCompFrac_dCompDens = m_dResCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & resPhaseVisc = m_resPhaseVisc;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & dResPhaseVisc_dPres = m_dResPhaseVisc_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & dResPhaseVisc_dComp = m_dResPhaseVisc_dComp;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & resPhaseCompFrac = m_resPhaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & dResPhaseCompFrac_dPres = m_dResPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 > > const & dResPhaseCompFrac_dComp = m_dResPhaseCompFrac_dComp;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 > > const & resPhaseRelPerm = m_resPhaseRelPerm;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 > > const & dResPhaseRelPerm_dPhaseVolFrac = m_dResPhaseRelPerm_dPhaseVolFrac;

  PerforationData * const perforationData = subRegion.GetPerforationData();

  // get depth
  arrayView1d< real64 const > const & wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // get well primary variables on well elements
  arrayView1d< real64 const > const & wellElemPressure =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dWellElemPressure =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  // get secondary well data on well elements
  arrayView1d< real64 const > const & wellElemMixtureDensity =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureDensityString );
  arrayView1d< real64 const > const & dWellElemMixtureDensity_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dMixtureDensity_dPressureString );
  arrayView2d< real64 const > const & dWellElemMixtureDensity_dComp =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

  arrayView2d< real64 const > const & wellElemCompFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
  arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // get well variables on perforations
  arrayView1d< real64 const > const & perfGravCoef =
    perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );
  arrayView1d< localIndex const > const & perfWellElemIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );
  arrayView1d< real64 const > const & perfTransmissibility =
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


  PerforationKernel::Launch< serialPolicy >( perforationData->size(),
                                             NumFluidComponents(),
                                             NumFluidPhases(),
                                             resPressure,
                                             dResPressure,
                                             resPhaseMob,
                                             dResPhaseMob_dPres,
                                             dResPhaseMob_dComp,
                                             dResPhaseVolFrac_dPres,
                                             dResPhaseVolFrac_dComp,
                                             dResCompFrac_dCompDens,
                                             resPhaseVisc,
                                             dResPhaseVisc_dPres,
                                             dResPhaseVisc_dComp,
                                             resPhaseCompFrac,
                                             dResPhaseCompFrac_dPres,
                                             dResPhaseCompFrac_dComp,
                                             resPhaseRelPerm,
                                             dResPhaseRelPerm_dPhaseVolFrac,
                                             wellElemGravCoef,
                                             wellElemPressure,
                                             dWellElemPressure,
                                             wellElemMixtureDensity,
                                             dWellElemMixtureDensity_dPres,
                                             dWellElemMixtureDensity_dComp,
                                             wellElemCompFrac,
                                             dWellElemCompFrac_dCompDens,
                                             perfGravCoef,
                                             perfWellElemIndex,
                                             perfTransmissibility,
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

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaGlobalCompDensityString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaMixtureConnRateString ) );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors() );

  // update properties
  UpdateStateAll( domain );

}

void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition & domain )
{

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

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

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const iwelem )
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

    WellControls const & wellControls = GetWellControls( subRegion );

    // get the degrees of freedom, depth info, next welem index
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
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

    PressureRelationKernel::Launch< serialPolicy >( subRegion.size(),
                                                    dofManager.rankOffset(),
                                                    NumFluidComponents(),
                                                    NumDofPerResElement(),
                                                    wellControls,
                                                    wellElemDofNumber,
                                                    wellElemGravCoef,
                                                    nextWellElemIndex,
                                                    wellElemPressure,
                                                    dWellElemPressure,
                                                    wellElemMixtureDensity,
                                                    dWellElemMixtureDensity_dPres,
                                                    dWellElemMixtureDensity_dComp,
                                                    localMatrix,
                                                    localRhs );

  } );
}

void CompositionalMultiphaseWell::FormControlEquation( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

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
    string const wellDofKey = dofManager.getKey( WellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get the index of the well element where the control is enforced
    localIndex const iwelemControl = wellControls.GetReferenceWellElementIndex();

    // get primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString );

    ControlEquationHelper::ComputeJacobianEntry( dofManager.rankOffset(),
                                                 m_numComponents,
                                                 wellControls,
                                                 wellElemDofNumber[iwelemControl],
                                                 wellElemPressure[iwelemControl],
                                                 dWellElemPressure[iwelemControl],
                                                 connRate[iwelemControl],
                                                 dConnRate[iwelemControl],
                                                 localMatrix,
                                                 localRhs );
  } );
}


void CompositionalMultiphaseWell::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition & domain )
{

  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

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


void CompositionalMultiphaseWell::CheckWellControlSwitch( DomainPartition & domain )
{
  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

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
