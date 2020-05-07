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
 * @file CompositionalMultiphaseFlow.cpp
 */

#include "CompositionalMultiphaseFlow.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "dataRepository/Group.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseFlowKernels;

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow( const string & name,
                                                          Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_capPressureFlag( 0 )
{
  // set the blockID for the block system interface
  // To generate the schema, multiple solvers of that use this command are constructed
  // Doing this can cause an error in the block setup, so move it to InitializePreSubGroups
  // getLinearSystemRepository()->SetBlockID(BlockIDs::compositionalBlock, this->getName());

//START_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::temperatureString, &m_temperature )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Temperature" );

  this->registerWrapper( viewKeyStruct::useMassFlagString, &m_useMass )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::relPermNamesString, &m_relPermModelNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setSizedFromParent( 0 )->
    setDescription( "Name of the relative permeability constitutive model to use" );

  this->registerWrapper( viewKeyStruct::capPressureNamesString, &m_capPressureModelNames )->
    setApplyDefaultValue( "" )->
    setSizedFromParent( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the capillary pressure constitutive model to use" );
}

void CompositionalMultiphaseFlow::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();
  CheckModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString );
  m_capPressureFlag = CheckModelNames( m_capPressureModelNames, viewKeyStruct::capPressureNamesString, true );
}

void CompositionalMultiphaseFlow::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    forTargetSubRegions( meshLevel, [&]( localIndex const,
                                         ElementSubRegionBase & elementSubRegion )
    {

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )->setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::bcPressureString );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompDensityString )->setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompFractionString )->setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString )->setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseMobilityString )->setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionOldString );
      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseDensityOldString );
      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::phaseComponentFractionOldString );
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::porosityOldString );
    } );
  }
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

void CompositionalMultiphaseFlow::ValidateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const
{
  MultiFluidBase const & fluid0 = *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );

  for( localIndex i = 1; i < m_fluidModelNames.size(); ++i )
  {
    MultiFluidBase const & fluid = *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[i] );
    CompareMultiphaseModels( fluid, fluid0 );
    CompareMulticomponentModels( fluid, fluid0 );
  }

  RelativePermeabilityBase const & relPerm0 = *cm.GetConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[0] );
  CompareMultiphaseModels( relPerm0, fluid0 );

  for( localIndex i = 1; i < m_relPermModelNames.size(); ++i )
  {
    RelativePermeabilityBase const & relPerm = *cm.GetConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[i] );
    CompareMultiphaseModels( relPerm, relPerm0 );
  }

  if( m_capPressureFlag )
  {
    CapillaryPressureBase const & capPres0 = *cm.GetConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[0] );
    CompareMultiphaseModels( capPres0, fluid0 );

    for( localIndex i = 1; i < m_capPressureModelNames.size(); ++i )
    {
      CapillaryPressureBase const & capPres = *cm.GetConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[i] );
      CompareMultiphaseModels( capPres, capPres0 );
    }
  }
}


void CompositionalMultiphaseFlow::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  ConstitutiveManager const & cm = *domain->getConstitutiveManager();

  // 1. Set key dimensions of the problem
  MultiFluidBase const & fluid0 = *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
  m_numPhases     = fluid0.numFluidPhases();
  m_numComponents = fluid0.numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  // 2. Validate various models against each other (must have same phases and components)
  ValidateConstitutiveModels( cm );

  // 3. Resize all fields as necessary, validate constitutive models in regions
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ResizeFields( meshLevel );

    ValidateModelMapping< MultiFluidBase >( *meshLevel.getElemManager(), m_fluidModelNames );
    ValidateModelMapping< RelativePermeabilityBase >( *meshLevel.getElemManager(), m_relPermModelNames );
    if( m_capPressureFlag )
    {
      ValidateModelMapping< CapillaryPressureBase >( *meshLevel.getElemManager(), m_capPressureModelNames );
    }
  }
}

void CompositionalMultiphaseFlow::ResizeFields( MeshLevel & meshLevel )
{
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  forTargetSubRegions( meshLevel, [&]( localIndex const,
                                       ElementSubRegionBase & subRegion )
  {
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString ).resizeDimension< 1 >( NC );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString ).resizeDimension< 1 >( NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString ).resizeDimension< 1 >( NC );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString ).resizeDimension< 1, 2 >( NC, NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString ).resizeDimension< 1 >( NP );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString ).resizeDimension< 1 >( NP );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString ).resizeDimension< 1, 2 >( NP, NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString ).resizeDimension< 1 >( NP );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString ).resizeDimension< 1 >( NP );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString ).resizeDimension< 1, 2 >( NP, NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionOldString ).resizeDimension< 1 >( NP );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseDensityOldString ).resizeDimension< 1 >( NP );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::phaseComponentFractionOldString ).resizeDimension< 1, 2 >( NP, NC );
  } );
}

void CompositionalMultiphaseFlow::UpdateComponentFraction( Group & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & compFrac =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );

  arrayView3d< real64 > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const & compDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );

  arrayView2d< real64 const > const & dCompDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  KernelLaunchSelector1< ComponentFractionKernel >( m_numComponents,
                                                    dataGroup.size(),
                                                    compDens,
                                                    dCompDens,
                                                    compFrac,
                                                    dCompFrac_dCompDens );
}

void CompositionalMultiphaseFlow::UpdatePhaseVolumeFraction( Group & dataGroup,
                                                             localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseVolFrac =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );

  arrayView2d< real64 > const & dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 > const & dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  // inputs

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  arrayView2d< real64 const > const & compDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );

  arrayView2d< real64 const > const & dCompDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  KernelLaunchSelector2< PhaseVolumeFractionKernel >( m_numComponents, m_numPhases,
                                                      dataGroup.size(),
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

void CompositionalMultiphaseFlow::UpdatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseMob =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString );

  arrayView2d< real64 > const & dPhaseMob_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString );

  arrayView3d< real64 > const & dPhaseMob_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 const > const & dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseVisc = fluid.phaseViscosity();
  arrayView3d< real64 const > const & dPhaseVisc_dPres = fluid.dPhaseViscosity_dPressure();
  arrayView4d< real64 const > const & dPhaseVisc_dComp = fluid.dPhaseViscosity_dGlobalCompFraction();

  RelativePermeabilityBase const & relperm = GetConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseRelPerm = relperm.phaseRelPerm();
  arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac = relperm.dPhaseRelPerm_dPhaseVolFraction();

  KernelLaunchSelector2< PhaseMobilityKernel >( m_numComponents, m_numPhases,
                                                dataGroup.size(),
                                                dCompFrac_dCompDens,
                                                phaseDens,
                                                dPhaseDens_dPres,
                                                dPhaseDens_dComp,
                                                phaseVisc,
                                                dPhaseVisc_dPres,
                                                dPhaseVisc_dComp,
                                                phaseRelPerm,
                                                dPhaseRelPerm_dPhaseVolFrac,
                                                dPhaseVolFrac_dPres,
                                                dPhaseVolFrac_dComp,
                                                phaseMob,
                                                dPhaseMob_dPres,
                                                dPhaseMob_dComp );
}

void CompositionalMultiphaseFlow::UpdateFluidModel( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView1d< real64 const > const & pres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
  arrayView2d< real64 const > const & compFrac = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );

  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forAll< serialPolicy >( dataGroup.size(), [&] ( localIndex const a )
  {
    fluid.PointUpdate( pres[a] + dPres[a], m_temperature, compFrac[a], a, 0 );
  } );
  //fluid->BatchUpdate( pres, temp, compFrac );
}

void CompositionalMultiphaseFlow::UpdateSolidModel( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase & solid = GetConstitutiveModel< ConstitutiveBase >( dataGroup, m_solidModelNames[targetIndex] );

  arrayView1d< real64 const > const & pres  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  forAll< serialPolicy >( dataGroup.size(), [&] ( localIndex const a )
  {
    solid.StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  } );
}

void CompositionalMultiphaseFlow::UpdateRelPermModel( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  RelativePermeabilityBase & relPerm =
    GetConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermModelNames[targetIndex] );

  arrayView2d< real64 > const & phaseVolFrac =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );

  relPerm.BatchUpdate( phaseVolFrac );
}

void CompositionalMultiphaseFlow::UpdateCapPressureModel( Group & dataGroup, localIndex const targetIndex )
{
  if( m_capPressureFlag )
  {
    CapillaryPressureBase & capPressure =
      GetConstitutiveModel< CapillaryPressureBase >( dataGroup, m_capPressureModelNames[targetIndex] );

    arrayView2d< real64 > const & phaseVolFrac =
      dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );

    capPressure.BatchUpdate( phaseVolFrac );
  }
}

void CompositionalMultiphaseFlow::UpdateState( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  UpdateComponentFraction( dataGroup );
  UpdateFluidModel( dataGroup, targetIndex );
  UpdatePhaseVolumeFraction( dataGroup, targetIndex );
  UpdateSolidModel( dataGroup, targetIndex );
  UpdateRelPermModel( dataGroup, targetIndex );
  UpdatePhaseMobility( dataGroup, targetIndex );
  UpdateCapPressureModel( dataGroup, targetIndex );
}

void CompositionalMultiphaseFlow::InitializeFluidState( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const targetIndex,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase &,
                                    ElementSubRegionBase & subRegion )
  {
    // 1. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    UpdateFluidModel( subRegion, targetIndex );

    // 2. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    arrayView2d< real64 const > const & compFrac = m_compFrac[er][esr];
    arrayView2d< real64 const > const & totalDens = m_totalDens[er][esr];
    arrayView2d< real64 > compDens = m_globalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    } );

    // 3. Update dependent state quantities
    UpdatePhaseVolumeFraction( subRegion, targetIndex );
    UpdateSolidModel( subRegion, targetIndex );
    UpdateRelPermModel( subRegion, targetIndex );
    UpdatePhaseMobility( subRegion, targetIndex );
    UpdateCapPressureModel( subRegion, targetIndex );

  } );
}

void CompositionalMultiphaseFlow::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::globalCompDensityString );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  // set mass fraction flag on fluid models
  forTargetSubRegions( *mesh, [&]( localIndex const targetIndex,
                                   ElementSubRegionBase & subRegion )
  {
    MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    fluid.setMassFlag( m_useMass );
  } );

  // Initialize primary variables from applied initial conditions
  ResetViews( domain );
  InitializeFluidState( domain );
}

real64 CompositionalMultiphaseFlow::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                integer const cycleNumber,
                                                DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return;

  ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // currently the only method is implicit time integration
  dt_return = NonlinearImplicitStep( time_n, dt, cycleNumber, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void CompositionalMultiphaseFlow::BackupFields( DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // backup some fields used in time derivative approximation
  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase &,
                                    ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d< real64 const > const & poroRef = m_porosityRef[er][esr];

    arrayView2d< real64 const > const & phaseVolFrac = m_phaseVolFrac[er][esr];

    arrayView3d< real64 const > const & phaseDens = m_phaseDens[er][esr];
    arrayView4d< real64 const > const & phaseCompFrac = m_phaseCompFrac[er][esr];
    arrayView2d< real64 const > const & pvMult = m_pvMult[er][esr];

    arrayView2d< real64 > const & phaseDensOld = m_phaseDensOld[er][esr];
    arrayView2d< real64 > const & phaseVolFracOld = m_phaseVolFracOld[er][esr];
    arrayView3d< real64 > const & phaseCompFracOld = m_phaseCompFracOld[er][esr];
    arrayView1d< real64 > const & poroOld = m_porosityOld[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        for( localIndex ip = 0; ip < m_numPhases; ++ip )
        {
          phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
          phaseVolFracOld[ei][ip] = phaseVolFrac[ei][ip];

          for( localIndex ic = 0; ic < m_numComponents; ++ic )
          {
            phaseCompFracOld[ei][ip][ic] = phaseCompFrac[ei][0][ip][ic];
          }
        }
        poroOld[ei] = poroRef[ei] * pvMult[ei][0];
      }
    } );
  } );
}

void
CompositionalMultiphaseFlow::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                DomainPartition * const domain,
                                                DofManager & dofManager,
                                                ParallelMatrix & matrix,
                                                ParallelVector & rhs,
                                                ParallelVector & solution )
{
  // bind the stored views to the current domain
  ResetViews( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain );

  // setup dof numbers and linear system
  if( !m_coupledWellsFlag )
  {
    SetupSystem( domain, dofManager, matrix, rhs, solution );
  }
}

void CompositionalMultiphaseFlow::SetupDofs( DomainPartition const * const domain,
                                             DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::dofFieldString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const * const fvManager =
    numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::dofFieldString, fluxApprox );
}

void CompositionalMultiphaseFlow::AssembleSystem( real64 const time_n,
                                                  real64 const dt,
                                                  DomainPartition * const domain,
                                                  DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.open();
  rhs.open();

  AssembleAccumulationTerms( time_n, dt, domain, &dofManager, &matrix, &rhs );
  AssembleFluxTerms( time_n, dt, domain, &dofManager, &matrix, &rhs );
  AssembleVolumeBalanceTerms( time_n, dt, domain, &dofManager, &matrix, &rhs );

  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseFlow::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, LAIOutputFormat::MATRIX_MARKET );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, LAIOutputFormat::MATRIX_MARKET );

    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseFlow::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void CompositionalMultiphaseFlow::AssembleAccumulationTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                             real64 const GEOSX_UNUSED_PARAM( dt ),
                                                             DomainPartition const * const domain,
                                                             DofManager const * const dofManager,
                                                             ParallelMatrix * const matrix,
                                                             ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  string const dofKey = dofManager->getKey( viewKeyStruct::dofFieldString );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & volume = m_volume[er][esr];
    arrayView1d< real64 const > const & porosityRef = m_porosityRef[er][esr];

    arrayView2d< real64 const > const & phaseVolFrac = m_phaseVolFrac[er][esr];
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];
    arrayView3d< real64 const > const & dCompFrac_dCompDens = m_dCompFrac_dCompDens[er][esr];

    arrayView1d< real64 const > const & porosityOld = m_porosityOld[er][esr];
    arrayView2d< real64 const > const & phaseVolFracOld = m_phaseVolFracOld[er][esr];
    arrayView2d< real64 const > const & phaseDensOld = m_phaseDensOld[er][esr];
    arrayView3d< real64 const > const & phaseCompFracOld = m_phaseCompFracOld[er][esr];

    arrayView2d< real64 const > const & pvMult = m_pvMult[er][esr];
    arrayView2d< real64 const > const & dPvMult_dPres = m_dPvMult_dPres[er][esr];
    arrayView3d< real64 const > const & phaseDens = m_phaseDens[er][esr];
    arrayView3d< real64 const > const & dPhaseDens_dPres = m_dPhaseDens_dPres[er][esr];
    arrayView4d< real64 const > const & dPhaseDens_dComp = m_dPhaseDens_dComp[er][esr];
    arrayView4d< real64 const > const & phaseCompFrac = m_phaseCompFrac[er][esr];
    arrayView4d< real64 const > const & dPhaseCompFrac_dPres = m_dPhaseCompFrac_dPres[er][esr];
    arrayView5d< real64 const > const & dPhaseCompFrac_dComp = m_dPhaseCompFrac_dComp[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        stackArray1d< globalIndex, maxNumDof > localAccumDOF( NDOF );
        stackArray1d< real64, maxNumComp > localAccum( NC );
        stackArray2d< real64, maxNumComp * maxNumDof > localAccumJacobian( NC, NDOF );

        AccumulationKernel::Compute( NC, NP,
                                     volume[ei],
                                     porosityOld[ei],
                                     porosityRef[ei],
                                     pvMult[ei][0],
                                     dPvMult_dPres[ei][0],
                                     dCompFrac_dCompDens[ei],
                                     phaseVolFracOld[ei],
                                     phaseVolFrac[ei],
                                     dPhaseVolFrac_dPres[ei],
                                     dPhaseVolFrac_dCompDens[ei],
                                     phaseDensOld[ei],
                                     phaseDens[ei][0],
                                     dPhaseDens_dPres[ei][0],
                                     dPhaseDens_dComp[ei][0],
                                     phaseCompFracOld[ei],
                                     phaseCompFrac[ei][0],
                                     dPhaseCompFrac_dPres[ei][0],
                                     dPhaseCompFrac_dComp[ei][0],
                                     localAccum,
                                     localAccumJacobian );

        // set DOF indices for this block
        for( localIndex idof = 0; idof < NDOF; ++idof )
        {
          localAccumDOF[idof] = dofNumber[ei] + idof;
        }

        // TODO: apply equation/variable change transformation(s)

        // add contribution to global residual and dRdP
        rhs->add( localAccumDOF.data(),
                  localAccum.data(),
                  NC );

        matrix->add( localAccumDOF.data(),
                     localAccumDOF.data(),
                     localAccumJacobian.data(),
                     NC, NDOF );
      }
    } );
  } );
}

void CompositionalMultiphaseFlow::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                     real64 const dt,
                                                     DomainPartition const * const domain,
                                                     DofManager const * const dofManager,
                                                     ParallelMatrix * const matrix,
                                                     ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup< NumericalMethodsManager >( keys::numericalMethodsManager );

  FiniteVolumeManager const * const fvManager =
    numericalMethodManager->GetGroup< FiniteVolumeManager >( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager->getKey( viewKeyStruct::dofFieldString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = elemManager->ConstructViewAccessor< array1d< globalIndex >,
                                                          arrayView1d< globalIndex const > >( dofKey );

  FluxKernel::ElementView< arrayView1d< globalIndex const > > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementView< arrayView1d< real64 const > > const & pres                = m_pressure.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & dPres               = m_deltaPressure.toViewConst();
  FluxKernel::ElementView< arrayView1d< real64 const > > const & gravCoef            = m_gravCoef.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & phaseMob            = m_phaseMob.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres     = m_dPhaseMob_dPres.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dPhaseMob_dComp     = m_dPhaseMob_dCompDens.toViewConst();
  FluxKernel::ElementView< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp = m_dPhaseVolFrac_dCompDens.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens = m_dCompFrac_dCompDens.toViewConst();

  FluxKernel::ElementView< arrayView3d< real64 const > > const & phaseDens                   = m_phaseDens.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres            = m_dPhaseDens_dPres.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & dPhaseDens_dComp            = m_dPhaseDens_dComp.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & phaseCompFrac               = m_phaseCompFrac.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres        = m_dPhaseCompFrac_dPres.toViewConst();
  FluxKernel::ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp        = m_dPhaseCompFrac_dComp.toViewConst();
  FluxKernel::ElementView< arrayView3d< real64 const > > const & phaseCapPres                = m_phaseCapPressure.toViewConst();
  FluxKernel::ElementView< arrayView4d< real64 const > > const & dPhaseCapPres_dPhaseVolFrac = m_dPhaseCapPressure_dPhaseVolFrac.toViewConst();

  localIndex constexpr numElems   = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex constexpr maxSize1 = numElems * maxNumComp;
  localIndex constexpr maxSize2 = maxStencil * maxNumDof;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  integer const capPressureFlag = m_capPressureFlag;

  fluxApprox->forAllStencils( [&] ( auto const & stencil )
  {
    typedef TYPEOFREF( stencil ) STENCIL_TYPE;
    typename STENCIL_TYPE::IndexContainerViewConstType const & eri = stencil.getElementRegionIndices();
    typename STENCIL_TYPE::IndexContainerViewConstType const & esri = stencil.getElementSubRegionIndices();
    typename STENCIL_TYPE::IndexContainerViewConstType const & ei = stencil.getElementIndices();
    typename STENCIL_TYPE::WeightContainerViewConstType const & weights = stencil.getWeights();

    forAll< serialPolicy >( stencil.size(), [=] ( localIndex iconn )
    {
      localIndex const stencilSize = stencil.stencilSize( iconn );

      // create local work arrays
      stackArray1d< globalIndex, maxSize1 > eqnRowIndices( numElems * NC );
      stackArray1d< globalIndex, maxSize2 > dofColIndices( stencilSize * NDOF );

      stackArray1d< real64, maxSize1 >            localFlux( numElems * NC );
      stackArray2d< real64, maxSize1 * maxSize2 > localFluxJacobian( numElems * NC, stencilSize * NDOF );

      FluxKernel::Compute( NC, NP,
                           stencilSize,
                           eri[iconn],
                           esri[iconn],
                           ei[iconn],
                           weights[iconn],
                           pres,
                           dPres,
                           gravCoef,
                           phaseMob,
                           dPhaseMob_dPres,
                           dPhaseMob_dComp,
                           dPhaseVolFrac_dPres,
                           dPhaseVolFrac_dComp,
                           dCompFrac_dCompDens,
                           phaseDens,
                           dPhaseDens_dPres,
                           dPhaseDens_dComp,
                           phaseCompFrac,
                           dPhaseCompFrac_dPres,
                           dPhaseCompFrac_dComp,
                           phaseCapPres,
                           dPhaseCapPres_dPhaseVolFrac,
                           capPressureFlag,
                           dt,
                           localFlux,
                           localFluxJacobian );

      // set equation indices for both connected cells
      for( localIndex i = 0; i < numElems; ++i )
      {
        globalIndex const offset = dofNumber[eri( iconn, i )][esri( iconn, i )][ei( iconn, i )];

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          eqnRowIndices[i * NC + ic] = offset + ic;
        }
      }

      for( localIndex i = 0; i < stencilSize; ++i )
      {
        globalIndex const offset = dofNumber[eri( iconn, i )][esri( iconn, i )][ei( iconn, i )];

        for( localIndex jdof = 0; jdof < NDOF; ++jdof )
        {
          dofColIndices[i * NDOF + jdof] = offset + jdof;
        }
      }

      // TODO: apply equation/variable change transformation(s)

      // Add to global residual/jacobian
      rhs->add( eqnRowIndices.data(),
                localFlux.data(),
                numElems * NC );

      matrix->add( eqnRowIndices.data(),
                   dofColIndices.data(),
                   localFluxJacobian.data(),
                   numElems * NC,
                   stencilSize * NDOF );

    } );
  } );
}

void CompositionalMultiphaseFlow::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const GEOSX_UNUSED_PARAM( dt ),
                                                              DomainPartition const * const domain,
                                                              DofManager const * const dofManager,
                                                              ParallelMatrix * const matrix,
                                                              ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  string const dofKey = dofManager->getKey( viewKeyStruct::dofFieldString );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d< real64 const > const & volume = m_volume[er][esr];
    arrayView1d< real64 const > const & porosityRef = m_porosityRef[er][esr];

    arrayView2d< real64 const > const & phaseVolFrac = m_phaseVolFrac[er][esr];
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];

    arrayView2d< real64 const > const & pvMult = m_pvMult[er][esr];
    arrayView2d< real64 const > const & dPvMult_dPres = m_dPvMult_dPres[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      if( elemGhostRank[ei] < 0 )
      {

        real64 localVolBalance;
        stackArray1d< real64, maxNumDof > localVolBalanceJacobian( NDOF );
        stackArray1d< globalIndex, maxNumDof > localVolBalanceDOF( NDOF );

        VolumeBalanceKernel::Compute( NC, NP,
                                      volume[ei],
                                      porosityRef[ei],
                                      pvMult[ei][0],
                                      dPvMult_dPres[ei][0],
                                      phaseVolFrac[ei],
                                      dPhaseVolFrac_dPres[ei],
                                      dPhaseVolFrac_dCompDens[ei],
                                      localVolBalance,
                                      localVolBalanceJacobian );

        // get equation/dof indices
        globalIndex const localVolBalanceEqnIndex = dofNumber[ei] + NC;
        for( localIndex jdof = 0; jdof < NDOF; ++jdof )
        {
          localVolBalanceDOF[jdof] = dofNumber[ei] + jdof;
        }

        // TODO: apply equation/variable change transformation(s)

        // add contribution to global residual and dRdP
        rhs->add( localVolBalanceEqnIndex,
                  localVolBalance );

        matrix->add( localVolBalanceEqnIndex,
                     localVolBalanceDOF.data(),
                     localVolBalanceJacobian.data(),
                     NDOF );
      }
    } );
  } );
}

void CompositionalMultiphaseFlow::ApplyBoundaryConditions( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition * const domain,
                                                           DofManager const & dofManager,
                                                           ParallelMatrix & matrix,
                                                           ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.open();
  rhs.open();

  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit( time_n, dt, &dofManager, domain, &matrix, &rhs );

  // apply flux boundary conditions
  ApplySourceFluxBC( time_n, dt, &dofManager, domain, &matrix, &rhs );


  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseFlow::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, LAIOutputFormat::MATRIX_MARKET );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, LAIOutputFormat::MATRIX_MARKET );

    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseFlow::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }


}

void
CompositionalMultiphaseFlow::ApplySourceFluxBC( real64 const time,
                                                real64 const dt,
                                                DofManager const * const dofManager,
                                                DomainPartition * const domain,
                                                ParallelMatrix * const matrix,
                                                ParallelVector * const rhs )
{

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  string const dofKey = dofManager->getKey( viewKeyStruct::dofFieldString );

  fsManager.Apply( time + dt, domain, "ElementRegions", FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & ) -> void
  {

    arrayView1d< globalIndex const > const &
    dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const &
    ghostRank = subRegion->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    SortedArray< localIndex > localSet;
    for( localIndex const a : lset )
    {
      if( ghostRank[a] < 0 )
      {
        localSet.insert( a );
      }
    }

    fs->ApplyBoundaryConditionToSystem< FieldSpecificationAdd, LAInterface >( localSet,
                                                                              time + dt,
                                                                              dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              integer_conversion< int >( m_numDofPerCell ),
                                                                              *matrix,
                                                                              *rhs,
                                                                              [&] ( localIndex const GEOSX_UNUSED_PARAM( a )) -> real64
    {
      return 0;
    } );

  } );
}


void
CompositionalMultiphaseFlow::ApplyDirichletBC_implicit( real64 const time,
                                                        real64 const dt,
                                                        DofManager const * const dofManager,
                                                        DomainPartition * const domain,
                                                        ParallelMatrix * const matrix,
                                                        ParallelVector * const rhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

  // 1. apply pressure Dirichlet BCs
  fsManager.Apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::pressureString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const & setName,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * subRegion,
                        string const & )
  {
    // 1.0. Check whether pressure has already been applied to this set
    string const & subRegionName = subRegion->getName();
    GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0, "Conflicting pressure boundary conditions on set " << setName );
    bcStatusMap[subRegionName][setName].resize( m_numComponents );
    bcStatusMap[subRegionName][setName] = false;

    // 1.1. Apply BC to set the field values
    fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                    time + dt,
                                                    subRegion,
                                                    viewKeyStruct::bcPressureString );
  } );

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  fsManager.Apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::globalCompFractionString,
                   [&] ( FieldSpecificationBase const * const fs,
                         string const & setName,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group * subRegion,
                         string const & )
  {
    // 2.0. Check pressure and record composition bc application
    string const & subRegionName = subRegion->getName();
    localIndex const comp = fs->GetComponent();
    GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0, "Pressure boundary condition not prescribed on set '" << setName << "'" );
    GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
    bcStatusMap[subRegionName][setName][comp] = true;

    // 2.1. Apply BC to set the field values
    fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                    time + dt,
                                                    subRegion,
                                                    viewKeyStruct::globalCompFractionString );
  } );

  // 2.3 Check consistency between composition BC applied to sets
  bool bcConsistent = true;
  for( auto const & bcStatusEntryOuter : bcStatusMap )
  {
    for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
    {
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        bcConsistent &= bcStatusEntryInner.second[ic];
        GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic <<
                          " on region " << bcStatusEntryOuter.first << ", set " << bcStatusEntryInner.first );
      }
    }
  }
  GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

  string const dofKey = dofManager->getKey( viewKeyStruct::dofFieldString );

  // 3. Call constitutive update, back-calculate target global component densities and apply to the system
  fsManager.Apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::pressureString,
                   [&] ( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( fs ),
                         string const & GEOSX_UNUSED_PARAM( setName ),
                         SortedArrayView< localIndex const > const & targetSet,
                         Group * subRegion,
                         string const & )
  {
    // TODO: hack! Find a better way to get the fluid
    Group const * const region = subRegion->getParent()->getParent();
    string const & fluidName = m_fluidModelNames[ targetRegionIndex( region->getName() ) ];
    MultiFluidBase & fluid = GetConstitutiveModel< MultiFluidBase >( *subRegion, fluidName );

    arrayView1d< globalIndex const > const & dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const & pres      = subRegion->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 const > const & dPres     = subRegion->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 const > const & bcPres    = subRegion->getReference< array1d< real64 > >( viewKeyStruct::bcPressureString );
    arrayView2d< real64 const > const & compFrac  = subRegion->getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );
    arrayView2d< real64 const > const & compDens  = subRegion->getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 const > const & dCompDens = subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );
    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    array1d< real64 > rhsContribution( targetSet.size() * m_numDofPerCell );
    array1d< globalIndex > dof( targetSet.size() * m_numDofPerCell );

    integer counter = 0;

    for( localIndex a : targetSet )
    {
      fluid.PointUpdate( bcPres[a], m_temperature, compFrac[a], a, 0 );

      dof[counter] = dofNumber[a];

      // 4.1. Apply pressure to the matrix
      FieldSpecificationEqual::SpecifyFieldValue< LAInterface >( dof[counter],
                                                                 *matrix,
                                                                 rhsContribution[counter],
                                                                 bcPres[a],
                                                                 pres[a] + dPres[a] );

      ++counter;

      // 4.2. For each component, apply target global density value
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        dof[counter] = dofNumber[a] + ic + 1;
        real64 const targetCompDens = totalDens[a][0] * compFrac[a][ic];

        FieldSpecificationEqual::SpecifyFieldValue< LAInterface >( dof[counter],
                                                                   *matrix,
                                                                   rhsContribution[counter],
                                                                   targetCompDens,
                                                                   compDens[a][ic] + dCompDens[a][ic] );

        ++counter;
      }

      // 4.3. Apply accumulated rhs values
      FieldSpecificationEqual::PrescribeRhsValues< LAInterface >( *rhs,
                                                                  counter,
                                                                  dof.data(),
                                                                  rhsContribution.data() );
    }
  } );

}

real64
CompositionalMultiphaseFlow::CalculateResidualNorm( DomainPartition const * const domain,
                                                    DofManager const & dofManager,
                                                    ParallelVector const & rhs )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  real64 const * localResidual = rhs.extractLocalVector();
  real64 localResidualNorm = 0.0;

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase const & subRegion )
  {
    MultiFluidBase const & fluid = GetConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

    arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & refPoro =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );
    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        globalIndex const offset = dofNumber[ei];
        for( localIndex idof = 0; idof < m_numDofPerCell; ++idof )
        {
          localIndex const lid = rhs.getLocalRowID( offset + idof );
          real64 const val = localResidual[lid] / ( totalDens[ei][0] * refPoro[ei] * volume[ei] );
          localResidualNorm += val * val;
        }
      }
    }
  } );

  // compute global residual norm
  realT globalResidualNorm;
  MpiWrapper::allReduce( &localResidualNorm, &globalResidualNorm, 1, MPI_SUM, MPI_COMM_GEOSX );

  return sqrt( globalResidualNorm );
}

void CompositionalMultiphaseFlow::SolveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

bool
CompositionalMultiphaseFlow::CheckSystemSolution( DomainPartition const * const domain,
                                                  DofManager const & dofManager,
                                                  ParallelVector const & solution,
                                                  real64 const scalingFactor )
{
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  real64 const * localSolution = solution.extractLocalVector();
  int localCheck = 1;

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase const &,
                                    ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d< real64 const > const & pres = m_pressure[er][esr];
    arrayView1d< real64 const > const & dPres = m_deltaPressure[er][esr];
    arrayView2d< real64 const > const & compDens = m_globalCompDensity[er][esr];
    arrayView2d< real64 const > const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        globalIndex const offset = dofNumber[ei];
        // extract solution and apply to dP
        {
          localIndex const lid = solution.getLocalRowID( offset );
          real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[lid];

          if( newPres < 0.0 )
          {
            localCheck = 0;
          }
        }

        for( localIndex ic = 0; ic < m_numComponents; ++ic )
        {
          localIndex const lid = solution.getLocalRowID( offset + ic + 1 );
          real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[lid];
          if( newDens < 0.0 )
          {
            localCheck = 0;
          }
        }
      }
    } );
  } );

  return MpiWrapper::Min( localCheck, MPI_COMM_GEOSX );
}

void
CompositionalMultiphaseFlow::ApplySystemSolution( DofManager const & dofManager,
                                                  ParallelVector const & solution,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.addVectorToField( solution,
                               viewKeyStruct::dofFieldString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( solution,
                               viewKeyStruct::dofFieldString,
                               viewKeyStruct::deltaGlobalCompDensityString,
                               scalingFactor,
                               1, m_numDofPerCell );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaGlobalCompDensityString );
  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain->getNeighbors() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateState( subRegion, targetIndex );
  } );
}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase &,
                                    ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dPres = m_deltaPressure[er][esr];
    arrayView2d< real64 > const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      dPres[ei] = 0.0;
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        dCompDens[ei][ic] = 0.0;
      }
    } );
  } );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    UpdateState( subRegion, targetIndex );
  } );
}

void CompositionalMultiphaseFlow::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase &,
                                    ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const & dPres = m_deltaPressure[er][esr];
    arrayView2d< real64 const > const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    arrayView1d< real64 > const & pres = m_pressure[er][esr];
    arrayView2d< real64 > const & compDens = m_globalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        compDens[ei][ic] += dCompDens[ei][ic];
      }
    } );
  } );
}

void CompositionalMultiphaseFlow::ResetViews( DomainPartition * const domain )
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_pressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaPressureString );
  m_globalCompDensity =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::globalCompDensityString );
  m_deltaGlobalCompDensity =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  m_compFrac =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::globalCompFractionString );
  m_dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor< array3d< real64 >, arrayView3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
  m_phaseVolFrac =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );
  m_dPhaseVolFrac_dPres =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  m_dPhaseVolFrac_dCompDens =
    elemManager->ConstructViewAccessor< array3d< real64 >, arrayView3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );
  m_phaseMob =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::phaseMobilityString );
  m_dPhaseMob_dPres =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString );
  m_dPhaseMob_dCompDens =
    elemManager->ConstructViewAccessor< array3d< real64 >, arrayView3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  m_porosityOld =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::porosityOldString );
  m_phaseVolFracOld =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::phaseVolumeFractionOldString );
  m_phaseDensOld =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::phaseDensityOldString );
  m_phaseCompFracOld =
    elemManager->ConstructViewAccessor< array3d< real64 >, arrayView3d< real64 > >( viewKeyStruct::phaseComponentFractionOldString );

  {
    using keys = ConstitutiveBase::viewKeyStruct;

    m_pvMult =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::poreVolumeMultiplierString,
                                                                                              targetRegionNames(),
                                                                                              solidModelNames() );
    m_dPvMult_dPres =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::dPVMult_dPresString,
                                                                                              targetRegionNames(),
                                                                                              solidModelNames() );
  }
  {
    using keys = MultiFluidBase::viewKeyStruct;

    m_phaseFrac =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseFractionString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseFrac_dPres =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dPhaseFraction_dPressureString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseFrac_dComp =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseFraction_dGlobalCompFractionString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_phaseDens =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseDensityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseDens_dPres =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dPhaseDensity_dPressureString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseDens_dComp =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseDensity_dGlobalCompFractionString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_phaseVisc =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseViscosityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseVisc_dPres =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::dPhaseViscosity_dPressureString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseVisc_dComp =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseViscosity_dGlobalCompFractionString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_phaseCompFrac =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::phaseCompFractionString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseCompFrac_dPres =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseCompFraction_dPressureString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_dPhaseCompFrac_dComp =
      elemManager->ConstructMaterialViewAccessor< array5d< real64 >, arrayView5d< real64 > >( keys::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
    m_totalDens =
      elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( keys::totalDensityString,
                                                                                              targetRegionNames(),
                                                                                              fluidModelNames() );
  }
  {
    using keys = RelativePermeabilityBase::viewKeyStruct;

    m_phaseRelPerm =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseRelPermString,
                                                                                              targetRegionNames(),
                                                                                              relPermModelNames() );
    m_dPhaseRelPerm_dPhaseVolFrac =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                              targetRegionNames(),
                                                                                              relPermModelNames() );
  }
  if( m_capPressureFlag )
  {
    using keys = CapillaryPressureBase::viewKeyStruct;

    m_phaseCapPressure =
      elemManager->ConstructMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( keys::phaseCapPressureString,
                                                                                              targetRegionNames(),
                                                                                              m_capPressureModelNames );
    m_dPhaseCapPressure_dPhaseVolFrac =
      elemManager->ConstructMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( keys::dPhaseCapPressure_dPhaseVolFractionString,
                                                                                              targetRegionNames(),
                                                                                              m_capPressureModelNames );
  }
}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseFlow, string const &, Group * const )
}// namespace geosx
