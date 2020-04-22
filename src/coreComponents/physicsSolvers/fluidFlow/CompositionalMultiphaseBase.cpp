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
 * @file CompositionalMultiphaseBase.cpp
 */

#include "CompositionalMultiphaseBase.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
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
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseBaseKernels;

CompositionalMultiphaseBase::CompositionalMultiphaseBase( const string & name,
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
  this->registerWrapper( viewKeyStruct::temperatureString, &m_temperature, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Temperature" );

  this->registerWrapper( viewKeyStruct::useMassFlagString, &m_useMass, false )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::relPermNameString, &m_relPermName, false )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the relative permeability constitutive model to use" );

  this->registerWrapper( viewKeyStruct::relPermIndexString, &m_relPermIndex, false );

  this->registerWrapper( viewKeyStruct::capPressureNameString, &m_capPressureName, false )->
    setApplyDefaultValue( "" )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the capillary pressure constitutive model to use" );

  this->registerWrapper( viewKeyStruct::capPressureIndexString, &m_capPressureIndex, false );
}

localIndex CompositionalMultiphaseBase::numFluidComponents() const
{
  return m_numComponents;
}

localIndex CompositionalMultiphaseBase::numFluidPhases() const
{
  return m_numPhases;
}

void CompositionalMultiphaseBase::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();

  if( !m_capPressureName.empty())
  {
    m_capPressureFlag = 1;
  }
  else
  {
    m_capPressureIndex = -1;
  }
}

void CompositionalMultiphaseBase::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase & elementSubRegion )
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

void CompositionalMultiphaseBase::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  MultiFluidBase const * fluid = cm->GetConstitutiveRelation< MultiFluidBase >( m_fluidName );
  m_numPhases     = fluid->numFluidPhases();
  m_numComponents = fluid->numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  RelativePermeabilityBase const * relPerm = cm->GetConstitutiveRelation< RelativePermeabilityBase >( m_relPermName );
  GEOSX_ERROR_IF( relPerm == nullptr, "Relative permeability model " + m_relPermName + " not found" );
  m_relPermIndex = relPerm->getIndexInParent();

  CapillaryPressureBase const * capPressure = cm->GetConstitutiveRelation< CapillaryPressureBase >( m_capPressureName );
  if( m_capPressureFlag )
  {
    GEOSX_ERROR_IF( capPressure == nullptr, "Capillary pressure model " + m_capPressureName + " not found" );
    m_capPressureIndex = capPressure->getIndexInParent();
  }

  // Consistency check between the models
  GEOSX_ERROR_IF( fluid->numFluidPhases() != relPerm->numFluidPhases(),
                  "Number of fluid phases differs between fluid model '" << m_fluidName
                                                                         << "' and relperm model '" << m_relPermName << "'" );
  if( m_capPressureFlag )
  {
    GEOSX_ERROR_IF( fluid->numFluidPhases() != capPressure->numFluidPhases(),
                    "Number of fluid phases differs between fluid model '" << m_fluidName
                                                                           << "' and capillary pressure model '" << m_capPressureName << "'" );
  }

  for( localIndex ip = 0; ip < m_numPhases; ++ip )
  {
    string const & phase_fl = fluid->phaseName( ip );
    string const & phase_rp = relPerm->phaseName( ip );
    GEOSX_ERROR_IF( phase_fl != phase_rp, "Phase '" << phase_fl << "' in fluid model '" << m_fluidName
                                                    << "' does not match phase '" << phase_rp << "' in relperm model '" << m_relPermName << "'" );

    if( m_capPressureFlag )
    {
      string const & phase_pc = capPressure->phaseName( ip );
      GEOSX_ERROR_IF( phase_fl != phase_pc, "Phase '" << phase_fl << "' in fluid model '" << m_fluidName
                                                      << "' does not match phase '" << phase_pc << "' in cap pressure model '" << m_capPressureName << "'" );
    }
  }

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ResizeFields( meshLevel );
  }
}

void CompositionalMultiphaseBase::ResizeFields( MeshLevel * const meshLevel )
{
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase & subRegion )
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

void CompositionalMultiphaseBase::UpdateComponentFraction( Group * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & compFrac =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );

  arrayView3d< real64 > const & dCompFrac_dCompDens =
    dataGroup->getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const & compDens =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );

  arrayView2d< real64 const > const & dCompDens =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  KernelLaunchSelector1< ComponentFractionKernel >( m_numComponents,
                                                    dataGroup->size(),
                                                    compDens,
                                                    dCompDens,
                                                    compFrac,
                                                    dCompFrac_dCompDens );
}

void CompositionalMultiphaseBase::UpdateComponentFraction( localIndex er, localIndex esr ) const
{
  GEOSX_MARK_FUNCTION;

  KernelLaunchSelector1< ComponentFractionKernel >( m_numComponents,
                                                    m_compFrac[er][esr].size( 0 ),
                                                    m_globalCompDensity[er][esr],
                                                    m_deltaGlobalCompDensity[er][esr],
                                                    m_compFrac[er][esr],
                                                    m_dCompFrac_dCompDens[er][esr] );
}

void CompositionalMultiphaseBase::UpdatePhaseVolumeFraction( Group * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseVolFrac =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );

  arrayView2d< real64 > const & dPhaseVolFrac_dPres =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 > const & dPhaseVolFrac_dComp =
    dataGroup->getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  // inputs

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    dataGroup->getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  arrayView2d< real64 const > const & compDens =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString );

  arrayView2d< real64 const > const & dCompDens =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString );

  MultiFluidBase * fluid = GetConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidName );

  arrayView3d< real64 const > const & phaseFrac =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::phaseFractionString );

  arrayView3d< real64 const > const & dPhaseFrac_dPres =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

  arrayView4d< real64 const > const & dPhaseFrac_dComp =
    fluid->getReference< array4d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

  arrayView3d< real64 const > const & phaseDens =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d< real64 const > const & dPhaseDens_dPres =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d< real64 const > const & dPhaseDens_dComp =
    fluid->getReference< array4d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  KernelLaunchSelector2< PhaseVolumeFractionKernel >( m_numComponents, m_numPhases,
                                                      dataGroup->size(),
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

void CompositionalMultiphaseBase::UpdatePhaseVolumeFraction( localIndex er, localIndex esr ) const
{
  GEOSX_MARK_FUNCTION;

  KernelLaunchSelector2< PhaseVolumeFractionKernel >( m_numComponents, m_numPhases,
                                                      m_phaseVolFrac[er][esr].size( 0 ),
                                                      m_globalCompDensity[er][esr],
                                                      m_deltaGlobalCompDensity[er][esr],
                                                      m_dCompFrac_dCompDens[er][esr],
                                                      m_phaseDens[er][esr][m_fluidIndex],
                                                      m_dPhaseDens_dPres[er][esr][m_fluidIndex],
                                                      m_dPhaseDens_dComp[er][esr][m_fluidIndex],
                                                      m_phaseFrac[er][esr][m_fluidIndex],
                                                      m_dPhaseFrac_dPres[er][esr][m_fluidIndex],
                                                      m_dPhaseFrac_dComp[er][esr][m_fluidIndex],
                                                      m_phaseVolFrac[er][esr],
                                                      m_dPhaseVolFrac_dPres[er][esr],
                                                      m_dPhaseVolFrac_dCompDens[er][esr] );
}

void CompositionalMultiphaseBase::UpdatePhaseMobility( Group * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseMob =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString );

  arrayView2d< real64 > const & dPhaseMob_dPres =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString );

  arrayView3d< real64 > const & dPhaseMob_dComp =
    dataGroup->getReference< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 const > const & dPhaseVolFrac_dComp =
    dataGroup->getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    dataGroup->getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  MultiFluidBase const * fluid = GetConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidName );

  arrayView3d< real64 const > const & phaseDens =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d< real64 const > const & dPhaseDens_dPres =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d< real64 const > const & dPhaseDens_dComp =
    fluid->getReference< array4d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  arrayView3d< real64 const > const & phaseVisc =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::phaseViscosityString );

  arrayView3d< real64 const > const & dPhaseVisc_dPres =
    fluid->getReference< array3d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

  arrayView4d< real64 const > const & dPhaseVisc_dComp =
    fluid->getReference< array4d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString );

  RelativePermeabilityBase const * relperm = GetConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermName );

  arrayView3d< real64 const > const & phaseRelPerm =
    relperm->getReference< array3d< real64 > >( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString );

  arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac =
    relperm->getReference< array4d< real64 > >( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString );

  KernelLaunchSelector2< PhaseMobilityKernel >( m_numComponents, m_numPhases,
                                                dataGroup->size(),
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

void CompositionalMultiphaseBase::UpdatePhaseMobility( localIndex er, localIndex esr ) const
{
  GEOSX_MARK_FUNCTION;

  KernelLaunchSelector2< PhaseMobilityKernel >( m_numComponents, m_numPhases,
                                                m_phaseMob[er][esr].size( 0 ),
                                                m_dCompFrac_dCompDens[er][esr],
                                                m_phaseDens[er][esr][m_fluidIndex],
                                                m_dPhaseDens_dPres[er][esr][m_fluidIndex],
                                                m_dPhaseDens_dComp[er][esr][m_fluidIndex],
                                                m_phaseVisc[er][esr][m_fluidIndex],
                                                m_dPhaseVisc_dPres[er][esr][m_fluidIndex],
                                                m_dPhaseVisc_dComp[er][esr][m_fluidIndex],
                                                m_phaseRelPerm[er][esr][m_relPermIndex],
                                                m_dPhaseRelPerm_dPhaseVolFrac[er][esr][m_relPermIndex],
                                                m_dPhaseVolFrac_dPres[er][esr],
                                                m_dPhaseVolFrac_dCompDens[er][esr],
                                                m_phaseMob[er][esr],
                                                m_dPhaseMob_dPres[er][esr],
                                                m_dPhaseMob_dCompDens[er][esr] );
}

void CompositionalMultiphaseBase::UpdateFluidModel( Group * const dataGroup )
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase * const fluid = GetConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidName );

  arrayView1d< real64 const > const & pres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
  arrayView2d< real64 const > const & compFrac = dataGroup->getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString );

  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forAll< serialPolicy >( dataGroup->size(), [=] ( localIndex const a )
  {
    fluid->PointUpdate( pres[a] + dPres[a], m_temperature, compFrac[a], a, 0 );
  } );
  //fluid->BatchUpdate( pres, temp, compFrac );
}

void CompositionalMultiphaseBase::UpdateSolidModel( Group * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase * const solid = GetConstitutiveModel< ConstitutiveBase >( dataGroup, m_solidName );

  arrayView1d< real64 const > const & pres  = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  forAll< serialPolicy >( dataGroup->size(), [=] ( localIndex const a )
  {
    solid->StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  } );
}

void CompositionalMultiphaseBase::UpdateRelPermModel( Group * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  RelativePermeabilityBase * const relPerm = GetConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermName );

  arrayView2d< real64 > const & phaseVolFrac =
    dataGroup->getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );

  relPerm->BatchUpdate( phaseVolFrac );
}

void CompositionalMultiphaseBase::UpdateCapPressureModel( Group * dataGroup )
{
  if( m_capPressureFlag )
  {
    CapillaryPressureBase * const capPressure = GetConstitutiveModel< CapillaryPressureBase >( dataGroup, m_capPressureName );

    arrayView2d< real64 > const & phaseVolFrac =
      dataGroup->getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString );

    capPressure->BatchUpdate( phaseVolFrac );
  }
}

void CompositionalMultiphaseBase::UpdateState( Group * const dataGroup )
{
  GEOSX_MARK_FUNCTION;

  UpdateComponentFraction( dataGroup );
  UpdateFluidModel( dataGroup );
  UpdatePhaseVolumeFraction( dataGroup );
  UpdateSolidModel( dataGroup );
  UpdateRelPermModel( dataGroup );
  UpdatePhaseMobility( dataGroup );
  UpdateCapPressureModel( dataGroup );
}

void CompositionalMultiphaseBase::InitializeFluidState( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  applyToSubRegionsComplete( mesh,
                             [&] ( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    // 1. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    UpdateFluidModel( &subRegion );

    // 2. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    arrayView2d< real64 const > const & compFrac  = m_compFrac[er][esr];
    arrayView2d< real64 const > const & totalDens = m_totalDens[er][esr][m_fluidIndex];
    arrayView2d< real64 > compDens = m_globalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex ei )
    {
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    } );

    // 3. Calculate phase saturations
    UpdatePhaseVolumeFraction( &subRegion );

    // 4. Initialize solid state
    UpdateSolidModel( &subRegion );

    // 5. Initialize rel perm state
    UpdateRelPermModel( &subRegion );

    // 6. Initialize mobility
    UpdatePhaseMobility( &subRegion );

    // 7. Initialize cap pressure state
    UpdateCapPressureModel( &subRegion );

  } );
}

void CompositionalMultiphaseBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::globalCompDensityString );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  // set mass fraction flag on main model
  // TODO find a way to set this before constitutive model is duplicated and attached to subregions?
  {
    MultiFluidBase * const fluid = constitutiveManager->GetConstitutiveRelation< MultiFluidBase >( m_fluidName );
    fluid->setMassFlag( static_cast< bool >(m_useMass) );
  }

  // set mass fraction flag on subregion models
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase & subRegion )
  {
    MultiFluidBase * const fluid = GetConstitutiveModel< MultiFluidBase >( &subRegion, m_fluidName );
    fluid->setMassFlag( static_cast< bool >(m_useMass) );
  } );

  // Initialize primary variables from applied initial conditions
  ResetViews( domain );
  InitializeFluidState( domain );
}

real64 CompositionalMultiphaseBase::SolverStep( real64 const & time_n,
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

void CompositionalMultiphaseBase::BackupFields( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  // backup some fields used in time derivative approximation
  applyToSubRegionsComplete( mesh,
                             [&] ( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d< real64 const > const & poroRef       = m_porosityRef[er][esr];

    arrayView2d< real64 const > const & phaseVolFrac = m_phaseVolFrac[er][esr];

    arrayView3d< real64 const > const & phaseDens     = m_phaseDens[er][esr][m_fluidIndex];
    arrayView4d< real64 const > const & phaseCompFrac = m_phaseCompFrac[er][esr][m_fluidIndex];
    arrayView2d< real64 const > const & pvMult        = m_pvMult[er][esr][m_solidIndex];

    arrayView2d< real64 > const & phaseDensOld     = m_phaseDensOld[er][esr];
    arrayView2d< real64 > const & phaseVolFracOld  = m_phaseVolFracOld[er][esr];
    arrayView3d< real64 > const & phaseCompFracOld = m_phaseCompFracOld[er][esr];
    arrayView1d< real64 > const & poroOld          = m_porosityOld[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex ei )
    {
      if( elemGhostRank[ei] >= 0 )
        return;

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
    } );
  } );
}

void
CompositionalMultiphaseBase::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
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


void CompositionalMultiphaseBase::AssembleSystem( real64 const time_n,
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
    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseBase::AssembleSystem" );
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

    GEOSX_LOG_RANK_0( "After CompositionalMultiphaseBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void CompositionalMultiphaseBase::AssembleAccumulationTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                             real64 const GEOSX_UNUSED_PARAM( dt ),
                                                             DomainPartition const * const domain,
                                                             DofManager const * const dofManager,
                                                             ParallelMatrix * const matrix,
                                                             ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  string const dofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );

  applyToSubRegionsComplete( mesh,
                             [&] ( localIndex const er, localIndex const esr, ElementRegionBase const &, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d< real64 const > const & volume      = m_volume[er][esr];
    arrayView1d< real64 const > const & porosityRef = m_porosityRef[er][esr];

    arrayView2d< real64 const > const & phaseVolFrac            = m_phaseVolFrac[er][esr];
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres     = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];
    arrayView3d< real64 const > const & dCompFrac_dCompDens     = m_dCompFrac_dCompDens[er][esr];

    arrayView1d< real64 const > const & porosityOld      = m_porosityOld[er][esr];
    arrayView2d< real64 const > const & phaseVolFracOld  = m_phaseVolFracOld[er][esr];
    arrayView2d< real64 const > const & phaseDensOld     = m_phaseDensOld[er][esr];
    arrayView3d< real64 const > const & phaseCompFracOld = m_phaseCompFracOld[er][esr];

    arrayView2d< real64 const > const & pvMult               = m_pvMult[er][esr][m_solidIndex];
    arrayView2d< real64 const > const & dPvMult_dPres        = m_dPvMult_dPres[er][esr][m_solidIndex];
    arrayView3d< real64 const > const & phaseDens            = m_phaseDens[er][esr][m_fluidIndex];
    arrayView3d< real64 const > const & dPhaseDens_dPres     = m_dPhaseDens_dPres[er][esr][m_fluidIndex];
    arrayView4d< real64 const > const & dPhaseDens_dComp     = m_dPhaseDens_dComp[er][esr][m_fluidIndex];
    arrayView4d< real64 const > const & phaseCompFrac        = m_phaseCompFrac[er][esr][m_fluidIndex];
    arrayView4d< real64 const > const & dPhaseCompFrac_dPres = m_dPhaseCompFrac_dPres[er][esr][m_fluidIndex];
    arrayView5d< real64 const > const & dPhaseCompFrac_dComp = m_dPhaseCompFrac_dComp[er][esr][m_fluidIndex];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        stackArray1d< globalIndex, maxNumDof >         localAccumDOF( NDOF );
        stackArray1d< real64, maxNumComp >             localAccum( NC );
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


void CompositionalMultiphaseBase::AssembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const GEOSX_UNUSED_PARAM( dt ),
                                                              DomainPartition const * const domain,
                                                              DofManager const * const dofManager,
                                                              ParallelMatrix * const matrix,
                                                              ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  string const dofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );

  applyToSubRegionsComplete( mesh,
                             [&] ( localIndex const er, localIndex const esr, ElementRegionBase const &, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const >     const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d< real64 const > const & volume      = m_volume[er][esr];
    arrayView1d< real64 const > const & porosityRef = m_porosityRef[er][esr];

    arrayView2d< real64 const > const & phaseVolFrac            = m_phaseVolFrac[er][esr];
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres     = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];

    arrayView2d< real64 const > const & pvMult        = m_pvMult[er][esr][m_solidIndex];
    arrayView2d< real64 const > const & dPvMult_dPres = m_dPvMult_dPres[er][esr][m_solidIndex];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex ei )
    {
      if( elemGhostRank[ei] >= 0 )
      {
        return;
      }

      real64 localVolBalance;
      stackArray1d< real64, maxNumDof >      localVolBalanceJacobian( NDOF );
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
    } );
  } );
}


void CompositionalMultiphaseBase::SolveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}


void CompositionalMultiphaseBase::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  applyToSubRegionsComplete( mesh,
                             [&] ( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dPres     = m_deltaPressure[er][esr];
    arrayView2d< real64 > const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex ei )
    {
      dPres[ei] = 0.0;
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
        dCompDens[ei][ic] = 0.0;
    } );
  } );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase & subRegion )
  {
    UpdateState( &subRegion );
  } );
}

void CompositionalMultiphaseBase::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  applyToSubRegionsComplete( mesh,
                             [&] ( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const & dPres     = m_deltaPressure[er][esr];
    arrayView2d< real64 const > const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    arrayView1d< real64 > const & pres     = m_pressure[er][esr];
    arrayView2d< real64 > const & compDens = m_globalCompDensity[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
        compDens[ei][ic] += dCompDens[ei][ic];
    } );
  } );
}

void CompositionalMultiphaseBase::ResetViews( DomainPartition * const domain )
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

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

  m_pvMult =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                                constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                                constitutiveManager );
  m_phaseFrac =
    elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( MultiFluidBase::viewKeyStruct::phaseFractionString,
                                                                                                constitutiveManager );
  m_dPhaseFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString,
                                                                                                constitutiveManager );
  m_dPhaseFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >(
      MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString,
      constitutiveManager );
  m_phaseDens =
    elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                                constitutiveManager );
  m_dPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                                constitutiveManager );
  m_dPhaseDens_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >(
      MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString,
      constitutiveManager );
  m_phaseVisc =
    elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( MultiFluidBase::viewKeyStruct::phaseViscosityString,
                                                                                                constitutiveManager );
  m_dPhaseVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString,
                                                                                                constitutiveManager );
  m_dPhaseVisc_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >(
      MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString,
      constitutiveManager );
  m_phaseCompFrac =
    elemManager->ConstructFullMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >( MultiFluidBase::viewKeyStruct::phaseCompFractionString,
                                                                                                constitutiveManager );
  m_dPhaseCompFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >(
      MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString,
      constitutiveManager );
  m_dPhaseCompFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array5d< real64 >, arrayView5d< real64 > >(
      MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString,
      constitutiveManager );
  m_totalDens =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( MultiFluidBase::viewKeyStruct::totalDensityString,
                                                                                                constitutiveManager );
  m_phaseRelPerm =
    elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString,
                                                                                                constitutiveManager );
  m_dPhaseRelPerm_dPhaseVolFrac =
    elemManager->ConstructFullMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >(
      RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
      constitutiveManager );
  if( m_capPressureFlag )
  {
    m_phaseCapPressure =
      elemManager->ConstructFullMaterialViewAccessor< array3d< real64 >, arrayView3d< real64 > >( CapillaryPressureBase::viewKeyStruct::phaseCapPressureString,
                                                                                                  constitutiveManager );
    m_dPhaseCapPressure_dPhaseVolFrac =
      elemManager->ConstructFullMaterialViewAccessor< array4d< real64 >, arrayView4d< real64 > >(
        CapillaryPressureBase::viewKeyStruct::dPhaseCapPressure_dPhaseVolFractionString,
        constitutiveManager );
  }
}

}// namespace geosx
