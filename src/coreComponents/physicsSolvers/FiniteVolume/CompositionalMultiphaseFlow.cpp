/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file CompositionalMultiphaseFlow.cpp
 */

#include "CompositionalMultiphaseFlow.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"
#include "constitutive/CapillaryPressure/CapillaryPressureBase.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

#include "physicsSolvers/FiniteVolume/CompositionalMultiphaseFlowKernels.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;
using namespace CompositionalMultiphaseFlowKernels;

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow( const string & name,
                                                          ManagedGroup * const parent )
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

  this->RegisterViewWrapper( viewKeyStruct::temperatureString, &m_temperature, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Temperature");

  this->RegisterViewWrapper( viewKeyStruct::useMassFlagString, &m_useMass, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Use mass formulation instead of molar");

  this->RegisterViewWrapper( viewKeyStruct::relPermNameString,  &m_relPermName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the relative permeability constitutive model to use");

  this->RegisterViewWrapper( viewKeyStruct::relPermIndexString, &m_relPermIndex, false );

  this->RegisterViewWrapper( viewKeyStruct::capPressureNameString,  &m_capPressureName,  false )->
    setApplyDefaultValue("")->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of the capillary pressure constitutive model to use");

  this->RegisterViewWrapper( viewKeyStruct::capPressureIndexString, &m_capPressureIndex, false );
}

localIndex CompositionalMultiphaseFlow::numFluidComponents() const
{
  return m_numComponents;
}

localIndex CompositionalMultiphaseFlow::numFluidPhases() const
{
  return m_numPhases;
}

void CompositionalMultiphaseFlow::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();

  if (!m_capPressureName.empty())
  {
    m_capPressureFlag = 1;
  }
  else
  {
    m_capPressureIndex = -1;
  }
}

void CompositionalMultiphaseFlow::RegisterDataOnMesh(ManagedGroup * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const elementSubRegion)
    {
      elementSubRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      elementSubRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );
      elementSubRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::bcPressureString );

      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::globalCompDensityString )->setPlotLevel(PlotLevel::LEVEL_0);
      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::deltaGlobalCompDensityString );

      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::globalCompFractionString )->setPlotLevel(PlotLevel::LEVEL_0);
      elementSubRegion->RegisterViewWrapper< array3d<real64> >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::phaseVolumeFractionString )->setPlotLevel(PlotLevel::LEVEL_0);
      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
      elementSubRegion->RegisterViewWrapper< array3d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::phaseMobilityString )->setPlotLevel(PlotLevel::LEVEL_0);
      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::dPhaseMobility_dPressureString );
      elementSubRegion->RegisterViewWrapper< array3d<real64> >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::phaseVolumeFractionOldString );
      elementSubRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::phaseDensityOldString );
      elementSubRegion->RegisterViewWrapper< array3d<real64> >( viewKeyStruct::phaseComponentFractionOldString );
      elementSubRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::porosityOldString );

      elementSubRegion->RegisterViewWrapper< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );
    });
  }
}

void CompositionalMultiphaseFlow::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID(BlockIDs::compositionalBlock, this->getName());

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  MultiFluidBase const * fluid = cm->GetConstitituveRelation<MultiFluidBase>( m_fluidName );
  m_numPhases     = fluid->numFluidPhases();
  m_numComponents = fluid->numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  RelativePermeabilityBase const * relPerm = cm->GetConstitituveRelation<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Relative permeability model " + m_relPermName + " not found" );
  m_relPermIndex = relPerm->getIndexInParent();

  CapillaryPressureBase const * capPressure = cm->GetConstitituveRelation<CapillaryPressureBase>( m_capPressureName );
  if (m_capPressureFlag)
  {
    GEOS_ERROR_IF( capPressure == nullptr, "Capillary pressure model " + m_capPressureName + " not found" );
    m_capPressureIndex = capPressure->getIndexInParent();
  }

  // Consistency check between the models
  GEOS_ERROR_IF( fluid->numFluidPhases() != relPerm->numFluidPhases(),
                 "Number of fluid phases differs between fluid model '" << m_fluidName
                 << "' and relperm model '" << m_relPermName << "'" );
  if (m_capPressureFlag)
  {
    GEOS_ERROR_IF( fluid->numFluidPhases() != capPressure->numFluidPhases(),
                   "Number of fluid phases differs between fluid model '" << m_fluidName
                   << "' and capillary pressure model '" << m_capPressureName << "'" );
  }

  for (localIndex ip = 0; ip < m_numPhases; ++ip)
  {
    string const & phase_fl = fluid->phaseName( ip );
    string const & phase_rp = relPerm->phaseName( ip );
    GEOS_ERROR_IF( phase_fl != phase_rp, "Phase '" << phase_fl << "' in fluid model '" << m_fluidName
                   << "' does not match phase '" << phase_rp << "' in relperm model '" << m_relPermName << "'" );

    if (m_capPressureFlag)
    {
      string const & phase_pc = capPressure->phaseName( ip );
      GEOS_ERROR_IF( phase_fl != phase_pc, "Phase '" << phase_fl << "' in fluid model '" << m_fluidName
                     << "' does not match phase '" << phase_pc << "' in cap pressure model '" << m_capPressureName << "'" );
    }
  }

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ResizeFields( meshLevel );
  }
}

void CompositionalMultiphaseFlow::ResizeFields( MeshLevel * const meshLevel )
{
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const subRegion )
  {
    subRegion->getReference< array2d<real64> >(viewKeyStruct::globalCompDensityString).resizeDimension<1>(NC);
    subRegion->getReference< array2d<real64> >(viewKeyStruct::deltaGlobalCompDensityString).resizeDimension<1>(NC);

    subRegion->getReference< array2d<real64> >(viewKeyStruct::globalCompFractionString).resizeDimension<1>(NC);
    subRegion->getReference< array3d<real64> >(viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString).resizeDimension<1,2>(NC, NC);

    subRegion->getReference< array2d<real64> >(viewKeyStruct::phaseVolumeFractionString).resizeDimension<1>(NP);
    subRegion->getReference< array2d<real64> >(viewKeyStruct::dPhaseVolumeFraction_dPressureString).resizeDimension<1>(NP);
    subRegion->getReference< array3d<real64> >(viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

    subRegion->getReference< array2d<real64> >(viewKeyStruct::phaseMobilityString).resizeDimension<1>(NP);
    subRegion->getReference< array2d<real64> >(viewKeyStruct::dPhaseMobility_dPressureString).resizeDimension<1>(NP);
    subRegion->getReference< array3d<real64> >(viewKeyStruct::dPhaseMobility_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

    subRegion->getReference< array2d<real64> >(viewKeyStruct::phaseVolumeFractionOldString).resizeDimension<1>(NP);
    subRegion->getReference< array2d<real64> >(viewKeyStruct::phaseDensityOldString).resizeDimension<1>(NP);
    subRegion->getReference< array3d<real64> >(viewKeyStruct::phaseComponentFractionOldString).resizeDimension<1,2>(NP, NC);
  });
}

void CompositionalMultiphaseFlow::UpdateComponentFraction( ManagedGroup * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d<real64> const & compFrac =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::globalCompFractionString );

  arrayView3d<real64> const & dCompFrac_dCompDens =
    dataGroup->getReference< array3d<real64> >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // inputs

  arrayView2d<real64 const> const & compDens =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dCompDens =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::deltaGlobalCompDensityString );

  KernelLaunchSelector1<ComponentFractionKernel>( m_numComponents,
                                                  0, dataGroup->size(),
                                                  compDens,
                                                  dCompDens,
                                                  compFrac,
                                                  dCompFrac_dCompDens );
}

void CompositionalMultiphaseFlow::UpdateComponentFraction( localIndex er, localIndex esr ) const
{
  GEOSX_MARK_FUNCTION;

  KernelLaunchSelector1<ComponentFractionKernel>( m_numComponents,
                                                  0, m_compFrac[er][esr].size(0),
                                                  m_globalCompDensity[er][esr],
                                                  m_deltaGlobalCompDensity[er][esr],
                                                  m_compFrac[er][esr],
                                                  m_dCompFrac_dCompDens[er][esr] );
}

void CompositionalMultiphaseFlow::UpdatePhaseVolumeFraction( ManagedGroup * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d<real64> const & phaseVolFrac =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::phaseVolumeFractionString );

  arrayView2d<real64> const & dPhaseVolFrac_dPres =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64> const & dPhaseVolFrac_dComp =
    dataGroup->getReference< array3d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  // inputs

  arrayView3d<real64 const> const & dCompFrac_dCompDens =
    dataGroup->getReference< array3d<real64> >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  arrayView2d<real64 const> const & compDens =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dCompDens =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::deltaGlobalCompDensityString );

  MultiFluidBase * fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  arrayView3d<real64 const> const & phaseFrac =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::phaseFractionString );

  arrayView3d<real64 const> const & dPhaseFrac_dPres =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

  arrayView4d<real64 const> const & dPhaseFrac_dComp =
    fluid->getReference< array4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

  arrayView3d<real64 const> const & phaseDens =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dPhaseDens_dPres =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d<real64 const> const & dPhaseDens_dComp =
    fluid->getReference< array4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  KernelLaunchSelector2<PhaseVolumeFractionKernel>( m_numComponents, m_numPhases,
                                                    0, dataGroup->size(),
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

void CompositionalMultiphaseFlow::UpdatePhaseVolumeFraction( localIndex er, localIndex esr ) const
{
  GEOSX_MARK_FUNCTION;

  KernelLaunchSelector2<PhaseVolumeFractionKernel>( m_numComponents, m_numPhases,
                                                    0, m_phaseVolFrac[er][esr].size(0),
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

void CompositionalMultiphaseFlow::UpdatePhaseMobility( ManagedGroup * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d<real64> const & phaseMob =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::phaseMobilityString );

  arrayView2d<real64> const & dPhaseMob_dPres =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::dPhaseMobility_dPressureString );

  arrayView3d<real64> const & dPhaseMob_dComp =
    dataGroup->getReference< array3d<real64> >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  // inputs

  arrayView2d<real64 const> const & dPhaseVolFrac_dPres =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64 const> const & dPhaseVolFrac_dComp =
    dataGroup->getReference< array3d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d<real64 const> const & dCompFrac_dCompDens =
    dataGroup->getReference< array3d<real64> >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  MultiFluidBase const * fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  arrayView3d<real64 const> const & phaseDens =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dPhaseDens_dPres =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d<real64 const> const & dPhaseDens_dComp =
    fluid->getReference< array4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  arrayView3d<real64 const> const & phaseVisc =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::phaseViscosityString );

  arrayView3d<real64 const> const & dPhaseVisc_dPres =
    fluid->getReference< array3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

  arrayView4d<real64 const> const & dPhaseVisc_dComp =
    fluid->getReference< array4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString );

  RelativePermeabilityBase const * relperm = GetConstitutiveModel<RelativePermeabilityBase>( dataGroup, m_relPermName );

  arrayView3d<real64 const> const & phaseRelPerm =
    relperm->getReference< array3d<real64> >( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString );

  arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac =
    relperm->getReference< array4d<real64> >( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString );

  KernelLaunchSelector2<PhaseMobilityKernel>( m_numComponents, m_numPhases,
                                              0, dataGroup->size(),
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

void CompositionalMultiphaseFlow::UpdatePhaseMobility( localIndex er, localIndex esr ) const
{
  GEOSX_MARK_FUNCTION;

  KernelLaunchSelector2<PhaseMobilityKernel>( m_numComponents, m_numPhases,
                                              0, m_phaseMob[er][esr].size(0),
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

void CompositionalMultiphaseFlow::UpdateFluidModel( ManagedGroup * const dataGroup )
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference< array1d<real64> >( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );
  arrayView2d<real64 const> const & compFrac = dataGroup->getReference< array2d<real64> >( viewKeyStruct::globalCompFractionString );

  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    fluid->PointUpdate( pres[a] + dPres[a], m_temperature, compFrac[a], a, 0 );
  });
  //fluid->BatchUpdate( pres, temp, compFrac );
}

void CompositionalMultiphaseFlow::UpdateSolidModel( ManagedGroup * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase * const solid = GetConstitutiveModel<ConstitutiveBase>( dataGroup, m_solidName );

  arrayView1d<real64 const> const & pres  = dataGroup->getReference< array1d<real64> >( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    solid->StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  });
}

void CompositionalMultiphaseFlow::UpdateRelPermModel( ManagedGroup * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  RelativePermeabilityBase * const relPerm = GetConstitutiveModel<RelativePermeabilityBase>( dataGroup, m_relPermName );

  arrayView2d<real64> const & phaseVolFrac =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::phaseVolumeFractionString );

  relPerm->BatchUpdate( phaseVolFrac );
}

void CompositionalMultiphaseFlow::UpdateCapPressureModel( ManagedGroup * dataGroup )
{
  if (m_capPressureFlag)
  {
    CapillaryPressureBase * const capPressure = GetConstitutiveModel<CapillaryPressureBase>( dataGroup, m_capPressureName );

    arrayView2d<real64> const & phaseVolFrac =
      dataGroup->getReference< array2d<real64> >( viewKeyStruct::phaseVolumeFractionString );

    capPressure->BatchUpdate( phaseVolFrac );
  }
}

void CompositionalMultiphaseFlow::UpdateState( ManagedGroup * const dataGroup )
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

void CompositionalMultiphaseFlow::InitializeFluidState( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    // 1. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    UpdateFluidModel( subRegion );

    // 2. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    arrayView2d<real64 const> const & compFrac  = m_compFrac[er][esr];
    arrayView2d<real64 const> const & totalDens = m_totalDens[er][esr][m_fluidIndex];
    arrayView2d<real64> compDens = m_globalCompDensity[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    });

    // 3. Calculate phase saturations
    UpdatePhaseVolumeFraction( subRegion );

    // 4. Initialize solid state
    UpdateSolidModel( subRegion );

    // 5. Initialize rel perm state
    UpdateRelPermModel( subRegion );

    // 6. Initialize mobility
    UpdatePhaseMobility( subRegion );

    // 7. Initialize cap pressure state
    UpdateCapPressureModel( subRegion );

  });
}

void CompositionalMultiphaseFlow::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array> fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::globalCompDensityString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  // set mass fraction flag on main model
  // TODO find a way to set this before constitutive model is duplicated and attached to subregions?
  {
    MultiFluidBase * const fluid = constitutiveManager->GetConstitituveRelation<MultiFluidBase>( m_fluidName );
    fluid->setMassFlag( static_cast<bool>(m_useMass) );
  }

  // set mass fraction flag on subregion models
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( subRegion, m_fluidName );
    fluid->setMassFlag( static_cast<bool>(m_useMass) );
  });

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

  real64 dt_return = dt;

  ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

  // currently the only method is implicit time integration
  dt_return= this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          getLinearSystemRepository() );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void CompositionalMultiphaseFlow::BackupFields( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  // backup some fields used in time derivative approximation
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<real64  const> const & poroRef       = m_porosityRef[er][esr];

    arrayView2d<real64  const> const & phaseVolFrac = m_phaseVolFrac[er][esr];

    arrayView3d<real64  const> const & phaseDens     = m_phaseDens[er][esr][m_fluidIndex];
    arrayView4d<real64  const> const & phaseCompFrac = m_phaseCompFrac[er][esr][m_fluidIndex];
    arrayView2d<real64  const> const & pvMult        = m_pvMult[er][esr][m_solidIndex];

    arrayView2d<real64> const & phaseDensOld     = m_phaseDensOld[er][esr];
    arrayView2d<real64> const & phaseVolFracOld  = m_phaseVolFracOld[er][esr];
    arrayView3d<real64> const & phaseCompFracOld = m_phaseCompFracOld[er][esr];
    arrayView1d<real64> const & poroOld          = m_porosityOld[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      for (localIndex ip = 0; ip < m_numPhases; ++ip)
      {
        phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
        phaseVolFracOld[ei][ip] = phaseVolFrac[ei][ip];

        for (localIndex ic = 0; ic < m_numComponents; ++ic)
        {
          phaseCompFracOld[ei][ip][ic] = phaseCompFrac[ei][0][ip][ic];
        }
      }

      poroOld[ei] = poroRef[ei] * pvMult[ei][0];
    });
  });
}

void
CompositionalMultiphaseFlow::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                DomainPartition * const domain,
                                                EpetraBlockSystem * const blockSystem )
{
  // bind the stored views to the current domain
  ResetViews( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain );

  // setup dof numbers and linear system
  SetupSystem( domain, blockSystem );
}

void CompositionalMultiphaseFlow::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                                localIndex & numLocalRows,
                                                                globalIndex & numGlobalRows,
                                                                localIndex offset )
{
  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > const & dofNumber     = m_dofNumber;
  ElementRegionManager::ElementViewAccessor< arrayView1d<integer> >     const & elemGhostRank = m_elemGhostRank;

  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather(numMpiProcesses);

  // communicate the number of local rows to each process
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                gather.data(),
                                                1 );

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p = 0 ; p < numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if (p < thisMpiProcess)
      firstLocalRow += gather[p];
  }

  // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
  for( localIndex er=0 ; er < elemGhostRank.size() ; ++er )
  {
    for( localIndex esr=0 ; esr < elemGhostRank[er].size() ; ++esr )
    {
      dofNumber[er][esr] = -1;
    }
  }

  // loop over all elements and set the dof number if the element is not a ghost
  ReduceSum< reducePolicy, localIndex  > localCount(0);
  forAllElemsInMesh<RAJA::seq_exec>( meshLevel, GEOSX_LAMBDA ( localIndex const er,
                                                               localIndex const esr,
                                                               localIndex const ei )
  {
    if( elemGhostRank[er][esr][ei] < 0 )
    {
      dofNumber[er][esr][ei] = firstLocalRow + localCount + offset;
      localCount += 1;
    }
  });

  GEOS_ERROR_IF( localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}

void CompositionalMultiphaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                                      Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > const & dofNumber =
    elementRegionManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

  ElementRegionManager::ElementViewAccessor< arrayView1d<integer> > const & elemGhostRank =
    elementRegionManager->ConstructViewAccessor< array1d<integer>, arrayView1d<integer> >( ObjectManagerBase::viewKeyStruct::ghostRankString );

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  localIndex constexpr numElems   = FluxApproximationBase::CellStencil::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = FluxApproximationBase::CellStencil::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NDOF = m_numDofPerCell;

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  fluxApprox->forCellStencils( [&]( FluxApproximationBase::CellStencil const & stencil )
  {
    ArrayOfArraysView<FluxApproximationBase::CellStencil::Entry const, true> const & connections = stencil.getConnections();

    forall_in_range<stencilPolicy>( 0, connections.size(), GEOSX_LAMBDA ( localIndex iconn )
    {
      localIndex const stencilSize = connections.sizeOfArray( iconn );
      stackArray1d<globalIndex, numElems   * maxNumDof> elementLocalDofIndexRow( numElems * NDOF );
      stackArray1d<globalIndex, maxStencil * maxNumDof> elementLocalDofIndexCol( stencilSize * NDOF );

      for (localIndex i = 0; i < numElems; ++i)
      {
        CellDescriptor const & cell = connections( iconn, i ).index;
        globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];

        for (localIndex idof = 0; idof < NDOF; ++idof)
        {
          elementLocalDofIndexRow[i * NDOF + idof] = offset + idof;
        }
      }

      for (localIndex i = 0; i < stencilSize; ++i)
      {
        CellDescriptor const & cell = connections( iconn, i ).index;
        globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];

        for (localIndex idof = 0; idof < NDOF; ++idof)
        {
          elementLocalDofIndexCol[i * NDOF + idof] = offset + idof;
        }
      }

      sparsity->InsertGlobalIndices( integer_conversion<int>(numElems * m_numDofPerCell),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>(stencilSize * m_numDofPerCell),
                                     elementLocalDofIndexCol.data() );
    } );
  } );

  // loop over all elements and add all locals just in case the above connector loop missed some
  forAllElemsInMesh( meshLevel, [&] ( localIndex const er,
                                      localIndex const esr,
                                      localIndex const ei )
  {
    if (elemGhostRank[er][esr][ei] < 0)
    {
      stackArray1d<globalIndex, maxNumDof> elementLocalDofIndexRow( NDOF );

      globalIndex const offset = NDOF * dofNumber[er][esr][ei];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        elementLocalDofIndexRow[idof] = offset + idof;
      }

      sparsity->InsertGlobalIndices( integer_conversion<int>( NDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( NDOF ),
                                     elementLocalDofIndexRow.data() );
    }
  });

  // add additional connectivity resulting from boundary stencils
  fluxApprox->forBoundaryStencils( [&] ( FluxApproximationBase::BoundaryStencil const & stencil )
  {
    ArrayOfArraysView<FluxApproximationBase::BoundaryStencil::Entry const, true> const & connections = stencil.getConnections();

    forall_in_range<stencilPolicy>( 0, connections.size(), GEOSX_LAMBDA ( localIndex iconn )
    {
      localIndex const stencilSize = connections.sizeOfArray( iconn );
      stackArray1d<globalIndex, numElems   * maxNumDof> elementLocalDofIndexRow( numElems * NDOF );
      stackArray1d<globalIndex, maxStencil * maxNumDof> elementLocalDofIndexCol( stencilSize * NDOF );

      for (localIndex i = 0; i < numElems; ++i)
      {
        PointDescriptor const & point = connections( iconn, i ).index;

        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & cell = point.cellIndex;
          globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
          for (localIndex idof = 0; idof < NDOF; ++idof)
          {
            elementLocalDofIndexRow[idof] = offset + idof;
          }
        }
      }

      integer counter = 0;
      for (localIndex i = 0; i < stencilSize; ++i)
      {
        PointDescriptor const & point = connections( iconn, i ).index;

        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & cell = point.cellIndex;
          globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
          for (localIndex idof = 0; idof < NDOF; ++idof)
          {
            elementLocalDofIndexCol[counter * NDOF + idof] = offset + idof;
          }
          ++counter;
        }
      }

      sparsity->InsertGlobalIndices( integer_conversion<int>( NDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( counter * NDOF ),
                                     elementLocalDofIndexCol.data() );
    });
  });
}

void CompositionalMultiphaseFlow::SetupSystem( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem )
{
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  // for this solver, the dof are on the cell center, and the block of rows corresponds to a cell
  localIndex numLocalRows  = 0;
  globalIndex numGlobalRows = 0;

  // get the number of local elements, and ghost elements...i.e. local rows and ghost rows
  elementRegionManager->forElementSubRegions( [&]( ObjectManagerBase * const subRegion )
  {
    numLocalRows += subRegion->size() - subRegion->GetNumberOfGhosts();
  });

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                0 );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::blockLocalDofNumberString );
  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


  // construct row map, and set a pointer to the row map
  Epetra_Map * const rowMap =
    blockSystem->SetRowMap( BlockIDs::compositionalBlock,
                            std::make_unique<Epetra_Map>( numGlobalRows * m_numDofPerCell,
                                                          numLocalRows * m_numDofPerCell,
                                                          0,
                                                          m_linearSolverWrapper.m_epetraComm ) );

  // construct sparsity matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const sparsity =
    blockSystem->SetSparsity( BlockIDs::compositionalBlock,
                              BlockIDs::compositionalBlock,
                              std::make_unique<Epetra_FECrsGraph>(Copy, *rowMap, 0) );



  // set the sparsity patter
  SetSparsityPattern( domain, sparsity );

  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( BlockIDs::compositionalBlock,
                          BlockIDs::compositionalBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  // construct solution vector
  blockSystem->SetSolutionVector( BlockIDs::compositionalBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( BlockIDs::compositionalBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );
}

void CompositionalMultiphaseFlow::AssembleSystem( DomainPartition * const domain,
                                                  EpetraBlockSystem * const blockSystem,
                                                  real64 const time_n, real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleVolumeBalanceTerms( domain, jacobian, residual, time_n, dt );

  jacobian->GlobalAssemble( true );
  residual->GlobalAssemble();

  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After CompositionalMultiphaseFlow::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

}

void CompositionalMultiphaseFlow::AssembleAccumulationTerms( DomainPartition const * const domain,
                                                             Epetra_FECrsMatrix * const jacobian,
                                                             Epetra_FEVector * const residual,
                                                             real64 const time_n,
                                                             real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer     const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64 const> const & volume      = m_volume[er][esr];
    arrayView1d<real64 const> const & porosityRef = m_porosityRef[er][esr];

    arrayView2d<real64 const> const & phaseVolFrac            = m_phaseVolFrac[er][esr];
    arrayView2d<real64 const> const & dPhaseVolFrac_dPres     = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d<real64 const> const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];
    arrayView3d<real64 const> const & dCompFrac_dCompDens     = m_dCompFrac_dCompDens[er][esr];

    arrayView1d<real64 const> const & porosityOld      = m_porosityOld[er][esr];
    arrayView2d<real64 const> const & phaseVolFracOld  = m_phaseVolFracOld[er][esr];
    arrayView2d<real64 const> const & phaseDensOld     = m_phaseDensOld[er][esr];
    arrayView3d<real64 const> const & phaseCompFracOld = m_phaseCompFracOld[er][esr];

    arrayView2d<real64 const> const & pvMult               = m_pvMult[er][esr][m_solidIndex];
    arrayView2d<real64 const> const & dPvMult_dPres        = m_dPvMult_dPres[er][esr][m_solidIndex];
    arrayView3d<real64 const> const & phaseDens            = m_phaseDens[er][esr][m_fluidIndex];
    arrayView3d<real64 const> const & dPhaseDens_dPres     = m_dPhaseDens_dPres[er][esr][m_fluidIndex];
    arrayView4d<real64 const> const & dPhaseDens_dComp     = m_dPhaseDens_dComp[er][esr][m_fluidIndex];
    arrayView4d<real64 const> const & phaseCompFrac        = m_phaseCompFrac[er][esr][m_fluidIndex];
    arrayView4d<real64 const> const & dPhaseCompFrac_dPres = m_dPhaseCompFrac_dPres[er][esr][m_fluidIndex];
    arrayView5d<real64 const> const & dPhaseCompFrac_dComp = m_dPhaseCompFrac_dComp[er][esr][m_fluidIndex];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        stackArray1d<long long, maxNumDof>           localAccumDOF( NDOF );
        stackArray1d<real64, maxNumComp>             localAccum( NC );
        stackArray2d<real64, maxNumComp * maxNumDof> localAccumJacobian( NC, NDOF );

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
        globalIndex const offset = NDOF * dofNumber[ei];
        for (localIndex idof = 0; idof < NDOF; ++idof)
        {
          localAccumDOF[idof] = integer_conversion<long long>(offset + idof);
        }

        // TODO: apply equation/variable change transformation(s)

        // add contribution to global residual and dRdP
        residual->SumIntoGlobalValues( integer_conversion<int>( NC ),
                                       localAccumDOF.data(),
                                       localAccum.data() );

        jacobian->SumIntoGlobalValues( integer_conversion<int>( NC ),
                                       localAccumDOF.data(),
                                       integer_conversion<int>( NDOF ),
                                       localAccumDOF.data(),
                                       localAccumJacobian.data(),
                                       Epetra_FECrsMatrix::ROW_MAJOR );
      }
    });
  });
}

void CompositionalMultiphaseFlow::AssembleFluxTerms( DomainPartition const * const domain,
                                                     Epetra_FECrsMatrix * const jacobian,
                                                     Epetra_FEVector * const residual,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  NumericalMethodsManager const * const numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * const fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  FluxKernel::ElementView< arrayView1d<globalIndex> > const & blockLocalDofNumber = m_dofNumber.toViewConst();

  FluxKernel::ElementView< arrayView1d<real64 const> > const & pres                = m_pressure.toViewConst();
  FluxKernel::ElementView< arrayView1d<real64 const> > const & dPres               = m_deltaPressure.toViewConst();
  FluxKernel::ElementView< arrayView1d<real64 const> > const & gravDepth           = m_gravDepth.toViewConst();
  FluxKernel::ElementView< arrayView2d<real64 const> > const & phaseMob            = m_phaseMob.toViewConst();
  FluxKernel::ElementView< arrayView2d<real64 const> > const & dPhaseMob_dPres     = m_dPhaseMob_dPres.toViewConst();
  FluxKernel::ElementView< arrayView3d<real64 const> > const & dPhaseMob_dComp     = m_dPhaseMob_dCompDens.toViewConst();
  FluxKernel::ElementView< arrayView2d<real64 const> > const & dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres.toViewConst();
  FluxKernel::ElementView< arrayView3d<real64 const> > const & dPhaseVolFrac_dComp = m_dPhaseVolFrac_dCompDens.toViewConst();
  FluxKernel::ElementView< arrayView3d<real64 const> > const & dCompFrac_dCompDens = m_dCompFrac_dCompDens.toViewConst();

  FluxKernel::MaterialView< arrayView3d<real64 const> > const & phaseDens                   = m_phaseDens.toViewConst();
  FluxKernel::MaterialView< arrayView3d<real64 const> > const & dPhaseDens_dPres            = m_dPhaseDens_dPres.toViewConst();
  FluxKernel::MaterialView< arrayView4d<real64 const> > const & dPhaseDens_dComp            = m_dPhaseDens_dComp.toViewConst();
  FluxKernel::MaterialView< arrayView4d<real64 const> > const & phaseCompFrac               = m_phaseCompFrac.toViewConst();
  FluxKernel::MaterialView< arrayView4d<real64 const> > const & dPhaseCompFrac_dPres        = m_dPhaseCompFrac_dPres.toViewConst();
  FluxKernel::MaterialView< arrayView5d<real64 const> > const & dPhaseCompFrac_dComp        = m_dPhaseCompFrac_dComp.toViewConst();
  FluxKernel::MaterialView< arrayView3d<real64 const> > const & phaseCapPres                = m_phaseCapPressure.toViewConst();
  FluxKernel::MaterialView< arrayView4d<real64 const> > const & dPhaseCapPres_dPhaseVolFrac = m_dPhaseCapPressure_dPhaseVolFrac.toViewConst();

  localIndex constexpr numElems   = FluxApproximationBase::CellStencil::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = FluxApproximationBase::CellStencil::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex constexpr maxSize1 = numElems * maxNumComp;
  localIndex constexpr maxSize2 = maxStencil * maxNumDof;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  localIndex const fluidIndex       = m_fluidIndex;
  localIndex const capPressureIndex = m_capPressureIndex;
  integer const gravityFlag     = m_gravityFlag;
  integer const capPressureFlag = m_capPressureFlag;

  fluxApprox->forCellStencils( [&] ( FluxApproximationBase::CellStencil const & stencil )
  {
    ArrayOfArraysView<FluxApproximationBase::CellStencil::Entry const, true> const & connections = stencil.getConnections();

    forall_in_range<stencilPolicy>( 0, connections.size(), GEOSX_LAMBDA ( localIndex iconn )
    {
      localIndex const stencilSize = connections.sizeOfArray(iconn);

      // create local work arrays
      stackArray1d<long long, maxSize1> eqnRowIndices( numElems * NC );
      stackArray1d<long long, maxSize2> dofColIndices( stencilSize * NDOF );

      stackArray1d<double, maxSize1>            localFlux( numElems * NC );
      stackArray2d<double, maxSize1 * maxSize2> localFluxJacobian( numElems * NC, stencilSize * NDOF );

      FluxKernel::Compute( NC, NP,
                           stencilSize,
                           connections[iconn],
                           pres,
                           dPres,
                           gravDepth,
                           phaseMob,
                           dPhaseMob_dPres,
                           dPhaseMob_dComp,
                           dPhaseVolFrac_dPres,
                           dPhaseVolFrac_dComp,
                           dCompFrac_dCompDens,
                           phaseDens ,
                           dPhaseDens_dPres,
                           dPhaseDens_dComp,
                           phaseCompFrac,
                           dPhaseCompFrac_dPres,
                           dPhaseCompFrac_dComp,
                           phaseCapPres,
                           dPhaseCapPres_dPhaseVolFrac,
                           fluidIndex,
                           capPressureIndex,
                           gravityFlag,
                           capPressureFlag,
                           dt,
                           localFlux,
                           localFluxJacobian );

      // set equation indices for both connected cells
      for (localIndex i = 0; i < numElems; ++i)
      {
        CellDescriptor const & cell = connections(iconn, i).index;
        globalIndex const offset = NDOF * blockLocalDofNumber[cell.region][cell.subRegion][cell.index];

        for (localIndex ic = 0; ic < NC; ++ic)
        {
          eqnRowIndices[i * NC + ic] = offset + ic;
        }
      }

      for (localIndex i = 0; i < stencilSize; ++i)
      {
        CellDescriptor const & cell = connections(iconn, i).index;
        globalIndex const offset = NDOF * blockLocalDofNumber[cell.region][cell.subRegion][cell.index];

        for (localIndex jdof = 0; jdof < NDOF; ++jdof)
        {
          dofColIndices[i * NDOF + jdof] = offset + jdof;
        }
      }

      // TODO: apply equation/variable change transformation(s)

      // Add to global residual/jacobian
      residual->SumIntoGlobalValues( integer_conversion<int>( numElems * NC ),
                                     eqnRowIndices.data(),
                                     localFlux.data() );

      jacobian->SumIntoGlobalValues( integer_conversion<int>( numElems * NC ),
                                     eqnRowIndices.data(),
                                     integer_conversion<int>( stencilSize * NDOF ),
                                     dofColIndices.data(),
                                     localFluxJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR );

    } );
  } );
}

void CompositionalMultiphaseFlow::AssembleVolumeBalanceTerms( DomainPartition const * const domain,
                                                              Epetra_FECrsMatrix * const jacobian,
                                                              Epetra_FEVector * const residual,
                                                              real64 const time_n,
                                                              real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64 const> const & volume      = m_volume[er][esr];
    arrayView1d<real64 const> const & porosityRef = m_porosityRef[er][esr];

    arrayView2d<real64 const> const & phaseVolFrac            = m_phaseVolFrac[er][esr];
    arrayView2d<real64 const> const & dPhaseVolFrac_dPres     = m_dPhaseVolFrac_dPres[er][esr];
    arrayView3d<real64 const> const & dPhaseVolFrac_dCompDens = m_dPhaseVolFrac_dCompDens[er][esr];

    arrayView2d<real64 const> const & pvMult        = m_pvMult[er][esr][m_solidIndex];
    arrayView2d<real64 const> const & dPvMult_dPres = m_dPvMult_dPres[er][esr][m_solidIndex];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] >= 0)
      {
        return;
      }

      real64                             localVolBalance;
      stackArray1d<real64, maxNumDof>    localVolBalanceJacobian( NDOF );
      stackArray1d<long long, maxNumDof> localVolBalanceDOF( NDOF );

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
      globalIndex const offset = NDOF * dofNumber[ei];
      globalIndex const localVolBalanceEqnIndex = offset + NC;
      for (localIndex jdof = 0; jdof < NDOF; ++jdof)
      {
        localVolBalanceDOF[jdof] = offset + jdof;
      }

      // TODO: apply equation/variable change transformation(s)

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues( 1,
                                     &localVolBalanceEqnIndex,
                                     &localVolBalance );

      jacobian->SumIntoGlobalValues( 1,
                                     &localVolBalanceEqnIndex,
                                     integer_conversion<int>( NDOF ),
                                     localVolBalanceDOF.data(),
                                     localVolBalanceJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR );
    } );
  } );
}

void CompositionalMultiphaseFlow::ApplyBoundaryConditions( DomainPartition * const domain,
                                                           EpetraBlockSystem * const blockSystem,
                                                           real64 const time_n, real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit(domain, time_n, dt, blockSystem);
  ApplyFaceDirichletBC_implicit(domain, time_n, dt, blockSystem);

  jacobian->GlobalAssemble();
  residual->GlobalAssemble();

  if (verboseLevel() >= 2)
  {
    GEOS_LOG_RANK( "After CompositionalMultiphaseFlow::ApplyBoundaryCondition" );
    GEOS_LOG_RANK( "\nJacobian\n" << *jacobian );
    GEOS_LOG_RANK( "\nResidual\n" << *residual );
  }
}

void
CompositionalMultiphaseFlow::ApplyDirichletBC_implicit( DomainPartition * const domain,
                                                        real64 const time_n, real64 const dt,
                                                        EpetraBlockSystem * const blockSystem )
{
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  unordered_map< string, array1d<bool> > bcStatusMap; // map to check consistent application of BC

  // 1. apply pressure Dirichlet BCs
  fsManager->Apply( time_n + dt,
                    domain,
                    "ElementRegions",
                    viewKeyStruct::pressureString,
                    [&]( FieldSpecificationBase const * const fs,
                    string const & setName,
                    set<localIndex> const & targetSet,
                    ManagedGroup * subRegion,
                    string const & )
  {
    // 1.0. Check whether pressure has already been applied to this set
    GEOS_ERROR_IF( bcStatusMap.count( setName ) > 0, "Conflicting pressure boundary conditions on set " << setName );
    bcStatusMap[setName].resize( m_numComponents );
    bcStatusMap[setName] = false;

    // 1.1. Apply BC to set the field values
    fs->ApplyFieldValue<FieldSpecificationEqual>( targetSet,
                                                  time_n + dt,
                                                  subRegion,
                                                  viewKeyStruct::bcPressureString );
  });

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  fsManager->Apply( time_n + dt,
                    domain,
                    "ElementRegions",
                    viewKeyStruct::globalCompFractionString,
                    [&] ( FieldSpecificationBase const * const fs,
                    string const & setName,
                    set<localIndex> const & targetSet,
                    ManagedGroup * subRegion,
                    string const & )
  {
    // 2.0. Check pressure and record composition bc application
    localIndex const comp = fs->GetComponent();
    GEOS_ERROR_IF( bcStatusMap.count( setName ) == 0, "Pressure boundary condition not prescribed on set '" << setName << "'" );
    GEOS_ERROR_IF( bcStatusMap[setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
    bcStatusMap[setName][comp] = true;

    // 2.1. Apply BC to set the field values
    fs->ApplyFieldValue<FieldSpecificationEqual>( targetSet,
                                                  time_n + dt,
                                                  subRegion,
                                                  viewKeyStruct::globalCompFractionString );
  });

  // 2.3 Check consistency between composition BC applied to sets
  bool bcConsistent = true;
  for (auto const & bcEntry : bcStatusMap)
  {
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      bcConsistent &= bcEntry.second[ic];
      GEOS_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component "
                                      << ic << " on set '" << bcEntry.first << "'" );
    }
  }
  GEOS_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

  // 3. Call constitutive update, back-calculate target global component densities and apply to the system
  fsManager->Apply( time_n + dt,
                    domain,
                    "ElementRegions",
                    viewKeyStruct::pressureString,
                    [&] ( FieldSpecificationBase const * const bc,
                    string const & setName,
                    set<localIndex> const & targetSet,
                    ManagedGroup * subRegion,
                    string const & )
  {
    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( subRegion, m_fluidName );

    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    arrayView1d<real64 const> const & pres      = subRegion->getReference< array1d<real64> >( viewKeyStruct::pressureString );
    arrayView1d<real64 const> const & dPres     = subRegion->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );
    arrayView1d<real64 const> const & bcPres    = subRegion->getReference< array1d<real64> >( viewKeyStruct::bcPressureString );
    arrayView2d<real64 const> const & compFrac  = subRegion->getReference< array2d<real64> >( viewKeyStruct::globalCompFractionString );
    arrayView2d<real64 const> const & compDens  = subRegion->getReference< array2d<real64> >( viewKeyStruct::globalCompDensityString );
    arrayView2d<real64 const> const & dCompDens = subRegion->getReference< array2d<real64> >( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d<real64> const & totalDens = fluid->getReference< array2d<real64> >( MultiFluidBase::viewKeyStruct::totalDensityString );

    Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

    array1d<real64> rhsContribution( targetSet.size() * m_numDofPerCell );
    array1d<globalIndex> dof( targetSet.size() * m_numDofPerCell );

    integer counter = 0;

    for (localIndex a : targetSet)
    {
      fluid->PointUpdate( bcPres[a], m_temperature, compFrac[a], a, 0 );

      globalIndex const offset = m_numDofPerCell * dofNumber[a];
      dof[counter] = offset;

      // 4.1. Apply pressure to the matrix
      FieldSpecificationEqual::SpecifyFieldValue( dof[counter],
                                                  blockSystem,
                                                  BlockIDs::compositionalBlock,
                                                  rhsContribution[counter],
                                                  bcPres[a],
                                                  pres[a] + dPres[a] );

      ++counter;

      // 4.2. For each component, apply target global density value
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        dof[counter] = offset + ic + 1;
        real64 const targetCompDens = totalDens[a][0] * compFrac[a][ic];

        FieldSpecificationEqual::SpecifyFieldValue( dof[counter],
                                                    blockSystem,
                                                    BlockIDs::compositionalBlock,
                                                    rhsContribution[counter],
                                                    targetCompDens,
                                                    compDens[a][ic] + dCompDens[a][ic] );

        ++counter;
      }

      // 4.3. Apply accumulated rhs values
      FieldSpecificationEqual::ReplaceGlobalValues( rhs,
                                                    counter,
                                                    dof.data(),
                                                    rhsContribution.data() );
    }
  });

}

void
CompositionalMultiphaseFlow::ApplyFaceDirichletBC_implicit( DomainPartition * const domain,
                                                            real64 const time_n, real64 const dt,
                                                            EpetraBlockSystem * const blockSystem )
{

}

real64
CompositionalMultiphaseFlow::CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                                                    DomainPartition * const domain )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::compositionalBlock );

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  real64 localResidualNorm = 0.0;
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];
    arrayView1d<real64 const>      const & refPoro       = m_porosityRef[er][esr];
    arrayView1d<real64 const>      const & volume        = m_volume[er][esr];

    localResidualNorm += sum_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        real64 cell_norm = 0.0;
        globalIndex const offset = m_numDofPerCell * dofNumber[ei];
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          int const lid = rowMap->LID(integer_conversion<int>(offset + idof));
          real64 const val = localResidual[lid] / (refPoro[ei] * volume[ei]);
          cell_norm += val * val;
        }
        return cell_norm;
      }
      return 0.0;
    } );
  } );

  // compute global residual norm
  realT globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

void CompositionalMultiphaseFlow::SolveSystem( EpetraBlockSystem * const blockSystem,
                                               SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
    solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  Epetra_FEVector * const
    residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  residual->Scale(-1.0);
  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                BlockIDs::compositionalBlock );

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("\nSolution:\n" << *solution);
  }
}

bool
CompositionalMultiphaseFlow::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  RAJA::ReduceMin<reducePolicy, integer> result(1);

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer     const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64 const> const & pres      = m_pressure[er][esr];
    arrayView1d<real64 const> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64 const> const & compDens  = m_globalCompDensity[er][esr];
    arrayView2d<real64 const> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      globalIndex const offset = m_numDofPerCell * dofNumber[ei];
      // extract solution and apply to dP
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset));
        real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * local_solution[lid];
        if (newPres < 0.0)
        {
          result.min(0);
        }
      }

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        int const lid = rowMap->LID(integer_conversion<int>(offset + ic + 1));
        real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * local_solution[lid];
        if (newDens < 0.0)
        {
          result.min(0);
        }
      }
    });
  });

  return result.get() > 0;
}

void
CompositionalMultiphaseFlow::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  int dummy;
  double * local_solution = nullptr;
  solution->ExtractView( &local_solution, &dummy );

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<integer     const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] >= 0)
        return;

      globalIndex const offset = m_numDofPerCell * dofNumber[ei];
      // extract solution and apply to dP
      {
        int const lid = rowMap->LID( integer_conversion<int>( offset ) );
        dPres[ei] += scalingFactor * local_solution[lid];
      }

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        int const lid = rowMap->LID( integer_conversion<int>( offset + ic + 1 ) );
        dCompDens[ei][ic] += scalingFactor * local_solution[lid];
      }
    });
  });

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaGlobalCompDensityString );
  CommunicationTools::SynchronizeFields( fieldNames,
                                         mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );
  } );
}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion * const elementRegion,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
        dCompDens[ei][ic] = 0.0;
    });
  });

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );
  } );
}

void CompositionalMultiphaseFlow::ImplicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion * const elementRegion,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64 const> const & dPres     = m_deltaPressure[er][esr];
    arrayView2d<real64 const> const & dCompDens = m_deltaGlobalCompDensity[er][esr];

    arrayView1d<real64> const & pres     = m_pressure[er][esr];
    arrayView2d<real64> const & compDens = m_globalCompDensity[er][esr];

    forall_in_range<elemPolicy>(0,subRegion->size(), GEOSX_LAMBDA ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
        compDens[ei][ic] += dCompDens[ei][ic];
    });
  });
}

void CompositionalMultiphaseFlow::ResetViews(DomainPartition * const domain)
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_dofNumber =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );
  m_pressure =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::deltaPressureString );
  m_globalCompDensity =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::globalCompDensityString );
  m_deltaGlobalCompDensity =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::deltaGlobalCompDensityString );

  m_compFrac =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::globalCompFractionString );
  m_dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor< array3d<real64>, arrayView3d<real64> >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
  m_phaseVolFrac =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::phaseVolumeFractionString );
  m_dPhaseVolFrac_dPres =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  m_dPhaseVolFrac_dCompDens =
    elemManager->ConstructViewAccessor< array3d<real64>, arrayView3d<real64> >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );
  m_phaseMob =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::phaseMobilityString );
  m_dPhaseMob_dPres =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::dPhaseMobility_dPressureString );
  m_dPhaseMob_dCompDens =
    elemManager->ConstructViewAccessor< array3d<real64>, arrayView3d<real64> >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  m_porosityOld =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::porosityOldString );
  m_phaseVolFracOld =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::phaseVolumeFractionOldString );
  m_phaseDensOld =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::phaseDensityOldString );
  m_phaseCompFracOld =
    elemManager->ConstructViewAccessor< array3d<real64>, arrayView3d<real64> >( viewKeyStruct::phaseComponentFractionOldString );

  m_pvMult =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                      constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                      constitutiveManager );
  m_phaseFrac =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::phaseFractionString,
                                                                                      constitutiveManager );
  m_dPhaseFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array4d<real64>, arrayView4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseDens =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                      constitutiveManager );
  m_dPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseDens_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array4d<real64>, arrayView4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseVisc =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::phaseViscosityString,
                                                                                      constitutiveManager );
  m_dPhaseVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseVisc_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array4d<real64>, arrayView4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseCompFrac =
    elemManager->ConstructFullMaterialViewAccessor< array4d<real64>, arrayView4d<real64> >( MultiFluidBase::viewKeyStruct::phaseCompFractionString,
                                                                                      constitutiveManager );
  m_dPhaseCompFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array4d<real64>, arrayView4d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseCompFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor< array5d<real64>, arrayView5d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_totalDens =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >( MultiFluidBase::viewKeyStruct::totalDensityString,
                                                                                      constitutiveManager );
  m_phaseRelPerm =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString,
                                                                                      constitutiveManager );
  m_dPhaseRelPerm_dPhaseVolFrac =
    elemManager->ConstructFullMaterialViewAccessor< array4d<real64>, arrayView4d<real64> >( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                      constitutiveManager );
  if (m_capPressureFlag)
  {
    m_phaseCapPressure =
      elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( CapillaryPressureBase::viewKeyStruct::phaseCapPressureString,
                                                                                        constitutiveManager );
    m_dPhaseCapPressure_dPhaseVolFrac =
      elemManager->ConstructFullMaterialViewAccessor< array4d<real64>, arrayView4d<real64> >( CapillaryPressureBase::viewKeyStruct::dPhaseCapPressure_dPhaseVolFractionString,
                                                                                        constitutiveManager );
  }
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseFlow, string const &, ManagedGroup * const)
}// namespace geosx
