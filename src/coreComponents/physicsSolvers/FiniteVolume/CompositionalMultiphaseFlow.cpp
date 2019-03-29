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

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

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
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementSubRegions([&]( ElementSubRegionBase * const elementSubRegion) -> void
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

  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  elemManager->forElementSubRegions( [&] ( ObjectManagerBase * const subRegion )
  {
    subRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompDensityString).resizeDimension<1>(NC);
    subRegion->getReference<array2d<real64>>(viewKeyStruct::deltaGlobalCompDensityString).resizeDimension<1>(NC);

    subRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompFractionString).resizeDimension<1>(NC);
    subRegion->getReference<array3d<real64>>(viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString).resizeDimension<1,2>(NC, NC);

    subRegion->getReference<array2d<real64>>(viewKeyStruct::phaseVolumeFractionString).resizeDimension<1>(NP);
    subRegion->getReference<array2d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dPressureString).resizeDimension<1>(NP);
    subRegion->getReference<array3d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

    subRegion->getReference<array2d<real64>>(viewKeyStruct::phaseMobilityString).resizeDimension<1>(NP);
    subRegion->getReference<array2d<real64>>(viewKeyStruct::dPhaseMobility_dPressureString).resizeDimension<1>(NP);
    subRegion->getReference<array3d<real64>>(viewKeyStruct::dPhaseMobility_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

    subRegion->getReference<array2d<real64>>(viewKeyStruct::phaseVolumeFractionOldString).resizeDimension<1>(NP);
    subRegion->getReference<array2d<real64>>(viewKeyStruct::phaseDensityOldString).resizeDimension<1>(NP);
    subRegion->getReference<array3d<real64>>(viewKeyStruct::phaseComponentFractionOldString).resizeDimension<1,2>(NP, NC);
  });
}

void CompositionalMultiphaseFlow::UpdateComponentFraction( ManagedGroup * const dataGroup )
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d<real64> const & compFrac =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  arrayView3d<real64> const & dCompFrac_dCompDens =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // inputs

  arrayView2d<real64 const> const & compDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dCompDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    real64 totalDensity = 0.0;

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      totalDensity += compDens[a][ic] + dCompDens[a][ic];
    }

    real64 const totalDensityInv = 1.0 / totalDensity;

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      compFrac[a][ic] = (compDens[a][ic] + dCompDens[a][ic]) * totalDensityInv;
      for (localIndex jc = 0; jc < m_numComponents; ++jc)
      {
        dCompFrac_dCompDens[a][ic][jc] = - compFrac[a][ic] * totalDensityInv;
      }
      dCompFrac_dCompDens[a][ic][ic] += totalDensityInv;
    }
  });
}

void CompositionalMultiphaseFlow::UpdatePhaseVolumeFraction( ManagedGroup * const dataGroup )
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d<real64> const & phaseVolFrac =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

  arrayView2d<real64> const & dPhaseVolFrac_dPres =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64> const & dPhaseVolFrac_dComp =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  // inputs

  arrayView3d<real64 const> const & dCompFrac_dCompDens =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  arrayView2d<real64 const> const & compDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dCompDens =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  MultiFluidBase * fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  arrayView3d<real64 const> const & phaseFrac =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString );

  arrayView3d<real64 const> const & dPhaseFrac_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

  arrayView4d<real64 const> const & dPhaseFrac_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

  arrayView3d<real64 const> const & phaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d<real64 const> const & dPhaseDens_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    stackArray1d<real64, maxNumComp> work( NC );

    // compute total density from component partial densities
    real64 totalDensity = 0.0;
    real64 const dTotalDens_dCompDens = 1.0;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      totalDensity += compDens[a][ic] + dCompDens[a][ic];
    }

    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
      real64 const phaseDensInv = 1.0 / phaseDens[a][0][ip];

      // compute saturation and derivatives except multiplying by the total density
      phaseVolFrac[a][ip] = phaseFrac[a][0][ip] * phaseDensInv;

      dPhaseVolFrac_dPres[a][ip] =
        (dPhaseFrac_dPres[a][0][ip] - phaseVolFrac[a][ip] * dPhaseDens_dPres[a][0][ip]) * phaseDensInv;

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseVolFrac_dComp[a][ip][jc] =
          (dPhaseFrac_dComp[a][0][ip][jc] - phaseVolFrac[a][ip] * dPhaseDens_dComp[a][0][ip][jc]) * phaseDensInv;
      }

      // apply chain rule to convert derivatives from global component fractions to densities
      applyChainRuleInPlace( NC, dCompFrac_dCompDens[a], dPhaseVolFrac_dComp[a][ip], work );

      // now finalize the computation by multiplying by total density
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseVolFrac_dComp[a][ip][jc] *= totalDensity;
        dPhaseVolFrac_dComp[a][ip][jc] += phaseVolFrac[a][ip] * dTotalDens_dCompDens;
      }

      phaseVolFrac[a][ip] *= totalDensity;
      dPhaseVolFrac_dPres[a][ip] *= totalDensity;
    }
  });
}

void CompositionalMultiphaseFlow::UpdatePhaseMobility( ManagedGroup * const dataGroup )
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d<real64> const & phaseMob =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::phaseMobilityString );

  arrayView2d<real64> const & dPhaseMob_dPres =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::dPhaseMobility_dPressureString );

  arrayView3d<real64> const & dPhaseMob_dComp =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  // inputs

  arrayView2d<real64 const> const & dPhaseVolFrac_dPres =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64 const> const & dPhaseVolFrac_dComp =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d<real64 const> const & dCompFrac_dCompDens =
    dataGroup->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  MultiFluidBase const * fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  arrayView3d<real64 const> const & phaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d<real64 const> const & dPhaseDens_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  arrayView3d<real64 const> const & phaseVisc =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString );

  arrayView3d<real64 const> const & dPhaseVisc_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

  arrayView4d<real64 const> const & dPhaseVisc_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString );

  RelativePermeabilityBase const * relperm = GetConstitutiveModel<RelativePermeabilityBase>( dataGroup, m_relPermName );

  arrayView3d<real64 const> const & phaseRelPerm =
    relperm->getReference<array3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString );

  arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac =
    relperm->getReference<array4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString );

  localIndex constexpr maxNumComp  = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    stackArray1d<real64, maxNumComp> dRelPerm_dC( NC );
    stackArray1d<real64, maxNumComp> dDens_dC( NC );
    stackArray1d<real64, maxNumComp> dVisc_dC( NC );

    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const density = phaseDens[a][0][ip];
      real64 const dDens_dP = dPhaseDens_dPres[a][0][ip];
      applyChainRule( NC, dCompFrac_dCompDens[a], dPhaseDens_dComp[a][0][ip], dDens_dC );

      real64 const viscosity = phaseVisc[a][0][ip];
      real64 const dVisc_dP = dPhaseVisc_dPres[a][0][ip];
      applyChainRule( NC, dCompFrac_dCompDens[a], dPhaseVisc_dComp[a][0][ip], dVisc_dC );

      real64 const relPerm = phaseRelPerm[a][0][ip];
      real64 dRelPerm_dP = 0.0;
      dRelPerm_dC = 0.0;

      for (localIndex jp = 0; jp < NP; ++jp)
      {
        real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[a][0][ip][jp];
        dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[a][jp];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[a][jp][jc];
        }
      }

      real64 const mobility = relPerm * density / viscosity;

      phaseMob[a][ip] = mobility;
      dPhaseMob_dPres[a][ip] = dRelPerm_dP * density / viscosity
                             + mobility * (dDens_dP / density - dVisc_dP / viscosity);

      // compositional derivatives
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseMob_dComp[a][ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                   + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
      }
    }
  } );
}

void CompositionalMultiphaseFlow::UpdateFluidModel( ManagedGroup * const dataGroup )
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
  arrayView2d<real64 const> const & compFrac = dataGroup->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

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

  arrayView1d<real64 const> const & pres  = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

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
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

  relPerm->BatchUpdate( phaseVolFrac );
}

void CompositionalMultiphaseFlow::UpdateCapPressureModel( ManagedGroup * dataGroup )
{
  if (m_capPressureFlag)
  {
    CapillaryPressureBase * const capPressure = GetConstitutiveModel<CapillaryPressureBase>( dataGroup, m_capPressureName );

    arrayView2d<real64> const & phaseVolFrac =
      dataGroup->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

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

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
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

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array> fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::globalCompDensityString );
  CommunicationTools::SynchronizeFields(fieldNames,
                                        domain->getMeshBody(0)->getMeshLevel(0),
                                        domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  ConstitutiveManager * constitutiveManager = domain->getConstitutiveManager();

  // set mass fraction flag on main model
  // TODO find a way to set this before constitutive model is duplicated and attached to subregions?
  {
    MultiFluidBase * fluid = constitutiveManager->GetConstitituveRelation<MultiFluidBase>( m_fluidName );
    fluid->setMassFlag( static_cast<bool>(m_useMass) );
  }

  // set mass fraction flag on subregion models
  applyToSubRegions( domain, [&] ( ManagedGroup * subRegion )
  {
    MultiFluidBase * fluid = GetConstitutiveModel<MultiFluidBase>( subRegion, m_fluidName );
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
  // backup some fields used in time derivative approximation
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
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
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber     = m_dofNumber;
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>>     const & elemGhostRank = m_elemGhostRank;

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

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber =
    elementRegionManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( viewKeyStruct::blockLocalDofNumberString );

  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> const & elemGhostRank =
    elementRegionManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  localIndex constexpr numElems   = StencilCollection<CellDescriptor, real64>::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NDOF = m_numDofPerCell;

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  stencilCollection.forAll( GEOSX_LAMBDA ( StencilCollection<CellDescriptor, real64>::Accessor stencil )
  {
    localIndex const stencilSize = stencil.size();
    stackArray1d<globalIndex, numElems   * maxNumDof> elementLocalDofIndexRow( numElems * NDOF );
    stackArray1d<globalIndex, maxStencil * maxNumDof> elementLocalDofIndexCol( stencilSize * NDOF );

    stencil.forConnected( [&] ( CellDescriptor const & cell, localIndex const i )
    {
      globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        elementLocalDofIndexRow[i * NDOF + idof] = offset + idof;
      }
    } );

    stencil.forAll( [&] ( CellDescriptor const & cell, real64 w, localIndex const i )
    {
      globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        elementLocalDofIndexCol[i * NDOF + idof] = offset + idof;
      }
    } );

    sparsity->InsertGlobalIndices( integer_conversion<int>(numElems * m_numDofPerCell),
                                   elementLocalDofIndexRow.data(),
                                   integer_conversion<int>(stencilSize * m_numDofPerCell),
                                   elementLocalDofIndexCol.data() );
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
  fluxApprox->forBoundaryStencils( [&] ( FluxApproximationBase::BoundaryStencil const & boundaryStencilCollection )
  {
    boundaryStencilCollection.forAll( GEOSX_LAMBDA ( StencilCollection<PointDescriptor, real64>::Accessor stencil )
    {
      localIndex const stencilSize = stencil.size();
      stackArray1d<globalIndex, numElems   * maxNumDof> elementLocalDofIndexRow( numElems * NDOF );
      stackArray1d<globalIndex, maxStencil * maxNumDof> elementLocalDofIndexCol( stencilSize * NDOF );

      stencil.forConnected( [&] ( PointDescriptor const & point, localIndex const i )
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & cell = point.cellIndex;
          globalIndex const offset = NDOF * dofNumber[cell.region][cell.subRegion][cell.index];
          for (localIndex idof = 0; idof < NDOF; ++idof)
          {
            elementLocalDofIndexRow[idof] = offset + idof;
          }
        }
      });

      integer counter = 0;
      stencil.forAll( [&] ( PointDescriptor const & point, real64 w, localIndex i )
      {
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
      });

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

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
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
      if (elemGhostRank[ei] >= 0)
        return;

      stackArray1d<long long, maxNumDof>           localAccumDOF( NDOF );
      stackArray1d<double, maxNumComp>             localAccum( NC );
      stackArray2d<double, maxNumComp * maxNumDof> localAccumJacobian( NC, NDOF );

      // temporary work arrays
      stackArray1d<double, maxNumComp> dPhaseAmount_dC( NC );
      stackArray1d<double, maxNumComp> dPhaseCompFrac_dC( NC );

      // reset the local values
      localAccum = 0.0;
      localAccumJacobian = 0.0;

      // set DOF indices for this block
      globalIndex const offset = NDOF * dofNumber[ei];
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        localAccumDOF[idof] = integer_conversion<long long>(offset + idof);
      }

      // compute fluid-independent (pore volume) part
      real64 const volNew   = volume[ei];
      real64 const volOld   = volume[ei];
      real64 const dVol_dP  = 0.0; // used in poroelastic solver

      real64 const poroNew  = porosityRef[ei] * pvMult[ei][0];
      real64 const poroOld  = porosityOld[ei];
      real64 const dPoro_dP = porosityRef[ei] * dPvMult_dPres[ei][0];

      real64 const poreVolNew = volNew * poroNew;
      real64 const poreVolOld = volOld * poroOld;
      real64 const dPoreVol_dP = dVol_dP * poroNew + volNew * dPoro_dP;

      // sum contributions to component accumulation from each phase
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip];
        real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ei][ip] * phaseDensOld[ei][ip];

        real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip]
                                      + poreVolNew * (dPhaseVolFrac_dPres[ei][ip] * phaseDens[ei][0][ip]
                                                           + phaseVolFrac[ei][ip] * dPhaseDens_dPres[ei][0][ip]);

        // assemble density dependence
        applyChainRule( NC, dCompFrac_dCompDens[ei], dPhaseDens_dComp[ei][0][ip], dPhaseAmount_dC );
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc]     * phaseVolFrac[ei][ip]
                              + phaseDens[ei][0][ip] * dPhaseVolFrac_dCompDens[ei][ip][jc];
          dPhaseAmount_dC[jc] *= poreVolNew;
        }

        // ic - index of component whose conservation equation is assembled
        // (i.e. row number in local matrix)
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ei][0][ip][ic];
          real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ei][ip][ic];

          real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ei][0][ip][ic]
                                           + phaseAmountNew  * dPhaseCompFrac_dPres[ei][0][ip][ic];

          localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
          localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

          // jc - index of component w.r.t. whose compositional var the derivative is being taken
          // (i.e. col number in local matrix)

          // assemble phase composition dependence
          applyChainRule( NC, dCompFrac_dCompDens[ei], dPhaseCompFrac_dComp[ei][0][ip][ic], dPhaseCompFrac_dC );
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            real64 const dPhaseCompAmount_dC =           dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + phaseCompFrac[ei][0][ip][ic] * dPhaseAmount_dC[jc];
            localAccumJacobian[ic][jc+1] += dPhaseCompAmount_dC;
          }
        }
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
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & blockLocalDofNumber = m_dofNumber;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & pres = m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dPres = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & gravDepth = m_gravDepth;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & phaseMob = m_phaseMob;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dPhaseMob_dPres = m_dPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dPhaseMob_dComp = m_dPhaseMob_dCompDens;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dPhaseVolFrac_dPres = m_dPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dPhaseVolFrac_dComp = m_dPhaseVolFrac_dCompDens;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dCompFrac_dCompDens = m_dCompFrac_dCompDens;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & phaseDens        = m_phaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dPhaseDens_dPres = m_dPhaseDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseDens_dComp = m_dPhaseDens_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & phaseCompFrac    = m_phaseCompFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseCompFrac_dPres = m_dPhaseCompFrac_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView5d<real64>> const & dPhaseCompFrac_dComp = m_dPhaseCompFrac_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & phaseCapPressure                = m_phaseCapPressure;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseCapPressure_dPhaseVolFrac = m_dPhaseCapPressure_dPhaseVolFrac;

  localIndex constexpr numElems   = StencilCollection<CellDescriptor, real64>::NUM_POINT_IN_FLUX;
  localIndex constexpr maxStencil = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  real64 constexpr densWeight[numElems] = { 0.5, 0.5 };

  stencilCollection.forAll<RAJA::seq_exec>(GEOSX_LAMBDA (StencilCollection<CellDescriptor, real64>::Accessor stencil)
  {
    localIndex const stencilSize = stencil.size();

    // create local work arrays
    stackArray1d<long long, numElems * maxNumComp>  eqnRowIndices( numElems * NC );
    stackArray1d<long long, maxStencil * maxNumDof> dofColIndices( stencilSize * NDOF );

    stackArray1d<double, numElems * maxNumComp>                          localFlux( numElems * NC );
    stackArray2d<double, numElems * maxNumComp * maxStencil * maxNumDof> localFluxJacobian( numElems * NC, stencilSize * NDOF );

    stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

    stackArray1d<real64, maxStencil>              dPhaseFlux_dP( stencilSize );
    stackArray2d<real64, maxStencil * maxNumComp> dPhaseFlux_dC( stencilSize, NC );

    stackArray1d<real64, maxNumComp>                           compFlux( NC );
    stackArray2d<real64, maxStencil * maxNumComp>              dCompFlux_dP( stencilSize, NC );
    stackArray3d<real64, maxStencil * maxNumComp * maxNumComp> dCompFlux_dC( stencilSize, NC, NC );

    stackArray1d<real64, maxNumComp> dCapPressure_dC( NC );
    stackArray1d<real64, maxNumComp> dDens_dC( NC );

    stackArray1d<real64, numElems>              dDensMean_dP( numElems );
    stackArray2d<real64, numElems * maxNumComp> dDensMean_dC( numElems, NC );

    stackArray1d<real64, maxStencil>              dPresGrad_dP( stencilSize );
    stackArray2d<real64, maxStencil * maxNumComp> dPresGrad_dC( stencilSize, NC );
    
    stackArray1d<real64, numElems>                dGravHead_dP( numElems );
    stackArray2d<real64, numElems * maxNumComp>   dGravHead_dC( numElems, NC );

    // reset the local values
    compFlux = 0.0;
    dCompFlux_dP = 0.0;
    dCompFlux_dC = 0.0;

    localFlux = 0.0;
    localFluxJacobian = 0.0;

    // set equation indices for both connected cells
    stencil.forConnected( [&] ( auto const & cell,
                                localIndex i )
    {
      localIndex const er  = cell.region;
      localIndex const esr = cell.subRegion;
      localIndex const ei  = cell.index;

      globalIndex const offset = NDOF * blockLocalDofNumber[er][esr][ei];
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        eqnRowIndices[i * NC + ic] = offset + ic;
      }
    });

    // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // clear working arrays
      real64 densMean = 0.0;
      dDensMean_dP = 0.0;
      dDensMean_dC = 0.0;

      real64 presGrad = 0.0;
      dPresGrad_dP = 0.0;
      dPresGrad_dC = 0.0;

      real64 gravHead = 0.0;
      dGravHead_dP = 0.0;
      dGravHead_dC = 0.0;

      real64 phaseFlux = 0.0;
      dPhaseFlux_dP = 0.0;
      dPhaseFlux_dC = 0.0;

      // calculate quantities on primary connected cells
      stencil.forConnected( [&] ( auto const & cell,
                                 localIndex i )
      {
        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        // density
        real64 const density  = phaseDens[er][esr][m_fluidIndex][ei][0][ip];
        real64 const dDens_dP = dPhaseDens_dPres[er][esr][m_fluidIndex][ei][0][ip];
        applyChainRule( NC, dCompFrac_dCompDens[er][esr][ei], dPhaseDens_dComp[er][esr][m_fluidIndex][ei][0][ip], dDens_dC );

        // average density and pressure derivative
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;

        // compositional derivatives
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dDensMean_dC[i][jc] = densWeight[i] * dDens_dC[jc];
        }
      });

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      stencil.forAll( [&] ( CellDescriptor cell,
                            real64 w,
                            localIndex i )
      {
        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        globalIndex const offset = NDOF * blockLocalDofNumber[er][esr][ei];
        for (localIndex jdof = 0; jdof < NDOF; ++jdof)
        {
          dofColIndices[i * NDOF + jdof] = offset + jdof;
        }

        //capillary pressure
	real64 capPressure     = 0.0;
	real64 dCapPressure_dP = 0.0;
	dCapPressure_dC = 0.0;

	if (m_capPressureFlag)
	{
          capPressure = phaseCapPressure[er][esr][m_capPressureIndex][ei][0][ip];

          for (localIndex jp = 0; jp < NP; ++jp)
          {
            real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][m_capPressureIndex][ei][0][ip][jp];
            dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];

            for (localIndex jc = 0; jc < NC; ++jc)
            {
              dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
            }
          }
	}

        presGrad += w * (pres[er][esr][ei] + dPres[er][esr][ei] - capPressure);
        dPresGrad_dP[i] += w * (1 - dCapPressure_dP);
        for (localIndex jc = 0; jc < NC; ++jc)
        {
	  dPresGrad_dC[i][jc] += - w * dCapPressure_dC[jc];
	}

        if (m_gravityFlag)
        {
          real64 const gravD = w * gravDepth[er][esr][ei];
          gravHead += densMean * gravD;

          // need to add contributions from both cells the mean density depends on
          stencil.forConnected( [&] ( CellDescriptor,
                                      localIndex j )
          {
            dGravHead_dP[j] += dDensMean_dP[j] * gravD;
            for (localIndex jc = 0; jc < NC; ++jc)
            {
              dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
            }
          });
        }
      });

      // *** upwinding ***

      // use PPU currently; advanced stuff like IHU would go here
      // TODO isolate into a kernel?

      // compute phase potential gradient
      real64 potGrad = presGrad + gravHead;

      // choose upstream cell
      localIndex const k_up = (potGrad >= 0) ? 0 : 1;

      CellDescriptor cell_up = stencil.connectedIndex( k_up );
      localIndex er_up  = cell_up.region;
      localIndex esr_up = cell_up.subRegion;
      localIndex ei_up  = cell_up.index;

      real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

      // skip the phase flux if phase not present or immobile upstream
      if (std::fabs(mobility) < 1e-20) // TODO better constant
      {
        continue;
      }

      // pressure gradient depends on all points in the stencil
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
        }

      }

      // gravitational head depends only on the two cells connected (same as mean density)
      for (localIndex ke = 0; ke < numElems; ++ke)
      {
        dPhaseFlux_dP[ke] += dGravHead_dP[ke];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] += dGravHead_dC[ke][jc];
        }
      }

      // compute the phase flux and derivatives using upstream cell mobility
      phaseFlux = mobility * potGrad;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] *= mobility;
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] *= mobility;
        }
      }

      real64 const dMob_dP  = dPhaseMob_dPres[er_up][esr_up][ei_up][ip];
      arraySlice1d<real64 const> dPhaseMob_dCompSub = dPhaseMob_dComp[er_up][esr_up][ei_up][ip];

      // add contribution from upstream cell mobility derivatives
      dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseFlux_dC[k_up][jc] += dPhaseMob_dCompSub[jc] * potGrad;
      }

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d<real64 const> phaseCompFracSub = phaseCompFrac[er_up][esr_up][m_fluidIndex][ei_up][0][ip];
      arraySlice1d<real64 const> dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up][esr_up][m_fluidIndex][ei_up][0][ip];
      arraySlice2d<real64 const> dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up][esr_up][m_fluidIndex][ei_up][0][ip];

      // compute component fluxes and derivatives using upstream cell composition
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const ycp = phaseCompFracSub[ic];
        compFlux[ic] += phaseFlux * ycp;

        // derivatives stemming from phase flux
        for (localIndex ke = 0; ke < stencilSize; ++ke)
        {
          dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

        // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component densities
        applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dPhaseCompFrac_dCompDens );
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompFlux_dC[k_up][ic][jc] += phaseFlux * dPhaseCompFrac_dCompDens[jc];
        }
      }
    }

    // *** end of upwinding

    // populate local flux vector and derivatives
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      localFlux[ic]      =  dt * compFlux[ic];
      localFlux[NC + ic] = -dt * compFlux[ic];

      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        localIndex const localDofIndexPres = ke * NDOF;
        localFluxJacobian[ic][localDofIndexPres] = dt * dCompFlux_dP[ke][ic];
        localFluxJacobian[NC + ic][localDofIndexPres] = -dt * dCompFlux_dP[ke][ic];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
          localFluxJacobian[ic][localDofIndexComp] = dt * dCompFlux_dC[ke][ic][jc];
          localFluxJacobian[NC + ic][localDofIndexComp] = -dt * dCompFlux_dC[ke][ic][jc];
        }
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

  });
}

void CompositionalMultiphaseFlow::AssembleVolumeBalanceTerms( DomainPartition const * const domain,
                                                              Epetra_FECrsMatrix * const jacobian,
                                                              Epetra_FEVector * const residual,
                                                              real64 const time_n,
                                                              real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC   = m_numComponents;
  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
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
        return;

      stackArray1d<long long, maxNumDof> localVolBalanceDOF( NDOF );
      stackArray1d<double, maxNumDof>    localVolBalanceJacobian( NDOF );

      // compute pore volume
      real64 const vol      = volume[ei];
      real64 const dVol_dP  = 0.0; // used in poroelastic solver

      real64 const poro     = porosityRef[ei] * pvMult[ei][0];
      real64 const dPoro_dP = porosityRef[ei] * dPvMult_dPres[ei][0];

      real64 const poreVol     = vol * poro;
      real64 const dPoreVol_dP = dVol_dP * poro + vol * dPoro_dP;

      // get equation/dof indices
      globalIndex const offset = NDOF * dofNumber[ei];
      globalIndex const localVolBalanceEqnIndex = offset + NC;
      for (localIndex jdof = 0; jdof < NDOF; ++jdof)
      {
        localVolBalanceDOF[jdof] = offset + jdof;
      }

      real64 localVolBalance = 1.0;
      localVolBalanceJacobian = 0.0;

      // sum contributions to component accumulation from each phase
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        localVolBalance -= phaseVolFrac[ei][ip];
        localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ei][ip];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ei][ip][jc];
        }
      }

      // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
      for (localIndex idof = 0; idof < NDOF; ++idof)
      {
        localVolBalanceJacobian[idof] *= poreVol;
      }
      localVolBalanceJacobian[0] += dPoreVol_dP * localVolBalance;
      localVolBalance *= poreVol;

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
    });
  });
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

  unordered_map<string, array1d<bool>> bcStatusMap; // map to check consistent application of BC

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

    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference<array1d<globalIndex>>( viewKeyStruct::blockLocalDofNumberString );

    arrayView1d<real64 const> const & pres      = subRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );
    arrayView1d<real64 const> const & dPres     = subRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
    arrayView1d<real64 const> const & bcPres    = subRegion->getReference<array1d<real64>>( viewKeyStruct::bcPressureString );
    arrayView2d<real64 const> const & compFrac  = subRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );
    arrayView2d<real64 const> const & compDens  = subRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );
    arrayView2d<real64 const> const & dCompDens = subRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d<real64> const & totalDens = fluid->getReference<array2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString );

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

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  real64 localResidualNorm = 0.0;
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
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

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  RAJA::ReduceMin<reducePolicy, integer> result(1);

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
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

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
                                   ElementSubRegionBase const * const subRegion )
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

  applyToSubRegions( domain, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );
}

void CompositionalMultiphaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * elementRegion,
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

  applyToSubRegions( domain, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );
}

void CompositionalMultiphaseFlow::ImplicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * elementRegion,
                                   ElementSubRegionBase const * const subRegion )
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
    elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( viewKeyStruct::blockLocalDofNumberString );
  m_pressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::deltaPressureString );
  m_globalCompDensity =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::globalCompDensityString );
  m_deltaGlobalCompDensity =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  m_compFrac =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::globalCompFractionString );
  m_dCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
  m_phaseVolFrac =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::phaseVolumeFractionString );
  m_dPhaseVolFrac_dPres =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  m_dPhaseVolFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );
  m_phaseMob =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::phaseMobilityString );
  m_dPhaseMob_dPres =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::dPhaseMobility_dPressureString );
  m_dPhaseMob_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  m_porosityOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::porosityOldString );
  m_phaseVolFracOld =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::phaseVolumeFractionOldString );
  m_phaseDensOld =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::phaseDensityOldString );
  m_phaseCompFracOld =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( viewKeyStruct::phaseComponentFractionOldString );

  m_pvMult =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                      constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                      constitutiveManager );
  m_phaseFrac =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString,
                                                                                      constitutiveManager );
  m_dPhaseFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseDens =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                      constitutiveManager );
  m_dPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseDens_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseVisc =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString,
                                                                                      constitutiveManager );
  m_dPhaseVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseVisc_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_phaseCompFrac =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString,
                                                                                      constitutiveManager );
  m_dPhaseCompFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dPhaseCompFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array5d<real64>, arrayView5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_totalDens =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString,
                                                                                      constitutiveManager );
  m_phaseRelPerm =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString,
                                                                                      constitutiveManager );
  m_dPhaseRelPerm_dPhaseVolFrac =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                      constitutiveManager );
  if (m_capPressureFlag)
  {
    m_phaseCapPressure =
      elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( CapillaryPressureBase::viewKeyStruct::phaseCapPressureString,
                                                                                        constitutiveManager );
    m_dPhaseCapPressure_dPhaseVolFrac =
      elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( CapillaryPressureBase::viewKeyStruct::dPhaseCapPressure_dPhaseVolFractionString,
                                                                                        constitutiveManager );
  }
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseFlow, string const &, ManagedGroup * const)
}// namespace geosx
