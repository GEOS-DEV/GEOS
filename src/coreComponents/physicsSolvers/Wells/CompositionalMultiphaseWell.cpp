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
 * @file CompositionalMultiphaseWell.cpp
 */

#include "CompositionalMultiphaseWell.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/FiniteVolume/CompositionalMultiphaseFlow.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "wells/WellManager.hpp"
#include "wells/Well.hpp"
#include "wells/PerforationData.hpp"
#include "wells/Perforation.hpp"
#include "wells/ConnectionData.hpp"
#include "wells/Connection.hpp"
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

CompositionalMultiphaseWell::CompositionalMultiphaseWell( const string & name,
                                                                      ManagedGroup * const parent )
  :
  WellSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 )
{
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
  
}
  
localIndex CompositionalMultiphaseWell::numFluidComponents() const
{
  return m_numComponents;
}

localIndex CompositionalMultiphaseWell::numFluidPhases() const
{
  return m_numPhases;
}

void CompositionalMultiphaseWell::RegisterDataOnMesh(ManagedGroup * const meshBodies)
{
  WellSolverBase::RegisterDataOnMesh(meshBodies);

  WellManager * const wellManager = meshBodies->getParent()->group_cast<DomainPartition *>()->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {

    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dofNumberString ); 
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalCompDensityString );
    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalCompFractionString );
    wellElementSubRegion->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );
    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
    wellElementSubRegion->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::sumCompPerforationRateString );
    
    ConnectionData * const connectionData = well->getConnections();
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::mixtureRateString );
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaMixtureRateString );

    PerforationData * const perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::compPerforationRateString );
    perforationData->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString );
    perforationData->RegisterViewWrapper<array4d<real64>>( viewKeyStruct::dCompPerforationRate_dCompString );
   
  });    
}
  
void CompositionalMultiphaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );
  
  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    MultiFluidBase const * const fluid = GetFluidModel( wellElementSubRegion );
    // TODO: put this out of the loop
    m_numPhases     = fluid->numFluidPhases();
    m_numComponents = fluid->numFluidComponents();

    localIndex const NC = m_numComponents;
    localIndex const NP = m_numPhases;
    
    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompDensityString).resizeDimension<1>(NC);
    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::deltaGlobalCompDensityString).resizeDimension<1>(NC);

    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompFractionString).resizeDimension<1>(NC);
    wellElementSubRegion->getReference<array3d<real64>>(viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString).resizeDimension<1,2>(NC, NC);

    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::phaseVolumeFractionString).resizeDimension<1>(NP);
    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dPressureString).resizeDimension<1>(NP);
    wellElementSubRegion->getReference<array3d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::sumCompPerforationRateString).resizeDimension<1>(NC);
    
    PerforationData * const perforationData = well->getPerforations();
    perforationData->getReference<array3d<real64>>( viewKeyStruct::compPerforationRateString ).resizeDimension<1,2>( NP, NC );
    perforationData->getReference<array3d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString ).resizeDimension<1,2>( NP, NC );
    perforationData->getReference<array4d<real64>>( viewKeyStruct::dCompPerforationRate_dCompString ).resizeDimension<1,2,3>( NP, NC, NC );
  });
  
}

MultiFluidBase * CompositionalMultiphaseWell::GetFluidModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels =
    dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  MultiFluidBase * const fluid = constitutiveModels->GetGroup<MultiFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

MultiFluidBase const * CompositionalMultiphaseWell::GetFluidModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels =
    dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  MultiFluidBase const * const fluid = constitutiveModels->GetGroup<MultiFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

RelativePermeabilityBase * CompositionalMultiphaseWell::GetRelPermModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels =
    dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  RelativePermeabilityBase * const relPerm = constitutiveModels->GetGroup<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Target group does not contain the relative permeability model" );

  return relPerm;
}

RelativePermeabilityBase const * CompositionalMultiphaseWell::GetRelPermModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels =
    dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  RelativePermeabilityBase const * const relPerm = constitutiveModels->GetGroup<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Target group does not contain the relative permeability model" );

  return relPerm;
}

void CompositionalMultiphaseWell::UpdateComponentFraction( Well * well )
{
  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

  arrayView2d<real64 const> const & wellElemCompDens =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dWellElemCompDens =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  arrayView2d<real64> const & wellElemCompFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  arrayView3d<real64> const & dWellElemCompFrac_dCompDens =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
  {
    real64 wellElemTotalDensity = 0.0;

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      wellElemTotalDensity += wellElemCompDens[iwelem][ic]
                           + dWellElemCompDens[iwelem][ic];
    }

    real64 const wellElemTotalDensityInv = 1.0 / wellElemTotalDensity;

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      wellElemCompFrac[iwelem][ic] = (wellElemCompDens[iwelem][ic]
				   + dWellElemCompDens[iwelem][ic]) * wellElemTotalDensityInv;

      for (localIndex jc = 0; jc < m_numComponents; ++jc)
      {
        dWellElemCompFrac_dCompDens[iwelem][ic][jc] = - wellElemCompFrac[iwelem][ic] * wellElemTotalDensityInv;
      }
      dWellElemCompFrac_dCompDens[iwelem][ic][ic] += wellElemTotalDensityInv;
    }
  }
}

void CompositionalMultiphaseWell::UpdatePhaseVolumeFraction( Well * well )
{
  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
  MultiFluidBase const * const fluid = GetFluidModel( wellElementSubRegion );

  arrayView2d<real64> const & wellElemPhaseVolFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

  arrayView2d<real64> const & dWellElemPhaseVolFrac_dPres =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64> const & dWellElemPhaseVolFrac_dComp =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d<real64> const & dWellElemCompFrac_dCompDens =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  arrayView2d<real64 const> const & wellElemCompDens =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

  arrayView2d<real64 const> const & dWellElemCompDens =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

  arrayView3d<real64 const> const & wellElemPhaseFrac =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString );

  arrayView3d<real64 const> const & dWellElemPhaseFrac_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

  arrayView4d<real64 const> const & dWellElemPhaseFrac_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

  arrayView3d<real64 const> const & wellElemPhaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dWellElemPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d<real64 const> const & dWellElemPhaseDens_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
  {
    stackArray1d<real64, maxNumComp> work( NC );

    // compute total density from component partial densities
    real64 wellElemTotalDensity = 0.0;
    real64 const dWellElemTotalDens_dCompDens = 1.0;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      wellElemTotalDensity += wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic];
    }

    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
      real64 const wellElemPhaseDensInv = 1.0 / wellElemPhaseDens[iwelem][0][ip];

      // compute saturation and derivatives except multiplying by the total density
      wellElemPhaseVolFrac[iwelem][ip] = wellElemPhaseFrac[iwelem][0][ip] * wellElemPhaseDensInv;

      dWellElemPhaseVolFrac_dPres[iwelem][ip] =
        (dWellElemPhaseFrac_dPres[iwelem][0][ip] - wellElemPhaseVolFrac[iwelem][ip] * dWellElemPhaseDens_dPres[iwelem][0][ip])
	* wellElemPhaseDensInv;

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dWellElemPhaseVolFrac_dComp[iwelem][ip][jc] =
          (dWellElemPhaseFrac_dComp[iwelem][0][ip][jc] - wellElemPhaseVolFrac[iwelem][ip] * dWellElemPhaseDens_dComp[iwelem][0][ip][jc])
	  * wellElemPhaseDensInv;
      }

      // apply chain rule to convert derivatives from global component fractions to densities
      applyChainRuleInPlace( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseVolFrac_dComp[iwelem][ip], work );

      // now finalize the computation by multiplying by total density
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dWellElemPhaseVolFrac_dComp[iwelem][ip][jc] *= wellElemTotalDensity;
        dWellElemPhaseVolFrac_dComp[iwelem][ip][jc] += wellElemPhaseVolFrac[iwelem][ip] * dWellElemTotalDens_dCompDens;
      }

      wellElemPhaseVolFrac[iwelem][ip] *= wellElemTotalDensity;
      dWellElemPhaseVolFrac_dPres[iwelem][ip] *= wellElemTotalDensity;
    }
  }
}

void CompositionalMultiphaseWell::UpdateFluidModel( Well * well )
{
  WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
  MultiFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

  arrayView1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & dWellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  arrayView2d<real64 const> const & wellElemCompFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
  {
    fluid->PointUpdate( wellElemPressure[iwelem] + dWellElemPressure[iwelem],
			m_temperature,
			wellElemCompFrac[iwelem],
			iwelem,
			0 );
  }
}
  
void CompositionalMultiphaseWell::UpdateRelPermModel( Well * well )
{
  WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
  RelativePermeabilityBase * const relPerm = GetRelPermModel( wellElementSubRegion );

  arrayView2d<real64> const & wellElemPhaseVolFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

  relPerm->BatchUpdate( wellElemPhaseVolFrac );
}


void CompositionalMultiphaseWell::UpdateState( Well * well )
{
  UpdateComponentFraction( well );
  UpdateFluidModel( well );
  UpdatePhaseVolumeFraction( well );
  UpdateRelPermModel( well );  
}

void CompositionalMultiphaseWell::UpdateStateAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    UpdateState( well );
  });
}

void CompositionalMultiphaseWell::InitializeWells( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & resPressure = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & resCompDens = m_resGlobalCompDensity;
  
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseDens     = m_resPhaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseVisc     = m_resPhaseVisc;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & resPhaseCompFrac = m_resPhaseCompFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseRelPerm  = m_resPhaseRelPerm;

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    ConnectionData const * const connectionData = well->getConnections();
    PerforationData const * const perforationData = well->getPerforations();

    // get well primary variables on segments
    array1d<real64> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array2d<real64> const & wellElemCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );
    
    // get the info stored on segments
    array2d<real64> const & wellElemCompFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

    array2d<real64> const & wellElemSumCompPerfRates =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::sumCompPerforationRateString );

    array1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );
    
    // get a reference to the primary variables on connections
    array1d<real64> const & connRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::mixtureRateString );

    // get the well element indices corresponding to each connect
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    arrayView1d<localIndex const> const & prevWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::prevWellElementIndexString );    
    
    // get well data on perforations
    array2d<real64> const & compPerfRate =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString );

    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );    
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // TODO: fix this in the case of multiple MPI ranks working on the same well
    // In particular, steps 1), 3) and 7) will need to be rewritten
    
    // 1) Loop over all perforations to compute an average mixture density and component fraction
    real64 avgTotalDensity   = 0.0;
    real64 avgMixtureDensity = 0.0;
    stackArray1d<real64, maxNumComp> avgCompFrac;
    
    avgCompFrac = 0.0;
    
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // compute the total mobility
      real64 resTotalMobility = 0.;
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        real64 const resViscosity = resPhaseVisc[er][esr][m_fluidIndex][ei][0][ip];
        real64 const resRelPerm   = resPhaseRelPerm[er][esr][m_relPermIndex][ei][0][ip];
	
	resTotalMobility += resRelPerm / resViscosity;
      }
      
      // increment the average mixture density
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        real64 const resDensity   = resPhaseDens[er][esr][m_fluidIndex][ei][0][ip];
        real64 const resViscosity = resPhaseVisc[er][esr][m_fluidIndex][ei][0][ip];
        real64 const resRelPerm   = resPhaseRelPerm[er][esr][m_relPermIndex][ei][0][ip];
	
        real64 const resMobility = resRelPerm / resViscosity;
	avgMixtureDensity += resMobility * resDensity / resTotalMobility;
      }


      // increment the average global component fraction
      real64 perfTotalDensity = 0.0;
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        perfTotalDensity += resCompDens[er][esr][ei][ic];
      }

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        avgCompFrac[ic] += resCompDens[er][esr][ei][ic] / perfTotalDensity;
      }
      avgTotalDensity += perfTotalDensity;
    }
    
    globalIndex const numPerforationsGlobal = perforationData->numPerforationsGlobal();

    // compute average densities
    avgMixtureDensity /= numPerforationsGlobal;
    avgTotalDensity   /= numPerforationsGlobal;

    // compute average component fraction
    if ( well->getType() == Well::Type::PRODUCER )
    {
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        avgCompFrac[ic] /= numPerforationsGlobal; // use average comp frac from reservoir
      }  
    }
    else // injector
    {
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        avgCompFrac[ic] = well->getInjectionStream( ic ); // use average comp frac from XML file
      }
    }

    // set the global component fractions to avgCompFrac
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); iwelem++)
    {
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        wellElemCompFrac[iwelem][ic] = avgCompFrac[ic];
	wellElemCompDens[iwelem][ic] = avgCompFrac[ic] * avgTotalDensity;
      }
    }

    // get the reference data for this well
    localIndex const iwelemRef = 0; // here, we assume that the first segment is on this MPI rank, which will NOT be the case in general
    real64 const refGravDepth = well->getReferenceDepth();
    
    // 2) Initialize the reference pressure
    real64 const targetBHP = well->getTargetBHP();
    if (well->getControl() == Well::Control::BHP)
    {
      // if pressure constraint, set the ref pressure at the constraint
      wellElemPressure[iwelemRef] = targetBHP;
    }
    else // rate control
    {
      real64 resPres = well->getTargetBHP();
      if (perforationData->numPerforationsLocal() > 0)
      {
        // TODO: again, this will need to be improved
        // we have to make sure that iperf = 0 is the first perforation
        localIndex const er  = resElementRegion[0];
        localIndex const esr = resElementSubRegion[0];
        localIndex const ei  = resElementIndex[0];
        resPres = resPressure[er][esr][ei];
      }
      // if rate constraint, set the ref pressure slightly above/below the target pressure
      wellElemPressure[iwelemRef] = (well->getType() == Well::Type::PRODUCER)
	                          ? 0.9 * resPres
		    	          : 1.1 * resPres;
    }

    // TODO: communicate ref pressure
    // TODO: communicate avgDens
    
    // 3) Estimate the pressures in the segments using this avgDensity
    // TODO: implement this in parallel with the communication of the pressures
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (iwelem != iwelemRef)
        wellElemPressure[iwelem] = wellElemPressure[iwelemRef]
 	  + ( m_gravityFlag ? avgMixtureDensity * ( wellElemGravDepth[iwelem] - refGravDepth ) : 0 );
    }

    // 4) Recompute the pressure-dependent properties
    UpdateState( well );

    // 5) Compute the perforation rates
    ComputeAllPerforationRates( well );

    // 6) Collect all the perforation rates
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      for (localIndex ic = 0; ic < NC; ++ic)
        wellElemSumCompPerfRates[iwelem][ic] = 0.0; // TODO: ask for a better way to do that!
    }
    
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {  
      localIndex const iwelem = perfWellElemIndex[iperf];
      for (localIndex ic = 0; ic < NC; ++ic)
        wellElemSumCompPerfRates[iwelem][ic] += compPerfRate[iperf][ic];
    }

    // 7) Estimate the connection rates
    // TODO: implement this in parallel with the communication of wellElemSumRates
    real64 prevConnRate = 0;
    localIndex const lastConnectionIndex = connectionData->numConnectionsLocal()-1;
    for (localIndex iconn = lastConnectionIndex; iconn >= 0; --iconn)
    {
      localIndex const iwelemNext = nextWellElemIndex[iconn];

      connRate[iconn] = prevConnRate;
      for (localIndex ic = 0; ic < m_numComponents; ic++)
        connRate[iconn] -= wellElemSumCompPerfRates[iwelemNext][ic] / wellElemCompDens[iwelemNext][ic];
      
      prevConnRate = connRate[iconn];
    }
  });
}

void CompositionalMultiphaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );

  // Initialize the primary and secondary variables
  InitializeWells( domain );
}

void CompositionalMultiphaseWell::BackupFields( DomainPartition * const domain )
{
  // TODO: backup a normalizer
}

void
CompositionalMultiphaseWell::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                DomainPartition * const domain,
                                                EpetraBlockSystem * const blockSystem )
{
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain );

  // assumes that the setup of dof numbers and linear system
  // is done in ReservoirWellSolver

}

void CompositionalMultiphaseWell::SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
                                                                localIndex & numLocalRows,
                                                                globalIndex & numGlobalRows,
                                                                localIndex offset )
{
  std::cout << "CompositionalMultiphaseWell::SetNumRowsAndTrilinosIndices started" << std::endl;
  
  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather(numMpiProcesses);

  // communicate the number of local rows to each process
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                &gather.front(),
                                                1 );

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<thisMpiProcess)
      firstLocalRow += gather[p];
  }

  // TODO: double check this for multiple MPI processes
  
  // get the well information
  WellManager const * const wellManager = domain->getWellManager();

  localIndex localCount = 0;
  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    
    arrayView1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
  
    // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
    for ( localIndex iwelem = 0; iwelem < wellElemGhostRank.size(); ++iwelem )
    {
      wellElemDofNumber[iwelem] = -1;
    }

    // loop over all well elements and set the dof number if the element is not a ghost
    for ( localIndex iwelem = 0; iwelem < wellElemGhostRank.size(); ++iwelem )
    {
      if (wellElemGhostRank[iwelem] < 0 )
      {
        wellElemDofNumber[iwelem] = firstLocalRow + localCount + offset;
        localCount += 1;
      }
      else
      {
        wellElemDofNumber[iwelem] = -1;
      }
    }
  });
  	    
  GEOS_ERROR_IF(localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}

void CompositionalMultiphaseWell::SetSparsityPattern( DomainPartition const * const domain,
                                                      Epetra_FECrsGraph * const sparsity,
						      globalIndex firstWellElemDofNumber,
						      localIndex numDofPerResElement )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();
  WellManager const * const wellManager = domain->getWellManager();
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber =
    elementRegionManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( CompositionalMultiphaseFlow::viewKeyStruct::blockLocalDofNumberString );

  // save these numbers (reused to compute the well elem offsets in multiple functions)
  m_firstWellElemDofNumber = firstWellElemDofNumber;
  m_numDofPerResElement    = numDofPerResElement;

  // set the number of degrees of freedom per element on both sides (reservoir and well)  
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS; 
  localIndex constexpr maxNumDof  = maxNumComp + 2; // dofs are 1 pressure, NC comp densities (reservoir and well), 1 connection rate (well only)

  // reservoir dofs are 1 pressure and NC comp densities
  localIndex const resNDOF  = numDofPerResElement;
  // well dofs are 1 pressure, NC comp densities, 1 connection rate
  localIndex const wellNDOF = numDofPerElement()
                            + numDofPerConnection(); 

  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    PerforationData const * const perforationData = well->getPerforations();
    ConnectionData const * const connectionData   = well->getConnections();
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the well degrees of freedom 
    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get the well element indices corresponding to each connect
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    arrayView1d<localIndex const> const & prevWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::prevWellElementIndexString );    
    
    // get the well element indices corresponding to each perforation
    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );    
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // 1) Insert the entries corresponding to reservoir-well perforations
    //    This will fill J_WW, J_WR, and J_RW
    //    That is all we need for single-segmented wells
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( resNDOF + wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( resNDOF + wellNDOF );

      // get the offset of the reservoir element equation
      globalIndex const resOffset  = resNDOF * resDofNumber[er][esr][ei];
      for (localIndex idof = 0; idof < resNDOF; ++idof)
      {
	// specify the reservoir equation number
        elementLocalDofIndexRow[SubRegionTag::RES * wellNDOF + idof] = resOffset + idof;
	// specify the reservoir variable number
	elementLocalDofIndexCol[SubRegionTag::RES * wellNDOF + idof] = resOffset + idof;
      }

      /*
       * firstWellElemDofNumber denotes the first DOF number of the well segments, for all the wells (i.e., first segment of first well)
       * currentElemDofNumber denotes the DOF number of the current segment
       *
       * The coordinates of this element's (NC+2)x(NC+2) block in J_WW can be accessed using:
       *
       * IndexRow = firstWellElemDofNumber * resNDOF ( = all the equations in J_RR)
       *          + (currentElemDofNumber - firstWellElemDofNumber ) * wellNDOF ( = offset of current segment in J_WW)
       *          + idof ( = local dofs for this segment, pressure and velocity)
       *	   
       * This is needed because resNDOF is not equal to wellNDOF
       */

      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );
      
      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
	// specify the well equation number
	elementLocalDofIndexRow[SubRegionTag::WELL * resNDOF + idof] = elemOffset + idof;
	// specify the reservoir variable number
	elementLocalDofIndexCol[SubRegionTag::WELL * resNDOF + idof] = elemOffset + idof;
      }      

      sparsity->InsertGlobalIndices( integer_conversion<int>( resNDOF + wellNDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( resNDOF + wellNDOF ),
                                     elementLocalDofIndexCol.data() );
    }

    // 2) Insert the entries corresponding to well-well connection between segments
    //    This will fill J_WW only
    //    That is not needed for single-segmented wells
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {

      // get previous segment index
      localIndex const iwelemPrev =  prevWellElemIndex[iconn];

      // get second segment index
      localIndex const iwelemNext = nextWellElemIndex[iconn];

      // check if this is not an entry or exit
      if (iwelemPrev < 0 || iwelemNext < 0)
	continue;
      
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( 2 * wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( 2 * wellNDOF );

      // get the offset of the well element equations
      globalIndex const prevElemOffset = getElementOffset( wellElemDofNumber[iwelemPrev] );
      globalIndex const nextElemOffset = getElementOffset( wellElemDofNumber[iwelemNext] );

      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
	elementLocalDofIndexRow[idof]            = prevElemOffset + idof;
	elementLocalDofIndexRow[wellNDOF + idof] = nextElemOffset + idof;
	elementLocalDofIndexCol[idof]            = prevElemOffset + idof;
	elementLocalDofIndexCol[wellNDOF + idof] = nextElemOffset + idof;
      }      

      sparsity->InsertGlobalIndices( integer_conversion<int>( 2 * wellNDOF ),
                                     elementLocalDofIndexRow.data(),
                                     integer_conversion<int>( 2 * wellNDOF ),
                                     elementLocalDofIndexCol.data() );
    }
    
  });
  
}


void CompositionalMultiphaseWell::AssembleSystem( DomainPartition * const domain,
                                                  EpetraBlockSystem * const blockSystem,
                                                  real64 const time_n, real64 const dt )
{ 
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  // first deal with the well control and switch if necessary
  CheckWellControlSwitch( domain );

  // then assemble the mass balance equations
  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleSourceTerms( domain, jacobian, residual, time_n, dt );

  // then assemble the volume balance equations
  AssembleVolumeBalanceTerms( domain, jacobian, residual, time_n, dt );

  // then assemble the connection equations
  FormMomentumEquations( domain, jacobian, residual );

  // finally assemble the well control equation
  FormControlEquation( domain, jacobian, residual );
  
  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After CompositionalMultiphaseWell::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

}

void CompositionalMultiphaseWell::AssembleAccumulationTerms( DomainPartition * const domain,
                                                             Epetra_FECrsMatrix * const jacobian,
                                                             Epetra_FEVector * const residual,
                                                             real64 const time_n,
                                                             real64 const dt )
{
  // will not be implemented
}

void CompositionalMultiphaseWell::AssembleFluxTerms( DomainPartition * const domain,
                                                     Epetra_FECrsMatrix * const jacobian,
                                                     Epetra_FEVector * const residual,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  // TODO
}

void CompositionalMultiphaseWell::AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                                              Epetra_FECrsMatrix * const jacobian,
                                                              Epetra_FEVector * const residual,
                                                              real64 const time_n,
                                                              real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    // loop over the segments
  });    
}


void CompositionalMultiphaseWell::AssembleSourceTerms( DomainPartition * const domain,
                                                       Epetra_FECrsMatrix * const jacobian,
                                                       Epetra_FEVector * const residual,
                                                       real64 const time_n,
                                                       real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber = m_resDofNumber;

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;
  
  localIndex const NC      = m_numComponents;
  localIndex const NP      = m_numPhases;
  localIndex const resNDOF = numDofPerResElement();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData const * const perforationData = well->getPerforations();
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // compute the local rates for this well
    ComputeAllPerforationRates( well );
    
    // get the degrees of freedom 
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // get well variables on perforations
    array2d<real64 const> const & compPerfRate =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString );

    array3d<real64 const> const & dCompPerfRate_dPres =
      perforationData->getReference<array3d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString );

    array4d<real64 const> const & dCompPerfRate_dComp =
      perforationData->getReference<array4d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString );
    
    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // loop over the perforations and add the rates to the residual and jacobian
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {

      // local working variables and arrays
      stackArray1d<long long, 2 * maxNumComp> eqnRowIndices( 2 * NC );
      stackArray1d<long long, 2 * maxNumDof>  dofColIndices( 2 * resNDOF );

      stackArray1d<double, 2 * maxNumComp>                 localFlux( 2 * NC );
      stackArray2d<double, 2 * maxNumComp * 2 * maxNumDof> localFluxJacobian( 2 * NC, 2 * resNDOF );

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const resOffset = resNDOF * resDofNumber[er][esr][ei];
      globalIndex const wellElemOffset = getElementOffset( wellElemDofNumber[iwelem] );

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        eqnRowIndices[SubRegionTag::RES  * NC + ic] = resOffset + ic;
        eqnRowIndices[SubRegionTag::WELL * NC + ic] = wellElemOffset + ic + RowOffset::MASSBAL;
      }
      for (localIndex jdof = 0; jdof < resNDOF; ++jdof)
      {
        dofColIndices[SubRegionTag::RES  * resNDOF + jdof] = resOffset  + jdof;
        dofColIndices[SubRegionTag::WELL * resNDOF + jdof] = wellElemOffset + jdof + ColOffset::DPRES;
      }
      
      // populate local flux vector and derivatives
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        localFlux[SubRegionTag::RES  * NC + ic] =   dt * compPerfRate[iperf][ic];
        localFlux[SubRegionTag::WELL * NC + ic] = - dt * compPerfRate[iperf][ic];

        for (localIndex ke = 0; ke < 2; ++ke)
        {
          localIndex const localDofIndexPres = ke * resNDOF;
          localFluxJacobian[SubRegionTag::RES  * NC + ic][localDofIndexPres] =   dt * dCompPerfRate_dPres[iperf][ke][ic];
          localFluxJacobian[SubRegionTag::WELL * NC + ic][localDofIndexPres] = - dt * dCompPerfRate_dPres[iperf][ke][ic];
	
  	  for (localIndex jc = 0; jc < NC; ++jc)
          {
	    localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
            localFluxJacobian[SubRegionTag::RES  * NC + ic][localDofIndexComp] =   dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
            localFluxJacobian[SubRegionTag::WELL * NC + ic][localDofIndexComp] = - dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
          }
	}
      }

      // Add to global residual/jacobian
      residual->SumIntoGlobalValues( integer_conversion<int>( 2 * NC ),
				     eqnRowIndices.data(),
				     localFlux.data() );

      jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 * NC ),
                                     eqnRowIndices.data(),
                                     integer_conversion<int>( 2 * resNDOF),
                                     dofColIndices.data(),
                                     localFluxJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR);
    }
    
  });  
}


real64
CompositionalMultiphaseWell::CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                                                    DomainPartition * const domain )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);
 
  WellManager * const wellManager = domain->getWellManager();

  real64 localResidualNorm = 0;
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers
    array1d<globalIndex const > const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

        localIndex const wellNDOF = numDofPerElement() + numDofPerConnection(); 
        for (localIndex idof = 0; idof < wellNDOF; ++idof)
        {
          int const lid = rowMap->LID(integer_conversion<int>( elemOffset + idof ));
          real64 const val = localResidual[lid];
          localResidualNorm += val * val;
        }
      }
    }
  });

  // compute global residual norm
  real64 globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

bool
CompositionalMultiphaseWell::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  return false;
}

void
CompositionalMultiphaseWell::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  // get the update
  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers on segments
    array1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // get a reference to the primary variables on segments
    array1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    // get a reference to the primary variables on connections
    array1d<real64> const & dConnRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureRateString );

    // get a reference to the next well element index for each connection
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    

    // get the ghosting information on segments and on connections
    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );


    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract solution and apply to dP
	globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

        int lid = rowMap->LID( integer_conversion<int>( elemOffset ) );
        dWellElemPressure[iwelem] += scalingFactor * local_solution[lid];

        for (localIndex ic = 0; ic < m_numComponents; ++ic)
        {
          lid = rowMap->LID( integer_conversion<int>( elemOffset + ic + 1 ) );
          dWellElemGlobalCompDensity[iwelem][ic] += scalingFactor * local_solution[lid];
        }
      }
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // as a temporary solution, the connection variable is looked up with the next element's DOF number
      localIndex iwelemNext = nextWellElemIndex[iconn];

      // check if there is a variable defined here
      if (iwelemNext < 0)
        continue;

      if (wellElemGhostRank[iwelemNext] < 0)
      {      
        // extract solution and apply to dQ
	globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelemNext] );
        int const lid = rowMap->LID( integer_conversion<int>( elemOffset + m_numComponents + 1 ) );
        dConnRate[iconn] += scalingFactor * local_solution[lid];
      }
    }
  });  

  // TODO: call CommunicationTools::SynchronizeFields

  // update properties
  UpdateStateAll( domain );
    
}

void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the ghosting information
    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );
    
    // get a reference to the primary variables on segments
    array1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    // get a reference to the primary variables on connections
    array1d<real64> const & dConnRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureRateString );

    // get a reference to the well element index for each connection
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    
    
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract solution and apply to dP
        dWellElemPressure[iwelem] = 0;
        for (localIndex ic = 0; ic < m_numComponents; ++ic)
          dWellElemGlobalCompDensity[iwelem][ic] = 0;
      }
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // as a temporary solution, the connection variable is looked up with the next element's DOF number
      localIndex iwelemNext = nextWellElemIndex[iconn];

      // check if there is a variable defined here
      if (iwelemNext < 0)
        continue;

      dConnRate[iconn] = 0;
    }
  });

  // call constitutive models
  UpdateStateAll( domain );
}

void CompositionalMultiphaseWell::ResetViews(DomainPartition * const domain)
{
  WellSolverBase::ResetViews(domain);

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_resDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( CompositionalMultiphaseFlow::viewKeyStruct::blockLocalDofNumberString );
  
  m_resPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

  m_deltaResPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

  m_resGlobalCompDensity =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString );

  m_deltaResGlobalCompDensity =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString );

  m_resCompFrac =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString );

  m_dResCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  m_resPhaseVolFrac =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::phaseVolumeFractionString );

  m_dResPhaseVolFrac_dPres =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  m_dResPhaseVolFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  m_resPhaseFrac =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString,
                                                                                      constitutiveManager );
  m_dResPhaseFrac_dPres =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dResPhaseFrac_dComp =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_resPhaseDens =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                      constitutiveManager );
  m_dResPhaseDens_dPres =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                      constitutiveManager );
  m_dResPhaseDens_dComp =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_resPhaseVisc =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString,
                                                                                      constitutiveManager );
  m_dResPhaseVisc_dPres =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString,
                                                                                      constitutiveManager );
  m_dResPhaseVisc_dComp =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString,
                                                                                     constitutiveManager );
  m_resPhaseCompFrac =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString,
                                                                                      constitutiveManager );
  m_dResPhaseCompFrac_dPres =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString,
                                                                                      constitutiveManager );
  m_dResPhaseCompFrac_dComp =
    elemManager->ConstructMaterialViewAccessor<array5d<real64>, arrayView5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                      constitutiveManager );
  m_resTotalDens =
    elemManager->ConstructMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString,
                                                                                      constitutiveManager );
  m_resPhaseRelPerm =
    elemManager->ConstructMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString,
                                                                                      constitutiveManager );
  m_dResPhaseRelPerm_dPhaseVolFrac =
    elemManager->ConstructMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                      constitutiveManager );
}


void CompositionalMultiphaseWell::FormMomentumEquations( DomainPartition * const domain,
							 Epetra_FECrsMatrix * const jacobian,
							 Epetra_FEVector * const residual )
{

}
  
void CompositionalMultiphaseWell::FormControlEquation( DomainPartition * const domain,
                                                       Epetra_FECrsMatrix * const jacobian,
                                                       Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData * const connectionData = well->getConnections();
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom 
    array1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // TODO: check that the first connection is on my rank
    
    // the well control equation is defined on the first connection
    localIndex const iconnControl = 0;
    // again we assume here that the first segment is on this MPI rank
    localIndex const iwelemControl = 0;

    // get well control
    Well::Control const control = well->getControl();

    // BHP control
    if (control == Well::Control::BHP)
    {
      // get pressure data
      array1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      array1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // form the control equation now
      real64 const currentBHP = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];
      real64 const targetBHP  = well->getTargetBHP();
      
      real64 const controlEqn = currentBHP - targetBHP;
      real64 const dControlEqn_dPres = 1;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      globalIndex const eqnRowIndex = elemOffset;
      globalIndex const dofColIndex = elemOffset;
      
      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex, &controlEqn );
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex, 1, &dofColIndex, &dControlEqn_dPres );
    }
    // rate control
    else
    {
    }
  });
}


void CompositionalMultiphaseWell::ImplicitStepComplete( real64 const & time,
                                                              real64 const & dt,
                                                              DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData const * const connectionData = well->getConnections();
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get a reference to the primary variables on segments
    array1d<real64> const & wellElemPressure  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array2d<real64> const & wellElemGlobalCompDensity  =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    array2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    // get a reference to the next well element index
    arrayView1d<localIndex const> const & nextWellElemIndex =
      connectionData->getReference<array1d<localIndex>>( ConnectionData::viewKeyStruct::nextWellElementIndexString );    
    
    // get a reference to the primary variables on connections
    array1d<real64> const & connRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::mixtureRateString );

    array1d<real64> const & dConnRate =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureRateString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
	wellElemGlobalCompDensity[iwelem][ic] += dWellElemGlobalCompDensity[iwelem][ic];
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // as a temporary solution, the connection variable is looked up with the next element's DOF number
      localIndex iwelemNext = nextWellElemIndex[iconn];

      // check if there is a variable defined here
      if (iwelemNext < 0)
        continue;

      connRate[iconn] += dConnRate[iconn];
    }    
  });  
}

void CompositionalMultiphaseWell::ComputeAllPerforationRates( Well * well )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS; 
  localIndex constexpr maxNumDof  = maxNumComp + 2; // dofs are 1 pressure, NC comp densities (reservoir and well), 1 connection rate (well only)

  localIndex const resNDOF  = numDofPerResElement();
  
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & resDofNumber = m_resDofNumber;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & resPressure  = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dResPressure = m_deltaResPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & resGravDepth = m_resGravDepth;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dResPhaseVolFrac_dPres = m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dResPhaseVolFrac_dComp = m_dResPhaseVolFrac_dCompDens;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dResCompFrac_dCompDens = m_dResCompFrac_dCompDens;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseDens        = m_resPhaseDens;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dResPhaseDens_dPres = m_dResPhaseDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dResPhaseDens_dComp = m_dResPhaseDens_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseVisc        = m_resPhaseVisc;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dResPhaseVisc_dPres = m_dResPhaseVisc_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dResPhaseVisc_dComp = m_dResPhaseVisc_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & resPhaseCompFrac    = m_resPhaseCompFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dResPhaseCompFrac_dPres = m_dResPhaseCompFrac_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView5d<real64>> const & dResPhaseCompFrac_dComp = m_dResPhaseCompFrac_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseRelPerm                = m_resPhaseRelPerm;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dResPhaseRelPerm_dPhaseVolFrac = m_dResPhaseRelPerm_dPhaseVolFrac;

  PerforationData const * const perforationData = well->getPerforations();
  WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

  // get the degrees of freedom
  array1d<globalIndex const> const & wellElemDofNumber =
    wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
  // get well primary variables on segments
  array1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  array1d<real64 const> const & dWellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  array2d<real64 const> const & dWellElemPhaseVolFrac_dPres =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  array3d<real64 const> const & dWellElemPhaseVolFrac_dComp =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  array3d<real64 const> const & dWellElemCompFrac_dCompDens =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  array1d<real64 const> const & wellElemGravDepth =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );
    
  // get well secondary variables on segments
  MultiFluidBase const * const fluid = GetFluidModel( wellElementSubRegion );
  RelativePermeabilityBase const * const relPerm = GetRelPermModel( wellElementSubRegion );
    
  array3d<real64 const> const & wellElemPhaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  array3d<real64 const> const & dWellElemPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

  array4d<real64 const> const & dWellElemPhaseDens_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

  array3d<real64 const> const & wellElemPhaseVisc =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString );

  array3d<real64 const> const & dWellElemPhaseVisc_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

  array4d<real64 const> const & dWellElemPhaseVisc_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );
    
  array4d<real64 const> const & wellElemPhaseCompFrac =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString );

  array4d<real64 const> const & dWellElemPhaseCompFrac_dPres =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString );

  array5d<real64 const> const & dWellElemPhaseCompFrac_dComp =
    fluid->getReference<array5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString );

  array3d<real64 const> const & wellElemPhaseRelPerm =
    relPerm->getReference<array3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString );

  array4d<real64> const & dWellElemPhaseRelPerm_dPhaseVolFrac =
    relPerm->getReference<array4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString );

  // get well variables on perforations
  array1d<real64 const> const & perfGravDepth =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::gravityDepthString );

  arrayView1d<localIndex const> const & perfWellElemIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d<real64 const> const & perfTransmissibility =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::transmissibilityString );

  // get the element region, subregion, index
  arrayView1d<localIndex const> const & resElementRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

  arrayView1d<localIndex const> const & resElementSubRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

  arrayView1d<localIndex const> const & resElementIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

  // TODO: rewrite this loop to match the SinglePhaseFluid implementations
    
  for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
  {
      
    // local working variables and arrays
    stackArray1d<long long, 2 * maxNumComp> eqnRowIndices( 2 * NC );
    stackArray1d<long long, 2 * maxNumDof>  dofColIndices( 2 * resNDOF );

    stackArray1d<double, 2 * maxNumComp>                 localFlux( 2 * NC );
    stackArray2d<double, 2 * maxNumComp * 2 * maxNumDof> localFluxJacobian( 2 * NC, 2 * resNDOF );
      
    stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

    stackArray1d<real64, 2> pressure( 2 );
      
    stackArray1d<real64, 2>              dPhaseFlux_dP( 2 );
    stackArray2d<real64, 2 * maxNumComp> dPhaseFlux_dC( 2, NC );

    stackArray1d<real64, maxNumComp>                  compFlux( NC );
    stackArray2d<real64, 2 * maxNumComp>              dCompFlux_dP( 2, NC );
    stackArray3d<real64, 2 * maxNumComp * maxNumComp> dCompFlux_dC( 2, NC, NC );

    stackArray1d<real64, maxNumComp> dRelPerm_dC( NC );
    stackArray1d<real64, maxNumComp> dDens_dC( NC );
    stackArray1d<real64, maxNumComp> dVisc_dC( NC );

    stackArray1d<real64, 2>              mobility( 2 );
    stackArray1d<real64, 2>              dMobility_dP( 2 );
    stackArray2d<real64, 2 * maxNumComp> dMobility_dC( 2, NC );

    stackArray1d<real64, 2>              dDensMean_dP( 2 );
    stackArray2d<real64, 2 * maxNumComp> dDensMean_dC( 2, NC );

    stackArray1d<real64, 2>              dPresGrad_dP( 2 );

    stackArray1d<real64, 2> multiplier( 2 );
      
    // reset the local values
    compFlux = 0.0;
    dCompFlux_dP = 0.0;
    dCompFlux_dC = 0.0;

    // 1) Reservoir side
      
    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];
      
    // loop over phases, compute and upwind phase flux
    // and sum contributions to each component's perforation rate
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // clear working arrays
      real64 presGrad = 0.0;
      dPresGrad_dP = 0.0;

      real64 phaseFlux = 0.0;
      dPhaseFlux_dP = 0.0;
      dPhaseFlux_dC = 0.0;

      // 1) get reservoir variables

      // pressure
      pressure[SubRegionTag::RES]  = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
	
      // TODO: check m_fluidIndex from reservoir and see if this is right
	
      // density
      real64 const resDensity  = resPhaseDens[er][esr][m_fluidIndex][ei][0][ip];
      real64 const dResDens_dP = dResPhaseDens_dPres[er][esr][m_fluidIndex][ei][0][ip];
      applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], dResPhaseDens_dComp[er][esr][m_fluidIndex][ei][0][ip], dDens_dC );

      // viscosity
      real64 const resViscosity = resPhaseVisc[er][esr][m_fluidIndex][ei][0][ip];
      real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][m_fluidIndex][ei][0][ip];
      applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], dResPhaseVisc_dComp[er][esr][m_fluidIndex][ei][0][ip], dVisc_dC );

      //relative permeability
      real64 const resRelPerm = resPhaseRelPerm[er][esr][m_relPermIndex][ei][0][ip];
      real64 dResRelPerm_dP = 0.0;
      dRelPerm_dC = 0.0;
      for (localIndex jp = 0; jp < NP; ++jp)
      {
        real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[er][esr][m_relPermIndex][ei][0][ip][jp];
        dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[er][esr][ei][jp];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[er][esr][ei][jp][jc];
        }
      }

      // mobility and pressure derivative
      mobility[SubRegionTag::RES] = resRelPerm * resDensity / resViscosity;
      dMobility_dP[SubRegionTag::RES] = dResRelPerm_dP * resDensity / resViscosity
        + mobility[SubRegionTag::RES] * (dResDens_dP / resDensity - dResVisc_dP / resViscosity);

      // compositional derivatives
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dDensMean_dC[SubRegionTag::RES][jc] = dDens_dC[jc];

        dMobility_dC[SubRegionTag::RES][jc] = dRelPerm_dC[jc] * resDensity / resViscosity
	  + mobility[SubRegionTag::RES] * (dDens_dC[jc] / resDensity - dVisc_dC[jc] / resViscosity);
      }

      multiplier[SubRegionTag::RES] = 1;

      // 2) get well variables

      localIndex const iwelem = perfWellElemIndex[iperf]; 
	
      // pressure
      pressure[SubRegionTag::WELL]  = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
	
      // density
      real64 const wellDensity  = wellElemPhaseDens[iwelem][0][ip];
      real64 const dWellDens_dP = dWellElemPhaseDens_dPres[iwelem][0][ip];
      applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseDens_dComp[iwelem][0][ip], dDens_dC );

      // viscosity
      real64 const wellViscosity = wellElemPhaseVisc[iwelem][0][ip];
      real64 const dWellVisc_dP  = dWellElemPhaseVisc_dPres[iwelem][0][ip];
      applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseVisc_dComp[iwelem][0][ip], dVisc_dC );

      //relative permeability
      real64 const wellRelPerm = wellElemPhaseRelPerm[iwelem][0][ip];
      real64 dWellRelPerm_dP = 0.0;
      dRelPerm_dC = 0.0;
      for (localIndex jp = 0; jp < NP; ++jp)
      {
        real64 const dWellRelPerm_dS = dWellElemPhaseRelPerm_dPhaseVolFrac[iwelem][0][ip][jp];
        dWellRelPerm_dP += dWellRelPerm_dS * dWellElemPhaseVolFrac_dPres[iwelem][jp];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dRelPerm_dC[jc] += dWellRelPerm_dS * dWellElemPhaseVolFrac_dComp[iwelem][jp][jc];
        }
      }

      // mobility and pressure derivative
      mobility[SubRegionTag::WELL] = wellRelPerm * wellDensity / wellViscosity;
      dMobility_dP[SubRegionTag::WELL] = dWellRelPerm_dP * wellDensity / wellViscosity
        + mobility[SubRegionTag::WELL] * (dWellDens_dP / wellDensity - dWellVisc_dP / wellViscosity);

      // compositional derivatives
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dMobility_dC[SubRegionTag::WELL][jc] = dRelPerm_dC[jc] * wellDensity / wellViscosity
	  + mobility[SubRegionTag::WELL] * (dDens_dC[jc] / wellDensity - dVisc_dC[jc] / wellViscosity);
      }

      multiplier[SubRegionTag::WELL] = -1;

      // @Sergey: below, this function is still in construction
      
      //***** calculation of flux *****
	
      // TODO: use distinct treatments for injector and producer
      // TODO: multiply mobility of injector by total reservoir mobility
      // TODO: account for depth difference between well element center and perforation
	
      // get transmissibility at the interface
      real64 const trans = perfTransmissibility[iperf]; 
	
      // compute potential difference
      for (localIndex i = 0; i < 2; ++i)
      {
        presGrad += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
        dPresGrad_dP[i] += multiplier[i] * trans;

        if (m_gravityFlag)
        {
	  // TODO: implement this the same way as in SinglePhaseWell
	}
      }

      // *** upwinding ***

      // compute phase potential gradient
      real64 potGrad = presGrad;

      // choose upstream cell
      localIndex const k_up = (potGrad >= 0) ? SubRegionTag::RES : SubRegionTag::WELL;

      // skip the phase flux if phase not present or immobile upstream
      if (std::fabs(mobility[k_up]) < 1e-20) // TODO better constant
      {
        break;
      }

      // pressure gradient depends on all points in the stencil
      for (localIndex ke = 0; ke < 2; ++ke)
      {
        dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
      }

      // compute the phase flux and derivatives using upstream cell mobility
      phaseFlux = mobility[k_up] * potGrad;
      for (localIndex ke = 0; ke < 2; ++ke)
      {
        dPhaseFlux_dP[ke] *= mobility[k_up];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] *= mobility[k_up];
        }
      }

      // add contribution from upstream cell mobility derivatives
      dPhaseFlux_dP[k_up] += dMobility_dP[k_up] * potGrad;
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseFlux_dC[k_up][jc] += dMobility_dC[k_up][jc] * potGrad;
      }

      // The code below currently assumes that the upstream cell is the reservoir cell, which will not always be true
      // TODO: change the syntax to support the case in which the well cell is upstream
	
      // get global identifiers of the upstream cell
      localIndex er_up  = er;
      localIndex esr_up = esr;
      localIndex ei_up  = ei;

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d<real64> phaseCompFracSub = resPhaseCompFrac[er_up][esr_up][m_fluidIndex][ei_up][0][ip];
      arraySlice1d<real64> dPhaseCompFrac_dPresSub = dResPhaseCompFrac_dPres[er_up][esr_up][m_fluidIndex][ei_up][0][ip];
      arraySlice2d<real64> dPhaseCompFrac_dCompSub = dResPhaseCompFrac_dComp[er_up][esr_up][m_fluidIndex][ei_up][0][ip];

      // compute component fluxes and derivatives using upstream cell composition
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const ycp = phaseCompFracSub[ic];
        compFlux[ic] += phaseFlux * ycp;

        // derivatives stemming from phase flux
        for (localIndex ke = 0; ke < 2; ++ke)
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
        applyChainRule( NC, dResCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dPhaseCompFrac_dCompDens );
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompFlux_dC[k_up][ic][jc] += phaseFlux * dPhaseCompFrac_dCompDens[jc];
        }
      }
    }
    
    //***** end upwinding *****

  }
}

void CompositionalMultiphaseWell::RecordWellData( Well * well )
{
  localIndex const NC = m_numComponents;
  
  ConnectionData const * const connectionData = well->getConnections();
  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
  PerforationData const * perforationData = well->getPerforations();
  
  array1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  array1d<real64 const> const & connRate =
    connectionData->getReference<array1d<real64>>( viewKeyStruct::mixtureRateString );

  array2d<real64 const> const & compPerfRate =
    perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString );

  // here, we will save the well info
  // for now, brute force: output to terminal

  std::cout << "Well : " << well->getName() << std::endl;
  if (well->getType() == Well::Type::PRODUCER)
    std::cout << "Type : PRODUCER" << std::endl;
  else 
    std::cout << "Type : INJECTOR" << std::endl;
  if (well->getControl() == Well::Control::BHP)
    std::cout << "Control : BHP" << std::endl;
  else
    std::cout << "Control : RATE" << std::endl;

  std::cout << "Below, positive perforation rate means flow from reservoir to well" << std::endl;
  std::cout << "Negative perforation rate means flow from well to reservoir" << std::endl;
  
  // output perforation rates
  for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
  {
    for (localIndex ic = 0; ic < NC; ++ic)
      std::cout << "Rate at perforation #" << iperf << " for component #" << ic << ": " << compPerfRate[iperf][ic] << std::endl;
  }

  // output the reference pressure
  localIndex const iwelemRef = 0; // be careful here for the parallel case
  real64 const pressure = wellElemPressure[iwelemRef];
  real64 const targetPressure = well->getTargetBHP();

  if (well->getControl() == Well::Control::BHP)
  {
    std::cout << "Current reference pressure = " << pressure
  	      << ", targetPressure = "           << targetPressure
	      << std::endl;
  }
  else
  {
    if (well->getType() == Well::Type::PRODUCER)
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

  std::cout << "Below, negative connection rate means production" << std::endl;
  std::cout << "Positive connection rate means injection" << std::endl;
  
  for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
  {
    if (well->getControl() == Well::Control::BHP || iconn > 0)
      std::cout << "Rate at connection #" << iconn << ": " << connRate[iconn] << std::endl;
    else
      std::cout << "Rate at connection #" << iconn << ": " << connRate[iconn]
		<< ", target rate : " << well->getTargetRate() << std::endl;
  }
}

  
void CompositionalMultiphaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {    
    // if isViable is true at the end of the following checks, no need to switch
    bool controlIsViable = false;

    // get well control and type
    Well::Control const currentControl = well->getControl();
    Well::Type const type = well->getType();
    
    // BHP control
    if (currentControl == Well::Control::BHP)
    {
      PerforationData * const perforationData = well->getPerforations();

      arrayView2d<real64 const> const & compPerfRate =
        perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString );

      // compute the local rates for this well
      // TODO: this is a waste of computations, we only need to know the sign of the potential difference
      ComputeAllPerforationRates( well );

      // the control is viable if at least one of the  perforations can inject / produce without cross flow
      // TODO: this loop will require communication
      for (localIndex iperf = 0; !controlIsViable && iperf < perforationData->numPerforationsLocal(); ++iperf)
      {	
	// producer should have positive difference; injector should have negative difference
	for (localIndex ic = 0; !controlIsViable && ic < NC; ++ic)
	{
   	  controlIsViable = ( ( type == Well::Type::PRODUCER && compPerfRate[iperf][ic] > 0 ) ||
                              ( type == Well::Type::INJECTOR && compPerfRate[iperf][ic] < 0 ) );
	}
      }
    }
    else // rate control
    {
      WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
      
      // get pressure data
      array1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      array1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // again we assume here that the first segment is on this MPI rank
      localIndex const iwelemRef = 0;
      
      // the control is viable if the reference pressure is below/above the max/min pressure
      real64 const refPressure = wellElemPressure[iwelemRef] + dWellElemPressure[iwelemRef];

      if ( type == Well::Type::PRODUCER )
      {
        real64 const minPressure = well->getTargetBHP(); // targetBHP specifies a min pressure here
	controlIsViable = ( refPressure >= minPressure );
      }
      else
      {
        real64 const maxPressure = well->getTargetBHP(); // targetBHP specifies a max pressure here
	controlIsViable = ( refPressure <= maxPressure );
      }
    }

    if (!controlIsViable)
    {
      if ( currentControl == Well::Control::BHP )
      {
	well->setControl( Well::Control::OILRATE );
        if ( m_verboseLevel >= 1 )
        {
          GEOS_LOG_RANK_0( "Control switch for well " << well->getName()
			   << " from BHP constraint to rate constraint" );
        }
      }
      else // rate control
      {
	well->setControl( Well::Control::BHP );
	if ( m_verboseLevel >= 1 )
        {
          GEOS_LOG_RANK_0( "Control switch for well " << well->getName()
			   << " from rate constraint to BHP constraint" );
        }
      }
    }
  });    
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseWell, string const &, ManagedGroup * const)
}// namespace geosx
