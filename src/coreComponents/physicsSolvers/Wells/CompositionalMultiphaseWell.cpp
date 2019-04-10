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
  m_numComponents( 0 ),
  m_resRelPermIndex()
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

  this->RegisterViewWrapper( viewKeyStruct::resRelPermIndexString,  &m_resRelPermIndex,  false );
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
    wellElementSubRegion->RegisterViewWrapper<array1d<globalIndex>>( viewKeyStruct::dofNumberString ); 
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalCompDensityString );
    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::mixtureConnRateString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );

    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalCompFractionString );
    wellElementSubRegion->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );
    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );
    wellElementSubRegion->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::sumCompPerforationRateString );
    
    PerforationData * const perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::compPerforationRateString );
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

    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );
    // TODO: put this out of the loop
    m_numPhases     = fluid->numFluidPhases();
    m_numComponents = fluid->numFluidComponents();

    localIndex const NC = m_numComponents;
    localIndex const NP = m_numPhases;

    m_numDofPerElement = NC + 2; // 1 pressure + NC compositions + 1 connectionRate
 
    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompDensityString).resizeDimension<1>(NC);
    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::deltaGlobalCompDensityString).resizeDimension<1>(NC);

    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompFractionString).resizeDimension<1>(NC);
    wellElementSubRegion->getReference<array3d<real64>>(viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString).resizeDimension<1,2>(NC, NC);

    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::phaseVolumeFractionString).resizeDimension<1>(NP);
    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dPressureString).resizeDimension<1>(NP);
    wellElementSubRegion->getReference<array3d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

    wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::sumCompPerforationRateString).resizeDimension<1>(NC);
    
    PerforationData * const perforationData = well->getPerforations();
    perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString ).resizeDimension<1>( NC );
    perforationData->getReference<array3d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString ).resizeDimension<1,2>( 2, NC );
    perforationData->getReference<array4d<real64>>( viewKeyStruct::dCompPerforationRate_dCompString ).resizeDimension<1,2,3>( 2, NC, NC );
    
  });

  ConstitutiveManager const * const cm = domain->getConstitutiveManager();
  
  RelativePermeabilityBase const * relPerm = cm->GetConstitituveRelation<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Relative permeability model " + m_relPermName + " not found" );
 
  m_resRelPermIndex = relPerm->getIndexInParent(); 

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
  MultiFluidBase const * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );

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
  MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );

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
  RelativePermeabilityBase * const relPerm = GetConstitutiveModel<RelativePermeabilityBase>( wellElementSubRegion, m_relPermName );

  arrayView2d<real64 const> const & wellElemPhaseVolFrac =
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
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseRelPerm  = m_resPhaseRelPerm;

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
    PerforationData const * const perforationData = well->getPerforations();

    // get well primary variables on well elements
    arrayView1d<real64> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView2d<real64> const & wellElemCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    arrayView1d<real64> const & connRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );
    
    // get the info stored on well elements
    arrayView2d<real64> const & wellElemCompFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

    arrayView2d<real64> const & wellElemSumCompPerfRates =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::sumCompPerforationRateString );

    arrayView1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );

    arrayView2d<real64> const & wellElemPhaseVolFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );
    
    // get the well element indices corresponding to each connect
    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    

    // get well data on perforations
    arrayView2d<real64> const & compPerfRate =
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

    // get well secondary variables on well elements
    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );

    arrayView2d<real64> const & totalDens =
      fluid->getReference<array2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString );
    
    arrayView3d<real64 const> const & wellElemPhaseDens =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

    arrayView4d<real64 const> const & wellElemPhaseCompFrac =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString );

    
    // TODO: fix this in the case of multiple MPI ranks working on the same well
    // In particular, steps 1), 3) and 7) will need to be rewritten
    
   /*
    * The methodology used here is based on Yifan Zhou's PhD dissertation, Stanford University
    */
    
    // 1) Loop over all perforations to compute an average mixture density and component fraction
    real64 avgTotalDensity   = 0.0;
    real64 avgMixtureDensity = 0.0;
    stackArray1d<real64, maxNumComp> avgCompFrac( NC );
    
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
        real64 const resViscosity = resPhaseVisc[er][esr][m_resFluidIndex][ei][0][ip];
        real64 const resRelPerm   = resPhaseRelPerm[er][esr][m_resRelPermIndex][ei][0][ip];
        
        resTotalMobility += resRelPerm / resViscosity;
      }

      // increment the average mixture density
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        real64 const resDensity   = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip];
        real64 const resViscosity = resPhaseVisc[er][esr][m_resFluidIndex][ei][0][ip];
        real64 const resRelPerm   = resPhaseRelPerm[er][esr][m_resRelPermIndex][ei][0][ip];
        
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
      }
    }

    // get the reference data for this well
    localIndex const iwelemRef = 0; // here, we assume that the first well element is on this MPI rank, which will NOT be the case in general
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
                                  ? 0.9 * resPres // hard-coded value comes from AD-GPRS
                                  : 1.1 * resPres; 
    }

    // TODO: communicate ref pressure
    // TODO: communicate avgDens
    
    // 3) Estimate the pressures in the well elements using this avgDensity
    // TODO: implement this in parallel with the communication of the pressures
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (iwelem != iwelemRef)
        wellElemPressure[iwelem] = wellElemPressure[iwelemRef]
          + ( m_gravityFlag ? avgMixtureDensity * ( wellElemGravDepth[iwelem] - refGravDepth ) : 0 );
    }
    
    // 4) Back calculate component densities
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        fluid->PointUpdate( wellElemPressure[iwelem],
                            m_temperature,
                            wellElemCompFrac[iwelem],
                            iwelem,
                            0 );
        wellElemCompDens[iwelem][ic] = avgCompFrac[ic] * totalDens[iwelem][0];
      }
    }

    // 5) Recompute all the pressure-dependent properties
    UpdateState( well );
    
    // 6) Compute the perforation rates
    ComputeAllPerforationRates( well );

    // 7) Collect all the perforation rates
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      for (localIndex ic = 0; ic < NC; ++ic)
        wellElemSumCompPerfRates[iwelem][ic] = 0.0; // TODO: ask for a better way to do that!
    }
    
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {  
      localIndex const iwelem = perfWellElemIndex[iperf];
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        wellElemSumCompPerfRates[iwelem][ic] += compPerfRate[iperf][ic];
      }
    }

    // 8) Estimate the connection rates
    // TODO: implement this in parallel with the communication of wellElemSumRates
    // WARNING FRANCOIS: this assumes that connections are ordered!!!
    real64 prevConnRate = 0;
    for (localIndex iwelem = wellElementSubRegion->numWellElementsLocal()-1; 
         iwelem >= 0; --iwelem)
    {
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      connRate[iwelem] = prevConnRate;
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        if ( wellElemPhaseVolFrac[iwelem][ip] > 0.1 &&
             wellElemPhaseDens[iwelem][0][ip] > std::numeric_limits<real64>::epsilon() )
        {
          real64 multiplier = wellElemPhaseDens[iwelem][0][ip]
                            * wellElemPhaseVolFrac[iwelem][ip];
          for (localIndex ic = 0; ic < NC; ++ic)
          {
            if ( wellElemPhaseCompFrac[iwelem][0][ip][ic] > 0.1 )
            {
              connRate[iwelem] -= wellElemSumCompPerfRates[iwelem][ic]
                / ( wellElemPhaseCompFrac[iwelem][0][ip][ic] * multiplier );
            }
          }
        }
      }
      std::cout << "Predicted connection rate at #" << iwelem
                << " = " << connRate[iwelem] << std::endl;
      prevConnRate = connRate[iwelem];
    }
  });
}

void CompositionalMultiphaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );
}

void CompositionalMultiphaseWell::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                     DomainPartition * const domain,
                                                     EpetraBlockSystem * const blockSystem )
{
  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // Initialize the primary and secondary variables
  InitializeWells( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );
  
  // assumes that the setup of dof numbers and linear system
  // is done in ReservoirWellSolver

}

void CompositionalMultiphaseWell::SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
                                                                localIndex & numLocalRows,
                                                                globalIndex & numGlobalRows,
                                                                localIndex offset )
{
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
  localIndex constexpr maxNumDof  = maxNumComp + 1; // dofs are 1 pressure, NC comp densities (reservoir and well), no connection rate needed

  // reservoir dofs are 1 pressure and NC comp densities
  localIndex const resNDOF  = numDofPerResElement;
  // well dofs are 1 pressure, NC comp densities, 1 connection rate
  localIndex const wellNDOF = numDofPerElement(); 

  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    PerforationData const * const perforationData = well->getPerforations();
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the well degrees of freedom 
    arrayView1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get the well element indices corresponding to each connect
    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    

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

    // 2) Insert the entries corresponding to well-well connection between well elements
    //    This will fill J_WW only
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      // get next well element index
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      // check if this is not an entry or exit
      if (iwelemNext < 0)
      {
        continue;
      }
      
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( 2 * wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( 2 * wellNDOF );

      // get the offset of the well element equations
      globalIndex const currentElemOffset = getElementOffset( wellElemDofNumber[iwelem] );
      globalIndex const nextElemOffset    = getElementOffset( wellElemDofNumber[iwelemNext] );

      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
        elementLocalDofIndexRow[idof]            = currentElemOffset + idof;
        elementLocalDofIndexRow[wellNDOF + idof] = nextElemOffset    + idof;
        elementLocalDofIndexCol[idof]            = currentElemOffset + idof;
        elementLocalDofIndexCol[wellNDOF + idof] = nextElemOffset    + idof;
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
  FormPressureRelations( domain, jacobian, residual );

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
  WellManager * const wellManager = domain->getWellManager();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;
  
  localIndex const NC      = m_numComponents;
  localIndex const NP      = m_numPhases;
  localIndex const resNDOF = numDofPerResElement();
  
  // loop over the wells
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get a reference to the degree-of-freedom numbers and to ghosting info
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get a reference to the primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64 const> const & wellElemCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    arrayView2d<real64> const & dWellElemCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64 const> const & connRate  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

    arrayView1d<real64 const> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );
    
    // get the info stored on well elements
    arrayView2d<real64> const & wellElemCompFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

    arrayView3d<real64 const> const & dWellElemCompFrac_dCompDens =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    arrayView2d<real64> const & wellElemPhaseVolFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

    arrayView2d<real64> const & dWellElemPhaseVolFrac_dPres =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    arrayView3d<real64> const & dWellElemPhaseVolFrac_dComp =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    
    // get well secondary variables on well elements
    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );
   
    arrayView2d<real64> const & wellElemTotalDens =
      fluid->getReference<array2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString );

    arrayView3d<real64 const> const & wellElemPhaseDens =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

    arrayView3d<real64 const> const & dWellElemPhaseDens_dPres =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

    arrayView4d<real64 const> const & dWellElemPhaseDens_dComp =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );
    
    arrayView4d<real64 const> const & wellElemPhaseCompFrac =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString );

    arrayView4d<real64 const> const & dWellElemPhaseCompFrac_dPres =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString );

    arrayView5d<real64 const> const & dWellElemPhaseCompFrac_dComp =
      fluid->getReference<array5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString );

    // create local work arrays
    stackArray1d<real64, maxNumComp> dDens_dC( NC );
    stackArray1d<real64, maxNumComp> dVolFrac_dC( NC );
    stackArray1d<real64, maxNumComp> dPhaseFlux_dCompDensUp( NC );

    stackArray1d<real64, maxNumComp>              phaseCompFrac( NC );
    stackArray1d<real64, maxNumComp>              dPhaseCompFrac_dP( NC );
    stackArray2d<real64, maxNumComp * maxNumComp> dPhaseCompFrac_dC( NC, NC );

    stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

    stackArray1d<real64, maxNumComp>              compFlux( NC );
    stackArray1d<real64, maxNumComp>              dCompFlux_dRate( NC );
    stackArray1d<real64, maxNumComp>              dCompFlux_dPresUp( NC );
    stackArray2d<real64, maxNumComp * maxNumComp> dCompFlux_dCompDensUp( NC, NC );
    
    stackArray1d<real64, 2 * maxNumComp>             localFlux( 2 * NC );
    stackArray1d<real64, 2 * maxNumComp>             localFluxJacobian_dRate( 2 * NC );
    stackArray2d<real64, 2 * maxNumComp * maxNumDof> localFluxJacobian_dPresCompUp( 2 * NC, resNDOF );

    stackArray1d<long long, 2 * maxNumComp> eqnRowIndices( 2 * NC );
    stackArray1d<long long, maxNumDof>      dofColIndices_dPresCompUp( resNDOF );
    
    // loop over the well elements to compute the fluxes between elements    
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      localIndex const iwelemNext = nextWellElemIndex[iwelem];
      
      // Step 1) prepare variables
        
      /*
      if (iwelemPrev < 0 && well->getType() == Well::Type::INJECTOR)
      {
          
        // get the global component fraction from XML file input
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          wellElemCompFrac[iwelemNext][ic] = well->getInjectionStream( ic );
        }

        // TODO: improve that later
          
        real64 const currentDeltaPressure = dWellElemPressure[iwelemNext];
        dWellElemPressure[iwelemNext] = well->getTargetBHP()
          - wellElemPressure[iwelemNext];
         
        stackArray1d<real64, maxNumComp> currentDeltaCompDens( NC );
        for (localIndex ic = 0; ic < NC; ++ic)
          currentDeltaCompDens[ic] = dWellElemCompDens[iwelemNext][ic];
          
        // update fluid model
          
        fluid->PointUpdate( wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext],
                            m_temperature,
                            wellElemCompFrac[iwelemNext],
                            iwelemNext,
                            0 );
          
          
        // back calculate the component densities
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dWellElemCompDens[iwelemNext][ic] = wellElemCompFrac[iwelemNext][ic] * wellElemTotalDens[iwelemNext][0]
                                            - wellElemCompDens[iwelemNext][ic];
        }

          
        // reupdate the constitutive variables to get the correct derivative
        // TODO: allow pointwise update in the well elements to reduce cost
        UpdateState( well );

        // restore to current delta values
        dWellElemPressure[iwelemNext] = currentDeltaPressure;
        for (localIndex ic = 0; ic < NC; ++ic)
          dWellElemCompDens[iwelemNext][ic] = currentDeltaCompDens[ic];
          
      }
*/
          
      // Step 2) decide the upwind well element

      real64 const currentConnRate = connRate[iwelem] + dConnRate[iwelem];

      localIndex iwelemUp = -1;

      /*  currentConnRate < 0 flow from iwelem to iwelemNext
       *  currentConnRate > 0 flow from iwelemNext to iwelem
       *  With this convention, currentConnRate < 0 at the last connection for a producer
       *                        currentConnRate > 0 at the last connection for a injector
       */

          
      if (iwelemNext < 0 || // exit connection
          currentConnRate < 0) // iwelemNext is upstream
      { 
        iwelemUp = iwelemNext;
      }
      else // iwelemNext is downstream
      {
        iwelemUp = iwelem;
      }

      // Step 3) compute upstream transport coefficient

      compFlux = 0.0;
      dCompFlux_dRate = 0.0;
      dCompFlux_dPresUp = 0.0;
      dCompFlux_dCompDensUp = 0.0;

      for (localIndex ip = 0; ip < NP; ++ip)
      { 
        dDens_dC = 0.0;
        dVolFrac_dC = 0.0;
        dPhaseCompFrac_dC = 0.0;
        dPhaseFlux_dCompDensUp = 0.0;
          
        // density
        real64 const density  = wellElemPhaseDens[iwelemUp][0][ip];
        real64 const dDens_dP = dWellElemPhaseDens_dPres[iwelemUp][0][ip];
        applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelemUp], dWellElemPhaseDens_dComp[iwelemUp][0][ip], dDens_dC );

        // volume fraction
        real64 const volFrac = wellElemPhaseVolFrac[iwelemUp][ip];
        real64 const dVolFrac_dP = dWellElemPhaseVolFrac_dPres[iwelemUp][ip];
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dVolFrac_dC[ic] = dWellElemPhaseVolFrac_dComp[iwelemUp][ip][ic];
        }
          
        // component phase fractions 
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          phaseCompFrac[ic] = wellElemPhaseCompFrac[iwelemUp][0][ip][ic];
          dPhaseCompFrac_dP[ic] = dWellElemPhaseCompFrac_dPres[iwelemUp][0][ip][ic];
            
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dPhaseCompFrac_dC[ic][jc] = dWellElemPhaseCompFrac_dComp[iwelemUp][0][ip][ic][jc];
          }
        }

        // compute phase flux
        real64 const phaseFlux          = volFrac * density * currentConnRate;
        real64 const dPhaseFlux_dRate   = volFrac * density;
        real64 const dPhaseFlux_dPresUp = (dVolFrac_dP * density + volFrac * dDens_dP ) * currentConnRate;
          
        for ( localIndex ic = 0; ic < NC; ++ic)
        {
          dPhaseFlux_dCompDensUp[ic] = (dVolFrac_dC[ic] * density + volFrac * dDens_dC[ic] )
            * currentConnRate;
        }
          
        // compute component flux
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          compFlux[ic]          += phaseCompFrac[ic] * phaseFlux;
          dCompFlux_dRate[ic]   += phaseCompFrac[ic] * dPhaseFlux_dRate;
          dCompFlux_dPresUp[ic] += dPhaseCompFrac_dP[ic] * phaseFlux + phaseCompFrac[ic] * dPhaseFlux_dPresUp;

          applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelemUp], dPhaseCompFrac_dC[ic], dPhaseCompFrac_dCompDens );
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dCompFlux_dCompDensUp[ic][jc] += dPhaseCompFrac_dCompDens[jc] * phaseFlux
              + phaseCompFrac[ic] * dPhaseFlux_dCompDensUp[jc];
          }
        }
      }

      // TODO: check for ghost well elements
        
      globalIndex const offsetUp      = getElementOffset( wellElemDofNumber[iwelemUp] );
      globalIndex const offsetCurrent = getElementOffset( wellElemDofNumber[iwelem] );
      globalIndex const offsetNext    = getElementOffset( wellElemDofNumber[iwelemNext] );

      if ( iwelemNext < 0 ) // exit connection
      {

        // for this case, we only need NC mass conservation equations
        // so we do not use the arrays initialized before the loop
        stackArray1d<real64, maxNumComp>             oneSidedFlux( NC );
        stackArray1d<real64, maxNumComp>             oneSidedFluxJacobian_dRate( NC );
        stackArray2d<real64, maxNumComp * maxNumDof> oneSidedFluxJacobian_dPresCompUp( NC, resNDOF );

        stackArray1d<long long, maxNumComp> oneSidedEqnRowIndices( NC );
        stackArray1d<long long, maxNumDof>  oneSidedDofColIndices_dPresCompUp( resNDOF );

        // flux terms
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          oneSidedFlux[ic] = - dt * compFlux[ic];

          // derivative with respect to rate
          oneSidedFluxJacobian_dRate[ic] = - dt * dCompFlux_dRate[ic];
          
          // derivative with respect to upstream pressure
          oneSidedFluxJacobian_dPresCompUp[ic][ColOffset::DPRES] = - dt * dCompFlux_dPresUp[ic];
          
          // derivatives with respect to upstream component densities
          for (localIndex jdof = 0; jdof < NC; ++jdof)
          {
            oneSidedFluxJacobian_dPresCompUp[ic][ColOffset::DPRES + jdof + 1] = - dt * dCompFlux_dCompDensUp[ic][jdof];
          }
            
        }

        // jacobian indices
        for (localIndex ic = 0; ic < NC; ++ic )
        {
          // mass balance equations for all components
          oneSidedEqnRowIndices[ic] = offsetUp + RowOffset::MASSBAL + ic;
        }
          
        // in the dof ordering used in this class, there are 1 pressure dofs 
        // and NC compDens dofs before the rate dof in this block
        globalIndex const oneSidedDofColIndices_dRate = offsetCurrent + NC + 1; 
          
        for (localIndex jdof = 0; jdof < resNDOF; ++jdof)
        {
          // compositions are given by injection stream,  hence no derivatives
          if ( well->getType() == Well::Type::INJECTOR &&
               jdof >= 1 )
          {
            continue; 
          }
            
          // dofs are the **upstream** pressure and component densities
          oneSidedDofColIndices_dPresCompUp[jdof] = offsetUp + ColOffset::DPRES + jdof;
        }

        // Add to global residual/jacobian
        residual->SumIntoGlobalValues( integer_conversion<int>( NC ),
                                       oneSidedEqnRowIndices.data(),
                                       oneSidedFlux.data() );

        // derivatives with respect to rate
        jacobian->SumIntoGlobalValues( integer_conversion<int>( NC ),
                                       oneSidedEqnRowIndices.data(),
                                       1,
                                       &oneSidedDofColIndices_dRate,
                                       oneSidedFluxJacobian_dRate.data(),
                                       Epetra_FECrsMatrix::ROW_MAJOR );

        // derivatives with respect to upstream pressure and component densities
        jacobian->SumIntoGlobalValues( integer_conversion<int>( NC ),
                                       oneSidedEqnRowIndices.data(),
                                       integer_conversion<int>( resNDOF ),
                                       oneSidedDofColIndices_dPresCompUp.data(),
                                       oneSidedFluxJacobian_dPresCompUp.data(),
                                       Epetra_FECrsMatrix::ROW_MAJOR );

        /*
        if ( well->getType() == Well::Type::INJECTOR )
        {
          // reupdate the constitutive variables to get the correct derivative
          // TODO: allow pointwise update in the well elements to reduce cost
          UpdateState( well );
        }*/
          
      }
      else // not an exit connection
      {

        // flux terms
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          localFlux[ElemTag::NEXT * NC + ic]    =   dt * compFlux[ic];
          localFlux[ElemTag::CURRENT * NC + ic] = - dt * compFlux[ic];

          // derivative with respect to rate
          localFluxJacobian_dRate[ElemTag::NEXT * NC + ic]    =   dt * dCompFlux_dRate[ic];
          localFluxJacobian_dRate[ElemTag::CURRENT * NC + ic] = - dt * dCompFlux_dRate[ic];
          
          // derivative with respect to upstream pressure
          localFluxJacobian_dPresCompUp[ElemTag::NEXT * NC + ic][ColOffset::DPRES]    =    dt * dCompFlux_dPresUp[ic];
          localFluxJacobian_dPresCompUp[ElemTag::CURRENT * NC + ic][ColOffset::DPRES] =  - dt * dCompFlux_dPresUp[ic];
          
          // derivatives with respect to upstream component densities
          for (localIndex jdof = 0; jdof < NC; ++jdof)
          {
            localFluxJacobian_dPresCompUp[ElemTag::NEXT * NC + ic][ColOffset::DPRES + jdof + 1] =   
                dt * dCompFlux_dCompDensUp[ic][jdof];
            localFluxJacobian_dPresCompUp[ElemTag::CURRENT * NC + ic][ColOffset::DPRES + jdof + 1] = 
              - dt * dCompFlux_dCompDensUp[ic][jdof];
          }
            
        }

        // jacobian indices
        for (localIndex ic = 0; ic < NC; ++ic )
        {
          // mass balance equations for all components
          eqnRowIndices[ElemTag::NEXT * NC + ic]    = offsetNext    + RowOffset::MASSBAL + ic;
          eqnRowIndices[ElemTag::CURRENT * NC + ic] = offsetCurrent + RowOffset::MASSBAL + ic;
        }
          
        // in the dof ordering used in this class, there are 1 pressure dofs 
        // and NC compDens dofs before the rate dof in this block
        globalIndex const dofColIndices_dRate = offsetCurrent + NC + 1;
          
        for (localIndex jdof = 0; jdof < resNDOF; ++jdof)
        {
          // dofs are the **upstream** pressure and component densities
          dofColIndices_dPresCompUp[jdof] = offsetUp + ColOffset::DPRES + jdof;
        }

        // Add to global residual/jacobian
        residual->SumIntoGlobalValues( integer_conversion<int>( 2 * NC ),
                                       eqnRowIndices.data(),
                                       localFlux.data() );

        // derivatives with respect to rate
        jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 * NC ),
                                       eqnRowIndices.data(),
                                       1,
                                       &dofColIndices_dRate,
                                       localFluxJacobian_dRate.data(),
                                       Epetra_FECrsMatrix::ROW_MAJOR );

        // derivatives with respect to upstream pressure and component densities
        jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 * NC ),
                                       eqnRowIndices.data(),
                                       integer_conversion<int>( resNDOF ),
                                       dofColIndices_dPresCompUp.data(),
                                       localFluxJacobian_dPresCompUp.data(),
                                       Epetra_FECrsMatrix::ROW_MAJOR );
      }
    }
  });
}

void CompositionalMultiphaseWell::AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                                              Epetra_FECrsMatrix * const jacobian,
                                                              Epetra_FEVector * const residual,
                                                              real64 const time_n,
                                                              real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();
  
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC        = m_numComponents;
  localIndex const NP        = m_numPhases;
  localIndex const welemNDOF = m_numDofPerElement;

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData const * const perforationData = well->getPerforations();
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
  
    // get the degrees of freedom and ghosting info
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    // get the properties on the well element
    arrayView2d<real64 const> const & wellElemPhaseVolFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

    arrayView2d<real64 const> const & dWellElemPhaseVolFrac_dPres =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    arrayView3d<real64 const> const & dWellElemPhaseVolFrac_dComp =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    arrayView1d<real64 const> const & wellElemVolume =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::wellElementVolumeString );

    stackArray1d<long long, maxNumDof> localVolBalanceDOF( welemNDOF );
    stackArray1d<double, maxNumDof>    localVolBalanceJacobian( welemNDOF );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->size(); ++iwelem)
    {

      if (wellElemGhostRank[iwelem] >= 0)
      {
        continue;
      }
      
      localVolBalanceDOF = 0;
      localVolBalanceJacobian = 0;

      // get volume (assume porosity is 1)
      real64 const volume = wellElemVolume[iwelem];

      // get equation/dof indices
      globalIndex const offset = getElementOffset( wellElemDofNumber[iwelem] );
      globalIndex const localVolBalanceEqnIndex = offset + RowOffset::MASSBAL + NC;
      for (localIndex jdof = 0; jdof < welemNDOF; ++jdof)
      {
        localVolBalanceDOF[jdof] = offset + ColOffset::DPRES + jdof;
      }

      real64 localVolBalance = 1.0;
      localVolBalanceJacobian = 0.0;

      // sum contributions to component accumulation from each phase
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        localVolBalance -= wellElemPhaseVolFrac[iwelem][ip];
        localVolBalanceJacobian[0] -= dWellElemPhaseVolFrac_dPres[iwelem][ip];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          localVolBalanceJacobian[jc+1] -= dWellElemPhaseVolFrac_dComp[iwelem][ip][jc];
        }
      }

      // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
      for (localIndex idof = 0; idof < welemNDOF; ++idof)
      {
        localVolBalanceJacobian[idof] *= volume;
      }
      localVolBalance *= volume;

      // add contribution to global residual and dRdP
      residual->SumIntoGlobalValues( 1,
                                     &localVolBalanceEqnIndex,
                                     &localVolBalance );

      jacobian->SumIntoGlobalValues( 1,
                                     &localVolBalanceEqnIndex,
                                     integer_conversion<int>( welemNDOF ),
                                     localVolBalanceDOF.data(),
                                     localVolBalanceJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR );
    }
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
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // get well variables on perforations
    arrayView2d<real64 const> const & compPerfRate =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString );

    arrayView3d<real64 const> const & dCompPerfRate_dPres =
      perforationData->getReference<array3d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString );

    arrayView4d<real64 const> const & dCompPerfRate_dComp =
      perforationData->getReference<array4d<real64>>( viewKeyStruct::dCompPerforationRate_dCompString );
    
    arrayView1d<localIndex const> const & perfWellElemIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // local working variables and arrays
    stackArray1d<long long, 2 * maxNumComp> eqnRowIndices( 2 * NC );
    stackArray1d<long long, 2 * maxNumDof>  dofColIndices( 2 * resNDOF );

    stackArray1d<double, 2 * maxNumComp>                 localFlux( 2 * NC );
    stackArray2d<double, 2 * maxNumComp * 2 * maxNumDof> localFluxJacobian( 2 * NC, 2 * resNDOF );

    // loop over the perforations and add the rates to the residual and jacobian
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {

      // local working variables and arrays
      eqnRowIndices = -1;
      dofColIndices = -1;

      localFlux = 0;
      localFluxJacobian = 0;

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const resOffset      = resNDOF * resDofNumber[er][esr][ei];
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
    arrayView1d<globalIndex const > const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer const> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

        localIndex const wellNDOF = numDofPerElement(); 
        for (localIndex idof = 0; idof < wellNDOF; ++idof)
        {
          int const lid = rowMap->LID(integer_conversion<int>( elemOffset + idof ));
          real64 const val = localResidual[lid];
          std::cout << "lid = " << lid << " localResidual = " << localResidual[lid] << std::endl;
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
  // not implemented yet

  return true;
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
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers on well elements, next elem index, ghosting info
    arrayView1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );
    
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract solution and apply to dP
        globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

        int lid = rowMap->LID( integer_conversion<int>( elemOffset ) );
        dWellElemPressure[iwelem] += scalingFactor * local_solution[lid];
        std::cout << "pressure: local_solution " << local_solution[lid] << std::endl;
        
        for (localIndex ic = 0; ic < m_numComponents; ++ic)
        {
          lid = rowMap->LID( integer_conversion<int>( elemOffset + ic + 1 ) );
          dWellElemGlobalCompDensity[iwelem][ic] += scalingFactor * local_solution[lid];
          std::cout << "compDens #" << ic << ": local_solution " << local_solution[lid] << std::endl;
        }

        lid = rowMap->LID( integer_conversion<int>( elemOffset + m_numComponents + 1 ) );
        dConnRate[iwelem] += scalingFactor * local_solution[lid];
        std::cout << "rate: local_solution " << local_solution[lid] << std::endl;

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
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the ghosting information and next well elem index
    arrayView1d<integer> const & wellElemGhostRank =
      wellElementSubRegion->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );
    
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      if (wellElemGhostRank[iwelem] < 0)
      {
        // extract solution and apply to dP
        dWellElemPressure[iwelem] = 0;
        dConnRate[iwelem] = 0;
        for (localIndex ic = 0; ic < m_numComponents; ++ic)
        {
          dWellElemGlobalCompDensity[iwelem][ic] = 0;
        }
      }
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
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString,
                                                                                          constitutiveManager );
  m_dResPhaseFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString,
                                                                                          constitutiveManager );
  m_dResPhaseFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString,
                                                                                          constitutiveManager );
  m_resPhaseDens =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                          constitutiveManager );
  m_dResPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                          constitutiveManager );
  m_dResPhaseDens_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString,
                                                                                          constitutiveManager );
  m_resPhaseVisc =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString,
                                                                                          constitutiveManager );
  m_dResPhaseVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString,
                                                                                          constitutiveManager );
  m_dResPhaseVisc_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString,
                                                                                          constitutiveManager );
  m_resPhaseCompFrac =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString,
                                                                                          constitutiveManager );
  m_dResPhaseCompFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString,
                                                                                          constitutiveManager );
  m_dResPhaseCompFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array5d<real64>, arrayView5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                          constitutiveManager );
  m_resTotalDens =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString,
                                                                                          constitutiveManager );
  m_resPhaseRelPerm =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString,
                                                                                          constitutiveManager );
  m_dResPhaseRelPerm_dPhaseVolFrac =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                          constitutiveManager );
}


void CompositionalMultiphaseWell::FormPressureRelations( DomainPartition * const domain,
                                                         Epetra_FECrsMatrix * const jacobian,
                                                         Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS; 
  localIndex constexpr maxNumDof  = maxNumComp + 1; // here, dofs are 1 pressure, NC comp densities (reservoir and well)
  
  localIndex const resNDOF  = numDofPerResElement();

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom, depth info, next welem index
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    // get secondary data on well elements    
    arrayView3d<real64 const> const & dWellElemCompFrac_dCompDens =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
   
    // get well constitutive data
    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );

    arrayView3d<real64 const> const & wellElemPhaseDens =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

    arrayView3d<real64 const> const & dWellElemPhaseDens_dPres =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

    arrayView4d<real64 const> const & dWellElemPhaseDens_dComp =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

    // local working variables and arrays
    stackArray1d<globalIndex, 2 * maxNumDof > dofColIndices( 2 * resNDOF );
    stackArray1d<real64, 2 * maxNumDof > localFluxJacobian( 2 * resNDOF );

    stackArray1d<real64, maxNumComp> dDensNext_dC( NC );
    stackArray1d<real64, maxNumComp> dDensCurrent_dC( NC );
    
    stackArray1d<real64, maxNumComp> dAvgDensity_dCompNext( NC );
    stackArray1d<real64, maxNumComp> dAvgDensity_dCompCurrent( NC );
    
    // loop over the well elements to compute the pressure relations between well elements
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      // TODO: check for ghost well elements
      
      if ( iwelemNext < 0 )  // if iwelemNext < 0, form control equation, not momentum
      {
        continue;
      }      

      // reset local working variables and arrays
      dofColIndices = -1;
      localFluxJacobian = 0.0;

      dDensNext_dC = 0.0;
      dDensCurrent_dC = 0.0;
    
      dAvgDensity_dCompNext    = 0.0;
      dAvgDensity_dCompCurrent = 0.0;
      real64 dAvgDensity_dPresNext    = 0.0;
      real64 dAvgDensity_dPresCurrent = 0.0;
        
      real64 avgDensity = 0;

      // TODO: precompte avg densities somewhere before this
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        avgDensity += 0.5 * ( wellElemPhaseDens[iwelemNext][0][ip]
                            + wellElemPhaseDens[iwelem][0][ip] ) / NP;
        dAvgDensity_dPresNext    += 0.5 * dWellElemPhaseDens_dPres[iwelemNext][0][ip] / NP;       
        dAvgDensity_dPresCurrent += 0.5 * dWellElemPhaseDens_dPres[iwelem][0][ip] / NP;

        applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelemNext], dWellElemPhaseDens_dComp[iwelemNext][0][ip], dDensNext_dC );
        applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseDens_dComp[iwelem][0][ip], dDensCurrent_dC );
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dAvgDensity_dCompNext[ic]    += 0.5 * dDensNext_dC[ic] / NP;
          dAvgDensity_dCompCurrent[ic] += 0.5 * dDensCurrent_dC[ic] / NP;
        }
      }
        
      // compute depth diff times acceleration
      real64 const gravD = ( wellElemGravDepth[iwelemNext] - wellElemGravDepth[iwelem] );

      // compute the current pressure in the two well elements
      real64 const pressureNext    = wellElemPressure[iwelemNext]    + dWellElemPressure[iwelemNext];
      real64 const pressureCurrent = wellElemPressure[iwelem] + dWellElemPressure[iwelem];

      // compute a coefficient to normalize the momentum equation
      real64 normalizer = (pressureNext > pressureCurrent)
                        ? pressureNext 
                        : pressureCurrent;
      normalizer = 1. / normalizer;

      // compute momentum flux and derivatives
      localIndex const localDofIndexPresNext    = ElemTag::NEXT * resNDOF;
      localIndex const localDofIndexPresCurrent = ElemTag::CURRENT * resNDOF;

      globalIndex const offsetNext    = getElementOffset( wellElemDofNumber[iwelemNext] ); 
      globalIndex const offsetCurrent = getElementOffset( wellElemDofNumber[iwelem] );

      globalIndex const eqnRowIndex = offsetCurrent + RowOffset::CONTROL; // according to convention used in this class
      dofColIndices[localDofIndexPresNext]    = offsetNext    + ColOffset::DPRES;
      dofColIndices[localDofIndexPresCurrent] = offsetCurrent + ColOffset::DPRES;

      real64 const localFlux = ( pressureNext - pressureCurrent - avgDensity * gravD ) * normalizer;
      localFluxJacobian[localDofIndexPresNext]    = ( 1 - dAvgDensity_dPresNext * gravD )    * normalizer;
      localFluxJacobian[localDofIndexPresCurrent] = (-1 - dAvgDensity_dPresCurrent * gravD ) * normalizer;

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        localIndex const localDofIndexCompNext    = localDofIndexPresNext    + ic + 1;
        localIndex const localDofIndexCompCurrent = localDofIndexPresCurrent + ic + 1;

        dofColIndices[localDofIndexCompNext]    = offsetNext    + ColOffset::DPRES + ic + 1;
        dofColIndices[localDofIndexCompCurrent] = offsetCurrent + ColOffset::DPRES + ic + 1;
          
        localFluxJacobian[localDofIndexCompNext]    = - dAvgDensity_dCompNext[ic]    * gravD;  
        localFluxJacobian[localDofIndexCompCurrent] = - dAvgDensity_dCompCurrent[ic] * gravD;
      }
          
      // TODO: add friction and acceleration terms
        
      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     &localFlux );
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     integer_conversion<int>( 2 * resNDOF ), dofColIndices.data(),
                                     localFluxJacobian.data() );
    }
  });
}
  
void CompositionalMultiphaseWell::FormControlEquation( DomainPartition * const domain,
                                                       Epetra_FECrsMatrix * const jacobian,
                                                       Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom 
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // TODO: check that the first connection is on my rank
    // if not, there is nothing to do here
    
    // for now, we assume here that the first well element is on this MPI rank
    localIndex const iwelemControl = 0;

    // get well control
    Well::Control const control = well->getControl();

    // BHP control
    if (control == Well::Control::BHP)
    {
      // get primary variables on well elements
      arrayView1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      arrayView1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // get the pressure and compute normalizer
      real64 const currentBHP = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];
      real64 const targetBHP  = well->getTargetBHP();
      real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
                              ? 1.0 / targetBHP
                              : 1.0;

      // control equation is a normalized difference between current pressure and target pressure
      real64 const controlEqn = ( currentBHP - targetBHP ) * normalizer;
      real64 const dControlEqn_dPres = normalizer;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + ColOffset::DPRES;

      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     &controlEqn );
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex, 
                                     1, &dofColIndex, 
                                     &dControlEqn_dPres );
    }
    else if (control == Well::Control::LIQUIDRATE) // liquid rate control
    {
      localIndex const NC = m_numComponents;
      
      // get a reference to the primary variables on well elements
      arrayView1d<real64 const> const & connRate  =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

      arrayView1d<real64 const> const & dConnRate =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );
      
      // get rates and compute normalizer
      real64 const currentConnRate = connRate[iwelemControl] + dConnRate[iwelemControl];
      real64 const targetConnRate  = well->getTargetRate();
      real64 const normalizer      = targetConnRate > std::numeric_limits<real64>::min()
                                   ? 1.0 / ( 1e-2 * targetConnRate ) // hard-coded value comes from AD-GPRS
                                   : 1.0;

      // control equation is a normalized difference between current rate and target rate
      real64 const controlEqn = ( currentConnRate - targetConnRate ) * normalizer;
      real64 const dControlEqn_dRate = normalizer;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + NC + 1;

      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     &controlEqn );
      
      jacobian->SumIntoGlobalValues( 1, &eqnRowIndex,
                                     1, &dofColIndex,
                                     &dControlEqn_dRate );

    }
    else
    {
      GEOS_ERROR_IF( (control != Well::Control::BHP) && (control != Well::Control::LIQUIDRATE),
                     "Phase rate contraints for CompositionalMultiphaseWell will be implemented later" );
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
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & wellElemPressure  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64> const & wellElemGlobalCompDensity  =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    arrayView2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64> const & connRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

    arrayView1d<real64> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );   

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        wellElemGlobalCompDensity[iwelem][ic] += dWellElemGlobalCompDensity[iwelem][ic];
      }
      connRate[iwelem] += dConnRate[iwelem];
    }

    RecordWellData( well );
  });  
}

void CompositionalMultiphaseWell::ComputeAllPerforationRates( Well * well )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS; 
  localIndex constexpr maxNumDof  = maxNumComp + 1; // here, dofs are 1 pressure, NC comp densities (reservoir and well), no connection rate needed

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

  // get the degrees of freedom and depth
  arrayView1d<globalIndex const> const & wellElemDofNumber =
    wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

  arrayView1d<real64 const> const & wellElemGravDepth =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );
    
  // get well primary variables on well elements
  arrayView1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & dWellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  // get well secondary variables on well elements
  arrayView2d<real64 const> const & wellElemPhaseVolFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );
  
  arrayView2d<real64 const> const & dWellElemPhaseVolFrac_dPres =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64 const> const & dWellElemPhaseVolFrac_dComp =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d<real64 const> const & dWellElemCompFrac_dCompDens =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // get well constitutive data
  MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );
  RelativePermeabilityBase const * relPerm = GetConstitutiveModel<RelativePermeabilityBase>( wellElementSubRegion, m_relPermName );
    
  arrayView3d<real64 const> const & wellElemPhaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dWellElemPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );
  
  arrayView4d<real64 const> const & dWellElemPhaseDens_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  arrayView3d<real64 const> const & wellElemPhaseVisc =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString );

  arrayView3d<real64 const> const & dWellElemPhaseVisc_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

  arrayView4d<real64 const> const & dWellElemPhaseVisc_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString );
    
  arrayView4d<real64 const> const & wellElemPhaseCompFrac =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString );

  arrayView4d<real64 const> const & dWellElemPhaseCompFrac_dPres =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString );

  arrayView5d<real64 const> const & dWellElemPhaseCompFrac_dComp =
    fluid->getReference<array5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString );

  arrayView3d<real64 const> const & wellElemPhaseRelPerm =
    relPerm->getReference<array3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString );

  arrayView4d<real64> const & dWellElemPhaseRelPerm_dPhaseVolFrac =
    relPerm->getReference<array4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString );

  // get well variables on perforations
  arrayView1d<real64 const> const & perfGravDepth =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::gravityDepthString );

  arrayView1d<localIndex const> const & perfWellElemIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d<real64 const> const & perfTransmissibility =
    perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::transmissibilityString );

  arrayView2d<real64> const & compPerfRate =
    perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString );

  arrayView3d<real64> const & dCompPerfRate_dPres =
    perforationData->getReference<array3d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString );

  arrayView4d<real64> const & dCompPerfRate_dComp =
    perforationData->getReference<array4d<real64>>( viewKeyStruct::dCompPerforationRate_dCompString );

  // get the element region, subregion, index
  arrayView1d<localIndex const> const & resElementRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

  arrayView1d<localIndex const> const & resElementSubRegion =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

  arrayView1d<localIndex const> const & resElementIndex =
    perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

  // local working variables and arrays
  stackArray1d<long long, 2 * maxNumComp> eqnRowIndices( 2 * NC );
  stackArray1d<long long, 2 * maxNumDof>  dofColIndices( 2 * resNDOF );

  stackArray1d<double, 2 * maxNumComp>                 localFlux( 2 * NC );
  stackArray2d<double, 2 * maxNumComp * 2 * maxNumDof> localFluxJacobian( 2 * NC, 2 * resNDOF );
     
  stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

  stackArray1d<real64, 2>              pressure( 2 );
  stackArray1d<real64, 2>              dPressure_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dPressure_dC( 2, NC );
      
  stackArray1d<real64, 2>              dPhaseFlux_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dPhaseFlux_dC( 2, NC );

  stackArray1d<real64, maxNumComp> dVolFrac_dC( NC );
  stackArray1d<real64, maxNumComp> dRelPerm_dC( NC );
  stackArray1d<real64, maxNumComp> dDens_dC( NC );
  stackArray1d<real64, maxNumComp> dVisc_dC( NC );

  stackArray1d<real64, maxNumComp> dWellElemAvgDensity_dC( NC );
  stackArray1d<real64, maxNumComp> dResTotalMobility_dC( NC );

  stackArray2d<real64, 2 * maxNumComp>              phaseCompFrac( 2, NC );
  stackArray2d<real64, 2 * maxNumComp>              dPhaseCompFrac_dP( 2, NC );
  stackArray3d<real64, 2 * maxNumComp * maxNumComp> dPhaseCompFrac_dC( 2, NC, NC );
    
  stackArray1d<real64, 2>              mobility( 2 );
  stackArray2d<real64, 4>              dMobility_dP( 2, 2 );
  stackArray3d<real64, 4 * maxNumComp> dMobility_dC( 2, 2, NC );

  stackArray1d<real64, 2>              dPotDiff_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dPotDiff_dC( 2, NC );

  stackArray1d<real64, 2> multiplier( 2 );

  // loop over the perforations to compute the perforation rates
  for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
  {
    eqnRowIndices = -1;
    dofColIndices = -1;

    localFlux = 0.0;
    localFluxJacobian = 0.0;
     
    dPhaseCompFrac_dCompDens = 0;

    // TODO: find a more compact way to do that
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      for (localIndex ke = 0; ke < 2; ++ke)
      {
        compPerfRate[iperf][ic] = 0.0;
        dCompPerfRate_dPres[iperf][ke][ic] = 0.0;
        for (localIndex jc = 0; jc < NC; ++jc)
          dCompPerfRate_dComp[iperf][ke][ic][jc] = 0.0;
      }
    }

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];
    // get the index of the well elem
    localIndex const iwelem = perfWellElemIndex[iperf]; 

    // 1) compute wellElemAvgDensity in well element and resTotalMobility in reservoir
    // TODO: move this to a different function

    real64 wellElemAvgDensity = 0.0;
    real64 resTotalMobility   = 0.0;
    real64 dWellElemAvgDensity_dP = 0.0;
    real64 dResTotalMobility_dP   = 0.0;
    dWellElemAvgDensity_dC = 0.0;
    dResTotalMobility_dC = 0.0;
    
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // increment avg density
      // this will be used to compute the gravity head
      wellElemAvgDensity += wellElemPhaseDens[iwelem][0][ip];      
      dWellElemAvgDensity_dP  += dWellElemPhaseDens_dPres[iwelem][0][ip];         
      
      applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseDens_dComp[iwelem][0][ip], dDens_dC );
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dWellElemAvgDensity_dC[ic] += dDens_dC[ic];
      }

      // increment total mobility
      // this will be used to compute the mobility of the injector
      
      // viscosity
      real64 const resViscosity = resPhaseVisc[er][esr][m_resFluidIndex][ei][0][ip];
      real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][m_resFluidIndex][ei][0][ip];
      applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], dResPhaseVisc_dComp[er][esr][m_resFluidIndex][ei][0][ip], dVisc_dC );
      
      //relative permeability
      real64 const resRelPerm = resPhaseRelPerm[er][esr][m_resRelPermIndex][ei][0][ip];
      real64 dResRelPerm_dP = 0.0;
      dRelPerm_dC = 0.0;
      for (localIndex jp = 0; jp < NP; ++jp)
      {
        real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[er][esr][m_resRelPermIndex][ei][0][ip][jp];
        dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[er][esr][ei][jp];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[er][esr][ei][jp][jc];
        }
      }

      resTotalMobility     += resRelPerm / resViscosity;
      dResTotalMobility_dP += ( dResRelPerm_dP * resViscosity - resRelPerm * dResVisc_dP )
        / ( resViscosity * resViscosity );
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dResTotalMobility_dC[ic] += ( dRelPerm_dC[ic] * resViscosity - resRelPerm * dVisc_dC[ic] )
          / ( resViscosity * resViscosity );
      }
    }
        
    // loop over phases, compute and upwind phase flux
    // and sum contributions to each component's perforation rate
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // clear working arrays
      real64 potDiff = 0.0;
      dPotDiff_dP = 0.0;
      dPotDiff_dC = 0.0;

      real64 phaseFlux = 0.0;
      dPhaseFlux_dP = 0.0;
      dPhaseFlux_dC = 0.0;

      pressure = 0.0;
      dPressure_dP = 0.0;
      dPressure_dC = 0.0;

      phaseCompFrac = 0.0;
      dPhaseCompFrac_dP = 0.0;
      dPhaseCompFrac_dC = 0.0;

      mobility = 0.0;
      dMobility_dP = 0.0;
      dMobility_dC = 0.0;

      dDens_dC = 0.0;
      dVisc_dC = 0.0;
      dVolFrac_dC = 0.0;

      multiplier = 0.0;
      
      // 2) copy the variables from the reservoir and well element

      // TODO: step 2) is not necessary, we should first do the upwinding and then only compute the necessary quantities
      
      // a) get reservoir variables

      // pressure
      pressure[SubRegionTag::RES]  = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
      dPressure_dP[SubRegionTag::RES] = 1.0;
      
      // density
      real64 const resDensity  = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip];
      real64 const dResDens_dP = dResPhaseDens_dPres[er][esr][m_resFluidIndex][ei][0][ip];
      applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], dResPhaseDens_dComp[er][esr][m_resFluidIndex][ei][0][ip], dDens_dC );
      
      // viscosity
      real64 const resViscosity = resPhaseVisc[er][esr][m_resFluidIndex][ei][0][ip];
      real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][m_resFluidIndex][ei][0][ip];
      applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], dResPhaseVisc_dComp[er][esr][m_resFluidIndex][ei][0][ip], dVisc_dC );
      
      //relative permeability
      real64 const resRelPerm = resPhaseRelPerm[er][esr][m_resRelPermIndex][ei][0][ip];
      real64 dResRelPerm_dP = 0.0;
      dRelPerm_dC = 0.0;
      for (localIndex jp = 0; jp < NP; ++jp)
      {
        real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[er][esr][m_resRelPermIndex][ei][0][ip][jp];
        dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[er][esr][ei][jp];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[er][esr][ei][jp][jc];
        }
      }

      // component phase fractions and pressure derivatives
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        phaseCompFrac[SubRegionTag::RES][ic] = resPhaseCompFrac[er][esr][m_resFluidIndex][ei][0][ip][ic];
        dPhaseCompFrac_dP[SubRegionTag::RES][ic] = dResPhaseCompFrac_dPres[er][esr][m_resFluidIndex][ei][0][ip][ic];

        // compositional derivatives
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseCompFrac_dC[SubRegionTag::RES][ic][jc] = dResPhaseCompFrac_dComp[er][esr][m_resFluidIndex][ei][0][ip][ic][jc];
        }
      }
      
      // mobility and pressure derivative
      mobility[SubRegionTag::RES] = resRelPerm * resDensity / resViscosity;
      dMobility_dP[SubRegionTag::RES][SubRegionTag::RES] = dResRelPerm_dP * resDensity / resViscosity
        + mobility[SubRegionTag::RES] * (dResDens_dP / resDensity - dResVisc_dP / resViscosity);

      // compositional derivatives
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dMobility_dC[SubRegionTag::RES][SubRegionTag::RES][jc] = dRelPerm_dC[jc] * resDensity / resViscosity
          + mobility[SubRegionTag::RES] * (dDens_dC[jc] / resDensity - dVisc_dC[jc] / resViscosity);
      }

      multiplier[SubRegionTag::RES] = 1;

      
      // b) get well variables

      // pressure
      pressure[SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
      dPressure_dP[SubRegionTag::WELL] = 1.0;

      if (m_gravityFlag)
      {
        real64 const gravD = ( perfGravDepth[iperf] - wellElemGravDepth[iwelem] );
        pressure[SubRegionTag::WELL]  += wellElemAvgDensity * gravD;
        dPressure_dP[SubRegionTag::WELL] += dWellElemAvgDensity_dP * gravD;
        for (localIndex ic = 0; ic < NC; ++ic)
          dPressure_dC[SubRegionTag::WELL][ic] += dWellElemAvgDensity_dC[ic] * gravD;
      }
      
      // density
      real64 const wellDensity  = wellElemPhaseDens[iwelem][0][ip];
      real64 const dWellDens_dP = dWellElemPhaseDens_dPres[iwelem][0][ip];
      applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseDens_dComp[iwelem][0][ip], dDens_dC );
      
      // viscosity
      real64 const wellViscosity = wellElemPhaseVisc[iwelem][0][ip];
      real64 const dWellVisc_dP  = dWellElemPhaseVisc_dPres[iwelem][0][ip];
      applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseVisc_dComp[iwelem][0][ip], dVisc_dC );

      // volume fraction
      real64 const wellVolFrac = wellElemPhaseVolFrac[iwelem][ip];
      real64 const dWellVolFrac_dP = dWellElemPhaseVolFrac_dPres[iwelem][ip];
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dVolFrac_dC[ic] = dWellElemPhaseVolFrac_dComp[iwelem][ip][ic];
      }
      
      // component phase fractions and pressure derivatives
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        phaseCompFrac[SubRegionTag::WELL][ic] = wellElemPhaseCompFrac[iwelem][0][ip][ic];
        dPhaseCompFrac_dP[SubRegionTag::WELL][ic] = dWellElemPhaseCompFrac_dPres[iwelem][0][ip][ic];

        // compositional derivatives
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseCompFrac_dC[SubRegionTag::WELL][ic][jc] = dWellElemPhaseCompFrac_dComp[iwelem][0][ip][ic][jc];
        }
      }
      
      // mobility and pressure derivative
      mobility[SubRegionTag::WELL] = wellVolFrac * wellDensity * resTotalMobility;
      dMobility_dP[SubRegionTag::WELL][SubRegionTag::RES]  = wellVolFrac * wellDensity * dResTotalMobility_dP;
      dMobility_dP[SubRegionTag::WELL][SubRegionTag::WELL] = ( dWellVolFrac_dP * wellDensity + wellVolFrac * dWellDens_dP )
        * resTotalMobility;

      // compositional derivatives
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dMobility_dC[SubRegionTag::WELL][SubRegionTag::RES][ic] = wellVolFrac * wellDensity * dResTotalMobility_dC[ic];
        dMobility_dC[SubRegionTag::WELL][SubRegionTag::WELL][ic] = ( dVolFrac_dC[ic] * wellDensity + wellVolFrac * dDens_dC[ic] )
                                                                    *  resTotalMobility;
      }
      
      multiplier[SubRegionTag::WELL] = -1;

      // 3) calculation of the potential difference
                
      // get transmissibility at the interface
      real64 const trans = perfTransmissibility[iperf]; 
        
      // compute potential difference
      for (localIndex i = 0; i < 2; ++i)
      {
        potDiff += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
        dPotDiff_dP[i] += multiplier[i] * trans * dPressure_dP[i];

        for (localIndex ic = 0; ic < NC; ++ic)
          dPotDiff_dC[i][ic] += multiplier[i] * trans * dPressure_dC[i][ic];
      }

      
      // 4) upwinding 

      // choose upstream cell
      localIndex const k_up = (potDiff >= 0) ? SubRegionTag::RES : SubRegionTag::WELL;

      // skip the phase flux if phase not present or immobile upstream
      if (std::fabs(mobility[k_up]) < 1e-20) // TODO better constant
      {
        continue;
      }

      // 5) computation of the phase flux
      
      // pressure gradient depends on all points in the stencil
      for (localIndex ke = 0; ke < 2; ++ke)
      {
        dPhaseFlux_dP[ke] += dPotDiff_dP[ke];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] += dPotDiff_dC[ke][jc];
        }
      }
      
      // compute the phase flux and derivatives using upstream cell mobility
      phaseFlux = mobility[k_up] * potDiff;
      for (localIndex ke = 0; ke < 2; ++ke)
      {
        dPhaseFlux_dP[ke] *= mobility[k_up];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] *= mobility[k_up];
        }
      }

      // add contribution from upstream cell mobility derivatives
      for (localIndex ke = 0; ke < 2; ++ke)
      {
        dPhaseFlux_dP[ke] += dMobility_dP[k_up][ke] * potDiff;
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] += dMobility_dC[k_up][ke][jc] * potDiff;
        }
      }
      
      // 6) increment the component fluxes

      // compute component fluxes and derivatives using upstream cell composition
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        compPerfRate[iperf][ic] += phaseFlux * phaseCompFrac[k_up][ic];

        // derivatives stemming from phase flux
        for (localIndex ke = 0; ke < 2; ++ke)
        {
          dCompPerfRate_dPres[iperf][ke][ic] += dPhaseFlux_dP[ke] * phaseCompFrac[k_up][ic];
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dCompPerfRate_dComp[iperf][ke][ic][jc] += dPhaseFlux_dC[ke][jc] * phaseCompFrac[k_up][ic];
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        dCompPerfRate_dPres[iperf][k_up][ic] += phaseFlux * dPhaseCompFrac_dP[k_up][ic];

        // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component densities
        if (k_up == SubRegionTag::WELL)
          applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dPhaseCompFrac_dC[k_up][ic], dPhaseCompFrac_dCompDens );
        else
          applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], dPhaseCompFrac_dC[k_up][ic], dPhaseCompFrac_dCompDens );
        
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompPerfRate_dComp[iperf][k_up][ic][jc] += phaseFlux * dPhaseCompFrac_dCompDens[jc];
        }
      }
    }
  }
}

void CompositionalMultiphaseWell::RecordWellData( Well * well )
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;

  localIndex const NP = m_numPhases;
  localIndex const NC = m_numComponents;
  
  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
  PerforationData const * perforationData = well->getPerforations();
  
  arrayView1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView2d<real64> const & wellElemPhaseVolFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

  arrayView1d<real64 const> const & connRate =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

  // get the index of the next well element
  arrayView1d<localIndex const> const & nextWellElemIndex =
    wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    
  
  arrayView2d<real64 const> const & compPerfRate =
    perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString );
  
  // get well secondary variables on well elements
  MultiFluidBase const * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );

  arrayView3d<real64 const> const & wellElemPhaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView4d<real64 const> const & wellElemPhaseCompFrac =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString );

  
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
      std::cout << "Mass rate at perforation #" << iperf << " for component #" << ic << ": " << compPerfRate[iperf][ic] << std::endl;
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
  
  for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
  {
    if (well->getControl() == Well::Control::BHP)
      std::cout << "Volumetric rate at connection #" << iwelem << ": " << connRate[iwelem] << std::endl;
    else
    {
      std::cout << "Volumetric at connection #" << iwelem << ": " << connRate[iwelem]
                << ", target rate : " << well->getTargetRate() << std::endl;
    }
  
    if (iwelem > 0)
      continue;
    
    stackArray1d<real64, maxNumComp> compStream( NC );
    compStream  = 0.0;

    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const phaseDens = wellElemPhaseDens[iwelem][0][ip];
      real64 const volFrac = wellElemPhaseVolFrac[iwelem][ip];
      real64 const phaseStream = volFrac * phaseDens * connRate[iwelem];
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const phaseCompFrac = wellElemPhaseCompFrac[iwelem][0][ip][ic];
        compStream[ic] += phaseCompFrac * phaseStream;
      }
    }

    if (well->getType() == Well::Type::INJECTOR)
    {
      for (localIndex ic = 0; ic < NC; ++ic)
        std::cout << "Mass injection rate for component #" << ic
                  << " at connection #" << iwelem
                  << ": " << compStream[ic] << std::endl;
    }
    else
    {
      for (localIndex ic = 0; ic < NC; ++ic)
        std::cout << "Mass production rate for component #" << ic
                  << " at connection #" << iwelem
                  << ": " << compStream[ic] << std::endl;
    }
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
      arrayView1d<real64 const> const & wellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

      arrayView1d<real64 const> const & dWellElemPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      // again we assume here that the first well element is on this MPI rank
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
        well->setControl( Well::Control::LIQUIDRATE );
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
