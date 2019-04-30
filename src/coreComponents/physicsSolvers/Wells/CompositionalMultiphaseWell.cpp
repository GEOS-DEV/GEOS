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

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"
#include "dataRepository/ManagedGroup.hpp"
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
  m_numComponents( 0 )
{
  this->RegisterViewWrapper( viewKeyStruct::temperatureString, &m_temperature, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Temperature");

  this->RegisterViewWrapper( viewKeyStruct::useMassFlagString, &m_useMass, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Use mass formulation instead of molar");

  this->RegisterViewWrapper( viewKeyStruct::resRelPermNameString,  &m_resRelPermName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the relative permeability constitutive model to use");

  this->RegisterViewWrapper( viewKeyStruct::resRelPermIndexString, &m_resRelPermIndex, false );

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

  wellManager->forSubGroups<Well>( [&] ( Well * const well ) -> void
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

    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::mixtureDensityString );
    wellElementSubRegion->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dMixtureDensity_dPressureString );
    wellElementSubRegion->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

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

  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  RelativePermeabilityBase const * relPerm = cm->GetConstitituveRelation<RelativePermeabilityBase>( m_resRelPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Relative permeability model " + m_resRelPermName + " not found" );
  m_resRelPermIndex = relPerm->getIndexInParent();

  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * const well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );

    // TODO: put this out of the loop
    m_numPhases     = fluid->numFluidPhases();
    m_numComponents = fluid->numFluidComponents();

    // check well injection stream for injectors
    if ( well->getType() == Well::Type::INJECTOR )
    { 
      real64 compFracSum = 0;
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        real64 const compFrac = well->getInjectionStream( ic );
        GEOS_ERROR_IF( compFrac < 0.0 || compFrac > 1.0, 
                       "Invalid injection stream for well " << well->getName() );
        compFracSum += compFrac;
      }
      GEOS_ERROR_IF( compFracSum < 1.0 - std::numeric_limits<real64>::epsilon() || 
                     compFracSum > 1.0 + std::numeric_limits<real64>::epsilon(), 
                     "Invalid injection stream for well " << well->getName() );
    }

    m_numDofPerElement = m_numComponents + 2; // 1 pressure + NC compositions + 1 connectionRate
 
    ResizeFields( well );
  });
}

void CompositionalMultiphaseWell::ResizeFields( Well * const well )
{
  WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
  PerforationData * const perforationData = well->getPerforations();

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompDensityString).resizeDimension<1>(NC);
  wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::deltaGlobalCompDensityString).resizeDimension<1>(NC);

  wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::globalCompFractionString).resizeDimension<1>(NC);
  wellElementSubRegion->getReference<array3d<real64>>(viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString).resizeDimension<1,2>(NC, NC);

  wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::phaseVolumeFractionString).resizeDimension<1>(NP);
  wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dPressureString).resizeDimension<1>(NP);
  wellElementSubRegion->getReference<array3d<real64>>(viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString).resizeDimension<1,2>(NP, NC);

  wellElementSubRegion->getReference<array2d<real64>>(viewKeyStruct::dMixtureDensity_dGlobalCompDensityString).resizeDimension<1>(NC);

  perforationData->getReference<array2d<real64>>( viewKeyStruct::compPerforationRateString ).resizeDimension<1>( NC );
  perforationData->getReference<array3d<real64>>( viewKeyStruct::dCompPerforationRate_dPresString ).resizeDimension<1,2>( 2, NC );
  perforationData->getReference<array4d<real64>>( viewKeyStruct::dCompPerforationRate_dCompString ).resizeDimension<1,2,3>( 2, NC, NC );

}

void CompositionalMultiphaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  WellManager * const wellManager = domain->getWellManager();

  // set mass fraction flag on subregion models
  wellManager->forSubGroups<Well>( [&] ( Well * const well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();

    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );
    fluid->setMassFlag( static_cast<bool>(m_useMass) );
  });
}

void CompositionalMultiphaseWell::UpdateComponentFraction( Well const * const well )
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

void CompositionalMultiphaseWell::UpdatePhaseVolumeFraction( Well const * const well )
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

void CompositionalMultiphaseWell::UpdateFluidModel( Well * const well )
{
  WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
  MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );

  arrayView1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & dWellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  arrayView2d<real64 const> const & wellElemCompFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  forall_in_range<RAJA::seq_exec>( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
  {
    fluid->PointUpdate( wellElemPressure[iwelem] + dWellElemPressure[iwelem],
                        m_temperature,
                        wellElemCompFrac[iwelem],
                        iwelem,
                        0 );
  });
}
  
void CompositionalMultiphaseWell::UpdateMixtureDensity( Well const * const well )
{
  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  // get well secondary variables on well elements
  arrayView1d<real64> const & wellElemMixtureDensity =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureDensityString );
  
  arrayView1d<real64> const & dWellElemMixtureDensity_dPres =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::dMixtureDensity_dPressureString );

  arrayView2d<real64> const & dWellElemMixtureDensity_dComp =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

  arrayView2d<real64 const> const & wellElemPhaseVolFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );
  
  arrayView2d<real64 const> const & dWellElemPhaseVolFrac_dPres =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d<real64 const> const & dWellElemPhaseVolFrac_dComp =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d<real64 const> const & dWellElemCompFrac_dCompDens =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // get constitutive data
  MultiFluidBase const * const fluid = GetConstitutiveModel<MultiFluidBase>( wellElementSubRegion, m_fluidName );
   
  arrayView3d<real64 const> const & wellElemPhaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

  arrayView3d<real64 const> const & dWellElemPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  arrayView4d<real64 const> const & dWellElemPhaseDens_dComp =
    fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

  stackArray1d<real64, maxNumComp> dDens_dC( NC );

  for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
  {

    // reset to zero
    wellElemMixtureDensity[iwelem]        = 0.0;
    dWellElemMixtureDensity_dPres[iwelem] = 0.0;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      dWellElemMixtureDensity_dComp[iwelem][ic] = 0.0;
    }

    // increment mixture velocity
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      wellElemMixtureDensity[iwelem] += wellElemPhaseVolFrac[iwelem][ip] 
                                      * wellElemPhaseDens[iwelem][0][ip];      
      dWellElemMixtureDensity_dPres[iwelem] += dWellElemPhaseVolFrac_dPres[iwelem][ip] * wellElemPhaseDens[iwelem][0][ip]
                                             + wellElemPhaseVolFrac[iwelem][ip] * dWellElemPhaseDens_dPres[iwelem][0][ip];
  
      dDens_dC = 0.0;
    
      applyChainRule( NC, dWellElemCompFrac_dCompDens[iwelem], dWellElemPhaseDens_dComp[iwelem][0][ip], dDens_dC );
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dWellElemMixtureDensity_dComp[iwelem][ic] += dWellElemPhaseVolFrac_dComp[iwelem][ip][ic] * wellElemPhaseDens[iwelem][0][ip]
                                                   + wellElemPhaseVolFrac[iwelem][ip] * dDens_dC[ic];
      }
    }
  }
}

void CompositionalMultiphaseWell::UpdateState( Well * const well )
{
  UpdateComponentFraction( well );
  UpdateFluidModel( well );
  UpdatePhaseVolumeFraction( well );
  UpdateMixtureDensity( well );
}

void CompositionalMultiphaseWell::InitializeWells( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & resPressure = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & resCompDens = m_resGlobalCompDensity;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & resPhaseVolFrac = m_resPhaseVolFrac;
  
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseDens    = m_resPhaseDens;

  wellManager->forSubGroups<Well>( [&] ( Well * const well ) -> void
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

    arrayView1d<real64 const> const & wellElemGravDepth =
      wellElementSubRegion->getReference<array1d<real64>>( WellElementSubRegion::viewKeyStruct::gravityDepthString );
    
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

    arrayView2d<real64 const> const & totalDens =
      fluid->getReference<array2d<real64>>( MultiFluidBase::viewKeyStruct::totalDensityString );
    
    // 1) Loop over all perforations to compute an average mixture density 
    //    and component fraction
    real64 avgTotalDensity   = 0.0;
    real64 avgMixtureDensity = 0.0;
    stackArray1d<real64, maxNumComp> avgCompFrac( NC );
    
    avgCompFrac = 0.0;
    
    // define a reservoir pressure used for initialization
    real64 resPres = well->getTargetBHP();

    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];
 
      // save min pressure for producer
      if ( well->getType() == Well::Type::PRODUCER &&
           resPres > resPressure[er][esr][ei] )
      { 
        resPres = resPressure[er][esr][ei];
      }
      // save max pressure for injector
      else if ( well->getType() == Well::Type::INJECTOR &&
                resPres < resPressure[er][esr][ei] )
      { 
        resPres = resPressure[er][esr][ei];
      }

      // increment the average mixture density
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        real64 const resDensity = resPhaseDens[er][esr][m_resFluidIndex][ei][0][ip];
        real64 const resVolFrac = resPhaseVolFrac[er][esr][ei][ip];
        avgMixtureDensity += resVolFrac * resDensity;
      }

      // increment the average total density
      real64 perfTotalDensity = 0.0;
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        perfTotalDensity += resCompDens[er][esr][ei][ic];
      }
      avgTotalDensity += perfTotalDensity;

      // increment the average component fraction
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        avgCompFrac[ic] += resCompDens[er][esr][ei][ic] / perfTotalDensity;
      }
    }

    // compute average densities
    avgMixtureDensity /= perforationData->numPerforationsGlobal();
    avgTotalDensity   /= perforationData->numPerforationsGlobal();

    // compute average component fraction
    if ( well->getType() == Well::Type::PRODUCER )
    {
      // use average comp frac from reservoir
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        avgCompFrac[ic] /= perforationData->numPerforationsGlobal(); 
      }  
    }
    else // injector
    {
      // use average comp frac from XML file
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        avgCompFrac[ic] = well->getInjectionStream( ic ); 
      }
    }

    // set the global component fractions to avgCompFrac
    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        wellElemCompFrac[iwelem][ic] = avgCompFrac[ic];
      }
    });

    // get the reference data for this well
    localIndex const iwelemControl = well->getReferenceWellElementIndex();
    real64 const gravDepthControl = wellElemGravDepth[iwelemControl];

    // 2) Initialize the reference pressure
    real64 const targetBHP = well->getTargetBHP();
    if (well->getControl() == Well::Control::BHP)
    {
      // if pressure constraint, set the ref pressure at the constraint
      wellElemPressure[iwelemControl] = targetBHP;
    }
    else // rate control
    {
      // if rate constraint, set the ref pressure slightly 
      // above/below the target pressure depending on well type
      wellElemPressure[iwelemControl] = (well->getType() == Well::Type::PRODUCER)
        ? 0.5 * resPres // hard-coded values come from personal communication with Hui 
        : 2.0 * resPres; 
    }

    // 3) Estimate the pressures in the well elements using this avgDensity
    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] = wellElemPressure[iwelemControl]
        + ( m_gravityFlag 
          ? avgMixtureDensity * ( wellElemGravDepth[iwelem] - gravDepthControl ) 
          : 0 );
    });
    
    // 4) Back calculate component densities
    forall_in_range<RAJA::seq_exec>( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      fluid->PointUpdate( wellElemPressure[iwelem],
                          m_temperature,
                          wellElemCompFrac[iwelem],
                          iwelem,
                          0 );
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        wellElemCompDens[iwelem][ic] = avgCompFrac[ic] * totalDens[iwelem][0];
      }
    });

    // 5) Recompute all the pressure-dependent properties
    UpdateState( well );
    
    // 6) Estimate the connection rates based on the min/max pressure
    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      real64 const targetRate = well->getTargetRate();
      if (well->getControl() == Well::Control::BHP)
      {
        // if BHP constraint set rate below the absolute max rate 
        // with the appropriate sign (negative for prod, positive for inj)
        connRate[iwelem] = (well->getType() == Well::Type::PRODUCER)
          ? std::max( 0.1 * targetRate, - 1e3 ) // hard-coded values come from personal communication with Hui
          : std::min( 0.1 * targetRate,   1e3 );
      }
      else
      {
        connRate[iwelem] = targetRate;
      }
    });
  });
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

  // find the first local row on this partition, 
  // and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<thisMpiProcess)
    {
      firstLocalRow += gather[p];
    }
  }

  // TODO: double check this for multiple MPI processes
  
  // get the well information
  WellManager const * const wellManager = domain->getWellManager();

  localIndex localCount = 0;
  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
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

  // save these numbers (reused to compute the well elem offsets)
  m_firstWellElemDofNumber = firstWellElemDofNumber;
  m_numDofPerResElement    = numDofPerResElement;

  // set the number of degrees of freedom per element 
  // on both sides (reservoir and well)  
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS; 
  localIndex constexpr maxNumDof  = maxNumComp + 1; // dofs are 1 pressure, NC comp densities (reservoir and well), no connection rate needed

  // reservoir dofs are 1 pressure and NC comp densities
  localIndex const resNDOF  = numDofPerResElement;
  // well dofs are 1 pressure, NC comp densities, 1 connection rate
  localIndex const wellNDOF = numDofPerElement(); 

  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
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
    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {

      // get next well element index
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      // check if this is not an entry or exit
      if (iwelemNext >= 0)
      {
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
  });
  
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
  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get a reference to the degree-of-freedom numbers
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );

    // get a reference to the primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView1d<real64 const> const & connRate  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

    arrayView1d<real64 const> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );
    
    // get the info stored on well elements
    arrayView2d<real64 const> const & wellElemCompFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

    arrayView3d<real64 const> const & dWellElemCompFrac_dCompDens =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    arrayView1d<real64 const> const & wellElemMixtureDensity =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureDensityString );
  
    arrayView1d<real64 const> const & dWellElemMixtureDensity_dPres =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::dMixtureDensity_dPressureString );

    arrayView2d<real64 const> const & dWellElemMixtureDensity_dCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );
    
    // create local work arrays
    stackArray1d<real64, maxNumComp>              compFracUp( NC );
    stackArray1d<real64, maxNumComp>              dCompFrac_dPresUp( NC );
    stackArray2d<real64, maxNumComp * maxNumComp> dCompFrac_dCompDensUp( NC, NC );

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
    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
 
      // Step 1) decide the upwind well element

      /*  currentConnRate < 0 flow from iwelem to iwelemNext
       *  currentConnRate > 0 flow from iwelemNext to iwelem
       *  With this convention, currentConnRate < 0 at the last connection for a producer
       *                        currentConnRate > 0 at the last connection for a injector
       */

      localIndex const iwelemNext = nextWellElemIndex[iwelem];
      real64 const currentConnRate = connRate[iwelem] + dConnRate[iwelem];
      localIndex iwelemUp = -1;
          
      if (iwelemNext < 0 && well->getType() == Well::Type::INJECTOR) // exit connection, injector
      {
        // we still need to define iwelemUp for Jacobian assembly
        iwelemUp = iwelem;

        // just copy the injection stream into compFrac
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          compFracUp[ic] = well->getInjectionStream( ic );
        }
        // zero out the derivatives wrt composition
        dCompFrac_dCompDensUp = 0.0;

      }
      else 
      {
        // first set iwelemUp to the upstream cell
        if ( (iwelemNext < 0 && well->getType() == Well::Type::PRODUCER) // exit connection, producer
             || currentConnRate < 0) // not an exit connection, iwelem is upstream
        { 
          iwelemUp = iwelem;
        }
        else // not an exit connection, iwelemNext is upstream
        {
          iwelemUp = iwelemNext;
        }

        // copy the vars of iwelemUp into compFrac
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          compFracUp[ic] = wellElemCompFrac[iwelemUp][ic];
          for (localIndex jc = 0; jc < NC; ++jc)
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

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        compFlux[ic]          = compFracUp[ic] * currentConnRate;
        dCompFlux_dRate[ic]   = compFracUp[ic];
        dCompFlux_dPresUp[ic] = 0.0; // none of these quantities depend on pressure
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompFlux_dCompDensUp[ic][jc] = dCompFrac_dCompDensUp[ic][jc] * currentConnRate;
        }
      }

      globalIndex const offsetUp = getElementOffset( wellElemDofNumber[iwelemUp] );
      globalIndex const offsetCurrent = getElementOffset( wellElemDofNumber[iwelem] );

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
          oneSidedFluxJacobian_dPresCompUp[ic][0] = - dt * dCompFlux_dPresUp[ic];
          
          // derivatives with respect to upstream component densities
          for (localIndex jdof = 0; jdof < NC; ++jdof)
          {
            oneSidedFluxJacobian_dPresCompUp[ic][jdof + 1] = - dt * dCompFlux_dCompDensUp[ic][jdof];
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
        localIndex const dRateColOffset = ColOffset::DCOMP + NC; 
        globalIndex const oneSidedDofColIndices_dRate = offsetCurrent + dRateColOffset; 
          
        for (localIndex jdof = 0; jdof < resNDOF; ++jdof)
        {
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
      }
      else // not an exit connection
      {
        globalIndex const offsetNext = getElementOffset( wellElemDofNumber[iwelemNext] );

        // flux terms
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          localFlux[ElemTag::NEXT * NC + ic]    =   dt * compFlux[ic];
          localFlux[ElemTag::CURRENT * NC + ic] = - dt * compFlux[ic];

          // derivative with respect to rate
          localFluxJacobian_dRate[ElemTag::NEXT * NC + ic]    =   dt * dCompFlux_dRate[ic];
          localFluxJacobian_dRate[ElemTag::CURRENT * NC + ic] = - dt * dCompFlux_dRate[ic];
          
          // derivative with respect to upstream pressure
          localFluxJacobian_dPresCompUp[ElemTag::NEXT * NC + ic][0]    =    dt * dCompFlux_dPresUp[ic];
          localFluxJacobian_dPresCompUp[ElemTag::CURRENT * NC + ic][0] =  - dt * dCompFlux_dPresUp[ic];
          
          // derivatives with respect to upstream component densities
          for (localIndex jdof = 0; jdof < NC; ++jdof)
          {
            localFluxJacobian_dPresCompUp[ElemTag::NEXT * NC + ic][jdof + 1] =   
                dt * dCompFlux_dCompDensUp[ic][jdof];
            localFluxJacobian_dPresCompUp[ElemTag::CURRENT * NC + ic][jdof + 1] = 
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
        localIndex const dRateColOffset = ColOffset::DCOMP + NC; 
        globalIndex const dofColIndices_dRate = offsetCurrent + dRateColOffset;
          
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
    });
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

  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    PerforationData const * const perforationData = well->getPerforations();
  
    // get the degrees of freedom
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get the properties on the well element
    arrayView2d<real64 const> const & wellElemPhaseVolFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

    arrayView2d<real64 const> const & dWellElemPhaseVolFrac_dPres =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    arrayView3d<real64 const> const & dWellElemPhaseVolFrac_dComp =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    arrayView1d<real64 const> const & wellElemVolume =
      wellElementSubRegion->getReference<array1d<real64>>( ElementSubRegionBase::viewKeyStruct::elementVolumeString );

    stackArray1d<long long, maxNumDof> localVolBalanceDOF( welemNDOF );
    stackArray1d<double, maxNumDof>    localVolBalanceJacobian( welemNDOF );

    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      localVolBalanceDOF = 0;
      localVolBalanceJacobian = 0;

      // get volume (assume porosity is 1)
      real64 const volume = wellElemVolume[iwelem];

      // get equation/dof indices
      globalIndex const offset = getElementOffset( wellElemDofNumber[iwelem] );
      localIndex const volBalRowOffset = RowOffset::MASSBAL + NC; 
      globalIndex const localVolBalanceEqnIndex = offset + volBalRowOffset;
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
    });
  });
}


void CompositionalMultiphaseWell::AssemblePerforationTerms( DomainPartition * const domain,
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
  
  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
    PerforationData const * const perforationData = well->getPerforations();

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

    // TODO: make this work if the wellElement and the reservoir element are on different ranks

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
      globalIndex const resOffset = resNDOF * resDofNumber[er][esr][ei];
      globalIndex const wellElemOffset = getElementOffset( wellElemDofNumber[iwelem] );

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        eqnRowIndices[SubRegionTag::RES  * NC + ic] = resOffset + ic;
        eqnRowIndices[SubRegionTag::WELL * NC + ic] = wellElemOffset + RowOffset::MASSBAL + ic;
      }
      for (localIndex jdof = 0; jdof < resNDOF; ++jdof)
      {
        dofColIndices[SubRegionTag::RES  * resNDOF + jdof] = resOffset  + jdof;
        dofColIndices[SubRegionTag::WELL * resNDOF + jdof] = wellElemOffset + ColOffset::DPRES + jdof;
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

  real64 residualNorm = 0;
  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers
    arrayView1d<globalIndex const > const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

      localIndex const wellNDOF = numDofPerElement(); 
      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
        int const lid = rowMap->LID(integer_conversion<int>( elemOffset + idof ));
        real64 const val = localResidual[lid];
        residualNorm += val * val;
      }
    }
  });

  return sqrt(residualNorm);
}

bool
CompositionalMultiphaseWell::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
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

  bool isValid = true;

  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers on well elements
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    // get a reference to the primary variables on well elements
    arrayView1d<real64 const> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64 const> const & wellElemCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    arrayView2d<real64 const> const & dWellElemCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );
    
    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

      // pressure 
      int lid = rowMap->LID( integer_conversion<int>( elemOffset + ColOffset::DPRES) );
      real64 const newPres = wellElemPressure[iwelem] + dWellElemPressure[iwelem]
                           + scalingFactor * local_solution[lid];

      if (newPres < 0.0)
      {
        isValid = false;
      }

      // comp densities
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        lid = rowMap->LID(integer_conversion<int>( elemOffset + ic + 1));
        real64 const newDens = wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic] 
                             + scalingFactor * local_solution[lid];
        if (newDens < 0.0)
        {
          isValid = false;
        }
      }
    }
  });  

  return isValid;
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

  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the degree of freedom numbers on well elements, next elem index
    arrayView1d<globalIndex> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );

    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );
    
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );

    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      globalIndex const elemOffset = getElementOffset( wellElemDofNumber[iwelem] );

      // pressure
      int lid = rowMap->LID( integer_conversion<int>( elemOffset ) );
      dWellElemPressure[iwelem] += scalingFactor * local_solution[lid];
        
      // comp densities
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        lid = rowMap->LID( integer_conversion<int>( elemOffset + ic + 1 ) );
        dWellElemGlobalCompDensity[iwelem][ic] += scalingFactor * local_solution[lid];
      }

      // conn rate
      lid = rowMap->LID( integer_conversion<int>( elemOffset + m_numComponents + 1 ) );
      dConnRate[iwelem] += scalingFactor * local_solution[lid];
    });
  });  

  // update properties
  UpdateStateAll( domain );
    
}

void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get next well elem index
    arrayView1d<localIndex const> const & nextWellElemIndex =
      wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );
    
    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );

    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      // extract solution and apply to dP
      dWellElemPressure[iwelem] = 0;
      dConnRate[iwelem] = 0;
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        dWellElemGlobalCompDensity[iwelem][ic] = 0;
      }
    });
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

  m_resPhaseMob =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::phaseMobilityString );

  m_dResPhaseMob_dPres =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dPressureString );

  m_dResPhaseMob_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

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
  
  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

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
    arrayView1d<real64 const> const & wellElemMixtureDensity =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureDensityString );
  
   arrayView1d<real64 const> const & dWellElemMixtureDensity_dPres =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::dMixtureDensity_dPressureString );

    arrayView2d<real64 const> const & dWellElemMixtureDensity_dComp =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

    // local working variables and arrays
    stackArray1d<globalIndex, 2 * maxNumDof > dofColIndices( 2 * resNDOF );
    stackArray1d<real64, 2 * maxNumDof > localFluxJacobian( 2 * resNDOF );

    stackArray1d<real64, maxNumComp> dAvgDensity_dCompCurrent( NC );
    stackArray1d<real64, maxNumComp> dAvgDensity_dCompNext( NC ); 

    // loop over the well elements to compute the pressure relations between well elements
    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      if ( iwelemNext >= 0 )  // if iwelemNext < 0, form control equation, not momentum
      {

        // reset local working variables and arrays
        dofColIndices = -1;
        localFluxJacobian = 0.0;

        // compute the average density at the interface between well elements
        real64 const avgDensity = 0.5 * ( wellElemMixtureDensity[iwelemNext]
                                        + wellElemMixtureDensity[iwelem] );
        real64 const dAvgDensity_dPresNext    = 0.5 * dWellElemMixtureDensity_dPres[iwelemNext];
        real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemMixtureDensity_dPres[iwelem];
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dAvgDensity_dCompNext[ic]    = 0.5 * dWellElemMixtureDensity_dComp[iwelemNext][ic];
          dAvgDensity_dCompCurrent[ic] = 0.5 * dWellElemMixtureDensity_dComp[iwelem][ic];
        }
        
        // compute depth diff times acceleration
        real64 const gravD = ( wellElemGravDepth[iwelemNext] - wellElemGravDepth[iwelem] );

        // compute the current pressure in the two well elements
        real64 const pressureNext    = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];
        real64 const pressureCurrent = wellElemPressure[iwelem] + dWellElemPressure[iwelem];

        // compute a coefficient to normalize the momentum equation
        real64 const targetBHP  = well->getTargetBHP();
        real64 const normalizer = targetBHP > std::numeric_limits<real64>::epsilon()
                                ? 1.0 / targetBHP
                                : 1.0;

        // compute momentum flux and derivatives
        localIndex const localDofIndexPresNext    = ElemTag::NEXT * resNDOF;
        localIndex const localDofIndexPresCurrent = ElemTag::CURRENT * resNDOF;

        globalIndex const offsetNext    = getElementOffset( wellElemDofNumber[iwelemNext] ); 
        globalIndex const offsetCurrent = getElementOffset( wellElemDofNumber[iwelem] );

        globalIndex const eqnRowIndex = offsetCurrent + RowOffset::CONTROL; 

        dofColIndices[localDofIndexPresNext]    = offsetNext    + ColOffset::DPRES;
        dofColIndices[localDofIndexPresCurrent] = offsetCurrent + ColOffset::DPRES;

        real64 const localFlux = ( pressureNext - pressureCurrent - avgDensity * gravD ) * normalizer;
        localFluxJacobian[localDofIndexPresNext]    = ( 1 - dAvgDensity_dPresNext * gravD )    * normalizer;
        localFluxJacobian[localDofIndexPresCurrent] = (-1 - dAvgDensity_dPresCurrent * gravD ) * normalizer;

        for (localIndex ic = 0; ic < NC; ++ic)
        {
          localIndex const localDofIndexCompNext    = localDofIndexPresNext    + ic + 1;
          localIndex const localDofIndexCompCurrent = localDofIndexPresCurrent + ic + 1;

          dofColIndices[localDofIndexCompNext]    = offsetNext    + ColOffset::DCOMP + ic;
          dofColIndices[localDofIndexCompCurrent] = offsetCurrent + ColOffset::DCOMP + ic;
          
          localFluxJacobian[localDofIndexCompNext]    = - dAvgDensity_dCompNext[ic]    * gravD * normalizer;  
          localFluxJacobian[localDofIndexCompCurrent] = - dAvgDensity_dCompCurrent[ic] * gravD * normalizer;
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
  });
}
  
void CompositionalMultiphaseWell::FormControlEquation( DomainPartition * const domain,
                                                       Epetra_FECrsMatrix * const jacobian,
                                                       Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom 
    arrayView1d<globalIndex const> const & wellElemDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // get the index of the well element where the control is enforced
    localIndex const iwelemControl = well->getReferenceWellElementIndex();

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

      // control equation is a normalized difference 
      // between current pressure and target pressure
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

      // control equation is a normalized difference 
      // between current rate and target rate
      real64 const controlEqn = ( currentConnRate - targetConnRate ) * normalizer;
      real64 const dControlEqn_dRate = normalizer;

      globalIndex const elemOffset  = getElementOffset( wellElemDofNumber[iwelemControl] );
      localIndex const dRateColOffset = ColOffset::DCOMP + NC; 
      globalIndex const eqnRowIndex = elemOffset + RowOffset::CONTROL;
      globalIndex const dofColIndex = elemOffset + dRateColOffset;

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

  wellManager->forSubGroups<Well>( [&] ( Well const * const well ) -> void
  {
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get a reference to the primary variables on well elements
    arrayView1d<real64> const & wellElemPressure  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64> const & wellElemGlobalCompDensity  =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    arrayView2d<real64 const> const & dWellElemGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d<real64> const & connRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

    arrayView1d<real64 const> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );   

    forall_in_range( 0, wellElementSubRegion->numWellElementsLocal(), GEOSX_LAMBDA ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        wellElemGlobalCompDensity[iwelem][ic] += dWellElemGlobalCompDensity[iwelem][ic];
      }
      connRate[iwelem] += dConnRate[iwelem];
    });

    RecordWellData( well );
  });  
}


void CompositionalMultiphaseWell::ComputeAllPerforationRates( Well const * const well )
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
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & resPhaseMob = m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dResPhaseMob_dPres = m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dResPhaseMob_dComp = m_dResPhaseMob_dCompDens;
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
  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

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

  // get secondary well data on well elements
  arrayView1d<real64 const> const & wellElemMixtureDensity =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureDensityString );
  
  arrayView1d<real64 const> const & dWellElemMixtureDensity_dPres =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::dMixtureDensity_dPressureString );

  arrayView2d<real64 const> const & dWellElemMixtureDensity_dComp =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

  arrayView2d<real64 const> const & wellElemCompFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  arrayView3d<real64 const> const & dWellElemCompFrac_dCompDens =
    wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

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
  stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

  stackArray1d<real64, 2>              pressure( 2 );
  stackArray1d<real64, 2>              dPressure_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dPressure_dC( 2, NC );
      
  stackArray1d<real64, 2>              dFlux_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dFlux_dC( 2, NC );

  stackArray1d<real64, 2>              dMult_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dMult_dC( 2, NC );

  stackArray1d<real64, maxNumComp> dMixtureDensity_dC( NC );
  stackArray1d<real64, maxNumComp> dResTotalMobility_dC( NC );

  stackArray2d<real64, 2 * maxNumComp>              phaseCompFrac( 2, NC );
  stackArray2d<real64, 2 * maxNumComp>              dPhaseCompFrac_dP( 2, NC );
  stackArray3d<real64, 2 * maxNumComp * maxNumComp> dPhaseCompFrac_dC( 2, NC, NC );
    
  stackArray1d<real64, maxNumComp> dVisc_dC( NC );
  stackArray1d<real64, maxNumComp> dRelPerm_dC( NC );

  stackArray1d<real64, 2>              dPotDiff_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dPotDiff_dC( 2, NC );

  stackArray1d<real64, 2> multiplier( 2 );


  // TODO: make this work if the wellElement and the reservoir element are on different ranks

  // loop over the perforations to compute the perforation rates
  for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
  {     

    // reset the perforation rates
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      compPerfRate[iperf][ic] = 0.0;
      for (localIndex ke = 0; ke < 2; ++ke) 
      {
        dCompPerfRate_dPres[iperf][ke][ic] = 0.0;
        for (localIndex jc = 0; jc < NC; ++jc) 
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

    pressure[SubRegionTag::RES]  = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[SubRegionTag::RES] = 1.0;

    // TODO: add a buoyancy term for the reservoir side here 

    multiplier[SubRegionTag::RES] = 1.0;

    // b) get well variables

    pressure[SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    dPressure_dP[SubRegionTag::WELL] = 1.0;

    multiplier[SubRegionTag::WELL] = -1.0;

    if (m_gravityFlag)
    {
      real64 const gravD = ( perfGravDepth[iperf] - wellElemGravDepth[iwelem] );
      pressure[SubRegionTag::WELL]  += wellElemMixtureDensity[iwelem] * gravD;
      dPressure_dP[SubRegionTag::WELL] += dWellElemMixtureDensity_dPres[iwelem] * gravD;
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dPressure_dC[SubRegionTag::WELL][ic] += dWellElemMixtureDensity_dComp[iwelem][ic] * gravD;
      }
    }

    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf]; 
        
    // 2) compute potential difference

    real64 potDiff = 0.0;
    dPotDiff_dP = 0.0;
    dPotDiff_dC = 0.0;

    for (localIndex i = 0; i < 2; ++i)
    {
      potDiff += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
      dPotDiff_dP[i] += multiplier[i] * trans * dPressure_dP[i];

      for (localIndex ic = 0; ic < NC; ++ic) 
      {
        dPotDiff_dC[i][ic] += multiplier[i] * trans * dPressure_dC[i][ic];
      }
    }

    real64 flux = 0.0;
    dFlux_dP = 0.0;
    dFlux_dC = 0.0;

    // 3) upwinding 

    if ( potDiff >= 0 ) // ** reservoir cell is upstream **
    {

      phaseCompFrac = 0.0;
      dPhaseCompFrac_dP = 0.0;
      dPhaseCompFrac_dC = 0.0;

      dPhaseCompFrac_dCompDens = 0;

      // loop over phases, compute and upwind phase flux
      // and sum contributions to each component's perforation rate
      for (localIndex ip = 0; ip < NP; ++ip)
      {
      
        // compute the phase flux and derivatives using upstream cell mobility
        flux = resPhaseMob[er][esr][ei][ip] * potDiff;
        dFlux_dP[SubRegionTag::RES]  = dResPhaseMob_dPres[er][esr][ei][ip] * potDiff 
                                     + resPhaseMob[er][esr][ei][ip]        * dPotDiff_dP[SubRegionTag::RES];
        dFlux_dP[SubRegionTag::WELL] = resPhaseMob[er][esr][ei][ip]        * dPotDiff_dP[SubRegionTag::WELL];

        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dFlux_dC[SubRegionTag::RES][ic]  = dResPhaseMob_dComp[er][esr][ei][ip][ic] * potDiff
                                           + resPhaseMob[er][esr][ei][ip]            * dPotDiff_dC[SubRegionTag::RES][ic];
          dFlux_dC[SubRegionTag::WELL][ic] = resPhaseMob[er][esr][ei][ip]            * dPotDiff_dC[SubRegionTag::RES][ic];
        }
        
        // increment component fluxes
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          compPerfRate[iperf][ic] += flux * resPhaseCompFrac[er][esr][m_resFluidIndex][ei][0][ip][ic];
        
          dCompPerfRate_dPres[iperf][SubRegionTag::RES][ic]  += resPhaseCompFrac[er][esr][m_resFluidIndex][ei][0][ip][ic]
                                                              * dFlux_dP[SubRegionTag::RES];
          dCompPerfRate_dPres[iperf][SubRegionTag::RES][ic]  += dResPhaseCompFrac_dPres[er][esr][m_resFluidIndex][ei][0][ip][ic]
                                                              * flux;
          dCompPerfRate_dPres[iperf][SubRegionTag::WELL][ic] += resPhaseCompFrac[er][esr][m_resFluidIndex][ei][0][ip][ic]
                                                              * dFlux_dP[SubRegionTag::WELL];

          applyChainRule( NC, 
                          dResCompFrac_dCompDens[er][esr][ei], 
                          dResPhaseCompFrac_dComp[er][esr][m_resFluidIndex][ei][0][ip][ic], 
                          dPhaseCompFrac_dCompDens );
           
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dCompPerfRate_dComp[iperf][SubRegionTag::RES][ic][jc]  += dFlux_dC[SubRegionTag::RES][jc] 
                                                                    * resPhaseCompFrac[er][esr][m_resFluidIndex][ei][0][ip][ic];
            dCompPerfRate_dComp[iperf][SubRegionTag::RES][ic][jc]  += flux 
                                                                    * dPhaseCompFrac_dCompDens[jc];
            dCompPerfRate_dComp[iperf][SubRegionTag::WELL][ic][jc] += dFlux_dC[SubRegionTag::WELL][jc] 
                                                                    * resPhaseCompFrac[er][esr][m_resFluidIndex][ei][0][ip][ic];
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
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        // viscosity
        real64 const resViscosity = resPhaseVisc[er][esr][m_resFluidIndex][ei][0][ip];
        real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][m_resFluidIndex][ei][0][ip];
        applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], 
                        dResPhaseVisc_dComp[er][esr][m_resFluidIndex][ei][0][ip], 
                        dVisc_dC );
      
        // relative permeability
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

        // increment total mobility
        resTotalMobility     += resRelPerm / resViscosity;
        dResTotalMobility_dP += ( dResRelPerm_dP * resViscosity - resRelPerm * dResVisc_dP )
                              / ( resViscosity * resViscosity );
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dResTotalMobility_dC[ic] += ( dRelPerm_dC[ic] * resViscosity - resRelPerm * dVisc_dC[ic] )
            / ( resViscosity * resViscosity );
        }
      }

      // compute a potdiff multiplier = wellElemMixtureDensity * resTotalMobility 
      real64 const mult = wellElemMixtureDensity[iwelem] * resTotalMobility;
      dMult_dP[SubRegionTag::RES]  = wellElemMixtureDensity[iwelem] * dResTotalMobility_dP;
      dMult_dP[SubRegionTag::WELL] = dWellElemMixtureDensity_dPres[iwelem] * resTotalMobility;
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dMult_dC[SubRegionTag::RES][ic]  = wellElemMixtureDensity[iwelem] * dResTotalMobility_dC[ic];
        dMult_dC[SubRegionTag::WELL][ic] = dWellElemMixtureDensity_dComp[iwelem][ic] * resTotalMobility;
      }

      // compute the volumetric flux and derivatives using upstream cell mobility
      flux = mult * potDiff;
      dFlux_dP[SubRegionTag::RES]  =  dMult_dP[SubRegionTag::RES] * potDiff 
                                   +  mult * dPotDiff_dP[SubRegionTag::RES];
      dFlux_dP[SubRegionTag::WELL] =  dMult_dP[SubRegionTag::WELL] * potDiff 
                                   +  mult * dPotDiff_dP[SubRegionTag::WELL];

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dFlux_dC[SubRegionTag::RES][ic]  = dMult_dC[SubRegionTag::RES][ic] * potDiff 
                                         + mult * dPotDiff_dC[SubRegionTag::RES][ic];
        dFlux_dC[SubRegionTag::WELL][ic] = dMult_dC[SubRegionTag::WELL][ic] * potDiff 
                                         + mult * dPotDiff_dC[SubRegionTag::RES][ic];
      }
      
      // compute component fluxes
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        compPerfRate[iperf][ic] += wellElemCompFrac[iwelem][ic] * flux;
        
        dCompPerfRate_dPres[iperf][SubRegionTag::RES][ic] = wellElemCompFrac[iwelem][ic]
                                                          * dFlux_dP[SubRegionTag::RES];
        dCompPerfRate_dPres[iperf][SubRegionTag::WELL][ic] = wellElemCompFrac[iwelem][ic]
                                                           * dFlux_dP[SubRegionTag::WELL];
           
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompPerfRate_dComp[iperf][SubRegionTag::RES][ic][jc]  += wellElemCompFrac[iwelem][ic]
                                                                  * dFlux_dC[SubRegionTag::RES][jc]; 
          dCompPerfRate_dComp[iperf][SubRegionTag::WELL][ic][jc] += wellElemCompFrac[iwelem][ic]
                                                                  * dFlux_dC[SubRegionTag::WELL][jc]; 
          dCompPerfRate_dComp[iperf][SubRegionTag::WELL][ic][jc] += dWellElemCompFrac_dCompDens[iwelem][ic][jc]
                                                                  * flux;
        }
      }
    }      
  }
}

void CompositionalMultiphaseWell::RecordWellData( Well const * const well )
{
  // note: this function is for debug and will go away

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;

  localIndex const NP = m_numPhases;
  localIndex const NC = m_numComponents;
  
  WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();
  PerforationData const * const perforationData = well->getPerforations();
  
  arrayView1d<real64 const> const & wellElemPressure =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & connRate =
    wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

  arrayView2d<real64 const> const & wellElemCompFrac =
    wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

  // get the index of the next well element
  arrayView1d<localIndex const> const & nextWellElemIndex =
    wellElementSubRegion->getReference<array1d<localIndex>>( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString );    
  
  arrayView2d<real64 const> const & compPerfRate =
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
      std::cout << "Mass rate at perforation #" << iperf << " for component #" << ic << ": " << compPerfRate[iperf][ic] << std::endl;
  }

  // output the reference pressure
  localIndex const iwelemControl = well->getReferenceWellElementIndex();
  real64 const pressure = wellElemPressure[iwelemControl];
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
    if (iwelem > 0 || well->getControl() == Well::Control::BHP)
      std::cout << "Mass rate at connection #" << iwelem << ": " << connRate[iwelem] << std::endl;
    else
    {
      std::cout << "Mass rate at connection #" << iwelem << ": " << connRate[iwelem]
                << ", target rate : " << well->getTargetRate() << std::endl;
    }
  
    if (iwelem > 0)
    {
      continue;
    }   
 
    stackArray1d<real64, maxNumComp> compStream( NC );

    for (localIndex ic = 0; ic < NC; ++ic)
    {
      real64 const compFrac = wellElemCompFrac[iwelem][ic];
      compStream[ic] = compFrac * connRate[iwelem];
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
  
  wellManager->forSubGroups<Well>( [&] ( Well * const well ) -> void
  {    
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the primary variables
    arrayView1d<real64 const> const & wellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellElemPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView1d<real64 const> const & connRate  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::mixtureConnRateString );

    arrayView1d<real64 const> const & dConnRate =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureConnRateString );

    // if isViable is true at the end of the following checks, no need to switch
    bool controlIsViable = false;

    // get well control and type
    Well::Control const currentControl = well->getControl();
    Well::Type const type = well->getType();
    
    // get the index of the well element where the control is enforced
    localIndex const iwelemControl = well->getReferenceWellElementIndex();

    real64 const refRate = connRate[iwelemControl] + dConnRate[iwelemControl];    
    real64 const refPressure = wellElemPressure[iwelemControl] + dWellElemPressure[iwelemControl];

    // TODO: check all inactive constraints (possibly more than one) and switch the one which is most violated
    // TODO: for the rate, use surface conditions (flash for compositional, easier for BO)

    // BHP control
    if (currentControl == Well::Control::BHP)
    {
      // the control is viable if the reference rate is below/above the max/min rate
      // targetRate specifies a max rate here
      real64 const maxRate = well->getTargetRate(); 
      controlIsViable = ( fabs(refRate) <= fabs(maxRate) );
    }
    else // rate control
    {

      // the control is viable if the reference pressure is below/above the max/min pressure
      if ( type == Well::Type::PRODUCER )
      {
        // targetBHP specifies a min pressure here
        real64 const minPressure = well->getTargetBHP(); 
        controlIsViable = ( refPressure >= minPressure );
      }
      else
      {
        // targetBHP specifies a max pressure here
        real64 const maxPressure = well->getTargetBHP(); 
        controlIsViable = ( refPressure <= maxPressure );
      }
    }

    if (!controlIsViable)
    {
      if ( currentControl == Well::Control::BHP )
      {
        well->setControl( Well::Control::LIQUIDRATE, well->getTargetRate() );
        if ( m_verboseLevel >= 1 )
        {
          GEOS_LOG_RANK_0( "Control switch for well " << well->getName()
                           << " from BHP constraint to rate constraint" );
        }
      }
      else // rate control
      {
        well->setControl( Well::Control::BHP, well->getTargetBHP() );
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
