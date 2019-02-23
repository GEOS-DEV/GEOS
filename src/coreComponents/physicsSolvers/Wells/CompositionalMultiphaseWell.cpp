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

  WellManager * wellManager = meshBodies->getParent()->group_cast<DomainPartition *>()->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {

    WellElementSubRegion * wellElementSubRegion = well->getWellElements();
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

    ConnectionData * connectionData = well->getConnections();
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::mixtureVelocityString );
    connectionData->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

    PerforationData * perforationData = well->getPerforations();
    perforationData->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::phaseFlowRateString );

  });    
}
  
void CompositionalMultiphaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();
  
  MultiFluidBase const * fluid = cm->GetConstitituveRelation<MultiFluidBase>( m_fluidName );
  
  m_numPhases        = fluid->numFluidPhases();
  m_numComponents    = fluid->numFluidComponents();

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

void CompositionalMultiphaseWell::UpdateComponentFractionAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    arrayView2d<real64 const> const & wellCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    arrayView2d<real64 const> const & dWellCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView2d<real64> const & wellCompFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

    arrayView3d<real64> const & dWellCompFrac_dCompDens =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      real64 wellTotalDensity = 0.0;

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        wellTotalDensity += wellCompDens[iwelem][ic]
                          + dWellCompDens[iwelem][ic];
      }

      real64 const wellTotalDensityInv = 1.0 / wellTotalDensity;

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        wellCompFrac[iwelem][ic] = (wellCompDens[iwelem][ic]
				 + dWellCompDens[iwelem][ic]) * wellTotalDensityInv;

        for (localIndex jc = 0; jc < m_numComponents; ++jc)
        {
          dWellCompFrac_dCompDens[iwelem][ic][jc] = - wellCompFrac[iwelem][ic] * wellTotalDensityInv;
        }
        dWellCompFrac_dCompDens[iwelem][ic][ic] += wellTotalDensityInv;
      }
    }
  });
}

void CompositionalMultiphaseWell::UpdatePhaseVolumeFractionAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();  
    MultiFluidBase * fluid = GetFluidModel( wellElementSubRegion );

    arrayView2d<real64> const & wellPhaseVolFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

    arrayView2d<real64> const & dWellPhaseVolFrac_dPres =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    arrayView3d<real64> const & dWellPhaseVolFrac_dComp =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    arrayView3d<real64> const & dWellCompFrac_dCompDens =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

    arrayView2d<real64 const> const & wellCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    arrayView2d<real64 const> const & dWellCompDens =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    arrayView3d<real64 const> const & wellPhaseFrac =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseFractionString );

    arrayView3d<real64 const> const & dWellPhaseFrac_dPres =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

    arrayView4d<real64 const> const & dWellPhaseFrac_dComp =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

    arrayView3d<real64 const> const & wellPhaseDens =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

    arrayView3d<real64 const> const & dWellPhaseDens_dPres =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

    arrayView4d<real64 const> const & dWellPhaseDens_dComp =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dGlobalCompFractionString );

    localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
    localIndex const NC = m_numComponents;
    localIndex const NP = m_numPhases;

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      stackArray1d<real64, maxNumComp> work( NC );

      // compute total density from component partial densities
      real64 wellTotalDensity = 0.0;
      real64 const dWellTotalDens_dCompDens = 1.0;
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        wellTotalDensity += wellCompDens[iwelem][ic] + dWellCompDens[iwelem][ic];
      }

      for (localIndex ip = 0; ip < NP; ++ip)
      {
        // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
        real64 const wellPhaseDensInv = 1.0 / wellPhaseDens[iwelem][0][ip];

        // compute saturation and derivatives except multiplying by the total density
        wellPhaseVolFrac[iwelem][ip] = wellPhaseFrac[iwelem][0][ip] * wellPhaseDensInv;

        dWellPhaseVolFrac_dPres[iwelem][ip] =
          (dWellPhaseFrac_dPres[iwelem][0][ip] - wellPhaseVolFrac[iwelem][ip] * dWellPhaseDens_dPres[iwelem][0][ip])
	  * wellPhaseDensInv;

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dWellPhaseVolFrac_dComp[iwelem][ip][jc] =
            (dWellPhaseFrac_dComp[iwelem][0][ip][jc] - wellPhaseVolFrac[iwelem][ip] * dWellPhaseDens_dComp[iwelem][0][ip][jc])
	    * wellPhaseDensInv;
        }

        // apply chain rule to convert derivatives from global component fractions to densities
        applyChainRuleInPlace( NC, dWellCompFrac_dCompDens[iwelem], dWellPhaseVolFrac_dComp[iwelem][ip], work );

        // now finalize the computation by multiplying by total density
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dWellPhaseVolFrac_dComp[iwelem][ip][jc] *= wellTotalDensity;
          dWellPhaseVolFrac_dComp[iwelem][ip][jc] += wellPhaseVolFrac[iwelem][ip] * dWellTotalDens_dCompDens;
        }

        wellPhaseVolFrac[iwelem][ip] *= wellTotalDensity;
        dWellPhaseVolFrac_dPres[iwelem][ip] *= wellTotalDensity;
      }
    }
  });
}

void CompositionalMultiphaseWell::UpdateFluidModelAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
    MultiFluidBase * const fluid = GetFluidModel( wellElementSubRegion );

    arrayView1d<real64 const> const & wellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    arrayView2d<real64 const> const & wellCompFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompFractionString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      fluid->PointUpdate( wellPressure[iwelem] + dWellPressure[iwelem], m_temperature, wellCompFrac[iwelem], iwelem, 0 );
    }
  });
}
  
void CompositionalMultiphaseWell::UpdateRelPermModelAll( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * const wellElementSubRegion = well->getWellElements();
    RelativePermeabilityBase * const relPerm = GetRelPermModel( wellElementSubRegion );

    arrayView2d<real64> const & wellPhaseVolFrac =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::phaseVolumeFractionString );

    relPerm->BatchUpdate( wellPhaseVolFrac );
  });
}

void CompositionalMultiphaseWell::UpdateStateAll( DomainPartition * const domain )
{
  UpdateComponentFractionAll( domain );
  UpdateFluidModelAll( domain );
  UpdatePhaseVolumeFractionAll( domain );
  UpdateRelPermModelAll( domain );
}

void CompositionalMultiphaseWell::InitializeWellState( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    // do something
  });
}

void CompositionalMultiphaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );  
}

void CompositionalMultiphaseWell::BackupFields( DomainPartition * const domain )
{
  // TODO
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

  // TODO: use that
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> const & resElemGhostRank =
    elementRegionManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS; 
  localIndex constexpr maxNumDof  = maxNumComp + 2; // dofs are 1 pressure, NC comp densities (reservoir and well), 1 velocity (well only)
  localIndex const resNDOF  = numDofPerResElement; // reservoir dofs are 1 pressure and NC comp densities
  localIndex const wellNDOF = numDofPerElement()
                            + numDofPerConnection(); // well dofs are 1 pressure, NC comp densities, 1 velocity

  //**** loop over all perforations. Fill in sparsity for all pairs of DOF/elem that are connected by a perforation
  wellManager->forSubGroups<Well>( [&] ( Well const * well ) -> void
  {
    PerforationData const * const perforationData = well->getPerforations();
    ConnectionData const * const connectionData   = well->getConnections();
    WellElementSubRegion const * const wellElementSubRegion = well->getWellElements();

    // get the well degrees of freedom 
    array1d<globalIndex> const & wellDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
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
        elementLocalDofIndexRow[ElemTag::RES * wellNDOF + idof] = resOffset + idof;
	// specify the reservoir variable number
	elementLocalDofIndexCol[ElemTag::RES * wellNDOF + idof] = resOffset + idof;
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

      localIndex const iwelem = 0; // this is a hack (we assume single-segmented wells)
      
      // get the offset of the well element equation
      globalIndex const currentElemDofNumber = wellDofNumber[iwelem];
      globalIndex const currentElemOffset    = firstWellElemDofNumber * resNDOF // number of eqns in J_RR
                                             + (currentElemDofNumber - firstWellElemDofNumber) * wellNDOF; // number of eqns in J_WW, before this element's equations
      for (localIndex idof = 0; idof < wellNDOF; ++idof)
      {
	elementLocalDofIndexRow[ElemTag::WELL * resNDOF + idof] = currentElemOffset + idof;
	elementLocalDofIndexCol[ElemTag::WELL * resNDOF + idof] = currentElemOffset + idof;
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
      localIndex iwelemPrev = -1; // TODO

      // get second segment index
      localIndex iwelemNext = -1; // TODO

      // check if this is not an entry or exit
      if (iwelemPrev < 0 || iwelemNext < 0)
	continue;
      
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexRow( 2 * wellNDOF );
      stackArray1d<globalIndex, 2 * maxNumDof > elementLocalDofIndexCol( 2 * wellNDOF );
      
      // get the offset of the well element equations
      globalIndex const prevElemDofNumber = wellDofNumber[iwelemPrev];
      globalIndex const prevElemOffset    = firstWellElemDofNumber * resNDOF // number of eqns in J_RR
                                          + (prevElemDofNumber - firstWellElemDofNumber) * wellNDOF; // number of eqns in J_WW, before this element's equations
      globalIndex const nextElemDofNumber = wellDofNumber[iwelemNext];
      globalIndex const nextElemOffset    = firstWellElemDofNumber * resNDOF // number of eqns in J_RR
                                          + (nextElemDofNumber - firstWellElemDofNumber) * wellNDOF; // number of eqns in J_WW, before this element's equations
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

  CheckWellControlSwitch( domain );
  
  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleVolumeBalanceTerms( domain, jacobian, residual, time_n, dt );
  AssembleSourceTerms( domain, jacobian, residual, time_n, dt );

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
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    // loop over the segments
  });  
}

void CompositionalMultiphaseWell::AssembleFluxTerms( DomainPartition * const domain,
                                                     Epetra_FECrsMatrix * const jacobian,
                                                     Epetra_FEVector * const residual,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

   // loop over the wells
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    ConnectionData const * const connectionData = well->getConnections();

    // for each well, loop over the connections
    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      Connection const * const connection = connectionData->getConnection( iconn );

      std::cout << "CompositionalMultiphaseWell: computing flux terms for connection "
		<< connection->getName()
		<< std::endl;
      // compute flux term and add to mass conservation in the two neighboring segments
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

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;
  
  localIndex const NC      = m_numComponents;
  localIndex const NP      = m_numPhases;
  localIndex const resNDOF = NC + 2; 
 
  real64 constexpr densWeight[2] = { 1.0, 0.0 };

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    PerforationData * perforationData = well->getPerforations();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get the degrees of freedom
    array1d<globalIndex> const & wellDofNumber =
      wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
    // get well primary variables on segments
    array1d<real64> const & wellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array2d<real64> const & dWellPhaseVolFrac_dPres =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dPressureString );

    array3d<real64> const & dWellPhaseVolFrac_dComp =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

    array3d<real64> const & dWellCompFrac_dCompDens =
      wellElementSubRegion->getReference<array3d<real64>>( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
    
    // get well secondary variables on segments
    MultiFluidBase * const fluid = GetFluidModel( wellElementSubRegion );
    RelativePermeabilityBase * const relPerm = GetRelPermModel( wellElementSubRegion );
    
    array3d<real64> const & wellPhaseDens =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );

    array3d<real64> const & dWellPhaseDens_dPres =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dPressureString );

    array4d<real64> const & dWellPhaseDens_dComp =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );

    array3d<real64> const & wellPhaseVisc =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString );

    array3d<real64> const & dWellPhaseVisc_dPres =
      fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

    array4d<real64> const & dWellPhaseVisc_dComp =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseFraction_dGlobalCompFractionString );
    
    array4d<real64> const & wellPhaseCompFrac =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString );

    array4d<real64> const & dWellPhaseCompFrac_dPres =
      fluid->getReference<array4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString );

    array5d<real64> const & dWellPhaseCompFrac_dComp =
      fluid->getReference<array5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString );

    array3d<real64> const & wellPhaseRelPerm =
      relPerm->getReference<array3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString );

    array4d<real64> const & dWellPhaseRelPerm_dPhaseVolFrac =
      relPerm->getReference<array4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString );

    // get well variables on perforations
    array2d<real64> const & wellPhaseFlowRate  =
      perforationData->getReference<array2d<real64>>( viewKeyStruct::phaseFlowRateString );

    array1d<real64> const & wellGravDepth =
      perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::gravityDepthString );

    // get the element region, subregion, index
    arrayView1d<localIndex const> const & resElementRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & resElementSubRegion =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & resElementIndex =
      perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );
    
    for (localIndex iperf = 0; iperf < perforationData->numPerforationsLocal(); ++iperf)
    {
      
      // local working variables and arrays
      stackArray1d<long long, 2 * maxNumComp> eqnRowIndices( 2 * NC );
      stackArray1d<long long, 2 * maxNumDof>  dofColIndices( 2 * resNDOF );

      stackArray1d<double, 2 * maxNumComp>                 localFlux( 2 * NC );
      stackArray2d<double, 2 * maxNumComp * 2 * maxNumDof> localFluxJacobian( 2 * NC, 2 * resNDOF );
      
      stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

      stackArray1d<real64, 2> pressure( 2 );
      stackArray1d<real64, 2> gravDepth( 2 );
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

      stackArray1d<real64, 2>              dGravHead_dP( 2 );
      stackArray2d<real64, 2 * maxNumComp> dGravHead_dC( 2, NC );

      stackArray1d<real64, 2> multiplier( 2 );
      
      // reset the local values
      compFlux = 0.0;
      dCompFlux_dP = 0.0;
      dCompFlux_dC = 0.0;

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      globalIndex const resOffset  = resNDOF * resDofNumber[er][esr][ei];
      // TODO: wellOffset = firstWellOffset + wellNDOF * iwell * wellElemLocalDofNumber[iwelem]
      globalIndex const wellOffset = 0; // temp
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        eqnRowIndices[ElemTag::RES  * NC + ic] = resOffset + ic;
        eqnRowIndices[ElemTag::WELL * NC + ic] = wellOffset + ic + 1;
      }
      for (localIndex jdof = 0; jdof < resNDOF; ++jdof)
      {
        // below, we insert the well mass balance equations at wellDofNumber + 1,
        // the +1 is needed because the control equations are at wellDofNumber
        dofColIndices[ElemTag::RES  * resNDOF + jdof] = resOffset  + jdof + 1;
        dofColIndices[ElemTag::WELL * resNDOF + jdof] = wellOffset + jdof + 1;
      }
      
      // loop over phases, compute and upwind phase flux
      // and sum contributions to each component's perforation rate
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        // clear working arrays
        real64 densMean = 0.0;
        dDensMean_dP = 0.0;
        dDensMean_dC = 0.0;

        real64 presGrad = 0.0;
        dPresGrad_dP = 0.0;

        real64 gravHead = 0.0;
        dGravHead_dP = 0.0;
        dGravHead_dC = 0.0;

        real64 phaseFlux = 0.0;
        dPhaseFlux_dP = 0.0;
        dPhaseFlux_dC = 0.0;

	// 1) get reservoir variables

	// pressure and gravDepth
	pressure[ElemTag::RES]  = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
	gravDepth[ElemTag::RES] = resGravDepth[er][esr][ei];
	
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
        mobility[ElemTag::RES] = resRelPerm * resDensity / resViscosity;
        dMobility_dP[ElemTag::RES] = dResRelPerm_dP * resDensity / resViscosity
                                   + mobility[ElemTag::RES] * (dResDens_dP / resDensity - dResVisc_dP / resViscosity);

        // average density and pressure derivative
        densMean += densWeight[ElemTag::RES] * resDensity;
        dDensMean_dP[ElemTag::RES] = densWeight[ElemTag::RES] * dResDens_dP;

        // compositional derivatives
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dDensMean_dC[ElemTag::RES][jc] = densWeight[ElemTag::RES] * dDens_dC[jc];

          dMobility_dC[ElemTag::RES][jc] = dRelPerm_dC[jc] * resDensity / resViscosity
	                                 + mobility[ElemTag::RES] * (dDens_dC[jc] / resDensity - dVisc_dC[jc] / resViscosity);
        }

	multiplier[ElemTag::RES] = 1;

        // 2) get well variables

	localIndex const iwelem = 0; // this is a hack
	
	// pressure
	pressure[ElemTag::WELL]  = wellPressure[iwelem] + dWellPressure[iwelem];
        gravDepth[ElemTag::WELL] = wellGravDepth[iwelem];
	
        // density
        real64 const wellDensity  = wellPhaseDens[iwelem][0][ip];
        real64 const dWellDens_dP = dWellPhaseDens_dPres[iwelem][0][ip];
        applyChainRule( NC, dWellCompFrac_dCompDens[iwelem], dWellPhaseDens_dComp[iwelem][0][ip], dDens_dC );

        // viscosity
        real64 const wellViscosity = wellPhaseVisc[iwelem][0][ip];
        real64 const dWellVisc_dP  = dWellPhaseVisc_dPres[iwelem][0][ip];
        applyChainRule( NC, dWellCompFrac_dCompDens[iwelem], dWellPhaseVisc_dComp[iwelem][0][ip], dVisc_dC );

        //relative permeability
        real64 const wellRelPerm = wellPhaseRelPerm[iwelem][0][ip];
        real64 dWellRelPerm_dP = 0.0;
        dRelPerm_dC = 0.0;
        for (localIndex jp = 0; jp < NP; ++jp)
        {
          real64 const dWellRelPerm_dS = dWellPhaseRelPerm_dPhaseVolFrac[iwelem][0][ip][jp];
          dWellRelPerm_dP += dWellRelPerm_dS * dWellPhaseVolFrac_dPres[iwelem][jp];

          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dRelPerm_dC[jc] += dWellRelPerm_dS * dWellPhaseVolFrac_dComp[iwelem][jp][jc];
          }
        }

        // mobility and pressure derivative
        mobility[ElemTag::WELL] = wellRelPerm * wellDensity / wellViscosity;
        dMobility_dP[ElemTag::WELL] = dWellRelPerm_dP * wellDensity / wellViscosity
                                    + mobility[ElemTag::WELL] * (dWellDens_dP / wellDensity - dWellVisc_dP / wellViscosity);

        // average density and pressure derivative
        densMean += densWeight[ElemTag::WELL] * wellDensity;
        dDensMean_dP[ElemTag::WELL] = densWeight[ElemTag::WELL] * dResDens_dP;

        // compositional derivatives
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dDensMean_dC[ElemTag::WELL][jc] = densWeight[ElemTag::WELL] * dDens_dC[jc];

          dMobility_dC[ElemTag::WELL][jc] = dRelPerm_dC[jc] * wellDensity / wellViscosity
	                                  + mobility[ElemTag::WELL] * (dDens_dC[jc] / wellDensity - dVisc_dC[jc] / wellViscosity);
        }

	multiplier[ElemTag::WELL] = -1;
	
        //***** calculation of flux *****
	
        // TODO: use distinct treatments for injector and producer
	// TODO: multiply mobility of injector by total reservoir mobility
        // TODO: account for depth difference between well element center and perforation
	
        // get transmissibility at the interface
        Perforation * perforation = perforationData->getPerforation( iperf ); 
        const real64 trans = perforation->getTransmissibility(); // TODO: change to an array of transmissibilities in PerforationData

	// compute potential difference
        for (localIndex i = 0; i < 2; ++i)
        {
          presGrad += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
          dPresGrad_dP[i] += multiplier[i] * trans;

          if (m_gravityFlag)
          {
            real64 const gravD = multiplier[i] * trans * gravDepth[i];
            gravHead += densMean * gravD;

            // need to add contributions from both cells the mean density depends on
	    for (localIndex j = 0; j < 2; ++j)
	    {
              dGravHead_dP[j] += dDensMean_dP[j] * gravD;
              for (localIndex jc = 0; jc < NC; ++jc)
              {
                dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
              }
	    }
	  }
        }

        // *** upwinding ***

        // compute phase potential gradient
        real64 potGrad = presGrad + gravHead;

        // choose upstream cell
        localIndex const k_up = (potGrad >= 0) ? ElemTag::RES : ElemTag::WELL;

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

        // gravitational head depends only on the two cells connected (same as mean density)
        for (localIndex ke = 0; ke < 2; ++ke)
        {
          dPhaseFlux_dP[ke] += dGravHead_dP[ke];
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dPhaseFlux_dC[ke][jc] += dGravHead_dC[ke][jc];
          }
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
 
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        localFlux[ElemTag::RES  * NC + ic] =  dt * compFlux[ic];
        localFlux[ElemTag::WELL * NC + ic] = -dt * compFlux[ic];
 
        for (localIndex ke = 0; ke < 2; ++ke)
        {
          localIndex const localDofIndexPres = ke * resNDOF;
          localFluxJacobian[ElemTag::RES  * NC + ic][localDofIndexPres] = dt * dCompFlux_dP[ke][ic];
          localFluxJacobian[ElemTag::WELL * NC + ic][localDofIndexPres] = -dt * dCompFlux_dP[ke][ic];

          for (localIndex jc = 0; jc < NC; ++jc)
          {
            localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
            localFluxJacobian[ElemTag::RES  * NC + ic][localDofIndexComp] = dt * dCompFlux_dC[ke][ic][jc];
            localFluxJacobian[ElemTag::WELL * NC + ic][localDofIndexComp] = -dt * dCompFlux_dC[ke][ic][jc];
          }
        }
      } 
 
      // Add to global residual/jacobian
      residual->SumIntoGlobalValues( integer_conversion<int>( 2 * NC ),
                                     eqnRowIndices.data(),
                                     localFlux.data() );

      jacobian->SumIntoGlobalValues( integer_conversion<int>( 2 * NC ),
                                     eqnRowIndices.data(),
                                     integer_conversion<int>( 2 * resNDOF ),
                                     dofColIndices.data(),
                                     localFluxJacobian.data(),
                                     Epetra_FECrsMatrix::ROW_MAJOR );

    }

    FormControlEquation( domain, jacobian, residual );
    
  });
}

void CompositionalMultiphaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    // check if the well control needs to be switched
  });
}


real64
CompositionalMultiphaseWell::CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                                                    DomainPartition * const domain )
{
  return 0.0;
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

    // get a reference to the primary variables on segments
    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array2d<real64> const & dWellGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    // get a reference to the primary variables on connections
    array1d<real64> const & dWellVelocity =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      // TODO: check for ghost segments

      // extract solution and apply to dP
      globalIndex const dummyOffset = 0;
      int lid = rowMap->LID( integer_conversion<int>( dummyOffset ) );
      dWellPressure[iwelem] += scalingFactor * local_solution[lid];

      for (localIndex ic = 0; ic < m_numComponents; ++ic)
      {
        lid = rowMap->LID( integer_conversion<int>( dummyOffset + ic + 1 ) );
        dWellGlobalCompDensity[iwelem][ic] += scalingFactor * local_solution[lid];
      }

    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {

      // TODO: check for ghost connections if needed
      // TODO: check if there is a primary var defined on this connection

      // extract solution and apply to dP
      globalIndex const dummyDofNumber = 0;
      int const lid = rowMap->LID( integer_conversion<int>( dummyDofNumber ) );
      dWellVelocity[iconn] += scalingFactor * local_solution[lid];
      
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

    // get a reference to the primary variables on segments
    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array2d<real64> const & dWellGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    // get a reference to the primary variables on connections
    array1d<real64> const & dWellVelocity =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {

      // TODO: check for ghost segments

      // extract solution and apply to dP
      dWellPressure[iwelem] = 0;
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
        dWellGlobalCompDensity[iwelem][ic] = 0;      
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {

      // TODO: check for ghost connections if needed
      // TODO: check if there is a primary var defined on this connection

      // extract solution and apply to dP
      dWellVelocity[iconn] = 0;

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


void CompositionalMultiphaseWell::FormControlEquation( DomainPartition * const domain,
                                                       Epetra_FECrsMatrix * const jacobian,
                                                       Epetra_FEVector * const residual )
{
  WellManager * const wellManager = domain->getWellManager();
  
  wellManager->forSubGroups<Well>( [&] ( Well * well ) -> void
  {
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get local index of the reference element where the control eqn will be formed
    localIndex const iwelemRef = wellElementSubRegion->getReferenceElem();
    if (iwelemRef >= 0)
    {
    
      // get well control
      Well::Control control = well->getWellControl();

      // get the degrees of freedom 
      array1d<globalIndex> const & wellDofNumber =
        wellElementSubRegion->getReference<array1d<globalIndex>>( viewKeyStruct::dofNumberString );
    
      // get pressure data
      array1d<real64> const & wellPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );
      array1d<real64> const & dWellPressure =
        wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

      real64 controlEqn        = 0;
      real64 dControlEqn_dPres = 0;
    
      // BHP control
      if (control == Well::Control::BHP)
      {

        real64 const currentBHP = wellPressure[iwelemRef] + dWellPressure[iwelemRef];
        real64 const targetBHP  = well->getTargetBHP();

        controlEqn = currentBHP - targetBHP;
        dControlEqn_dPres = 1;

      }
      // rate control
      else 
      {
        // TODO: implement this
        controlEqn = 0;
        dControlEqn_dPres = 0;
      }

      // the control equation is the first equation of the element block 
      globalIndex welemRefDOF = wellDofNumber[iwelemRef];
    
      // add contribution to global residual and jacobian
      residual->SumIntoGlobalValues( 1, &welemRefDOF, &controlEqn );
      jacobian->SumIntoGlobalValues( 1, &welemRefDOF, 1, &welemRefDOF, &dControlEqn_dPres );
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
    ConnectionData * connectionData = well->getConnections();
    WellElementSubRegion * wellElementSubRegion = well->getWellElements();

    // get a reference to the primary variables on segments
    array1d<real64> const & wellPressure  =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::pressureString );

    array1d<real64> const & dWellPressure =
      wellElementSubRegion->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

    array2d<real64> const & wellGlobalCompDensity  =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::globalCompDensityString );

    array2d<real64> const & dWellGlobalCompDensity =
      wellElementSubRegion->getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );

    // get a reference to the primary variables on connections
    array1d<real64> const & wellVelocity  =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::mixtureVelocityString );

    array1d<real64> const & dWellVelocity =
      connectionData->getReference<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

    for (localIndex iwelem = 0; iwelem < wellElementSubRegion->numWellElementsLocal(); ++iwelem)
    {
      wellPressure[iwelem] += dWellPressure[iwelem];
      for (localIndex ic = 0; ic < m_numComponents; ++ic)
	wellGlobalCompDensity[iwelem][ic] += dWellGlobalCompDensity[iwelem][ic];
    }

    for (localIndex iconn = 0; iconn < connectionData->numConnectionsLocal(); ++iconn)
    {
      // TODO: check if there is a variable on this connection
      wellVelocity[iconn] += dWellVelocity[iconn];
    }    
  }); 
 
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseWell, string const &, ManagedGroup * const)
}// namespace geosx
