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
 * @file SinglePhaseReservoir.cpp
 *
 */


#include "SinglePhaseReservoir.hpp"

#include "common/TimingMacros.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseWell.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhaseReservoir::SinglePhaseReservoir( const std::string & name,
                                            Group * const parent ):
  ReservoirSolverBase( name, parent )
{}

SinglePhaseReservoir::~SinglePhaseReservoir()
{}

void SinglePhaseReservoir::SetupDofs( DomainPartition const * const domain,
                                      DofManager & dofManager ) const
{
  m_flowSolver->SetupDofs( domain, dofManager );
  m_wellSolver->SetupDofs( domain, dofManager );

  // TODO: add coupling when dofManager can support perforation connectors
}

void SinglePhaseReservoir::SetupSystem( DomainPartition * const domain,
                                        DofManager & dofManager,
                                        ParallelMatrix & matrix,
                                        ParallelVector & rhs,
                                        ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalDof = dofManager.numLocalDofs();

  matrix.createWithLocalSize( numLocalDof, numLocalDof, 8, MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );

  dofManager.setSparsityPattern( matrix, false ); // don't close the matrix

  // TODO: remove this and just call SolverBase::SetupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  string const resDofKey  = dofManager.getKey( m_wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager.getKey( m_wellSolver->WellElementDofName() );

  localIndex const resNDOF = m_wellSolver->NumDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->NumDofPerWellElement();

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resDofNumber =
    elemManager->ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( resDofKey );

  elemManager->forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion const & subRegion )
  {
    PerforationData const * const perforationData = subRegion.GetPerforationData();

    // get the well degrees of freedom and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get the well element indices corresponding to each perforation
    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    stackArray1d< globalIndex, 1 > dofIndexRes( resNDOF );
    stackArray1d< globalIndex, 2 > dofIndexWell( wellNDOF );
    stackArray2d< real64, 2 > values( resNDOF, wellNDOF );
    values = 1.0;

    // Insert the entries corresponding to reservoir-well perforations
    // This will fill J_WR, and J_RW
    for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er     = resElementRegion[iperf];
      localIndex const esr    = resElementSubRegion[iperf];
      localIndex const ei     = resElementIndex[iperf];
      localIndex const iwelem = perfWellElemIndex[iperf];

      for( localIndex idof = 0; idof < resNDOF; ++idof )
      {
        dofIndexRes[idof] = resDofNumber[er][esr][ei] + idof;
      }

      for( localIndex idof = 0; idof < wellNDOF; ++idof )
      {
        dofIndexWell[idof] = wellElemDofNumber[iwelem] + idof;
      }

      // fill J_RW
      matrix.insert( dofIndexRes.data(),
                     dofIndexWell.data(),
                     values.data(),
                     resNDOF,
                     wellNDOF );

      // fill J_WR
      matrix.insert( dofIndexWell.data(),
                     dofIndexRes.data(),
                     values.data(),
                     wellNDOF,
                     resNDOF );
    }

  } );

  matrix.close();
}

void SinglePhaseReservoir::AssembleCouplingTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition * const domain,
                                                  DofManager const * const dofManager,
                                                  ParallelMatrix * const matrix,
                                                  ParallelVector * const rhs )
{
  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  string const wellDofKey = dofManager->getKey( m_wellSolver->WellElementDofName() );
  string const resDofKey  = dofManager->getKey( m_wellSolver->ResElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  resDofNumberAccessor = elemManager->ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( resDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst
    resDofNumber = resDofNumberAccessor.toViewConst();

  // loop over the wells
  elemManager->forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion & subRegion )
  {
    PerforationData const * const perforationData = subRegion.GetPerforationData();

    // compute the local rates for this well
    ComputeAllPerforationRates( &subRegion );

    // get the degrees of freedom
    arrayView1d< globalIndex const > const &
    wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get well variables on perforations
    arrayView1d< real64 const > const &
    perfRate = perforationData->getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::perforationRateString );
    arrayView2d< real64 const > const &
    dPerfRate_dPres = perforationData->getReference< array2d< real64 > >( SinglePhaseWell::viewKeyStruct::dPerforationRate_dPresString );

    arrayView1d< localIndex const > const &
    perfWellElemIndex = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const &
    resElementRegion = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
    arrayView1d< localIndex const > const &
    resElementSubRegion = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
    arrayView1d< localIndex const > const &
    resElementIndex = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // local working variables and arrays
    stackArray1d< globalIndex, 2 > eqnRowIndices( 2 );
    stackArray1d< globalIndex, 2 > dofColIndices( 2 );

    stackArray1d< real64, 2 > localPerf( 2 );
    stackArray2d< real64, 4 > localPerfJacobian( 2, 2 );

    // loop over the perforations and add the rates to the residual and jacobian
    for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
    {
      eqnRowIndices = -1;
      dofColIndices = -1;

      localPerf = 0;
      localPerfJacobian = 0;

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem      = perfWellElemIndex[iperf];
      globalIndex const elemOffset = wellElemDofNumber[iwelem];

      // row index on reservoir side
      eqnRowIndices[WellSolverBase::SubRegionTag::RES] = resDofNumber[er][esr][ei];

      // column index on reservoir side
      dofColIndices[WellSolverBase::SubRegionTag::RES] = resDofNumber[er][esr][ei];

      // row index on well side
      eqnRowIndices[WellSolverBase::SubRegionTag::WELL] = elemOffset
                                                          + SinglePhaseWell::RowOffset::MASSBAL;

      // column index on well side
      dofColIndices[WellSolverBase::SubRegionTag::WELL] = elemOffset
                                                          + SinglePhaseWell::ColOffset::DPRES;

      // populate local flux vector and derivatives
      localPerf[WellSolverBase::SubRegionTag::RES]  =  dt * perfRate[iperf];
      localPerf[WellSolverBase::SubRegionTag::WELL] = -localPerf[WellSolverBase::SubRegionTag::RES];

      for( localIndex ke = 0; ke < 2; ++ke )
      {
        localPerfJacobian[WellSolverBase::SubRegionTag::RES][ke]  = dt * dPerfRate_dPres[iperf][ke];
        localPerfJacobian[WellSolverBase::SubRegionTag::WELL][ke] = -localPerfJacobian[WellSolverBase::SubRegionTag::RES][ke];
      }

      rhs->add( eqnRowIndices,
                localPerf );

      matrix->add( eqnRowIndices,
                   dofColIndices,
                   localPerfJacobian );
    }

  } );

}

void SinglePhaseReservoir::ComputeAllPerforationRates( WellElementSubRegion * const subRegion )
{

  // get the reservoir data
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  const & resPressure         = m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  const & dResPressure        = m_deltaResPressure;

  ElementRegionManager::MaterialViewAccessor< arrayView2d< real64 > > const & resDensity          = m_resDensity;
  ElementRegionManager::MaterialViewAccessor< arrayView2d< real64 > > const & dResDensity_dPres   = m_dResDens_dPres;
  ElementRegionManager::MaterialViewAccessor< arrayView2d< real64 > > const & resViscosity        = m_resViscosity;
  ElementRegionManager::MaterialViewAccessor< arrayView2d< real64 > > const & dResViscosity_dPres = m_dResVisc_dPres;

  // get the well data
  PerforationData * const perforationData = subRegion->GetPerforationData();

  // get the degrees of freedom and depth
  arrayView1d< real64 const > const &
  wellElemGravCoef = subRegion->getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::gravityCoefString );

  // get well primary variables on well elements
  arrayView1d< real64 const > const &
  wellElemPressure = subRegion->getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::pressureString );
  arrayView1d< real64 const > const &
  dWellElemPressure = subRegion->getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::deltaPressureString );

  // get well constitutive data
  string const & fluidName = m_wellSolver->getReference< string >( WellSolverBase::viewKeyStruct::fluidNameString );
  SingleFluidBase const * const fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

  arrayView2d< real64 const > const &
  wellElemDensity = fluid->getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString );
  arrayView2d< real64 const > const &
  dWellElemDensity_dPres = fluid->getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::dDens_dPresString );

  arrayView2d< real64 const > const &
  wellElemViscosity = fluid->getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::viscosityString );
  arrayView2d< real64 const > const &
  dWellElemViscosity_dPres = fluid->getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::dVisc_dPresString );

  // get well variables on perforations
  arrayView1d< real64 const > const &
  perfGravCoef = perforationData->getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::gravityCoefString );

  arrayView1d< localIndex const > const &
  perfWellElemIndex = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d< real64 const > const &
  perfTransmissibility = perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::transmissibilityString );

  arrayView1d< real64 > const &
  perfRate = perforationData->getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::perforationRateString );
  arrayView2d< real64 > const &
  dPerfRate_dPres = perforationData->getReference< array2d< real64 > >( SinglePhaseWell::viewKeyStruct::dPerforationRate_dPresString );

  // get the element region, subregion, index
  arrayView1d< localIndex const > const &
  resElementRegion = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString );
  arrayView1d< localIndex const > const &
  resElementSubRegion = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString );
  arrayView1d< localIndex const > const &
  resElementIndex = perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString );

  // local working variables and arrays
  stackArray1d< globalIndex, 2 > eqnRowIndices( 2 );
  stackArray1d< globalIndex, 2 > dofColIndices( 2 );

  stackArray1d< real64, 2 > pressure( 2 );
  stackArray1d< real64, 2 > dPressure_dP( 2 );

  stackArray1d< localIndex, 2 > multiplier( 2 );

  // loop over the perforations to compute the perforation rates
  for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
  {
    eqnRowIndices = -1;
    dofColIndices = -1;

    pressure = 0;
    dPressure_dP = 0;

    multiplier = 0;

    // 1) Reservoir side

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    // get reservoir variables
    pressure[WellSolverBase::SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[WellSolverBase::SubRegionTag::RES] = 1;

    // TODO: add a buoyancy term for the reservoir side here

    // multiplier for reservoir side in the flux
    multiplier[WellSolverBase::SubRegionTag::RES] = 1;

    // 2) Well side

    // get the local index of the well element
    localIndex const iwelem = perfWellElemIndex[iperf];

    // get well variables
    pressure[WellSolverBase::SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    dPressure_dP[WellSolverBase::SubRegionTag::WELL] = 1.0;

    real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );
    pressure[WellSolverBase::SubRegionTag::WELL]     += wellElemDensity[iwelem][0] * gravD;
    dPressure_dP[WellSolverBase::SubRegionTag::WELL] += dWellElemDensity_dPres[iwelem][0] * gravD;

    // multiplier for well side in the flux
    multiplier[WellSolverBase::SubRegionTag::WELL] = -1;

    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf];

    // compute potential difference
    real64 potDif = 0.0;
    for( localIndex i = 0; i < 2; ++i )
    {
      potDif += multiplier[i] * trans * pressure[i];
      dPerfRate_dPres[iperf][i] = multiplier[i] * trans * dPressure_dP[i];
    }

    // choose upstream cell based on potential difference
    localIndex const k_up = (potDif >= 0)
                          ? WellSolverBase::SubRegionTag::RES
                          : WellSolverBase::SubRegionTag::WELL;

    // compute upstream density, viscosity, and mobility
    real64 densityUp       = 0.0;
    real64 dDensityUp_dP   = 0.0;
    real64 viscosityUp     = 0.0;
    real64 dViscosityUp_dP = 0.0;

    // upwinding the variables
    if( k_up == WellSolverBase::SubRegionTag::RES ) // use reservoir vars
    {
      localIndex const & resFluidIndex = m_wellSolver->getReference< localIndex >( WellSolverBase::viewKeyStruct::resFluidIndexString );

      densityUp     = resDensity[er][esr][resFluidIndex][ei][0];
      dDensityUp_dP = dResDensity_dPres[er][esr][resFluidIndex][ei][0];

      viscosityUp     = resViscosity[er][esr][resFluidIndex][ei][0];
      dViscosityUp_dP = dResViscosity_dPres[er][esr][resFluidIndex][ei][0];
    }
    else // use well vars
    {
      densityUp = wellElemDensity[iwelem][0];
      dDensityUp_dP = dWellElemDensity_dPres[iwelem][0];

      viscosityUp = wellElemViscosity[iwelem][0];
      dViscosityUp_dP = dWellElemViscosity_dPres[iwelem][0];
    }

    // compute mobility
    real64 const mobilityUp     = densityUp / viscosityUp;
    real64 const dMobilityUp_dP = dDensityUp_dP / viscosityUp
                                  - mobilityUp / viscosityUp * dViscosityUp_dP;

    perfRate[iperf] = mobilityUp * potDif;
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      dPerfRate_dPres[iperf][ke] *= mobilityUp;
    }
    dPerfRate_dPres[iperf][k_up] += dMobilityUp_dP * potDif;
  }
}


void SinglePhaseReservoir::ResetViews( DomainPartition * const domain )
{
  ReservoirSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_resPressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString );

  m_deltaResPressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString );

  m_resDensity =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString,
                                                                                                constitutiveManager );
  m_dResDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::dDens_dPresString,
                                                                                                constitutiveManager );
  m_resViscosity =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::viscosityString,
                                                                                                constitutiveManager );
  m_dResVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::dVisc_dPresString,
                                                                                                constitutiveManager );
}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseReservoir, std::string const &, Group * const )

} /* namespace geosx */
