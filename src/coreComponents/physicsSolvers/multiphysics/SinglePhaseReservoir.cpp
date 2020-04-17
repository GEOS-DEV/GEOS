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

void SinglePhaseReservoir::AddCouplingSparsityPattern( DomainPartition * const domain,
                                                       DofManager & dofManager,
                                                       ParallelMatrix & matrix )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = *domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *meshLevel.getElemManager();

  // TODO: remove this and just call SolverBase::SetupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  string const resDofKey  = dofManager.getKey( m_wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager.getKey( m_wellSolver->WellElementDofName() );

  localIndex const resNDOF = m_wellSolver->NumDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->NumDofPerWellElement();

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resDofNumber =
    elemManager.ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( resDofKey );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
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
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];
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
      matrix.insert( dofIndexRes,
                     dofIndexWell,
                     values );

      // fill J_WR
      matrix.insert( dofIndexWell,
                     dofIndexRes,
                     values );
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
  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager & elemManager = *meshLevel.getElemManager();

  string const wellDofKey = dofManager->getKey( m_wellSolver->WellElementDofName() );
  string const resDofKey  = dofManager->getKey( m_wellSolver->ResElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  resDofNumberAccessor = elemManager.ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( resDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst
    resDofNumber = resDofNumberAccessor.toViewConst();

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    PerforationData const * const perforationData = subRegion.GetPerforationData();

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
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
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
      localPerf[WellSolverBase::SubRegionTag::RES] = dt * perfRate[iperf];
      localPerf[WellSolverBase::SubRegionTag::WELL] = -localPerf[WellSolverBase::SubRegionTag::RES];

      for( localIndex ke = 0; ke < 2; ++ke )
      {
        localPerfJacobian[WellSolverBase::SubRegionTag::RES][ke] = dt * dPerfRate_dPres[iperf][ke];
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


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseReservoir, std::string const &, Group * const )

} /* namespace geosx */
