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
 * @file CompositionalMultiphaseReservoir.cpp
 *
 */


#include "CompositionalMultiphaseReservoir.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlow.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseWell.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

CompositionalMultiphaseReservoir::CompositionalMultiphaseReservoir( const std::string & name,
                                                                    Group * const parent ):
  ReservoirSolverBase( name, parent )
{ m_linearSolverParameters.get().mgr.strategy = "CompositionalMultiphaseReservoir"; }

CompositionalMultiphaseReservoir::~CompositionalMultiphaseReservoir()
{}

void CompositionalMultiphaseReservoir::AddCouplingSparsityPattern( DomainPartition * const domain,
                                                                   DofManager & dofManager,
                                                                   ParallelMatrix & matrix )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  // TODO: remove this and just call SolverBase::SetupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  string const resDofKey  = dofManager.getKey( m_wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager.getKey( m_wellSolver->WellElementDofName() );

  localIndex const resNDOF = m_wellSolver->NumDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->NumDofPerWellElement();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

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

    stackArray1d< globalIndex, maxNumDof > dofIndexRes( resNDOF );
    stackArray1d< globalIndex, maxNumDof > dofIndexWell( wellNDOF );
    stackArray2d< real64, maxNumDof * maxNumDof > valuesResWell( resNDOF, wellNDOF );
    stackArray2d< real64, maxNumDof * maxNumDof > valuesWellRes( wellNDOF, resNDOF );
    valuesResWell = 1.0;
    valuesWellRes = 1.0;

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
                     valuesResWell );

      // fill J_WR
      matrix.insert( dofIndexWell,
                     dofIndexRes,
                     valuesWellRes );
    }

  } );

  matrix.close();
}

void CompositionalMultiphaseReservoir::AssembleCouplingTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const dt,
                                                              DomainPartition * const domain,
                                                              DofManager const * const dofManager,
                                                              ParallelMatrix * const matrix,
                                                              ParallelVector * const rhs )
{
  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  localIndex const NC      = m_wellSolver->NumFluidComponents();
  localIndex const resNDOF = m_wellSolver->NumDofPerResElement();

  string const resDofKey  = dofManager->getKey( m_wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager->getKey( m_wellSolver->WellElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  resDofNumberAccessor = elemManager->ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( resDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >::ViewTypeConst
    resDofNumber = resDofNumberAccessor.toViewConst();

  elemManager->forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion & subRegion )
  {

    PerforationData const * const perforationData = subRegion.GetPerforationData();

    // get the degrees of freedom
    arrayView1d< globalIndex const > const &
    wellElemDofNumber = subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get well variables on perforations
    arrayView2d< real64 const > const &
    compPerfRate = perforationData->getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::compPerforationRateString );
    arrayView3d< real64 const > const &
    dCompPerfRate_dPres = perforationData->getReference< array3d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dPresString );
    arrayView4d< real64 const > const &
    dCompPerfRate_dComp = perforationData->getReference< array4d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dCompString );

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
    stackArray1d< long long, 2 * maxNumComp > eqnRowIndices( 2 * NC );
    stackArray1d< long long, 2 * maxNumDof >  dofColIndices( 2 * resNDOF );

    stackArray1d< double, 2 * maxNumComp >                 localPerf( 2 * NC );
    stackArray2d< double, 2 * maxNumComp * 2 * maxNumDof > localPerfJacobian( 2 * NC, 2 * resNDOF );

    // loop over the perforations and add the rates to the residual and jacobian
    for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
    {

      // local working variables and arrays
      eqnRowIndices = -1;
      dofColIndices = -1;

      localPerf = 0;
      localPerfJacobian = 0;

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const resOffset = resDofNumber[er][esr][ei];
      globalIndex const wellElemOffset = wellElemDofNumber[iwelem];

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        eqnRowIndices[WellSolverBase::SubRegionTag::RES * NC + ic] = resOffset + ic;
        eqnRowIndices[WellSolverBase::SubRegionTag::WELL * NC + ic] =
          wellElemOffset + CompositionalMultiphaseWell::RowOffset::MASSBAL + ic;
      }
      for( localIndex jdof = 0; jdof < resNDOF; ++jdof )
      {
        dofColIndices[WellSolverBase::SubRegionTag::RES * resNDOF + jdof] = resOffset + jdof;
        dofColIndices[WellSolverBase::SubRegionTag::WELL * resNDOF + jdof] =
          wellElemOffset + CompositionalMultiphaseWell::ColOffset::DPRES + jdof;
      }

      // populate local flux vector and derivatives
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        localPerf[WellSolverBase::SubRegionTag::RES * NC + ic] = dt * compPerfRate[iperf][ic];
        localPerf[WellSolverBase::SubRegionTag::WELL * NC + ic] = -dt * compPerfRate[iperf][ic];

        for( localIndex ke = 0; ke < 2; ++ke )
        {
          localIndex const localDofIndexPres = ke * resNDOF;
          localPerfJacobian[WellSolverBase::SubRegionTag::RES * NC + ic][localDofIndexPres] =
            dt * dCompPerfRate_dPres[iperf][ke][ic];
          localPerfJacobian[WellSolverBase::SubRegionTag::WELL * NC + ic][localDofIndexPres] =
            -dt * dCompPerfRate_dPres[iperf][ke][ic];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
            localPerfJacobian[WellSolverBase::SubRegionTag::RES * NC + ic][localDofIndexComp] =
              dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
            localPerfJacobian[WellSolverBase::SubRegionTag::WELL * NC + ic][localDofIndexComp] =
              -dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
          }
        }
      }

      rhs->add( eqnRowIndices,
                localPerf );

      matrix->add( eqnRowIndices,
                   dofColIndices,
                   localPerfJacobian );
    }
  } );

}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseReservoir, std::string const &, Group * const )

} /* namespace geosx */
