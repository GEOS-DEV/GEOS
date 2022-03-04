/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseReservoir.cpp
 *
 */

#include "SinglePhaseReservoir.hpp"

#include "common/TimingMacros.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhaseReservoir::SinglePhaseReservoir( const string & name,
                                            Group * const parent ):
  ReservoirSolverBase( name, parent )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirFVM;
}

SinglePhaseReservoir::~SinglePhaseReservoir()
{}

void SinglePhaseReservoir::initializePostInitialConditionsPreSubGroups()
{
  ReservoirSolverBase::initializePostInitialConditionsPreSubGroups();

  if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirHybridFVM;
  }
}

void SinglePhaseReservoir::addCouplingSparsityPattern( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    // TODO: remove this and just call SolverBase::setupSystem when DofManager can handle the coupling

    // Populate off-diagonal sparsity between well and reservoir

    string const resDofKey  = dofManager.getKey( m_wellSolver->resElementDofName() );
    string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );

    localIndex const wellNDOF = m_wellSolver->numDofPerWellElement();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resDofNumber =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

    globalIndex const rankOffset = dofManager.rankOffset();

    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion const & subRegion )
    {
      PerforationData const * const perforationData = subRegion.getPerforationData();

      // get the well degrees of freedom and ghosting info
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the well element indices corresponding to each perforation
      arrayView1d< localIndex const > const & perfWellElemIndex =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

      // Insert the entries corresponding to reservoir-well perforations
      // This will fill J_WR, and J_RW
      forAll< serialPolicy >( perforationData->size(), [=] ( localIndex const iperf )
      {
        // Get the reservoir (sub)region and element indices
        localIndex const er = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei = resElementIndex[iperf];
        localIndex const iwelem = perfWellElemIndex[iperf];

        globalIndex const eqnRowIndexRes = resDofNumber[er][esr][ei] - rankOffset;
        globalIndex const dofColIndexRes = resDofNumber[er][esr][ei];

        // working arrays
        stackArray1d< globalIndex, 2 > eqnRowIndicesWell( wellNDOF );
        stackArray1d< globalIndex, 2 > dofColIndicesWell( wellNDOF );

        for( localIndex idof = 0; idof < wellNDOF; ++idof )
        {
          eqnRowIndicesWell[idof] = wellElemDofNumber[iwelem] + idof - rankOffset;
          dofColIndicesWell[idof] = wellElemDofNumber[iwelem] + idof;
        }

        if( eqnRowIndexRes >= 0 && eqnRowIndexRes < pattern.numRows() )
        {
          for( localIndex j = 0; j < dofColIndicesWell.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndexRes, dofColIndicesWell[j] );
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesWell.size(); ++i )
        {
          if( eqnRowIndicesWell[i] >= 0 && eqnRowIndicesWell[i] < pattern.numRows() )
          {
            pattern.insertNonZero( eqnRowIndicesWell[i], dofColIndexRes );
          }
        }
      } );
    } );
  } );
}

void SinglePhaseReservoir::assembleCouplingTerms( real64 const time_n,
                                                  real64 const dt,
                                                  DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  using TAG = singlePhaseWellKernels::SubRegionTag;
  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const resDofKey = dofManager.getKey( m_wellSolver->resElementDofName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const resDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );
    ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber =
      resDofNumberAccessor.toNestedViewConst();
    globalIndex const rankOffset = dofManager.rankOffset();

    // loop over the wells
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion const & subRegion )
    {

      // if the well is shut, we neglect reservoir-well flow that may occur despite the zero rate
      // therefore, we do not want to compute perforation rates and we simply assume they are zero
      WellControls const & wellControls = m_wellSolver->getWellControls( subRegion );
      if( !wellControls.isWellOpen( time_n + dt ) )
      {
        return;
      }

      PerforationData const * const perforationData = subRegion.getPerforationData();

      // get the degrees of freedom
      string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );
      arrayView1d< globalIndex const > const wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get well variables on perforations
      arrayView1d< real64 const > const perfRate =
        perforationData->getExtrinsicData< extrinsicMeshData::well::perforationRate >();
      arrayView2d< real64 const > const dPerfRate_dPres =
        perforationData->getExtrinsicData< extrinsicMeshData::well::dPerforationRate_dPres >();

      arrayView1d< localIndex const > const perfWellElemIndex =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

      // get the element region, subregion, index
      arrayView1d< localIndex const > const resElementRegion =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
      arrayView1d< localIndex const > const resElementSubRegion =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
      arrayView1d< localIndex const > const resElementIndex =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

      // loop over the perforations and add the rates to the residual and jacobian
      forAll< parallelDevicePolicy<> >( perforationData->size(), [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
      {
        // local working variables and arrays
        localIndex eqnRowIndices[ 2 ] = { -1 };
        globalIndex dofColIndices[ 2 ] = { -1 };


        real64 localPerf[ 2 ] = { 0.0 };
        real64 localPerfJacobian[ 2 ][ 2 ] = {{ 0.0 }};

        // get the reservoir (sub)region and element indices
        localIndex const er = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei = resElementIndex[iperf];

        // get the well element index for this perforation
        localIndex const iwelem = perfWellElemIndex[iperf];
        globalIndex const elemOffset = wellElemDofNumber[iwelem];

        // row index on reservoir side
        eqnRowIndices[TAG::RES] = resDofNumber[er][esr][ei] - rankOffset;
        // column index on reservoir side
        dofColIndices[TAG::RES] = resDofNumber[er][esr][ei];

        // row index on well side
        eqnRowIndices[TAG::WELL] = LvArray::integerConversion< localIndex >( elemOffset - rankOffset ) + ROFFSET::MASSBAL;
        // column index on well side
        dofColIndices[TAG::WELL] = elemOffset + COFFSET::DPRES;

        // populate local flux vector and derivatives
        localPerf[TAG::RES] = dt * perfRate[iperf];
        localPerf[TAG::WELL] = -localPerf[TAG::RES];

        for( localIndex ke = 0; ke < 2; ++ke )
        {
          localPerfJacobian[TAG::RES][ke] = dt * dPerfRate_dPres[iperf][ke];
          localPerfJacobian[TAG::WELL][ke] = -localPerfJacobian[TAG::RES][ke];
        }

        for( localIndex i = 0; i < 2; ++i )
        {
          if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                              &dofColIndices[0],
                                                                              &localPerfJacobian[0][0] + 2 * i,
                                                                              2 );
            atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localPerf[i] );
          }
        }
      } );
    } );
  } );

}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseReservoir, string const &, Group * const )

} /* namespace geosx */
