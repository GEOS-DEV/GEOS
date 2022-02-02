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
 * @file CompositionalMultiphaseReservoir.cpp
 *
 */


#include "CompositionalMultiphaseReservoir.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

CompositionalMultiphaseReservoir::CompositionalMultiphaseReservoir( const string & name,
                                                                    Group * const parent ):
  ReservoirSolverBase( name, parent )
{
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirFVM;
}

CompositionalMultiphaseReservoir::~CompositionalMultiphaseReservoir()
{}

void CompositionalMultiphaseReservoir::postProcessInput()
{
  ReservoirSolverBase::postProcessInput();

  integer const & useMassFlow = m_flowSolver->getReference< integer >( CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString() );
  integer const & useMassWell = m_wellSolver->getReference< integer >( CompositionalMultiphaseWell::viewKeyStruct::useMassFlagString() );
  GEOSX_THROW_IF( useMassFlow != useMassWell,
                  GEOSX_FMT( "CompositionalMultiphaseReservoir '{}': the input flag {} must be the same in the flow and well solvers, respectively '{}' and '{}'",
                             getName(), CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString(),
                             m_flowSolver->getName(), m_wellSolver->getName() ),
                  InputError );
}

void CompositionalMultiphaseReservoir::initializePostInitialConditionsPreSubGroups()
{
  ReservoirSolverBase::initializePostInitialConditionsPreSubGroups();

  if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseReservoirHybridFVM;
  }
}

void CompositionalMultiphaseReservoir::computeStatistics( real64 const & time,
                                                          real64 const & dt,
                                                          integer cycleNumber,
                                                          DomainPartition & domain,
                                                          bool outputStatisticsToScreen )
{
  // output the number of Newton iterations if this is the main solver
  if( outputStatisticsToScreen && m_nonlinearSolverParameters.m_totalSuccessfulNewtonNumIterations > 0 )
  {
    GEOSX_LOG_LEVEL_RANK_0( 1, getName()
                            << ": Total number of time steps = " << cycleNumber+1
                            << ", successful nonlinear iterations = " << m_nonlinearSolverParameters.m_totalSuccessfulNewtonNumIterations
                            << ", wasted nonlinear iterations = " << m_nonlinearSolverParameters.m_totalWastedNewtonNumIterations );
  }

  m_flowSolver->computeStatistics( time, dt, cycleNumber, domain, outputStatisticsToScreen );
}

void CompositionalMultiphaseReservoir::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                   DofManager const & dofManager,
                                                                   SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  // TODO: remove this and just call SolverBase::setupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  localIndex const resNDOF = m_wellSolver->numDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->numDofPerWellElement();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );
  string const resDofKey  = dofManager.getKey( m_wellSolver->resElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resDofNumber =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

  globalIndex const rankOffset = dofManager.rankOffset();

  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion const & subRegion )
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
      stackArray1d< globalIndex, maxNumDof > eqnRowIndicesRes( resNDOF );
      stackArray1d< globalIndex, maxNumDof > eqnRowIndicesWell( wellNDOF );
      stackArray1d< globalIndex, maxNumDof > dofColIndicesRes( resNDOF );
      stackArray1d< globalIndex, maxNumDof > dofColIndicesWell( wellNDOF );

      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];
      localIndex const iwelem = perfWellElemIndex[iperf];

      for( localIndex idof = 0; idof < resNDOF; ++idof )
      {
        eqnRowIndicesRes[idof] = resDofNumber[er][esr][ei] + idof - rankOffset;
        dofColIndicesRes[idof] = resDofNumber[er][esr][ei] + idof;
      }

      for( localIndex idof = 0; idof < wellNDOF; ++idof )
      {
        eqnRowIndicesWell[idof] = wellElemDofNumber[iwelem] + idof - rankOffset;
        dofColIndicesWell[idof] = wellElemDofNumber[iwelem] + idof;
      }

      for( localIndex i = 0; i < eqnRowIndicesRes.size(); ++i )
      {
        if( eqnRowIndicesRes[i] >= 0 && eqnRowIndicesRes[i] < pattern.numRows() )
        {
          for( localIndex j = 0; j < dofColIndicesWell.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesRes[i], dofColIndicesWell[j] );
          }
        }
      }

      for( localIndex i = 0; i < eqnRowIndicesWell.size(); ++i )
      {
        if( eqnRowIndicesWell[i] >= 0 && eqnRowIndicesWell[i] < pattern.numRows() )
        {
          for( localIndex j = 0; j < dofColIndicesRes.size(); ++j )
          {
            pattern.insertNonZero( eqnRowIndicesWell[i], dofColIndicesRes[j] );
          }
        }
      }
    } );
  } );
}

void CompositionalMultiphaseReservoir::assembleCouplingTerms( real64 const time_n,
                                                              real64 const dt,
                                                              DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  using namespace CompositionalMultiphaseUtilities;

  using TAG = CompositionalMultiphaseWellKernels::SubRegionTag;
  using ROFFSET = CompositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = CompositionalMultiphaseWellKernels::ColOffset;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  localIndex constexpr MAX_NUM_COMP = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr MAX_NUM_DOF  = MAX_NUM_COMP + 1;

  localIndex const numComps = m_wellSolver->numFluidComponents();
  localIndex const resNumDofs  = m_wellSolver->numDofPerResElement();

  string const resDofKey = dofManager.getKey( m_wellSolver->resElementDofName() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const resDofNumberAccessor =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );
  ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber =
    resDofNumberAccessor.toNestedViewConst();
  globalIndex const rankOffset = dofManager.rankOffset();

  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion const & subRegion )
  {

    // if the well is shut, we neglect reservoir-well flow that may occur despite the zero rate
    // therefore, we do not want to compute perforation rates and we simply assume they are zero
    WellControls const & wellControls = m_wellSolver->getWellControls( subRegion );
    if( !wellControls.wellIsOpen( time_n + dt ) )
    {
      return;
    }

    PerforationData const * const perforationData = subRegion.getPerforationData();

    // get the degrees of freedom
    string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get well variables on perforations
    arrayView2d< real64 const > const & compPerfRate =
      perforationData->getExtrinsicData< extrinsicMeshData::well::compPerforationRate >();
    arrayView3d< real64 const > const & dCompPerfRate_dPres =
      perforationData->getExtrinsicData< extrinsicMeshData::well::dCompPerforationRate_dPres >();
    arrayView4d< real64 const > const & dCompPerfRate_dComp =
      perforationData->getExtrinsicData< extrinsicMeshData::well::dCompPerforationRate_dComp >();

    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

    // loop over the perforations and add the rates to the residual and jacobian
    forAll< parallelDevicePolicy<> >( perforationData->size(), [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
    {
      // local working variables and arrays
      stackArray1d< localIndex, 2 * MAX_NUM_COMP > eqnRowIndices( 2 * numComps );
      stackArray1d< globalIndex, 2 * MAX_NUM_DOF > dofColIndices( 2 * resNumDofs );

      stackArray1d< real64, 2 * MAX_NUM_COMP > localPerf( 2 * numComps );
      stackArray2d< real64, 2 * MAX_NUM_COMP * 2 * MAX_NUM_DOF > localPerfJacobian( 2 * numComps, 2 * resNumDofs );

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const resOffset = resDofNumber[er][esr][ei];
      globalIndex const wellElemOffset = wellElemDofNumber[iwelem];

      for( localIndex ic = 0; ic < numComps; ++ic )
      {
        eqnRowIndices[TAG::RES * numComps + ic] = LvArray::integerConversion< localIndex >( resOffset - rankOffset ) + ic;
        eqnRowIndices[TAG::WELL * numComps + ic] = LvArray::integerConversion< localIndex >( wellElemOffset - rankOffset ) + ROFFSET::MASSBAL + ic;
      }
      for( localIndex jdof = 0; jdof < resNumDofs; ++jdof )
      {
        dofColIndices[TAG::RES * resNumDofs + jdof] = resOffset + jdof;
        dofColIndices[TAG::WELL * resNumDofs + jdof] = wellElemOffset + COFFSET::DPRES + jdof;
      }

      // populate local flux vector and derivatives
      for( localIndex ic = 0; ic < numComps; ++ic )
      {
        localPerf[TAG::RES * numComps + ic] = dt * compPerfRate[iperf][ic];
        localPerf[TAG::WELL * numComps + ic] = -dt * compPerfRate[iperf][ic];

        for( localIndex ke = 0; ke < 2; ++ke )
        {
          localIndex const localDofIndexPres = ke * resNumDofs;
          localPerfJacobian[TAG::RES * numComps + ic][localDofIndexPres] = dt * dCompPerfRate_dPres[iperf][ke][ic];
          localPerfJacobian[TAG::WELL * numComps + ic][localDofIndexPres] = -dt * dCompPerfRate_dPres[iperf][ke][ic];

          for( localIndex jc = 0; jc < numComps; ++jc )
          {
            localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
            localPerfJacobian[TAG::RES * numComps + ic][localDofIndexComp] = dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
            localPerfJacobian[TAG::WELL * numComps + ic][localDofIndexComp] = -dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
          }
        }
      }

      // Apply equation/variable change transformation(s)
      stackArray1d< real64, 2 * MAX_NUM_DOF > work( 2 * resNumDofs );
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComps, resNumDofs*2, 2, localPerfJacobian, work );
      shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComps, 2, localPerf );

      for( localIndex i = 0; i < localPerf.size(); ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                            dofColIndices.data(),
                                                                            localPerfJacobian[i].dataIfContiguous(),
                                                                            2 * resNumDofs );
          atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localPerf[i] );
        }
      }
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseReservoir, string const &, Group * const )

} /* namespace geosx */
