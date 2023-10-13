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
 * @file SinglePhaseReservoirAndWells.cpp
 *
 */

#include "SinglePhaseReservoirAndWells.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/PerforationFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

namespace
{

// This is meant to be specialized to work, see below
template< typename SINGLEPHASE_RESERVOIR_SOLVER > class
  SinglePhaseCatalogNames {};

// Class specialization for a RESERVOIR_SOLVER set to SinglePhaseFlow
template<> class SinglePhaseCatalogNames< SinglePhaseBase >
{
public:
  // TODO: find a way to use the catalog name here
  static string name() { return "SinglePhaseReservoir"; }
};
// Class specialization for a RESERVOIR_SOLVER set to SinglePhasePoromechanics
template<> class SinglePhaseCatalogNames< SinglePhasePoromechanics >
{
public:
  static string name() { return SinglePhasePoromechanics::catalogName()+"Reservoir"; }
};
}

// provide a definition for catalogName()
template< typename SINGLEPHASE_RESERVOIR_SOLVER >
string
SinglePhaseReservoirAndWells< SINGLEPHASE_RESERVOIR_SOLVER >::
catalogName()
{
  return SinglePhaseCatalogNames< SINGLEPHASE_RESERVOIR_SOLVER >::name();
}

template< typename SINGLEPHASE_RESERVOIR_SOLVER >
SinglePhaseReservoirAndWells< SINGLEPHASE_RESERVOIR_SOLVER >::
SinglePhaseReservoirAndWells( const string & name,
                              Group * const parent )
  : Base( name, parent )
{}

template< typename SINGLEPHASE_RESERVOIR_SOLVER >
SinglePhaseReservoirAndWells< SINGLEPHASE_RESERVOIR_SOLVER >::
~SinglePhaseReservoirAndWells()
{}

template<>
SinglePhaseBase const *
SinglePhaseReservoirAndWells< SinglePhaseBase >::
flowSolver() const
{
  return this->reservoirSolver();
}

template<>
SinglePhaseBase const *
SinglePhaseReservoirAndWells< SinglePhasePoromechanics >::
flowSolver() const
{
  return this->reservoirSolver()->flowSolver();
}

template<>
void
SinglePhaseReservoirAndWells< SinglePhaseBase >::
setMGRStrategy()
{
  if( flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirHybridFVM;
  }
  else
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirFVM;
  }
}

template<>
void
SinglePhaseReservoirAndWells< SinglePhasePoromechanics >::
setMGRStrategy()
{
  if( flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseReservoirHybridFVM )
  {
    GEOS_LOG_RANK_0( "The MGR strategy for hybrid FVM is not implemented" );
  }
  else
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsReservoirFVM;
  }
}

template< typename SINGLEPHASE_RESERVOIR_SOLVER >
void
SinglePhaseReservoirAndWells< SINGLEPHASE_RESERVOIR_SOLVER >::
initializePreSubGroups()
{
  Base::initializePreSubGroups();
  SinglePhaseBase const * const flowSolver = this->flowSolver();
  Base::wellSolver()->setFlowSolverName( flowSolver->getName() );
}

template< typename SINGLEPHASE_RESERVOIR_SOLVER >
void
SinglePhaseReservoirAndWells< SINGLEPHASE_RESERVOIR_SOLVER >::
initializePostInitialConditionsPreSubGroups()
{
  Base::initializePostInitialConditionsPreSubGroups();
  setMGRStrategy();
}

template< typename SINGLEPHASE_RESERVOIR_SOLVER >
void
SinglePhaseReservoirAndWells< SINGLEPHASE_RESERVOIR_SOLVER >::
addCouplingSparsityPattern( DomainPartition const & domain,
                            DofManager const & dofManager,
                            SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel const & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();

    // TODO: remove this and just call SolverBase::setupSystem when DofManager can handle the coupling

    // Populate off-diagonal sparsity between well and reservoir

    string const resDofKey  = dofManager.getKey( Base::wellSolver()->resElementDofName() );
    string const wellDofKey = dofManager.getKey( Base::wellSolver()->wellElementDofName() );

    integer const wellNDOF = Base::wellSolver()->numDofPerWellElement();

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
        perforationData->getField< fields::perforation::wellElementIndex >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getField< fields::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData->getField< fields::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData->getField< fields::perforation::reservoirElementIndex >();

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

        for( integer idof = 0; idof < wellNDOF; ++idof )
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

template< typename SINGLEPHASE_RESERVOIR_SOLVER >
void
SinglePhaseReservoirAndWells< SINGLEPHASE_RESERVOIR_SOLVER >::
assembleCouplingTerms( real64 const time_n,
                       real64 const dt,
                       DomainPartition const & domain,
                       DofManager const & dofManager,
                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                       arrayView1d< real64 > const & localRhs )
{
  using TAG = singlePhaseWellKernels::SubRegionTag;
  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  GEOS_THROW_IF( !Base::m_isWellTransmissibilityComputed,
                 GEOS_FMT( "{} `{}`: The well transmissibility has not been computed yet",
                           catalogName(), this->getName() ),
                 std::runtime_error );

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel const & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    integer areWellsShut = 1;

    ElementRegionManager const & elemManager = mesh.getElemManager();

    string const resDofKey = dofManager.getKey( Base::wellSolver()->resElementDofName() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const resDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );
    ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber =
      resDofNumberAccessor.toNestedViewConst();
    globalIndex const rankOffset = dofManager.rankOffset();

    // loop over the wells
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion const & subRegion )
    {
      PerforationData const * const perforationData = subRegion.getPerforationData();

      WellControls const & wellControls = Base::wellSolver()->getWellControls( subRegion );
      if( !wellControls.isWellOpen( time_n + dt ) )
      {
        return;
      }

      areWellsShut = 0;

      // get the degrees of freedom
      string const wellDofKey = dofManager.getKey( Base::wellSolver()->wellElementDofName() );
      arrayView1d< globalIndex const > const wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get well variables on perforations
      arrayView1d< real64 const > const perfRate =
        perforationData->getField< fields::well::perforationRate >();
      arrayView2d< real64 const > const dPerfRate_dPres =
        perforationData->getField< fields::well::dPerforationRate_dPres >();

      arrayView1d< localIndex const > const perfWellElemIndex =
        perforationData->getField< fields::perforation::wellElementIndex >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const resElementRegion =
        perforationData->getField< fields::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const resElementSubRegion =
        perforationData->getField< fields::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const resElementIndex =
        perforationData->getField< fields::perforation::reservoirElementIndex >();

      // loop over the perforations and add the rates to the residual and jacobian
      forAll< parallelDevicePolicy<> >( perforationData->size(), [=] GEOS_HOST_DEVICE ( localIndex const iperf )
      {
        // local working variables and arrays
        localIndex eqnRowIndices[ 2 ] = { -1 };
        globalIndex dofColIndices[ 2 ] = { -1 };


        real64 localPerf[ 2 ]{};
        real64 localPerfJacobian[ 2 ][ 2 ]{};

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

        for( integer ke = 0; ke < 2; ++ke )
        {
          localPerfJacobian[TAG::RES][ke] = dt * dPerfRate_dPres[iperf][ke];
          localPerfJacobian[TAG::WELL][ke] = -localPerfJacobian[TAG::RES][ke];
        }

        for( integer i = 0; i < 2; ++i )
        {
          if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                              &dofColIndices[0],
                                                                              &localPerfJacobian[0][0] + 2 * i,
                                                                              2 );
            RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localPerf[i] );
          }
        }
      } );
    } );

    // update dynamically the MGR recipe to optimize the linear solve if all wells are shut
    areWellsShut = MpiWrapper::min( areWellsShut );
    m_linearSolverParameters.get().mgr.areWellsShut = areWellsShut;


  } );

}

namespace
{
typedef SinglePhaseReservoirAndWells< SinglePhaseBase > SinglePhaseFlowAndWells;
typedef SinglePhaseReservoirAndWells< SinglePhasePoromechanics > SinglePhasePoromechanicsAndWells;
REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlowAndWells, string const &, Group * const )
REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsAndWells, string const &, Group * const )
}

} /* namespace geos */
