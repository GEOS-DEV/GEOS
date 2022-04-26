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
 * @file CoupledReservoirAndWellsBase.cpp
 *
 */

#include "CoupledReservoirAndWellsBase.hpp"

#include "common/TimingMacros.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
CoupledReservoirAndWellsBase( const string & name,
                              Group * const parent )
  : AbstractBase( name, parent ),
  m_reservoirSolverName(),
  m_wellSolverName()
{
  this->template registerWrapper( viewKeyStruct::reservoirSolverNameString(), &m_reservoirSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the flow solver to use in the reservoir-well system solver" );

  this->template registerWrapper( viewKeyStruct::wellSolverNameString(), &m_wellSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the well solver to use in the reservoir-well system solver" );

  this->template getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
~CoupledReservoirAndWellsBase()
{}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
RESERVOIR_SOLVER *
CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
getReservoirSolver() const { return std::get< toUnderlying< SolverType >( SolverType::Reservoir ) >( m_solvers ); }

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
WELL_SOLVER *
CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
getWellSolver() const { return std::get< toUnderlying< SolverType >( SolverType::Well ) >( m_solvers ); }

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
void CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
initializePostInitialConditionsPreSubGroups()
{
  AbstractBase::initializePostInitialConditionsPreSubGroups( );

  DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );

  this->template forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                               MeshLevel & meshLevel,
                                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    // loop over the wells
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion & subRegion )
    {
      array1d< array1d< arrayView3d< real64 const > > > const permeability =
        elemManager.constructMaterialExtrinsicAccessor< PermeabilityBase, extrinsicMeshData::permeability::permeability >();

      PerforationData * const perforationData = subRegion.getPerforationData();

      // compute the Peaceman index (if not read from XML)
      perforationData->computeWellTransmissibility( meshLevel, subRegion, permeability );
    } );
  } );
}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
void CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
addCouplingNumNonzeros( DomainPartition & domain,
                        DofManager & dofManager,
                        arrayView1d< localIndex > const & rowLengths ) const
{
  integer const resNDOF = getWellSolver()->numDofPerResElement();
  integer const wellNDOF = getWellSolver()->numDofPerWellElement();

  this->template forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                               MeshLevel const & meshLevel,
                                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = meshLevel.getElemManager();

    string const wellDofKey = dofManager.getKey( getWellSolver()->wellElementDofName() );
    string const resDofKey = dofManager.getKey( getWellSolver()->resElementDofName() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resElemDofNumber =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const & resElemGhostRank =
      elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

    globalIndex const rankOffset = dofManager.rankOffset();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const, WellElementSubRegion const & subRegion )
    {
      PerforationData const * const perforationData = subRegion.getPerforationData();

      arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

      // get the well degrees of freedom and ghosting info
      arrayView1d< globalIndex const > const & wellElemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( wellDofKey );

      // get the well element indices corresponding to each perforation
      arrayView1d< localIndex const > const & perfWellElemIndex =
        perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion = perforationData->getMeshElements().m_toElementRegion;
      arrayView1d< localIndex const > const & resElementSubRegion = perforationData->getMeshElements().m_toElementSubRegion;
      arrayView1d< localIndex const > const & resElementIndex = perforationData->getMeshElements().m_toElementIndex;

      // Loop over perforations and increase row lengths for reservoir and well elements accordingly
      forAll< serialPolicy >( perforationData->size(), [=] ( localIndex const iperf )
      {
        // get the reservoir (sub)region and element indices
        localIndex const er = resElementRegion[iperf];
        localIndex const esr = resElementSubRegion[iperf];
        localIndex const ei = resElementIndex[iperf];
        localIndex const iwelem = perfWellElemIndex[iperf];

        if( resElemGhostRank[er][esr][ei] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( resElemDofNumber[er][esr][ei] - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GE( rowLengths.size(), localRow + resNDOF );

          for( integer idof = 0; idof < resNDOF; ++idof )
          {
            rowLengths[localRow + idof] += wellNDOF;
          }
        }

        if( wellElemGhostRank[iwelem] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( wellElemDofNumber[iwelem] - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GE( rowLengths.size(), localRow + wellNDOF );

          for( integer idof = 0; idof < wellNDOF; ++idof )
          {
            rowLengths[localRow + idof] += resNDOF;
          }
        }
      } );
    } );
  } );
}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
void CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
setupSystem( DomainPartition & domain,
             DofManager & dofManager,
             CRSMatrix< real64, globalIndex > & localMatrix,
             ParallelVector & rhs,
             ParallelVector & solution,
             bool const )
{
  GEOSX_MARK_FUNCTION;

  dofManager.setDomain( domain );

  AbstractBase::setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  // Set the sparsity pattern without reservoir-well coupling
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling on perforations
  addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/localMatrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOSX );
}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
void CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
assembleSystem( real64 const time_n,
                real64 const dt,
                DomainPartition & domain,
                DofManager const & dofManager,
                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                arrayView1d< real64 > const & localRhs )
{
  // assemble J_RR (excluding perforation rates)
  getReservoirSolver()->assembleSystem( time_n, dt,
                                        domain,
                                        dofManager,
                                        localMatrix,
                                        localRhs );

  // assemble J_WW (excluding perforation rates)
  getWellSolver()->assembleSystem( time_n, dt,
                                   domain,
                                   dofManager,
                                   localMatrix,
                                   localRhs );

  // assemble perforation rates in J_WR, J_RW, J_RR and J_WW
  assembleCouplingTerms( time_n, dt,
                         domain,
                         dofManager,
                         localMatrix,
                         localRhs );
}

template< typename RESERVOIR_SOLVER, typename WELL_SOLVER >
void CoupledReservoirAndWellsBase< RESERVOIR_SOLVER, WELL_SOLVER >::
solveLinearSystem( DofManager const & dofManager,
                   ParallelMatrix & matrix,
                   ParallelVector & rhs,
                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();
  SolverBase::solveLinearSystem( dofManager, matrix, rhs, solution );
}

template class CoupledReservoirAndWellsBase< SinglePhaseFVM< SinglePhaseBase >, SinglePhaseWell >;
template class CoupledReservoirAndWellsBase< SinglePhaseHybridFVM, SinglePhaseWell >;
template class CoupledReservoirAndWellsBase< CompositionalMultiphaseFVM, CompositionalMultiphaseWell >;
template class CoupledReservoirAndWellsBase< CompositionalMultiphaseHybridFVM, CompositionalMultiphaseWell >;


} /* namespace geosx */
