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

#include "mesh/PerforationExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geosx
{

namespace coupledReservoirAndWellsInternal
{

void
initializePostInitialConditionsPreSubGroups( SolverBase * const solver )
{

  DomainPartition & domain = solver->getGroupByPath< DomainPartition >( "/Problem/domain" );

  solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                        MeshLevel & meshLevel,
                                                                        arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    // loop over the wells
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion & subRegion )
    {
      array1d< array1d< arrayView3d< real64 const > > > const permeability =
        elemManager.constructMaterialExtrinsicAccessor< constitutive::PermeabilityBase, extrinsicMeshData::permeability::permeability >();

      PerforationData * const perforationData = subRegion.getPerforationData();

      // compute the Peaceman index (if not read from XML)
      perforationData->computeWellTransmissibility( meshLevel, subRegion, permeability );
    } );
  } );
}

void
addCouplingNumNonzeros( SolverBase const * const solver,
                        DomainPartition & domain,
                        DofManager & dofManager,
                        arrayView1d< localIndex > const & rowLengths,
                        integer const resNumDof,
                        integer const wellNumDof,
                        string const & resElemDofName,
                        string const & wellElemDofName )
{
  solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                        MeshLevel const & meshLevel,
                                                                        arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = meshLevel.getElemManager();

    string const wellDofKey = dofManager.getKey( wellElemDofName );
    string const resDofKey = dofManager.getKey( resElemDofName );

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
        perforationData->getExtrinsicData< extrinsicMeshData::perforation::wellElementIndex >();

      // get the element region, subregion, index
      arrayView1d< localIndex const > const & resElementRegion =
        perforationData->getExtrinsicData< extrinsicMeshData::perforation::reservoirElementRegion >();
      arrayView1d< localIndex const > const & resElementSubRegion =
        perforationData->getExtrinsicData< extrinsicMeshData::perforation::reservoirElementSubRegion >();
      arrayView1d< localIndex const > const & resElementIndex =
        perforationData->getExtrinsicData< extrinsicMeshData::perforation::reservoirElementIndex >();

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
          GEOSX_ASSERT_GE( rowLengths.size(), localRow + resNumDof );

          for( integer idof = 0; idof < resNumDof; ++idof )
          {
            rowLengths[localRow + idof] += wellNumDof;
          }
        }

        if( wellElemGhostRank[iwelem] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( wellElemDofNumber[iwelem] - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GE( rowLengths.size(), localRow + wellNumDof );

          for( integer idof = 0; idof < wellNumDof; ++idof )
          {
            rowLengths[localRow + idof] += resNumDof;
          }
        }
      } );
    } );
  } );
}

}

} /* namespace geosx */
