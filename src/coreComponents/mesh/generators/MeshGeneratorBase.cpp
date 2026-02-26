/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "MeshGeneratorBase.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/LogLevelsInfo.hpp"
#include "mesh/generators/ParticleBlockManager.hpp"
#include "mesh/generators/MeshComponentBase.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/logger/Logger.hpp"
#include "common/MpiWrapper.hpp"
namespace geos
{
using namespace dataRepository;

MeshGeneratorBase::MeshGeneratorBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::checkEulerCharacteristicString(), &m_checkEulerCharacteristic ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "When set to 1, compute the Euler-Poincaré characteristic χ = V − E + F − C "
                    "after the mesh is loaded and issue a warning if χ ≠ 1. "
                    "χ = 1 is expected for a single connected solid mesh without interior voids. "
                    "Non-unity values may indicate multiple disconnected bodies, interior voids, "
                    "or non-manifold topology. Default: 0 (disabled)." );
}

Group * MeshGeneratorBase::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< MeshComponentBase > meshComp =
    MeshComponentBase::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &this->registerGroup< MeshComponentBase >( childName, std::move( meshComp ) );
}

void MeshGeneratorBase::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from MeshComponentBase here
  for( auto & catalogIter: MeshComponentBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

MeshGeneratorBase::CatalogInterface::CatalogType & MeshGeneratorBase::getCatalog()
{
  static MeshGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void MeshGeneratorBase::generateMesh( Group & parent, SpatialPartition & partition )
{
  MeshBody & meshBody = dynamic_cast< MeshBody & >( parent );
  if( meshBody.hasParticles() )
  {
    ParticleBlockManager & particleBlockManager = parent.registerGroup< ParticleBlockManager >( keys::particleManager );

    MeshLevel & meshLevel0 = meshBody.getBaseDiscretization();
    ParticleManager & particleManager = meshLevel0.getParticleManager();

    fillParticleBlockManager( particleBlockManager, particleManager, partition );
  }
  else
  {
    CellBlockManager & cellBlockManager = parent.registerGroup< CellBlockManager >( keys::cellManager );

    fillCellBlockManager( cellBlockManager, partition );

    this->attachWellInfo( cellBlockManager );
  }
}

void MeshGeneratorBase::attachWellInfo( CellBlockManager & cellBlockManager )
{
  forSubGroups< WellGeneratorBase >( [&]( WellGeneratorBase & wellGen ) {
    wellGen.generateWellGeometry( );
    LineBlock & lb = cellBlockManager.registerLineBlock( wellGen.getWellRegionName() );
    lb.setNumElements( wellGen.numElements() );
    lb.setElemCoords( wellGen.getElemCoords() );
    lb.setNextElemIndex( wellGen.getNextElemIndex() );
    lb.setPrevElemIndices( wellGen.getPrevElemIndices() );
    lb.setElemToNodesMap( wellGen.getElemToNodesMap() );
    lb.setElemVolume( wellGen.getElemVolume() );
    lb.setElementRadius( wellGen.getElementRadius() );
    lb.setNumNodes( wellGen.numNodes() );
    lb.setNodeCoords( wellGen.getNodeCoords() );
    lb.setNumPerforations( wellGen.numPerforations() );
    lb.setPerfCoords( wellGen.getPerfCoords() );
    lb.setPerfTransmissibility( wellGen.getPerfTransmissibility() );
    lb.setPerfSkinFactor( wellGen.getPerfSkinFactor() );
    lb.setPerfTargetRegion( wellGen.getPerfTargetRegion() );
    lb.setPerfElemIndex( wellGen.getPerfElemIndex() );
    lb.setPerfStatusTableName( wellGen.getPerfStatusTableName());
    lb.setPerfName( wellGen.getPerfName());
    lb.setWellControlsName( wellGen.getWellControlsName() );
    lb.setWellGeneratorName( wellGen.getName() );

  } );
}

integer computeEulerCharacteristic( NodeManager const & nodeManager,
                                    EdgeManager const & edgeManager,
                                    FaceManager const & faceManager,
                                    ElementRegionManager const & elemManager )
{
  // Obtain global unique entity counts for the Euler characteristic
  // χ = V − E + F − C.
  //
  // This function may be called before ghost exchange (setupCommunications),
  // so ghost ranks are not necessarily set for every entity type.
  //
  // Strategy per entity type:
  //
  //  Nodes – ghost ranks are NOT set at this point (assignGlobalIndices is
  //          not called for nodes in setupBaseLevelMeshGlobalInfo).  However,
  //          node global IDs originate from the VTK mesh and are a contiguous
  //          0-based sequence [0 .. N-1].  Therefore maxGlobalIndex()+1 is the
  //          true global node count.
  //
  //  Edges / Faces – assignGlobalIndices has been called, which properly sets
  //          ghost ranks for duplicate copies.  getNumberOfLocalIndices()
  //          counts only owned entities (ghostRank <= -1), so the MPI sum
  //          gives the correct global unique count.
  //
  //  Cells – each cell is assigned to exactly one MPI rank by the partitioner
  //          (no duplication before ghost exchange).  size() on each rank is
  //          the owned count and the MPI sum is correct.

  // --- Nodes: use maxGlobalIndex (contiguous 0-based VTK global IDs) ---
  globalIndex const V = nodeManager.maxGlobalIndex() + 1;

  // --- Edges & Faces: ghost ranks correctly set by assignGlobalIndices ---
  localIndex const localE = edgeManager.getNumberOfLocalIndices();
  localIndex const localF = faceManager.getNumberOfLocalIndices();
  globalIndex const E = MpiWrapper::sum< globalIndex >( localE );
  globalIndex const F = MpiWrapper::sum< globalIndex >( localF );

  // --- Cells: uniquely partitioned before ghost exchange; use
  //     getNumberOfLocalIndices() so this also works after ghost exchange. ---
  localIndex localC = 0;
  elemManager.forElementSubRegions< CellElementSubRegion >(
    [&]( CellElementSubRegion const & subRegion )
  {
    localC += subRegion.getNumberOfLocalIndices();
  } );
  globalIndex const C = MpiWrapper::sum< globalIndex >( localC );

  return static_cast< integer >( V - E + F - C );
}

}
