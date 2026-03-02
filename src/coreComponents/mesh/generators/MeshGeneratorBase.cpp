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
#include "mesh/CellElementSubRegion.hpp"

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

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
                                    EdgeManager const & GEOS_UNUSED_PARAM( edgeManager ),
                                    FaceManager const & GEOS_UNUSED_PARAM( faceManager ),
                                    ElementRegionManager const & elemManager )
{
  // Compute the Euler-Poincaré characteristic χ = V − E + F − C for the
  // **bulk** mesh only (CellElementSubRegion entities).
  //
  // We reconstruct V, E, F directly from the cell-to-node maps so that
  // fracture-induced duplicates in EdgeManager / FaceManager are excluded.
  //
  // MPI correctness: each rank iterates over its *owned* cells only
  // (ghostRank < 0).  Shared faces / edges / nodes on partition boundaries
  // are seen by multiple ranks.  To avoid double-counting we assign
  // ownership of every entity to the rank that owns the node with the
  // smallest global ID in that entity.  The MPI sum of owned counts then
  // gives the true global unique count.

  using EdgeKey = std::pair< globalIndex, globalIndex >;
  using FaceKey = std::vector< globalIndex >;

  arrayView1d< globalIndex const > const localToGlobal = nodeManager.localToGlobalMap();
  arrayView1d< integer const > const nodeGhostRank = nodeManager.ghostRank();

  std::set< globalIndex > allNodes;
  std::set< EdgeKey >     allEdges;
  std::set< FaceKey >     allFaces;
  localIndex              numOwnedCells = 0;

  // --- helpers -----------------------------------------------------------
  auto makeEdge = []( globalIndex a, globalIndex b ) -> EdgeKey
  {
    return a < b ? std::make_pair( a, b ) : std::make_pair( b, a );
  };

  auto makeFace = []( std::initializer_list< globalIndex > nodes ) -> FaceKey
  {
    FaceKey f( nodes );
    std::sort( f.begin(), f.end() );
    return f;
  };

  // --- iterate owned cells, collecting global-ID-based entities ----------
  elemManager.forElementSubRegions< CellElementSubRegion >(
    [&]( CellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap =
      subRegion.nodeList();
    arrayView1d< integer const > const cellGhostRank = subRegion.ghostRank();

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      if( cellGhostRank[ei] >= 0 )
        continue;                       // skip ghost cells

      ++numOwnedCells;

      localIndex const npe = elemToNodeMap.size( 1 );
      std::vector< globalIndex > gn( npe );
      for( localIndex ni = 0; ni < npe; ++ni )
      {
        gn[ni] = localToGlobal[ elemToNodeMap[ei][ni] ];
        allNodes.insert( gn[ni] );
      }

      if( npe == 4 )  // ---- Tetrahedron --------------------------------
      {
        allEdges.insert( makeEdge( gn[0], gn[1] ) );
        allEdges.insert( makeEdge( gn[0], gn[2] ) );
        allEdges.insert( makeEdge( gn[0], gn[3] ) );
        allEdges.insert( makeEdge( gn[1], gn[2] ) );
        allEdges.insert( makeEdge( gn[1], gn[3] ) );
        allEdges.insert( makeEdge( gn[2], gn[3] ) );

        allFaces.insert( makeFace( {gn[0], gn[1], gn[2]} ) );
        allFaces.insert( makeFace( {gn[0], gn[1], gn[3]} ) );
        allFaces.insert( makeFace( {gn[0], gn[2], gn[3]} ) );
        allFaces.insert( makeFace( {gn[1], gn[2], gn[3]} ) );
      }
      else if( npe == 8 )  // ---- Hexahedron (GEOS ordering 0,1,3,2,4,5,7,6) ---
      {
        allEdges.insert( makeEdge( gn[0], gn[1] ) );
        allEdges.insert( makeEdge( gn[1], gn[3] ) );
        allEdges.insert( makeEdge( gn[3], gn[2] ) );
        allEdges.insert( makeEdge( gn[2], gn[0] ) );
        allEdges.insert( makeEdge( gn[4], gn[5] ) );
        allEdges.insert( makeEdge( gn[5], gn[7] ) );
        allEdges.insert( makeEdge( gn[7], gn[6] ) );
        allEdges.insert( makeEdge( gn[6], gn[4] ) );
        allEdges.insert( makeEdge( gn[0], gn[4] ) );
        allEdges.insert( makeEdge( gn[1], gn[5] ) );
        allEdges.insert( makeEdge( gn[3], gn[7] ) );
        allEdges.insert( makeEdge( gn[2], gn[6] ) );

        allFaces.insert( makeFace( {gn[0], gn[1], gn[3], gn[2]} ) );
        allFaces.insert( makeFace( {gn[4], gn[5], gn[7], gn[6]} ) );
        allFaces.insert( makeFace( {gn[0], gn[1], gn[5], gn[4]} ) );
        allFaces.insert( makeFace( {gn[2], gn[3], gn[7], gn[6]} ) );
        allFaces.insert( makeFace( {gn[1], gn[3], gn[7], gn[5]} ) );
        allFaces.insert( makeFace( {gn[0], gn[2], gn[6], gn[4]} ) );
      }
      else if( npe == 6 )  // ---- Wedge / Prism ---------------------------
      {
        allEdges.insert( makeEdge( gn[0], gn[1] ) );
        allEdges.insert( makeEdge( gn[1], gn[2] ) );
        allEdges.insert( makeEdge( gn[2], gn[0] ) );
        allEdges.insert( makeEdge( gn[3], gn[4] ) );
        allEdges.insert( makeEdge( gn[4], gn[5] ) );
        allEdges.insert( makeEdge( gn[5], gn[3] ) );
        allEdges.insert( makeEdge( gn[0], gn[3] ) );
        allEdges.insert( makeEdge( gn[1], gn[4] ) );
        allEdges.insert( makeEdge( gn[2], gn[5] ) );

        allFaces.insert( makeFace( {gn[0], gn[1], gn[2]} ) );
        allFaces.insert( makeFace( {gn[3], gn[4], gn[5]} ) );
        allFaces.insert( makeFace( {gn[0], gn[1], gn[4], gn[3]} ) );
        allFaces.insert( makeFace( {gn[1], gn[2], gn[5], gn[4]} ) );
        allFaces.insert( makeFace( {gn[2], gn[0], gn[3], gn[5]} ) );
      }
      else if( npe == 5 )  // ---- Pyramid ---------------------------------
      {
        allEdges.insert( makeEdge( gn[0], gn[1] ) );
        allEdges.insert( makeEdge( gn[1], gn[2] ) );
        allEdges.insert( makeEdge( gn[2], gn[3] ) );
        allEdges.insert( makeEdge( gn[3], gn[0] ) );
        allEdges.insert( makeEdge( gn[0], gn[4] ) );
        allEdges.insert( makeEdge( gn[1], gn[4] ) );
        allEdges.insert( makeEdge( gn[2], gn[4] ) );
        allEdges.insert( makeEdge( gn[3], gn[4] ) );

        allFaces.insert( makeFace( {gn[0], gn[1], gn[2], gn[3]} ) );
        allFaces.insert( makeFace( {gn[0], gn[1], gn[4]} ) );
        allFaces.insert( makeFace( {gn[1], gn[2], gn[4]} ) );
        allFaces.insert( makeFace( {gn[2], gn[3], gn[4]} ) );
        allFaces.insert( makeFace( {gn[3], gn[0], gn[4]} ) );
      }
    }
  } );

  // --- ownership: entity owned by this rank if the node with the smallest
  //     global ID in the entity is owned (ghostRank < 0) on this rank. ----

  auto nodeIsOwned = [&]( globalIndex gid ) -> bool
  {
    auto it = nodeManager.globalToLocalMap().find( gid );
    if( it == nodeManager.globalToLocalMap().end() )
      return false;
    return nodeGhostRank[ it->second ] < 0;
  };

  localIndex ownedV = 0;
  for( globalIndex const & gid : allNodes )
  {
    if( nodeIsOwned( gid ) )
      ++ownedV;
  }

  localIndex ownedE = 0;
  for( EdgeKey const & e : allEdges )
  {
    if( nodeIsOwned( e.first ) )   // e.first is the min global ID (sorted)
      ++ownedE;
  }

  localIndex ownedF = 0;
  for( FaceKey const & f : allFaces )
  {
    if( nodeIsOwned( f[0] ) )      // f[0] is the min global ID (sorted)
      ++ownedF;
  }

  // --- MPI reduce --------------------------------------------------------
  globalIndex const V = MpiWrapper::sum< globalIndex >( ownedV );
  globalIndex const E = MpiWrapper::sum< globalIndex >( ownedE );
  globalIndex const F = MpiWrapper::sum< globalIndex >( ownedF );
  globalIndex const C = MpiWrapper::sum< globalIndex >( numOwnedCells );

  GEOS_LOG_RANK_0( "computeEulerCharacteristic (bulk):"
                   << " V=" << V << " E=" << E << " F=" << F << " C=" << C
                   << " chi=" << (V - E + F - C) );

  return static_cast< integer >( V - E + F - C );
}

}
