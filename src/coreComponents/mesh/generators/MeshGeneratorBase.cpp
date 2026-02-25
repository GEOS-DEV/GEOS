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
  // Count bulk volumetric cells only (exclude fracture / surface sub-regions)
  localIndex numCells = 0;
  elemManager.forElementSubRegions< CellElementSubRegion >(
    [&]( CellElementSubRegion const & subRegion )
  {
    numCells += subRegion.size();
  } );

  localIndex const V = nodeManager.size();
  localIndex const E = edgeManager.size();
  localIndex const F = faceManager.size();
  localIndex const C = numCells;

  return static_cast< integer >( V - E + F - C );
}

}
