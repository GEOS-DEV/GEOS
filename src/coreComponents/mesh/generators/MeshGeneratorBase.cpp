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

#include "MeshGeneratorBase.hpp"
#include "InternalWellGenerator.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/generators/CellBlockManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshGeneratorBase::MeshGeneratorBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

Group * MeshGeneratorBase::createChild( string const & childKey, string const & childName )
{
  GEOSX_LOG_RANK_0( "Adding Mesh attribute: " << childKey << ", " << childName );
  std::unique_ptr< InternalWellGenerator > wellGen = InternalWellGenerator::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< InternalWellGenerator >( childName, std::move( wellGen ) );
}

MeshGeneratorBase::CatalogInterface::CatalogType & MeshGeneratorBase::getCatalog()
{
  static MeshGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

CellBlockManagerABC & MeshGeneratorBase::generateMesh( Group & parent )
{
  CellBlockManager & cellBlockManager = parent.registerGroup< CellBlockManager >( keys::cellManager );

  generateCellBlockManager( cellBlockManager );

  this->generateWells( cellBlockManager );

  return cellBlockManager;
}

void MeshGeneratorBase::generateWells( CellBlockManager & cellBlockManager )
{
  forSubGroups< InternalWellGenerator >( [&]( InternalWellGenerator & wellGen ) {
    LineBlock & wb = cellBlockManager.registerLineBlock( wellGen.getName() );
    wb.setNumElements( wellGen.getNumElements() );
    wb.setElemCoords( wellGen.getElemCoords() );
    wb.setNextElemIndex( wellGen.getNextElemIndex() );
    wb.setPrevElemIndices( wellGen.getPrevElemIndices() );
    wb.setElemToNodesMap( wellGen.getElemToNodesMap() );
    wb.setElemVolume( wellGen.getElemVolume() );
    wb.setElementRadius( wellGen.getElementRadius() );
    wb.setNumNodes( wellGen.getNumNodes() );
    wb.setNodeCoords( wellGen.getNodeCoords() );
    wb.setNumPerforations( wellGen.getNumPerforations() );
    wb.setPerfCoords( wellGen.getPerfCoords() );
    wb.setPerfTransmissibility( wellGen.getPerfTransmissibility() );
    wb.setPerfElemIndex( wellGen.getPerfElemIndex() );
  } );
}
}
