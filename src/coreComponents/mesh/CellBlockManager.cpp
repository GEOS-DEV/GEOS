/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ElementManagerT.cpp
 */

#include "CellBlockManager.hpp"

#include "FaceManager.hpp"

namespace geosx
{
using namespace dataRepository;

CellBlockManager::CellBlockManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  this->registerGroup< Group >( keys::cellBlocks );
}

CellBlockManager::~CellBlockManager()
{
  // TODO Auto-generated destructor stub
}

void CellBlockManager::resize( integer_array const & numElements,
                               string_array const & regionNames,
                               string_array const & GEOSX_UNUSED_PARAM( elementTypes ) )
{
  localIndex const numRegions = LvArray::integerConversion< localIndex >( regionNames.size());
  for( localIndex reg=0; reg<numRegions; ++reg )
  {
    this->getRegion( regionNames[reg] ).resize( numElements[reg] );
  }
}

Group * CellBlockManager::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

std::map< localIndex, std::vector< localIndex > > CellBlockManager::getNodeToElem() const
{
  // TODO check that we can use the [0, nNodes[ size for first dimension (mainly parallel issues),
  // and use a std::vector< std::vector< localIndex > > instead;
  std::map< localIndex, std::vector< localIndex > > result;

  // This is a dummy non parallel implementation...
  for( localIndex iCellBlock = 0; iCellBlock < numCellBlocks(); ++iCellBlock)// Do not need index
  {
    const CellBlockABC & cb = this->getCellBlocks().getGroup< const CellBlockABC >( iCellBlock );
    CellBlockABC::NodeMapType const elemToNode = cb.getElemToNode();
    for( localIndex iElem = 0; iElem < cb.numElements(); ++iElem )
    {
      for( localIndex a = 0; a <  cb.numNodesPerElement(); ++a )
      {
        localIndex const nodeIndex = elemToNode( iElem, a );
        result[nodeIndex].push_back( iElem );
      }
    }
  }

  return result;
}

const Group & CellBlockManager::getCellBlocks() const
{
  return this->getGroup( keys::cellBlocks );
}

localIndex CellBlockManager::numNodes() const
{
  localIndex numNodes = 0;
  const Group & cellBlocks = this->getCellBlocks();
  cellBlocks.forSubGroups< const CellBlockABC >([&](const CellBlockABC & cb ){
    numNodes += cb.numNodesPerElement() * cb.numElements();
  } );
  return numNodes;
}

localIndex CellBlockManager::numCellBlocks() const
{
  return this->getCellBlocks().numSubGroups();
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlockManager, string const &, Group * const )
}
