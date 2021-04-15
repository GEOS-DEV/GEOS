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
#include "NodeManager.hpp"

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

/**
 * @brief Convenience programming structure that holds the nodes for a given Face. No information about which cell though.
 *
 * This structure holds `<` and `==` such that equal instances in a sorted array
 * are meant to be consecutive. (see `std::unique` for example).
 * Two instances with same nodes in different order are considered equal.
 * Such a case often happen since the same surface is shared by two adjacent elements
 * in conformal meshes.
 */
struct NodesOfFace
{
  NodesOfFace( std::vector< localIndex > nodes_ ):
    nodes( nodes_ ),
    sortedNodes( nodes_ )
  {
    std::sort( sortedNodes.begin(), sortedNodes.end() );
  }

  /**
   * @brief Imposes an ordering on NodesOfFaces.
   * @param [in] rhs the NodesOfFace to compare against.
   * @return a boolean.
   */
  bool operator<( NodesOfFace const & rhs ) const
  {
    return sortedNodes < rhs.sortedNodes;
  }

  /**
   * @brief Two NodesOfFace instances are considered if they share the same node set, whatever the order.
   * @param [in] rhs the NodesOfFace to compare against.
   * @return a boolean.
   */
  bool operator==( NodesOfFace const & rhs ) const
  {
    return sortedNodes == rhs.sortedNodes;
  }

  /// The list of nodes describing the face.
  std::vector< localIndex > nodes;

private:
  /// Sorted nodes describing the face; mainly for comparison reasons.
  std::vector< localIndex > sortedNodes;
};

// FIXME unmodified, copied for convenience.
localIndex calculateTotalNumberOfFaces2( ArrayOfArraysView< NodesOfFace const > const & facesByLowestNode,
                                        arrayView1d< localIndex > const & uniqueFaceOffsets )
{
  localIndex const numNodes = facesByLowestNode.size();
  GEOSX_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );

  uniqueFaceOffsets[0] = 0;

  // Loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );

    // If there are no faces associated with this node we can skip it.
    if( numFaces == 0 )
      return;

    localIndex & numUniqueFaces = uniqueFaceOffsets[ nodeID + 1 ];
    numUniqueFaces = 0;

    // Otherwise since facesByLowestNode[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries. Since each face can appear at most twice if we find a match we
    // can skip the next entry as well.
    localIndex j = 0;
    for(; j < numFaces - 1; ++j )
    {
      ++numUniqueFaces;
      j += facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 );
    }

    numUniqueFaces += j == numFaces - 1;
  } );

  // At this point uniqueFaceOffsets[ i ] holds the number of unique face associated with node i - 1.
  // Perform an inplace prefix-sum to get the unique face offset.
  RAJA::inclusive_scan_inplace< parallelHostPolicy >( uniqueFaceOffsets.begin(), uniqueFaceOffsets.end() );
  return uniqueFaceOffsets.back();
}

/**
 * @brief Fills the appropriate subpart @p faceToNodes for face @p faceID and nodes contained in @p fb.
 * @param [in] faceID The face index.
 * @param [in] fb The nodes for @p faceID.
 * @param [out] faceToNodes No input data is used.
 */
void addFace( localIndex const faceID,
              NodesOfFace const & fb,
              ArrayOfArrays< localIndex > & faceToNodes )
{
  localIndex const numFaceNodes = fb.nodes.size();
  // FIXME The size should be OK because it's been allocated previously.
  for (localIndex i = 0; i < numFaceNodes; ++i ){
    faceToNodes[ faceID ][i] = fb.nodes[i];
  }
  GEOSX_ASSERT_EQ( numFaceNodes, faceToNodes.sizeOfArray( faceID ) );
  GEOSX_DEBUG_VAR( numFaceNodes );
}

/**
 * @brief Fills the face to nodes map.
 */
void populateFaceToNodesMap( ArrayOfArraysView< NodesOfFace const > const & facesByLowestNode,
                             arrayView1d< localIndex const > const & uniqueFaceOffsets,
                             ArrayOfArrays< localIndex > & faceToNodes )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  GEOSX_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );
  GEOSX_ERROR_IF_NE( numUniqueFaces, faceToNodes.size() );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );
    // TODO same loop concern than for resizeFaceToNodesMap below.
    // loop over all the FaceBuilders associated with the node
    localIndex j, curFaceID = uniqueFaceOffsets[ nodeID ];
    for( j = 0; j < numFaces - 1; ++j, ++curFaceID )
    {
      addFace( curFaceID, facesByLowestNode( nodeID, j ), faceToNodes );
      // If two subsequent FaceBuilders compare equal then they describe an interior face.
      if( facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 ) ) { ++j; }
      // Otherwise it's a boundary face.
    }

    if( j == numFaces - 1 )
    {
      addFace( curFaceID, facesByLowestNode( nodeID, j ), faceToNodes );
    }
  } );
}

// We can get rid of all the Region stuff if `NodesOfFace` holds its nodes.
void resizeFaceToNodesMap( ArrayOfArraysView< NodesOfFace const > const & facesByLowestNode,
                          arrayView1d< localIndex const > const & uniqueFaceOffsets,
                          ArrayOfArrays< localIndex > & faceToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  array1d< localIndex > numNodesPerFace( numUniqueFaces );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalFaceNodes( 0.0 );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex curFaceID = uniqueFaceOffsets[ nodeID ];
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );

    // If there are no faces associated with this node we can skip it.
    // TODO The loop should be bullet proof w.r.t. this kind of cases.
    if( numFaces == 0 )
      return;

    // TODO This kind of loops with a side effect after (`j == numFaces - 1`) can surely be improved.
    // loop over all the FaceBuilders associated with the node
    localIndex j = 0;
    for(; j < numFaces - 1; ++j, ++curFaceID )
    {
      const NodesOfFace & fb = facesByLowestNode( nodeID, j );
      numNodesPerFace[ curFaceID ] = fb.nodes.size();
      totalFaceNodes += numNodesPerFace[ curFaceID ];

      if( facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 ) ) { ++j; }
    }

    if( j == numFaces - 1 )
    {
      const NodesOfFace & fb = facesByLowestNode( nodeID, j );
      numNodesPerFace[ curFaceID ] = fb.nodes.size();
      totalFaceNodes += numNodesPerFace[ curFaceID ];
    }
  } );

  // Resize the face to node map.
  faceToNodeMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numUniqueFaces;
  faceToNodeMap.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalFaceNodes.get() + numUniqueFaces * FaceManager::nodeMapExtraSpacePerFace() * ( 1 + 2 * overAllocationFactor );
  faceToNodeMap.reserveValues( valuesToReserve );

  // Append the individual arrays.
  for( localIndex faceID = 0; faceID < numUniqueFaces; ++faceID )
  {
    faceToNodeMap.appendArray( numNodesPerFace[ faceID ] );
    faceToNodeMap.setCapacityOfArray( faceToNodeMap.size() - 1,
                                      numNodesPerFace[ faceID ] + FaceManager::nodeMapExtraSpacePerFace() );
  }
}

ArrayOfArrays< localIndex > CellBlockManager::getFaceToNodes( localIndex numNodes ) const
{
  const localIndex maxFacesPerNode = 200; // TODO deal with this differently
  ArrayOfArrays< NodesOfFace > facesByLowestNode( numNodes, 2 * maxFacesPerNode );
  // Begin of `createFacesByLowestNode`
  // The function is a bit simplified and is not run in parallel anymore.
  // Can be improved.
  const Group & cellBlocks = this->getCellBlocks();
  cellBlocks.forSubGroups< const CellBlockABC >( [&]( const CellBlockABC & cb ) {
    localIndex const numFacesPerElement = cb.numFacesPerElement();
    localIndex const numElements = cb.numElements();

    for( localIndex iElement = 0; iElement < numElements; ++iElement )
    {
      // Looping on the faces of the cell
      for( localIndex iFace = 0; iFace < numFacesPerElement; ++iFace )
      {
        // Get all the nodes of the cell
        std::vector< localIndex > nodesInFace = cb.getFaceNodes( iElement, iFace );
        localIndex const & lowestNode = *std::min_element( nodesInFace.cbegin(), nodesInFace.cend() );
        facesByLowestNode.emplaceBack( lowestNode, nodesInFace );
      }
    }
  } );

  // Loop over all the nodes and sort the associated faces.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    NodesOfFace * const faces = facesByLowestNode[ nodeID ];
    std::sort( faces, faces + facesByLowestNode.sizeOfArray( nodeID ) );
  } );
  // End of `createFacesByLowestNode`

  array1d< localIndex > uniqueFaceOffsets( numNodes + 1 );
  // FIXME We do not use the `numFaces` returned by `calculateTotalNumberOfFaces2` anymore...
  calculateTotalNumberOfFaces2( facesByLowestNode.toViewConst(), uniqueFaceOffsets );

  ArrayOfArrays< localIndex > faceToNodeMap;

  resizeFaceToNodesMap( facesByLowestNode.toViewConst(),
                       uniqueFaceOffsets,
                       faceToNodeMap );

  populateFaceToNodesMap( facesByLowestNode.toViewConst(),
                          uniqueFaceOffsets,
                          faceToNodeMap );

  return faceToNodeMap;
}

std::map< localIndex, std::vector< localIndex > > CellBlockManager::getFaceToElem() const
{
  std::map< localIndex, std::vector< localIndex > > result;
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

localIndex CellBlockManager::numUniqueNodes() const
{
  std::set< localIndex > uniqueNodes;
  const Group & cellBlocks = this->getCellBlocks();
  cellBlocks.forSubGroups< const CellBlockABC >([&](const CellBlockABC & cb ){
    CellBlockABC::NodeMapType const elemToNode = cb.getElemToNode();
    for( localIndex iElem = 0; iElem < cb.numElements(); ++iElem )
    {
      for( localIndex a = 0; a <  cb.numNodesPerElement(); ++a )
      {
        uniqueNodes.insert(elemToNode( iElem, a ));
      }
    }
  } );
  return uniqueNodes.size();
}

localIndex CellBlockManager::numCellBlocks() const
{
  return this->getCellBlocks().numSubGroups();
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlockManager, string const &, Group * const )
}
