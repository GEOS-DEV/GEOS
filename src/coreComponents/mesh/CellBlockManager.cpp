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

std::map< localIndex, std::vector< localIndex > > CellBlockManager::getNodeToElements() const
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
struct NodesAndElementOfFace
{
  NodesAndElementOfFace( std::vector< localIndex > nodes_, localIndex element_, localIndex iCellBlock_, localIndex iFace_ ):
    nodes( nodes_ ),
    element( element_ ),
    iCellBlock( iCellBlock_ ),
    iFace( iFace_ ),
    sortedNodes( nodes_ )
  {
    std::sort( sortedNodes.begin(), sortedNodes.end() );
  }

  /**
   * @brief Imposes an ordering on NodesAndElementOfFace.
   * @param [in] rhs the NodesAndElementOfFace to compare against.
   * @return a boolean.
   */
  bool operator<( NodesAndElementOfFace const & rhs ) const
  {
    return sortedNodes < rhs.sortedNodes;
  }

  /**
   * @brief Two NodesAndElementOfFace instances are considered if they share the same node set, whatever the order.
   * @param [in] rhs the NodesAndElementOfFace to compare against.
   * @return a boolean.
   */
  bool operator==( NodesAndElementOfFace const & rhs ) const
  {
    return sortedNodes == rhs.sortedNodes;
  }

  /// The list of nodes describing the face.
  std::vector< localIndex > nodes;

  /**
   * @brief The element to which this face belongs.
   *
   * Each face may belong to multiple elements.
   * But during the identification process (we loop on each face of each element),
   * we store the the bovious cell we are iterating on.
   * The we'll be able to identify the duplicated faces because we also have the nodes.
   */
  localIndex element;
  localIndex iCellBlock;
  localIndex iFace;

private:
  /// Sorted nodes describing the face; mainly for comparison reasons.
  std::vector< localIndex > sortedNodes;
};

// FIXME unmodified, copied for convenience.
localIndex calculateTotalNumberOfFaces2( ArrayOfArraysView< NodesAndElementOfFace const > const & facesByLowestNode,
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
 * @brief Fills the appropriate subpart @p faceToNodes for face @p faceID and nodes contained in @p nodesAndElementOfFace.
 * @param [in] faceID The face index.
 * @param [in] nodesAndElementOfFace The nodes and element for @p faceID.
 * @param [out] faceToNodes No input data is used.
 */
void insertFaceToNodesEntry( localIndex const faceID,
                             NodesAndElementOfFace const & nodesAndElementOfFace,
                             ArrayOfArrays< localIndex > & faceToNodes )
{
  localIndex const numFaceNodes = nodesAndElementOfFace.nodes.size();
  // FIXME The size should be OK because it's been allocated previously.
  for( localIndex i = 0; i < numFaceNodes; ++i )
  {
    faceToNodes[faceID][i] = nodesAndElementOfFace.nodes[i];
  }
  GEOSX_ASSERT_EQ( numFaceNodes, faceToNodes.sizeOfArray( faceID ) );
  GEOSX_DEBUG_VAR( numFaceNodes );
}

/**
 * @brief Fills the face to nodes map and element to node maps
 */
void populateMaps( ArrayOfArraysView< NodesAndElementOfFace const > const & facesByLowestNode,
                   arrayView1d< localIndex const > const & uniqueFaceOffsets,
                   ArrayOfArrays< localIndex > & faceToNodes,
                   arrayView2d< localIndex > const & faceToElements )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  GEOSX_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );
  GEOSX_ERROR_IF_NE( numUniqueFaces, faceToNodes.size() );
  GEOSX_ERROR_IF_NE( numUniqueFaces, faceToElements.size( 0 ) );
  GEOSX_ERROR_IF_NE( 2, faceToElements.size( 1 ) );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );
    // TODO same loop concern than for resizeFaceToNodesMap below.
    // loop over all the `NodesAndElementOfFace` associated with the node
    localIndex j, curFaceID = uniqueFaceOffsets[nodeID];
    for( j = 0; j < numFaces - 1; ++j, ++curFaceID )
    {
      insertFaceToNodesEntry( curFaceID, facesByLowestNode( nodeID, j ), faceToNodes );
      // If two subsequent `NodesAndElementOfFace` compare equal then they describe an interior face.
      // It must therefore be considered "twice".
      // Otherwise it's a boundary face, and "once" is enough.
      NodesAndElementOfFace const & f0 = facesByLowestNode( nodeID, j );
      NodesAndElementOfFace const & f1 = facesByLowestNode( nodeID, j + 1 );
      faceToElements( curFaceID, 0 ) = f0.element;
      if( f0 == f1 )
      {
        faceToElements( curFaceID, 1 ) = f1.element;
        ++j;
      }
      else
      {
        faceToElements( curFaceID, 1 ) = -1;
      }
    }

    if( j == numFaces - 1 )
    {
      insertFaceToNodesEntry( curFaceID, facesByLowestNode( nodeID, j ), faceToNodes );
      faceToElements( curFaceID, 0 ) = facesByLowestNode( nodeID, j ).element;
      faceToElements( curFaceID, 1 ) = -1;
    }
  } );
}

/**
 * @brief Resize the Maps
 * @param [in] facesByLowestNode and array of size numNodes of arrays of FaceBuilders associated with each node.
 * @param [in] uniqueFaceOffsets an containing the unique face IDs for each node in facesByLowestNode.
 * @param [out] faceToNodeMap the map from faces to nodes. This function resizes the array appropriately.
 * @param [out] faceToElemMap the map from faces to elements. This function resizes the array appropriately.
 */
void resizeMaps( ArrayOfArraysView< NodesAndElementOfFace const > const & facesByLowestNode,
                arrayView1d< localIndex const > const & uniqueFaceOffsets,
                ArrayOfArrays< localIndex > & faceToNodeMap,
                array2d< localIndex > & faceToElemMap )
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
    localIndex j;
    for( j = 0; j < numFaces - 1; ++j, ++curFaceID )
    {
      const NodesAndElementOfFace & f0 = facesByLowestNode( nodeID, j );
      const NodesAndElementOfFace & f1 = facesByLowestNode( nodeID, j + 1 );

      numNodesPerFace[curFaceID] = f0.nodes.size();
      totalFaceNodes += numNodesPerFace[curFaceID];

      if( f0 == f1 )
      { ++j; }
    }

    if( j == numFaces - 1 )
    {
      const NodesAndElementOfFace & f = facesByLowestNode( nodeID, j );
      numNodesPerFace[curFaceID] = f.nodes.size();
      totalFaceNodes += numNodesPerFace[curFaceID];
    }
  } );

  // Resize the face to node map.
  faceToNodeMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numUniqueFaces; //TODO why this allocation factor
  faceToNodeMap.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalFaceNodes.get() + numUniqueFaces * FaceManager::nodeMapExtraSpacePerFace() * ( 1 + 2 * overAllocationFactor );
  faceToNodeMap.reserveValues( valuesToReserve );
  faceToNodeMap.reserveValues( 2 * entriesToReserve ); //TODO I don"t undertand anythin about LVARRAY :@

  // Append the individual arrays.
  for( localIndex faceID = 0; faceID < numUniqueFaces; ++faceID )
  {
    faceToNodeMap.appendArray( numNodesPerFace[ faceID ] );
    faceToNodeMap.setCapacityOfArray( faceToNodeMap.size() - 1,
                                      numNodesPerFace[ faceID ] + FaceManager::nodeMapExtraSpacePerFace() );
                                      //TODO THIS LINE IS DAMN STRANGE
  }

  // Each face may belong to _maximum_ 2 elements.
  // If it belongs to only one, we put `-1` for the undefined value.
  faceToElemMap.resize( numUniqueFaces, 2 );
}

void CellBlockManager::buildMaps( localIndex numNodes )
{
  const localIndex maxFacesPerNode = 200; // TODO deal with this differently
  ArrayOfArrays< NodesAndElementOfFace > facesByLowestNode( numNodes, 2 * maxFacesPerNode );
  // Begin of `createFacesByLowestNode`
  // The function is a bit simplified and is not run in parallel anymore.
  // Can be improved.
  const Group & cellBlocks = this->getCellBlocks();
  for( localIndex iCellBlock = 0; iCellBlock < cellBlocks.numSubGroups(); ++iCellBlock )
  {
    const CellBlockABC & cb = cellBlocks.getGroup< CellBlockABC >( iCellBlock );
    localIndex const numFacesPerElement = cb.numFacesPerElement();
    localIndex const numElements = cb.numElements();

    for( localIndex iElement = 0; iElement < numElements; ++iElement )
    {
      // Looping on the faces of the cell
      for( localIndex iFace = 0; iFace < numFacesPerElement; ++iFace )
      {
        // Get all the nodes of the cell
        std::vector< localIndex > const nodesInFace = cb.getFaceNodes( iElement, iFace );
        localIndex const & lowestNode = *std::min_element( nodesInFace.cbegin(), nodesInFace.cend() );
        facesByLowestNode.emplaceBack( lowestNode, nodesInFace, iElement, iCellBlock, iFace );
      }
    }
  }

  // Loop over all the nodes and sort the associated faces.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    NodesAndElementOfFace * const faces = facesByLowestNode[ nodeID ];
    std::sort( faces, faces + facesByLowestNode.sizeOfArray( nodeID ) );
  } );
  // End of `createFacesByLowestNode`

  array1d< localIndex > uniqueFaceOffsets( numNodes + 1 );
  // FIXME We do not use the `numFaces` returned by `calculateTotalNumberOfFaces2` anymore...
  calculateTotalNumberOfFaces2( facesByLowestNode.toViewConst(), uniqueFaceOffsets );


  resizeMaps( facesByLowestNode.toViewConst(),
              uniqueFaceOffsets,
              m_faceToNodes,
              m_faceToElements );

  populateMaps( facesByLowestNode.toViewConst(),
                uniqueFaceOffsets,
                m_faceToNodes,
                m_faceToElements );

  // Filling the elements to faces maps in the cell blocks.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );
    localIndex curFaceID = uniqueFaceOffsets[nodeID];
    localIndex j;
    for( j = 0; j < numFaces - 1; ++j, ++curFaceID )
    {
      const NodesAndElementOfFace & f0 = facesByLowestNode( nodeID, j );
      CellBlock & cb0 = this->getGroup( keys::cellBlocks ).getGroup< CellBlock >( f0.iCellBlock );
      cb0.setElementToFaces( f0.element, f0.iFace, curFaceID );
      const NodesAndElementOfFace & f1 = facesByLowestNode( nodeID, j + 1 );
      CellBlock & cb1 = this->getGroup( keys::cellBlocks ).getGroup< CellBlock >( f1.iCellBlock );
      cb1.setElementToFaces( f1.element, f1.iFace, curFaceID );

      if( f0 == f1 ) { ++j; }
    }
    if( j == numFaces - 1 )
    {
      const NodesAndElementOfFace & f0 = facesByLowestNode( nodeID, j );
      CellBlock & cb0 = this->getGroup( keys::cellBlocks ).getGroup< CellBlock >( f0.iCellBlock );
      cb0.setElementToFaces( f0.element, f0.iFace, curFaceID );
    }
  }
}

//TODO return views
ArrayOfArrays< localIndex > CellBlockManager::getFaceToNodes() const
{
  return m_faceToNodes;
}

// TODO return views
array2d< localIndex > CellBlockManager::getFaceToElements() const
{
  return m_faceToElements;
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
