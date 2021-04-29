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

#include "CellBlockManager.hpp"
#include "NodeManager.hpp"

#include "FaceManager.hpp"

namespace geosx
{
using namespace dataRepository;

string const CellBlockManager::cellBlocksKey = "cellBlocks";

const int maxEdgesPerNode = 200;// TODO deal with this.

CellBlockManager::CellBlockManager( string const & name, Group * const parent ):
  ObjectManagerBase( name, parent )
{
  this->registerGroup< Group >( cellBlocksKey );
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

/**
 * @class EdgeBuilder
 * @brief This class stores the data necessary to construct the various edge maps.
 */
struct EdgeBuilder
{

  /**
   * @brief Constructor.
   * @param [in] n1_ the greater of the two node indices that comprise the edge.
   * @param [in] faceID_ the ID of the face this edge came from.
   * @param [in] faceLocalEdgeIndex_ the face local index of this edge.
   */
  EdgeBuilder( localIndex const n1_,
               localIndex const faceID_,
               localIndex const faceLocalEdgeIndex_ ):
    n1( int32_t( n1_ ) ),
    faceID( int32_t( faceID_ ) ),
    faceLocalEdgeIndex( int32_t( faceLocalEdgeIndex_ ) )
  {}

  /**
   * @brief Imposes an ordering on EdgeBuilders. First compares n1 and then the faceID.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator<( EdgeBuilder const & rhs ) const
  {
    if( n1 < rhs.n1 ) return true;
    if( n1 > rhs.n1 ) return false;
    return faceID < rhs.faceID;
  }

  /**
   * @brief Return true if the two EdgeBuilders share the same greatest node index.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator==( EdgeBuilder const & rhs ) const
  { return n1 == rhs.n1; }

  /**
   * @brief Return true if the two EdgeBuilders don't share the same greatest node index.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator!=( EdgeBuilder const & rhs ) const
  { return n1 != rhs.n1; }

  int32_t n1;                  // The larger of the two node indices that comprise the edge.
  int32_t faceID;              // The face the edge came from.
  int32_t faceLocalEdgeIndex;  // The face local index of the edge.
};

/**
 * @brief Return the total number of unique faces and fill in the uniqueFaceOffsets array.
 * @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
 * @param [out] uniqueFaceOffsets an array of size numNodes + 1. After this function returns node i contains
 *              faces with IDs ranging from uniqueFaceOffsets[ i ] to uniqueFaceOffsets[ i + 1 ] - 1.
 * @return return total number of faces
 */
localIndex calculateTotalNumberOfFaces( ArrayOfArraysView< NodesAndElementOfFace const > const & lowestNodeToFaces,
                                        arrayView1d< localIndex > const & uniqueFaceOffsets )
{
  localIndex const numNodes = lowestNodeToFaces.size();
  GEOSX_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );
  uniqueFaceOffsets.setValues< serialPolicy >( 0. );

  // Loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
    // Since lowestNodeToFaces[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries. Since each face can appear at most twice if we find a match we
    // can skip the next entry as well.
    localIndex & numUniqueFaces = uniqueFaceOffsets[ nodeID + 1 ];
    for( localIndex j = 0; j < numFaces; ++j )
    {
      ++numUniqueFaces;
      if ( j < numFaces - 1 )
      {
        if( lowestNodeToFaces( nodeID, j ) == lowestNodeToFaces( nodeID, j + 1 ) )
        {
          ++j;
        }
      }
    }
  } );
  // At this point uniqueFaceOffsets[ i ] holds the number of unique face associated with node i - 1.
  // Perform an inplace prefix-sum to get the unique face offset.
  RAJA::inclusive_scan_inplace< parallelHostPolicy >( uniqueFaceOffsets.begin(), uniqueFaceOffsets.end() );
  return uniqueFaceOffsets.back();
}

/**
 * @brief Return the total number of unique edges and fill in the uniqueEdgeOffsets array.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [out] uniqueEdgeOffsets an array of size numNodes + 1. After this function returns node i contains
 * edges with IDs ranging from uniqueEdgeOffsets[ i ] to uniqueEdgeOffsets[ i + 1 ] - 1.
 */
localIndex calculateTotalNumberOfEdges( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                                        arrayView1d< localIndex > const & uniqueEdgeOffsets )
{
  localIndex const numNodes = edgesByLowestNode.size();
  GEOSX_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );

  uniqueEdgeOffsets[0] = 0;

  // Loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // If there are no edges associated with this node we can skip it.
    if( numEdges == 0 )
      return;

    localIndex & numUniqueEdges = uniqueEdgeOffsets[ nodeID + 1 ];
    numUniqueEdges = 0;

    // Otherwise since edgesByLowestNode[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries.
    localIndex j = 0;
    for(; j < numEdges - 1; ++j )
    {
      numUniqueEdges += edgesByLowestNode( nodeID, j ) != edgesByLowestNode( nodeID, j + 1 );
    }

    numUniqueEdges += j == numEdges - 1;
  } );

  // At this point uniqueEdgeOffsets[ i ] holds the number of unique edges associated with node i - 1.
  // Perform an inplace prefix-sum to get the unique edge offset.
  RAJA::inclusive_scan_inplace< parallelHostPolicy >( uniqueEdgeOffsets.begin(), uniqueEdgeOffsets.end() );

  return uniqueEdgeOffsets.back();
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
 * @brief Add an edge to the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in/out] faceToEdgeMap the map from face IDs to edge IDs.
 * @param [in/out] edgeToFacemap the map from edgeIDs to faceIDs.
 * @param [in/out] edgeToNodeMap the map from edgeIDs to nodeIDs.
 * @param [in] edgeID the ID of the edge to add.
 * @param [in] firstNodeID the ID of the first node of the edge.
 * @param [in] firstMatch the index of the first EdgeBuilder that describes this edge in edgesByLowestNode[ firstNodeID ].
 * @param [in] numMatches the number of EdgeBuilders that describe this edge in edgesByLowestNode[ firstNodeID ].
 */
void addEdge( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
              ArrayOfArraysView< localIndex > const & faceToEdgeMap,
              ArrayOfSetsView< localIndex > const & edgeToFaceMap,
              arrayView2d< localIndex > const & edgeToNodeMap,
              localIndex const edgeID,
              localIndex const firstNodeID,
              localIndex const firstMatch,
              localIndex const numMatches )
{
  GEOSX_ASSERT_GE( edgeToFaceMap.capacityOfSet( edgeID ), numMatches );

  // Populate the edge to node map.
  edgeToNodeMap( edgeID, 0 ) = firstNodeID;
  edgeToNodeMap( edgeID, 1 ) = edgesByLowestNode( firstNodeID, firstMatch ).n1;

  // Loop through all the matches and fill in the face to edge and edge to face maps.
  for( localIndex i = 0; i < numMatches; ++i )
  {
    localIndex const faceID = edgesByLowestNode( firstNodeID, firstMatch + i ).faceID;
    localIndex const faceLocalEdgeIndex = edgesByLowestNode( firstNodeID, firstMatch + i ).faceLocalEdgeIndex;

    faceToEdgeMap( faceID, faceLocalEdgeIndex ) = edgeID;
    edgeToFaceMap.insertIntoSet( edgeID, faceID );
  }
}

/**
 * @brief Fills the face to nodes map and face to element maps
 * @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
 * @param [in] uniqueFaceOffsets an array containing the unique ID of the first face associated with each node.
 * @param [inout] faceToElements the face to element map.
 * @param [inout] faceToNodes the face to node map.
 */
void populateFaceMaps( ArrayOfArraysView< NodesAndElementOfFace const > const & lowestNodeToFaces,
                       arrayView1d< localIndex const > const & uniqueFaceOffsets,
                       ArrayOfArrays< localIndex > & faceToNodes,
                       arrayView2d< localIndex > & faceToElements )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = lowestNodeToFaces.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  GEOSX_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );
  GEOSX_ERROR_IF_NE( numUniqueFaces, faceToNodes.size() );
  GEOSX_ERROR_IF_NE( numUniqueFaces, faceToElements.size( 0 ) );
  GEOSX_ERROR_IF_NE( 2, faceToElements.size( 1 ) );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
    // loop over all the `NodesAndElementOfFace` associated with the node.
    for( localIndex j = 0, curFaceID = uniqueFaceOffsets[nodeID]; j < numFaces ; ++j, ++curFaceID )
    {
      // If two subsequent `NodesAndElementOfFace` compare equal then they describe an interior face.
      // It must therefore be considered "twice".
      // Otherwise it's a boundary face, and "once" is enough.
      NodesAndElementOfFace const & f0 = lowestNodeToFaces( nodeID, j );
      insertFaceToNodesEntry( curFaceID, f0, faceToNodes );
      faceToElements( curFaceID, 0 ) = f0.element;
      faceToElements( curFaceID, 1 ) = -1; // TODO Make a constant

      // This is where we check for the two subsequent faces (when they exist).
      if( j < numFaces - 1 )
      {
        NodesAndElementOfFace const & f1 = lowestNodeToFaces( nodeID, j + 1 );
        if( f0 == f1 )
        {
          faceToElements( curFaceID, 1 ) = f1.element;
          ++j;
        }
      }
    }
  } );
}

/**
 * @brief Populate the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an array containing the unique ID of the first edge associated with each node.
 * @param [in] faceToNodeMap the map from faces to nodes.
 * @param [in/out] faceToEdgeMap the map from face IDs to edge IDs.
 * @param [in/out] edgeToFacemap the map from edgeIDs to faceIDs.
 * @param [in/out] edgeToNodeMap the map from edgeIDs to nodeIDs.
 */
void populateEdgeMaps( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                       arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                       ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                       ArrayOfArrays< localIndex > & faceToEdgeMap,
                       ArrayOfSets< localIndex > & edgeToFaceMap,
                       arrayView2d< localIndex > const & edgeToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  const localIndex edgeMapExtraSpacePerFace = 4; // TODO

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numFaces = faceToNodeMap.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  GEOSX_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );
  GEOSX_ERROR_IF_NE( numFaces, faceToEdgeMap.size() );
  GEOSX_ERROR_IF_NE( numUniqueEdges, edgeToFaceMap.size() );
  GEOSX_ERROR_IF_NE( numUniqueEdges, edgeToNodeMap.size( 0 ) );

  // The face to edge map has the same shape as the face to node map, so we can resize appropriately.
  localIndex totalSize = 0;
  for( localIndex faceID = 0; faceID < numFaces; ++faceID )
  {
    totalSize += faceToNodeMap.sizeOfArray( faceID );
  }

  // Resize the face to edge map
  faceToEdgeMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numFaces;
  faceToEdgeMap.reserve( entriesToReserve );

  // Reserve space for the total number of face edges + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalSize + numFaces * edgeMapExtraSpacePerFace * ( 1 + 2 * overAllocationFactor );
  faceToEdgeMap.reserveValues( valuesToReserve );
  for( localIndex faceID = 0; faceID < numFaces; ++faceID )
  {
    faceToEdgeMap.appendArray( faceToNodeMap.sizeOfArray( faceID ) );
    faceToEdgeMap.setCapacityOfArray( faceToEdgeMap.size() - 1,
                                      faceToNodeMap.sizeOfArray( faceID ) + edgeMapExtraSpacePerFace );
  }

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeID ];
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // loop over all the EdgeBuilders associated with the node
    localIndex j = 0;
    while( j < numEdges - 1 )
    {
      // Find the number of EdgeBuilders that describe the same edge
      localIndex numMatches = 1;
      while( edgesByLowestNode( nodeID, j ) == edgesByLowestNode( nodeID, j + numMatches ) )
      {
        ++numMatches;
        if( j + numMatches == numEdges )
          break;
      }
      // Then add the edge.
      addEdge( edgesByLowestNode, faceToEdgeMap.toView(), edgeToFaceMap.toView(), edgeToNodeMap, curEdgeID, nodeID, j, numMatches );
      ++curEdgeID;
      j += numMatches;
    }

    if( j == numEdges - 1 )
    {
      addEdge( edgesByLowestNode, faceToEdgeMap.toView(), edgeToFaceMap.toView(), edgeToNodeMap, curEdgeID, nodeID, j, 1 );
    }
  } );
}


/**
 * @brief Resize the face maps
 * @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
 * @param [in] uniqueFaceOffsets an containing the unique face IDs for each node in lowestNodeToFaces.
 * @param [out] faceToNodeMap the map from faces to nodes. This function resizes the array appropriately.
 * @param [out] faceToElemMap the map from faces to elements. This function resizes the array appropriately.
 */
void resizeFaceMaps( ArrayOfArraysView< NodesAndElementOfFace const > const & lowestNodeToFaces,
                     arrayView1d< localIndex const > const & uniqueFaceOffsets,
                     ArrayOfArrays< localIndex > & faceToNodeMap,
                     ArrayOfArrays< localIndex > & faceToEdgesMap,
                     array2d< localIndex > & faceToElemMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = lowestNodeToFaces.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  array1d< localIndex > numNodesPerFace( numUniqueFaces );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalFaceNodes( 0.0 );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
    // loop over all the NodesAndElementOfFace associated with the node
    for( localIndex j = 0, curFaceID = uniqueFaceOffsets[ nodeID ]; j < numFaces; ++j, ++curFaceID )
    {
      const NodesAndElementOfFace & f0 = lowestNodeToFaces( nodeID, j );
      numNodesPerFace[curFaceID] = f0.nodes.size();
      totalFaceNodes += numNodesPerFace[curFaceID];

      if( ( j < numFaces - 1 ) and ( f0 == lowestNodeToFaces( nodeID, j + 1 ) ) )
      { ++j; }
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

  // TODO I'm really not sure.
  const localIndex edgeMapExtraSpacePerFace = 4;// TODO
  faceToEdgesMap.resize( numUniqueFaces, 2 * edgeMapExtraSpacePerFace );
}

/**
 * @brief Resize the edge to face map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an containing the unique edge IDs for each node in edgesByLowestNode.
 * param [out] edgeToFaceMap the map from edges to faces. This function resizes the array appropriately.
 */
void resizeEdgeToFaceMap( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                          arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                          ArrayOfSets< localIndex > & edgeToFaceMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  array1d< localIndex > numFacesPerEdge( numUniqueEdges );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalEdgeFaces( 0.0 );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeID ];
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // loop over all the EdgeBuilders associated with the node
    localIndex j = 0;
    while( j < numEdges - 1 )
    {
      // Find the number of EdgeBuilders that describe the same edge
      localIndex numMatches = 1;
      while( edgesByLowestNode( nodeID, j ) == edgesByLowestNode( nodeID, j + numMatches ) )
      {
        ++numMatches;
        if( j + numMatches == numEdges )
          break;
      }

      // The number of matches is the number of faces associated with this edge.
      numFacesPerEdge( curEdgeID ) = numMatches;
      totalEdgeFaces += numFacesPerEdge( curEdgeID );
      ++curEdgeID;
      j += numMatches;
    }

    if( j == numEdges - 1 )
    {
      numFacesPerEdge( curEdgeID ) = 1;
      totalEdgeFaces += numFacesPerEdge( curEdgeID );
    }
  } );

  // Resize the edge to face map
  edgeToFaceMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numUniqueEdges;
  edgeToFaceMap.reserve( entriesToReserve );

  const localIndex faceMapExtraSpacePerEdge = 4;// TODO
  // Reserve space for the total number of edge faces + extra space for existing edges + even more space for new edges.
  localIndex const valuesToReserve = totalEdgeFaces.get() + numUniqueEdges * faceMapExtraSpacePerEdge * ( 1 + 2 * overAllocationFactor );
  edgeToFaceMap.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex faceID = 0; faceID < numUniqueEdges; ++faceID )
  {
    edgeToFaceMap.appendSet( numFacesPerEdge[ faceID ] + faceMapExtraSpacePerEdge );
  }
}

/**
 * @brief Resize the edge to face map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an containing the unique edge IDs for each node in edgesByLowestNode.
 * param [out] edgeToFaceMap the map from edges to faces. This function resizes the array appropriately.
 */
void resizeEdgeMaps( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                           arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                           array2d< localIndex > & edgeToNodesMap,
                           ArrayOfSets< localIndex > & edgeToFaceMap )
{
  GEOSX_MARK_FUNCTION;

  resizeEdgeToFaceMap( edgesByLowestNode, uniqueEdgeOffsets, edgeToFaceMap );

  // Each face may belong to _maximum_ 2 elements.
  // If it belongs to only one, we put `-1` for the undefined value.
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  edgeToNodesMap.resize( numUniqueEdges, 2 );

}

/**
 * @brief Populate the lowestNodeToFaces map.
 * @param [in] numNodes Number of nodes
 * @param [in] cellBlocks The cell blocks on which we need to operate.
 *
 * For each face of each element of each cell blocks,
 * this function stores some information of the faces (@see NodesAndElementOfFace)
 * in a node to faces (information) mapping.
 * The key of this mapping is the lowest node index of the face.
 * E.g. faces {3, 5, 6, 2} and {4, 2, 9, 7} will both be stored in "bucket" of node 2.
 * Also, bucket of faces information are sorted (@see NodesAndElementOfFace) to make specific computations possible.
 */
ArrayOfArrays< NodesAndElementOfFace > createLowestNodeToFaces( localIndex numNodes, const Group & cellBlocks )
{
  const localIndex maxFacesPerNode = 200; // TODO deal with this differently
  ArrayOfArrays< NodesAndElementOfFace > lowestNodeToFaces( numNodes, 2 * maxFacesPerNode );

  // The function is a bit simplified and is not run in parallel anymore.
  // Can be improved.
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
        lowestNodeToFaces.emplaceBack( lowestNode, nodesInFace, iElement, iCellBlock, iFace );
      }
    }
  }

  // Loop over all the nodes and sort the associated faces.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    NodesAndElementOfFace * const faces = lowestNodeToFaces[ nodeID ];
    std::sort( faces, faces + lowestNodeToFaces.sizeOfArray( nodeID ) );
  } );

  return lowestNodeToFaces;
}

/**
 * @brief Populate the edgesByLowestNode map.
 * @param [in] numNodes The number of nodes.
 * @param [in] faceToNodeMap a map that associates an ordered list of nodes with each face.
 * @param [in/out] edgesByLowestNode of size numNodes, where each sub array has been preallocated to hold
 *        *enough* space.
 * For each edge of each face, this function gets the lowest node in the edge n0, creates an EdgeBuilder
 * associated with the edge and then appends the EdgeBuilder to edgesByLowestNode[ n0 ]. Finally it sorts
 * the contents of each sub-array of edgesByLowestNode from least to greatest.
 */
ArrayOfArrays< EdgeBuilder > createEdgesByLowestNode( localIndex numNodes,
                                                      ArrayOfArraysView< localIndex const > const & faceToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< EdgeBuilder > edgesByLowestNode( numNodes, 2 * maxEdgesPerNode );

//  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numFaces = faceToNodeMap.size();

  // loop over all the faces.
  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const faceID )
  {
    localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceID );

    // loop over all the nodes in the face. there will be an edge for each node.
    for( localIndex a=0; a< numNodesInFace; ++a )
    {
      // sort the nodes in order of index value.
      localIndex node0 = faceToNodeMap( faceID, a );
      localIndex node1 = faceToNodeMap( faceID, ( a + 1 ) % numNodesInFace );
      if( node0 > node1 )
        std::swap( node0, node1 );

      // And append the edge to edgesByLowestNode.
      edgesByLowestNode.emplaceBackAtomic< parallelHostAtomic >( node0, node1, faceID, a );
    }
  } );

  // Loop over all the nodes and sort the associated edges.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    EdgeBuilder * const edges = edgesByLowestNode[ nodeID ];
    std::sort( edges, edges + edgesByLowestNode.sizeOfArray( nodeID ) );
  } );

  return edgesByLowestNode;
}

/**
 * @brief Filling the elements to faces maps in the cell blocks.
 * @param lowestNodeToFaces The lowest node to faces information array.
 * @param uniqueFaceOffsets The unique face offsets.
 * @param cellBlocks The cell blocks for which we need to compute the element to faces mappings.
 *
 * @note @p lowestNodeToFaces and @p uniqueFaceOffsets are better described in the documentations of the functions that build them.
 */
void fillElementToFacesOfCellBlocks( ArrayOfArrays< NodesAndElementOfFace > const & lowestNodeToFaces,
                                     array1d <localIndex> const & uniqueFaceOffsets,
                                     Group & cellBlocks )
{
  localIndex const numNodes = lowestNodeToFaces.size();
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
    for( localIndex j = 0, curFaceID = uniqueFaceOffsets[nodeID]; j < numFaces; ++j, ++curFaceID )
    {
      const NodesAndElementOfFace & f0 = lowestNodeToFaces( nodeID, j );
      CellBlock & cb0 = cellBlocks.getGroup< CellBlock >( f0.iCellBlock );
      cb0.setElementToFaces( f0.element, f0.iFace, curFaceID );

      // If the following face exists and is identical to the current one,
      // Then we insert the following face with the same face id
      // and remove it from the next iteration (it's already inserted).
      if( j < numFaces - 1 )
      {
        const NodesAndElementOfFace & f1 = lowestNodeToFaces( nodeID, j + 1 );
        if( f0 == f1 )
        {
          CellBlock & cb1 = cellBlocks.getGroup< CellBlock >( f1.iCellBlock );
          cb1.setElementToFaces( f1.element, f1.iFace, curFaceID );
          ++j;
        }
      }
    }
  }
}

/**
 * @brief Fills the element to edges mappings of all the cells provided through @p cellBlocks.
 * @param faceToEdges We need the face to edges mapping to get some edge index.
 * @param cellBlocks The cell blocks for which the mappings will be constructed.
 */
void fillElementToEdgesOfCellBlocks( ArrayOfArrays< localIndex > const & faceToEdges,
                                     Group & cellBlocks )
{
  for( localIndex iCellBlock = 0; iCellBlock < cellBlocks.numSubGroups(); ++iCellBlock )
  {
    CellBlock & cellBlock = cellBlocks.getGroup< CellBlock >( iCellBlock );
    array2d< localIndex > const & cellToFaces = cellBlock.getElemToFaces();

    // We build the edges of each face of each cell,
    // so we can construct the cell to edges mapping.
    // Some specific care is required not to insert edges twice (faces share edges).
    // Another implementation (used in other contexts) would use some edge signature
    // to remove the duplicates.

    // Loop over the cells
    for( localIndex kc = 0; kc < cellBlock.size(); kc++ )
    {
      int count = 0;
      for( localIndex kf = 0; kf < cellBlock.numFacesPerElement(); kf++ )
      {
        // Loop over edges of each face
        localIndex faceIndex = cellToFaces[kc][kf];
        for( localIndex ke = 0; ke < faceToEdges.sizeOfArray( faceIndex ); ke++ )
        {
          bool isUnique = true;
          localIndex edgeIndex = faceToEdges[faceIndex][ke];

          // Loop over edges that have already been added to the element.
          for( localIndex kec = 0; kec < count + 1; kec++ )
          {
            // make sure that the edge has not been counted yet
            if( cellBlock.hasElementToEdges( kc, kec, edgeIndex ) )
            {
              isUnique = false;
              break;
            }
          }
          if( isUnique )
          {
            cellBlock.setElementToEdges( kc, count, edgeIndex );
            count++;
          }
        } // end edge loop
      } // end face loop
    } // end cell loop
  }
}

void CellBlockManager::buildFaceMaps( localIndex numNodes )
{
  const ArrayOfArrays< NodesAndElementOfFace > lowestNodeToFaces = createLowestNodeToFaces( numNodes, this->getCellBlocks() );

  array1d< localIndex > uniqueFaceOffsets( numNodes + 1 );
  m_numFaces = calculateTotalNumberOfFaces( lowestNodeToFaces.toViewConst(), uniqueFaceOffsets );

  resizeFaceMaps( lowestNodeToFaces.toViewConst(),
                  uniqueFaceOffsets,
                  m_faceToNodes,
                  m_faceToEdges,
                  m_faceToElements );

  populateFaceMaps( lowestNodeToFaces.toViewConst(),
                    uniqueFaceOffsets,
                    m_faceToNodes,
                    m_faceToElements );

  fillElementToFacesOfCellBlocks( lowestNodeToFaces, uniqueFaceOffsets, this->getCellBlocks() );
}

localIndex buildEdgeMaps( localIndex numNodes,
                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                          ArrayOfArrays< localIndex > & faceToEdgeMap,
                          ArrayOfSets< localIndex > & edgeToFaceMap,
                          array2d< localIndex > & edgeToNodeMap )
{
  ArrayOfArrays< EdgeBuilder > edgesByLowestNode = createEdgesByLowestNode( numNodes, faceToNodeMap );

  array1d< localIndex > uniqueEdgeOffsets( numNodes + 1 );
  localIndex const numEdges = calculateTotalNumberOfEdges( edgesByLowestNode.toViewConst(), uniqueEdgeOffsets );

  resizeEdgeMaps( edgesByLowestNode.toViewConst(),
                  uniqueEdgeOffsets,
                  edgeToNodeMap,
                  edgeToFaceMap );

  populateEdgeMaps( edgesByLowestNode.toViewConst(),
                    uniqueEdgeOffsets,
                    faceToNodeMap,
                    faceToEdgeMap,
                    edgeToFaceMap,
                    edgeToNodeMap );

  return numEdges;
}

void CellBlockManager::buildNodeToEdges( localIndex numNodes )
{
  ArrayOfArrays< localIndex > toEdgesTemp( numNodes, maxEdgesPerNode );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeEdges = 0;

  forAll< parallelHostPolicy >( m_numEdges, [&]( localIndex const edgeID )
  {
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( m_edgeToNodes( edgeID, 0 ), edgeID );
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( m_edgeToNodes( edgeID, 1 ), edgeID );
    totalNodeEdges += 2;
  } );

  // Resize the node to edge map.
  m_nodeToEdges.resize( 0 );

  // Reserve space for the number of current nodes plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  m_nodeToEdges.reserve( entriesToReserve );

  localIndex getEdgeMapOverallocation = 8; // TODO deal with this.
  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeEdges.get() + numNodes * getEdgeMapOverallocation * ( 1 + 2 * overAllocationFactor );
  m_nodeToEdges.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    m_nodeToEdges.appendSet( toEdgesTemp.sizeOfArray( nodeID ) + getEdgeMapOverallocation );
  }

  ArrayOfSetsView< localIndex > const & toEdgesView = m_nodeToEdges.toView(); // FIXME why a view here?
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const edges = toEdgesTemp[ nodeID ];
    localIndex const numNodeEdges = toEdgesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueEdges = LvArray::sortedArrayManipulation::makeSortedUnique( edges, edges + numNodeEdges );
    toEdgesView.insertIntoSet( nodeID, edges, edges + numUniqueEdges );
  } );
}

void CellBlockManager::buildMaps( localIndex numNodes )
{
  buildFaceMaps( numNodes );
  m_numEdges = buildEdgeMaps( numNodes,
                              m_faceToNodes.toViewConst(),
                              m_faceToEdges,
                              m_edgeToFaces,
                              m_edgeToNodes );
  buildNodeToEdges( numNodes );

  fillElementToEdgesOfCellBlocks( m_faceToEdges, this->getCellBlocks() );
}

//TODO return views
ArrayOfArrays< localIndex > CellBlockManager::getFaceToNodes() const
{
  return m_faceToNodes;
}

ArrayOfSets< localIndex > CellBlockManager::getNodeToFaces(localIndex numNodes) const // TODO remove numNodes
{
  localIndex const numFaces = m_faceToNodes.size();
  const localIndex maxFacesPerNode = 200; // TODO deal with this differently
  ArrayOfArrays< localIndex > nodeToFacesTemp( numNodes, maxFacesPerNode );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeFaces = 0;

  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const faceID )
  {
    localIndex const numFaceNodes = m_faceToNodes.sizeOfArray( faceID );
    totalNodeFaces += numFaceNodes;
    for( localIndex a = 0; a < numFaceNodes; ++a )
    {
      nodeToFacesTemp.emplaceBackAtomic< parallelHostAtomic >( m_faceToNodes( faceID, a ), faceID );
    }
  } );

  ArrayOfSets< localIndex > result;

  const localIndex faceMapOverallocation = 8; // TODO deal with this differently
  const localIndex nodeMapExtraSpacePerFace = 4; // TODO deal with this differently
  // Resize the node to face map.
  result.resize( 0 );

  // Reserve space for the number of nodes faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  result.reserve( entriesToReserve );

  // Reserve space for the total number of node faces + extra space for existing nodes + even more space for new nodes.
  localIndex const valuesToReserve = totalNodeFaces.get() + numNodes * nodeMapExtraSpacePerFace * ( 1 + 2 * overAllocationFactor );
  result.reserveValues( valuesToReserve );

  // Append the individual arrays.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    result.appendSet( nodeToFacesTemp.sizeOfArray( nodeID ) + faceMapOverallocation );
  }

  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const faces = nodeToFacesTemp[ nodeID ];
    localIndex const numNodeFaces = nodeToFacesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueFaces = LvArray::sortedArrayManipulation::makeSortedUnique( faces, faces + numNodeFaces );
    result.insertIntoSet( nodeID, faces, faces + numUniqueFaces );
  } );

  return result;
}

// TODO return views
array2d< localIndex > CellBlockManager::getFaceToElements() const
{
  return m_faceToElements;
}

const Group & CellBlockManager::getCellBlocks() const
{
  return this->getGroup( cellBlocksKey );
}

Group & CellBlockManager::getCellBlocks()
{
  return this->getGroup( cellBlocksKey );
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

localIndex CellBlockManager::numFaces() const
{
  return m_numFaces;
}

ArrayOfSets< geosx::localIndex > const & CellBlockManager::getEdgeToFaces() const
{
  return m_edgeToFaces;
}

array2d< geosx::localIndex > const & CellBlockManager::getEdgeToNodes() const
{
  return m_edgeToNodes;
}

ArrayOfArrays< geosx::localIndex > const & CellBlockManager::getFaceToEdges() const
{
  return m_faceToEdges;
}

ArrayOfSets< localIndex > CellBlockManager::getNodeToEdges() const
{
  return m_nodeToEdges;
}

localIndex CellBlockManager::numEdges() const
{
  return m_numEdges;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlockManager, string const &, Group * const )
}
