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

#include "CellBlockManager.hpp"

#include "CellBlockUtilities.hpp"

#include <algorithm>

namespace geosx
{
using namespace dataRepository;

CellBlockManager::CellBlockManager( string const & name, Group * const parent ):
  CellBlockManagerABC( name, parent ),
  m_nodesPositions( 0, 3 )
{
  this->registerGroup< Group >( viewKeyStruct::cellBlocks() );
}

void CellBlockManager::resize( integer_array const & numElements,
                               string_array const & regionNames )
{
  localIndex const numRegions = LvArray::integerConversion< localIndex >( regionNames.size());
  for( localIndex reg=0; reg<numRegions; ++reg )
  {
    this->getCellBlock( regionNames[reg] ).resize( numElements[reg] );
  }
}

Group * CellBlockManager::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

ArrayOfArrays< localIndex > CellBlockManager::getNodeToElements() const
{
  // The function works in three steps.
  // First finds the number of elements attached to each node.
  // Then second step allocates the vectors, including extra allocations.
  // Last, the output vector is filled.

  // First step: how many elements for each node, stored in the elemsPerNode array.
  array1d< localIndex > elemsPerNode( m_numNodes );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeElems = 0;

  for( localIndex iCellBlock = 0; iCellBlock < numCellBlocks(); ++iCellBlock )
  {
    const CellBlockABC & cb = this->getCellBlock( iCellBlock );
    array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & elemToNode = cb.getElemToNodes();
    forAll< parallelHostPolicy >( cb.numElements(), [&elemsPerNode, totalNodeElems, &elemToNode, &cb] ( localIndex const k )
    {
      localIndex const numNodesPerElement = cb.numNodesPerElement();
      totalNodeElems += numNodesPerElement;
      for( localIndex a = 0; a < numNodesPerElement; ++a )
      {
        localIndex const nodeIndex = elemToNode( k, a );
        RAJA::atomicInc< parallelHostAtomic >( &elemsPerNode[ nodeIndex ] );
      }
    } );
  }

  // Second, allocation.
  ArrayOfArrays< localIndex > result;

  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * m_numNodes;
  result.reserve( entriesToReserve );

  localIndex const valuesToReserve = totalNodeElems.get() + m_numNodes * getElemMapOverAllocation() * ( 1 + 2 * overAllocationFactor );
  result.reserveValues( valuesToReserve );

  // Append an array for each node with capacity to hold the appropriate number of elements plus some wiggle room.
  for( localIndex nodeID = 0; nodeID < m_numNodes; ++nodeID )
  {
    result.appendArray( 0 );
    result.setCapacityOfArray( nodeID, elemsPerNode[ nodeID ] + getElemMapOverAllocation() );
  }

  // Third, filling the result.
  for( localIndex iCellBlock = 0; iCellBlock < numCellBlocks(); ++iCellBlock )// Do not need index
  {
    const CellBlockABC & cb = this->getCellBlock( iCellBlock );
    array2d< localIndex, cells::NODE_MAP_PERMUTATION > const elemToNode = cb.getElemToNodes();
    for( localIndex iElem = 0; iElem < cb.numElements(); ++iElem )
    {
      for( localIndex iNode = 0; iNode < cb.numNodesPerElement(); ++iNode )
      {
        localIndex const nodeIndex = elemToNode( iElem, iNode );
        result.emplaceBack( nodeIndex, iElem );
      }
    }
  }

  return result;
}

/**
 * @brief Convenience programming structure that holds the nodes for a given Face. No information about which cell though.
 *
 * This structure holds `<` and `==` operators such that equal instances in a sorted array
 * are meant to be consecutive. (see `std::unique` for example).
 * Two instances with same nodes in different order are considered equal.
 * Such a case often happen since the same surface is shared by two adjacent elements
 * in conformal meshes.
 */

// TODO Evaluate the cost of this structure 
// The number of nodes is either 3 or 4 - using array1d is likely an overkill 
struct NodesAndElementOfFace
{
  NodesAndElementOfFace( localIndex* nodes_, localIndex element_, localIndex iCellBlock_, localIndex iFace_ ):
    nodes{ nodes_[0],nodes_[1],nodes_[2],nodes_[3] },
    element( element_ ),
    iCellBlock( iCellBlock_ ),
    iFace( iFace_ ),
    sortedNodes { nodes_[0],nodes_[1],nodes_[2],nodes_[3] }
  {
    
    if( sortedNodes[0] > sortedNodes[1] ) std::swap( sortedNodes[0], sortedNodes[1] );
    if( sortedNodes[2] > sortedNodes[3] ) std::swap( sortedNodes[2], sortedNodes[3] );
    if( sortedNodes[0] > sortedNodes[2] ) std::swap( sortedNodes[0], sortedNodes[2] );
    if( sortedNodes[1] > sortedNodes[3] ) std::swap( sortedNodes[1], sortedNodes[3] );
    if( sortedNodes[1] > sortedNodes[2] ) std::swap( sortedNodes[1], sortedNodes[2] );


    // TODO Replace by tests - There are 3 or 4 values
    //std::sort( sortedNodes.begin(), sortedNodes.end() );
  }

  /**
   * @brief Imposes an ordering on NodesAndElementOfFace.
   * @param [in] rhs the NodesAndElementOfFace to compare against.
   * @return a boolean.
   */
  bool operator<( NodesAndElementOfFace const & rhs ) const
  {
    // Using some standard comparison like vector::operator<.
    // Two subsequent NodesAndElementOfFace may still be equal
    // We are consistent with operator==, which is what we require.
    
    if(sortedNodes[0] != rhs.sortedNodes[0]) return sortedNodes[0] < rhs.sortedNodes[0];
    if(sortedNodes[1] != rhs.sortedNodes[1]) return sortedNodes[1] < rhs.sortedNodes[1];
    if(sortedNodes[2] != rhs.sortedNodes[2]) return sortedNodes[2] < rhs.sortedNodes[2];
    if(sortedNodes[3] != rhs.sortedNodes[3]) return sortedNodes[3] < rhs.sortedNodes[3];
    return false;
    
    //return std::lexicographical_compare( sortedNodes.begin(), sortedNodes.end(),
                                        // rhs.sortedNodes.begin(), rhs.sortedNodes.end() );
  }

  /**
   * @brief Two NodesAndElementOfFace instances are considered equal if they share the same node set, whatever the order.
   * @param [in] rhs the NodesAndElementOfFace to compare against.
   * @return a boolean.
   */
  bool operator==( NodesAndElementOfFace const & rhs ) const
  {
    // Comparing term by term like STL does.

   return (sortedNodes[0] == rhs.sortedNodes[0]) &&
          (sortedNodes[1] == rhs.sortedNodes[1]) &&
          (sortedNodes[2] == rhs.sortedNodes[2]) &&
          (sortedNodes[3] == rhs.sortedNodes[3]);

    //return ( sortedNodes.size() == rhs.sortedNodes.size() && std::equal( sortedNodes.begin(), sortedNodes.end(), rhs.sortedNodes.begin() ) );
  }

  /// The list of nodes describing the face.
  localIndex nodes[4] = {-1, -1, -1, -1};
  // array1d< localIndex > nodes;

  /**
   * @brief The element to which this face belongs.
   *
   * Each face may belong to multiple elements.
   * But during the identification process (we loop on each face of each element),
   * we store the cell we are iterating on.
   * The we'll be able to identify the duplicated faces because we also have the nodes.
   */
  localIndex element;

  /**
   * @brief Cell block index
   *
   * During the process, we need to know form which cell block the instance was created.
   */
  localIndex iCellBlock;

  /**
   * @brief Face index
   *
   * During the process, we need to know what was the face index when this instance was created.
   */
  localIndex iFace;

private:
  /// Sorted nodes describing the face; only for comparison/sorting reasons.
  
  // TODO Think about storing the permutation of the input nodes 
  localIndex sortedNodes[4]={ -1, -1, -1, -1 };
  // array1d< localIndex > sortedNodes;
};

/**
 * @brief Return the total number of unique faces and fill in the uniqueFaceOffsets array.
 * @param [in] lowestNodeToFaces An array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
 * @param [out] uniqueFaceOffsets An array of size numNodes + 1. After this function returns node i contains
 *              faces with IDs ranging from uniqueFaceOffsets[ i ] to uniqueFaceOffsets[ i + 1 ] - 1.
 * @return return total number of faces
 */
// TODO Why is this separated from createLowestNodeToFaces ? 
// TODO Flag interior - boundary faces - compute areas - normals - could also be done 
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
      if( j < numFaces - 1 )
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
 * @brief Copies the nodes from @p nodesAndElementOfFace into @p faceToNodes[@p faceID ].
 * @param [in] faceID The face index.
 * @param [in] nodesAndElementOfFace The nodes and element for @p faceID.
 * @param [out] faceToNodes No input data is used.
 */
void insertFaceToNodesEntry( localIndex const faceID,
                             NodesAndElementOfFace const & nodesAndElementOfFace,
                             ArrayOfArrays< localIndex > & faceToNodes )
{
  localIndex const numFaceNodes = 4; // nodesAndElementOfFace.nodes.size();

  // FIXME The size should be OK because it's been allocated previously.
  for( localIndex i = 0; i < numFaceNodes; ++i )
  {
    faceToNodes[faceID][i] = nodesAndElementOfFace.nodes[i];
  }
  GEOSX_ASSERT_EQ( numFaceNodes, faceToNodes.sizeOfArray( faceID ) );
  GEOSX_DEBUG_VAR( numFaceNodes );
}

/**
 * @brief Fills the face to nodes map and face to element maps
 * @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
 * @param [in] uniqueFaceOffsets an array containing the unique ID of the first face associated with each node.
 * @param [inout] faceToElements the face to element map.
 * @param [inout] faceToNodes the face to node map.
 */
// TODO  The faceToNode and nodeToFace maps should be the first filled and should have been 
// before so that we do not carry around this annoying lowestNodeToFaces and uniqueFaceOffsets
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
    for( localIndex j = 0, curFaceID = uniqueFaceOffsets[nodeID]; j < numFaces; ++j, ++curFaceID )
    {
      // If two subsequent `NodesAndElementOfFace` compare equal then they describe an interior face.
      // It must therefore be considered "twice".
      // Otherwise it's a boundary face, and "once" is enough.
      NodesAndElementOfFace const & f0 = lowestNodeToFaces( nodeID, j );
      
      //insertFaceToNodesEntry( curFaceID, f0, faceToNodes );

      faceToNodes[curFaceID][0] = f0.nodes[0];
      faceToNodes[curFaceID][1] = f0.nodes[1];
      faceToNodes[curFaceID][2] = f0.nodes[2];
      faceToNodes[curFaceID][3] = f0.nodes[3];  

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
  GEOSX_UNUSED_VAR( lowestNodeToFaces );

  //localIndex const numNodes = lowestNodeToFaces.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();

  // TODO  The faceToNode and nodeToFace maps should be the first filled and should have been 
  // before so that we do not carry around this annoying lowestNodeToFaces and uniqueFaceOffsets

//  array1d< localIndex > numNodesPerFace( numUniqueFaces );
//  RAJA::ReduceSum< parallelHostReduce, localIndex > totalFaceNodes( 0.0 );

  // Compute the size needed to store the faceToNodeMap - hexahedron nbFaces x 4
  
  // loop over all the nodes.
  // TODO optimize
/*  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
    // loop over all the NodesAndElementOfFace associated with the node
    for( localIndex j = 0, curFaceID = uniqueFaceOffsets[ nodeID ]; j < numFaces; ++j, ++curFaceID )
    {
      const NodesAndElementOfFace & f0 = lowestNodeToFaces( nodeID, j );
      numNodesPerFace[curFaceID] = 4; // f0.nodes.size();
      totalFaceNodes += numNodesPerFace[curFaceID];

      if( ( j < numFaces - 1 ) and ( f0 == lowestNodeToFaces( nodeID, j + 1 ) ) )
      {
        ++j;
      }
    }
  } );
*/

  // TODO The same allocation should be done for faceToNode and faceToEdge


  // There are 3 or 4 nodes per face
  // A code that works for general polygons is not efficient for hexahedra
  // TODO Optimize

  // Resize the face to node map.
  // faceToNodeMap.resize( numUniqueFaces, 4 ); // TODO How can we allocate what we need ? // This does not do it
  
    // Reserve space for the number of current faces plus some extra.
//  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = numUniqueFaces; //TODO why this allocation factor
  
  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = 4 * numUniqueFaces; // totalFaceNodes.get() + numUniqueFaces * CellBlockManager::nodeMapExtraSpacePerFace() * ( 1 + 2 * overAllocationFactor );


  faceToNodeMap.reserve( entriesToReserve );
  faceToNodeMap.reserveValues( valuesToReserve );
  // faceToNodeMap.reserveValues( 2 * entriesToReserve ); //TODO I don"t undertand anythin about LVARRAY :@ 

  // I do not understand why do we have to put it the arrays 
  // Should be the same for Elem and Edges ? I do not get it ....
  // Append the individual arrays.
  for( localIndex faceID = 0; faceID < numUniqueFaces; ++faceID )
  {
    faceToNodeMap.appendArray( 4 ); //  numNodesPerFace[ faceID ] );
    //faceToNodeMap.setCapacityOfArray( faceToNodeMap.size() - 1, 4); 
                                      //numNodesPerFace[ faceID ] + CellBlockManager::nodeMapExtraSpacePerFace() );
  } 

  // Each face may belong to _maximum_ 2 elements.
  // If it belongs to only one, we put `-1` for the undefined value.
  faceToElemMap.resize( numUniqueFaces, 2 );

  // The number of edges per face is 3 or 4 - let's go with 4 
  // TODO adequate strategy for triangles
  faceToEdgesMap.resize( numUniqueFaces, 4 );
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
  // TODO - Use information on the mesh to avoid this overkill allocation
  // For meshes of the InternalMeshGenerator - regular cartesian - we are good with 6 faces per node - no big gain compared to 24 
  ArrayOfArrays< NodesAndElementOfFace > lowestNodeToFaces( numNodes, 6 );
  //ArrayOfArrays< NodesAndElementOfFace > lowestNodeToFaces( numNodes, 2 * CellBlockManager::maxFacesPerNode() );

  // The function is a bit simplified and is not run in parallel anymore.
  // Can be improved.
  for( localIndex iCellBlock = 0; iCellBlock < cellBlocks.numSubGroups(); ++iCellBlock )
  {
    const CellBlock & cb = cellBlocks.getGroup< CellBlock >( iCellBlock );
    localIndex const numFacesPerElement = cb.numFacesPerElement();
    localIndex const numElements = cb.numElements();

    for( localIndex iElement = 0; iElement < numElements; ++iElement )
    {
      // Looping on the faces of the cell
      for( localIndex iFace = 0; iFace < numFacesPerElement; ++iFace )
      {
        // Get all the nodes of the face
        localIndex nodesInFace[4];
        // For now suppose everything goes well - all cells having 4 nodes
        cb.getFaceNodes( iElement, iFace, nodesInFace );
        // Fill the result with the collected data.       
        // localIndex const & lowestNode = *std::min_element( nodesInFace.begin(), nodesInFace.end() );
        localIndex lowestNode = std::min( std::min(nodesInFace[0],nodesInFace[1]), std::min(nodesInFace[2], nodesInFace[3]) );
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
 * @brief Filling the elements to faces maps in the cell blocks.
 * @param lowestNodeToFaces The lowest node to faces information array.
 * @param uniqueFaceOffsets The unique face offsets.
 * @param cellBlocks The cell blocks for which we need to compute the element to faces mappings.
 *
 * @note @p lowestNodeToFaces and @p uniqueFaceOffsets are better described in the documentations of the functions that build them.
 */
void fillElementToFacesOfCellBlocks( ArrayOfArrays< NodesAndElementOfFace > const & lowestNodeToFaces,
                                     array1d< localIndex > const & uniqueFaceOffsets,
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
// What is the usefulness of the CellBlock concept ?
void fillElementToEdgesOfCellBlocks( ArrayOfArrays< localIndex > const & faceToEdges,
                                     Group & cellBlocks )
{
  for( localIndex iCellBlock = 0; iCellBlock < cellBlocks.numSubGroups(); ++iCellBlock )
  {
    CellBlock & cellBlock = cellBlocks.getGroup< CellBlock >( iCellBlock );
    arrayView2d< localIndex const > const cellToFaces = cellBlock.getElemToFacesConstView();

    // We build the edges of each face of each cell,
    // so we can construct the cell to edges mapping.
    // Some specific care is required not to insert edges twice (faces share edges).
    // Another implementation (used in other contexts) would use some edge signature
    // to remove the duplicates.

    // TODO Optimize the list of edges is known for the types of cells
    // managed by the CellBlocks
    // Access to the global edge index could be gotten through the node

    // Loop over the cells
    for( localIndex kc = 0; kc < cellBlock.numElements(); kc++ )
    {
      int count = 0;
      // TODO Remove this loop 
      for( localIndex kf = 0; kf < cellBlock.numFacesPerElement(); kf++ )
      {
        // Loop over edges of each face
        localIndex const faceIndex = cellToFaces[kc][kf];
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

void CellBlockManager::buildFaceMaps()
{
  const ArrayOfArrays< NodesAndElementOfFace > lowestNodeToFaces = createLowestNodeToFaces( m_numNodes, this->getCellBlocks() );

  array1d< localIndex > uniqueFaceOffsets( m_numNodes + 1 );
  m_numFaces = calculateTotalNumberOfFaces( lowestNodeToFaces.toViewConst(), uniqueFaceOffsets );

  std::cout << " Number of faces  " << m_numFaces << std::endl;

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

void CellBlockManager::buildNodeToEdges()
{
  ArrayOfArrays< localIndex > toEdgesTemp( m_numNodes, maxEdgesPerNode() );
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
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * m_numNodes;
  m_nodeToEdges.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeEdges.get() + m_numNodes * getEdgeMapOverallocation() * ( 1 + 2 * overAllocationFactor );
  m_nodeToEdges.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex nodeID = 0; nodeID < m_numNodes; ++nodeID )
  {
    m_nodeToEdges.appendSet( toEdgesTemp.sizeOfArray( nodeID ) + getEdgeMapOverallocation() );
  }

  ArrayOfSetsView< localIndex > const & toEdgesView = m_nodeToEdges.toView(); // FIXME why a view here?
  forAll< parallelHostPolicy >( m_numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const edges = toEdgesTemp[ nodeID ];
    localIndex const numNodeEdges = toEdgesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueEdges = LvArray::sortedArrayManipulation::makeSortedUnique( edges, edges + numNodeEdges );
    toEdgesView.insertIntoSet( nodeID, edges, edges + numUniqueEdges );
  } );
}

void CellBlockManager::buildMaps()
{
  buildFaceMaps();
  m_numEdges = buildEdgeMaps( m_numNodes,
                              m_faceToNodes.toViewConst(),
                              m_faceToEdges,
                              m_edgeToFaces,
                              m_edgeToNodes );
  buildNodeToEdges();

  fillElementToEdgesOfCellBlocks( m_faceToEdges, this->getCellBlocks() );
}

ArrayOfArrays< localIndex > CellBlockManager::getFaceToNodes() const
{
  std::cout << " This is a face " <<
             m_faceToNodes[0][0] << "  " <<
             m_faceToNodes[0][1] << "  " <<
             m_faceToNodes[0][2] << "  " <<
             m_faceToNodes[0][3] << "  " << std::endl;
  return m_faceToNodes;
}

ArrayOfSets< localIndex > CellBlockManager::getNodeToFaces() const
{
  localIndex const numFaces = m_faceToNodes.size();
  ArrayOfArrays< localIndex > nodeToFacesTemp( m_numNodes, maxFacesPerNode() );
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

  // Resize the node to face map.
  result.resize( 0 );

  // Reserve space for the number of nodes faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * m_numNodes;
  result.reserve( entriesToReserve );

  // Reserve space for the total number of node faces + extra space for existing nodes + even more space for new nodes.
  localIndex const valuesToReserve = totalNodeFaces.get() + m_numNodes * nodeMapExtraSpacePerFace() * ( 1 + 2 * overAllocationFactor );
  result.reserveValues( valuesToReserve );

  // Append the individual arrays.
  for( localIndex nodeID = 0; nodeID < m_numNodes; ++nodeID )
  {
    result.appendSet( nodeToFacesTemp.sizeOfArray( nodeID ) + getFaceMapOverallocation() );
  }

  forAll< parallelHostPolicy >( m_numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const faces = nodeToFacesTemp[ nodeID ];
    localIndex const numNodeFaces = nodeToFacesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueFaces = LvArray::sortedArrayManipulation::makeSortedUnique( faces, faces + numNodeFaces );
    result.insertIntoSet( nodeID, faces, faces + numUniqueFaces );
  } );

  return result;
}

array2d< localIndex > CellBlockManager::getFaceToElements() const
{
  return m_faceToElements;
}

const Group & CellBlockManager::getCellBlocks() const
{
  return this->getGroup( viewKeyStruct::cellBlocks() );
}

Group & CellBlockManager::getCellBlocks()
{
  return this->getGroup( viewKeyStruct::cellBlocks() );
}

localIndex CellBlockManager::numNodes() const
{
  return m_numNodes;
}

localIndex CellBlockManager::numCellBlocks() const
{
  return this->getCellBlocks().numSubGroups();
}

const CellBlockABC & CellBlockManager::getCellBlock( localIndex iCellBlock ) const
{
  return this->getCellBlocks().getGroup< const CellBlockABC >( iCellBlock );
}

localIndex CellBlockManager::numFaces() const
{
  return m_numFaces;
}

ArrayOfSets< geosx::localIndex > CellBlockManager::getEdgeToFaces() const
{
  return m_edgeToFaces;
}

array2d< geosx::localIndex > CellBlockManager::getEdgeToNodes() const
{
  return m_edgeToNodes;
}

ArrayOfArrays< geosx::localIndex > CellBlockManager::getFaceToEdges() const
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

CellBlock & CellBlockManager::registerCellBlock( string name )
{
  return this->getCellBlocks().registerGroup< CellBlock >( name );
}

array2d< real64, nodes::REFERENCE_POSITION_PERM > CellBlockManager::getNodesPositions() const
{
  return m_nodesPositions;
}

arrayView2d< real64, nodes::REFERENCE_POSITION_USD > CellBlockManager::getNodesPositions()
{
  return m_nodesPositions.toView();
}

void CellBlockManager::setNumNodes( localIndex numNodes )
{
  m_numNodes = numNodes;
  m_nodesPositions.resize( m_numNodes );
  m_nodeLocalToGlobal.resize( m_numNodes );
  m_nodeLocalToGlobal.setValues< serialPolicy >( -1 );
}

array1d< globalIndex > CellBlockManager::getNodeLocalToGlobal() const
{
  return m_nodeLocalToGlobal;
}

arrayView1d< globalIndex > CellBlockManager::getNodeLocalToGlobal()
{
  return m_nodeLocalToGlobal.toView();
}

std::map< string, SortedArray< localIndex > > const & CellBlockManager::getNodeSets() const
{
  return m_nodeSets;
}

std::map< string, SortedArray< localIndex > > & CellBlockManager::getNodeSets()
{
  return m_nodeSets;
}

}
