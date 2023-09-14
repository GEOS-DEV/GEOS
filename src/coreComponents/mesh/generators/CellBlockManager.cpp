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

#include "mesh/generators/CellBlockUtilities.hpp"
#include "mesh/generators/LineBlock.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"
#include "RAJA/util/mutex.hpp"

#include <tuple>
#include <algorithm>

namespace geos
{
using namespace dataRepository; 
CellBlockManager::CellBlockManager( string const & name, Group * const parent ):
  CellBlockManagerABC( name, parent ),
  m_nodesPositions( 0, 3 )
{
  this->registerGroup< Group >( viewKeyStruct::cellBlocks() );
  this->registerGroup< Group >( viewKeyStruct::faceBlocks() );
  this->registerGroup< Group >( viewKeyStruct::lineBlocks() );
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

Group * CellBlockManager::createChild( string const & GEOS_UNUSED_PARAM( childKey ), string const & GEOS_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

/// Element identifier containing (block index, cell index).
using CellBlockIndexPair = std::pair< localIndex, localIndex >;

template< typename POLICY >
void convertFromCellBlockPairMap( ArrayOfArraysView< CellBlockIndexPair const > const & srcMap,
                                  ToCellRelation< ArrayOfArrays< localIndex > > & dstMap )
{
  ArrayOfArrays< localIndex > & toBlock = dstMap.toBlockIndex;
  ArrayOfArrays< localIndex > & toCell = dstMap.toCellIndex;

  localIndex const numObjects = srcMap.size();

  toBlock.resizeFromOffsets( numObjects, srcMap.toViewConst().getOffsets() );
  toCell.resizeFromOffsets( numObjects, srcMap.toViewConst().getOffsets() );

  forAll< parallelHostPolicy >( numObjects, [toBlock = toBlock.toView(),
                                             toCell = toCell.toView(),
                                             srcMap]( localIndex const objIndex )
  {
    arraySlice1d< CellBlockIndexPair const > const cells = srcMap[ objIndex ];
    for( CellBlockIndexPair const & e : cells )
    {
      toBlock.emplaceBack( objIndex, std::get< 0 >( e ) );
      toCell.emplaceBack( objIndex, std::get< 1 >( e ) );
    }
  } );
}

template< typename POLICY >
void convertFromCellBlockPairMap( ArrayOfArraysView< CellBlockIndexPair const > const & srcMap,
                                  ToCellRelation< array2d< localIndex > > & dstMap )
{
  array2d< localIndex > & toBlock = dstMap.toBlockIndex;
  array2d< localIndex > & toCell = dstMap.toCellIndex;

  localIndex const numObjects = srcMap.size();
  localIndex const maxNumElem = toCell.size( 1 );

  GEOS_ERROR_IF_NE( toBlock.size( 1 ), maxNumElem );
  GEOS_ERROR_IF_NE( toCell.size( 1 ), maxNumElem );

  toBlock.resizeDimension< 0 >( numObjects );
  toCell.resizeDimension< 0 >( numObjects );

  // We allow a fixed-size map to represent a variable relationship, as long as
  // the number of elements does not exceed the fixed size (set by the caller).
  // In this case, a dummy value "-1" is used to represent non-present elements.
  toBlock.setValues< POLICY >( -1 );
  toCell.setValues< POLICY >( -1 );

  forAll< parallelHostPolicy >( numObjects, [=, // needed to optionally capture maxNumElem in Debug
                                             toBlock = toBlock.toView(),
                                             toCell = toCell.toView()]( localIndex const objIndex )
  {
    arraySlice1d< CellBlockIndexPair const > const cells = srcMap[ objIndex ];
    GEOS_ASSERT_GE( maxNumElem, cells.size() );
    localIndex count = 0;
    for( CellBlockIndexPair const & e : cells )
    {
      toBlock( objIndex, count ) = std::get< 0 >( e );
      toCell( objIndex, count ) = std::get< 1 >( e );
      ++count;
    }
  } );
}

/**
 * @brief Build to-cell map by inverting existing maps in cell blocks.
 * @tparam BASEMAP underlying type of to-cell map
 * @tparam FUNC type of @p cellToObjectGetter
 * @param cellIndex number of objects (nodes, faces, etc.) for which maps are built
 * @param toCells container for to-cell (block, index) maps
 * @param cellToObjectGetter function used to extract maps from subregions
 * @param overAlloc overallocation for the resulting maps
 *                  (extra capacity per row, only meaningful for variable-secondary-size containers)
 */
template< typename BASEMAP, typename FUNC >
void CellBlockManager::buildToCellMap( localIndex const numObjects,
                                       ToCellRelation< BASEMAP > & toCells,
                                       FUNC cellToObjectGetter,
                                       localIndex const overAlloc ) const
{
  // Calculate the number of entries in each sub-array
  array1d< localIndex > counts( numObjects );
  counts.setValues< serialPolicy >( overAlloc );

  for( localIndex blockIndex = 0; blockIndex < numCellBlocks(); ++blockIndex )
  {
    CellBlock const & cb = getCellBlock( blockIndex );

    auto const cellToObject = cellToObjectGetter( cb );
    forAll< parallelHostPolicy >( cb.size(), [counts = counts.toView(),
                                              cellToObject]( localIndex const ei )
    {
      auto const objects = cellToObject[ ei ];
      // can't use range-based for loop when slice is not contiguous
      for( localIndex i = 0; i < objects.size(); ++i )
      {
        RAJA::atomicInc< parallelHostAtomic >( &counts[ objects[i] ] );
      }
    } );
  }
  ;

  // Allocate memory
  ArrayOfArrays< CellBlockIndexPair > cellBlockPairList;
  cellBlockPairList.resizeFromCapacities< parallelHostPolicy >( numObjects, counts.data() );

  // Populate map of tuples in parallel
  for( localIndex blockIndex = 0; blockIndex < numCellBlocks(); ++blockIndex )
  {
    CellBlock const & cb = getCellBlock( blockIndex );

    auto const cellToObject = cellToObjectGetter( cb );
    forAll< parallelHostPolicy >( cb.size(), [cellBlockPairList = cellBlockPairList.toView(),
                                              cellToObject,
                                              blockIndex]( localIndex const cellIndex )
    {
      auto const objects = cellToObject[ cellIndex ];
      // can't use range-based for loop when slice is not contiguous
      for( localIndex i = 0; i < objects.size(); ++i )
      {
        cellBlockPairList.emplaceBackAtomic< parallelHostAtomic >( objects[i], blockIndex, cellIndex );
      }
    } );
  }
  ;

  // Sort each element list to ensure unique race-condition-free map order
  forAll< parallelHostPolicy >( numObjects, [cellBlockPairList = cellBlockPairList.toView()]( localIndex const objIndex )
  {
    arraySlice1d< CellBlockIndexPair > const cells = cellBlockPairList[ objIndex ];
    LvArray::sortedArrayManipulation::makeSorted( cells.begin(), cells.end() );
  } );

  // Finally, split arrays-of-tuples into separate arrays
  convertFromCellBlockPairMap< parallelHostPolicy >( cellBlockPairList.toViewConst(), toCells );
}

ToCellRelation< ArrayOfArrays< localIndex > > CellBlockManager::getNodeToElements() const
{
  ToCellRelation< ArrayOfArrays< localIndex > > result;
  buildToCellMap( m_numNodes,
                  result,
                  []( CellBlock const & cb ) { return cb.getElemToNode(); },
                  elemMapExtraSpacePerNode() );
  return result;
}

/**
 * @brief Holds information about face used in face map construction.
 *
 * These records are created and stored when visiting faces through cells.
 * Organizing and sorting them by node lists enables identification of
 * matching faces of different cells. Cell and block indices are preserved
 * so that cell-face maps may be constructed.
 */
struct NodesAndElementOfFace
{
  NodesAndElementOfFace( localIndex const duplicateFaceNodesIdx,
                         localIndex const cellIdx,
                         localIndex const blockIdx,
                         localIndex const faceNum ):
    duplicateFaceNodesIndex( duplicateFaceNodesIdx ),
    cellIndex( cellIdx ),
    blockIndex( blockIdx ),
    faceNumber( faceNum )
  {}

  localIndex duplicateFaceNodesIndex;

  /**
   * @brief Index of the cell to which this face belongs.
   *
   * Each face may belong to multiple elements.
   * But during the identification process (we loop on each face of each element),
   * we store the cell we are iterating on.
   * The we'll be able to identify the duplicated faces because we also have the nodes.
   */
  localIndex cellIndex;

  /// Cell block index of the cell to which this face belongs.
  localIndex blockIndex;

  /// Face number within a cell
  localIndex faceNumber;
};

/**
 * @brief Holds all the information used to build up the face maps.
 */
struct FaceBuilder
{
  /**
   * @brief Return a functor used to compare if two duplicate faces are the same.
   * @return A comparator that takes in two @c NodesAndElementOfFace and compares
   *   the nodes that make up the two faces.
   */
  auto duplicateFaceEquality() const
  {
    return [duplicateFaces = duplicateFaces.toViewConst()]
             ( NodesAndElementOfFace const & lhs, NodesAndElementOfFace const & rhs )
    {
      return std::equal( duplicateFaces[ lhs.duplicateFaceNodesIndex ].begin(),
                         duplicateFaces[ lhs.duplicateFaceNodesIndex ].end(),
                         duplicateFaces[ rhs.duplicateFaceNodesIndex ].begin(),
                         duplicateFaces[ rhs.duplicateFaceNodesIndex ].end() );
    };
  }

  /// A map from the lowest node index in each face to a list of duplicate faces.
  ArrayOfArrays< NodesAndElementOfFace > lowestNodeToFaces;

  /// Contains the entire node list of each duplicate face.
  ArrayOfArrays< localIndex > duplicateFaces;
};

/**
 * @brief Fills the face to nodes map and face to element maps
 * @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
 * @param [in] uniqueFaceOffsets an array containing the unique ID of the first face associated with each node.
 * @param [inout] faceToCells the face to element map.
 * @param [inout] faceToNodes the face to node map.
 */
void populateFaceMaps( Group const & cellBlocks,
                       FaceBuilder const & faceBuilder,
                       arrayView1d< localIndex const > const & uniqueFaceOffsets,
                       ArrayOfArraysView< localIndex > const & faceToNodes,
                       arrayView2d< localIndex > const & faceToCells,
                       arrayView2d< localIndex > const & faceToBlocks )
{
  ArrayOfArraysView< NodesAndElementOfFace const > const & lowestNodeToFaces =
    faceBuilder.lowestNodeToFaces.toViewConst();

  localIndex const numNodes = lowestNodeToFaces.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  GEOS_ERROR_IF_NE( uniqueFaceOffsets.size() - 1, numNodes );
  GEOS_ERROR_IF_NE( faceToNodes.size(), numUniqueFaces );
  GEOS_ERROR_IF_NE( faceToCells.size( 0 ), numUniqueFaces );
  GEOS_ERROR_IF_NE( faceToCells.size( 1 ), 2 );
  GEOS_ERROR_IF_NE( faceToBlocks.size( 0 ), numUniqueFaces );
  GEOS_ERROR_IF_NE( faceToBlocks.size( 1 ), 2 );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [ numNodes,
                                            uniqueFaceOffsets,
                                            lowestNodeToFaces,
                                            faceToNodes,
                                            faceToCells,
                                            faceToBlocks,
                                            &cellBlocks,
                                            &faceBuilder ]( localIndex const nodeIndex )
  {
    localIndex nodesInFace[ CellBlockManager::maxNodesPerFace() ];
    localIndex curFaceID = uniqueFaceOffsets[nodeIndex];
    arraySlice1d< NodesAndElementOfFace const > const faces = lowestNodeToFaces[ nodeIndex ];
    forEqualRanges( faces.begin(), faces.end(), [&]( auto first, auto last )
    {
      NodesAndElementOfFace const & f0 = *first;
      CellBlock const & cb = cellBlocks.getGroup< CellBlock >( f0.blockIndex );
      localIndex const numNodesInFace = cb.getFaceNodes( f0.cellIndex, f0.faceNumber, nodesInFace );
      GEOS_ASSERT_EQ( numNodesInFace, numNodes );
      GEOS_UNUSED_VAR( numNodes );

      for( localIndex i = 0; i < numNodesInFace; ++i )
      {
        faceToNodes.emplaceBack( curFaceID, nodesInFace[i] );
      }

      faceToCells( curFaceID, 0 ) = f0.cellIndex;
      faceToBlocks( curFaceID, 0 ) = f0.blockIndex;

      if( ++first != last )
      {
        NodesAndElementOfFace const & f1 = *first++;
        faceToCells( curFaceID, 1 ) = f1.cellIndex;
        faceToBlocks( curFaceID, 1 ) = f1.blockIndex;
      }
      else
      {
        faceToCells( curFaceID, 1 ) = -1;
        faceToBlocks( curFaceID, 1 ) = -1;
      }
      GEOS_ASSERT( first == last ); // Should not be more than 2 faces

      ++curFaceID;
    }, faceBuilder.duplicateFaceEquality() );
  } );
}


/**
 * @brief Resize the face maps
 * @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
 * @param [in] uniqueFaceOffsets an containing the unique face IDs for each node in lowestNodeToFaces.
 * @param [out] faceToNodeMap the map from faces to nodes. This function resizes the array appropriately.
 * @param [out] faceToCellMap the map from faces to elements. This function resizes the array appropriately.
 */
void resizeFaceMaps( FaceBuilder const & faceBuilder,
                     arrayView1d< localIndex const > const & uniqueFaceOffsets,
                     ArrayOfArrays< localIndex > & faceToNodeMap,
                     ArrayOfArrays< localIndex > & faceToEdgesMap,
                     array2d< localIndex > & faceToCellMap,
                     array2d< localIndex > & faceToBlockMap )
{
  localIndex const numNodes = faceBuilder.lowestNodeToFaces.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  array1d< localIndex > numNodesPerFace( numUniqueFaces );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [uniqueFaceOffsets,
                                           numNodesPerFace = numNodesPerFace.toView(),
                                           lowestNodeToFaces = faceBuilder.lowestNodeToFaces.toViewConst(),
                                           duplicateFaces = faceBuilder.duplicateFaces.toViewConst(),
                                           &faceBuilder]( localIndex const nodeIndex )
  {
    localIndex curFaceID = uniqueFaceOffsets[ nodeIndex ];
    arraySlice1d< NodesAndElementOfFace const > const faces = lowestNodeToFaces[ nodeIndex ];
    forUniqueValues( faces.begin(), faces.end(), [&]( NodesAndElementOfFace const & f, localIndex )
    {
      numNodesPerFace[ curFaceID++ ] = duplicateFaces.sizeOfArray( f.duplicateFaceNodesIndex ) + CellBlockManager::nodeMapExtraSpacePerFace();
    }, faceBuilder.duplicateFaceEquality() );
  } );

  faceToNodeMap.resizeFromCapacities< parallelHostPolicy >( numUniqueFaces, numNodesPerFace.data() );

  // Each face may belong to _maximum_ 2 elements.
  // If it belongs to only one, we put `-1` for the undefined value.
  faceToCellMap.resize( numUniqueFaces, 2 );
  faceToBlockMap.resize( numUniqueFaces, 2 );

  // TODO I'm really not sure.
  faceToEdgesMap.resize( numUniqueFaces, 2 * CellBlockManager::edgeMapExtraSpacePerFace() );
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
FaceBuilder createLowestNodeToFaces( localIndex const numNodes, const Group & cellBlocks )
{
  array1d< localIndex > faceCounts( numNodes );
  localIndex totalDuplicateFaces = 0;
  localIndex totalDuplicateFaceNodes = 0;
  for( localIndex blockIndex = 0; blockIndex < cellBlocks.numSubGroups(); ++blockIndex )
  {
    CellBlock const & cb = cellBlocks.getGroup< CellBlock >( blockIndex );
    localIndex const numFacesPerElement = cb.numFacesPerElement();
    localIndex const numElements = cb.numElements();

    forAll< parallelHostPolicy >( numElements, [&cb, numFacesPerElement,
                                                counts = faceCounts.toView()]( localIndex const elemID )
    {
      localIndex nodesInFace[ CellBlockManager::maxNodesPerFace() ];
      for( localIndex faceNum = 0; faceNum < numFacesPerElement; ++faceNum )
      {
        // Get all the nodes of the face
        localIndex const numNodesInFace = cb.getFaceNodes( elemID, faceNum, nodesInFace );
        localIndex const lowestNode = *std::min_element( nodesInFace, nodesInFace + numNodesInFace );
        RAJA::atomicInc< parallelHostAtomic >( &counts[ lowestNode ] );
      }
    } );

    totalDuplicateFaces += numFacesPerElement * numElements;
    totalDuplicateFaceNodes += numElements * cb.maxNodesPerFace() * cb.numFacesPerElement();
  }

  FaceBuilder faceBuilder;
  faceBuilder.lowestNodeToFaces.resizeFromCapacities< parallelHostPolicy >( numNodes, faceCounts.data() );

  faceBuilder.duplicateFaces.reserve( totalDuplicateFaces );
  faceBuilder.duplicateFaces.reserveValues( totalDuplicateFaceNodes );

  {

    for( localIndex blockIndex = 0; blockIndex < cellBlocks.numSubGroups(); ++blockIndex )
    {
      CellBlock const & cb = cellBlocks.getGroup< CellBlock >( blockIndex );
      localIndex const numFacesPerElement = cb.numFacesPerElement();
      localIndex const numElements = cb.numElements();

      localIndex const prevFaceOffset = faceBuilder.duplicateFaces.size();
      faceBuilder.duplicateFaces.resize( faceBuilder.duplicateFaces.size() + numFacesPerElement * numElements, cb.maxNodesPerFace() );

      forAll< parallelHostPolicy >( numElements, [&cb, numFacesPerElement, blockIndex, prevFaceOffset,
                                                  lowestNodeToFaces = faceBuilder.lowestNodeToFaces.toView(),
                                                  duplicateFaces = faceBuilder.duplicateFaces.toView()]( localIndex const elemID )
      {
        localIndex nodesInFace[ CellBlockManager::maxNodesPerFace() ];
        for( localIndex faceNum = 0; faceNum < numFacesPerElement; ++faceNum )
        {
          localIndex const duplicateFaceIndex = prevFaceOffset + elemID * numFacesPerElement + faceNum;

          localIndex const numNodesInFace = cb.getFaceNodes( elemID, faceNum, nodesInFace );
          std::sort( nodesInFace, nodesInFace + numNodesInFace );

          duplicateFaces.appendToArray( duplicateFaceIndex, nodesInFace, nodesInFace + numNodesInFace );

          lowestNodeToFaces.emplaceBackAtomic< parallelHostAtomic >( nodesInFace[ 0 ],
                                                                     duplicateFaceIndex,
                                                                     elemID,
                                                                     blockIndex,
                                                                     faceNum );
        }
      } );
    }
  }

  // Loop over all the nodes and sort the associated faces.
  forAll< parallelHostPolicy >( numNodes, [lowestNodeToFaces = faceBuilder.lowestNodeToFaces.toView(),
                                           duplicateFaces = faceBuilder.duplicateFaces.toViewConst()]( localIndex const nodeIndex )
  {
    arraySlice1d< NodesAndElementOfFace > const faces = lowestNodeToFaces[ nodeIndex ];
    std::sort( faces.begin(), faces.end(), [&]( NodesAndElementOfFace const & lhs, NodesAndElementOfFace const & rhs )
    {
      // With C++20 this can all be replaced with std::lexicographical_compare_three_way
      auto const pairOfIters = std::mismatch( duplicateFaces[ lhs.duplicateFaceNodesIndex ].begin(),
                                              duplicateFaces[ lhs.duplicateFaceNodesIndex ].end(),
                                              duplicateFaces[ rhs.duplicateFaceNodesIndex ].begin(),
                                              duplicateFaces[ rhs.duplicateFaceNodesIndex ].end() );

      // If the ranges are equal
      if( pairOfIters.first == duplicateFaces[ lhs.duplicateFaceNodesIndex ].end() &&
          pairOfIters.second == duplicateFaces[ rhs.duplicateFaceNodesIndex ].end() )
      {
        return std::tie( lhs.blockIndex, lhs.cellIndex ) < std::tie( rhs.blockIndex, rhs.cellIndex );
      }

      // If the second range is a prefix of the first, the first is greater than the secondd.
      if( pairOfIters.first == duplicateFaces[ lhs.duplicateFaceNodesIndex ].end() )
      {
        return false;
      }

      // If the first range is a prefix of the second, the first is less than the second.
      if( pairOfIters.second == duplicateFaces[ rhs.duplicateFaceNodesIndex ].end() )
      {
        return true;
      }

      // Otherwise simply compare the first non-equal values.
      return *pairOfIters.first < *pairOfIters.second;
    } );
  } );

  return faceBuilder;
}


/**
 * @brief Filling the elements to faces maps in the cell blocks.
 * @param lowestNodeToFaces The lowest node to faces information array.
 * @param uniqueFaceOffsets The unique face offsets.
 * @param cellBlocks The cell blocks for which we need to compute the element to faces mappings.
 *
 * @note @p lowestNodeToFaces and @p uniqueFaceOffsets are better described in the documentations of the functions that build them.
 */
void fillElementToFacesOfCellBlocks( FaceBuilder const & faceBuilder,
                                     arrayView1d< localIndex const > const & uniqueFaceOffsets,
                                     Group & cellBlocks )
{
  ArrayOfArraysView< NodesAndElementOfFace const > const & lowestNodeToFaces
    = faceBuilder.lowestNodeToFaces.toViewConst();

  localIndex const numNodes = lowestNodeToFaces.size();
  forAll< parallelHostPolicy >( numNodes, [lowestNodeToFaces,
                                           uniqueFaceOffsets,
                                           &cellBlocks,
                                           &faceBuilder]( localIndex const nodeIndex )
  {
    arraySlice1d< NodesAndElementOfFace const > const faces = lowestNodeToFaces[ nodeIndex ];
    localIndex curFaceID = uniqueFaceOffsets[nodeIndex];

    forEqualRanges( faces.begin(), faces.end(), [&]( auto first, auto last )
    {
      while( first != last )
      {
        NodesAndElementOfFace const & f = *first++;
        CellBlock & cb0 = cellBlocks.getGroup< CellBlock >( f.blockIndex );
        cb0.setElementToFaces( f.cellIndex, f.faceNumber, curFaceID );
      }
      ++curFaceID;
    }, faceBuilder.duplicateFaceEquality() );
  } );
}

/**
 * @brief Fills the element to edges mappings of all the cells provided through @p cellBlocks.
 * @param faceToEdges We need the face to edges mapping to get some edge index.
 * @param cellBlocks The cell blocks for which the mappings will be constructed.
 */
void fillElementToEdgesOfCellBlocks( ArrayOfArraysView< localIndex const > const & faceToEdges,
                                     Group & cellBlocks )
{
  for( localIndex blockIndex = 0; blockIndex < cellBlocks.numSubGroups(); ++blockIndex )
  {
    CellBlock & cellBlock = cellBlocks.getGroup< CellBlock >( blockIndex );
    arrayView2d< localIndex const > const cellToFaces = cellBlock.getElemToFacesConstView();

    // We build the edges of each face of each cell,
    // so we can construct the cell to edges mapping.
    // Some specific care is required not to insert edges twice (faces share edges).
    // Another implementation (used in other contexts) would use some edge signature
    // to remove the duplicates.

    // Loop over the cells
    forAll< parallelHostPolicy >( cellBlock.numElements(), [cellToFaces,
                                                            faceToEdges,
                                                            &cellBlock]( localIndex const kc )
    {
      int count = 0;
      for( localIndex kf = 0; kf < cellBlock.numFacesPerElement(); ++kf )
      {
        // Loop over edges of each face
        localIndex const faceIndex = cellToFaces[kc][kf];
        for( localIndex ke = 0; ke < faceToEdges.sizeOfArray( faceIndex ); ++ke )
        {
          bool isUnique = true;
          localIndex const edgeIndex = faceToEdges[faceIndex][ke];

          // Loop over edges that have already been added to the element.
          for( localIndex kec = 0; kec < count + 1; ++kec )
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
    } ); // end cell loop
  }
}

void CellBlockManager::buildFaceMaps()
{
  GEOS_MARK_FUNCTION;

  FaceBuilder const faceBuilder =
    createLowestNodeToFaces( m_numNodes, this->getCellBlocks() );

  array1d< localIndex > const uniqueFaceOffsets =
    computeUniqueValueOffsets< parallelHostPolicy >( faceBuilder.lowestNodeToFaces.toViewConst(),
                                                     faceBuilder.duplicateFaceEquality() );


  m_numFaces = uniqueFaceOffsets.back();

  resizeFaceMaps( faceBuilder,
                  uniqueFaceOffsets,
                  m_faceToNodes,
                  m_faceToEdges,
                  m_faceToCells.toCellIndex,
                  m_faceToCells.toBlockIndex );

  populateFaceMaps( getCellBlocks(),
                    faceBuilder,
                    uniqueFaceOffsets,
                    m_faceToNodes.toView(),
                    m_faceToCells.toCellIndex,
                    m_faceToCells.toBlockIndex );

  fillElementToFacesOfCellBlocks( faceBuilder,
                                  uniqueFaceOffsets,
                                  this->getCellBlocks() );
}

void CellBlockManager::buildNodeToEdges()
{
  m_nodeToEdges = meshMapUtilities::transposeIndexMap< parallelHostPolicy >( m_edgeToNodes.toViewConst(),
                                                                             m_numNodes,
                                                                             edgeMapExtraSpacePerNode() );
}

void CellBlockManager::buildMaps()
{
  GEOS_MARK_FUNCTION;

  buildFaceMaps();
  m_numEdges = buildEdgeMaps( m_numNodes,
                              m_faceToNodes.toViewConst(),
                              m_faceToEdges,
                              m_edgeToFaces,
                              m_edgeToNodes );
  buildNodeToEdges();

  fillElementToEdgesOfCellBlocks( m_faceToEdges.toViewConst(), this->getCellBlocks() );
}

ArrayOfArrays< localIndex > CellBlockManager::getFaceToNodes() const
{
  return m_faceToNodes;
}

ArrayOfArrays< localIndex > CellBlockManager::getNodeToFaces() const
{
  return meshMapUtilities::transposeIndexMap< parallelHostPolicy >( m_faceToNodes.toViewConst(),
                                                                    m_numNodes,
                                                                    faceMapExtraSpacePerNode() );
}

ToCellRelation< array2d< localIndex > > CellBlockManager::getFaceToElements() const
{
  return m_faceToCells;
}

const Group & CellBlockManager::getCellBlocks() const
{
  return this->getGroup( viewKeyStruct::cellBlocks() );
}

Group & CellBlockManager::getCellBlocks()
{
  return this->getGroup( viewKeyStruct::cellBlocks() );
}

Group const & CellBlockManager::getFaceBlocks() const
{
  return this->getGroup( viewKeyStruct::faceBlocks() );
}

Group & CellBlockManager::getFaceBlocks()
{
  return this->getGroup( viewKeyStruct::faceBlocks() );
}

Group & CellBlockManager::getLineBlocks()
{
  return this->getGroup( viewKeyStruct::lineBlocks() );
}

LineBlockABC const & CellBlockManager::getLineBlock( string name ) const
{
  return this->getGroup( viewKeyStruct::lineBlocks() ).getGroup< LineBlockABC >( name );
}

localIndex CellBlockManager::numNodes() const
{
  return m_numNodes;
}

localIndex CellBlockManager::numCellBlocks() const
{
  return this->getCellBlocks().numSubGroups();
}

CellBlock const & CellBlockManager::getCellBlock( localIndex const blockIndex ) const
{
  return this->getCellBlocks().getGroup< CellBlock >( blockIndex );
}

localIndex CellBlockManager::numFaces() const
{
  return m_numFaces;
}

ArrayOfArrays< localIndex > CellBlockManager::getEdgeToFaces() const
{
  return m_edgeToFaces;
}

array2d< localIndex > CellBlockManager::getEdgeToNodes() const
{
  return m_edgeToNodes;
}

ArrayOfArrays< localIndex > CellBlockManager::getFaceToEdges() const
{
  return m_faceToEdges;
}

ArrayOfArrays< localIndex > CellBlockManager::getNodeToEdges() const
{
  return meshMapUtilities::transposeIndexMap< parallelHostPolicy >( m_edgeToNodes.toViewConst(),
                                                                    m_numNodes,
                                                                    edgeMapExtraSpacePerNode() );
}

localIndex CellBlockManager::numEdges() const
{
  return m_numEdges;
}

CellBlock & CellBlockManager::registerCellBlock( string const & name )
{
  return this->getCellBlocks().registerGroup< CellBlock >( name );
}

FaceBlock & CellBlockManager::registerFaceBlock( string const & name )
{
  return this->getFaceBlocks().registerGroup< FaceBlock >( name );
}

LineBlock & CellBlockManager::registerLineBlock( string const & name )
{
  return this->getLineBlocks().registerGroup< LineBlock >( name );
}

array2d< real64, nodes::REFERENCE_POSITION_PERM > CellBlockManager::getNodePositions() const
{
  return m_nodesPositions;
}

arrayView2d< real64, nodes::REFERENCE_POSITION_USD > CellBlockManager::getNodePositions()
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

static array1d< real64 > gaussLobattoPoints( int order )
{
  array1d< real64 > GaussLobattoPts( order+1 );

  switch( order )
  {
    case 2:
      GaussLobattoPts[0] = 0.0;
      break;
    case 3:
      static constexpr real64 sqrt5 = 2.2360679774997897;
      GaussLobattoPts[0] = -1./sqrt5;
      GaussLobattoPts[1] = 1./sqrt5;
      break;
    case 4:
      static constexpr real64 sqrt3_7 = 0.6546536707079771;
      GaussLobattoPts[0] = -sqrt3_7;
      GaussLobattoPts[1] = 0.0;
      GaussLobattoPts[2] = sqrt3_7;
      break;
    case 5:
      static constexpr real64 sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
      static constexpr real64 sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
      static constexpr real64 sqrt_inv21 = 0.218217890235992381;
      GaussLobattoPts[0] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
      GaussLobattoPts[1] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
      GaussLobattoPts[2] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
      GaussLobattoPts[3] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
      break;
  }
  return GaussLobattoPts;
}

inline void trilinearInterp( real64 const alpha,
                             real64 const beta,
                             real64 const gamma,
                             real64 const (&X)[8][3],
                             real64 (& coords)[3] )
{
  for( int i=0; i<3; i++ )
  {
    coords[i] = X[0][i]*( 1.0-alpha )*( 1.0-beta )*( 1.0-gamma )+
                X[1][i]*    alpha    *( 1.0-beta )*( 1.0-gamma )+
                X[2][i]*( 1.0-alpha )*    beta    *( 1.0-gamma )+
                X[3][i]*    alpha    *    beta    *( 1.0-gamma )+
                X[4][i]*( 1.0-alpha )*( 1.0-beta )*  gamma+
                X[5][i]*    alpha    *( 1.0-beta )*  gamma+
                X[6][i]*( 1.0-alpha )*    beta    *  gamma+
                X[7][i]*    alpha    *    beta    *  gamma;
  }
}

void CellBlockManager::generateHighOrderMaps( localIndex const order,
                                              globalIndex const maxVertexGlobalID,
                                              globalIndex const maxEdgeGlobalID,
                                              globalIndex const maxFaceGlobalID,
                                              arrayView1d< globalIndex const > const edgeLocalToGlobal,
                                              arrayView1d< globalIndex const > const faceLocalToGlobal )
{

  GEOS_MARK_FUNCTION;

  // constants for hex mesh
  localIndex const numEdgesPerFace = 4;
  localIndex const numFacesPerCell = 6;
  localIndex const numEdgesPerCell = 12;
  localIndex const numVerticesPerFace = 4;
  localIndex const numVerticesPerCell = 8;

  localIndex const numNodesPerEdge = ( order+1 );
  localIndex const numNodesPerFace = ( order+1 )*( order+1 );
  localIndex const numNodesPerCell = ( order+1 )*( order+1 )*( order+1 );

  localIndex const numInternalNodesPerEdge = ( order-1 );
  localIndex const numInternalNodesPerFace = ( order-1 )*( order-1 );
  localIndex const numInternalNodesPerCell = ( order-1 )*( order-1 )*( order-1 );

  localIndex const numLocalVertices = this->numNodes();
  localIndex const numLocalEdges = this->numEdges();
  localIndex const numLocalFaces = this->numFaces();

  GEOS_LOG_RANK_0("TEST TEST: numLocalEdges = " << numLocalEdges);
  GEOS_LOG_RANK_0("TEST TEST: numLocalFaces = " << numLocalFaces);
  GEOS_LOG_RANK_0("TEST TEST: numLocalVertices = " << numLocalVertices);

  localIndex numLocalCells = 0;
  this->getCellBlocks().forSubGroups< CellBlock >( [&]( CellBlock & cellBlock )
  {
    numLocalCells += cellBlock.numElements();
  } );

  ////////////////////////////////
  // Get the new number of nodes
  ////////////////////////////////

  localIndex numLocalNodes = numLocalVertices
                             + numLocalEdges * numInternalNodesPerEdge
                             + numLocalFaces * numInternalNodesPerFace
                             + numLocalCells * numInternalNodesPerCell;

  array1d< globalIndex > const nodeLocalToGlobalSource ( m_nodeLocalToGlobal );
  array2d< localIndex > const edgeToNodesMapSource( m_edgeToNodes );
  ArrayOfArrays< localIndex > const faceToNodesMapSource( m_faceToNodes );
  array2d< real64, nodes::REFERENCE_POSITION_PERM > const refPosSource ( m_nodesPositions );

  m_numNodes = numLocalNodes;
  m_nodeLocalToGlobal.resize( m_numNodes );
  m_edgeToNodes.resize( m_numEdges, order+1 );
  m_nodesPositions.resize( m_numNodes );

  // also assign node coordinates using trilinear interpolation in th elements
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > refPosNew = this->getNodePositions();
  refPosNew.setValues< parallelHostPolicy >( -1.0 );

  array1d< real64 > glCoords = gaussLobattoPoints( order );

  // ---------------------------------------
  // initialize the node local to global map
  // ---------------------------------------

  // The general idea for building nodes local-to-global map for the high-order meshes are:
  // 1. for all the nodes that already in the base mesh-level, we keep their globalIDs the same
  // 2. for the new nodes on the edges, we assign their global ID according to the global information from edgeLocalToGlobal map
  // 3. for the new nodes on the faces, we assign their global ID according to the global information from faceLocalToGlobal map
  // 4. for the new nodes on the internal of each elements, we assign the global ID for the new nodes
  //    according to the global information for the elements from elementLocalToGlobal map.
  arrayView1d< globalIndex > nodeLocalToGlobalNew = m_nodeLocalToGlobal.toView();

  // nodeIDs is a hash map that contains unique vertices and their indices for shared nodes:
  // - A vertex is identified by 6 integers [i1 i2 i3 i4 a b], as follows:
  // - nodes on a vertex v are identified by the vector [v -1 -1 -1 -1 -1]
  // - nodes on a edge are given by a linear interpolation between vertices v1 and v2, 'a' steps away from v1
  //   (we assume that v1 < v2 and identify these nodes with [v1 v2 -1 -1 a -1]).
  // - nodes on a face are given by a bilinear interpolation between edges v1-v2 and v3-v4
  //   (v1-v4 and v2-v3 are the diagonals), with interpolation parameters 'a' and 'b'.
  //   (we assume that v1 is the smallest, and that v2 < v3) Then these nodes are identified with [v1 v2 v3 v4 a b]
  // - nodes within the internal of a cell are encountered only once, and thus do not need to be put in the hash map

  // Create new nodes, with local and global IDs
  localIndex localNodeOffset = 0;
  
  GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Vertices --- 11111");
  forAll< parallelHostPolicy >( numLocalVertices, 
                                [ nodeLocalToGlobalSource=nodeLocalToGlobalSource.toView(), 
                                  nodeLocalToGlobalNew=nodeLocalToGlobalNew.toView(), 
                                  refPosSource=refPosSource.toView(),
                                  refPosNew=refPosNew.toView() ]( localIndex const iter_vertex)
  {
    nodeLocalToGlobalNew[ iter_vertex ] = nodeLocalToGlobalSource[ iter_vertex ];
    refPosNew[ iter_vertex ][0] = refPosSource[ iter_vertex ][0];
    refPosNew[ iter_vertex ][1] = refPosSource[ iter_vertex ][1];
    refPosNew[ iter_vertex ][2] = refPosSource[ iter_vertex ][2];
  });
  localNodeOffset = numLocalVertices;

  GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Vertices --- 11111");

  //GEOS_LOG_RANK_0("TEST TEST: nodeLocalToGlobalSource = " << nodeLocalToGlobalSource);
  //GEOS_LOG_RANK_0("TEST TEST: nodeLocalToGlobalNew = " << nodeLocalToGlobalNew);
  
  //GEOS_LOG_RANK_0("TEST TEST: refPosSource = " << refPosSource);
  //GEOS_LOG_RANK_0("TEST TEST: refPosNew = " << refPosNew);
  /*
  std::unordered_map< std::array< localIndex, 6 >, localIndex, NodeKeyHasher< localIndex > > nodeIDs;
  localIndex localNodeID = 0;

  GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Vertices --- 22222");
  for( localIndex iter_vertex=0; iter_vertex < numLocalVertices; iter_vertex++ )
  {
    nodeLocalToGlobalNew[ iter_vertex ] = nodeLocalToGlobalSource[ iter_vertex ];
    nodeIDs[ createNodeKey( iter_vertex ) ] = localNodeID;
    localNodeID++;
  }
  GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Vertices --- 22222");
  */
  
  /*
  forAll< parallelHostPolicy >( numLocalVertices, 
                                [ nodeLocalToGlobalSource=nodeLocalToGlobalSource.toView(), &mutexNodeID, 
                                  &nodeLocalToGlobalNew, &nodeIDs, &localNodeID]( localIndex const iter_vertex)
  {
    nodeLocalToGlobalNew[ localNodeID ] = nodeLocalToGlobalSource[ iter_vertex ];
    mutexNodeID.lock();
    nodeIDs[ createNodeKey( iter_vertex ) ] = localNodeID;
    mutexNodeID.unlock();
    localNodeID++;
  } );
  //}

  */

  //////////////////////////
  // Edges
  //////////////////////////

  // -------------------------------------
  // ---- initialize edge-to-node map ----
  // -------------------------------------

  arrayView2d< localIndex > edgeToNodeMapNew = m_edgeToNodes.toView();
  // create / retrieve nodes on edges
  localIndex globalNodeOffset = maxVertexGlobalID;

  GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Edges -- 111111");
  forAll< parallelHostPolicy >( numLocalEdges, 
                                [ edgeToNodesMapSource=edgeToNodesMapSource.toView(),
                                  edgeToNodeMapNew=edgeToNodeMapNew.toView(),
                                  refPosSrc=refPosSource.toView(),
                                  refPosNew=refPosNew.toView(),
                                  edgeLocalToGlobal=edgeLocalToGlobal.toView(),
                                  nodeLocalToGlobalNew=nodeLocalToGlobalNew.toView(),
                                  numInternalNodesPerEdge, globalNodeOffset, glCoords, localNodeOffset]( localIndex const iter_edge)
  {
      localIndex edgeHeadNode = edgeToNodesMapSource[ iter_edge ][ 0 ];
      localIndex edgeEndNode = edgeToNodesMapSource[ iter_edge ][ 1 ];

      edgeToNodeMapNew[ iter_edge ][ 0 ] = edgeHeadNode;
      edgeToNodeMapNew[ iter_edge ][ numInternalNodesPerEdge + 1] = edgeEndNode;

      for( localIndex iter_node = 0; iter_node < numInternalNodesPerEdge; iter_node++ )
      {
        real64 alpha = ( glCoords[ iter_node % numInternalNodesPerEdge ] + 1.0 ) / 2.0;
        localIndex nodeLocalID = localNodeOffset + iter_edge * ( numInternalNodesPerEdge ) + iter_node; 
        edgeToNodeMapNew[ iter_edge ][ iter_node + 1 ] = nodeLocalID;
        nodeLocalToGlobalNew[ nodeLocalID ] = globalNodeOffset + edgeLocalToGlobal[ iter_edge ] * numInternalNodesPerEdge + iter_node;

        refPosNew[nodeLocalID][0] = refPosSrc[edgeEndNode][0] * alpha + refPosSrc[edgeHeadNode][0] * ( 1.0-alpha ) ;
        refPosNew[nodeLocalID][1] = refPosSrc[edgeEndNode][1] * alpha + refPosSrc[edgeHeadNode][1] * ( 1.0-alpha ) ;
        refPosNew[nodeLocalID][2] = refPosSrc[edgeEndNode][2] * alpha + refPosSrc[edgeHeadNode][2] * ( 1.0-alpha ) ;

      }
  });
  localNodeOffset += numInternalNodesPerEdge * numLocalEdges;
  //GEOS_LOG_RANK_0("TEST TEST: refPosNew = " << refPosNew); 
  GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Edges -- 111111");

  //GEOS_LOG_RANK_0("TEST TEST: edgeToNodesMapSource = " << edgeToNodesMapSource);
  //GEOS_LOG_RANK_0("TEST TEST: edgeToNodeMapNew = " << edgeToNodeMapNew);

  //GEOS_LOG_RANK_0("TEST TEST: refPosSource = " << refPosSource);
  //GEOS_LOG_RANK_0("TEST TEST: refPosNew = " << refPosNew);

  //GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Edges -- 222222");
  //for( localIndex iter_edge = 0; iter_edge < numLocalEdges; iter_edge++ )
  //RAJA::omp::mutex mutexNodeID;
  /*
  forAll< parallelHostPolicy >( numLocalEdges, 
                                [ nodeLocalToGlobalSource=nodeLocalToGlobalSource.toView(),
                                  edgeToNodesMapSource=edgeToNodesMapSource.toView(),
                                  edgeLocalToGlobal=edgeLocalToGlobal.toView(),
                                  numNodesPerEdge, order, globalNodeOffset, numInternalNodesPerEdge, &mutexNodeID,
                                  &edgeToNodeMapNew,&nodeLocalToGlobalNew, &nodeIDs, &localNodeID]( localIndex const iter_edge)
  */
  /*
  {
    localIndex v1 = edgeToNodesMapSource[ iter_edge ][ 0 ];
    localIndex v2 = edgeToNodesMapSource[ iter_edge ][ 1 ];
    globalIndex gv1 = nodeLocalToGlobalSource[v1];
    globalIndex gv2 = nodeLocalToGlobalSource[v2];

    for( int q=0; q<numNodesPerEdge; q++ )
    {
      localIndex nodeID;
      std::array< localIndex, 6 > nodeKey = createNodeKey( v1, v2, q, order );
      if( nodeIDs.count( nodeKey ) == 0 )
      {
        // this is an internal edge node: create it
        nodeID = localNodeID;
        mutexNodeID.lock();
        nodeIDs[ nodeKey ] = nodeID;
        mutexNodeID.unlock();
        std::array< globalIndex, 6 > referenceOrientation = createNodeKey( gv1, gv2, q, order );
        int gq = referenceOrientation[4] - 1;
        nodeLocalToGlobalNew[ nodeID ] = globalNodeOffset + edgeLocalToGlobal[ iter_edge ] * numInternalNodesPerEdge + gq;
        localNodeID++;
      }
      else
      {
        nodeID = nodeIDs[ nodeKey ];
      }
      edgeToNodeMapNew[ iter_edge ][ q ] = nodeID;
    }
  }
  //} );
  GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Edges -- 222222");
  */

  //GEOS_LOG_RANK_0("TEST TEST: edgeToNodeMapNew = " << edgeToNodeMapNew);
  //GEOS_LOG_RANK_0("TEST TEST: nodeLocalToGlobalSource = " << nodeLocalToGlobalSource);
  //GEOS_LOG_RANK_0("TEST TEST: nodeLocalToGlobalNew = " << nodeLocalToGlobalNew);

  /////////////////////////
  // Faces
  //////////////////////////

  /*
  // initialize faceToNodeMap for the high-order mesh-level
  ArrayOfArrays< localIndex > & faceToNodeMapNew = m_faceToNodes;
  // number of elements in each row of the map as capacity
  array1d< localIndex > counts( faceToNodeMapNew.size());
  counts.setValues< parallelHostPolicy >( numNodesPerFace );
  //  reconstructs the faceToNodeMap with the provided capacity in counts
  faceToNodeMapNew.resizeFromCapacities< parallelHostPolicy >( faceToNodeMapNew.size(), counts.data() );
  // setup initial values of the faceToNodeMap using emplaceBack
  forAll< parallelHostPolicy >( faceToNodeMapNew.size(),
                                [ faceToNodeMapNew = faceToNodeMapNew.toView() ]( localIndex const faceIndex )
  {
    for( localIndex i = 0; i < faceToNodeMapNew.capacityOfArray( faceIndex ); ++i )
    {
      faceToNodeMapNew.emplaceBack( faceIndex, -1 );
    }
  } );
  */

  // initialize faceToNodeMap for the high-order mesh-level
  ArrayOfArrays< localIndex > faceToNodeMapNew2 = m_faceToNodes;
  // number of elements in each row of the map as capacity
  array1d< localIndex > counts( faceToNodeMapNew2.size());
  counts.setValues< parallelHostPolicy >( numNodesPerFace );
  //  reconstructs the faceToNodeMap with the provided capacity in counts
  faceToNodeMapNew2.resizeFromCapacities< parallelHostPolicy >( faceToNodeMapNew2.size(), counts.data() );
  // setup initial values of the faceToNodeMap using emplaceBack
  forAll< parallelHostPolicy >( faceToNodeMapNew2.size(),
                                [ faceToNodeMapNew2 = faceToNodeMapNew2.toView() ]( localIndex const faceIndex )
  {
    for( localIndex i = 0; i < faceToNodeMapNew2.capacityOfArray( faceIndex ); ++i )
    {
      faceToNodeMapNew2.emplaceBack( faceIndex, -1 );
    }
  } );

  //GEOS_LOG_RANK_0("TEST TEST: faceToNodesMapSource = " << faceToNodesMapSource);
  globalNodeOffset = maxVertexGlobalID + maxEdgeGlobalID * numInternalNodesPerEdge;

  GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Faces -- 111111");
  //for( localIndex iter_edge = 0; iter_edge < numLocalEdges; iter_edge++ )
  forAll< parallelHostPolicy >( numLocalFaces, 
                                [ faceToNodesMapSource=faceToNodesMapSource.toView(),
                                  edgeToNodeMapNew=edgeToNodeMapNew.toView(),
                                  faceToNodeMapNew=faceToNodeMapNew2.toView(),
                                  m_faceToEdges=m_faceToEdges.toView(),
                                  refPosSrc=refPosSource.toView(),
                                  refPosNew=refPosNew.toView(),
                                  faceLocalToGlobal=faceLocalToGlobal.toView(),
                                  nodeLocalToGlobalNew=nodeLocalToGlobalNew.toView(),
                                  numNodesPerEdge, numEdgesPerFace, numVerticesPerFace,
                                  numInternalNodesPerEdge, numInternalNodesPerFace,
                                  globalNodeOffset, glCoords, localNodeOffset ]( localIndex const iter_face)
  {
      localIndex faceVertID[ numVerticesPerFace ];
      for( localIndex iter_node=0; iter_node<numVerticesPerFace; iter_node++ )
      {
        faceVertID[ iter_node ] = faceToNodesMapSource[ iter_face ][ iter_node ];
        faceToNodeMapNew[ iter_face ][ iter_node ] = faceVertID[ iter_node ];
      }

      for( localIndex iter_edge=0; iter_edge<numEdgesPerFace; iter_edge++ )
        for( localIndex iter_node=0; iter_node<numInternalNodesPerEdge; iter_node++ )
          faceToNodeMapNew[ iter_face ][ numVerticesPerFace + iter_edge * numInternalNodesPerEdge + iter_node ] 
            = edgeToNodeMapNew[m_faceToEdges[iter_face][iter_edge]][iter_node+1];

      for( localIndex iter_node=0; iter_node<numInternalNodesPerFace; iter_node++ )
      {
        real64 alpha = ( glCoords[ iter_node % numInternalNodesPerEdge ] + 1.0 ) / 2.0;
        real64 beta = ( glCoords[  ( iter_node / numInternalNodesPerEdge ) % numInternalNodesPerEdge ] + 1.0 ) / 2.0;

        localIndex nodeLocalID = localNodeOffset + iter_face * numInternalNodesPerFace + iter_node;
        faceToNodeMapNew[ iter_face ][ numVerticesPerFace + numEdgesPerFace * numInternalNodesPerEdge + iter_node ] = nodeLocalID;
        nodeLocalToGlobalNew[ nodeLocalID ] = globalNodeOffset + faceLocalToGlobal[ iter_face ] * numInternalNodesPerFace + iter_node;

        refPosNew[nodeLocalID][0] = refPosSrc[faceVertID[0]][0]*( 1.0-alpha )*( 1.0-beta ) + refPosSrc[faceVertID[1]][0]*alpha*( 1.0-beta )+
                                    refPosSrc[faceVertID[2]][0]*( 1.0-alpha )*beta + refPosSrc[faceVertID[3]][0]*alpha*beta;

        refPosNew[nodeLocalID][1] = refPosSrc[faceVertID[0]][1]*( 1.0-alpha )*( 1.0-beta ) + refPosSrc[faceVertID[1]][1]*alpha*( 1.0-beta )+
                                    refPosSrc[faceVertID[2]][1]*( 1.0-alpha )*beta + refPosSrc[faceVertID[3]][1]*alpha*beta;

        refPosNew[nodeLocalID][2] = refPosSrc[faceVertID[0]][2]*( 1.0-alpha )*( 1.0-beta ) + refPosSrc[faceVertID[1]][2]*alpha*( 1.0-beta )+
                                    refPosSrc[faceVertID[2]][2]*( 1.0-alpha )*beta + refPosSrc[faceVertID[3]][2]*alpha*beta;
      }

  });

  localNodeOffset += numInternalNodesPerFace * numLocalFaces ;
  //}
  GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Faces -- 111111");
  //GEOS_LOG_RANK_0("TEST TEST: faceToNodeMapNew2 = " << faceToNodeMapNew2);
  //GEOS_LOG_RANK_0("TEST TEST: refPosNew = " << refPosNew); 


  // create / retrieve nodes on faces
  // by default, faces are oriented so that the normal has an outward orientation for the current rank.
  // For this reason :
  // - The 3rd and 4th node need to be swapped, to be consistent with the GL ordering
  // - The global IDs of the internal nodes must be referred to a global "reference" orientation (using the createNodeKey method)

  /*
  GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Faces -- 222222");
  for( localIndex iter_face = 0; iter_face < numLocalFaces; iter_face++ )
  {
    localIndex v1 = faceToNodesMapSource[ iter_face ][ 0 ];
    localIndex v2 = faceToNodesMapSource[ iter_face ][ 1 ];
    localIndex v3 = faceToNodesMapSource[ iter_face ][ 2 ];
    localIndex v4 = faceToNodesMapSource[ iter_face ][ 3 ];
    std::swap( v3, v4 );
    globalIndex gv1 = nodeLocalToGlobalSource.toView()[v1];
    globalIndex gv2 = nodeLocalToGlobalSource.toView()[v2];
    globalIndex gv3 = nodeLocalToGlobalSource.toView()[v3];
    globalIndex gv4 = nodeLocalToGlobalSource.toView()[v4];
    for( int q1=0; q1<numNodesPerEdge; q1++ )
    {
      for( int q2=0; q2<numNodesPerEdge; q2++ )
      {
        localIndex nodeID;
        std::array< localIndex, 6 > nodeKey = createNodeKey( v1, v2, v3, v4, q1, q2, order );
        if( nodeIDs.count( nodeKey ) == 0 )
        {
          // this is an internal face node: create it
          nodeID = localNodeID;
          nodeIDs[ nodeKey ] = nodeID;
          std::array< globalIndex, 6 > referenceOrientation = createNodeKey( gv1, gv2, gv3, gv4, q1, q2, order );
          int gq1 = referenceOrientation[4] - 1;
          int gq2 = referenceOrientation[5] - 1;
          nodeLocalToGlobalNew[ nodeID ] = globalNodeOffset + faceLocalToGlobal[ iter_face ] * numInternalNodesPerFace + gq1* numInternalNodesPerEdge + gq2;
          localNodeID++;
        }
        else
        {
          nodeID = nodeIDs[ nodeKey ];
        }
        faceToNodeMapNew[ iter_face ][ q2 + q1*numNodesPerEdge ] = nodeID;
      }
    }
  }
  GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Faces -- 222222");
  // GEOS_LOG_RANK_0("TEST TEST: faceToNodeMapNew = " << faceToNodeMapNew);
  */

  // add all nodes to the target set "all"
  SortedArray< localIndex > & allNodesSet = this->getNodeSets()[ "all" ];
  allNodesSet.reserve( numLocalNodes );

  for( localIndex iter_nodes=0; iter_nodes< numLocalNodes; ++iter_nodes )
  {
    allNodesSet.insert( iter_nodes );
  }

  /////////////////////////
  // Elements
  //////////////////////////

  globalNodeOffset = maxVertexGlobalID + maxEdgeGlobalID * numInternalNodesPerEdge + maxFaceGlobalID * numInternalNodesPerFace;
  //std::array< localIndex, 6 > const nullKey = std::array< localIndex, 6 >{ -1, -1, -1, -1, -1, -1 };

  // initialize the elements-to-nodes map
  arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodesNew;

  this->getCellBlocks().forSubGroups< CellBlock >( [&]( CellBlock & cellBlock )
  {
    arrayView1d< globalIndex > elementLocalToGlobal( cellBlock.localToGlobalMap() );
    array2d< localIndex, cells::NODE_MAP_PERMUTATION > elemsToNodesSource ( cellBlock.getElemToNodes() );

    array2d< localIndex > elemsToEdges ( cellBlock.getElemToEdges() );
    array2d< localIndex > elemsToFaces ( cellBlock.getElemToFaces() );

    cellBlock.resizeNumNodes( numNodesPerCell );
    elemsToNodesNew = cellBlock.getElemToNode();
    localIndex const numCellElements = cellBlock.numElements();

    // then loop through all the elements and assign the globalID according to the globalID of the Element
    // and insert the new local to global ID ( for the internal nodes of elements ) into the nodeLocalToGlobal
    // retrieve finite element type

    //GEOS_LOG_RANK_0("TEST TEST: elemsToNodesSource = " << elemsToNodesSource);
    //GEOS_LOG_RANK_0("TEST TEST: elemsToNodesNew = " << elemsToNodesNew);
    //GEOS_LOG_RANK_0("TEST TEST: elemsToEdges = " << elemsToEdges);
    //GEOS_LOG_RANK_0("TEST TEST: elemsToFaces = " << elemsToFaces);

    GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Elements -- 111111 ");
    //GEOS_LOG_RANK_0("TEST TEST: faceToNodeMapNew2 = " << faceToNodeMapNew2);

    forAll< RAJA::omp_parallel_for_exec >( numCellElements, 
                                  [ elemsToNodesSource=elemsToNodesSource.toView(), 
                                    elemsToNodesNew=elemsToNodesNew.toView(),
                                    elemsToEdges=elemsToEdges.toView(), 
                                    elemsToFaces=elemsToFaces.toView(), 
                                    edgeToNodeMapNew=edgeToNodeMapNew.toView(), 
                                    faceToNodeMapNew2=faceToNodeMapNew2.toView(), 
                                    refPosSrc=refPosSource.toView(),
                                    refPosNew=refPosNew.toView(),
                                    elementLocalToGlobal=elementLocalToGlobal.toView(),
                                    nodeLocalToGlobalNew=nodeLocalToGlobalNew.toView(),
                                    numEdgesPerFace, numVerticesPerFace,
                                    numVerticesPerCell, numEdgesPerCell, numFacesPerCell, numNodesPerCell,
                                    numInternalNodesPerCell, numInternalNodesPerEdge, numInternalNodesPerFace,
                                    globalNodeOffset, glCoords, localNodeOffset ]( localIndex const iter_elem )
    {
      localIndex elemVertID[ numVerticesPerCell];
      for( localIndex iter_node=0; iter_node<numVerticesPerCell; iter_node++ )
      {
        elemVertID[ iter_node ] = elemsToNodesSource[ iter_elem ][ iter_node ];
        elemsToNodesNew[ iter_elem ][ iter_node ] = elemVertID[ iter_node ];
      }

      for( localIndex iter_edge=0; iter_edge<numEdgesPerCell; iter_edge++ )
        for( localIndex iter_node=0; iter_node<numInternalNodesPerEdge; iter_node++ )
          elemsToNodesNew[ iter_elem ][ numVerticesPerCell + iter_edge * numInternalNodesPerEdge + iter_node ] 
            = edgeToNodeMapNew[elemsToEdges[iter_elem][iter_edge]][iter_node+1];

      for( localIndex iter_face=0; iter_face<numFacesPerCell; iter_face++ )
        for( localIndex iter_node=0; iter_node<numInternalNodesPerFace; iter_node++ )
        {
          elemsToNodesNew[ iter_elem ][ numVerticesPerCell + numEdgesPerCell * numInternalNodesPerEdge + iter_face * numInternalNodesPerFace + iter_node ] 
            = faceToNodeMapNew2[elemsToFaces[iter_elem][iter_face]][numVerticesPerFace + numEdgesPerFace * numInternalNodesPerEdge + iter_node ];
        }

      for( localIndex iter_node=0; iter_node<numInternalNodesPerCell; iter_node++ )
      {
        real64 alpha = ( glCoords[ iter_node % numInternalNodesPerEdge ] + 1.0 ) / 2.0;
        real64 beta = ( glCoords[ ( iter_node / numInternalNodesPerEdge ) % numInternalNodesPerEdge ] + 1.0 ) / 2.0;
        real64 gamma = ( glCoords[ ( iter_node / numInternalNodesPerEdge / numInternalNodesPerEdge ) % numInternalNodesPerEdge ] + 1.0 ) / 2.0;

        //GEOS_LOG_RANK_0("TEST TEST: alpha = " << alpha << "; beta = "<< beta << "; gamma= "<< gamma ); 

        localIndex nodeLocalID = localNodeOffset + iter_elem * numInternalNodesPerCell + iter_node;
        elemsToNodesNew[ iter_elem ][ numVerticesPerCell + numEdgesPerCell * numInternalNodesPerEdge + numFacesPerCell * numInternalNodesPerFace + iter_node ] = nodeLocalID;
        nodeLocalToGlobalNew[ nodeLocalID ] = globalNodeOffset + elementLocalToGlobal[ iter_elem ] * numInternalNodesPerCell + iter_node;

        refPosNew[nodeLocalID][0] = refPosSrc[elemVertID[0]][0]*( 1.0-alpha )*( 1.0-beta )*( 1.0-gamma ) + refPosSrc[elemVertID[1]][0]*alpha*( 1.0-beta )*( 1.0-gamma ) +
                                    refPosSrc[elemVertID[2]][0]*( 1.0-alpha )*beta*( 1.0-gamma )         + refPosSrc[elemVertID[3]][0]*alpha*beta*( 1.0-gamma ) +
                                    refPosSrc[elemVertID[4]][0]*( 1.0-alpha )*( 1.0-beta )*gamma         + refPosSrc[elemVertID[5]][0]*alpha*( 1.0-beta )*gamma +
                                    refPosSrc[elemVertID[6]][0]*( 1.0-alpha )*beta*gamma                 + refPosSrc[elemVertID[7]][0]*alpha*beta*gamma;

        refPosNew[nodeLocalID][1] = refPosSrc[elemVertID[0]][1]*( 1.0-alpha )*( 1.0-beta )*( 1.0-gamma ) + refPosSrc[elemVertID[1]][1]*alpha*( 1.0-beta )*( 1.0-gamma )+
                                    refPosSrc[elemVertID[2]][1]*( 1.0-alpha )*beta*( 1.0-gamma )         + refPosSrc[elemVertID[3]][1]*alpha*beta*( 1.0-gamma )+
                                    refPosSrc[elemVertID[4]][1]*( 1.0-alpha )*( 1.0-beta )*gamma         + refPosSrc[elemVertID[5]][1]*alpha*( 1.0-beta )*gamma+
                                    refPosSrc[elemVertID[6]][1]*( 1.0-alpha )*beta*gamma                 + refPosSrc[elemVertID[7]][1]*alpha*beta*gamma;

        refPosNew[nodeLocalID][2] = refPosSrc[elemVertID[0]][2]*( 1.0-alpha )*( 1.0-beta )*( 1.0-gamma ) + refPosSrc[elemVertID[1]][2]*alpha*( 1.0-beta )*( 1.0-gamma )+
                                    refPosSrc[elemVertID[2]][2]*( 1.0-alpha )*beta*( 1.0-gamma )         + refPosSrc[elemVertID[3]][2]*alpha*beta*( 1.0-gamma )+
                                    refPosSrc[elemVertID[4]][2]*( 1.0-alpha )*( 1.0-beta )*gamma         + refPosSrc[elemVertID[5]][2]*alpha*( 1.0-beta )*gamma+
                                    refPosSrc[elemVertID[6]][2]*( 1.0-alpha )*beta*gamma                 + refPosSrc[elemVertID[7]][2]*alpha*beta*gamma;
    
        //GEOS_LOG_RANK_0("TEST TEST: refPosNew["<<nodeLocalID<<"][0]="<<refPosNew[nodeLocalID][0]);
        //GEOS_LOG_RANK_0("TEST TEST: refPosNew["<<nodeLocalID<<"][1]="<<refPosNew[nodeLocalID][1]);
        //GEOS_LOG_RANK_0("TEST TEST: refPosNew["<<nodeLocalID<<"][2]="<<refPosNew[nodeLocalID][2]);

      }

      std::vector<std::tuple<real64, real64, real64, localIndex> > elemRefPos;
      for( localIndex iter_node=0; iter_node<numNodesPerCell; iter_node++ )
      {
        localIndex elemNodeID = elemsToNodesNew[ iter_elem ][iter_node];
        elemRefPos.push_back(std::make_tuple(refPosNew[elemNodeID][2], refPosNew[elemNodeID][1], refPosNew[elemNodeID][0], elemNodeID));
      }

      std::sort(elemRefPos.begin(), elemRefPos.end());
      for( localIndex iter_node=0; iter_node<numNodesPerCell; iter_node++ )
        elemsToNodesNew[ iter_elem ][iter_node] = std::get<3>(elemRefPos[iter_node]);

      /*
      std::cout << std::endl << "ElementID = "<<iter_elem<< std::endl;
      for (int i = 0; i < elemRefPos.size(); i++)
        std::cout << std::get<2>(elemRefPos[i]) << " "
             << std::get<1>(elemRefPos[i]) << " "
             << std::get<0>(elemRefPos[i]) << " "
             << std::get<3>(elemRefPos[i]) << "\n";
      */

    } );

    GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Elements -- 111111 ");
    //GEOS_LOG_RANK_0("TEST TEST: elemsToNodesNew = " << elemsToNodesNew);

    //GEOS_LOG_RANK_0("TEST TEST: refPosSource = " << refPosSource);
    //GEOS_LOG_RANK_0("TEST TEST: refPosNew = " << refPosNew);

    /*
    GEOS_MARK_BEGIN("geos::CellBlockManager::generateHighOrderMaps -- Elements -- 222222");
    //for( localIndex iter_elem = 0; iter_elem < numCellElements; ++iter_elem )
    forAll< RAJA::omp_parallel_for_exec >( numCellElements, 
                                  [ elemsToNodesSource=elemsToNodesSource.toView(), 
                                    elementLocalToGlobal=elementLocalToGlobal.toView(),
                                    elemsToNodesNew=elemsToNodesNew.toView(),
                                    refPosSource=refPosSource.toView(),
                                    refPosNew=refPosNew.toView(),
                                    numVerticesPerCell, numNodesPerCell, numNodesPerEdge, glCoords, order, globalNodeOffset, numInternalNodesPerCell, nullKey,
                                    nodeLocalToGlobalNew=nodeLocalToGlobalNew.toView(), 
                                    &localNodeID, &nodeIDs
                                  ]( localIndex const iter_elem )
    {
      localIndex newCellNodes = 0;
      real64 Xmesh[ numVerticesPerCell ][ 3 ] = { { } };
      real64 X[ 3 ] = { { } };
      localIndex elemMeshVertices[ numVerticesPerCell ] = { };

      for( localIndex iter_vertex = 0; iter_vertex < numVerticesPerCell; iter_vertex++ )
      {
        elemMeshVertices[ iter_vertex ] = elemsToNodesSource[ iter_elem ][ iter_vertex ];
        for( int i =0; i < 3; i++ )
        {
          Xmesh[ iter_vertex ][ i ] = refPosSource[ elemMeshVertices[ iter_vertex ] ][ i ];
        }
      }

      for( int q = 0; q < numNodesPerCell; q++ )
      {
        localIndex nodeID;
        int dof = q;
        int q1 = dof % numNodesPerEdge;
        //GEOS_LOG_RANK_0("TEST TEST: dof = "<<dof<<"; q1 = "<<q1<<"; q="<<q);
        dof /= ( numNodesPerEdge );
        int q2 = dof % ( numNodesPerEdge );
        //GEOS_LOG_RANK_0("TEST TEST: dof = "<<dof<<"; q2 = "<<q2<<"; q="<<q);
        dof /= ( numNodesPerEdge );
        int q3 = dof % ( numNodesPerEdge );
        //GEOS_LOG_RANK_0("TEST TEST: dof = "<<dof<<"; q3 = "<<q3<<"; q="<<q);

        // compute node coords
        real64 alpha = ( glCoords[ q1 ] + 1.0 ) / 2.0;
        real64 beta = ( glCoords[ q2 ] + 1.0 ) / 2.0;
        real64 gamma = ( glCoords[ q3 ] + 1.0 ) / 2.0;
        //GEOS_LOG_RANK_0("TEST TEST: alpha = "<<alpha<<"; beta = "<<beta<<"; gamma="<<gamma);

        // find node ID
        std::array< localIndex, 6 > nodeKey = createNodeKey( elemMeshVertices, q1, q2, q3, order );
        if( nodeKey == nullKey )
        {
          // the node is internal to a cell -- create it
          nodeID = localNodeID;
          nodeLocalToGlobalNew[ nodeID ] = globalNodeOffset + elementLocalToGlobal[ iter_elem ] * numInternalNodesPerCell + newCellNodes;
          localNodeID++;
          newCellNodes++;
        }
        else
        {
          nodeID = nodeIDs.at( nodeKey );
        }

        trilinearInterp( alpha, beta, gamma, Xmesh, X );


        for( int i=0; i<3; i++ )
        {
          refPosNew( nodeID, i ) = X[ i ];
        }

        elemsToNodesNew[ iter_elem ][ q ] = nodeID;
      }
    //}
    } );
    GEOS_MARK_END("geos::CellBlockManager::generateHighOrderMaps -- Elements -- 222222");
    */
  } );
  //GEOS_LOG_RANK_0("TEST TEST: elemsToNodesNew = " << elemsToNodesNew);
  //GEOS_LOG_RANK_0("TEST TEST: nodeLocalToGlobalNew = " << nodeLocalToGlobalNew);
}

}
