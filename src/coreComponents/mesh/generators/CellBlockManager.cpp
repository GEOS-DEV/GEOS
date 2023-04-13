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
#include "mesh/utilities/MeshMapUtilities.hpp"

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
  forAll< parallelHostPolicy >( numNodes, [uniqueFaceOffsets,
                                           lowestNodeToFaces,
                                           faceToNodes,
                                           faceToCells,
                                           faceToBlocks,
                                           &cellBlocks,
                                           &faceBuilder]( localIndex const nodeIndex )
  {
    localIndex nodesInFace[ CellBlockManager::maxNodesPerFace() ];
    localIndex curFaceID = uniqueFaceOffsets[nodeIndex];
    arraySlice1d< NodesAndElementOfFace const > const faces = lowestNodeToFaces[ nodeIndex ];
    forEqualRanges( faces.begin(), faces.end(), [&]( auto first, auto last )
    {
      NodesAndElementOfFace const & f0 = *first;
      CellBlock const & cb = cellBlocks.getGroup< CellBlock >( f0.blockIndex );
      localIndex const numNodesInFace = cb.getFaceNodes( f0.cellIndex, f0.faceNumber, nodesInFace );
      GEOS_ASSERT_EQ( numNodesInFace, f0.numNodes );

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
  GEOSX_MARK_FUNCTION;

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
  GEOSX_MARK_FUNCTION;

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

Group & CellBlockManager::getFaceBlocks()
{
  return this->getGroup( viewKeyStruct::faceBlocks() );
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

}
