/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_MESH_GENERATORS_CELLBLOCKUTILITIES_HPP_
#define GEOSX_MESH_GENERATORS_CELLBLOCKUTILITIES_HPP_

#include "mesh/ElementType.hpp"
#include "common/DataTypes.hpp"
#include "common/Span.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

/**
 * @brief Free function that generates face to edges, edge to faces and edge to nodes mappings.
 * @param[in] numNodes The number of nodes.
 * @param[in] faceToNodeMap Face to node mappings as an input.
 * @param[out] faceToEdgeMap Face to edges will be resized and filled.
 * @param[out] edgeToFaceMap Ege to faces will be resized and filled.
 * @param[out] edgeToNodeMap Edge to nodes will be resized and filled.
 * @return The number of edges.
 */
localIndex buildEdgeMaps( localIndex numNodes,
                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                          ArrayOfArrays< localIndex > & faceToEdgeMap,
                          ArrayOfArrays< localIndex > & edgeToFaceMap,
                          array2d< localIndex > & edgeToNodeMap );

/**
 * @brief Get the local indices of the nodes in a face of the element.
 * @param[in] elementType Type of the element
 * @param[in] elemIdx the local index of the target element
 * @param[in] faceNumber the local index of the target face in the element  (this will be 0-numFacesInElement)
 * @param[in] elementToNodes Element to nodes mapping.
 * @param[out] faceNodes space to write node indices to, must be of sufficient size
 * @return number of nodes in the face
 */
localIndex getFaceNodes( ElementType const elementType,
                         localIndex const elemIdx,
                         localIndex const faceNumber,
                         arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elementToNodes,
                         Span< localIndex > const faceNodes );

/**
 * @brief Find and count ranges of repeated values in an array of sorted arrays and compute offsets.
 * @tparam POLICY execution policy
 * @tparam T value type of input arrays
 * @param [in] sortedLists the input array of sorted arrays of values
 * @return an array of size @p sortedLists.size() + 1, where element i contains starting index of
 *         unique values from sub-array i in a global list of unique values, while the last value
 *         contains the total number of unique edges (i.e. an exclusive scan + total reduction).
 */
template< typename POLICY, typename T >
array1d< localIndex >
computeUniqueValueOffsets( ArrayOfArraysView< T const > const & sortedLists )
{
  localIndex const numNodes = sortedLists.size();
  array1d< localIndex > uniqueValueOffsets( numNodes + 1 );
  uniqueValueOffsets[0] = 0;

  // For each node, count number of unique edges that have the node as its lowest
  arrayView1d< localIndex > const numUniqueValuesView = uniqueValueOffsets.toView();
  forAll< POLICY >( numNodes, [sortedLists, numUniqueValuesView]( localIndex const i )
  {
    numUniqueValuesView[ i + 1 ] = 0;
    arraySlice1d< T const > const list = sortedLists[ i ];
    forEqualRanges( list.begin(), list.end(), [&]( auto, auto )
    {
      ++numUniqueValuesView[ i + 1 ];
    } );
  } );

  // Perform an inplace prefix-sum to get the unique edge offset.
  RAJA::inclusive_scan_inplace< POLICY >( uniqueValueOffsets );
  return uniqueValueOffsets;
}

}

#endif // GEOSX_MESH_GENERATORS_CELLBLOCKUTILITIES_HPP_
