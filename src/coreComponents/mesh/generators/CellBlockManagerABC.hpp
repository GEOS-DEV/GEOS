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

#ifndef GEOS_CELLBLOCKMANAGERABC_HPP
#define GEOS_CELLBLOCKMANAGERABC_HPP

#include "CellBlockUtilities.hpp"
#include "dataRepository/Group.hpp"

#include <map>

namespace geos
{

/**
 * @brief Abstract base class for CellBlockManager.
 */
class CellBlockManagerABC : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  CellBlockManagerABC( string const & name, Group * const parent ):
    Group( name, parent )
  {
    // Left blank
  }

  /**
   * @brief Extra space for node to edges mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex edgeMapExtraSpacePerNode()
  { return 8; }

  /**
   * @brief Extra space for node to faces mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex faceMapExtraSpacePerNode()
  { return 8; }

  /**
   * @brief Extra space for node to elements mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex elemMapExtraSpacePerNode()
  { return 8; }

  /**
   * @brief Extra space for extra nodes.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex nodeMapExtraSpacePerFace()
  { return 4; }

  /**
   * @brief Extra space for extra faces.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex edgeMapExtraSpacePerFace()
  { return 4; }

  /**
   * @brief Extra space for extra edges.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex faceMapExtraSpacePerEdge()
  { return 4; }

  /**
   * @brief Returns a group containing the cell blocks as @p CellBlockABC instances.
   * @return Mutable reference to the cell blocks group.
   *
   * @note It should probably be better not to expose a non-const accessor here.
   */
  virtual Group & getCellBlocks() = 0;

  /**
   * @brief Returns a group containing the face blocks as @p FaceBlockABC instances.
   * @return Mutable reference to the face blocks group.
   *
   * @note It should probably be better not to expose a non-const accessor here.
   */
  virtual Group & getFaceBlocks() = 0;

  /**
   * @brief Returns a group containing the cell blocks as CellBlockABC instances
   * @return Const reference to the Group instance.
   */
  virtual const Group & getCellBlocks() const = 0;

  /**
   * @brief Total number of nodes across all the cell blocks.
   * @return The total number of nodes.
   *
   * Nodes shared by multiple cell blocks are counted only once.
   */
  virtual localIndex numNodes() const = 0;

  /**
   * @brief Total number of edges across all the cell blocks.
   * @return The total number of edges.
   */
  virtual localIndex numEdges() const = 0;

  /**
   * @brief Total number of faces across all the cell blocks.
   * @return The total number of faces.
   */
  virtual localIndex numFaces() const = 0;

  /**
   * @brief Returns the node coordinates in a (numNodes, 3) 2d array.
   * @return A const view to the array.
   */
  virtual array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const = 0;

  /**
   * @brief Returns a mutable view of the node coordinates in a (numNodes, 3) 2d array.
   * @return A view to the array.
   */
  virtual arrayView2d< real64, nodes::REFERENCE_POSITION_USD > getNodePositions() = 0;

  /**
   * @brief Returns the node to edges mapping.
   * @return The one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > getNodeToEdges() const = 0;

  /**
   * @brief Returns the face to nodes mappings.
   * @return The one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > getNodeToFaces() const = 0;

  /**
   * @brief Returns the node to elements mapping.
   * @return A one to many relationship.
   */
  virtual ToCellRelation< ArrayOfArrays< localIndex > > getNodeToElements() const = 0;

  /**
   * @brief Returns the edge to nodes mapping.
   * @return A 1 to 2 relationship. The result is meant to have size (numEdges, 2).
   */
  virtual array2d< localIndex > getEdgeToNodes() const = 0;

  /**
   * @brief Returns a mutable view of the edge to nodes mapping.
   * @return A 1 to 2 relationship. The result is meant to have size (numEdges, 2).
   */
  virtual arrayView2d< localIndex > getEdgeToNodes() = 0;

  /**
   * @brief Returns the edge to faces mapping.
   * @return A one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > getEdgeToFaces() const = 0;

  /**
   * @brief Returns the face to nodes mapping.
   * @return The one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > getFaceToNodes() const = 0;

  /**
   * @brief Returns a mutable view of the face to nodes mapping.
   * @return The one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > & getFaceToNodes() = 0;

  /**
   * @brief Returns the face to edges mapping.
   * @return A one to many relationship.
   */
  virtual ArrayOfArrays< localIndex > getFaceToEdges() const = 0;

  /**
   * @brief Returns the face to elements mapping.
   * @return A 1 to 2 relationship. The result is meant to have size (numFaces, 2).
   *
   * In case the face only belongs to one single element, the second value of the table is -1.
   */
  virtual ToCellRelation< array2d< localIndex > > getFaceToElements() const = 0;

  /**
   * @brief Returns a mutable view of the element to nodes mapping, and resizes it.
   * @param name The name of the cell block
   * @param numNodesPartition The number of nodes in the partition, used for resizing
   * @return A one to many relationship.
   */
  virtual arrayView2d< localIndex, cells::NODE_MAP_USD > getElemToNodes( string const & name, localIndex const numNodesPartition ) = 0;

  /**
   * @brief The node to global mapping for nodes.
   * @return The mapping as an array of size numNodes.
   */
  virtual array1d< globalIndex > getNodeLocalToGlobal() const = 0;

  /**
   * @brief Returns a mutable view of the node to global mapping for nodes.
   * @return The mapping as an array of size numNodes.
   */
  virtual arrayView1d< globalIndex > getNodeLocalToGlobal()  = 0;

  /**
   * @brief Sets the number of nodes and the order, and initializes some arrays
   * @param numNodes The number of nodes in the cell block
   * @param order The order of the cell block
   */
  virtual void setNumNodes( localIndex numNodes, localIndex const order ) = 0;

  /**
   * @brief Returns the node sets. Key of the map is the name of the set.
   * @return A reference to constant map.
   */
  virtual std::map< string, SortedArray< localIndex > > const & getNodeSets() const = 0;

  /**
   * @brief Returns a mutable reference to the node sets. Key of the map is the name of the set.
   * @return A mutable reference to the map.
   */
  virtual std::map< string, SortedArray< localIndex > > & getNodeSets() = 0;
};

}
#endif // include guard
