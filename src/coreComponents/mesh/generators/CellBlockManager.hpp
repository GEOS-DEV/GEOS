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

/**
 * @file CellBlockManager.hpp
 */

#ifndef GEOS_MESH_CELLBLOCKMANAGER_H_
#define GEOS_MESH_CELLBLOCKMANAGER_H_

#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlock.hpp"
#include "mesh/generators/FaceBlock.hpp"

namespace geos
{

/**
 * @class CellBlockManager
 * @brief The CellBlockManager class provides an interface to ObjectManagerBase in order to manage CellBlock data.
 */
class CellBlockManager : public CellBlockManagerABC
{
public:

  /**
   * @brief Constructor for CellBlockManager object.
   * @param name name of this instantiation of CellBlockManager
   * @param parent pointer to the parent Group of this instantiation of CellBlockManager
   */
  CellBlockManager( string const & name, Group * const parent );

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Maximum number of nodes allowed (in memory) per each face.
   * @return The number as an integer.
   */
  static constexpr int maxNodesPerFace()
  { return 64; }

  array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const override;

  /**
   * @brief Returns a view to the vector holding the nodes coordinates
   * @return The reference
   *
   * @note This is meant to be used as a values setter.
   */
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > getNodePositions() override;

  ArrayOfArrays< localIndex > getNodeToEdges() const override;

  ArrayOfArrays< localIndex > getNodeToFaces() const override;

  ToCellRelation< ArrayOfArrays< localIndex > > getNodeToElements() const override;

  array2d< localIndex > getEdgeToNodes() const override;

  /**
   * @brief Returns a view to the vector holding the edges to nodes map
   * @return The reference
   *
   * @note This is meant to be used as a values setter.
   */
  arrayView2d< localIndex > getEdgeToNodes() override;

  ArrayOfArrays< localIndex > getEdgeToFaces() const override;

  ArrayOfArrays< localIndex > getFaceToNodes() const override;

  /**
   * @brief Returns a view to the vector holding the faces to nodes map
   * @return The reference
   *
   * @note This is meant to be used as a values setter.
   */
  ArrayOfArrays< localIndex > & getFaceToNodes() override;

  ArrayOfArrays< localIndex > getFaceToEdges() const override;

  ToCellRelation< array2d< localIndex > > getFaceToElements() const override;

  /**
   * @brief Returns a mutable view to the vector holding the elements to nodes map, after resizing it.
   * @param name The name of the cell block
   * @param numNodesPartition The number of nodes in the partition, used for resizing
   * @return The reference
   *
   * @note This is meant to be used as a values setter.
   */
  arrayView2d< localIndex, cells::NODE_MAP_USD > getElemToNodes( string const & name, localIndex const numNodesPartition ) override;

  array1d< globalIndex > getNodeLocalToGlobal() const override;

  /**
   * @brief Returns a view to the vector holding the node to global mapping.
   * @return The reference
   *
   * @note This is meant to be used as a values setter.
   */
  arrayView1d< globalIndex > getNodeLocalToGlobal() override;

  std::map< string, SortedArray< localIndex > > const & getNodeSets() const override;

  /**
   * @brief Returns a mutable reference to the node sets.
   * @return A reference to the mapping.
   *
   * The key of the map is the name of the set.
   * While the values are sorted arrays which sizes are meant to be managed by the client code.
   * This member function is meant to be used like a setter.
   */
  std::map< string, SortedArray< localIndex > > & getNodeSets() override;

  /**
   * @brief Defines the number of nodes and resizes some underlying arrays appropriately.
   * @param[in] numNodes The number of nodes.
   *
   * The nodes coordinates and nodes local to global mappings get resized to @p numNodes.
   */
  void setNumNodes( localIndex numNodes ); // TODO Improve doc. Is it per domain, are there duplicated nodes because of subregions?

  /**
   * @brief Defines the number of nodes and the order, and resizes some underlying arrays appropriately.
   * @param[in] numNodes The number of nodes.
   * @param[in] order The order of the cell block
   *
   * The nodes coordinates and nodes local to global mappings get resized to @p numNodes.
   */
  void setNumNodes( localIndex numNodes, localIndex const order ) override;

  localIndex numNodes() const override;

  localIndex numEdges() const override;

  localIndex numFaces() const override;

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param[in] numElements list of the new element numbers
   * @param[in] regionNames list of the element region names
   */
  void resize( integer_array const & numElements,
               string_array const & regionNames );

  /**
   * @brief Trigger the computation of all the mappings.
   *
   * Call this member function to compute all the mappings.
   * Computations could be done lazily when calling getters.
   * But this is not yet implemented.
   */
  void buildMaps();

  /**
   * @brief Get cell block by name.
   * @param[in] name Name of the cell block.
   * @return Reference to the cell block instance.
   */
  CellBlock & getCellBlock( string const & name )
  {
    return this->getGroup( viewKeyStruct::cellBlocks() ).getGroup< CellBlock >( name );
  }

  const Group & getCellBlocks() const override;

  Group & getCellBlocks() override;

  Group & getFaceBlocks() override;

  /**
   * @brief Registers and returns a cell block of name @p name.
   * @param name The name of the created cell block.
   * @return A reference to the new cell block. The CellBlockManager owns this new instance.
   */
  CellBlock & registerCellBlock( string const & name );

  /**
   * @brief Registers and returns a face block of name @p name.
   * @param name The name of the created face block.
   * @return A reference to the new face block. The CellBlockManager owns this new instance.
   */
  FaceBlock & registerFaceBlock( string const & name );

  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA lambda )
  {
    this->getGroup( viewKeyStruct::cellBlocks() ).forSubGroups< CellBlock >( lambda );
  }

private:

  struct viewKeyStruct
  {
    /// Cell blocks key
    static constexpr char const * cellBlocks()
    { return "cellBlocks"; }

    /// Face blocks key
    static constexpr char const * faceBlocks()
    { return "faceBlocks"; }
  };

  /**
   * @brief Get cell block at index @p blockIndex.
   * @param[in] blockIndex The cell block index.
   * @return Const reference to the instance.
   *
   * @note Mainly useful for iteration purposes.
   */
  CellBlock const & getCellBlock( localIndex const blockIndex ) const;

  /**
   * @brief Returns the number of cells blocks
   * @return Number of cell blocks
   */
  localIndex numCellBlocks() const;

  /**
   * @brief Trigger the node to edges mapping computation.
   */
  void buildNodeToEdges();

  /**
   * @brief Trigger the face to nodes, edges and elements mappings.
   */
  void buildFaceMaps();

  template< typename BASEMAP, typename FUNC >
  void buildToCellMap( localIndex const cellIndex,
                       ToCellRelation< BASEMAP > & toCells,
                       FUNC cellToObjectGetter,
                       localIndex const overAlloc = 0 ) const;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_nodesPositions;

  ArrayOfArrays< localIndex > m_nodeToEdges;
  ArrayOfArrays< localIndex > m_edgeToFaces;
  array2d< localIndex > m_edgeToNodes;
  ArrayOfArrays< localIndex > m_faceToNodes;
  ArrayOfArrays< localIndex > m_faceToEdges;
  ToCellRelation< array2d< localIndex > > m_faceToCells;

  array1d< globalIndex > m_nodeLocalToGlobal;

  std::map< string, SortedArray< localIndex > > m_nodeSets;

  localIndex m_numNodes;
  localIndex m_numFaces;
  localIndex m_numEdges;
};

}
#endif /* GEOS_MESH_CELLBLOCKMANAGER_H_ */
