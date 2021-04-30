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
 * @file CellBlockManager.hpp
 */

#ifndef GEOSX_MESH_CELLBLOCKMANAGER_H_
#define GEOSX_MESH_CELLBLOCKMANAGER_H_

#include "mesh/ObjectManagerBase.hpp"
#include "mesh/generators/CellBlock.hpp"

#include <map>
#include <vector>

namespace geosx
{

/**
 * @class CellBlockManager
 * @brief The CellBlockManager class provides an interface to ObjectManagerBase in order to manage CellBlock data.
 */
class CellBlockManager : public ObjectManagerBase
{
public:

  /**
   * @brief The function is to return the name of the CellBlockManager in the object catalog
   * @return string that contains the catalog name used to register/lookup this class in the object catalog
   */
  static string catalogName()
  {
    return "CellBlockManager";
  }

  virtual const string getCatalogName() const override final
  { return CellBlockManager::catalogName(); }


  /**
   * @brief Constructor for CellBlockManager object.
   * @param name name of this instantiation of CellBlockManager
   * @param parent pointer to the parent Group of this instantiation of CellBlockManager
   */
  CellBlockManager( string const & name, Group * const parent );

  /**
   * @brief Destructor
   */
  ~CellBlockManager() override = default;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  static constexpr int maxEdgesPerNode()
  { return 200; }

  static constexpr int maxFacesPerNode()
  { return 200; }

  static constexpr localIndex getEdgeMapOverallocation()
  { return 8; }

  static constexpr localIndex getFaceMapOverallocation()
  { return 8; }

  static constexpr localIndex getElemMapOverAllocation()
  { return 8; }

  static constexpr localIndex nodeMapExtraSpacePerFace()
  { return 4; }

  static constexpr localIndex edgeMapExtraSpacePerFace()
  { return 4; }

  static constexpr localIndex faceMapExtraSpacePerEdge()
  { return 4; }

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param numElements list of the new element numbers
   * @param regionNames list of the element region names
   * @param elementTypes list of the element types
   */
  void resize( integer_array const & numElements,
               string_array const & regionNames,
               string_array const & elementTypes );

  /**
   * @brief Get cell block by name.
   * @param name Name of the cell block.
   * @return Reference to the cell block instance.
   */
  CellBlock & getCellBlock( string const & name )
  {
    return this->getGroup( cellBlocksKey ).getGroup< CellBlock >( name );
  }

  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA lambda )
  {
    this->getGroup( cellBlocksKey ).forSubGroups< CellBlock >( lambda );
  }

  /**
   * @brief Returns the node to elements mappings.
   * @return A one to many relationship.
   *
   * @note The mapping is computed on the fly and returned. It is not stored in the instance.
   */
  ArrayOfArrays< localIndex > getNodeToElements() const;

  /**
   * @brief Returns the face to elements mappings.
   * @return A one to many relationship.
   *
   * In case the face only belongs to one single element, the second value of the table is -1.
   */
  array2d< localIndex > getFaceToElements() const;

  /**
   * @brief Returns the face to nodes mappings.
   * @param numNodes This should not be here, needs to be removed TODO
   * @return The one to many relationship.
   */
  ArrayOfArrays< localIndex > getFaceToNodes() const;

  /**
   * @brief Returns the face to nodes mappings.
   * @return The one to many relationship.
   */
  ArrayOfSets< localIndex > getNodeToFaces() const;

  ArrayOfSets< localIndex > getNodeToEdges() const;

  ArrayOfSets< geosx::localIndex > const & getEdgeToFaces() const;

  array2d< geosx::localIndex > const & getEdgeToNodes() const;

  ArrayOfArrays< geosx::localIndex > const & getFaceToEdges() const;

  void buildMaps();

  /**
   * @brief Total number of nodes across all the cell blocks.
   * @return
   */
  localIndex numNodes() const; // TODO Improve doc

  /**
   * @brief Defines the number of nodes.
   * @param numNodes The number of nodes.
   *
   * @note This has to be defined consistent w.r.t. the points positions defined in the node manager.
   * @deprecated To be more consistent, define the points positions instead and fill the node manager with it.
   */
  void setNumNodes( localIndex numNodes ) // TODO Improve doc. Is it per domain, are there duplicated nodes because of subregions?
  { m_numNodes = numNodes; }

  /**
   * @brief Total number of faces across all the cell blocks.
   * @return
   */
  localIndex numFaces() const; // TODO Improve doc

  localIndex numEdges() const; // TODO Improve doc

  /**
   * @brief Returns a group containing the cell blocks as CellBlockABC instances
   * @return
   */
  Group & getCellBlocks();

private:

  /**
   * @brief Cell key
   */
  static string const cellBlocksKey;

  /**
   * @brief Copy constructor.
   */
  CellBlockManager( const CellBlockManager & ) = delete;

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  CellBlockManager & operator=( const CellBlockManager & ) = delete;

  /**
   * @brief Returns a group containing the cell blocks as CellBlockABC instances
   * @return Reference to the Group instance.
   */
  const Group & getCellBlocks() const;

  /**
   * @brief Get cell block at index @p iCellBlock.
   * @param iCellBlock The cell block index.
   * @return Const reference to the instance.
   *
   * @note Mainly useful for iteration purposes.
   */
  const CellBlockABC & getCellBlock( localIndex iCellBlock ) const;

  localIndex numCellBlocks() const;

  void buildNodeToEdges();
  void buildFaceMaps();

  ArrayOfArrays< localIndex >  m_faceToNodes;
  array2d< localIndex >  m_faceToElements;
  ArrayOfSets< localIndex > m_edgeToFaces;
  array2d< localIndex > m_edgeToNodes;
  ArrayOfSets< localIndex > m_nodeToEdges;
  ArrayOfArrays< localIndex > m_faceToEdges;
  localIndex m_numNodes;
  localIndex m_numFaces;
  localIndex m_numEdges;
};

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
                          ArrayOfSets< localIndex > & edgeToFaceMap,
                          array2d< localIndex > & edgeToNodeMap );

}
#endif /* GEOSX_MESH_CELLBLOCKMANAGER_H_ */
