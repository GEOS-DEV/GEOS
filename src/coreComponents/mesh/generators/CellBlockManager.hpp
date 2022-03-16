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

#ifndef GEOSX_MESH_CELLBLOCKMANAGER_H_
#define GEOSX_MESH_CELLBLOCKMANAGER_H_

#include "mesh/generators/CellBlockManagerBase.hpp"
#include "mesh/generators/CellBlock.hpp"

namespace geosx
{

/**
 * @class CellBlockManager
 * @brief The CellBlockManager class provides an interface to ObjectManagerBase in order to manage CellBlock data.
 */
class CellBlockManager : public CellBlockManagerBase
{
public:

  /**
   * @brief Constructor for CellBlockManager object.
   * @param name name of this instantiation of CellBlockManager
   * @param parent pointer to the parent Group of this instantiation of CellBlockManager
   */
  CellBlockManager( string const & name, Group * const parent );

  CellBlockManager( const CellBlockManager & ) = delete;

  CellBlockManager & operator=( const CellBlockManager & ) = delete;

  ~CellBlockManager() override = default;

  /**
   * @brief Maximum number of faces allowed (in memory) per each node.
   * @return The number as an integer.
   */
  static constexpr int maxFacesPerNode()
  { return 200; }

  /**
   * @brief Extra space for node to faces mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex getFaceMapOverallocation()
  { return 8; }


  localIndex numEdges() const override;

  localIndex numFaces() const override;

  ArrayOfSets< localIndex > getNodeToEdges() const override;

  ArrayOfSets< localIndex > getNodeToFaces() const override;

  ArrayOfArrays< localIndex > getNodeToElements() const override;

  array2d< geosx::localIndex > getEdgeToNodes() const override;

  ArrayOfSets< geosx::localIndex > getEdgeToFaces() const override;

  ArrayOfArrays< localIndex > getFaceToNodes() const override;

  ArrayOfArrays< geosx::localIndex > getFaceToEdges() const override;

  array2d< localIndex > getFaceToElements() const override;

  
  /**
   * @brief Trigger the computation of all the mappings.
   *
   * Call this member function to compute all the mappings.
   * Computations could be done lazily when calling getters.
   * But this is not yet implemented.
   */
  void buildMaps() override;

private:

  /**
   * @brief Trigger the node to edges mapping computation.
   */
  void buildNodeToEdges();

  /**
   * @brief Trigger the face to nodes, edges and elements mappings.
   */
  void buildFaceMaps();

  localIndex m_numFaces;
  localIndex m_numEdges;

  ArrayOfSets< localIndex > m_nodeToEdges;
  ArrayOfSets< localIndex > m_edgeToFaces;
  array2d< localIndex > m_edgeToNodes;
  ArrayOfArrays< localIndex >  m_faceToNodes;
  ArrayOfArrays< localIndex > m_faceToEdges;
  array2d< localIndex >  m_faceToElements;

};

}
#endif /* GEOSX_MESH_CELLBLOCKMANAGER_H_ */
