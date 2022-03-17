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
 * @file HexCellBlockManager.hpp
 */

#ifndef GEOSX_MESH_HEXCELLBLOCKMANAGER_H_
#define GEOSX_MESH_HEXCELLBLOCKMANAGER_H_

#include "mesh/generators/CellBlockManagerBase.hpp"


namespace geosx
{
  
class MeshConnectivityBuilder;

/**
 * @class HexCellBlockManager
 * @brief The HexCellBlockManager specializes CellBlockManagerBase 
 * with a lazy computation strategy and NO overallocation.
 * 
 * Only implemented for hexahedral meshes.
 * The hexahedral mesh may be structured or unstructured.
 * 
 * TODO Implement for other type of cells the MeshConnectivityBuilder
 * It should not be necessary to modify this class.
 * 
 * POTENTIAL ISSUE: Where are the Element indices valid in the maps NodeToElements
 * and FaceToElements? Should it be the index in the CellBlock? 
 */
class HexCellBlockManager : public CellBlockManagerBase
{
public:
  /**
   * @brief Constructor for HexCellBlockManager object.
   * @param name name of this instantiation of CellBlockManagerBase
   * @param parent pointer to the parent Group of this instantiation of CellBlockManagerBase
   */
  HexCellBlockManager( string const & name, Group * const parent );
  HexCellBlockManager( const HexCellBlockManager & ) = delete;
  HexCellBlockManager & operator=( const HexCellBlockManager & ) = delete;
  ~HexCellBlockManager() override;

  localIndex numEdges() const override;
  localIndex numFaces() const override;
  /**
   * @brief Initialize the mapping computations 
   * @details Does not build the maps.
   * Computations are done lazily when calling getters.
   * 
   * @warning MUST be called before any other access to the number of edges / faces
   * or any of the mappings
   */
  void buildMaps() override;

  array2d<geosx::localIndex> getEdgeToNodes() const override;
  ArrayOfSets<geosx::localIndex> getEdgeToFaces() const override;
  ArrayOfArrays<localIndex> getFaceToNodes() const override;
  ArrayOfArrays<geosx::localIndex> getFaceToEdges() const override;
  array2d<localIndex> getFaceToElements() const override;
  ArrayOfSets<localIndex> getNodeToEdges() const override;
  ArrayOfSets<localIndex> getNodeToFaces() const override;
  ArrayOfArrays<localIndex> getNodeToElements() const override;

private:
  /// Instance of the class that build the mappings
  MeshConnectivityBuilder * m_theOneWhoDoesTheJob;
};

}
#endif
