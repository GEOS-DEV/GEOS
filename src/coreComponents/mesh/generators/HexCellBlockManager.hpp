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

#include "common/DataTypes.hpp"

#include <cassert>

namespace geosx
{
  
class HexMeshConnectivityBuilder;

/**
 * @class HexCellBlockManager
 * @brief The HexCellBlockManager specializes CellBlockManagerBase for full hexahedral meshes  
 * The hexahedral mesh may be structured or unstructured.
 */
class HexCellBlockManager : public CellBlockManagerBase
{
public:

  /**
   * @brief Constructor for HexCellBlockManager object.
   * @param name name of this instantiation of HexCellBlockManager
   * @param parent pointer to the parent Group of this instantiation of HexCellBlockManager
   */
  HexCellBlockManager( string const & name, Group * const parent );
  HexCellBlockManager( const HexCellBlockManager & ) = delete;
  HexCellBlockManager & operator=( const HexCellBlockManager & ) = delete;
  ~HexCellBlockManager() override;

  localIndex numEdges() const override
  { return m_numEdges; }

  localIndex numFaces() const override
  { return m_numFaces; }

  localIndex numElements() const
  { return m_numElements; }

  array2d<geosx::localIndex> getEdgeToNodes() override;
  ArrayOfSets<geosx::localIndex> getEdgeToFaces() override;

  ArrayOfArrays<localIndex> getFaceToNodes() override;
  ArrayOfArrays<geosx::localIndex> getFaceToEdges() override;

  // TODO We have a problem - where are the Element index valid ?
  array2d<localIndex> getFaceToElements() override;

  ArrayOfSets<localIndex> getNodeToEdges() override;
  ArrayOfSets<localIndex> getNodeToFaces() override;

  // TODO We have a problem - where are the Element index valid ?
  ArrayOfArrays<localIndex> getNodeToElements() override;

  /**
   * @brief Compute all possible maps and more
   */
  void buildMaps();

private:
  HexMeshConnectivityBuilder * m_theOneWhoDoesTheJob;

  // The numbers of things we are dealing with 
  localIndex m_numEdges = 0;
  localIndex m_numFaces = 0;
  localIndex m_numElements = 0;
};

}
#endif
