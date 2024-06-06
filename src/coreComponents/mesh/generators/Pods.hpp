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

#ifndef GEOS_PODS_HPP
#define GEOS_PODS_HPP

#include "Indices.hpp"

#include "include/NodeMgr.hpp"
#include "include/EdgeMgr.hpp"
#include "include/FaceMgr.hpp"
#include "include/CellMgr.hpp"
#include "include/CellBlk.hpp"

namespace geos
{

class NodeMgrImpl : public generators::NodeMgr
{
public:
  NodeMgrImpl( localIndex numNodes,
               array2d< real64, nodes::REFERENCE_POSITION_PERM > && positions,
               array1d< integer > && ghostRank,
               ArrayOfArrays< localIndex > const & n2e,
               ArrayOfArrays< localIndex > const & n2f,
               ArrayOfArrays< localIndex > const & n2c,
               array1d< globalIndex > const & l2g );

  [[nodiscard]] localIndex numNodes() const override;

  [[nodiscard]] array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getNodeToEdges() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getNodeToFaces() const override;

  [[nodiscard]] ToCellRelation< ArrayOfArrays< localIndex>> getNodeToElements() const override;

  [[nodiscard]] array1d< globalIndex > getLocalToGlobal() const override;

  [[nodiscard]] std::map< string, SortedArray< localIndex > > const & getNodeSets() const override;

private:
  localIndex m_numNodes;
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_positions;
  array1d< integer > m_ghostRank;
  ArrayOfArrays< localIndex > m_n2e;
  ArrayOfArrays< localIndex > m_n2f;
  ArrayOfArrays< localIndex > m_n2c;
  array1d< globalIndex > m_l2g;
  std::map< string, SortedArray< localIndex > > m_todo;
};

class EdgeMgrImpl : public generators::EdgeMgr
{
public:
  EdgeMgrImpl( std::size_t numEdges,
               array1d< integer > && ghostRank,
               array2d< localIndex > && e2n,
               ArrayOfArrays< localIndex > && e2f,
               unordered_map< globalIndex, localIndex > && g2l,
               array1d< globalIndex > && l2g );

  [[nodiscard]] localIndex numEdges() const override;

  [[nodiscard]] array2d< localIndex > getEdgeToNodes() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getEdgeToFaces() const override;

  [[nodiscard]] array1d< integer > getGhostRank() const override;

  [[nodiscard]] array1d< globalIndex > getLocalToGlobal() const override;

private:
  localIndex m_numEdges;
  array1d< integer > m_ghostRank;
  array2d< localIndex > m_e2n;
  ArrayOfArrays< localIndex > m_e2f;
  unordered_map< globalIndex, localIndex > m_g2l;
  array1d< globalIndex > m_l2g;
};

class FaceMgrImpl : public generators::FaceMgr
{
public:
  FaceMgrImpl( std::size_t numFaces,
               array1d< integer > && ghostRank,
               ArrayOfArrays< localIndex > && f2n,
               ArrayOfArrays< localIndex > && f2e,
               array2d< localIndex > && f2c,
               unordered_map< globalIndex, localIndex > && g2l,
               array1d< globalIndex > && l2g );

  [[nodiscard]] localIndex numFaces() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getFaceToNodes() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getFaceToEdges() const override;

  [[nodiscard]] ToCellRelation< array2d< localIndex > > getFaceToElements() const override;

  [[nodiscard]] array1d< integer > getGhostRank() const override;

  [[nodiscard]] array1d< globalIndex > getLocalToGlobal() const override;

private:
  localIndex m_numFaces;
  array1d< integer > m_ghostRank;
  ArrayOfArrays< localIndex > m_f2n;
  ArrayOfArrays< localIndex > m_f2e;
  array2d< localIndex > m_f2c;
  unordered_map< globalIndex, localIndex > m_g2l;
  array1d< globalIndex > m_l2g;
};

class CellBlkImpl : public CellBlk
{
public:
  CellBlkImpl( localIndex numCells,
               array1d< integer > const & ghostRank,
               array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & c2n,
               array2d< localIndex > const & c2e,
               array2d< localIndex > const & c2f,
               array1d< globalIndex > const & l2g );

  [[nodiscard]] ElementType getElementType() const override;

  [[nodiscard]] localIndex numNodesPerElement() const override;

  [[nodiscard]] localIndex numEdgesPerElement() const override;

  [[nodiscard]] localIndex numFacesPerElement() const override;

  [[nodiscard]] localIndex numElements() const override;

  [[nodiscard]] array2d< localIndex, cells::NODE_MAP_PERMUTATION > getElemToNodes() const override;

  [[nodiscard]] array2d< localIndex > getElemToEdges() const override;

  [[nodiscard]] array2d< localIndex > getElemToFaces() const override;

  [[nodiscard]] array1d< globalIndex > localToGlobalMap() const override;

private:
  localIndex m_numCells;
  array1d< integer > m_ghostRank;
  array2d< localIndex, cells::NODE_MAP_PERMUTATION > m_c2n;
  array2d< localIndex > m_c2e;
  array2d< localIndex > m_c2f;
  array1d< globalIndex > m_l2g;
};

} // geos

#endif //GEOS_PODS_HPP
