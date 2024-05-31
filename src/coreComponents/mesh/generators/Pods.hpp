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

namespace geos
{

class NodeMgrImpl : public generators::NodeMgr
{
public:
  explicit NodeMgrImpl( NodeLocIdx const & numNodes );

  [[nodiscard]] localIndex numNodes() const override;

  [[nodiscard]] array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getNodeToEdges() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getNodeToFaces() const override;

  [[nodiscard]] ToCellRelation< ArrayOfArrays< localIndex>> getNodeToElements() const override;

  [[nodiscard]] array1d< globalIndex > getLocalToGlobal() const override;

  [[nodiscard]] std::map< string, SortedArray< localIndex > > const & getNodeSets() const override;

private:
  NodeLocIdx m_numNodes;
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
  explicit FaceMgrImpl( FaceLocIdx const & numFaces );

  [[nodiscard]] localIndex numFaces() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getFaceToNodes() const override;

  [[nodiscard]] ArrayOfArrays< localIndex > getFaceToEdges() const override;

  [[nodiscard]] ToCellRelation< array2d< localIndex > > getFaceToElements() const override;

private:
  FaceLocIdx m_numFaces;
};

} // geos

#endif //GEOS_PODS_HPP
