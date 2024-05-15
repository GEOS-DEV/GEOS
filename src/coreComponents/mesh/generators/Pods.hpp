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

  localIndex numNodes() const override;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const override;

  ArrayOfArrays< localIndex > getNodeToEdges() const override;

  ArrayOfArrays< localIndex > getNodeToFaces() const override;

  ToCellRelation< ArrayOfArrays< localIndex>> getNodeToElements() const override;

  array1d< globalIndex > getLocalToGlobal() const override;

  std::map< string, SortedArray< localIndex > > const & getNodeSets() const override;

private:
  NodeLocIdx m_numNodes;
  std::map< string, SortedArray< localIndex > > m_todo;
};

class EdgeMgrImpl : public generators::EdgeMgr
{
public:
  explicit EdgeMgrImpl( EdgeLocIdx const & numEdges );

  localIndex numEdges() const override;

  array2d< localIndex > getEdgeToNodes() const override;

  ArrayOfArrays< localIndex > getEdgeToFaces() const override;

  array1d< integer > getGhostRank() const override;

private:
  EdgeLocIdx m_numEdges;
};

class FaceMgrImpl : public generators::FaceMgr
{
public:
  explicit FaceMgrImpl( FaceLocIdx const & numFaces );

  localIndex numFaces() const override;

  ArrayOfArrays< localIndex > getFaceToNodes() const override;

  ArrayOfArrays< localIndex > getFaceToEdges() const override;

  ToCellRelation< array2d< localIndex > > getFaceToElements() const override;

private:
  FaceLocIdx m_numFaces;
};

} // geos

#endif //GEOS_PODS_HPP
