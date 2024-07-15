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

#include "include/GhostExchange.hpp"
#include "include/NodeMgr.hpp"
#include "include/EdgeMgr.hpp"
#include "include/FaceMgr.hpp"
#include "include/CellMgr.hpp"
#include "include/CellBlk.hpp"
#include "include/MeshMappings.hpp"

namespace geos
{

struct GhostMapping
{
  unordered_map< globalIndex, localIndex > m_g2l;
  std::map< integer, array1d< localIndex > > m_send;
  std::map< integer, array1d< localIndex > > m_recv;
};

class NodeMgrImpl : public generators::NodeMgr
{
public:
  NodeMgrImpl() = default;

  NodeMgrImpl( localIndex numNodes,
               array2d< real64, nodes::REFERENCE_POSITION_PERM > && positions,
               ArrayOfArrays< localIndex > const & n2e,
               ArrayOfArrays< localIndex > const & n2f,
               ArrayOfArrays< localIndex > const & n2c,
               unordered_map< globalIndex, localIndex > && g2l,
               std::map< integer, array1d< localIndex > > && send,
               std::map< integer, array1d< localIndex > > && recv );

  [[nodiscard]] localIndex numNodes() const override
  {
    return m_numNodes;
  }

  [[nodiscard]] array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const override
  {
    return m_positions;
  }

  [[nodiscard]] ArrayOfArrays< localIndex > getNodeToEdges() const override
  {
    return m_n2e;
  }

  [[nodiscard]] ArrayOfArrays< localIndex > getNodeToFaces() const override
  {
    return m_n2f;
  }

  [[nodiscard]] ToCellRelation< ArrayOfArrays< localIndex > > getNodeToElements() const override;

  [[nodiscard]] std::map< string, SortedArray< localIndex > > const & getNodeSets() const override
  {
    return m_todo;
  }

  // Diamond
  [[nodiscard]] unordered_map< globalIndex, localIndex > getGlobalToLocal() const override
  {
    return m_ghost.m_g2l;
  }

  [[nodiscard]] std::map< integer, array1d< localIndex > > getSend() const override
  {
    return m_ghost.m_send;
  }

  [[nodiscard]] std::map< integer, array1d< localIndex > > getRecv() const override
  {
    return m_ghost.m_recv;
  }
private:
  GhostMapping m_ghost; // Diamond
  localIndex m_numNodes;
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_positions;
  ArrayOfArrays< localIndex > m_n2e;
  ArrayOfArrays< localIndex > m_n2f;
  ArrayOfArrays< localIndex > m_n2c;
  std::map< string, SortedArray< localIndex > > m_todo;
};

class EdgeMgrImpl : public generators::EdgeMgr
{
public:
  EdgeMgrImpl() = default;

  EdgeMgrImpl( std::size_t numEdges,
               array2d< localIndex > && e2n,
               ArrayOfArrays< localIndex > && e2f,
               unordered_map< globalIndex, localIndex > && g2l,
               std::map< integer, array1d< localIndex > > && send,
               std::map< integer, array1d< localIndex > > && recv );

  [[nodiscard]] localIndex numEdges() const override
  {
    return m_numEdges;
  }

  [[nodiscard]] array2d< localIndex > getEdgeToNodes() const override
  {
    return m_e2n;
  }

  [[nodiscard]] ArrayOfArrays< localIndex > getEdgeToFaces() const override
  {
    return m_e2f;
  }

  // Diamond
  [[nodiscard]] unordered_map< globalIndex, localIndex > getGlobalToLocal() const override
  {
    return m_ghost.m_g2l;
  }

  [[nodiscard]] std::map< integer, array1d< localIndex > > getSend() const override
  {
    return m_ghost.m_send;
  }

  [[nodiscard]] std::map< integer, array1d< localIndex > > getRecv() const override
  {
    return m_ghost.m_recv;
  }

private:
  GhostMapping m_ghost; // Diamond
  localIndex m_numEdges;
  array2d< localIndex > m_e2n;
  ArrayOfArrays< localIndex > m_e2f;
};

class FaceMgrImpl : public generators::FaceMgr
{
public:
  FaceMgrImpl() = default;

  FaceMgrImpl( std::size_t numFaces,
               ArrayOfArrays< localIndex > && f2n,
               ArrayOfArrays< localIndex > && f2e,
               array2d< localIndex > && f2c,
               unordered_map< globalIndex, localIndex > && g2l,
               std::map< integer, array1d< localIndex > > && send,
               std::map< integer, array1d< localIndex > > && recv );

  [[nodiscard]] localIndex numFaces() const override
  {
    return intConv< localIndex >( m_numFaces );
  }

  [[nodiscard]] ArrayOfArrays< localIndex > getFaceToNodes() const override
  {
    return m_f2n;
  }

  [[nodiscard]] ArrayOfArrays< localIndex > getFaceToEdges() const override
  {
    return m_f2e;
  }

  [[nodiscard]] ToCellRelation< array2d< localIndex > > getFaceToElements() const override;

  // Diamond
  [[nodiscard]] unordered_map< globalIndex, localIndex > getGlobalToLocal() const override
  {
    return m_ghost.m_g2l;
  }
  [[nodiscard]] std::map< integer, array1d< localIndex > > getSend() const override
  {
    return m_ghost.m_send;
  }

  [[nodiscard]] std::map< integer, array1d< localIndex > > getRecv() const override
  {
    return m_ghost.m_recv;
  }

private:
  GhostMapping m_ghost; // Diamond
  localIndex m_numFaces;
  ArrayOfArrays< localIndex > m_f2n;
  ArrayOfArrays< localIndex > m_f2e;
  array2d< localIndex > m_f2c;
};

class CellBlkImpl : public generators::CellBlk
{
public:
  CellBlkImpl() = default;

  CellBlkImpl( ElementType cellType,
               localIndex numCells,
               array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & c2n,
               array2d< localIndex > const & c2e,
               array2d< localIndex > const & c2f,
               unordered_map< globalIndex, localIndex > && m_g2l,
               std::map< integer, array1d< localIndex > > && send,
               std::map< integer, array1d< localIndex > > && recv );

  [[nodiscard]] ElementType getElementType() const override;

  [[nodiscard]] localIndex numNodesPerElement() const override;

  [[nodiscard]] localIndex numEdgesPerElement() const override;

  [[nodiscard]] localIndex numFacesPerElement() const override;

  [[nodiscard]] localIndex numElements() const override
  {
    return m_numCells;
  }

  [[nodiscard]] array2d< localIndex, cells::NODE_MAP_PERMUTATION > getElemToNodes() const override
  {
    return m_c2n;
  }

  [[nodiscard]] array2d< localIndex > getElemToEdges() const override
  {
    return m_c2e;
  }

  [[nodiscard]] array2d< localIndex > getElemToFaces() const override
  {
    return m_c2f;
  }

  // Diamond
  [[nodiscard]] unordered_map< globalIndex, localIndex > getGlobalToLocal() const override
  {
    return m_ghost.m_g2l;
  }
  [[nodiscard]] std::map< integer, array1d< localIndex > > getSend() const override
  {
    return m_ghost.m_send;
  }

  [[nodiscard]] std::map< integer, array1d< localIndex > > getRecv() const override
  {
    return m_ghost.m_recv;
  }

private:
  GhostMapping m_ghost; // Diamond
  ElementType m_cellType;
  localIndex m_numCells;
  array2d< localIndex, cells::NODE_MAP_PERMUTATION > m_c2n;
  array2d< localIndex > m_c2e;
  array2d< localIndex > m_c2f;
};

class CellMgrImpl : public generators::CellMgr
{
public:
  CellMgrImpl() = default;

  CellMgrImpl( std::vector< CellBlkImpl > && cellBlks )
    :
    m_cellBlk( cellBlks )
  { 
    std::cout << "CellMgr has blocks of:" << std::endl;
    for (size_t i = 0; i < m_cellBlk.size(); ++i)
    {
      std::cout << m_cellBlk[i].getElementType() << std::endl;
    }  
  }

  [[nodiscard]] std::map< string, generators::CellBlk const * > getCellBlks() const override;

private:
  std::vector< CellBlkImpl > m_cellBlk;
};

class MeshMappingImpl : public generators::MeshMappings
{
public:
  MeshMappingImpl( string const & name,
                   Group * const parent )
    : MeshMappings( name, parent )
  { }

  [[nodiscard]] generators::CellMgr const & getCellMgr() const override
  {
    return m_cellMgr;
  }

  [[nodiscard]] generators::EdgeMgr const & getEdgeMgr() const override
  {
    return m_edgeMgr;
  }

  [[nodiscard]] generators::FaceMgr const & getFaceMgr() const override
  {
    return m_faceMgr;
  }

  [[nodiscard]] generators::NodeMgr const & getNodeMgr() const override
  {
    return m_nodeMgr;
  }

  [[nodiscard]] std::set< integer > const & getNeighbors() const override
  {
    return m_neighbors;
  }

  void setNeighbors( std::set< integer > && neighbors )
  {
    m_neighbors = std::move( neighbors );
  }

  void setCellMgr( CellMgrImpl && cellMgr )
  {
    m_cellMgr = cellMgr;
  }

  void setEdgeMgr( EdgeMgrImpl && edgeMgr )
  {
    m_edgeMgr = edgeMgr;
  }

  void setFaceMgr( FaceMgrImpl && faceMgr )
  {
    m_faceMgr = faceMgr;
  }

  void setNodeMgr( NodeMgrImpl && nodeMgr )
  {
    m_nodeMgr = nodeMgr;
  }

private:
  CellMgrImpl m_cellMgr;
  EdgeMgrImpl m_edgeMgr;
  FaceMgrImpl m_faceMgr;
  NodeMgrImpl m_nodeMgr;

  std::set< integer > m_neighbors;
};

} // geos

#endif //GEOS_PODS_HPP
