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

#include "Pods.hpp"

namespace geos
{

NodeMgrImpl::NodeMgrImpl( localIndex numNodes,
                          array2d< real64, nodes::REFERENCE_POSITION_PERM > && positions,
                          array1d< integer > && ghostRank,
                          ArrayOfArrays< localIndex > const & n2e,
                          ArrayOfArrays< localIndex > const & n2f,
                          ArrayOfArrays< localIndex > const & n2c,
                          array1d< globalIndex > const & l2g )
  : m_numNodes( numNodes ),
    m_positions( positions ),
    m_ghostRank( ghostRank ),
    m_n2e( n2e ),
    m_n2f( n2f ),
    m_n2c( n2c ),
    m_l2g( l2g )
{ }

localIndex NodeMgrImpl::numNodes() const
{
  return m_numNodes;
}

array2d< real64, nodes::REFERENCE_POSITION_PERM > NodeMgrImpl::getNodePositions() const
{
  return m_positions;
}

ArrayOfArrays< localIndex > NodeMgrImpl::getNodeToEdges() const
{
  return m_n2e;
}

ArrayOfArrays< localIndex > NodeMgrImpl::getNodeToFaces() const
{
  return m_n2f;
}

ToCellRelation< ArrayOfArrays< localIndex > > NodeMgrImpl::getNodeToElements() const
{
  ArrayOfArrays< localIndex > toBlockIndex( m_n2c );  // TODO cheat in case there's one unique cell block!
  for( int i = 0; i < toBlockIndex.size(); ++i )
  {
    for( int j = 0; j < toBlockIndex.sizeOfArray(i); ++j )
    {
      if( m_n2c( i, j ) > -1 )
      {
        toBlockIndex( i, j ) = 0;
      }
    }
  }
  return { std::move( toBlockIndex ), m_n2c };
}

array1d< globalIndex > NodeMgrImpl::getLocalToGlobal() const
{
  return m_l2g;
}

std::map< string, SortedArray< localIndex > > const & NodeMgrImpl::getNodeSets() const
{
  return m_todo;
}

EdgeMgrImpl::EdgeMgrImpl( std::size_t numEdges,
                          array1d< integer > && ghostRank,
                          array2d< localIndex > && e2n,
                          ArrayOfArrays< localIndex > && e2f,
                          unordered_map< globalIndex, localIndex > && g2l,
                          array1d< globalIndex > && l2g )
  : m_numEdges( numEdges ),
    m_ghostRank( ghostRank ),
    m_e2n( e2n ),
    m_e2f( e2f ),
    m_g2l( g2l ),
    m_l2g( l2g )
{ }

localIndex EdgeMgrImpl::numEdges() const
{
  return m_numEdges;
}

array2d< localIndex > EdgeMgrImpl::getEdgeToNodes() const
{
  return m_e2n;
}

ArrayOfArrays< localIndex > EdgeMgrImpl::getEdgeToFaces() const
{
  return m_e2f;
}

array1d< integer > EdgeMgrImpl::getGhostRank() const
{
  return m_ghostRank;
}

array1d< globalIndex > EdgeMgrImpl::getLocalToGlobal() const
{
  return m_l2g;
}

FaceMgrImpl::FaceMgrImpl( std::size_t numFaces,
                          array1d< integer > && ghostRank,
                          ArrayOfArrays< localIndex > && f2n,
                          ArrayOfArrays< localIndex > && f2e,
                          array2d< localIndex > && f2c,
                          unordered_map< globalIndex, localIndex > && g2l,
                          array1d< globalIndex > && l2g )
  : m_numFaces( numFaces ),
    m_ghostRank( ghostRank ),
    m_f2n( f2n ),
    m_f2e( f2e ),
    m_f2c( f2c ),
    m_g2l( g2l ),
    m_l2g( l2g )
{ }

localIndex FaceMgrImpl::numFaces() const
{
  return intConv< localIndex >( m_numFaces );
}

ArrayOfArrays< localIndex > FaceMgrImpl::getFaceToNodes() const
{
  return m_f2n;
}

ArrayOfArrays< localIndex > FaceMgrImpl::getFaceToEdges() const
{
  return m_f2e;
}

ToCellRelation< array2d< localIndex > > FaceMgrImpl::getFaceToElements() const
{
  array2d< localIndex > toBlockIndex( m_f2c );  // TODO cheat in case there's one unique cell block!
  for( int i = 0; i < toBlockIndex.size( 0 ); ++i )
  {
    for( int j = 0; j < toBlockIndex.size( 1 ); ++j )
    {
      if( m_f2c( i, j ) > -1 )
      {
        toBlockIndex( i, j ) = 0;
      }
    }
  }
  return { std::move( toBlockIndex ), m_f2c };
}

array1d< integer > FaceMgrImpl::getGhostRank() const
{
  return m_ghostRank;
}

array1d< globalIndex > FaceMgrImpl::getLocalToGlobal() const
{
  return m_l2g;
}

CellBlkImpl::CellBlkImpl( localIndex numCells,
                          array1d< integer > const & ghostRank,
                          array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & c2n,
                          array2d< localIndex > const & c2e,
                          array2d< localIndex > const & c2f,
                          array1d< globalIndex > const & l2g )
  : m_numCells( numCells ),
    m_ghostRank( ghostRank ),
    m_c2n( c2n ),
    m_c2e( c2e ),
    m_c2f( c2f ),
    m_l2g( l2g )
{ }

ElementType CellBlkImpl::getElementType() const
{
  return ElementType::Hexahedron;
}

localIndex CellBlkImpl::numNodesPerElement() const
{
  return 8;
}

localIndex CellBlkImpl::numEdgesPerElement() const
{
  return 12;
}

localIndex CellBlkImpl::numFacesPerElement() const
{
  return 6;
}

localIndex CellBlkImpl::numElements() const
{
  return m_numCells;
}

array2d< localIndex, cells::NODE_MAP_PERMUTATION > CellBlkImpl::getElemToNodes() const
{
  return m_c2n;
}

array2d< localIndex > CellBlkImpl::getElemToEdges() const
{
  return m_c2e;
}

array2d< localIndex > CellBlkImpl::getElemToFaces() const
{
  return m_c2f;
}

array1d< globalIndex > CellBlkImpl::localToGlobalMap() const
{
  return m_l2g;
}

//std::list< dataRepository::WrapperBase const * > CellBlkImpl::getExternalProperties() const
//{
//  return {};
//}

generators::CellMgr const & MeshMappingImpl::getCellMgr() const
{
  return m_cellMgr;
}

generators::EdgeMgr const & MeshMappingImpl::getEdgeMgr() const
{
  return m_edgeMgr;
}

generators::FaceMgr const & MeshMappingImpl::getFaceMgr() const
{
  return m_faceMgr;
}

generators::NodeMgr const & MeshMappingImpl::getNodeMgr() const
{
  return m_nodeMgr;
}

std::list< CellBlk const * > CellMgrImpl::getCellBlks() const
{
  return { &m_cellBlk };
}

} // end of namespace