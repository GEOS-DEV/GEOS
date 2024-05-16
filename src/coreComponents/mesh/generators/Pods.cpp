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

NodeMgrImpl::NodeMgrImpl( NodeLocIdx const & numNodes )
  : m_numNodes( numNodes )
{ }

localIndex NodeMgrImpl::numNodes() const
{
  return intConv< localIndex >( m_numNodes.get() );
}

array2d< real64, nodes::REFERENCE_POSITION_PERM > NodeMgrImpl::getNodePositions() const
{
  return  {};
}

ArrayOfArrays< localIndex > NodeMgrImpl::getNodeToEdges() const
{
  return  {};
}

ArrayOfArrays< localIndex > NodeMgrImpl::getNodeToFaces() const
{
  return  {};
}

ToCellRelation< ArrayOfArrays< localIndex > > NodeMgrImpl::getNodeToElements() const
{
  return  {};
}

array1d< globalIndex > NodeMgrImpl::getLocalToGlobal() const
{
  return  {};
}

std::map< string, SortedArray< localIndex > > const & NodeMgrImpl::getNodeSets() const
{
  return m_todo;
}

EdgeMgrImpl::EdgeMgrImpl( EdgeLocIdx const & numEdges,
                          std::vector< MpiRank > && ghostRank,
                          std::vector< EdgeGlbIdx > && l2g )
  : m_numEdges( numEdges ),
    m_ghostRank( ghostRank ),
    m_l2g( l2g )
{ }

localIndex EdgeMgrImpl::numEdges() const
{
  return intConv< localIndex >( m_numEdges.get() );
}

array2d< localIndex > EdgeMgrImpl::getEdgeToNodes() const
{
  return {};
}

ArrayOfArrays< localIndex > EdgeMgrImpl::getEdgeToFaces() const
{
  return {};
}

array1d< integer > EdgeMgrImpl::getGhostRank() const
{
  array1d< integer > result;
  result.reserve( std::size( m_ghostRank ) );
  for( MpiRank const & rank: m_ghostRank )
  {
    result.emplace_back( rank.get() );
  }
  return result;
}

array1d< globalIndex > EdgeMgrImpl::getLocalToGlobal() const
{
  array1d< globalIndex > result;
  result.reserve( std::size( m_l2g ) );
  for( EdgeGlbIdx const & edge: m_l2g )
  {
    result.emplace_back( edge.get() );
  }
  return result;
}

FaceMgrImpl::FaceMgrImpl( FaceLocIdx const & numFaces )
  : m_numFaces( numFaces )
{ }

localIndex FaceMgrImpl::numFaces() const
{
  return intConv< localIndex >( m_numFaces.get() );
}

ArrayOfArrays< localIndex > FaceMgrImpl::getFaceToNodes() const
{
  return {};
}

ArrayOfArrays< localIndex > FaceMgrImpl::getFaceToEdges() const
{
  return {};
}

ToCellRelation< array2d< localIndex > > FaceMgrImpl::getFaceToElements() const
{
  return {};
}

} // geos