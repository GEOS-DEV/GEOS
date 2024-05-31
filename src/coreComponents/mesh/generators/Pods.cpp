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