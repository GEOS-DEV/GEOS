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
                          ArrayOfArrays< localIndex > const & n2e,
                          ArrayOfArrays< localIndex > const & n2f,
                          ArrayOfArrays< localIndex > const & n2c,
                          unordered_map< globalIndex, localIndex > && g2l,
                          std::map< integer, array1d< localIndex > > && send,
                          std::map< integer, array1d< localIndex > > && recv )
  : m_ghost{ std::move( g2l ), std::move( send ), std::move( recv ) },
    m_numNodes( numNodes ),
    m_positions( positions ),
    m_n2e( n2e ),
    m_n2f( n2f ),
    m_n2c( n2c )
{ }

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

EdgeMgrImpl::EdgeMgrImpl( std::size_t numEdges,
                          array2d< localIndex > && e2n,
                          ArrayOfArrays< localIndex > && e2f,
                          unordered_map< globalIndex, localIndex > && g2l,
                          std::map< integer, array1d< localIndex > > && send,
                          std::map< integer, array1d< localIndex > > && recv )
  : m_ghost{ std::move( g2l ), std::move( send ), std::move( recv ) },
    m_numEdges( numEdges ),
    m_e2n( e2n ),
    m_e2f( e2f )
{ }

FaceMgrImpl::FaceMgrImpl( std::size_t numFaces,
                          ArrayOfArrays< localIndex > && f2n,
                          ArrayOfArrays< localIndex > && f2e,
                          array2d< localIndex > && f2c,
                          unordered_map< globalIndex, localIndex > && g2l,
                          std::map< integer, array1d< localIndex > > && send,
                          std::map< integer, array1d< localIndex > > && recv )
  : m_ghost{ std::move( g2l ), std::move( send ), std::move( recv ) },
    m_numFaces( numFaces ),
    m_f2n( f2n ),
    m_f2e( f2e ),
    m_f2c( f2c )
{ }

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

CellBlkImpl::CellBlkImpl( localIndex numCells,
                          array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & c2n,
                          array2d< localIndex > const & c2e,
                          array2d< localIndex > const & c2f,
                          unordered_map< globalIndex, localIndex > && g2l,
                          std::map< integer, array1d< localIndex > > && send,
                          std::map< integer, array1d< localIndex > > && recv )
  : m_ghost{ std::move( g2l ), std::move( send ), std::move( recv ) },
    m_numCells( numCells ),
    m_c2n( c2n ),
    m_c2e( c2e ),
    m_c2f( c2f )
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

std::map< string, generators::CellBlk const * > CellMgrImpl::getCellBlks() const
{
  return { { string( "hexahedra" ), &m_cellBlk } };  // TODO hard coded values.
}

} // end of namespace