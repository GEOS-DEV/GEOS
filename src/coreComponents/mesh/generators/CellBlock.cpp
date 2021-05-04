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

#include "CellBlock.hpp"

namespace geosx
{
using namespace dataRepository;

CellBlock::CellBlock( string const & name, Group * const parent ):
  CellBlockABC( name, parent )
{}

void CellBlock::setElementType( string const & elementType )
{
  m_elementTypeString = elementType;

  if( m_elementTypeString =="C3D8" )
  {
    // Hexahedron
    m_numNodesPerElement = 8;
    m_numEdgesPerElement = 12;
    m_numFacesPerElement = 6;
  }
  else if( m_elementTypeString == "C3D4" )
  {
    // Tetrahedron
    m_numNodesPerElement = 4;
    m_numEdgesPerElement = 6;
    m_numFacesPerElement = 4;
  }
  else if( m_elementTypeString =="C3D6" )
  {
    // Triangular prism
    m_numNodesPerElement = 6;
    m_numEdgesPerElement = 9;
    m_numFacesPerElement = 5;
  }
  else if( m_elementTypeString == "C3D5" )
  {
    // Pyramid
    m_numNodesPerElement = 5;
    m_numEdgesPerElement = 8;
    m_numFacesPerElement = 5;
  }
  else
  {
    GEOSX_ERROR( "Error.  Don't know what kind of element this is." );
  }

  // If `setElementType` is called after the resize, the first dimension would be removed.
  // We do not want that so we try to keep it.
  m_elementsToNodes.resize( this->numCells(), m_numNodesPerElement );
  m_elementsToEdges.resize( this->numCells(), m_numEdgesPerElement );
  m_elementsToFaces.resize( this->numCells(), m_numFacesPerElement );
}

void CellBlock::resize( dataRepository::indexType const numCells )
{
  Group::resize( numCells );

  // Those members are not registered as wrappers because I do not want them
  // to be exposed though the `Group` public interface.
  m_localToGlobalMap.resize( numCells );
  m_elementsToNodes.resize( numCells );
  m_elementsToEdges.resize( numCells );
  m_elementsToFaces.resize( numCells );
}

void CellBlock::getFaceNodes( localIndex iElement,
                              localIndex iFace,
                              array1d< localIndex > & nodesInFaces ) const
{
  geosx::getFaceNodes( m_elementTypeString, iElement, iFace, m_elementsToNodes, nodesInFaces );
}

}
