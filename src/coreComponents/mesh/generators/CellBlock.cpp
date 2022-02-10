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

#include "CellBlock.hpp"

#include "mesh/generators/CellBlockUtilities.hpp"

namespace geosx
{
using namespace dataRepository;

CellBlock::CellBlock( string const & name, Group * const parent ):
  CellBlockABC( name, parent )
{}

void CellBlock::setElementType( ElementType elementType )
{
  m_elementType = elementType;

  switch( m_elementType )
  {
    case ElementType::Hexahedron:
    {
      m_numNodesPerElement = 8;
      m_numEdgesPerElement = 12;
      m_numFacesPerElement = 6;
      break;
    }
    case ElementType::Tetrahedron:
    {
      m_numNodesPerElement = 4;
      m_numEdgesPerElement = 6;
      m_numFacesPerElement = 4;
      break;
    }
    case ElementType::Prism:
    {
      m_numNodesPerElement = 6;
      m_numEdgesPerElement = 9;
      m_numFacesPerElement = 5;
      break;
    }
    case ElementType::Pyramid:
    {
      m_numNodesPerElement = 5;
      m_numEdgesPerElement = 8;
      m_numFacesPerElement = 5;
      break;
    }
    default:
    {
      GEOSX_ERROR( "Invalid element type " << m_elementType << " for CellBlock " << getName() );
    }
  }

  // If `setElementType` is called after the resize, the first dimension would be removed.
  // We do not want that so we try to keep it.
  m_elementsToNodes.resize( this->numElements(), m_numNodesPerElement );
  m_elementsToEdges.resize( this->numElements(), m_numEdgesPerElement );
  m_elementsToFaces.resize( this->numElements(), m_numFacesPerElement );
}

void CellBlock::resize( dataRepository::indexType const numElements )
{
  Group::resize( numElements );

  // Those members are not registered as wrappers because I do not want them
  // to be exposed though the `Group` public interface.
  m_localToGlobalMap.resize( numElements );
  m_elementsToNodes.resize( numElements );
  m_elementsToEdges.resize( numElements );
  m_elementsToFaces.resize( numElements );
}

void CellBlock::getFaceNodes( localIndex iElement,
                              localIndex iFace,
                              localIndex* nodesInFaces ) const
{
  // TODO Change name -- this is another function
  geosx::getFaceNodes( m_elementType, iElement, iFace, m_elementsToNodes, nodesInFaces );
}

}
