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

/**
 * @file CellBlock.cpp
 *
 */

#include "CellBlock.hpp"
#include "mesh/MeshLevel.hpp"

#include "mesh/NodeManager.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

CellBlock::CellBlock( string const & name, Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_toFacesRelation(),
  m_externalPropertyNames()
{
  registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation );
  registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation );
}

CellBlock::~CellBlock()
{}

void CellBlock::setElementType( string const & elementType )
{
  m_elementTypeString = elementType;

  if( !m_elementTypeString.compare( 0, 4, "C3D8" ))
  {
    // Hexahedron
    this->setNumNodesPerElement( 8 );
    this->setNumEdgesPerElement( 12 );
    this->setNumFacesPerElement( 6 );
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D4" ))
  {
    // Tetrahedron
    this->setNumNodesPerElement( 4 );
    this->setNumEdgesPerElement( 6 );
    this->setNumFacesPerElement( 4 );
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D6" ))
  {
    // Triangular prism
    this->setNumNodesPerElement( 6 );
    this->setNumEdgesPerElement( 9 );
    this->setNumFacesPerElement( 5 );
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D5" ))
  {
    // Pyramid
    this->setNumNodesPerElement( 5 );
    this->setNumEdgesPerElement( 8 );
    this->setNumFacesPerElement( 5 );
  }
  else
  {
    GEOSX_ERROR( "Error.  Don't know what kind of element this is." );
  }

  m_toNodesRelation.resize( 0, m_numNodesPerElement );
  m_toEdgesRelation.resize( 0, m_numEdgesPerElement );
  m_toFacesRelation.resize( 0, m_numFacesPerElement );

}

void CellBlock::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

void CellBlock::calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                     FaceManager const & GEOSX_UNUSED_PARAM( faceManager ) )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  forAll< serialPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateCellVolumesKernel( k, X );
  } );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlock, string const &, Group * const )

}
