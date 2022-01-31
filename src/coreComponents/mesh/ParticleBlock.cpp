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

/**
 * @file ParticleBlock.cpp
 *
 */

#include "ParticleBlock.hpp"
#include "MeshLevel.hpp"

#include "NodeManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

ParticleBlock::ParticleBlock( string const & name, Group * const parent ):
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

ParticleBlock::~ParticleBlock()
{}

localIndex ParticleBlock::getNumFaceNodes( localIndex const GEOSX_UNUSED_PARAM( elementIndex ),
                                       localIndex const localFaceIndex ) const
{
  switch( m_elementType )
  {
    case ElementType::Hexahedron:
    {
      return 4;
    }
    case ElementType::Prism:
    {
      switch( localFaceIndex )
      {
        case 0:
        case 1:
        case 4: return 4;
        case 2:
        case 3: return 3;
        default:
        {
          GEOSX_ERROR( "Invalid local face index: " << localFaceIndex );
          return -1;
        }
      }
    }
    case ElementType::Tetrahedron:
    {
      return 3;
    }
    case ElementType::Pyramid:
    {
      switch( localFaceIndex )
      {
        case 0: return 4;
        case 1:
        case 2:
        case 3:
        case 4: return 3;
        default:
        {
          GEOSX_ERROR( "Invalid local face index: " << localFaceIndex );
          return -1;
        }
      }
    }
    default:
    {
      GEOSX_ERROR( "Invalid element type: " << m_elementType );
      return -1;
    }
  }
}

localIndex ParticleBlock::getFaceNodes( localIndex const elementIndex,
                                    localIndex const localFaceIndex,
                                    localIndex * const nodeIndicies ) const
{
  switch( m_elementType )
  {
    case ElementType::Hexahedron:
    {
      switch( localFaceIndex )
      {
        case 0:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
          break;
        }
        case 1:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][1];
          break;
        }
        case 2:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][4];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][6];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][2];
          break;
        }
        case 3:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][7];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][5];
          break;
        }
        case 4:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][6];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][7];
          break;
        }
        case 5:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][4];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][5];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][7];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][6];
          break;
        }
        default:
        {
          GEOSX_ERROR( "Invalid local face index: " << localFaceIndex );
          return -1;
        }
      }
      return 4;
    }
    case ElementType::Prism:
    {
      switch( localFaceIndex )
      {
        case 0:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
          return 4;
        }
        case 1:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][1];
          return 4;
        }
        case 2:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
          return 3;
        }
        case 3:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
          return 3;
        }
        case 4:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
          return 4;
        }
        default:
        {
          GEOSX_ERROR( "Invalid local face index: " << localFaceIndex );
          return -1;
        }
      }
    }
    case ElementType::Tetrahedron:
    {
      switch( localFaceIndex )
      {
        case 0:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][1];
          break;
        }
        case 1:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
          break;
        }
        case 2:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
          break;
        }
        case 3:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
          break;
        }
        default:
        {
          GEOSX_ERROR( "Invalid local face index: " << localFaceIndex );
          return -1;
        }
      }
      return 3;
    }
    case ElementType::Pyramid:
    {
      switch( localFaceIndex )
      {
        case 0:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[3] = m_toNodesRelation[elementIndex][3];
          return 4;
        }
        case 1:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
          return 3;
        }
        case 2:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
          return 3;
        }
        case 3:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
          return 3;
        }
        case 4:
        {
          nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
          nodeIndicies[1] = m_toNodesRelation[elementIndex][0];
          nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
          return 3;
        }
        default:
        {
          GEOSX_ERROR( "Invalid local face index: " << localFaceIndex );
          return -1;
        }
      }
    }
    default:
    {
      GEOSX_ERROR( "Invalid element type: " << m_elementType );
      return -1;
    }
  }
}

void ParticleBlock::getFaceNodes( localIndex const elementIndex,
                              localIndex const localFaceIndex,
                              localIndex_array & nodeIndicies ) const
{
  nodeIndicies.resize( getNumFaceNodes( elementIndex, localFaceIndex ) );
  localIndex const numNodes = getFaceNodes( elementIndex, localFaceIndex, nodeIndicies.data() );
  GEOSX_DEBUG_VAR( numNodes );
  GEOSX_ASSERT_EQ( numNodes, nodeIndicies.size() );
}

void ParticleBlock::setElementType( ElementType const elementType )
{
  m_elementType = elementType;

  switch( m_elementType )
  {
    case ElementType::Hexahedron:
    {
      // Hexahedron
      this->setNumNodesPerElement( 8 );
      this->setNumIndependentNodesPerElement( 8 );
      this->setNumEdgesPerElement( 12 );
      this->setNumFacesPerElement( 6 );
      break;
    }
    case ElementType::Tetrahedron:
    {
      // Tetrahedron
      this->setNumNodesPerElement( 4 );
      this->setNumIndependentNodesPerElement( 4 );
      this->setNumEdgesPerElement( 6 );
      this->setNumFacesPerElement( 4 );
      break;
    }
    case ElementType::Prism:
    {
      // Triangular prism
      this->setNumNodesPerElement( 6 );
      this->setNumIndependentNodesPerElement( 6 );
      this->setNumEdgesPerElement( 9 );
      this->setNumFacesPerElement( 5 );
      break;
    }
    case ElementType::Pyramid:
    {
      // Pyramid
      this->setNumNodesPerElement( 5 );
      this->setNumIndependentNodesPerElement( 5 );
      this->setNumEdgesPerElement( 8 );
      this->setNumFacesPerElement( 5 );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Invalid element type: " << m_elementType );
    }
  }

  m_toNodesRelation.resize( 0, m_numNodesPerElement );
  m_toEdgesRelation.resize( 0, m_numEdgesPerElement );
  m_toFacesRelation.resize( 0, m_numFacesPerElement );

}

void ParticleBlock::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

void ParticleBlock::calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                     FaceManager const & GEOSX_UNUSED_PARAM( faceManager ) )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  forAll< serialPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateCellVolumesKernel( k, X );
  } );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleBlock, string const &, Group * const )

}
