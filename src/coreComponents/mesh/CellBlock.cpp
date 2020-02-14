/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CellBlock.cpp
 *
 */

#include "CellBlock.hpp"
#include "MeshLevel.hpp"

#include "NodeManager.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

CellBlock::CellBlock( string const & name, Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_toNodesRelation(),
  m_toFacesRelation()
{
  registerWrapper(viewKeyStruct::nodeListString, &m_toNodesRelation, 0 );
  registerWrapper(viewKeyStruct::faceListString, &m_toFacesRelation, 0 );
  registerWrapper(viewKeyStruct::numNodesPerElementString, &m_numNodesPerElement, 0 );
  registerWrapper(viewKeyStruct::numEdgesPerElementString, &m_numEdgesPerElement, 0 );
  registerWrapper(viewKeyStruct::numFacesPerElementString, &m_numFacesPerElement, 0 );
  registerWrapper(viewKeyStruct::elementCenterString, &m_elementCenter, 0 );
  registerWrapper(viewKeyStruct::elementVolumeString, &m_elementVolume, 0 );
}

CellBlock::~CellBlock()
{}

localIndex CellBlock::GetNumFaceNodes( localIndex const GEOSX_UNUSED_PARAM( elementIndex ),
                                       localIndex const localFaceIndex ) const
{
  if (!m_elementTypeString.compare(0, 4, "C3D8")) return 4;
  if (!m_elementTypeString.compare(0, 4, "C3D6"))
  {
    if (localFaceIndex == 0) return 4;
    if (localFaceIndex == 1) return 4;
    if (localFaceIndex == 2) return 3;
    if (localFaceIndex == 3) return 3;
    if (localFaceIndex == 4) return 4;
  }
  if (!m_elementTypeString.compare(0, 4, "C3D4")) return 3;
  if (!m_elementTypeString.compare(0, 4, "C3D5"))
  {
    if (localFaceIndex == 0) return 4;
    if (localFaceIndex == 1) return 3;
    if (localFaceIndex == 2) return 3;
    if (localFaceIndex == 3) return 3;
    if (localFaceIndex == 4) return 3;
  }

  GEOSX_ERROR("Error. Don't know what kind of element this is and cannot build faces.");
  return -1;
}

localIndex CellBlock::GetFaceNodes( localIndex const elementIndex,
                                    localIndex const localFaceIndex,
                                    localIndex * const nodeIndicies) const
{
  if (!m_elementTypeString.compare(0, 4, "C3D8"))
  {
    if (localFaceIndex == 0)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
    }
    else if (localFaceIndex == 1)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][1];
    }
    else if (localFaceIndex == 2)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][4];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][6];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][2];
    }
    else if (localFaceIndex == 3)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][7];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][5];
    }
    else if (localFaceIndex == 4)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][6];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][7];
    }
    else if (localFaceIndex == 5)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][4];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][7];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][6];
    }
    return 4;
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D6"))
  {
    if (localFaceIndex == 0)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
      return 4;
    }
    else if (localFaceIndex == 1)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][1];
      return 4;
    }
    else if (localFaceIndex == 2)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if (localFaceIndex == 3)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      return 3;
    }
    else if (localFaceIndex == 4)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
      return 4;
    }
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D4"))
  {
    if (localFaceIndex == 0)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][1];
    }
    else if (localFaceIndex == 1)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
    }
    else if (localFaceIndex == 2)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
    }
    else if (localFaceIndex == 3)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
    }

    return 3;
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D5"))
  {
    if (localFaceIndex == 0)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][3];
      return 4;
    }
    else if (localFaceIndex == 1)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if (localFaceIndex == 2)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if (localFaceIndex == 3)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if (localFaceIndex == 4)
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
  }

  GEOSX_ERROR("Error. Don't know what kind of element this is and cannot build faces.");
  return -1;
}

void CellBlock::GetFaceNodes( localIndex const elementIndex,
                              localIndex const localFaceIndex,
                              localIndex_array & nodeIndicies) const
{
  nodeIndicies.resize( GetNumFaceNodes( elementIndex, localFaceIndex ) );
  localIndex const numNodes = GetFaceNodes( elementIndex, localFaceIndex, nodeIndicies.data() );
  GEOSX_DEBUG_VAR( numNodes);
  GEOSX_ASSERT_EQ( numNodes, nodeIndicies.size() );
}

void CellBlock::SetElementType( string const & elementType)
{
  m_elementTypeString = elementType;

  if (!m_elementTypeString.compare(0, 4, "C3D8"))
  {
    m_toNodesRelation.resize(0,8);
    m_toFacesRelation.resize(0,6);
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D4"))
  {
    m_toNodesRelation.resize(0,4);
    m_toFacesRelation.resize(0,4);
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D6"))
  {
    m_toNodesRelation.resize(0,8);
    m_toFacesRelation.resize(0,5);
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D5"))
  {
    m_toNodesRelation.resize(0,5);
    m_toFacesRelation.resize(0,5);
  }
  else
  {
    GEOSX_ERROR("Error.  Don't know what kind of element this is.");
  }

  if (!m_elementTypeString.compare(0, 4, "C3D8"))
  {
    this->numNodesPerElement() = 8;
    this->numFacesPerElement() = 6;
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D4"))
  {
    this->numNodesPerElement() = 4;
    this->numFacesPerElement() = 4;
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D6"))
  {
    this->numNodesPerElement() = 8;
    this->numFacesPerElement() = 5;
  }
  else if (!m_elementTypeString.compare(0, 4, "C3D5"))
  {
    this->numNodesPerElement() = 5;
    this->numFacesPerElement() = 5;
  }

}

void CellBlock::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
  this->m_toFacesRelation.SetRelatedObject( mesh->getFaceManager() );
}

void CellBlock::CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                     FaceManager const & GEOSX_UNUSED_PARAM( facemanager ) )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  forall_in_range<serialPolicy>( 0, this->size(), GEOSX_LAMBDA ( localIndex const k )
  {
    CalculateCellVolumesKernel( k, X );
  });
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlock, std::string const &, Group * const )

}
