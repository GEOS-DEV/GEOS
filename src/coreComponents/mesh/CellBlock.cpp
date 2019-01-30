/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file CellBlock.cpp
 *
 */

#include "CellBlock.hpp"
#include "MeshLevel.hpp"

#include "NodeManager.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
namespace geosx
{
using namespace dataRepository;
//using namespace constitutive;


CellBlock::CellBlock( string const & name, ManagedGroup * const parent ):
  CellBase( name, parent ),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_toFacesRelation()
{
  RegisterViewWrapper(viewKeyStruct::nodeListString, &m_toNodesRelation, 0 );
  RegisterViewWrapper(viewKeyStruct::edgeListString, &m_toEdgesRelation, 0 );
  RegisterViewWrapper(viewKeyStruct::faceListString, &m_toFacesRelation, 0 );
  RegisterViewWrapper(viewKeyStruct::numNodesPerElementString, &m_numNodesPerElement, 0 );
  RegisterViewWrapper(viewKeyStruct::numEdgesPerElementString, &m_numEdgesPerElement, 0 );
  RegisterViewWrapper(viewKeyStruct::numFacesPerElementString, &m_numFacesPerElement, 0 );
  RegisterViewWrapper(viewKeyStruct::elementCenterString, &m_elementCenter, 0 );
  RegisterViewWrapper(viewKeyStruct::elementVolumeString, &m_elementVolume, 0 );
}


CellBlock::~CellBlock()
{}

void CellBlock::GetFaceNodes( const localIndex elementIndex,
                              const localIndex localFaceIndex,
                              localIndex_array& nodeIndicies) const
{
  if (!m_elementType.compare(0, 4, "C3D8"))
  {
    nodeIndicies.resize(4);

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

  }
  else if (!m_elementType.compare(0, 4, "C3D6"))
  {
    if (localFaceIndex == 0)
    {
      nodeIndicies.resize(3);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
    }
    else if (localFaceIndex == 1)
    {
      nodeIndicies.resize(3);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][4];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
    }
    else if (localFaceIndex == 2)
    {
      nodeIndicies.resize(4);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
    }
    else if (localFaceIndex == 3)
    {
      nodeIndicies.resize(4);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][5];
    }
    else if (localFaceIndex == 4)
    {
      nodeIndicies.resize(4);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][5];
    }
  }
  else if (!m_elementType.compare(0, 4, "C3D4"))
  {
    nodeIndicies.resize(3);

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
  }
  else if (!m_elementType.compare(0, 4, "C3D5"))
  {
    if (localFaceIndex == 0)
    {
      nodeIndicies.resize(4);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][3];
    }
    else if (localFaceIndex == 1)
    {
      nodeIndicies.resize(3);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
    }
    else if (localFaceIndex == 2)
    {
      nodeIndicies.resize(3);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
    }
    else if (localFaceIndex == 3)
    {
      nodeIndicies.resize(3);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
    }
    else if (localFaceIndex == 4)
    {
      nodeIndicies.resize(3);
      nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
    }
  }

//  else if ( !m_elementGeometryID.compare(0,4,"CPE2") )
//  {
//    if( localFaceIndex == 0 )
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
//    }
//  }
//
//  else if ( !m_elementGeometryID.compare(0,4,"CPE3") )
//  {
//    if( localFaceIndex == 0 )
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
//    }
//    else if( localFaceIndex == 1 )
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
//    }
//    else if( localFaceIndex == 2 )
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][0];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 4, "CPE4"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
//    }
//    else if (localFaceIndex == 1)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
//    }
//    else if (localFaceIndex == 2)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
//    }
//    else if (localFaceIndex == 3)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][0];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 4, "STRI"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
//    }
//    else if (localFaceIndex == 1)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
//    }
//    else if (localFaceIndex == 2)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][0];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 3, "S4R"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
//      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
//      nodeIndicies[3] = m_toNodesRelation[elementIndex][3];
//    }
//  }
//
//  else if (!m_elementGeometryID.compare(0, 4, "TRSH"))
//  {
//    if (localFaceIndex == 0)
//    {
//      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
//      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
//      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
//    }
//  }
//
  else
  {
    GEOS_ERROR("Error.  Don't know what kind of element this is and cannot build faces.");
  }

}

R1Tensor const & CellBlock::calculateElementCenter( localIndex k,
                                                    const NodeManager& nodeManager,
                                                    const bool ) const
{

  r1_array const & X = nodeManager.referencePosition();
  m_elementCenter[k] = 0;
  for ( localIndex a = 0 ; a < numNodesPerElement() ; ++a)
  {
    const localIndex b = m_toNodesRelation[k][a];
    m_elementCenter[k] += X[b];
  }
  m_elementCenter[k] /= numNodesPerElement();

  return m_elementCenter[k];
}

void CellBlock::SetElementType( string const & elementType)
{
  m_elementType = elementType;

  if (!m_elementType.compare(0, 4, "C3D8"))
  {
    m_toNodesRelation.resize(0,8);
    m_toEdgesRelation.resize(0,12);
    m_toFacesRelation.resize(0,6);
  }
  else if (!m_elementType.compare(0, 4, "C3D4"))
  {
    m_toNodesRelation.resize(0,4);
    m_toEdgesRelation.resize(0,6);
    m_toFacesRelation.resize(0,4);
  }
  else if (!m_elementType.compare(0, 4, "C3D6"))
  {
    m_toNodesRelation.resize(0,6);
    m_toEdgesRelation.resize(0,9);
    m_toFacesRelation.resize(0,5);
  }
  else if (!m_elementType.compare(0, 4, "C3D5"))
  {
    m_toNodesRelation.resize(0,5);
    m_toEdgesRelation.resize(0,8);
    m_toFacesRelation.resize(0,5);
  }
  else
  {
    GEOS_ERROR("Error.  Don't know what kind of element this is.");
  }

  if (!m_elementType.compare(0, 4, "C3D8"))
  {
    this->numNodesPerElement() = 8;
    this->numFacesPerElement() = 6;
  }
  else if (!m_elementType.compare(0, 4, "C3D4"))
  {
    this->numNodesPerElement() = 4;
    this->numFacesPerElement() = 4;
  }
  else if (!m_elementType.compare(0, 4, "C3D6"))
  {
    this->numNodesPerElement() = 6;
    this->numFacesPerElement() = 5;
  }
  else if (!m_elementType.compare(0, 4, "C3D5"))
  {
    this->numNodesPerElement() = 5;
    this->numFacesPerElement() = 5;
  }

}

void CellBlock::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
  this->m_toEdgesRelation.SetRelatedObject( mesh->getEdgeManager() );
  this->m_toFacesRelation.SetRelatedObject( mesh->getFaceManager() );
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellBlock, std::string const &, ManagedGroup * const )

}
