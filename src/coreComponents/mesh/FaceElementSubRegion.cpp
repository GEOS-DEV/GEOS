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
 * @file FaceElementSubRegion.cpp
 */

#include "FaceElementSubRegion.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"

namespace geosx
{

FaceElementSubRegion::FaceElementSubRegion( string const & name,
                                      dataRepository::ManagedGroup * const parent ):
  ElementSubRegionBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::nodeListString, &m_toNodesRelation, false )->
    setDescription("Map to the nodes attached to each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation, false )->
    setDescription("Map to the edges attached to each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::faceListString, &m_toFacesRelation, false )->
    setDescription("Map to the faces attached to each FaceElement.")->
    reference().resize(0,2);

  RegisterViewWrapper( viewKeyStruct::elementApertureString, &m_elementAperture, false )->
    setApplyDefaultValue(1.0e-5)->
    setDescription("The aperture of each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::elementAreaString, &m_elementArea, false )->
    setApplyDefaultValue(-1.0)->
    setDescription("The area of each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::elementCenterString, &m_elementCenter, false )->
    setDescription("The center of each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::elementVolumeString, &m_elementVolume, false )->
    setApplyDefaultValue(-1.0)->
    setDescription("The volume of each FaceElement.");

  m_numNodesPerElement = 8;
}

FaceElementSubRegion::~FaceElementSubRegion()
{}


R1Tensor const & FaceElementSubRegion::calculateElementCenter( localIndex k,
                                                               const NodeManager& nodeManager,
                                                               const bool useReferencePos ) const
{
  r1_array const & X = nodeManager.referencePosition();
  m_elementCenter[k] = 0;
  localIndex const numNodes = numNodesPerElement( k );
  for ( localIndex a = 0 ; a < numNodes ; ++a)
  {
    const localIndex b = m_toNodesRelation[k][a];
    m_elementCenter[k] += X[b];
    if(!useReferencePos)
      m_elementCenter[k] += X[b];
  }
  m_elementCenter[k] /= numNodes;

  return m_elementCenter[k];

}

void FaceElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
  this->m_toEdgesRelation.SetRelatedObject( mesh->getEdgeManager() );
  this->m_toFacesRelation.SetRelatedObject( mesh->getFaceManager() );
}

void FaceElementSubRegion::CalculateElementGeometricQuantities( localIndex const k,
                                                                arrayView1d<real64 const> const & faceArea )
{
  m_elementArea[k] = faceArea[ m_toFacesRelation[k][0] ];
  m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
}

void FaceElementSubRegion::CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                                FaceManager const & faceManager )
{
  arrayView1d<real64 const> const & faceArea = faceManager.faceArea();

  forall_in_range<elemPolicy>( 0, this->size(), GEOSX_LAMBDA ( localIndex const k )
  {
    m_elementArea[k] = faceArea[ m_toFacesRelation[k][0] ];
    m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
  });
}



} /* namespace geosx */
