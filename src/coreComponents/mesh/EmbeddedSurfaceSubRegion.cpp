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
 * @file EmbeddedSurfaceSubRegion.cpp
 */

#include "EmbeddedSurfaceSubRegion.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"

namespace geosx
{
using namespace dataRepository;


EmbeddedSurfaceSubRegion::EmbeddedSurfaceSubRegion( string const & name,
                                      dataRepository::Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_normalVector(),
  m_embeddedSurfaceToCell(),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_toFacesRelation(),
  m_elementAperture(),
  m_elementArea()
{
  registerWrapper( viewKeyStruct::nodeListString, &m_toNodesRelation, false )->
    setDescription("Map to the nodes attached to each EmbeddedSurface.");

  registerWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation, false )->
      setDescription("Map to the edges attached to each FaceElement.");

  registerWrapper( viewKeyStruct::faceListString, &m_toFacesRelation, false )->
      setDescription("Map to the faces attached to each FaceElement.")->
      reference().resize(0,2);

  registerWrapper( viewKeyStruct::cellListString, &m_embeddedSurfaceToCell, false )->
        setDescription(".");

  registerWrapper( viewKeyStruct::normalVectorString, &m_normalVector, false )->
          setDescription(".");

  registerWrapper( viewKeyStruct::elementApertureString, &m_elementAperture, false )->
    setApplyDefaultValue(1.0e-5)->
    setPlotLevel(dataRepository::PlotLevel::LEVEL_0)->
    setDescription("The aperture of each EmbeddedSurface.");

  registerWrapper( viewKeyStruct::elementAreaString, &m_elementArea, false )->
    setApplyDefaultValue(-1.0)->
    setDescription("The area of each EmbeddedSurface.");

  registerWrapper( viewKeyStruct::elementCenterString, &m_elementCenter, false )->
    setDescription("The center of each EmbeddedSurface.");

  registerWrapper( viewKeyStruct::elementVolumeString, &m_elementVolume, false )->
    setApplyDefaultValue(-1.0)->
    setPlotLevel(dataRepository::PlotLevel::LEVEL_0)->
    setDescription("The volume of each EmbeddedSurface.");

  m_numNodesPerElement = 4; // Let s assume it's a plane for now
}


EmbeddedSurfaceSubRegion::~EmbeddedSurfaceSubRegion()
{}


void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( localIndex const k)
{
  m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
}

void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( NodeManager const & GEOSX_UNUSED_ARG( nodeManager ),
                                                                    FaceManager const & GEOSX_UNUSED_ARG(facemanager) )
{
  // loop over the elements
  forall_in_range<serialPolicy>( 0, this->size(), GEOSX_LAMBDA ( localIndex const k )
    {
      m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
    });
}

void EmbeddedSurfaceSubRegion::AddNewEmbeddedSurface (localIndex const cellIndex,
                                                      R1Tensor normalVector)
{
  m_embeddedSurfaceToCell.push_back(cellIndex);
  m_normalVector.push_back(normalVector);

  // resize
  this->resize(this->size() + 1);

}

void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities(NodeManager const & nodeManager,
                                                                   EdgeManager const & edgeManager,
                                                                   FixedOneToManyRelation const & cellToEdges,
                                                                   R1Tensor origin)
{
  /*
   * Compute area of each embedded surface element:
   *
   * To compute the area of each embedded surface I need to find the intersection point
   * between the plane and each edge of each fractured element. So, I perform a loop over
   * all the embeddedSurfaceElements, then, for each of them I loop over the edges of the
   * relative fractured element (cut by the ) embedded surface) contained in m_embeddedSurfaceToCell.
   *
   * For each edge:
   *
   * 1. I check if it is cut by the plane using the Dot product between the distance of each node
   * from the origin and the normal vector. If an edgde is cut by the plane it is just a
   * matter of finding the intersection between a line and a plane.
   *
   * 2. Once I have the intersection points computing the area is easy as long as
   * they are ordered either CW or CCW.
   */

  array1d<R1Tensor> const & nodesCoord = nodeManager.referencePosition();
  array2d<localIndex> const & edgeToNodes = edgeManager.nodeList();
  localIndex count, edgeIndex;
  R1Tensor lineDir, dist;
  real64 prodScalarProd;

  // loop over the embedded surface elements
  for ( localIndex k=0; k < this->size(); k++ )
   {
    count = 0;
    array1d<R1Tensor> intersectionPoints;
    // loop over the edges
    std::cout << "embedded elem " << k << std::endl;
    for (localIndex ke = 0; ke < 12; ke++)
    {
      edgeIndex = cellToEdges[m_embeddedSurfaceToCell[k]][ke];
      dist = nodesCoord[edgeToNodes[edgeIndex][0]];
      dist -= origin;
      prodScalarProd = Dot(dist, m_normalVector[k]);
      dist = nodesCoord[edgeToNodes[edgeIndex][1]];
      dist -= origin;
      prodScalarProd *= Dot(dist, m_normalVector[k]);

      if (prodScalarProd < 0)
      {
        lineDir  = nodesCoord[edgeToNodes[edgeIndex][0]];
        lineDir -= nodesCoord[edgeToNodes[edgeIndex][1]];
        lineDir.Normalize();
        intersectionPoints.push_back(computationalGeometry::LinePlaneIntersection(lineDir, nodesCoord[edgeToNodes[edgeIndex][0]],
            m_normalVector[k], origin));

        m_elementCenter[k] += intersectionPoints[count];
        count++;
      }
    } //end of edge loop

    m_elementArea[k]   = computationalGeometry::ComputeSurfaceArea(intersectionPoints, count, m_normalVector[k]);
    m_elementCenter[k] /= count;
    this->CalculateElementGeometricQuantities(k);

    // std::cout << "Area = " << m_elementArea[k] << std::endl;
   } // end of embedded element loop
}


void EmbeddedSurfaceSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
}

} /* namespace geosx */



