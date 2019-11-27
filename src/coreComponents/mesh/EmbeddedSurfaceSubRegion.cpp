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
  m_elementAperture(),
  m_elementArea()
{
  registerWrapper( viewKeyStruct::nodeListString, &m_toNodesRelation, false )->
    setDescription("Map to the nodes attached to each EmbeddedSurface.");

  registerWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation, false )->
    setDescription("Map to the edges.");

  registerWrapper( viewKeyStruct::cellListString, &m_embeddedSurfaceToCell, false )->
    setDescription("Map to the cells.");

  registerWrapper( viewKeyStruct::normalVectorString, &m_normalVector, false )->
    setDescription("Unit normal vector to the embedded surface.");

  registerWrapper( viewKeyStruct::elementApertureString, &m_elementAperture, false )->
    setApplyDefaultValue(1.0e-5)->
    setPlotLevel(dataRepository::PlotLevel::LEVEL_0)->
    setDescription("The aperture of each EmbeddedSurface.");

  registerWrapper( viewKeyStruct::elementAreaString, &m_elementArea, false )->
    setApplyDefaultValue(-1.0)->
    setDescription("The area of each EmbeddedSurface element.");

  registerWrapper( viewKeyStruct::elementCenterString, &m_elementCenter, false )->
    setDescription("The center of each EmbeddedSurface element.");

  registerWrapper( viewKeyStruct::elementVolumeString, &m_elementVolume, false )->
    setApplyDefaultValue(-1.0)->
    setPlotLevel(dataRepository::PlotLevel::LEVEL_0)->
    setDescription("The volume of each EmbeddedSurface element.");

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

void EmbeddedSurfaceSubRegion::AddNewEmbeddedSurface (localIndex const cellIndex,
                                                      R1Tensor normalVector,
                                                      NodeManager const & nodeManager,
                                                      EdgeManager const & edgeManager,
                                                      FixedOneToManyRelation const & cellToEdges,
                                                      BoundedPlane const * fracture)
{
  /* The goal is to add an embeddedSurfaceElem if it is contained within the BoundedPlane
   *
   * A. Identify whether the cell falls within the bounded plane or not
   *
   * we loop over each edge:
   *
   *   1. check if it is cut by the plane using the Dot product between the distance of each node
   *   from the origin and the normal vector.
   *   2. If an edge is cut by the plane we look for the intersection between a line and a plane.
   *   3. Once we have the intersection we check if it falls inside the bounded plane.
   *
   * Only elements for which all intersection points fall within the fracture plane limits will be added.
   * If the frac does not cut through the entire element we will just chop it (it's a discretization error).
   *
   * B. Once we know the element has to be added we compute it's geometric properties:
   * - Surface Area: this is trivial given the intersection points as long as they are ordered either CW or CCW.
   * - centre:
   * - Volume:
   * - Heaviside:
   */

  std::cout << "cell " << cellIndex  << std::endl;

  bool addEmbeddedElem = true;
  array1d<R1Tensor> const & nodesCoord = nodeManager.referencePosition();
  array2d<localIndex> const & edgeToNodes = edgeManager.nodeList();
  R1Tensor origin  = fracture->getCenter();
  localIndex edgeIndex;
  R1Tensor lineDir, dist, point;
  real64 prodScalarProd;

  array1d<R1Tensor> intersectionPoints;
  for (localIndex ke = 0; ke < 12; ke++)
  {
    edgeIndex = cellToEdges[cellIndex][ke];
    dist = nodesCoord[edgeToNodes[edgeIndex][0]];
    dist -= origin;
    prodScalarProd = Dot(dist, normalVector);
    dist = nodesCoord[edgeToNodes[edgeIndex][1]];
    dist -= origin;
    prodScalarProd *= Dot(dist, normalVector);

    if (prodScalarProd < 0)
    {
      std::cout << "node 1: " << nodesCoord[edgeToNodes[edgeIndex][0]] <<  " node 2: " << nodesCoord[edgeToNodes[edgeIndex][1]] << std::endl;
      lineDir  = nodesCoord[edgeToNodes[edgeIndex][0]];
      lineDir -= nodesCoord[edgeToNodes[edgeIndex][1]];
      lineDir.Normalize();
      point = computationalGeometry::LinePlaneIntersection(lineDir,
                                                           nodesCoord[edgeToNodes[edgeIndex][0]],
                                                           normalVector,
                                                           origin);
      std::cout << "origin " << origin <<  " point " << point << std::endl;
      // Check if the point is inside the fracture (bounded plane)
      if ( !fracture->IsCoordInObject(point) )
      {
        addEmbeddedElem = false;
      }
      intersectionPoints.push_back(point);
    }
  } //end of edge loop

  if (addEmbeddedElem)
  {
    std::cout << "adding cell " << cellIndex  << std::endl;
    m_embeddedSurfaceToCell.push_back(cellIndex);
    m_normalVector.push_back(normalVector);
    // resize
    this->resize(this->size() + 1);
    this->CalculateElementGeometricQuantities(intersectionPoints, this->size()-1);
  }
}

void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( array1d<R1Tensor> const  intersectionPoints,
                                                                    localIndex k )
{
      for (localIndex p = 0; p < intersectionPoints.size(); p++)
      {
        m_elementCenter[k] += intersectionPoints[p];
      }

    m_elementArea[k]   = computationalGeometry::ComputeSurfaceArea(intersectionPoints, intersectionPoints.size(), m_normalVector[k]);
    m_elementCenter[k] /= intersectionPoints.size();
    this->CalculateElementGeometricQuantities(k);
}

void EmbeddedSurfaceSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
}

} /* namespace geosx */



