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
  m_tangentVector1(),
  m_tangentVector2(),
  m_embeddedSurfaceToRegion(),
  m_embeddedSurfaceToSubRegion(),
  m_embeddedSurfaceToCell(),
  m_toNodesRelation(),
  m_elementAperture(),
  m_elementArea(),
  m_numNodesPerSurface(),
  m_numOfJumpEnrichments( 3 ),
  m_connectivityIndex()
{
  registerWrapper( viewKeyStruct::regionListString, &m_embeddedSurfaceToRegion )->
    setDescription( "Map to the region cut by each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::subregionListString, &m_embeddedSurfaceToSubRegion )->
    setDescription( "Map to the subregion cut by each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::nodeListString, &m_toNodesRelation )->
    setDescription( "Map to the nodes attached to each EmbeddedSurface." );
//
//  registerWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation, false )->
//    setDescription( "Map to the edges." );

  registerWrapper( viewKeyStruct::cellListString, &m_embeddedSurfaceToCell )->
    setDescription( "Map to the cells." );

  registerWrapper( viewKeyStruct::normalVectorString, &m_normalVector )->
    setDescription( "Unit normal vector to the embedded surface." );

  registerWrapper( viewKeyStruct::t1VectorString, &m_tangentVector1 )->
    setDescription( "Unit vector in the first tangent direction to the embedded surface." );

  registerWrapper( viewKeyStruct::t2VectorString, &m_tangentVector2 )->
    setDescription( "Unit vector in the second tangent direction to the embedded surface." );

  registerWrapper( viewKeyStruct::elementApertureString, &m_elementAperture )->
    setApplyDefaultValue( 1.0e-5 )->
    setPlotLevel( dataRepository::PlotLevel::LEVEL_0 )->
    setDescription( "The aperture of each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::elementAreaString, &m_elementArea )->
    setApplyDefaultValue( -1.0 )->
    setDescription( "The area of each EmbeddedSurface element." );

  registerWrapper( viewKeyStruct::elementCenterString, &m_elementCenter )->
    setDescription( "The center of each EmbeddedSurface element." );

  registerWrapper( viewKeyStruct::elementVolumeString, &m_elementVolume )->
    setApplyDefaultValue( -1.0 )->
    setPlotLevel( dataRepository::PlotLevel::LEVEL_0 )->
    setDescription( "The volume of each EmbeddedSurface element." );

  registerWrapper( viewKeyStruct::numNodesString, &m_numNodesPerSurface )->
    setApplyDefaultValue( 0 )->
    setDescription( "Number of nodes of each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::connectivityIndexString, &m_connectivityIndex )->
    setApplyDefaultValue( 1 )->
    setDescription( "Connectivity index of each EmbeddedSurface." );
}


EmbeddedSurfaceSubRegion::~EmbeddedSurfaceSubRegion()
{}


void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( localIndex const k )
{
  m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
}

void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
                                                                    FaceManager const & GEOSX_UNUSED_PARAM( facemanager ) )
{
  // loop over the elements
  forAll< serialPolicy >( this->size(), [=] ( localIndex const k )
  {
    m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
  } );
}

void EmbeddedSurfaceSubRegion::AddNewEmbeddedSurface ( localIndex const cellIndex,
                                                       R1Tensor normalVector )
{
  m_embeddedSurfaceToCell.emplace_back( cellIndex );
  m_normalVector.emplace_back( normalVector );

  // resize
  this->resize( this->size() + 1 );

}

bool EmbeddedSurfaceSubRegion::AddNewEmbeddedSurface ( localIndex const cellIndex,
                                                       localIndex const regionIndex,
                                                       localIndex const subRegionIndex,
                                                       NodeManager & nodeManager,
                                                       EdgeManager const & edgeManager,
                                                       FixedOneToManyRelation const & cellToEdges,
                                                       BoundedPlane const * fracture )
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
   * If the fracture does not cut through the entire element we will just chop it (it's a discretization error).
   *
   * B. Once we know the element has to be added we compute it's geometric properties:
   * - Surface Area: this is trivial given the intersection points as long as they are ordered either CW or CCW.
   * - centre:
   * - Volume:
   */

  bool addEmbeddedElem = true;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager.referencePosition();
  EdgeManager::NodeMapType::ViewTypeConst const & edgeToNodes = edgeManager.nodeList();
  R1Tensor origin        = fracture->getCenter();
  R1Tensor normalVector  = fracture->getNormal();
  localIndex edgeIndex;
  R1Tensor lineDir, dist, point, distance;
  real64 prodScalarProd;

  array1d< R1Tensor > intersectionPoints;
  for( localIndex ke = 0; ke < cellToEdges.size( 1 ); ke++ )
  {
    edgeIndex = cellToEdges[cellIndex][ke];
    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][0]] );
    LvArray::tensorOps::subtract< 3 >( dist, origin );
    prodScalarProd = LvArray::tensorOps::AiBi< 3 >( dist, normalVector );
    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][1]] );
    LvArray::tensorOps::subtract< 3 >( dist, origin );
    prodScalarProd *= LvArray::tensorOps::AiBi< 3 >( dist, normalVector );

    // check if the plane intersects the edge
    if( prodScalarProd < 0 )
    {
      lineDir = LVARRAY_TENSOROPS_INIT_LOCAL_3( nodesCoord[edgeToNodes[edgeIndex][0]] );
      LvArray::tensorOps::subtract< 3 >( lineDir, nodesCoord[edgeToNodes[edgeIndex][1]] );
      LvArray::tensorOps::normalize< 3 >( lineDir );
      //find the intersection point
      point = computationalGeometry::LinePlaneIntersection( lineDir,
                                                            nodesCoord[edgeToNodes[edgeIndex][0]],
                                                            normalVector,
                                                            origin );

      // Check if the point is inside the fracture (bounded plane)
      if( !fracture->IsCoordInObject( point ) )
      {
        addEmbeddedElem = false;
      }
      intersectionPoints.emplace_back( point );
    }
  } //end of edge loop

  if( addEmbeddedElem )
  {

    // resize
    localIndex surfaceIndex = this->size();
    this->resize( surfaceIndex + 1 );
    m_numNodesPerSurface[surfaceIndex] = intersectionPoints.size();

    // Reorder the points CCW and then add the point to the list in the nodeManager if it is a new one.
    intersectionPoints = computationalGeometry::orderPointsCCW( intersectionPoints, intersectionPoints.size(), normalVector );
    array2d< real64, nodes::REFERENCE_POSITION_PERM > & embSurfNodesPos = nodeManager.embSurfNodesPosition();

    bool isNew;
    localIndex NodeIndex;
    array1d< localIndex > elemNodes( m_numNodesPerSurface[surfaceIndex] );

    for( localIndex j=0; j < m_numNodesPerSurface[surfaceIndex]; j++ )
    {
      isNew = true;
      for( localIndex h=0; h < embSurfNodesPos.size( 0 ); h++ )
      {
        LvArray::tensorOps::copy< 3 >( distance, intersectionPoints[j] );
        LvArray::tensorOps::subtract< 3 >( distance, embSurfNodesPos.toViewConst()[h] );
        if( distance.L2_Norm() < 1e-9 )
        {
          isNew = false;
          NodeIndex = h;
          break;
        }
      }
      if( isNew )
      {
        // Add the point to the
        NodeIndex = embSurfNodesPos.size( 0 );
        embSurfNodesPos.resize( NodeIndex + 1 );
        LvArray::tensorOps::copy< 3 >( embSurfNodesPos[NodeIndex], intersectionPoints[j] );
      }
      elemNodes[j] =  NodeIndex;
    }

    m_toNodesRelation.resizeArray( surfaceIndex, m_numNodesPerSurface[surfaceIndex] );
    for( localIndex inode = 0; inode < m_numNodesPerSurface[surfaceIndex]; inode++ )
    {
      m_toNodesRelation( surfaceIndex, inode ) = elemNodes[inode];
    }

    m_embeddedSurfaceToCell[ surfaceIndex ]      = cellIndex;
    m_embeddedSurfaceToRegion[ surfaceIndex ]    =  regionIndex;
    m_embeddedSurfaceToSubRegion[ surfaceIndex ] =  subRegionIndex;
    LvArray::tensorOps::copy< 3 >( m_normalVector[ surfaceIndex ], normalVector );
    LvArray::tensorOps::copy< 3 >( m_tangentVector1[ surfaceIndex ], fracture->getWidthVector());
    LvArray::tensorOps::copy< 3 >( m_tangentVector2[ surfaceIndex ], fracture->getLengthVector());
    this->CalculateElementGeometricQuantities( intersectionPoints, this->size()-1 );
  }
  return addEmbeddedElem;
}

void EmbeddedSurfaceSubRegion::inheritGhostRank( array1d< array1d< arrayView1d< integer const > > > const & cellGhostRank )
{
  arrayView1d< integer > const & ghostRank = this->ghostRank();
  for( localIndex k=0; k < size(); ++k )
  {
    ghostRank[k] = cellGhostRank[m_embeddedSurfaceToRegion[k]][m_embeddedSurfaceToSubRegion[k]][m_embeddedSurfaceToCell[k]];
  }
}

void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( array1d< R1Tensor > const intersectionPoints,
                                                                    localIndex const k )
{
  for( localIndex p = 0; p < intersectionPoints.size(); p++ )
  {
    LvArray::tensorOps::add< 3 >( m_elementCenter[k], intersectionPoints[ p ] );
  }

  m_elementArea[ k ] = computationalGeometry::ComputeSurfaceArea( intersectionPoints, intersectionPoints.size(), m_normalVector[k] );

  LvArray::tensorOps::scale< 3 >( m_elementCenter[ k ], 1.0 / intersectionPoints.size() );
  this->CalculateElementGeometricQuantities( k );
}

void EmbeddedSurfaceSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
}

int EmbeddedSurfaceSubRegion::totalNumberOfNodes() const
{
  int totalNumNodes = 0;
  for( localIndex esi=0; esi<size(); esi++ )
  {
    totalNumNodes += m_numNodesPerSurface[esi];
  }
  return totalNumNodes;
}

real64 EmbeddedSurfaceSubRegion::ComputeHeavisideFunction( ArraySlice< real64 const, 1, nodes::REFERENCE_POSITION_USD - 1 > const nodeCoord,
                                                           localIndex const k ) const
{
  real64 heaviside;
  R1Tensor distanceVector;
  LvArray::tensorOps::copy< 3 >( distanceVector, nodeCoord );
  LvArray::tensorOps::subtract< 3 >( distanceVector, m_elementCenter[k] );

  heaviside = LvArray::tensorOps::AiBi< 3 >( distanceVector, m_normalVector[k] ) > 0 ? 1 : 0;

  return heaviside;
}

} /* namespace geosx */
