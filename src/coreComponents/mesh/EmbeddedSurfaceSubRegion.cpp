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
 * @file EmbeddedSurfaceSubRegion.cpp
 */

#include "EmbeddedSurfaceSubRegion.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "common/MpiWrapper.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
using namespace dataRepository;


void surfaceWithGhostNodes::insert( globalIndex const & edgeIndex )
{
  parentEdgeIndex.push_back( edgeIndex );
  numOfNodes += 1;
}

EmbeddedSurfaceSubRegion::EmbeddedSurfaceSubRegion( string const & name,
                                                    dataRepository::Group * const parent ):
  SurfaceElementSubRegion( name, parent ),
  m_normalVector(),
  m_tangentVector1(),
  m_tangentVector2(),
  m_numOfJumpEnrichments( 3 ),
  m_connectivityIndex(),
  m_parentPlaneName()
{
  m_elementType = ElementType::Polygon;

  registerWrapper( viewKeyStruct::normalVectorString(), &m_normalVector ).
    setDescription( "Unit normal vector to the embedded surface." );

  registerWrapper( viewKeyStruct::t1VectorString(), &m_tangentVector1 ).
    setDescription( "Unit vector in the first tangent direction to the embedded surface." );

  registerWrapper( viewKeyStruct::t2VectorString(), &m_tangentVector2 ).
    setDescription( "Unit vector in the second tangent direction to the embedded surface." );

  registerWrapper( viewKeyStruct::elementCenterString(), &m_elementCenter ).
    setDescription( "The center of each EmbeddedSurface element." );

  registerWrapper( viewKeyStruct::elementVolumeString(), &m_elementVolume ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_0 ).
    setDescription( "The volume of each EmbeddedSurface element." );

  registerWrapper( viewKeyStruct::connectivityIndexString(), &m_connectivityIndex ).
    setApplyDefaultValue( 1 ).
    setDescription( "Connectivity index of each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::surfaceElementToParentPlaneString(), &m_parentPlaneName ).
    setDescription( "A map of surface element to the parent fracture name" );

  m_normalVector.resizeDimension< 1 >( 3 );
  m_tangentVector1.resizeDimension< 1 >( 3 );
  m_tangentVector2.resizeDimension< 1 >( 3 );
  m_surfaceElementsToCells.resize( 0, 1 );
}


EmbeddedSurfaceSubRegion::~EmbeddedSurfaceSubRegion()
{}

void EmbeddedSurfaceSubRegion::calculateElementGeometricQuantities( NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
                                                                    FaceManager const & GEOSX_UNUSED_PARAM( facemanager ) )
{
  // loop over the elements
  forAll< serialPolicy >( this->size(), [=] ( localIndex const k )
  {
    m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
  } );
}

void EmbeddedSurfaceSubRegion::calculateElementGeometricQuantities( arrayView2d< real64 const > const intersectionPoints,
                                                                    localIndex const k )
{
  for( localIndex p = 0; p < intersectionPoints.size( 0 ); p++ )
  {
    LvArray::tensorOps::add< 3 >( m_elementCenter[k], intersectionPoints[ p ] );
  }

  // update area
  m_elementArea[ k ] = computationalGeometry::ComputeSurfaceArea( intersectionPoints, m_normalVector[k] );

  LvArray::tensorOps::scale< 3 >( m_elementCenter[ k ], 1.0 / intersectionPoints.size( 0 ) );

  // update volume
  m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
}

void EmbeddedSurfaceSubRegion::computeConnectivityIndex( localIndex const k,
                                                         arrayView2d< localIndex const, cells::NODE_MAP_USD > const cellToNodes,
                                                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodesCoord )
{
  // 1. Compute average distance
  // TODO: This is a pretty bad approximation of the average distance and proper numerical integration should
  // be implemented.
  real64 averageDistance = 0.0;

  localIndex const cellIndex = m_surfaceElementsToCells.m_toElementIndex[k][0];
  localIndex const numOfNodes = cellToNodes.size( 1 );

  real64 nodeToFracCenter[3];
  for( localIndex a=0; a < numOfNodes; a++ )
  {
    localIndex const nodeIndex = cellToNodes[cellIndex][a];
    LvArray::tensorOps::copy< 3 >( nodeToFracCenter, nodesCoord[nodeIndex] );
    LvArray::tensorOps::subtract< 3 >( nodeToFracCenter, m_elementCenter[k] );
    real64 distance = LvArray::tensorOps::AiBi< 3 >( nodeToFracCenter, m_normalVector[k] );
    averageDistance += std::sqrt( distance * distance );
  }
  averageDistance /= numOfNodes;

  //2. Compute connectivity index
  m_connectivityIndex[k] = m_elementArea[ k ] / averageDistance;
}

bool EmbeddedSurfaceSubRegion::addNewEmbeddedSurface ( localIndex const cellIndex,
                                                       localIndex const subRegionIndex,
                                                       localIndex const regionIndex,
                                                       NodeManager const & nodeManager,
                                                       EmbeddedSurfaceNodeManager & embSurfNodeManager,
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
  arrayView2d< localIndex const > const edgeToNodes = edgeManager.nodeList();
  arrayView1d< integer const > const edgeGhostRank = edgeManager.ghostRank();
  arrayView1d< globalIndex const > const & edgeLocalToGlobal = edgeManager.localToGlobalMap();

  R1Tensor origin        = fracture->getCenter();
  R1Tensor normalVector  = fracture->getNormal();
  localIndex edgeIndex;
  real64 lineDir[3], dist[3], point[3], distance[3], prodScalarProd;

  array2d< real64 > intersectionPoints( 0, 3 );
  array1d< integer > pointGhostRank;
  array1d< localIndex > pointParentIndex;

  localIndex numPoints = 0;
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
      LvArray::tensorOps::copy< 3 >( lineDir, nodesCoord[edgeToNodes[edgeIndex][0]] );
      LvArray::tensorOps::subtract< 3 >( lineDir, nodesCoord[edgeToNodes[edgeIndex][1]] );
      LvArray::tensorOps::normalize< 3 >( lineDir );
      //find the intersection point
      computationalGeometry::LinePlaneIntersection( lineDir,
                                                    nodesCoord[edgeToNodes[edgeIndex][0]],
                                                    normalVector,
                                                    origin,
                                                    point );

      // Check if the point is inside the fracture (bounded plane)
      if( !(fracture->isCoordInObject( point )) )
      {
        addEmbeddedElem = false;
      }
      intersectionPoints.resizeDimension< 0 >( numPoints+1 );
      pointGhostRank.resize( numPoints+1 );
      pointParentIndex.resize( numPoints+1 );
      pointGhostRank[numPoints] = edgeGhostRank[edgeIndex];
      pointParentIndex[numPoints] = edgeIndex;
      LvArray::tensorOps::copy< 3 >( intersectionPoints[numPoints], point );
      numPoints++;
    }
  } //end of edge loop


  if( addEmbeddedElem && intersectionPoints.size( 0 ) > 0 )
  {
    // resize
    localIndex surfaceIndex = this->size();
    this->resize( surfaceIndex + 1 );

    // Reorder the points CCW and then add the point to the list in the nodeManager if it is a new one.
    array1d< int > originalIndices = computationalGeometry::orderPointsCCW( intersectionPoints, normalVector );

    // Get location of embedded surfaces nodes.
    array2d< real64, nodes::REFERENCE_POSITION_PERM > & embSurfNodesPos =
      embSurfNodeManager.referencePosition();

    // fill out elemNodes array with the previously found intersection points
    // add new nodes to embSurfNodes
    bool isNew;
    bool hasGhostNode = false;
    localIndex nodeIndex;
    array1d< localIndex > elemNodes( intersectionPoints.size( 0 ) );

    for( localIndex j=0; j < intersectionPoints.size( 0 ); j++ )
    {
      isNew = true;
      for( localIndex h=0; h < embSurfNodeManager.size(); h++ )
      {
        LvArray::tensorOps::copy< 3 >( distance, intersectionPoints[ j ] );
        LvArray::tensorOps::subtract< 3 >( distance, embSurfNodesPos[ h ] );
        if( LvArray::tensorOps::l2Norm< 3 >( distance ) < 1e-9 )
        {
          isNew = false;
          nodeIndex = h;
          break;
        }
      }
      if( isNew )
      {
        // Add the point to the node Manager if it's not a ghost
        globalIndex parentEdgeID = edgeLocalToGlobal[ pointParentIndex[ originalIndices[ j ] ] ];
        nodeIndex = embSurfNodeManager.size();

        if( pointGhostRank[ originalIndices[ j ] ] < 0 )
        {
          embSurfNodeManager.appendNode( intersectionPoints[ j ],
                                         pointGhostRank[ originalIndices[ j ] ] );

          arrayView1d< localIndex > const & parentIndex =
            embSurfNodeManager.getExtrinsicData< extrinsicMeshData::ParentEdgeIndex >();

          parentIndex[nodeIndex] = pointParentIndex[ originalIndices[ j ] ];

          array1d< globalIndex > & parentEdgeGlobalIndex = embSurfNodeManager.getParentEdgeGlobalIndex();
          parentEdgeGlobalIndex[nodeIndex] = parentEdgeID;
        }
        else
        {
          hasGhostNode = true;
        }
      }
      elemNodes[ j ] =  nodeIndex;
    }

    if( hasGhostNode )
    {
      localIndex surfWithGhostsIndex = m_surfaceWithGhostNodes.size();
      m_surfaceWithGhostNodes.resize( surfWithGhostsIndex + 1 );
      m_surfaceWithGhostNodes[ surfWithGhostsIndex ].surfaceIndex = surfaceIndex;
      for( int ii = 0; ii < elemNodes.size(); ii++ )
      {
        globalIndex parentEdgeID = edgeLocalToGlobal[ pointParentIndex[ originalIndices[ ii ] ] ];
        m_surfaceWithGhostNodes[ surfWithGhostsIndex ].insert( parentEdgeID );
      }
    }

    m_toNodesRelation.resizeArray( surfaceIndex, elemNodes.size() );
    for( localIndex inode = 0; inode < elemNodes.size(); inode++ )
    {
      m_toNodesRelation( surfaceIndex, inode ) = elemNodes[ inode ];
    }

    m_surfaceElementsToCells.m_toElementIndex[ surfaceIndex ][0]     = cellIndex;
    m_surfaceElementsToCells.m_toElementSubRegion[ surfaceIndex ][0] =  subRegionIndex;
    m_surfaceElementsToCells.m_toElementRegion[ surfaceIndex ][0]    =  regionIndex;
    m_parentPlaneName[ surfaceIndex ] = fracture->getName();
    LvArray::tensorOps::copy< 3 >( m_normalVector[ surfaceIndex ], normalVector );
    LvArray::tensorOps::copy< 3 >( m_tangentVector1[ surfaceIndex ], fracture->getWidthVector());
    LvArray::tensorOps::copy< 3 >( m_tangentVector2[ surfaceIndex ], fracture->getLengthVector());
    this->calculateElementGeometricQuantities( intersectionPoints.toViewConst(), this->size()-1 );
  }
  return addEmbeddedElem;
}

void EmbeddedSurfaceSubRegion::inheritGhostRank( array1d< array1d< arrayView1d< integer const > > > const & cellGhostRank )
{
  arrayView1d< integer > const & ghostRank = this->ghostRank();
  for( localIndex k=0; k < size(); ++k )
  {
    localIndex regionIndex    = m_surfaceElementsToCells.m_toElementRegion[k][0];
    localIndex subRegionIndex = m_surfaceElementsToCells.m_toElementSubRegion[k][0];
    localIndex cellIndex      = m_surfaceElementsToCells.m_toElementIndex[k][0];

    ghostRank[k] = cellGhostRank[regionIndex][subRegionIndex][cellIndex];
  }
}

void EmbeddedSurfaceSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getEmbSurfNodeManager() );
}

localIndex EmbeddedSurfaceSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}

localIndex EmbeddedSurfaceSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex EmbeddedSurfaceSubRegion::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                            arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           nodeList().base().toViewConst(),
                                           m_unmappedGlobalIndicesInToNodes,
                                           packList,
                                           localToGlobal,
                                           nodeLocalToGlobal );


  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::surfaceElementsToCellRegionsString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           this->m_surfaceElementsToCells,
                                           packList,
                                           m_surfaceElementsToCells.getElementRegionManager() );

  return packedSize;
}


localIndex EmbeddedSurfaceSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                       localIndex_array & packList,
                                                       bool const overwriteUpMaps,
                                                       bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOSX_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOSX_ERROR_IF_NE( elementListString, viewKeyStruct::surfaceElementsToCellRegionsString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_surfaceElementsToCells,
                                     packList.toViewConst(),
                                     m_surfaceElementsToCells.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void EmbeddedSurfaceSubRegion::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::viewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::surfaceElementsToCellRegionsString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::surfaceElementsToCellSubRegionsString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::surfaceElementsToCellIndexString() ));
  //exclusionList.insert( this->getWrapperIndex( viewKeyStruct::surfaceWithGhostNodesString() ));
}

} /* namespace geosx */
