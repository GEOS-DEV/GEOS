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
  m_unmappedGlobalIndicesInToNodes(),
  m_unmappedGlobalIndicesInToEdges(),
  m_unmappedGlobalIndicesInToFaces(),
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


R1Tensor const & EmbeddedSurfaceSubRegion::calculateElementCenter( localIndex k,
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

void EmbeddedSurfaceSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
}


void EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( localIndex const k)
{
  // Matteo: needs to be filled in with the proper computation.
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

  //for (localIndex k =0; k < this->size(); k++){
  //std::cout << m_embeddedSurfaceToCell[k] << std::endl;
  //}
  // add the cellIndex and the normalVector

  // m_embeddedSurfaceToCell[this->size() - 1] = cellIndex;
  // m_normalVector[this->size() - 1]          = normalVector;
}

void EmbeddedSurfaceSubRegion::ComputeElementArea(NodeManager const & nodeManager,
                                                  EdgeManager const & edgeManager,
                                                  FixedOneToManyRelation const & cellToEdges,
                                                  R1Tensor origin)
{
  array1d<R1Tensor> const & nodesCoord = nodeManager.referencePosition();
  array2d<localIndex> const & edgeToNodes = edgeManager.nodeList();
  localIndex count;
  localIndex edgeIndex;
  R1Tensor lineDir;

  // loop over the elements
  for ( localIndex k=0; k < this->size(); k++ )
   {
    count = 0;
    array1d<R1Tensor> intersectionPoints; // I am assuming there are 4 intersections (needs to be fixed)
    for (localIndex ke = 0; ke < 12; ke++)
    {
      edgeIndex = cellToEdges[k][ke];
      real64 prodScalarProd = Dot(nodesCoord[edgeToNodes[edgeIndex][0]], m_normalVector[k])
                    * Dot(nodesCoord[edgeToNodes[edgeIndex][1]], m_normalVector[k]);
      if (prodScalarProd < 0)
      {
        lineDir  = nodesCoord[edgeToNodes[edgeIndex][0]];
        lineDir -= nodesCoord[edgeToNodes[edgeIndex][1]];
        lineDir.Normalize();
        intersectionPoints.push_back(computationalGeometry::LinePlaneIntersection(lineDir, nodesCoord[edgeToNodes[edgeIndex][0]],
            m_normalVector[k], origin));

        R1Tensor a;
        a = intersectionPoints[count];
        a.print(std::cout);
        count += 1;
      }
    }
   }
}


localIndex EmbeddedSurfaceSubRegion::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex EmbeddedSurfaceSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template<bool DOPACK>
localIndex EmbeddedSurfaceSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                        arrayView1d<localIndex const> const & packList ) const
{

  localIndex packedSize = 0;
  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::nodeListString) );

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toNodesRelation.Base(),
                                         m_unmappedGlobalIndicesInToNodes,
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_toNodesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::edgeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toEdgesRelation.Base(),
                                         m_unmappedGlobalIndicesInToEdges,
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_toEdgesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::faceListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toFacesRelation.Base(),
                                         m_unmappedGlobalIndicesInToFaces,
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_toFacesRelation.RelatedObjectLocalToGlobal() );
  return packedSize;
}



localIndex EmbeddedSurfaceSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const GEOSX_UNUSED_ARG( overwriteUpMaps   ),
                                                   bool const GEOSX_UNUSED_ARG( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;
  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->m_globalToLocalMap,
                                     m_toNodesRelation.RelatedObjectGlobalToLocal() );


  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->m_globalToLocalMap,
                                     m_toEdgesRelation.RelatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation.Base(),
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->m_globalToLocalMap,
                                     m_toFacesRelation.RelatedObjectGlobalToLocal() );
  return unPackedSize;
}

void EmbeddedSurfaceSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );

}

void EmbeddedSurfaceSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
  //exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceElementsToCellRegionsString));
  //exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceElementsToCellSubRegionsString));
  //exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceElementsToCellIndexString));
}

} /* namespace geosx */



