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
#include "BufferOps.hpp"

namespace geosx
{
using namespace dataRepository;


FaceElementSubRegion::FaceElementSubRegion( string const & name,
                                      dataRepository::ManagedGroup * const parent ):
  ElementSubRegionBase( name, parent ),
  m_unmappedGlobalIndicesInToNodes(),
  m_unmappedGlobalIndicesInToEdges(),
  m_unmappedGlobalIndicesInToFaces(),
  m_faceElementsToCells(),
  m_newFaceElements(),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_toFacesRelation(),
  m_elementAperture(),
  m_elementArea()
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
    setPlotLevel(dataRepository::PlotLevel::LEVEL_0)->
    setDescription("The aperture of each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::elementAreaString, &m_elementArea, false )->
    setApplyDefaultValue(-1.0)->
    setDescription("The area of each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::elementCenterString, &m_elementCenter, false )->
    setDescription("The center of each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::elementVolumeString, &m_elementVolume, false )->
    setApplyDefaultValue(-1.0)->
    setPlotLevel(dataRepository::PlotLevel::LEVEL_0)->
    setDescription("The volume of each FaceElement.");

  RegisterViewWrapper( viewKeyStruct::faceElementsToCellRegionsString,
                       &(m_faceElementsToCells.m_toElementRegion), 0 )->
    setApplyDefaultValue(-1)->
    setPlotLevel(PlotLevel::NOPLOT)->
    setDescription( "A map of face element local indices to the cell local indices");

  RegisterViewWrapper( viewKeyStruct::faceElementsToCellSubRegionsString,
                       &(m_faceElementsToCells.m_toElementSubRegion), 0 )->
    setApplyDefaultValue(-1)->
    setPlotLevel(PlotLevel::NOPLOT)->
    setDescription( "A map of face element local indices to the cell local indices");

  RegisterViewWrapper( viewKeyStruct::faceElementsToCellIndexString,
                       &(m_faceElementsToCells.m_toElementIndex), 0 )->
    setApplyDefaultValue(-1)->
    setPlotLevel(PlotLevel::NOPLOT)->
    setDescription( "A map of face element local indices to the cell local indices");

  m_faceElementsToCells.resize(0,2);
  m_faceElementsToCells.setElementRegionManager( getParent()->getParent()->getParent()->getParent()->group_cast<ElementRegionManager*>() );

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

  forall_in_range<serialPolicy>( 0, this->size(), GEOSX_LAMBDA ( localIndex const k )
  {
    m_elementArea[k] = faceArea[ m_toFacesRelation[k][0] ];
    m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
  });
}




localIndex FaceElementSubRegion::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex FaceElementSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template<bool DOPACK>
localIndex FaceElementSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
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


  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::faceElementsToCellRegionsString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         this->m_faceElementsToCells,
                                         packList,
                                         m_faceElementsToCells.getElementRegionManager() );

  return packedSize;
}



localIndex FaceElementSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const overwriteUpMaps,
                                                   bool const overwriteDownMaps )
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

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::faceElementsToCellRegionsString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_faceElementsToCells,
                                     packList,
                                     m_faceElementsToCells.getElementRegionManager(),
                                     overwriteUpMaps );



  return unPackedSize;
}

void FaceElementSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
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

void FaceElementSubRegion::inheritGhostRankFromParentFace( FaceManager const * const faceManager,
                                                           std::set<localIndex> const & indices )
{
  for( localIndex const & index : indices )
  {
    m_ghostRank[index] = faceManager->m_ghostRank[ m_toFacesRelation[index][0] ];
  }
}

void FaceElementSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceElementsToCellRegionsString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceElementsToCellSubRegionsString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceElementsToCellIndexString));
}

} /* namespace geosx */
