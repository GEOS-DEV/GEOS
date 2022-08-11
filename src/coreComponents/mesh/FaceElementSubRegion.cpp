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
 * @file FaceElementSubRegion.cpp
 */

#include "FaceElementSubRegion.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/MpiWrapper.hpp"

namespace geosx
{
using namespace dataRepository;


FaceElementSubRegion::FaceElementSubRegion( string const & name,
                                            dataRepository::Group * const parent ):
  SurfaceElementSubRegion( name, parent ),
  m_unmappedGlobalIndicesInToEdges(),
  m_unmappedGlobalIndicesInToFaces(),
  m_newFaceElements(),
  m_toFacesRelation()
{
  m_elementType = ElementType::Hexahedron;

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();

  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation ).
    setDescription( "Map to the faces attached to each FaceElement." ).
    reference().resize( 0, 2 );

  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString(), &m_edgesToFractureConnectorsEdges ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of edge local indices to the fracture connector local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString(), &m_fractureConnectorsEdgesToEdges ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices to edge local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString(), &m_fractureConnectorEdgesToFaceElements ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices face element local indices" ).
    setSizedFromParent( 0 );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  registerWrapper( viewKeyStruct::separationCoeffString(), &m_separationCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
    setDescription( "Scalar indicator of level of separation for a fracturing face." );
#endif

  excludeWrappersFromPacking( { viewKeyStruct::faceListString() } );

  m_surfaceElementsToCells.resize( 0, 2 );

  m_numNodesPerElement = 8;
}

ArrayOfArrays< localIndex > convert( array2d< localIndex > const & vv )
{
  ArrayOfArrays< localIndex > res;

  for( localIndex i = 0; i < vv.size(0); ++i ){
    auto const & vvv = vv[i];
    res.appendArray( vvv.begin(), vvv.end() );
  }

  return res;
}

//array1d< localIndex > convert( std::vector< localIndex > const & v )
//{
//  array1d< localIndex > res;
//  for( auto const & val: v )
//  {
//    res.emplace_back( val );
//  }
//
//  return res;
//}

ArrayOfArrays< localIndex > convert( std::vector< std::vector< localIndex > > const & vv )
{
  ArrayOfArrays< localIndex > res;

  for( std::size_t i = 0; i < vv.size(); ++i )
  {
    auto const & vvv = vv[i];
    res.appendArray( vvv.begin(), vvv.end() );
  }

  return res;
}


void FaceElementSubRegion::copyFromCellBlock( CellBlockABC const & cellBlock, CellBlockManagerABC const & cellBlockManager )
{
  // TODO

  this->resize( cellBlock.numElements() );

  this->m_toNodesRelation.base() = convert( cellBlock.getElemToNodes() ); // Inconsistent with the dimensions of the line below: it's 10x8
//  this->m_toEdgesRelation.base() = convert( cellBlock.getElemToFaces() ); // Warning, not the proper dimension, it should be 10x4! Ou plutôt il faut cbm.getFaceToEdges() ?
//  std::vector< std::vector< localIndex > > tmp;
  auto & toEdges = this->m_toEdgesRelation.base();
  toEdges.resize( cellBlock.numElements(), 4 );
//  auto const & e2f = cellBlock.getElemToFaces();
  auto const & f2ed= cellBlockManager.getFaceToEdges();
  std::vector< localIndex > const leftCommonFaces{ 205, 209, 213, 217, 221, 225, 229, 233, 237, 241 };
  std::vector< localIndex > const rightCommonFaces{ 246, 250, 254, 258, 262, 266, 270, 274, 278, 282 };
  for( std::size_t i = 0; i < 10; ++i )
  {
    toEdges.resizeArray( i, 4 );
    for( std::size_t j = 0; j < 4; ++j )
    {
      toEdges[i][j] = f2ed[leftCommonFaces[i]][j];
    }
  }
//  for( std::size_t i = 0; i < 10; ++i )
//  {
//    localIndex const f = e2f[i][0];
//    auto const & edges = f2ed[f];
//    std::vector< localIndex > tmp2;
//    for( auto const & val: edges )
//    {
//      tmp2.push_back( val );
//    }
//    tmp.push_back( tmp2 );
//  }
//  this->m_toEdgesRelation.base() = convert( tmp );
  // TODO manipulate `faceBlock.getFaceToElements();` to feed `this->m_surfaceElementsToCells`.
  // Use `transformCellBlockToRegionMap` instead of the wrong code below.
//  ToCellRelation< array2d< localIndex > > const & f2e = cellBlock.getElemToFaces();
//  this->m_surfaceElementsToCells.m_toElementIndex = f2e.toCellIndex;
//  this->m_surfaceElementsToCells.m_toElementSubRegion = f2e.toBlockIndex;
//  this->m_surfaceElementsToCells.m_toElementRegion.setValues<serialPolicy>(0);

  auto & e = this->m_surfaceElementsToCells.m_toElementIndex ;
  auto & esr= this->m_surfaceElementsToCells.m_toElementSubRegion;
  auto & er= this->m_surfaceElementsToCells.m_toElementRegion;
//  this->m_surfaceElementsToCells.m_toElementRegion.setValues<serialPolicy>(0);
  for( localIndex i = 0; i < 10; ++i )
  {
    e[i][0] = 40 + i ;
    e[i][1] = 50 + i ;
    esr[i][0] = 0 ;
    esr[i][1] = 0 ;
    er[i][0] = 0 ;
    er[i][1] = 0 ;
  }

  for( localIndex i = 0; i < 10; ++i )
  {
//    this->m_toFacesRelation.base()[i][0] = 205 + 4 * i;
//    this->m_toFacesRelation.base()[i][1] = 246 + 4 * i;
    this->m_toFacesRelation.base()[i][0] = leftCommonFaces[i];
    this->m_toFacesRelation.base()[i][1] = rightCommonFaces[i];
  }
  for( localIndex i = 0; i < 10; ++i )
  {
    this->m_newFaceElements.insert(i);
  }

  // TODO what about the external fields?
}


void FaceElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

void FaceElementSubRegion::calculateSingleElementGeometricQuantities( localIndex const k,
                                                                      arrayView1d< real64 const > const & faceArea )
{
  m_elementArea[k] = faceArea[ m_toFacesRelation[k][0] ];
  m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
//  m_elementArea[k] = 1.;
//  m_elementVolume[k] = 1.e-5;
}

void FaceElementSubRegion::calculateElementGeometricQuantities( NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
                                                                FaceManager const & faceManager )
{
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();

  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateSingleElementGeometricQuantities( k, faceArea );
  } );
}



localIndex FaceElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}

localIndex FaceElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex FaceElementSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > const nodeLocalToGlobal = m_toNodesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const edgeLocalToGlobal = m_toEdgesRelation.relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > const faceLocalToGlobal = m_toFacesRelation.relatedObjectLocalToGlobal();

  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::nodeListString() ) );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toNodesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToNodes,
                                               packList,
                                               localToGlobal,
                                               nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toEdgesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToEdges,
                                               packList,
                                               localToGlobal,
                                               edgeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_toFacesRelation.base().toViewConst(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::surfaceElementsToCellRegionsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               this->m_surfaceElementsToCells,
                                               packList,
                                               m_surfaceElementsToCells.getElementRegionManager() );

  return packedSize;
}



localIndex FaceElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
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


  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOSX_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOSX_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation.base(),
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal() );

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

void FaceElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );

}

void FaceElementSubRegion::inheritGhostRankFromParentFace( FaceManager const & faceManager,
                                                           std::set< localIndex > const & indices )
{
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  for( localIndex const & index : indices )
  {
    m_ghostRank[index] = faceGhostRank[ m_toFacesRelation[index][0] ];
  }
}

} /* namespace geosx */
