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

#include "FaceElementSubRegion.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
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

  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString(), &m_edgesTo2dFaces ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of edge local indices to the fracture connector local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString(), &m_2dFaceToEdge ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices to edge local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString(), &m_2dFaceTo2dElems ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices face element local indices" ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::elem2dToCollocatedNodesBucketsString(), &m_2dElemToCollocatedNodesBuckets ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map eventually containing all the collocated nodes." ).
    setSizedFromParent( 1 );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  registerWrapper( viewKeyStruct::separationCoeffString(), &m_separationCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
    setDescription( "Scalar indicator of level of separation for a fracturing face." );
#endif

  excludeWrappersFromPacking( { viewKeyStruct::faceListString() } );

  m_2dElemToElems.resize( 0, 2 );

  m_numNodesPerElement = 8;
}

void FaceElementSubRegion::copyFromCellBlock( FaceBlockABC const & faceBlock )
{
  localIndex const num2dElements = faceBlock.num2dElements();
  resize( faceBlock.num2dElements() );

  m_toNodesRelation.base() = faceBlock.get2dElemToNodes();
  m_toEdgesRelation.base() = faceBlock.get2dElemToEdges();

  // `FaceBlockABC` is designed to be heterogeneous.
  // `FaceElementSubRegion` inherits from `ElementSubRegionBase` which is meant to be homogeneous.
  // But `FaceElementSubRegion` is sometimes used as an heterogeneous sub region,
  // which emphasizes the need of a refactoring.
  // In the meantime, we try to fill the face block into the sub region and hope for the best...
  {
    auto const deduce3dElemType = [&]( integer maxSize ) -> ElementType
    {
      switch( maxSize )
      {
        case 3:
          return ElementType::Wedge;
        case 4:
          return ElementType::Hexahedron;
        case 5:
          return ElementType::Prism5;
        case 6:
          return ElementType::Prism6;
        case 7:
          return ElementType::Prism7;
        case 8:
          return ElementType::Prism8;
        case 9:
          return ElementType::Prism9;
        case 10:
          return ElementType::Prism10;
        case 11:
          return ElementType::Prism11;
        case 0:
          // In the case the fracture is empty (on this rank), then we default to hexahedron. Otherwise, there's something wrong
          GEOS_ERROR_IF_NE_MSG( num2dElements, 0, "Could not determine the element type of the fracture \"" << getName() << "\"." );
          return ElementType::Hexahedron;
        default:
          GEOS_ERROR( "Unsupported type of elements during the face element sub region creation." );
          return {};
      }
    };

    m_2dElemToCollocatedNodesBuckets = faceBlock.get2dElemsToCollocatedNodesBuckets();
    // Checking if all the 2d elements are homogeneous.
    // We rely on the number of nodes for each element to find out.
    std::vector< integer > numNodesPerElement( num2dElements );
    for( int i = 0; i < num2dElements; ++i )
    {
      numNodesPerElement[i] = m_2dElemToCollocatedNodesBuckets[i].size();
    }
    std::set< integer > const sizes( numNodesPerElement.cbegin(), numNodesPerElement.cend() );

    if( sizes.size() > 1 )
    {
      // If we have found that the input face block contains 2d elements of different types,
      // we inform the used that the situation may be at risk.
      // (We're storing the face block in a homogeneous container while it's actually heterogeneous).
      GEOS_WARNING( "Heterogeneous face element sub region found and stored as homogeneous. Use at your own risk." );
    }

    auto const it = std::max_element( sizes.cbegin(), sizes.cend() );
    integer const maxSize = it != sizes.cend() ? *it : 0;
    m_elementType = deduce3dElemType( maxSize );
    m_numNodesPerElement = maxSize;
  }

  // The `m_2dElemToElems` mappings involves element, sub regions and regions indices.
  // We store the element indices that are correct.
  // But we only have access to the cell block indices, not the sub regions indices.
  // Temporarily, and also because they share the same dimensions,
  // we store the cell block mapping at the sub region mapping location.
  // It will later be transformed into a sub regions mapping.
  // Last, we fill the regions mapping with dummy -1 values that should all be replaced eventually.
  auto const elem2dToElems = faceBlock.get2dElemToElems();
  m_2dElemToElems.resize( num2dElements, 2 );
  for( int i = 0; i < num2dElements; ++i )
  {
    for( localIndex const & j: elem2dToElems.toCellIndex[i] )
    {
      m_2dElemToElems.m_toElementIndex.emplaceBack( i, j );
    }
    for( localIndex const & j: elem2dToElems.toBlockIndex[i] )
    {
      m_2dElemToElems.m_toElementSubRegion.emplaceBack( i, j );
    }
  }

  m_toFacesRelation.base() = faceBlock.get2dElemToFaces();

  m_2dFaceToEdge = faceBlock.get2dFaceToEdge();
  m_2dFaceTo2dElems = faceBlock.get2dFaceTo2dElems();

  m_localToGlobalMap = faceBlock.localToGlobalMap();
  this->constructGlobalToLocalMap();

  for( int i = 0; i < faceBlock.num2dFaces(); ++i )
  {
    m_recalculateConnectionsFor2dFaces.insert( i );
  }

  for( localIndex i = 0; i < faceBlock.num2dElements(); ++i )
  {
    m_newFaceElements.insert( i );
  }

  // TODO We still need to be able to import fields on the FaceElementSubRegion.
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
}

void FaceElementSubRegion::calculateElementGeometricQuantities( NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                                                                FaceManager const & faceManager )
{
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();

  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateSingleElementGeometricQuantities( k, faceArea );
  } );
}

ElementType FaceElementSubRegion::getElementType( localIndex ei ) const
{
  // We try to get the element type information from the bucket.
  // If this information does not appear to be reliable,
  // let's fall back on the 2d elem to nodes.
  auto const sizeFromBuckets = m_2dElemToCollocatedNodesBuckets[ei].size();
  auto const size = sizeFromBuckets > 0 ? sizeFromBuckets : m_toNodesRelation[ei].size() / 2;
  switch( size )
  {
    case 3:
      return ElementType::Triangle;
    case 4:
      return ElementType::Quadrilateral;
    default:
      return ElementType::Polygon;
  }
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
                                               m_toFacesRelation.toViewConst(),
                                               m_unmappedGlobalIndicesInToFaces,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::surfaceElementsToCellRegionsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_2dElemToElems,
                                               packList,
                                               m_2dElemToElems.getElementRegionManager() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::elem2dToCollocatedNodesBucketsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_2dElemToCollocatedNodesBuckets.toViewConst(),
                                               packList,
                                               localToGlobal );

  return packedSize;
}


localIndex FaceElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const overwriteUpMaps,
                                                   bool const GEOS_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOS_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOS_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOS_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOS_ERROR_IF_NE( elementListString, viewKeyStruct::surfaceElementsToCellRegionsString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToElems,
                                     packList.toViewConst(),
                                     m_2dElemToElems.getElementRegionManager(),
                                     overwriteUpMaps );

  string elem2dToCollocatedNodesBucketsString;
  unPackedSize += bufferOps::Unpack( buffer, elem2dToCollocatedNodesBucketsString );
  GEOS_ERROR_IF_NE( elem2dToCollocatedNodesBucketsString, viewKeyStruct::elem2dToCollocatedNodesBucketsString() );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_2dElemToCollocatedNodesBuckets,
                                     packList,
                                     this->globalToLocalMap() );

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
  for( localIndex const & index: indices )
  {
    m_ghostRank[index] = faceGhostRank[m_toFacesRelation[index][0]];
  }
}

std::set< std::set< globalIndex > > FaceElementSubRegion::getCollocatedNodes() const
{
  std::set< std::set< globalIndex > > result;

  for( localIndex e2d = 0; e2d < m_2dElemToCollocatedNodesBuckets.size(); ++e2d )
  {
    for( integer ni = 0; ni < m_2dElemToCollocatedNodesBuckets[e2d].size(); ++ni )
    {
      array1d< globalIndex > const & bucket = m_2dElemToCollocatedNodesBuckets( e2d, ni );
      std::set< globalIndex > const s( bucket.begin(), bucket.end() );
      result.insert( s );
    }
  }

  return result;
}

} /* namespace geos */
