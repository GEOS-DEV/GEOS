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


#include "CellElementSubRegion.hpp"

#include "common/TypeDispatch.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/generators/CellBlockUtilities.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellElementSubRegion::CellElementSubRegion( string const & name, Group * const parent ):
  ElementSubRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation );
  registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation );

  registerWrapper( viewKeyStruct::constitutiveGroupingString(), &m_constitutiveGrouping ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::constitutivePointVolumeFractionString(), &m_constitutivePointVolumeFraction );

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();

  registerWrapper( viewKeyStruct::toEmbSurfString(), &m_toEmbeddedSurfaces ).setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::fracturedCellsString(), &m_fracturedCells ).setSizedFromParent( 1 );

  excludeWrappersFromPacking( { viewKeyStruct::nodeListString(),
                                viewKeyStruct::edgeListString(),
                                viewKeyStruct::faceListString(),
                                viewKeyStruct::fracturedCellsString(),
                                viewKeyStruct::toEmbSurfString() } );
}

void CellElementSubRegion::copyFromCellBlock( CellBlockABC & cellBlock )
{
  // Defines the (unique) element type of this cell element region,
  // and its associated number of nodes, edges, faces.
  m_elementType = cellBlock.getElementType();
  m_numNodesPerElement = cellBlock.numNodesPerElement();
  m_numEdgesPerElement = cellBlock.numEdgesPerElement();
  m_numFacesPerElement = cellBlock.numFacesPerElement();

  // We call the `resize` member function of the cell to (nodes, edges, faces) relations,
  // before calling the `CellElementSubRegion::resize` in order to keep the first dimension.
  // Be careful when refactoring.
  m_toNodesRelation.resize( this->size(), m_numNodesPerElement );
  m_toEdgesRelation.resize( this->size(), m_numEdgesPerElement );
  m_toFacesRelation.resize( this->size(), m_numFacesPerElement );
  this->resize( cellBlock.numElements() );

  this->nodeList() = cellBlock.getElemToNodes();
  this->edgeList() = cellBlock.getElemToEdges();
  this->faceList() = cellBlock.getElemToFaces();

  this->m_localToGlobalMap = cellBlock.localToGlobalMap();

  this->constructGlobalToLocalMap();
  cellBlock.forExternalProperties( [&]( WrapperBase & wrapper )
  {
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
    {
      using ArrayType = decltype( array );
      Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
      this->registerWrapper( wrapper.getName(), std::make_unique< ArrayType >( wrapperT.reference() ) );
    } );
  } );
}

void CellElementSubRegion::addFracturedElement( localIndex const cellElemIndex,
                                                localIndex const embSurfIndex )
{
  // add the connection between the element and the embedded surface to the map
  m_toEmbeddedSurfaces.emplaceBack( cellElemIndex, embSurfIndex );
  // add the element to the fractured elements list
  m_fracturedCells.insert( cellElemIndex );
}


localIndex CellElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsImpl< false >( junk, packList );
}


localIndex CellElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsImpl< true >( buffer, packList );
}

template< bool DO_PACKING >
localIndex CellElementSubRegion::packUpDownMapsImpl( buffer_unit_type * & buffer,
                                                     arrayView1d< localIndex const > const & packList ) const
{

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > edgeLocalToGlobal = edgeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();


  localIndex packedSize = bufferOps::Pack< DO_PACKING >( buffer,
                                                         nodeList().base().toViewConst(),
                                                         m_unmappedGlobalIndicesInNodelist,
                                                         packList,
                                                         localToGlobal,
                                                         nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               edgeList().base().toViewConst(),
                                               m_unmappedGlobalIndicesInEdgelist,
                                               packList,
                                               localToGlobal,
                                               edgeLocalToGlobal );


  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               faceList().base().toViewConst(),
                                               m_unmappedGlobalIndicesInFacelist,
                                               packList,
                                               localToGlobal,
                                               faceLocalToGlobal );

  return packedSize;
}

localIndex CellElementSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                   localIndex_array & packList,
                                                   bool const GEOSX_UNUSED_PARAM( overwriteUpMaps ),
                                                   bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;
  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().base().toView(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->globalToLocalMap(),
                                     nodeList().relatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     edgeList().base(),
                                     packList,
                                     m_unmappedGlobalIndicesInEdgelist,
                                     this->globalToLocalMap(),
                                     edgeList().relatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     faceList().base(),
                                     packList,
                                     m_unmappedGlobalIndicesInFacelist,
                                     this->globalToLocalMap(),
                                     faceList().relatedObjectGlobalToLocal() );

  return unPackedSize;
}

localIndex CellElementSubRegion::packFracturedElementsSize( arrayView1d< localIndex const > const & packList,
                                                            arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  buffer_unit_type * junk = nullptr;
  return packFracturedElementsImpl< false >( junk, packList, embeddedSurfacesLocalToGlobal );
}


localIndex CellElementSubRegion::packFracturedElements( buffer_unit_type * & buffer,
                                                        arrayView1d< localIndex const > const & packList,
                                                        arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  return packFracturedElementsImpl< true >( buffer, packList, embeddedSurfacesLocalToGlobal );
}

template< bool DO_PACKING >
localIndex CellElementSubRegion::packFracturedElementsImpl( buffer_unit_type * & buffer,
                                                            arrayView1d< localIndex const > const & packList,
                                                            arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  localIndex packedSize = 0;

  // only here to use that packing function
  map< localIndex, array1d< globalIndex > > unmappedGlobalIndices;

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::toEmbSurfString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               embeddedSurfacesList().base().toViewConst(),
                                               unmappedGlobalIndices,
                                               packList,
                                               localToGlobal,
                                               embeddedSurfacesLocalToGlobal );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( viewKeyStruct::fracturedCellsString() ) );
  packedSize += bufferOps::Pack< DO_PACKING >( buffer,
                                               m_fracturedCells.toViewConst(),
                                               packList,
                                               localToGlobal );

  return packedSize;
}


localIndex CellElementSubRegion::unpackFracturedElements( buffer_unit_type const * & buffer,
                                                          localIndex_array & packList,
                                                          unordered_map< globalIndex, localIndex > const & embeddedSurfacesGlobalToLocal )
{
  localIndex unPackedSize = 0;

  string toEmbSurfString;
  unPackedSize += bufferOps::Unpack( buffer, toEmbSurfString );
  GEOSX_ERROR_IF_NE( toEmbSurfString, viewKeyStruct::toEmbSurfString() );

  // only here to use that packing function
  map< localIndex, array1d< globalIndex > > unmappedGlobalIndices;

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEmbeddedSurfaces,
                                     packList,
                                     unmappedGlobalIndices,
                                     this->globalToLocalMap(),
                                     embeddedSurfacesGlobalToLocal );

  string fracturedCellsString;
  unPackedSize += bufferOps::Unpack( buffer, fracturedCellsString );
  GEOSX_ERROR_IF_NE( fracturedCellsString, viewKeyStruct::fracturedCellsString() );

  SortedArray< globalIndex > junk;
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_fracturedCells,
                                     junk,
                                     this->globalToLocalMap(),
                                     false );

  return unPackedSize;
}



void CellElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( edgeList(),
                                    m_unmappedGlobalIndicesInEdgelist,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( faceList(),
                                    m_unmappedGlobalIndicesInFacelist,
                                    clearIfUnmapped );
}

localIndex CellElementSubRegion::getFaceNodes( localIndex const elementIndex,
                                               localIndex const localFaceIndex,
                                               Span< localIndex > const nodeIndices ) const
{
  return geosx::getFaceNodes( m_elementType, elementIndex, localFaceIndex, m_toNodesRelation, nodeIndices );
}

void CellElementSubRegion::
  calculateElementCenterAndVolume( localIndex const k,
                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
{
  LvArray::tensorOps::fill< 3 >( m_elementCenter[ k ], 0 );

  switch( m_elementType )
  {
    case ElementType::Hexahedron:
    {
      real64 Xlocal[8][3];
      for( localIndex a = 0; a < m_numNodesPerElement; ++a )
      {
        LvArray::tensorOps::copy< 3 >( Xlocal[ a ], X[ m_toNodesRelation( k, a ) ] );
        LvArray::tensorOps::add< 3 >( m_elementCenter[ k ], Xlocal[ a ] );
      }
      m_elementVolume[k] = computationalGeometry::hexVolume( Xlocal );
      break;
    }
    case ElementType::Tetrahedron:
    {
      real64 Xlocal[4][3];
      for( localIndex a = 0; a < m_numNodesPerElement; ++a )
      {
        LvArray::tensorOps::copy< 3 >( Xlocal[ a ], X[ m_toNodesRelation( k, a ) ] );
        LvArray::tensorOps::add< 3 >( m_elementCenter[ k ], Xlocal[ a ] );
      }
      m_elementVolume[k] = computationalGeometry::tetVolume( Xlocal );
      break;
    }
    case ElementType::Prism:
    {
      real64 Xlocal[6][3];
      for( localIndex a = 0; a < m_numNodesPerElement; ++a )
      {
        LvArray::tensorOps::copy< 3 >( Xlocal[ a ], X[ m_toNodesRelation( k, a ) ] );
        LvArray::tensorOps::add< 3 >( m_elementCenter[ k ], Xlocal[ a ] );
      }
      m_elementVolume[k] = computationalGeometry::wedgeVolume( Xlocal );
      break;
    }
    case ElementType::Pyramid:
    {
      real64 Xlocal[5][3];
      for( localIndex a = 0; a < m_numNodesPerElement; ++a )
      {
        LvArray::tensorOps::copy< 3 >( Xlocal[ a ], X[ m_toNodesRelation( k, a ) ] );
        LvArray::tensorOps::add< 3 >( m_elementCenter[ k ], Xlocal[ a ] );
      }
      m_elementVolume[k] = computationalGeometry::pyramidVolume( Xlocal );
      break;
    }
    default:
    {
      GEOSX_ERROR( GEOSX_FMT( "Volume calculation not supported for element type {} in subregion {}",
                              m_elementType, getName() ) );
    }
  }

  LvArray::tensorOps::scale< 3 >( m_elementCenter[ k ], 1.0 / m_numNodesPerElement );
}

void CellElementSubRegion::calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                                FaceManager const & GEOSX_UNUSED_PARAM( faceManager ) )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();

  forAll< parallelHostPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateElementCenterAndVolume( k, X );
  } );
}

void CellElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const & mesh )
{
  this->m_toNodesRelation.setRelatedObject( mesh.getNodeManager() );
  this->m_toEdgesRelation.setRelatedObject( mesh.getEdgeManager() );
  this->m_toFacesRelation.setRelatedObject( mesh.getFaceManager() );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementSubRegion, string const &, Group * const )

} /* namespace geosx */
