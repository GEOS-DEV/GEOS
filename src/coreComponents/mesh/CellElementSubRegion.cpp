/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#include "CellElementSubRegion.hpp"

#include "mesh/MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellElementSubRegion::CellElementSubRegion( string const & name, Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_toFacesRelation(),
  m_externalPropertyNames()
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
}

CellElementSubRegion::~CellElementSubRegion()
{
  // Left blank
}

void CellElementSubRegion::setElementType( string const & elementType )
{
  m_elementTypeString = elementType;

  if( !m_elementTypeString.compare( 0, 4, "C3D8" ))
  {
    // Hexahedron
    m_numNodesPerElement = 8;
    m_numEdgesPerElement = 12;
    m_numFacesPerElement = 6;
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D4" ))
  {
    // Tetrahedron
    m_numNodesPerElement = 4;
    m_numEdgesPerElement = 6;
    m_numFacesPerElement = 4;

  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D6" ))
  {
    // Triangular prism
    m_numNodesPerElement = 6;
    m_numEdgesPerElement = 9;
    m_numFacesPerElement = 5;
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D5" ))
  {
    // Pyramid
    m_numNodesPerElement = 5;
    m_numEdgesPerElement = 8;
    m_numFacesPerElement = 5;
  }
  else
  {
    GEOSX_ERROR( "Error.  Don't know what kind of element this is." );
  }

  // If `setElementType` is called after the resize, the first dimension would be removed.
  // We do not want that so we try to keep it.
  m_toNodesRelation.resize( this->size(), m_numNodesPerElement );
  m_toEdgesRelation.resize( this->size(), m_numEdgesPerElement );
  m_toFacesRelation.resize( this->size(), m_numFacesPerElement );
}

void CellElementSubRegion::copyFromCellBlock( CellBlockABC & source )
{
  this->setElementType( source.getElementTypeString() );
  m_numNodesPerElement = source.numNodesPerElement();
  m_numFacesPerElement = source.numFacesPerElement();
  this->resize( source.size() );
  this->nodeList() = source.getElemToNode();
  this->edgeList() = source.getElemToEdges();
  this->faceList() = source.getElemToFaces();

  arrayView1d< globalIndex const > const sourceLocalToGlobal = source.localToGlobalMap();
  this->m_localToGlobalMap.resize( sourceLocalToGlobal.size() );
  for( localIndex i = 0; i < localToGlobalMap().size(); ++i )
  {
    this->m_localToGlobalMap[ i ] = sourceLocalToGlobal[ i ];
  }

  this->constructGlobalToLocalMap();
  source.forExternalProperties( [&]( dataRepository::WrapperBase & wrapper )
  {
    std::type_index typeIndex = std::type_index( wrapper.getTypeId());
    rtTypes::applyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                    true,
                                    [&]( auto type, auto GEOSX_UNUSED_PARAM( baseType ) )
    {
      using fieldType = decltype(type);
      dataRepository::Wrapper< fieldType > & field = dynamicCast< dataRepository::Wrapper< fieldType > & >( wrapper );
      const fieldType & fieldref = field.reference();
      this->registerWrapper( wrapper.getName(), &const_cast< fieldType & >( fieldref ) ); //TODO remove const_cast
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

void CellElementSubRegion::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::viewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString()  ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString()  ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString()  ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::toEmbSurfString() ));
}


localIndex CellElementSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}


localIndex CellElementSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                                 arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex CellElementSubRegion::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                        arrayView1d< localIndex const > const & packList ) const
{

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > edgeLocalToGlobal = edgeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();


  localIndex packedSize = bufferOps::Pack< DOPACK >( buffer,
                                                     nodeList().base().toViewConst(),
                                                     m_unmappedGlobalIndicesInNodelist,
                                                     packList,
                                                     localToGlobal,
                                                     nodeLocalToGlobal );

  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           edgeList().base().toViewConst(),
                                           m_unmappedGlobalIndicesInEdgelist,
                                           packList,
                                           localToGlobal,
                                           edgeLocalToGlobal );


  packedSize += bufferOps::Pack< DOPACK >( buffer,
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
  return packFracturedElementsPrivate< false >( junk, packList, embeddedSurfacesLocalToGlobal );
}


localIndex CellElementSubRegion::packFracturedElements( buffer_unit_type * & buffer,
                                                        arrayView1d< localIndex const > const & packList,
                                                        arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  return packFracturedElementsPrivate< true >( buffer, packList, embeddedSurfacesLocalToGlobal );
}

template< bool DOPACK >
localIndex CellElementSubRegion::packFracturedElementsPrivate( buffer_unit_type * & buffer,
                                                               arrayView1d< localIndex const > const & packList,
                                                               arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const
{
  localIndex packedSize = 0;

  // only here to use that packing function
  map< localIndex, array1d< globalIndex > > unmappedGlobalIndices;

  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::toEmbSurfString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           embeddedSurfacesList().base().toViewConst(),
                                           unmappedGlobalIndices,
                                           packList,
                                           localToGlobal,
                                           embeddedSurfacesLocalToGlobal );

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::fracturedCellsString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
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

localIndex CellElementSubRegion::getNumFaceNodes( localIndex const GEOSX_UNUSED_PARAM( elementIndex ),
                                       localIndex const localFaceIndex ) const
{
  if( !m_elementTypeString.compare( 0, 4, "C3D8" ))
    return 4;
  if( !m_elementTypeString.compare( 0, 4, "C3D6" ))
  {
    if( localFaceIndex == 0 )
      return 4;
    if( localFaceIndex == 1 )
      return 4;
    if( localFaceIndex == 2 )
      return 3;
    if( localFaceIndex == 3 )
      return 3;
    if( localFaceIndex == 4 )
      return 4;
  }
  if( !m_elementTypeString.compare( 0, 4, "C3D4" ))
    return 3;
  if( !m_elementTypeString.compare( 0, 4, "C3D5" ))
  {
    if( localFaceIndex == 0 )
      return 4;
    if( localFaceIndex == 1 )
      return 3;
    if( localFaceIndex == 2 )
      return 3;
    if( localFaceIndex == 3 )
      return 3;
    if( localFaceIndex == 4 )
      return 3;
  }

  GEOSX_ERROR( "Error. Don't know what kind of element this is and cannot build faces." );
  return -1;
}

// TODO this should not be there... Use the standard map instead.
localIndex CellElementSubRegion::getFaceNodes( localIndex const elementIndex,
                                    localIndex const localFaceIndex,
                                    localIndex * const nodeIndicies ) const
{
  if( !m_elementTypeString.compare( 0, 4, "C3D8" ))
  {
    if( localFaceIndex == 0 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
    }
    else if( localFaceIndex == 1 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][1];
    }
    else if( localFaceIndex == 2 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][4];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][6];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][2];
    }
    else if( localFaceIndex == 3 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][7];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][5];
    }
    else if( localFaceIndex == 4 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][6];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][7];
    }
    else if( localFaceIndex == 5 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][4];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][7];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][6];
    }
    return 4;
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D6" ))
  {
    if( localFaceIndex == 0 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
      return 4;
    }
    else if( localFaceIndex == 1 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][1];
      return 4;
    }
    else if( localFaceIndex == 2 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if( localFaceIndex == 3 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      return 3;
    }
    else if( localFaceIndex == 4 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][5];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][4];
      return 4;
    }
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D4" ))
  {
    if( localFaceIndex == 0 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][1];
    }
    else if( localFaceIndex == 1 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
    }
    else if( localFaceIndex == 2 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
    }
    else if( localFaceIndex == 3 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][3];
    }

    return 3;
  }
  else if( !m_elementTypeString.compare( 0, 4, "C3D5" ))
  {
    if( localFaceIndex == 0 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[3] = m_toNodesRelation[elementIndex][3];
      return 4;
    }
    else if( localFaceIndex == 1 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if( localFaceIndex == 2 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][1];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if( localFaceIndex == 3 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][2];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
    else if( localFaceIndex == 4 )
    {
      nodeIndicies[0] = m_toNodesRelation[elementIndex][3];
      nodeIndicies[1] = m_toNodesRelation[elementIndex][0];
      nodeIndicies[2] = m_toNodesRelation[elementIndex][4];
      return 3;
    }
  }

  GEOSX_ERROR( "Error. Don't know what kind of element this is and cannot build faces." );
  return -1;
}

void CellElementSubRegion::getFaceNodes( localIndex const elementIndex,
                              localIndex const localFaceIndex,
                              localIndex_array & nodeIndicies ) const
{
  nodeIndicies.resize( getNumFaceNodes( elementIndex, localFaceIndex ) );
  localIndex const numNodes = getFaceNodes( elementIndex, localFaceIndex, nodeIndicies.data() );
  GEOSX_DEBUG_VAR( numNodes );
  GEOSX_ASSERT_EQ( numNodes, nodeIndicies.size() );
}

void CellElementSubRegion::calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                     FaceManager const & GEOSX_UNUSED_PARAM( faceManager ) )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  forAll< serialPolicy >( this->size(), [=] ( localIndex const k )
  {
    calculateCellVolumesKernel( k, X );
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
