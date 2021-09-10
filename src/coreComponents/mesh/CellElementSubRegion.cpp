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
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellElementSubRegion::CellElementSubRegion( string const & name, Group * const parent ):
  CellBlock( name, parent )
{
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
  // TODO Auto-generated destructor stub
}

void CellElementSubRegion::copyFromCellBlock( CellBlock & source )
{
  this->setElementType( source.getElementType());
  this->setNumNodesPerElement( source.numNodesPerElement() );
  this->setNumFacesPerElement( source.numFacesPerElement() );
  this->resize( source.size());
  this->nodeList() = source.nodeList();

  arrayView1d< globalIndex const > const sourceLocalToGlobal = source.localToGlobalMap();
  this->m_localToGlobalMap.resize( sourceLocalToGlobal.size() );
  for( localIndex i = 0; i < localToGlobalMap().size(); ++i )
  {
    this->m_localToGlobalMap[ i ] = sourceLocalToGlobal[ i ];
  }

  this->constructGlobalToLocalMap();
  source.forExternalProperties( [&]( WrapperBase & wrapper )
  {
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
    {
      using ArrayType = decltype( array );
      Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
      this->registerWrapper( wrapper.getName(), &wrapperT.reference() );
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


} /* namespace geosx */
