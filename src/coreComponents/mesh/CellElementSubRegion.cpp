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

#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellElementSubRegion::CellElementSubRegion( string const & name, Group * const parent ):
  CellBlock( name, parent )
{
  registerWrapper( viewKeyStruct::constitutiveGroupingString, &m_constitutiveGrouping )->
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::constitutivePointVolumeFraction, &m_constitutivePointVolumeFraction );

  registerWrapper( viewKeyStruct::dNdXString, &m_dNdX )->setSizedFromParent( 1 )->reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString, &m_detJ )->setSizedFromParent( 1 )->reference();

  registerWrapper( viewKeyStruct::toEmbSurfString, &m_toEmbeddedSurfaces )->setSizedFromParent( 1 );
}

CellElementSubRegion::~CellElementSubRegion()
{
  // TODO Auto-generated destructor stub
}

void CellElementSubRegion::copyFromCellBlock( CellBlock * source )
{
  this->setElementType( source->getElementTypeString());
  this->setNumNodesPerElement( source->numNodesPerElement() );
  this->setNumFacesPerElement( source->numFacesPerElement() );
  this->resize( source->size());
  this->nodeList() = source->nodeList();

  arrayView1d< globalIndex const > const sourceLocalToGlobal = source->localToGlobalMap();
  this->m_localToGlobalMap.resize( sourceLocalToGlobal.size() );
  for( localIndex i = 0; i < localToGlobalMap().size(); ++i )
  {
    this->m_localToGlobalMap[ i ] = sourceLocalToGlobal[ i ];
  }

  this->constructGlobalToLocalMap();
  source->forExternalProperties( [&]( dataRepository::WrapperBase * const wrapper )
  {
    std::type_index typeIndex = std::type_index( wrapper->getTypeId());
    rtTypes::applyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                    true,
                                    [&]( auto type, auto GEOSX_UNUSED_PARAM( baseType ) )
    {
      using fieldType = decltype(type);
      dataRepository::Wrapper< fieldType > & field = dataRepository::Wrapper< fieldType >::cast( *wrapper );
      const fieldType & fieldref = field.reference();
      this->registerWrapper( wrapper->getName(), &const_cast< fieldType & >( fieldref ) ); //TODO remove const_cast
//      auto const & origFieldRef = field.reference();
//      fieldType & fieldRef = this->registerWrapper<fieldType>( wrapper->getName() )->reference();
//      fieldRef.resize( origFieldRef.size() );
    } );
  } );
}

void CellElementSubRegion::constructSubRegionFromFaceSet( FaceManager const * const faceManager,
                                                          string const & setName )
{
  SortedArrayView< localIndex const > const & targetSet = faceManager->sets().getReference< SortedArray< localIndex > >( setName );
  m_toFacesRelation.resize( 0, 2 );
  this->resize( targetSet.size() );
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
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString ));
//  exclusionList.insert(this->getWrapperIndex(this->viewKeys.edgeListString));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString ));
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
  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();


  localIndex packedSize = bufferOps::Pack< DOPACK >( buffer,
                                                     nodeList().base().toViewConst(),
                                                     m_unmappedGlobalIndicesInNodelist,
                                                     packList,
                                                     localToGlobal,
                                                     nodeLocalToGlobal );

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
                                     faceList().base(),
                                     packList,
                                     m_unmappedGlobalIndicesInFacelist,
                                     this->globalToLocalMap(),
                                     faceList().relatedObjectGlobalToLocal() );

  return unPackedSize;
}

void CellElementSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( faceList(),
                                    m_unmappedGlobalIndicesInFacelist,
                                    clearIfUnmapped );
}


} /* namespace geosx */
