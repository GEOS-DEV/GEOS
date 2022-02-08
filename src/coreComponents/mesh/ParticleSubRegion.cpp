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


#include "ParticleSubRegion.hpp"

#include "common/TypeDispatch.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

ParticleSubRegion::ParticleSubRegion( string const & name, Group * const parent ):
  ParticleBlock( name, parent )
{
  registerWrapper( viewKeyStruct::constitutiveGroupingString(), &m_constitutiveGrouping ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();
}

ParticleSubRegion::~ParticleSubRegion()
{
  // TODO Auto-generated destructor stub
}

void ParticleSubRegion::copyFromParticleBlock( ParticleBlock & source )
{
  this->setParticleType( source.getParticleType());
  this->resize( source.size());

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

//void ParticleSubRegion::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
//{
//  ObjectManagerBase::viewPackingExclusionList( exclusionList );
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString()  ));
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString()  ));
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString()  ));
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::toEmbSurfString() ));
//}


//localIndex ParticleSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
//{
//  buffer_unit_type * junk = nullptr;
//  return packUpDownMapsPrivate< false >( junk, packList );
//}


//localIndex ParticleSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
//                                                 arrayView1d< localIndex const > const & packList ) const
//{
//  return packUpDownMapsPrivate< true >( buffer, packList );
//}

//template< bool DOPACK >
//localIndex ParticleSubRegion::packUpDownMapsPrivate( buffer_unit_type * & buffer,
//                                                        arrayView1d< localIndex const > const & packList ) const
//{
//
//  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
//  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
//  arrayView1d< globalIndex const > edgeLocalToGlobal = edgeList().relatedObjectLocalToGlobal();
//  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();
//
//
//  localIndex packedSize = bufferOps::Pack< DOPACK >( buffer,
//                                                     nodeList().base().toViewConst(),
//                                                     m_unmappedGlobalIndicesInNodelist,
//                                                     packList,
//                                                     localToGlobal,
//                                                     nodeLocalToGlobal );
//
//  packedSize += bufferOps::Pack< DOPACK >( buffer,
//                                           edgeList().base().toViewConst(),
//                                           m_unmappedGlobalIndicesInEdgelist,
//                                           packList,
//                                           localToGlobal,
//                                           edgeLocalToGlobal );
//
//
//  packedSize += bufferOps::Pack< DOPACK >( buffer,
//                                           faceList().base().toViewConst(),
//                                           m_unmappedGlobalIndicesInFacelist,
//                                           packList,
//                                           localToGlobal,
//                                           faceLocalToGlobal );
//
//  return packedSize;
//}

//localIndex ParticleSubRegion::unpackUpDownMaps( buffer_unit_type const * & buffer,
//                                                   localIndex_array & packList,
//                                                   bool const GEOSX_UNUSED_PARAM( overwriteUpMaps ),
//                                                   bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
//{
//  localIndex unPackedSize = 0;
//  unPackedSize += bufferOps::Unpack( buffer,
//                                     nodeList().base().toView(),
//                                     packList,
//                                     m_unmappedGlobalIndicesInNodelist,
//                                     this->globalToLocalMap(),
//                                     nodeList().relatedObjectGlobalToLocal() );
//
//  unPackedSize += bufferOps::Unpack( buffer,
//                                     edgeList().base(),
//                                     packList,
//                                     m_unmappedGlobalIndicesInEdgelist,
//                                     this->globalToLocalMap(),
//                                     edgeList().relatedObjectGlobalToLocal() );
//
//  unPackedSize += bufferOps::Unpack( buffer,
//                                     faceList().base(),
//                                     packList,
//                                     m_unmappedGlobalIndicesInFacelist,
//                                     this->globalToLocalMap(),
//                                     faceList().relatedObjectGlobalToLocal() );
//
//  return unPackedSize;
//}


//void ParticleSubRegion::fixUpDownMaps( bool const clearIfUnmapped )
//{
//  ObjectManagerBase::fixUpDownMaps( nodeList(),
//                                    m_unmappedGlobalIndicesInNodelist,
//                                    clearIfUnmapped );
//
//  ObjectManagerBase::fixUpDownMaps( edgeList(),
//                                    m_unmappedGlobalIndicesInEdgelist,
//                                    clearIfUnmapped );
//
//  ObjectManagerBase::fixUpDownMaps( faceList(),
//                                    m_unmappedGlobalIndicesInFacelist,
//                                    clearIfUnmapped );
//}


} /* namespace geosx */
