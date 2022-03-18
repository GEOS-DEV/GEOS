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
#include "mesh/MeshLevel.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

ParticleSubRegion::ParticleSubRegion( string const & name, Group * const parent ):
  ParticleSubRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::constitutiveGroupingString(), &m_constitutiveGrouping ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::constitutivePointVolumeFractionString(), &m_constitutivePointVolumeFraction );

  registerWrapper( viewKeyStruct::dNdXString(), &m_dNdX ).setSizedFromParent( 1 ).reference().resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString(), &m_detJ ).setSizedFromParent( 1 ).reference();
}

ParticleSubRegion::~ParticleSubRegion()
{
  // Left blank
}

void ParticleSubRegion::copyFromParticleBlock( ParticleBlockABC & particleBlock )
{
  // Defines the (unique) particle type of this cell particle region,
  // and its associated number of nodes, edges, faces.
  m_particleType = particleBlock.getParticleType();
  m_particleCenter = particleBlock.getParticleCenter();
  m_particleVelocity = particleBlock.getParticleVelocity();
  m_particleVolume = particleBlock.getParticleVolume();
  m_particleVolume0 = particleBlock.getParticleVolume();
  m_particleDeformationGradient.resize(particleBlock.size(),3,3);
  m_particleRVectors = particleBlock.getParticleRVectors();
  m_particleRVectors0 = particleBlock.getParticleRVectors();
  //m_particleMass.resize(particleBlock.size()); // handled by ParticleSubRegionBase constructor?
  m_hasRVectors = particleBlock.hasRVectors();


  // We call the `resize` member function of the cell to (nodes, edges, faces) relations,
  // before calling the `ParticleSubRegion::resize` in order to keep the first dimension.
  // Be careful when refactoring.
  this->resize( particleBlock.numParticles() );

  this->m_localToGlobalMap = particleBlock.localToGlobalMap();

  this->constructGlobalToLocalMap();
  particleBlock.forExternalProperties( [&]( WrapperBase & wrapper )
  {
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
    {
      using ArrayType = decltype( array );
      Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
      this->registerWrapper( wrapper.getName(), &wrapperT.reference() );
    } );
  } );
}

void ParticleSubRegion::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::viewPackingExclusionList( exclusionList );
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString()  ));
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString()  ));
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString()  ));
//  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::toEmbSurfString() ));
}


localIndex ParticleSubRegion::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}


localIndex ParticleSubRegion::packUpDownMaps( buffer_unit_type * & buffer,
                                              arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex ParticleSubRegion::packUpDownMapsPrivate( buffer_unit_type * & GEOSX_UNUSED_PARAM(buffer),
                                                     arrayView1d< localIndex const > const & GEOSX_UNUSED_PARAM(packList) ) const
{

//  arrayView1d< globalIndex const > const localToGlobal = this->localToGlobalMap();
//  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
//  arrayView1d< globalIndex const > edgeLocalToGlobal = edgeList().relatedObjectLocalToGlobal();
//  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();


  localIndex packedSize = 0;

//  packedSize += bufferOps::Pack< DOPACK >( buffer,
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

  return packedSize;
}

localIndex ParticleSubRegion::unpackUpDownMaps( buffer_unit_type const * & GEOSX_UNUSED_PARAM(buffer),
                                                   localIndex_array & GEOSX_UNUSED_PARAM(packList),
                                                   bool const GEOSX_UNUSED_PARAM( overwriteUpMaps ),
                                                   bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;
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

  return unPackedSize;
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleSubRegion, string const &, Group * const )

} /* namespace geosx */
