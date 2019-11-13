/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
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
  registerWrapper( viewKeyStruct::constitutiveGroupingString, &m_constitutiveGrouping, 0)->
    setSizedFromParent(0);

  registerWrapper( viewKeyStruct::constitutiveMapString,
                       &m_constitutiveMapView, 0);

  registerWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0);


  registerWrapper( viewKeyStruct::constitutivePointVolumeFraction,
                       &m_constitutivePointVolumeFraction, 0);

  registerWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0)->setSizedFromParent(1);

}

CellElementSubRegion::~CellElementSubRegion()
{
  // TODO Auto-generated destructor stub
}

void CellElementSubRegion::CopyFromCellBlock( CellBlock const * source )
{
  this->SetElementType(source->GetElementTypeString());
  this->numNodesPerElement() = source->numNodesPerElement();
  this->numFacesPerElement() = source->numFacesPerElement();
  this->resize(source->size());
  this->nodeList() = source->nodeList();
  this->m_localToGlobalMap = source->m_localToGlobalMap;
  this->ConstructGlobalToLocalMap();
  source->forExternalProperties([&]( const dataRepository::WrapperBase * wrapper )->void
  {
    std::type_index typeIndex = std::type_index( wrapper->get_typeid());
    rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                    true,
                                    [&]( auto type, auto GEOSX_UNUSED_ARG( baseType ) ) -> void
    {
      using fieldType = decltype(type);
      const dataRepository::Wrapper<fieldType> & field = dataRepository::Wrapper< fieldType >::cast( *wrapper );
      const fieldType & fieldref = field.reference();
      this->registerWrapper( wrapper->getName(), &const_cast< fieldType & >( fieldref ), 0 ); //TODO remove const_cast
    });
  });
}

void CellElementSubRegion::ConstructSubRegionFromFaceSet( FaceManager const * const faceManager,
                                                        string const & setName )
{
  set<localIndex> const & targetSet = faceManager->sets()->getReference<set<localIndex> >(setName);
  m_toFacesRelation.resize(0,2);
  this->resize( targetSet.size() );


}

void CellElementSubRegion::MaterialPassThru( string const & GEOSX_UNUSED_ARG( matName ),
                                             string const & GEOSX_UNUSED_ARG( setName ),
                                             set<localIndex> & GEOSX_UNUSED_ARG( materialSet ),
                                             Group * GEOSX_UNUSED_ARG( material ) )
{}





void CellElementSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
//  exclusionList.insert(this->getWrapperIndex(this->viewKeys.edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
}


localIndex CellElementSubRegion::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}


localIndex CellElementSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                                               arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
localIndex CellElementSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                      arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         nodeList().Base().toViewConst(),
                                         m_unmappedGlobalIndicesInNodelist,
                                         packList,
                                         this->m_localToGlobalMap,
                                         nodeList().RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         faceList().Base().toViewConst(),
                                         m_unmappedGlobalIndicesInFacelist,
                                         packList,
                                         this->m_localToGlobalMap,
                                         faceList().RelatedObjectLocalToGlobal() );

  return packedSize;
}


localIndex CellElementSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                                 localIndex_array & packList,
                                                 bool const GEOSX_UNUSED_ARG( overwriteUpMaps ),
                                                 bool const GEOSX_UNUSED_ARG( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;
  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().Base().toView(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->m_globalToLocalMap,
                                     nodeList().RelatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     faceList().Base(),
                                     packList,
                                     m_unmappedGlobalIndicesInFacelist,
                                     this->m_globalToLocalMap,
                                     faceList().RelatedObjectGlobalToLocal() );

  return unPackedSize;
}

void CellElementSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( faceList(),
                                    m_unmappedGlobalIndicesInFacelist,
                                    clearIfUnmapped );
}


} /* namespace geosx */
