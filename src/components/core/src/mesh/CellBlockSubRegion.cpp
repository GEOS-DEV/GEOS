/*
 * CellBlockSubRegion.cpp
 *
 *  Created on: May 11, 2017
 *      Author: rrsettgast
 */

#include "CellBlockSubRegion.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

CellBlockSubRegion::CellBlockSubRegion( string const & name, ManagedGroup * const parent ):
  CellBlock( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::constitutiveGroupingString, &m_constitutiveGrouping, 0)->
      setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::constitutiveMapString, &m_constitutiveMapView, 0);

  RegisterViewWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0)->setSizedFromParent(1);



}

CellBlockSubRegion::~CellBlockSubRegion()
{
  // TODO Auto-generated destructor stub
}



void CellBlockSubRegion::FillDocumentationNode()
{
  CellBlock::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->getCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "an element region" );

  CellBlock::FillDocumentationNode();

//  docNode->AllocateChildNode( viewKeys.numNodesPerElement.Key(),
//                              viewKeys.numNodesPerElement.Key(),
//                              -1,
//                              "integer",
//                              "integer",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              0,
//                              1,
//                              0 );

//  docNode->AllocateChildNode( keys::numNodesPerElement,
//                              keys::numNodesPerElement,
//                              -1,
//                              "integer",
//                              "integer",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              0,
//                              1,
//                              0 );



//  docNode->AllocateChildNode( keys::constitutiveMap,
//                              keys::constitutiveMap,
//                              -1,
//                              "mapPair_array",
//                              "mapPair_array",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              1,
//                              0,
//                              0 );

//  docNode->AllocateChildNode( keys::constitutiveMap,
//                              keys::constitutiveMap,
//                              -1,
//                              "mapPair_array",
//                              "mapPair_array",
//                              "Number of Nodes Per Element",
//                              "Number of Nodes Per Element",
//                              "1",
//                              "",
//                              1,
//                              0,
//                              0 );



}

void CellBlockSubRegion::ReadXML_PostProcess()
{
//  integer & numNodesPerElem = numNodesPerElement();
//  numNodesPerElem = 8;

}

void CellBlockSubRegion::InitializePreSubGroups( ManagedGroup * const )
{
//  auto const & elementRegion = static_cast<ElementRegion const&>(
// *(this->getParent()) );
//  auto const & numMethod = elementRegion.getNumericalMethod();

}

void CellBlockSubRegion::InitializePostSubGroups( ManagedGroup * const )
{
  ObjectManagerBase::InitializePostSubGroups(nullptr);
  this->numNodesPerElement() = 8;
  this->numFacesPerElement() = 6;
}

void CellBlockSubRegion::CopyFromCellBlock( CellBlock const * source )
{
  this->resize(source->size());
  this->m_toNodesRelation = source->m_toNodesRelation;
  this->m_localToGlobalMap = source->m_localToGlobalMap;
  this->ConstructGlobalToLocalMap();
}

void CellBlockSubRegion::MaterialPassThru( string const & matName,
                                           string const & setName,
                                           lSet & materialSet,
                                           ManagedGroup * material )
{}





void CellBlockSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
//  exclusionList.insert(this->getWrapperIndex(this->viewKeys.edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
}


int CellBlockSubRegion::PackUpDownMapsSize( localIndex_array const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}


int CellBlockSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                               localIndex_array const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
int CellBlockSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               localIndex_array const & packList ) const
{
  int packedSize = 0;

  packedSize += CommBufferOps::Pack<DOPACK>( buffer,
                                             m_toNodesRelation,
                                             packList,
                                             this->m_localToGlobalMap,
                                             m_toNodesRelation.RelatedObjectLocalToGlobal() );

  packedSize += CommBufferOps::Pack<DOPACK>( buffer,
                                             m_toFacesRelation,
                                             packList,
                                             this->m_localToGlobalMap,
                                             m_toFacesRelation.RelatedObjectLocalToGlobal() );

  return packedSize;
}


int CellBlockSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                 localIndex_array const & packList )
{
  int unPackedSize = 0;

  unPackedSize += CommBufferOps::Unpack( buffer,
                                         m_toNodesRelation,
                                         packList,
                                         this->m_globalToLocalMap,
                                         m_toNodesRelation.RelatedObjectGlobalToLocal() );

  unPackedSize += CommBufferOps::Unpack( buffer,
                                         m_toFacesRelation,
                                         packList,
                                         this->m_globalToLocalMap,
                                         m_toFacesRelation.RelatedObjectGlobalToLocal() );

  return unPackedSize;
}


} /* namespace geosx */
