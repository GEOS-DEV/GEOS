/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

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
  CellBlock( name, parent ),
  m_constitutiveModels(groupKeyStruct::constitutiveModelsString,this)
{
  RegisterViewWrapper( viewKeyStruct::constitutiveGroupingString, &m_constitutiveGrouping, 0)->
      setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::constitutiveMapString,
                       &m_constitutiveMapView, 0);

  RegisterViewWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0);

//  RegisterViewWrapper( viewKeyStruct::constitutiveRelationIndexString,
//                       &m_constitutiveRelationIndex, 0);
//
//  RegisterViewWrapper( viewKeyStruct::constitutivePointIndexString,
//                       &m_constitutivePointIndex, 0);

  RegisterViewWrapper( viewKeyStruct::constitutivePointVolumeFraction,
                       &m_constitutivePointVolumeFraction, 0);

  RegisterViewWrapper( viewKeyStruct::dNdXString, &m_dNdX, 0)->setSizedFromParent(1);


  RegisterGroup( groupKeyStruct::constitutiveModelsString, &m_constitutiveModels, 0 );
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
  this->nodeList() = source->nodeList();
  this->m_localToGlobalMap = source->m_localToGlobalMap;
  this->ConstructGlobalToLocalMap();
}

void CellBlockSubRegion::MaterialPassThru( string const & matName,
                                           string const & setName,
                                           set<localIndex> & materialSet,
                                           ManagedGroup * material )
{}





void CellBlockSubRegion::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::nodeListString));
//  exclusionList.insert(this->getWrapperIndex(this->viewKeys.edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
}


localIndex CellBlockSubRegion::PackUpDownMapsSize( localIndex_array const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}


localIndex CellBlockSubRegion::PackUpDownMaps( buffer_unit_type * & buffer,
                               localIndex_array const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
localIndex CellBlockSubRegion::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               localIndex_array const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                             nodeList(),
                                             packList,
                                             this->m_localToGlobalMap,
                                             nodeList().RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                             faceList(),
                                             packList,
                                             this->m_localToGlobalMap,
                                             faceList().RelatedObjectLocalToGlobal() );

  return packedSize;
}


localIndex CellBlockSubRegion::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                 localIndex_array const & packList )
{
  localIndex unPackedSize = 0;

  unPackedSize += bufferOps::Unpack( buffer,
                                         nodeList(),
                                         packList,
                                         this->m_globalToLocalMap,
                                         nodeList().RelatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer,
                                         faceList(),
                                         packList,
                                         this->m_globalToLocalMap,
                                         faceList().RelatedObjectGlobalToLocal() );

  return unPackedSize;
}


} /* namespace geosx */
