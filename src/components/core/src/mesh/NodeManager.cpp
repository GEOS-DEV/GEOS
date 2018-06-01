// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * NodeManagerT.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "NodeManager.hpp"
//#include "managers/DomainPartition.hpp"
//#include "ObjectManagers/FaceManagerT.h"
//#include "ObjectManagers/EdgeManagerT.h"
//#include "ObjectManagers/ElementManagerT.h"
//#include "Utilities/Utilities.h"
//#include <fstream>
//#include "ElementRegionT.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace multidimensionalArray;

// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager( std::string const & name,
                          ManagedGroup * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition()
{
  RegisterViewWrapper(viewKeyStruct::referencePositionString, &m_referencePosition, false );

  this->RegisterViewWrapper< array<localIndex_array> >(viewKeyStruct::elementRegionListString);
  this->RegisterViewWrapper< array<localIndex_array> >(viewKeyStruct::elementSubRegionListString);
  this->RegisterViewWrapper< array<localIndex_array> >(viewKeyStruct::elementListString);

}



// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
/*
   NodeManagerT::NodeManagerT( const NodeManagerT& init ):
   ObjectDataStructureBaseT(init),
   DataLengths()(this->m_DataLengths),
   m_refposition(NULL),
   m_displacement(NULL),
   m_incrementalDisplacement(NULL),
   m_velocity(NULL),
   m_acceleration(NULL),
   m_force(NULL),
   m_mass(NULL),
   m_toElementsRelation(init.m_toElementsRelation),
   m_nodeToFaceMap(m_UnorderedVariableOneToManyMaps["nodeToFaceMap"]),
   m_nodeToEdgeMap(m_UnorderedVariableOneToManyMaps["nodeToEdgeMap"])
   {}
 */

// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::~NodeManager()
{}


//void NodeManager::Initialize()
//{
//  this->AddKeyedDataField<FieldInfo::referencePosition>();
//
//}



void NodeManager::FillDocumentationNode()
{
  ObjectManagerBase::FillDocumentationNode();
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName( this->getCatalogName() );
  docNode->setSchemaType( "Node" );
  docNode->setShortDescription( "a node manager" );

//  docNode->AllocateChildNode( keys::elementRegionMap,
//                              keys::elementRegionMap,
//                              -1,
//                              "integer_array",
//                              "integer_array",
//                              "map to element region",
//                              "map to element region",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );
//
//  docNode->AllocateChildNode( keys::elementSubRegionMap,
//                              keys::elementSubRegionMap,
//                              -1,
//                              "integer_array",
//                              "integer_array",
//                              "map to element sub regions",
//                              "map to element sub regions",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );
//
//  docNode->AllocateChildNode( keys::elementMap,
//                              keys::elementMap,
//                              -1,
//                              "localIndex_array",
//                              "localIndex_array",
//                              "map to element in a subregion",
//                              "map to element in a subregion",
//                              "",
//                              "",
//                              1,
//                              0,
//                              0 );

  docNode->AllocateChildNode( keys::referencePositionString,
                              keys::referencePositionString,
                              -1,
                              "r1_array",
                              "r1_array",
                              "reference position of nodes",
                              "reference position of nodes",
                              "",
                              "",
                              1,
                              0,
                              0 );
}

void NodeManager::SetElementMaps( ElementRegionManager const * const elementRegionManager )
{
  array<localIndex_array> & elementRegionList = this->getReference< array<localIndex_array> >(viewKeyStruct::elementRegionListString);
  array<localIndex_array> & elementSubRegionList = this->getReference< array<localIndex_array> >(viewKeyStruct::elementSubRegionListString);
  array<localIndex_array> & elementList = this->getReference< array<localIndex_array> >(viewKeyStruct::elementListString);

  for( localIndex a=0 ; a<size() ; ++a )
  {
    elementRegionList[a].clear();
    elementSubRegionList[a].clear();
    elementList[a].clear();
  }

  for( typename dataRepository::indexType kReg=0 ; kReg<elementRegionManager->numRegions() ; ++kReg  )
  {
    ElementRegion const * const elemRegion = elementRegionManager->GetRegion(kReg);

    for( typename dataRepository::indexType kSubReg=0 ; kSubReg<elemRegion->numSubRegions() ; ++kSubReg  )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(kSubReg);

      FixedOneToManyRelation const & elemsToNodes = subRegion->getReference<FixedOneToManyRelation>(subRegion->viewKeys().nodeList);

      for( localIndex ke=0 ; ke<subRegion->size() ; ++ke )
      {
        arrayView1d<localIndex const> const nodeList = elemsToNodes[ke];
        for( localIndex a=0 ; a<elemsToNodes.size(1) ; ++a )
        {
          localIndex nodeIndex = nodeList[a];
          elementRegionList[nodeIndex].push_back( kReg );
          elementSubRegionList[nodeIndex].push_back( kSubReg );
          elementList[nodeIndex].push_back( ke );
        }
      }
    }
  }
}


void NodeManager::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementSubRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementListString));
}


localIndex NodeManager::PackUpDownMapsSize( localIndex_array const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

localIndex NodeManager::PackUpDownMaps( buffer_unit_type * & buffer,
                             localIndex_array const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

template< bool DOPACK >
localIndex NodeManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                        localIndex_array const & packList ) const
{
  localIndex packedSize = 0;


  return packedSize;
}


localIndex NodeManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                               localIndex_array const & packList )
{
  localIndex unPackedSize = 0;

  return unPackedSize;
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, std::string const &, ManagedGroup * const )

}
