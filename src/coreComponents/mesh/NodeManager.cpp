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

/**
 * @file NodeManager.hpp
 */

#include "NodeManager.hpp"
//#include "managers/DomainPartition.hpp"
#include "FaceManager.hpp"
#include "EdgeManager.hpp"
//#include "ObjectManagers/ElementManagerT.h"
//#include "Utilities/Utilities.h"
//#include <fstream>
//#include "ElementRegionT.hpp"
#include "ElementRegionManager.hpp"
#include "ToElementRelation.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"

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


  this->RegisterViewWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation, false );

  this->RegisterViewWrapper( viewKeyStruct::faceListString, &m_toFacesRelation, false );

  this->RegisterViewWrapper( viewKeyStruct::elementRegionListString,
                             &(elementRegionList()),
                             false );

  this->RegisterViewWrapper( viewKeyStruct::elementSubRegionListString,
                             &(elementSubRegionList()),
                             false );

  this->RegisterViewWrapper( viewKeyStruct::elementListString,
                             &(elementList()),
                             false );

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
                              1 );
}

//**************************************************************************************************
void NodeManager::SetEdgeMaps( EdgeManager const * const edgeManager )
{
  /// flow of function:
  /// <ol>
  /// <li> Extract edgeToNode map from the edge manager
  FixedOneToManyRelation const & edgeToNodes = edgeManager->nodeList();

  /// <li> Get number of edges
  localIndex const numEdges = edgeManager->size();

  /// <li> Loop over all edges
  /// <ol>
  for( localIndex ke=0 ; ke<numEdges ; ++ke )
  {
    localIndex const numNodes = edgeToNodes.size(1);
    /// <li> Loop over all nodes for each edge
    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      /// <ul>
      /// <li> Set the nodeToEdge relation using the current node and edge indices
      m_toEdgesRelation[a].insert(ke);
      /// </ul>
    }
  }
  /// </ol>

  /// <li> Set the related object pointer in the edge manager
  m_toEdgesRelation.SetRelatedObject( edgeManager );
  /// </ol>
}

//**************************************************************************************************
void NodeManager::SetFaceMaps( FaceManager const * const faceManager )
{
  OrderedVariableOneToManyRelation const & faceToNodes = faceManager->nodeList();
  localIndex const numFaces = faceManager->size();
  for( localIndex ke=0 ; ke<numFaces ; ++ke )
  {
    localIndex const numNodes = faceToNodes[ke].size();
    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      m_toFacesRelation[a].insert(ke);
    }
  }
  m_toFacesRelation.SetRelatedObject( faceManager );
}


//**************************************************************************************************
void NodeManager::SetElementMaps( ElementRegionManager const * const elementRegionManager )
{
  array1d<set<localIndex>> & toElementRegionList = elementRegionList();
  array1d<set<localIndex>> & toElementSubRegionList = elementSubRegionList();
  array1d<set<localIndex>> & toElementList = elementList();

  for( localIndex a=0 ; a<size() ; ++a )
  {
    toElementRegionList[a].clear();
    toElementSubRegionList[a].clear();
    toElementList[a].clear();
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
        localIndex const * const nodeList = elemsToNodes[ke];
        for( localIndex a=0 ; a<elemsToNodes.size(1) ; ++a )
        {
          localIndex nodeIndex = nodeList[a];
          toElementRegionList[nodeIndex].insert( kReg );
          toElementSubRegionList[nodeIndex].insert( kSubReg );
          toElementList[nodeIndex].insert( ke );
        }
      }
    }
  }

  this->m_toElements.setElementRegionManager( elementRegionManager );
}

//**************************************************************************************************
void NodeManager::ViewPackingExclusionList( set<localIndex> & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList(exclusionList);
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::edgeListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::faceListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementSubRegionListString));
  exclusionList.insert(this->getWrapperIndex(viewKeyStruct::elementListString));
}

//**************************************************************************************************
localIndex NodeManager::PackUpDownMapsSize( localIndex_array const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

//**************************************************************************************************
localIndex NodeManager::PackUpDownMaps( buffer_unit_type * & buffer,
                             localIndex_array const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

//**************************************************************************************************
template< bool DOPACK >
localIndex NodeManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               localIndex_array const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::edgeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toEdgesRelation,
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_toEdgesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::faceListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         m_toFacesRelation,
                                         packList,
                                         this->m_localToGlobalMap,
                                         m_toFacesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::elementListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                         this->m_toElements,
                                         packList,
                                         m_toElements.getElementRegionManager() );
  return packedSize;
}

//**************************************************************************************************
localIndex NodeManager::UnpackUpDownMaps( buffer_unit_type const * & buffer,
                               localIndex_array const & packList )
{
  localIndex unPackedSize = 0;

  string temp;
  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ASSERT( temp==viewKeyStruct::edgeListString, "")
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     this->m_globalToLocalMap,
                                     m_toEdgesRelation.RelatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ASSERT( temp==viewKeyStruct::faceListString, "")
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     this->m_globalToLocalMap,
                                     m_toFacesRelation.RelatedObjectGlobalToLocal() );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ASSERT( temp==viewKeyStruct::elementListString, "")
  unPackedSize += bufferOps::Unpack( buffer,
                                     this->m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager() );

  return unPackedSize;
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, std::string const &, ManagedGroup * const )

}
