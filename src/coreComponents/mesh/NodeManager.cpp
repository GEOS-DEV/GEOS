/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::~NodeManager()
{}


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
      m_toEdgesRelation[edgeToNodes[ke][a]].insert(ke);
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
      m_toFacesRelation[faceToNodes[ke][a]].insert(ke);
    }
  }
  m_toFacesRelation.SetRelatedObject( faceManager );
}


//**************************************************************************************************
void NodeManager::SetElementMaps( ElementRegionManager const * const elementRegionManager )
{
  array1d<localIndex_array> & toElementRegionList = elementRegionList();
  array1d<localIndex_array> & toElementSubRegionList = elementSubRegionList();
  array1d<localIndex_array> & toElementList = elementList();

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
      ElementSubRegionBase const * const subRegion = elemRegion->GetGroup(ElementRegion::viewKeyStruct::elementSubRegions)->GetGroup<ElementSubRegionBase>(kSubReg);


      for( localIndex ke=0 ; ke<subRegion->size() ; ++ke )
      {
        arraySlice1d<localIndex const> const elemToNodes = subRegion->nodeList(ke);
        for( localIndex a=0 ; a<subRegion->numNodesPerElement() ; ++a )
        {
          localIndex nodeIndex = elemToNodes[a];
          toElementRegionList[nodeIndex].push_back( kReg );
          toElementSubRegionList[nodeIndex].push_back( kSubReg );
          toElementList[nodeIndex].push_back( ke );
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

  if( this->hasView( "usedFaces" ) )
  {
    exclusionList.insert(this->getWrapperIndex("usedFaces"));
  }
}

//**************************************************************************************************
localIndex NodeManager::PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate<false>( junk, packList );
}

//**************************************************************************************************
localIndex NodeManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d<localIndex const> const & packList ) const
{
  return PackUpDownMapsPrivate<true>( buffer, packList );
}

//**************************************************************************************************
template< bool DOPACK >
localIndex NodeManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d<localIndex const> const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::edgeListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                       m_toEdgesRelation,
                                       m_unmappedGlobalIndicesInToEdges,
                                       packList,
                                       this->m_localToGlobalMap,
                                       m_toEdgesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack<DOPACK>( buffer, string(viewKeyStruct::faceListString) );
  packedSize += bufferOps::Pack<DOPACK>( buffer,
                                       m_toFacesRelation,
                                       m_unmappedGlobalIndicesInToFaces,
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
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const )
{
  localIndex unPackedSize = 0;

  string temp;
  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ERROR_IF( temp != viewKeyStruct::edgeListString, "");
  unPackedSize += bufferOps::Unpack( buffer,
                                   m_toEdgesRelation,
                                   packList,
                                   m_unmappedGlobalIndicesInToEdges,
                                   this->m_globalToLocalMap,
                                   m_toEdgesRelation.RelatedObjectGlobalToLocal(),
                                   overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ERROR_IF( temp != viewKeyStruct::faceListString, "");
  unPackedSize += bufferOps::Unpack( buffer,
                                   m_toFacesRelation,
                                   packList,
                                   m_unmappedGlobalIndicesInToFaces,
                                   this->m_globalToLocalMap,
                                   m_toFacesRelation.RelatedObjectGlobalToLocal(),
                                   overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOS_ERROR_IF( temp != viewKeyStruct::elementListString, "");
  unPackedSize += bufferOps::Unpack( buffer,
                                   this->m_toElements,
                                   packList,
                                   m_toElements.getElementRegionManager(),
                                   overwriteUpMaps );

  return unPackedSize;
}

void NodeManager::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );

}

void NodeManager::depopulateUpMaps( std::set<localIndex> const & receivedNodes,
                                    array2d< localIndex > const & edgesToNodes,
                                    array1d< array1d< localIndex > > const & facesToNodes,
                                    ElementRegionManager const & elemRegionManager )
{

  ObjectManagerBase::CleanUpMap( receivedNodes, m_toEdgesRelation, edgesToNodes );
  ObjectManagerBase::CleanUpMap( receivedNodes, m_toFacesRelation, facesToNodes );

  for( auto const & targetIndex : receivedNodes )
  {
    set<std::tuple<localIndex,localIndex,localIndex> > eraseList;
    for( localIndex k=0 ; k<m_toElements.m_toElementRegion[targetIndex].size() ; ++k )
    {
      localIndex const elemRegionIndex    = m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex = m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex          = m_toElements.m_toElementIndex[targetIndex][k];

      CellElementSubRegion const * subRegion = elemRegionManager.GetRegion(elemRegionIndex)->
                                               GetSubRegion<CellElementSubRegion>(elemSubRegionIndex);
      array2d<localIndex> const & downmap = subRegion->nodeList();
      bool hasTargetIndex = false;

      for( localIndex a=0 ; a<downmap.size(1) ; ++a )
      {
        localIndex const compositeLocalIndex = downmap[elemIndex][a];
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert(std::make_tuple(elemRegionIndex, elemSubRegionIndex, elemIndex) );
      }
    }
    for( auto const & val : eraseList )
    {
      erase( m_toElements, targetIndex, std::get<0>(val), std::get<1>(val), std::get<2>(val) );
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, std::string const &, ManagedGroup * const )

}
