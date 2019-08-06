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
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager->nodeList();
  localIndex const totalNumEdges = edgeToNodeMap.size( 0 );
  localIndex const totalNumNodes = size();

  constexpr int MAX_EDGES_PER_NODE = 10;
  ArrayOfArrays< localIndex > toEdgesTemp( totalNumNodes, MAX_EDGES_PER_NODE );

  forall_in_range< parallelHostPolicy >( 0, totalNumEdges, [&]( localIndex const edgeID )
  {
    toEdgesTemp.atomicAppendToArray( RAJA::atomic::auto_atomic{}, edgeToNodeMap( edgeID, 0 ), edgeID );
    toEdgesTemp.atomicAppendToArray( RAJA::atomic::auto_atomic{}, edgeToNodeMap( edgeID, 1 ), edgeID );
  } );

  GEOSX_MARK_BEGIN("Reserving space in m_toEdgesRelation");
  for ( localIndex i = 0; i < totalNumNodes; ++i )
  {
    m_toEdgesRelation[ i ].reserve( toEdgesTemp.sizeOfArray( i ) );
  }
  GEOSX_MARK_END("Reserving space in m_toEdgesRelation");

  forall_in_range< parallelHostPolicy >( 0, totalNumNodes, [&]( localIndex const nodeID )
  {
    localIndex * const edges = toEdgesTemp[ nodeID ];
    localIndex const numEdges = toEdgesTemp.sizeOfArray( nodeID );
    std::sort( edges, edges + numEdges );
    m_toEdgesRelation[ nodeID ].insertSorted( edges, numEdges );
  } );

  m_toEdgesRelation.SetRelatedObject( edgeManager );
}

//**************************************************************************************************
void NodeManager::SetFaceMaps( FaceManager const * const faceManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< arrayView1d< localIndex const > const > const & faceToNodes = faceManager->nodeList().toViewConst();
  localIndex const totalNumFaces = faceToNodes.size();
  localIndex const totalNumNodes = size();

  constexpr int MAX_FACES_PER_NODE = 20;
  ArrayOfArrays< localIndex > toFacesTemp( totalNumNodes, MAX_FACES_PER_NODE );

  forall_in_range< parallelHostPolicy >( 0, totalNumFaces, [&]( localIndex const faceID )
  {
    localIndex const numFaceNodes = faceToNodes[ faceID ].size();
    for ( localIndex a = 0; a < numFaceNodes; ++a )
    {
      toFacesTemp.atomicAppendToArray( RAJA::atomic::auto_atomic{}, faceToNodes[ faceID ][ a ], faceID );
    }
  } );

  GEOSX_MARK_BEGIN("Reserving space in m_toFacesRelation");
  for ( localIndex i = 0; i < totalNumNodes; ++i )
  {
    m_toFacesRelation[ i ].reserve( toFacesTemp.sizeOfArray( i ) );
  }
  GEOSX_MARK_END("Reserving space in m_toFacesRelation");

  forall_in_range< parallelHostPolicy >( 0, totalNumNodes, [&]( localIndex const nodeID )
  {
    localIndex * const faces = toFacesTemp[ nodeID ];
    localIndex const numFaces = toFacesTemp.sizeOfArray( nodeID );
    std::sort( faces, faces + numFaces );
    m_toFacesRelation[ nodeID ].insertSorted( faces, numFaces );
  } );

  m_toFacesRelation.SetRelatedObject( faceManager );
}


//**************************************************************************************************
void NodeManager::SetElementMaps( ElementRegionManager const * const elementRegionManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays<localIndex> & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays<localIndex> & toElementSubRegionList = m_toElements.m_toElementSubRegion;
  ArrayOfArrays<localIndex> & toElementList = m_toElements.m_toElementIndex;

  // This sets the capacity of each sub-array to 10. If this is using a bunch of memory
  // add a compress + shrink method to ArrayOfArrays that we can call afterwards.
  toElementRegionList.resize(0);
  toElementRegionList.resize(size(), 10);

  toElementSubRegionList.resize(0);
  toElementSubRegionList.resize(size(), 10);

  toElementList.resize(0);
  toElementList.resize(size(), 10);

  for( typename dataRepository::indexType kReg=0 ; kReg<elementRegionManager->numRegions() ; ++kReg )
  {
    ElementRegion const * const elemRegion = elementRegionManager->GetRegion(kReg);

    elemRegion->forElementSubRegionsIndex<CellElementSubRegion>( [&]( localIndex const kSubReg,
                                                                      CellElementSubRegion const * const subRegion )
    {
      for( localIndex ke=0 ; ke<subRegion->size() ; ++ke )
      {
        arraySlice1d<localIndex const> const elemToNodes = subRegion->nodeList(ke);
        
        for( localIndex a=0 ; a<subRegion->numNodesPerElement() ; ++a )
        {
          localIndex nodeIndex = elemToNodes[a];
          toElementRegionList.appendToArray( nodeIndex, kReg );
          toElementSubRegionList.appendToArray( nodeIndex, kSubReg );
          toElementList.appendToArray( nodeIndex, ke );
        }
      }
    });
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
    std::set<std::tuple<localIndex,localIndex,localIndex> > eraseList;
    for( localIndex k=0 ; k<m_toElements.m_toElementRegion.sizeOfArray(targetIndex) ; ++k )
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
