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
#include "ToElementRelation.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{
using namespace dataRepository;

// *********************************************************************************************************************
/**
 * @return
 */
NodeManager::NodeManager( std::string const & name,
                          Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString, &m_referencePosition );

  this->registerWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation );

  this->registerWrapper( viewKeyStruct::faceListString, &m_toFacesRelation );

  this->registerWrapper( viewKeyStruct::elementRegionListString, &elementRegionList() );

  this->registerWrapper( viewKeyStruct::elementSubRegionListString, &elementSubRegionList() );

  this->registerWrapper( viewKeyStruct::elementListString, &elementList() );

}



// *********************************************************************************************************************

// *********************************************************************************************************************
/**
 * @return
 */
NodeManager::~NodeManager()
{}


//**************************************************************************************************
void NodeManager::SetEdgeMaps( EdgeManager const * const edgeManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager->nodeList();
  localIndex const numEdges = edgeToNodeMap.size( 0 );
  localIndex const numNodes = size();

  ArrayOfArrays< localIndex > toEdgesTemp( numNodes, edgeManager->maxEdgesPerNode() );

  forAll< parallelHostPolicy >( numEdges, [&]( localIndex const edgeID )
  {
    toEdgesTemp.atomicAppendToArray( RAJA::auto_atomic{}, edgeToNodeMap( edgeID, 0 ), edgeID );
    toEdgesTemp.atomicAppendToArray( RAJA::auto_atomic{}, edgeToNodeMap( edgeID, 1 ), edgeID );
  } );

  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeEdges( 0 );
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    totalNodeEdges += toEdgesTemp.sizeOfArray( nodeID );
  } );

  m_toEdgesRelation.resize( 0 );
  m_toEdgesRelation.reserve( numNodes );
  m_toEdgesRelation.reserveValues( totalNodeEdges.get() + numNodes * GetEdgeMapOverallocation() );
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    m_toEdgesRelation.appendSet( toEdgesTemp.sizeOfArray( nodeID ) + GetEdgeMapOverallocation() );
  }

  ArrayOfSetsView< localIndex > const & toEdgesView = m_toEdgesRelation;
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const edges = toEdgesTemp[ nodeID ];
    localIndex const numNodeEdges = toEdgesTemp.sizeOfArray( nodeID );
    std::sort( edges, edges + numNodeEdges );
    toEdgesView.insertSortedIntoSet( nodeID, edges, numNodeEdges );
  } );

  m_toEdgesRelation.SetRelatedObject( edgeManager );
}

//**************************************************************************************************
void NodeManager::SetFaceMaps( FaceManager const * const faceManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager->nodeList();
  localIndex const numFaces = faceToNodes.size();
  localIndex const numNodes = size();

  ArrayOfArrays< localIndex > toFacesTemp( numNodes, faceManager->maxFacesPerNode() );

  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const faceID )
  {
    localIndex const numFaceNodes = faceToNodes.sizeOfArray( faceID );
    for( localIndex a = 0; a < numFaceNodes; ++a )
    {
      toFacesTemp.atomicAppendToArray( RAJA::auto_atomic{}, faceToNodes( faceID, a ), faceID );
    }
  } );

  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeFaces( 0 );
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    totalNodeFaces += toFacesTemp.sizeOfArray( nodeID );
  } );

  m_toFacesRelation.resize( 0 );
  m_toFacesRelation.reserve( numNodes );
  m_toFacesRelation.reserveValues( totalNodeFaces.get() + numNodes * GetFaceMapOverallocation() );
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    m_toFacesRelation.appendSet( toFacesTemp.sizeOfArray( nodeID ) + GetFaceMapOverallocation() );
  }

  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const faces = toFacesTemp[ nodeID ];
    localIndex const numNodeFaces = toFacesTemp.sizeOfArray( nodeID );
    std::sort( faces, faces + numNodeFaces );
    m_toFacesRelation.insertSortedIntoSet( nodeID, faces, numNodeFaces );
  } );

  m_toFacesRelation.SetRelatedObject( faceManager );
}


//**************************************************************************************************
void NodeManager::SetElementMaps( ElementRegionManager const * const elementRegionManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< localIndex > & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays< localIndex > & toElementSubRegionList = m_toElements.m_toElementSubRegion;
  ArrayOfArrays< localIndex > & toElementList = m_toElements.m_toElementIndex;


  // This sets the capacity of each sub-array based on the maximum number of
  // elements attached to a node. This is an over-allocation, so we need to
  // compress the arrays as part of the initialization.

  array1d< localIndex > sizeOfArrays( size());

  elementRegionManager->
    forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemToNodeMap = subRegion.nodeList();
    for( localIndex k=0; k<subRegion.size(); ++k )
    {
      for( localIndex a=0; a<subRegion.numIndependentNodesPerElement(); ++a )
      {
        localIndex nodeIndex = elemToNodeMap( k, a );
        ++sizeOfArrays[nodeIndex];
      }
    }
  } );
  localIndex const arrayCapacity = *(std::max_element( sizeOfArrays.begin(), sizeOfArrays.end() ));

  toElementRegionList.resize( 0 );
  toElementRegionList.resize( size(), arrayCapacity );
  toElementSubRegionList.resize( 0 );
  toElementSubRegionList.resize( size(), arrayCapacity );
  toElementList.resize( 0 );
  toElementList.resize( size(), arrayCapacity );

  for( typename dataRepository::indexType kReg=0; kReg<elementRegionManager->numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = elementRegionManager->GetRegion( kReg );

    elemRegion->forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const kSubReg,
                                                                        CellElementSubRegion const & subRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemToNodeMap = subRegion.nodeList();

      for( localIndex k=0; k<subRegion.size(); ++k )
      {
        for( localIndex a=0; a<subRegion.numIndependentNodesPerElement(); ++a )
        {
          localIndex nodeIndex = elemToNodeMap( k, a );

          toElementRegionList.appendToArray( nodeIndex, kReg );
          toElementSubRegionList.appendToArray( nodeIndex, kSubReg );
          toElementList.appendToArray( nodeIndex, k );
        }
      }
    } );
  }

  this->m_toElements.setElementRegionManager( elementRegionManager );
}

void NodeManager::CompressRelationMaps()
{
  //GEOSX_MARK_FUNCTION;
  m_toEdgesRelation.compress();
  m_toFacesRelation.compress();
  m_toElements.m_toElementRegion.compress();
  m_toElements.m_toElementSubRegion.compress();
  m_toElements.m_toElementIndex.compress();
}

//**************************************************************************************************
void NodeManager::ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementRegionListString ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementSubRegionListString ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementListString ));

  if( this->hasWrapper( "usedFaces" ) )
  {
    exclusionList.insert( this->getWrapperIndex( "usedFaces" ));
  }
}

//**************************************************************************************************
localIndex NodeManager::PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}

//**************************************************************************************************
localIndex NodeManager::PackUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}

//**************************************************************************************************
template< bool DOPACK >
localIndex NodeManager::PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::edgeListString ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toEdgesRelation.toArrayOfArraysView(),
                                           m_unmappedGlobalIndicesInToEdges,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_toEdgesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::faceListString ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toFacesRelation.toArrayOfArraysView(),
                                           m_unmappedGlobalIndicesInToFaces,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_toFacesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::elementListString ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
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
  GEOSX_ERROR_IF( temp != viewKeyStruct::edgeListString, "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.RelatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::faceListString, "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.RelatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::elementListString, "" );
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
                                    m_toEdgesRelation.RelatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toFacesRelation,
                                    m_toFacesRelation.RelatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );

}

void NodeManager::depopulateUpMaps( std::set< localIndex > const & receivedNodes,
                                    array2d< localIndex > const & edgesToNodes,
                                    ArrayOfArraysView< localIndex const > const & facesToNodes,
                                    ElementRegionManager const & elemRegionManager )
{

  ObjectManagerBase::CleanUpMap( receivedNodes, m_toEdgesRelation, edgesToNodes );
  ObjectManagerBase::CleanUpMap( receivedNodes, m_toFacesRelation, facesToNodes );

  for( auto const & targetIndex : receivedNodes )
  {
    std::set< std::tuple< localIndex, localIndex, localIndex > > eraseList;
    for( localIndex k=0; k<m_toElements.m_toElementRegion.sizeOfArray( targetIndex ); ++k )
    {
      localIndex const elemRegionIndex    = m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex = m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex          = m_toElements.m_toElementIndex[targetIndex][k];

      CellElementSubRegion const * subRegion = elemRegionManager.GetRegion( elemRegionIndex )->
                                                 GetSubRegion< CellElementSubRegion >( elemSubRegionIndex );
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & downmap = subRegion->nodeList();
      bool hasTargetIndex = false;

      for( localIndex a=0; a<downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap( elemIndex, a );
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert( std::make_tuple( elemRegionIndex, elemSubRegionIndex, elemIndex ) );
      }
    }
    for( auto const & val : eraseList )
    {
      erase( m_toElements, targetIndex, std::get< 0 >( val ), std::get< 1 >( val ), std::get< 2 >( val ) );
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, std::string const &, Group * const )

}
