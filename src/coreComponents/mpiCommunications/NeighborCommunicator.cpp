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
 * @file NeighborCommunicator.cpp
 *
 */

#include "mpiCommunications/NeighborCommunicator.hpp"

#include "common/TimingMacros.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "mesh/MeshLevel.hpp"
#include <sys/time.h>

namespace geosx
{

using namespace dataRepository;

NeighborCommunicator::NeighborCommunicator():
  m_neighborRank( -1 ),
  m_sendBufferSize(),
  m_receiveBufferSize(),
  m_sendBuffer(),
  m_receiveBuffer(),
  m_mpiSendBufferRequest(),
  m_mpiRecvBufferRequest(),
  m_mpiSendBufferStatus(),
  m_mpiRecvBufferStatus()
{}

//NeighborCommunicator::~NeighborCommunicator()
//{
//}

void NeighborCommunicator::MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                                             int const sendSize,
                                             MPI_Request& sendRequest,
                                             buffer_unit_type * const receiveBuffer,
                                             int const receiveSize,
                                             MPI_Request& receiveRequest,
                                             int const commID,
                                             MPI_Comm mpiComm )
{
  int const sendTag = CommTag( Rank(), m_neighborRank, commID );
  //m_rank * m_size + m_neighborRank + m_size*m_size*commID;
  MpiWrapper::iSend( const_cast<buffer_unit_type*>(sendBuffer),
             sendSize,
             m_neighborRank,
             sendTag,
             mpiComm,
             &sendRequest );

  int const receiveTag = CommTag( m_neighborRank, Rank(), commID );
  //m_neighborRank * m_size + m_rank + m_size*m_size*commID;
  MpiWrapper::iRecv( receiveBuffer,
             receiveSize,
             m_neighborRank,
             receiveTag,
             mpiComm,
             &receiveRequest );
}

void NeighborCommunicator::MPI_iSendReceiveBufferSizes( int const commID,
                                                        MPI_Comm mpiComm )
{
  MPI_iSendReceiveBufferSizes( commID,
                               m_mpiRecvBufferRequest[commID],
                               m_mpiSendBufferRequest[commID],
                               mpiComm );
}

void NeighborCommunicator::MPI_iSendReceiveBufferSizes( int const commID,
                                                        MPI_Request& mpiSendRequest,
                                                        MPI_Request& mpiRecvRequest,
                                                        MPI_Comm mpiComm )
{
//  m_sendBufferSize[commID] = integer_conversion<int>( m_sendBuffer[commID].size());
  MPI_iSendReceive( &m_sendBufferSize[commID], 1, mpiSendRequest,
                    &m_receiveBufferSize[commID],
                    1, mpiRecvRequest,
                    commID,
                    mpiComm );
}

void NeighborCommunicator::MPI_iSendReceiveBuffers( int const commID,
                                                    MPI_Comm mpiComm )
{
  MPI_iSendReceiveBuffers( commID,
                           m_mpiSendBufferRequest[commID],
                           m_mpiRecvBufferRequest[commID],
                           mpiComm );
}

void NeighborCommunicator::MPI_iSendReceiveBuffers( int const commID,
                                                    MPI_Request& mpiSendRequest,
                                                    MPI_Request& mpiRecvRequest,
                                                    MPI_Comm mpiComm )
{
  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( m_sendBuffer[commID].data(),
                    integer_conversion<int>( m_sendBuffer[commID].size()),
                    mpiSendRequest,
                    m_receiveBuffer[commID].data(),
                    integer_conversion<int>( m_receiveBuffer[commID].size()),
                    mpiRecvRequest,
                    commID,
                    mpiComm );

}

void NeighborCommunicator::MPI_iSendReceive( int const commID,
                                             MPI_Comm mpiComm )
{
  MPI_iSendReceive( commID,
                    m_mpiSendBufferRequest[commID],
                    m_mpiRecvBufferRequest[commID],
                    mpiComm );
}

void NeighborCommunicator::MPI_iSendReceive( int const commID,
                                             MPI_Request& mpiSendRequest,
                                             MPI_Request& mpiRecvRequest,
                                             MPI_Comm mpiComm )
{
  MPI_iSendReceiveBufferSizes( commID, mpiComm );

  MpiWrapper::Waitall( 1, &( m_mpiRecvBufferRequest[commID] ), &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Waitall( 1, &( m_mpiSendBufferRequest[commID] ), &( m_mpiSendBufferStatus[commID] ) );

  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( m_sendBuffer[commID].data(),
                    m_sendBufferSize[commID],
                    mpiSendRequest,
                    m_receiveBuffer[commID].data(),
                    m_receiveBufferSize[commID],
                    mpiRecvRequest,
                    commID,
                    mpiComm );
}

void NeighborCommunicator::MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                                             int const sendSize,
                                             int const commID,
                                             MPI_Comm mpiComm )
{
  MPI_iSendReceive( &sendSize,
                    1,
                    m_mpiSendBufferRequest[commID],
                    &m_receiveBufferSize[commID],
                    1,
                    m_mpiRecvBufferRequest[commID],
                    commID,
                    mpiComm );

  MpiWrapper::Waitall( 1, &( m_mpiRecvBufferRequest[commID] ), &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Waitall( 1, &( m_mpiSendBufferRequest[commID] ), &( m_mpiSendBufferStatus[commID] ) );

  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( sendBuffer,
                    sendSize,
                    m_mpiSendBufferRequest[commID],
                    m_receiveBuffer[commID].data(),
                    m_receiveBufferSize[commID],
                    m_mpiRecvBufferRequest[commID],
                    commID,
                    mpiComm );
}

void NeighborCommunicator::MPI_WaitAll( int const GEOSX_UNUSED_ARG( commID ),
                                        MPI_Request& mpiSendRequest,
                                        MPI_Status& mpiSendStatus,
                                        MPI_Request& mpiRecvRequest,
                                        MPI_Status& mpiReceiveStatus )

{
  MpiWrapper::Waitall( 1, &mpiRecvRequest, &mpiReceiveStatus );
  MpiWrapper::Waitall( 1, &mpiSendRequest, &mpiSendStatus );
}

void NeighborCommunicator::MPI_WaitAll( int const commID )
{
  MpiWrapper::Waitall( 1, &( m_mpiRecvBufferRequest[commID] ), &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Waitall( 1, &( m_mpiSendBufferRequest[commID] ), &( m_mpiSendBufferStatus[commID] ) );
}

int NeighborCommunicator::Rank()
{
  return MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
}

int NeighborCommunicator::MPISize()
{
  return MpiWrapper::Comm_size( MPI_COMM_GEOSX );
}

void NeighborCommunicator::Clear()
{
  for( int i = 0 ; i < maxComm ; ++i )
  {
    m_sendBuffer[i].clear();
    m_receiveBuffer[i].clear();
  }
}

void NeighborCommunicator::AddNeighborGroupToMesh( MeshLevel * const mesh ) const
{
  Group * neighborGroups[100];
  localIndex numNeighborGroups = 0;

  ObjectManagerBase * const nodeManager = mesh->getNodeManager();
  neighborGroups[numNeighborGroups++] = nodeManager->
                                        GetGroup( nodeManager->m_ObjectManagerBaseGroupKeys.neighborData )->
                                        RegisterGroup( std::to_string( this->m_neighborRank ));

  ObjectManagerBase * const edgeManager = mesh->getEdgeManager();
  neighborGroups[numNeighborGroups++] = edgeManager->
                                        GetGroup( edgeManager->m_ObjectManagerBaseGroupKeys.neighborData )->
                                        RegisterGroup( std::to_string( this->m_neighborRank ));

  ObjectManagerBase * const faceManager = mesh->getFaceManager();
  neighborGroups[numNeighborGroups++] = faceManager->
                                        GetGroup( faceManager->m_ObjectManagerBaseGroupKeys.neighborData )->
                                        RegisterGroup( std::to_string( this->m_neighborRank ));

  ElementRegionManager * const elemManager = mesh->getElemManager();
  elemManager->forElementSubRegions( [&]( Group * const elementSubRegion ) -> void
  {
    neighborGroups[numNeighborGroups++] = elementSubRegion->
                                          GetGroup( faceManager->m_ObjectManagerBaseGroupKeys.neighborData )->
                                          RegisterGroup( std::to_string( this->m_neighborRank ));
  } );

  for( localIndex a=0 ; a<numNeighborGroups ; ++a )
  {
    neighborGroups[a]->
    registerWrapper<localIndex_array>( ObjectManagerBase::viewKeyStruct::matchedPartitionBoundaryObjectsString )->
    setSizedFromParent( 0 );

    neighborGroups[a]->
    registerWrapper<localIndex_array>( ObjectManagerBase::viewKeyStruct::ghostsToSendString )->
    setSizedFromParent( 0 );

    neighborGroups[a]->
    registerWrapper<localIndex_array>( ObjectManagerBase::viewKeyStruct::ghostsToReceiveString )->
    setSizedFromParent( 0 );

    neighborGroups[a]->
    registerWrapper<localIndex_array>( ObjectManagerBase::viewKeyStruct::adjacencyListString )->
    setSizedFromParent( 0 );

  }



}

void NeighborCommunicator::FindAndPackGhosts( bool const GEOSX_UNUSED_ARG( contactActive ),
                                              integer const depth,
                                              MeshLevel * const mesh,
                                              int const commID )
{
  GEOSX_MARK_FUNCTION;
  NodeManager & nodeManager = *(mesh->getNodeManager());
  EdgeManager & edgeManager = *(mesh->getEdgeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());

  Group * const nodeNeighborData = nodeManager.
                                          GetGroup( nodeManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const edgeNeighborData = edgeManager.
                                          GetGroup( edgeManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const faceNeighborData = faceManager.
                                          GetGroup( faceManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  localIndex_array & nodeAdjacencyList = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.adjacencyList );
  localIndex_array & edgeAdjacencyList = edgeNeighborData->getReference<localIndex_array>( edgeManager.viewKeys.adjacencyList );
  localIndex_array & faceAdjacencyList = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.adjacencyList );

  {
    ElementRegionManager::ElementViewAccessor<ReferenceWrapper<localIndex_array>>
    elementAdjacencyList =
      elemManager.ConstructReferenceAccessor<localIndex_array>( ObjectManagerBase::viewKeyStruct::adjacencyListString,
                                                           std::to_string( this->m_neighborRank ) );

    mesh->GenerateAdjacencyLists( nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.matchedPartitionBoundaryObjects ),
                                  nodeAdjacencyList,
                                  edgeAdjacencyList,
                                  faceAdjacencyList,
                                  elementAdjacencyList,
                                  depth );
  }

  ElementRegionManager::ElementViewAccessor<arrayView1d<localIndex>> const elementAdjacencyList =
    elemManager.ConstructViewAccessor<array1d<localIndex>, arrayView1d<localIndex>>( ObjectManagerBase::viewKeyStruct::adjacencyListString,
                                                           std::to_string( this->m_neighborRank ) );

  int bufferSize = 0;

  bufferSize += nodeManager.PackGlobalMapsSize( nodeAdjacencyList, 0 );
  bufferSize += edgeManager.PackGlobalMapsSize( edgeAdjacencyList, 0 );
  bufferSize += faceManager.PackGlobalMapsSize( faceAdjacencyList, 0 );
  bufferSize += elemManager.PackGlobalMapsSize( elementAdjacencyList );

  bufferSize += nodeManager.PackUpDownMapsSize( nodeAdjacencyList );
  bufferSize += edgeManager.PackUpDownMapsSize( edgeAdjacencyList );
  bufferSize += faceManager.PackUpDownMapsSize( faceAdjacencyList );
  bufferSize += elemManager.PackUpDownMapsSize( elementAdjacencyList );

  bufferSize += nodeManager.PackSize( {}, nodeAdjacencyList, 0 );
  bufferSize += edgeManager.PackSize( {}, edgeAdjacencyList, 0 );
  bufferSize += faceManager.PackSize( {}, faceAdjacencyList, 0 );
  bufferSize += elemManager.PackSize( {}, elementAdjacencyList );

  resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = SendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = 0;

  packedSize += nodeManager.PackGlobalMaps( sendBufferPtr, nodeAdjacencyList, 0 );
  packedSize += edgeManager.PackGlobalMaps( sendBufferPtr, edgeAdjacencyList, 0 );
  packedSize += faceManager.PackGlobalMaps( sendBufferPtr, faceAdjacencyList, 0 );
  packedSize += elemManager.PackGlobalMaps( sendBufferPtr, elementAdjacencyList );

  packedSize += nodeManager.PackUpDownMaps( sendBufferPtr, nodeAdjacencyList );
  packedSize += edgeManager.PackUpDownMaps( sendBufferPtr, edgeAdjacencyList );
  packedSize += faceManager.PackUpDownMaps( sendBufferPtr, faceAdjacencyList );
  packedSize += elemManager.PackUpDownMaps( sendBufferPtr, elementAdjacencyList );

  packedSize += nodeManager.Pack( sendBufferPtr, {}, nodeAdjacencyList, 0 );
  packedSize += edgeManager.Pack( sendBufferPtr, {}, edgeAdjacencyList, 0 );
  packedSize += faceManager.Pack( sendBufferPtr, {}, faceAdjacencyList, 0 );
  packedSize += elemManager.Pack( sendBufferPtr, {}, elementAdjacencyList );


  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );

  this->MPI_iSendReceive( commID, MPI_COMM_GEOSX );
}

void NeighborCommunicator::UnpackGhosts( MeshLevel * const mesh,
                                         int const commID )
{
  GEOSX_MARK_FUNCTION;
  NodeManager & nodeManager = *(mesh->getNodeManager());
  EdgeManager & edgeManager = *(mesh->getEdgeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());

  buffer_type const & receiveBuffer = ReceiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  int unpackedSize = 0;

  localIndex_array nodeUnpackList;
  unpackedSize += nodeManager.UnpackGlobalMaps( receiveBufferPtr,
                                                nodeUnpackList,
                                                0 );

  localIndex_array edgeUnpackList;
  unpackedSize += edgeManager.UnpackGlobalMaps( receiveBufferPtr,
                                                edgeUnpackList,
                                                0 );

  localIndex_array faceUnpackList;
  unpackedSize += faceManager.UnpackGlobalMaps( receiveBufferPtr,
                                                faceUnpackList,
                                                0 );

    ElementRegionManager::ElementReferenceAccessor<localIndex_array>
    elementAdjacencyReceiveListArray =
      elemManager.ConstructReferenceAccessor<localIndex_array>( ObjectManagerBase::viewKeyStruct::ghostsToReceiveString,
                                                                std::to_string( this->m_neighborRank ) );
    unpackedSize += elemManager.UnpackGlobalMaps( receiveBufferPtr,
                                                  elementAdjacencyReceiveListArray );

  ElementRegionManager::ElementViewAccessor<arrayView1d<localIndex>>
  elementAdjacencyReceiveList =
    elemManager.ConstructViewAccessor<array1d<localIndex>, arrayView1d<localIndex>>( ObjectManagerBase::viewKeyStruct::ghostsToReceiveString,
                                                                                     std::to_string( this->m_neighborRank ) );

  unpackedSize += nodeManager.UnpackUpDownMaps( receiveBufferPtr, nodeUnpackList, false, false );
  unpackedSize += edgeManager.UnpackUpDownMaps( receiveBufferPtr, edgeUnpackList, false, false );
  unpackedSize += faceManager.UnpackUpDownMaps( receiveBufferPtr, faceUnpackList, false, false );
  unpackedSize += elemManager.UnpackUpDownMaps( receiveBufferPtr, elementAdjacencyReceiveListArray, false );


  unpackedSize += nodeManager.Unpack( receiveBufferPtr, nodeUnpackList, 0 );
  unpackedSize += edgeManager.Unpack( receiveBufferPtr, edgeUnpackList, 0 );
  unpackedSize += faceManager.Unpack( receiveBufferPtr, faceUnpackList, 0 );
  unpackedSize += elemManager.Unpack( receiveBufferPtr, elementAdjacencyReceiveList );

}




void NeighborCommunicator::RebuildSyncLists( MeshLevel * const mesh,
                                             int const commID )
{
  NodeManager & nodeManager = *(mesh->getNodeManager());
  EdgeManager & edgeManager = *(mesh->getEdgeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());

  Group * const nodeNeighborData = nodeManager.GetGroup( nodeManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const edgeNeighborData = edgeManager.GetGroup( edgeManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const faceNeighborData = faceManager.
                                          GetGroup( faceManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );


  localIndex_array & nodeGhostsToSend = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );
  localIndex_array & edgeGhostsToSend = edgeNeighborData->getReference<localIndex_array>( edgeManager.viewKeys.ghostsToSend );
  localIndex_array & faceGhostsToSend = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToSend );

  localIndex_array const & nodeGhostsToReceive = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToReceive );
  localIndex_array const & edgeGhostsToReceive = edgeNeighborData->getReference<localIndex_array>( edgeManager.viewKeys.ghostsToReceive );
  localIndex_array const & faceGhostsToReceive = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToReceive );

  ElementRegionManager::ElementViewAccessor<arrayView1d<localIndex>>
  elementGhostToReceive =
    elemManager.ConstructViewAccessor<array1d<localIndex>, arrayView1d<localIndex>>( ObjectManagerBase::
                                                         viewKeyStruct::ghostsToReceiveString,
                                                         std::to_string( this->m_neighborRank ) );


  buffer_type & sendBuffer = SendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int bufferSize = 0;
  bufferSize += bufferOps::Pack<false>( sendBufferPtr,
                                        nodeGhostsToReceive.toSliceConst(),
                                        nullptr,
                                        nodeGhostsToReceive.size(),
                                        nodeManager.m_localToGlobalMap.toSliceConst() );

  bufferSize += bufferOps::Pack<false>( sendBufferPtr,
                                        edgeGhostsToReceive.toSliceConst(),
                                        nullptr,
                                        edgeGhostsToReceive.size(),
                                        edgeManager.m_localToGlobalMap.toSliceConst() );

  bufferSize += bufferOps::Pack<false>( sendBufferPtr,
                                        faceGhostsToReceive.toSliceConst(),
                                        nullptr,
                                        faceGhostsToReceive.size(),
                                        faceManager.m_localToGlobalMap.toSliceConst() );

  for( localIndex er=0 ; er<elemManager.numRegions() ; ++er )
  {
    ElementRegionBase const * const elemRegion = elemManager.GetRegion( er );
    elemRegion->forElementSubRegionsIndex( [&]( localIndex const esr,
                                                ElementSubRegionBase const * const subRegion )
    {
      bufferSize+= bufferOps::Pack<false>( sendBufferPtr,
                                           elementGhostToReceive[er][esr].toSliceConst(),
                                           nullptr,
                                           elementGhostToReceive[er][esr].size(),
                                           subRegion->m_localToGlobalMap.toSliceConst() );
    });
  }

  this->resizeSendBuffer( commID, bufferSize );

  int packedSize = 0;
  packedSize += bufferOps::Pack<true>( sendBufferPtr,
                                       nodeGhostsToReceive.toSliceConst(),
                                       nullptr,
                                       nodeGhostsToReceive.size(),
                                       nodeManager.m_localToGlobalMap.toSliceConst() );

  packedSize += bufferOps::Pack<true>( sendBufferPtr,
                                       edgeGhostsToReceive.toSliceConst(),
                                       nullptr,
                                       edgeGhostsToReceive.size(),
                                       edgeManager.m_localToGlobalMap.toSliceConst() );

  packedSize += bufferOps::Pack<true>( sendBufferPtr,
                                       faceGhostsToReceive.toSliceConst(),
                                       nullptr,
                                       faceGhostsToReceive.size(),
                                       faceManager.m_localToGlobalMap.toSliceConst() );



  for( localIndex er=0 ; er<elemManager.numRegions() ; ++er )
  {
    ElementRegionBase const * const elemRegion = elemManager.GetRegion( er );
    elemRegion->forElementSubRegionsIndex([&]( localIndex const esr, ElementSubRegionBase const * const subRegion )
    {
      packedSize+= bufferOps::Pack<true>( sendBufferPtr,
                                          elementGhostToReceive[er][esr].toSliceConst(),
                                          nullptr,
                                          elementGhostToReceive[er][esr].size(),
                                          subRegion->m_localToGlobalMap.toSliceConst() );
    });
  }

  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );

  this->MPI_iSendReceive( commID, MPI_COMM_GEOSX );
  this->MPI_WaitAll( commID );

  buffer_type const & receiveBuffer = ReceiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  int unpackedSize = 0;
  array1d<globalIndex> unmappedNodeIndices;
  array1d<globalIndex> unmappedEdgeIndices;
  array1d<globalIndex> unmappedFaceIndices;
  unpackedSize += bufferOps::Unpack( receiveBufferPtr,
                                     nodeGhostsToSend,
                                     unmappedNodeIndices,
                                     nodeManager.m_globalToLocalMap );
  GEOS_ERROR_IF( unmappedNodeIndices.size()!=0,
                 "Some node global indices were not mappable to a localIndex" );

  unpackedSize += bufferOps::Unpack( receiveBufferPtr,
                                     edgeGhostsToSend,
                                     unmappedEdgeIndices,
                                     edgeManager.m_globalToLocalMap );
  GEOS_ERROR_IF( unmappedEdgeIndices.size()!=0,
                 "Some edge global indices were not mappable to a localIndex" );

  unpackedSize += bufferOps::Unpack( receiveBufferPtr,
                                     faceGhostsToSend,
                                     unmappedFaceIndices,
                                     faceManager.m_globalToLocalMap );
  GEOS_ERROR_IF( unmappedFaceIndices.size()!=0,
                 "Some face global indices were not mappable to a localIndex" );


  nodeManager.SetGhostRankForSenders( nodeGhostsToSend );
  edgeManager.SetGhostRankForSenders( edgeGhostsToSend );
  faceManager.SetGhostRankForSenders( faceGhostsToSend );

  ElementRegionManager::ElementViewAccessor<ReferenceWrapper<localIndex_array>> elementGhostToSend =
    elemManager.ConstructReferenceAccessor<localIndex_array>( ObjectManagerBase:: viewKeyStruct::ghostsToSendString,
                                                              std::to_string( this->m_neighborRank ) );

  for( localIndex er=0 ; er<elemManager.numRegions() ; ++er )
  {
    ElementRegionBase * const elemRegion = elemManager.GetRegion( er );
    elemRegion->forElementSubRegionsIndex([&]( localIndex const esr, ElementSubRegionBase * const subRegion )
    {
      array1d<globalIndex> unmappedIndices;
      unpackedSize+= bufferOps::Unpack( receiveBufferPtr,
                                        elementGhostToSend[er][esr].get(),
                                        unmappedIndices,
                                        subRegion->m_globalToLocalMap );
      GEOS_ERROR_IF( unmappedIndices.size()!=0,
                     "Some element global indices were not mappable to a localIndex" );

      subRegion->SetGhostRankForSenders( elementGhostToSend[er][esr].get() );
    });
  }
}


int NeighborCommunicator::PackCommSizeForSync( std::map<string, string_array > const & fieldNames,
                                               MeshLevel * const mesh,
                                               int const commID )
{
  GEOSX_MARK_FUNCTION;

  NodeManager & nodeManager = *(mesh->getNodeManager());
  EdgeManager & edgeManager = *(mesh->getEdgeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());
  Group * const
  nodeNeighborData = nodeManager.GetGroup( nodeManager.groupKeys.neighborData )->
                     GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const
  edgeNeighborData = edgeManager.GetGroup( edgeManager.groupKeys.neighborData )->
                     GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const
  faceNeighborData = faceManager.GetGroup( faceManager.groupKeys.neighborData )->
                     GetGroup( std::to_string( this->m_neighborRank ) );

  localIndex_array const &
  nodeGhostsToSend = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );

  localIndex_array const &
  edgeGhostsToSend = edgeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );

  localIndex_array const &
  faceGhostsToSend = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToSend );

  ElementRegionManager::ElementViewAccessor<arrayView1d<localIndex>> const elementGhostToSend =
    elemManager.ConstructViewAccessor<array1d<localIndex>, arrayView1d<localIndex>>( ObjectManagerBase::
                                                         viewKeyStruct::
                                                         ghostsToSendString,
                                                         std::to_string( this->m_neighborRank ) );

  int bufferSize = 0;

  if( fieldNames.count( "node" ) > 0 )
  {
    bufferSize += nodeManager.PackSize( fieldNames.at( "node" ), nodeGhostsToSend, 0 );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    bufferSize += edgeManager.PackSize( fieldNames.at( "edge" ), edgeGhostsToSend, 0 );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    bufferSize += faceManager.PackSize( fieldNames.at( "face" ), faceGhostsToSend, 0 );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    for( localIndex er=0 ; er<elemManager.numRegions() ; ++er )
    {
      ElementRegionBase const * const elemRegion = elemManager.GetRegion( er );
      elemRegion->forElementSubRegionsIndex( [&]( localIndex const esr,
                                                  auto const * const subRegion )
      {
        bufferSize += subRegion->PackSize( fieldNames.at( "elems" ), elementGhostToSend[er][esr], 0 );
      });
    }
  }

  this->m_sendBufferSize[commID] = bufferSize;
  return bufferSize;
}


void NeighborCommunicator::PackCommBufferForSync( std::map<string, string_array > const & fieldNames,
                                                  MeshLevel * const mesh,
                                                  int const commID )
{
  GEOSX_MARK_FUNCTION;

  NodeManager & nodeManager = *(mesh->getNodeManager());
  EdgeManager & edgeManager = *(mesh->getEdgeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());
  Group * const
  nodeNeighborData = nodeManager.GetGroup( nodeManager.groupKeys.neighborData )->
                     GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const
  edgeNeighborData = edgeManager.GetGroup( edgeManager.groupKeys.neighborData )->
                     GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const
  faceNeighborData = faceManager.GetGroup( faceManager.groupKeys.neighborData )->
                     GetGroup( std::to_string( this->m_neighborRank ) );

  localIndex_array const &
  nodeGhostsToSend = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );

  localIndex_array const &
  edgeGhostsToSend = edgeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );

  localIndex_array const &
  faceGhostsToSend = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToSend );


  ElementRegionManager::ElementViewAccessor<arrayView1d<localIndex>> const elementGhostToSend =
    elemManager.ConstructViewAccessor<array1d<localIndex>, arrayView1d<localIndex>>( ObjectManagerBase::
                                                         viewKeyStruct::
                                                         ghostsToSendString,
                                                         std::to_string( this->m_neighborRank ) );

  buffer_type & sendBuffer = SendBuffer( commID );
  int const bufferSize =  integer_conversion<int>(sendBuffer.size());
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = 0;
  if( fieldNames.count( "node" ) > 0 )
  {
    packedSize += nodeManager.Pack( sendBufferPtr, fieldNames.at( "node" ), nodeGhostsToSend, 0 );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    packedSize += edgeManager.Pack( sendBufferPtr, fieldNames.at( "edge" ), edgeGhostsToSend, 0 );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    packedSize += faceManager.Pack( sendBufferPtr, fieldNames.at( "face" ), faceGhostsToSend, 0 );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    for( localIndex er=0 ; er<elemManager.numRegions() ; ++er )
    {
      ElementRegionBase const * const elemRegion = elemManager.GetRegion( er );
      elemRegion->forElementSubRegionsIndex([&]( localIndex const esr, auto const * const subRegion )
      {
        packedSize += subRegion->Pack( sendBufferPtr, fieldNames.at( "elems" ), elementGhostToSend[er][esr], 0 );
      });
    }
  }

  GEOS_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );
}


void NeighborCommunicator::SendRecvBuffers( int const commID )
{
  this->MPI_iSendReceive( commID, MPI_COMM_GEOSX );
}


void NeighborCommunicator::UnpackBufferForSync( std::map<string, string_array > const & fieldNames,
                                                MeshLevel * const mesh,
                                                int const commID )
{
  GEOSX_MARK_FUNCTION;

  buffer_type const & receiveBuffer = ReceiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  NodeManager & nodeManager = *(mesh->getNodeManager());
  EdgeManager & edgeManager = *(mesh->getEdgeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());

  Group * const nodeNeighborData = nodeManager.GetGroup( nodeManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const
  edgeNeighborData = edgeManager.GetGroup( edgeManager.groupKeys.neighborData )->
                     GetGroup( std::to_string( this->m_neighborRank ) );

  Group * const faceNeighborData = faceManager.
                                          GetGroup( faceManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  localIndex_array & nodeGhostsToReceive =
    nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToReceive );

  localIndex_array & edgeGhostsToReceive =
    edgeNeighborData->getReference<localIndex_array>( edgeManager.viewKeys.ghostsToReceive );

  localIndex_array & faceGhostsToReceive =
    faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToReceive );

  ElementRegionManager::ElementViewAccessor<arrayView1d<localIndex>> elementGhostToReceive =
    elemManager.ConstructViewAccessor<array1d<localIndex>, arrayView1d<localIndex>>( ObjectManagerBase::
                                                         viewKeyStruct::ghostsToReceiveString,
                                                         std::to_string( this->m_neighborRank ) );

  int unpackedSize = 0;

  if( fieldNames.count( "node" ) > 0 )
  {
    unpackedSize += nodeManager.Unpack( receiveBufferPtr, nodeGhostsToReceive, 0 );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    unpackedSize += edgeManager.Unpack( receiveBufferPtr, edgeGhostsToReceive, 0 );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    unpackedSize += faceManager.Unpack( receiveBufferPtr, faceGhostsToReceive, 0 );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    for( localIndex er=0 ; er<elemManager.numRegions() ; ++er )
    {
      ElementRegionBase * const elemRegion = elemManager.GetRegion( er );
      elemRegion->forElementSubRegionsIndex([&]( localIndex const esr,
                                                 auto * const subRegion )
      {
        unpackedSize += subRegion->Unpack( receiveBufferPtr, elementGhostToReceive[er][esr], 0 );
      });
    }
  }
}


} /* namespace geosx */
