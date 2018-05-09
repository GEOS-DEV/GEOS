/*
 * NeighborCommunicator.cpp
 *
 *  Created on: Jan 2, 2018
 *      Author: settgast
 */

#include "NeighborCommunicator.hpp"
#include "mesh/MeshLevel.hpp"

namespace geosx
{

NeighborCommunicator::NeighborCommunicator():
  m_neighborRank(-1),
  m_sendBufferSize(),
  m_receiveBufferSize(),
  m_sendBuffer(),
  m_receiveBuffer(),
  m_mpiSendBufferRequest(),
  m_mpiRecvBufferRequest(),
  m_mpiSendBufferStatus(),
  m_mpiRecvBufferStatus()
{
}

//NeighborCommunicator::~NeighborCommunicator()
//{
//}

void NeighborCommunicator::MPI_iSendReceive( char const * const sendBuffer,
                                             int const sendSize,
                                             MPI_Request& sendRequest,
                                             char * const receiveBuffer,
                                             int const receiveSize,
                                             MPI_Request& receiveRequest,
                                             int const commID,
                                             MPI_Comm mpiComm )
{
  int const sendTag = CommTag( Rank(), m_neighborRank, commID );
  //m_rank * m_size + m_neighborRank + m_size*m_size*commID;
  MPI_Isend( const_cast<char*>(sendBuffer),
             sendSize,
             MPI_CHAR,
             m_neighborRank,
             sendTag,
             mpiComm,
             &sendRequest );

  int const receiveTag = CommTag( m_neighborRank, Rank(), commID );
  //m_neighborRank * m_size + m_rank + m_size*m_size*commID;
  MPI_Irecv( receiveBuffer,
             receiveSize,
             MPI_CHAR,
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
  m_sendBufferSize[commID] = integer_conversion<int>(m_sendBuffer[commID].size());
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
                    integer_conversion<int>(m_sendBuffer[commID].size()),
                    mpiSendRequest,
                    m_receiveBuffer[commID].data(),
                    integer_conversion<int>(m_receiveBuffer[commID].size()),
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

  MPI_Waitall( 1, &( m_mpiRecvBufferRequest[commID] ), &( m_mpiRecvBufferStatus[commID] ) );
  MPI_Waitall( 1, &( m_mpiSendBufferRequest[commID] ), &( m_mpiSendBufferStatus[commID] ) );

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

void NeighborCommunicator::MPI_iSendReceive( char const * const sendBuffer,
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

  MPI_Waitall( 1, &( m_mpiRecvBufferRequest[commID] ), &( m_mpiRecvBufferStatus[commID] ) );
  MPI_Waitall( 1, &( m_mpiSendBufferRequest[commID] ), &( m_mpiSendBufferStatus[commID] ) );

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

void NeighborCommunicator::MPI_WaitAll( int const commID,
                                        MPI_Request& mpiSendRequest,
                                        MPI_Status& mpiSendStatus,
                                        MPI_Request& mpiRecvRequest,
                                        MPI_Status& mpiReceiveStatus )

{
  MPI_Waitall( 1, &mpiRecvRequest, &mpiReceiveStatus );
  MPI_Waitall( 1, &mpiSendRequest, &mpiSendStatus );
}

void NeighborCommunicator::MPI_WaitAll( int const commID )
{
  MPI_Waitall( 1, &( m_mpiRecvBufferRequest[commID] ), &( m_mpiRecvBufferStatus[commID] ) );
  MPI_Waitall( 1, &( m_mpiSendBufferRequest[commID] ), &( m_mpiSendBufferStatus[commID] ) );
}

int NeighborCommunicator::Rank()
{
  int rank = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  return rank;
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
  ManagedGroup * neighborGroups[100];
  localIndex numNeighborGroups = 0;

  ObjectManagerBase * const nodeManager = mesh->getNodeManager();
  neighborGroups[numNeighborGroups++] = nodeManager->
      GetGroup(nodeManager->m_ObjectManagerBaseGroupKeys.neighborData)->
      RegisterGroup( std::to_string( this->m_neighborRank ));

  ObjectManagerBase * const faceManager = mesh->getFaceManager();
  neighborGroups[numNeighborGroups++] = faceManager->
      GetGroup(faceManager->m_ObjectManagerBaseGroupKeys.neighborData)->
      RegisterGroup( std::to_string( this->m_neighborRank ));

  ElementRegionManager * const elemManager = mesh->getElemManager();
  elemManager->forCellBlocks( [&]( ManagedGroup * const cellBlock ) -> void
  {
    neighborGroups[numNeighborGroups++] = cellBlock->
        GetGroup(faceManager->m_ObjectManagerBaseGroupKeys.neighborData)->
        RegisterGroup( std::to_string( this->m_neighborRank ));
  });


  for( localIndex a=0 ; a<numNeighborGroups ; ++a )
  {
    neighborGroups[a]->
    RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::matchedPartitionBoundaryObjectsString)->
    setSizedFromParent(0);

    neighborGroups[a]->
    RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::ghostsToSendString)->
    setSizedFromParent(0);

    neighborGroups[a]->
    RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::ghostsToReceiveString)->
    setSizedFromParent(0);

    neighborGroups[a]->
    RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::adjacencyListString)->
    setSizedFromParent(0);

  }



}

void NeighborCommunicator::FindAndPackGhosts( bool const contactActive,
                                       integer const depth,
                                       MeshLevel * const mesh,
                                       int const commID )
{
  NodeManager & nodeManager = *(mesh->getNodeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());

  ManagedGroup * const nodeNeighborData = nodeManager.
                                          GetGroup(nodeManager.groupKeys.neighborData)->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  ManagedGroup * const faceNeighborData = faceManager.
                                          GetGroup( faceManager.groupKeys.neighborData )->
                                          GetGroup( std::to_string( this->m_neighborRank ) );

  localIndex_array & nodeAdjacencyList = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.adjacencyList );
  localIndex_array  edgeAdjacencyList ;//= edgeNeighborData->getReference<localIndex_array>( edgeManager.viewKeys.ghostsToSend );;
  localIndex_array & faceAdjacencyList = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.adjacencyList );
  ElementRegionManager::ElementViewAccessor<localIndex_array> elementAdjacencyList =
  elemManager.ConstructViewAccessor<localIndex_array>( ObjectManagerBase::viewKeyStruct::adjacencyListString,
                                                       std::to_string( this->m_neighborRank ) );

  mesh->GenerateAdjacencyLists( nodeNeighborData->getReference<localIndex_array>(nodeManager.viewKeys.matchedPartitionBoundaryObjects),
                                nodeAdjacencyList,
                                edgeAdjacencyList,
                                faceAdjacencyList,
                                elementAdjacencyList,
                                2 );

  int bufferSize = 0;

  bufferSize += nodeManager.PackGlobalMapsSize( nodeAdjacencyList, 0 );
  bufferSize += faceManager.PackGlobalMapsSize( faceAdjacencyList, 0 );
  bufferSize += elemManager.PackGlobalMapsSize( elementAdjacencyList );

  bufferSize += nodeManager.PackUpDownMapsSize( nodeAdjacencyList );
  bufferSize += faceManager.PackUpDownMapsSize( faceAdjacencyList );
  bufferSize += elemManager.PackUpDownMapsSize( elementAdjacencyList );

  bufferSize += nodeManager.PackSize( {}, nodeAdjacencyList, 0 );
  bufferSize += faceManager.PackSize( {}, faceAdjacencyList, 0 );
  bufferSize += elemManager.PackSize( {}, elementAdjacencyList );




  buffer_type & sendBuffer = SendBuffer(commID);
  sendBuffer.resize(bufferSize);

  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = 0;

  packedSize += nodeManager.PackGlobalMaps( sendBufferPtr, nodeAdjacencyList, 0 );
  packedSize += faceManager.PackGlobalMaps( sendBufferPtr, faceAdjacencyList, 0 );
  packedSize += elemManager.PackGlobalMaps( sendBufferPtr, elementAdjacencyList );

  packedSize += nodeManager.PackUpDownMaps( sendBufferPtr, nodeAdjacencyList );
  packedSize += faceManager.PackUpDownMaps( sendBufferPtr, faceAdjacencyList );
  packedSize += elemManager.PackUpDownMaps( sendBufferPtr, elementAdjacencyList );

  packedSize += nodeManager.Pack( sendBufferPtr, {} ,nodeAdjacencyList, 0 );
  packedSize += faceManager.Pack( sendBufferPtr, {}, faceAdjacencyList, 0 );
  packedSize += elemManager.Pack( sendBufferPtr, {}, elementAdjacencyList );


  GEOS_ASSERT( bufferSize == packedSize, "Allocated Buffer Size is not equal to packed buffer size")

  this->MPI_iSendReceive( commID, MPI_COMM_WORLD );

}


void NeighborCommunicator::UnpackGhosts( MeshLevel * const mesh,
                                         int const commID )
{
  NodeManager & nodeManager = *(mesh->getNodeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());

//  ManagedGroup * const nodeNeighborData = nodeManager.
//                                          GetGroup(nodeManager.groupKeys.neighborData)->
//                                          GetGroup( std::to_string( this->m_neighborRank ) );
//
//  ManagedGroup * const faceNeighborData = faceManager.
//                                          GetGroup( faceManager.groupKeys.neighborData )->
//                                          GetGroup( std::to_string( this->m_neighborRank ) );


  buffer_type const & receiveBuffer = ReceiveBuffer(commID);
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();


  int unpackedSize = 0;

//  localIndex_array & nodeGhostsToSend = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );
//  localIndex_array & faceGhostsToSend = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToSend );
//
//  localIndex_array & nodeGhostsToReceive = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToReceive );
//  localIndex_array & faceGhostsToReceive = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToReceive );
  ElementRegionManager::ElementViewAccessor<localIndex_array> elementAdjacencyReceiveList =
      elemManager.ConstructViewAccessor<localIndex_array>( ObjectManagerBase::viewKeyStruct::ghostsToReceiveString,
                                                           std::to_string( this->m_neighborRank ) );


  localIndex_array nodeUpackList;
  unpackedSize += nodeManager.UnpackGlobalMaps( receiveBufferPtr, nodeUpackList, 0);

  localIndex_array faceUnpackList;
  unpackedSize += faceManager.UnpackGlobalMaps( receiveBufferPtr, faceUnpackList, 0);

  unpackedSize += elemManager.UnpackGlobalMaps( receiveBufferPtr, elementAdjacencyReceiveList);


  unpackedSize += nodeManager.UnpackUpDownMaps( receiveBufferPtr, nodeUpackList);
  unpackedSize += faceManager.UnpackUpDownMaps( receiveBufferPtr, faceUnpackList);
  unpackedSize += elemManager.UnpackUpDownMaps( receiveBufferPtr, elementAdjacencyReceiveList);


  unpackedSize += nodeManager.Unpack( receiveBufferPtr, nodeUpackList, 0);
  unpackedSize += faceManager.Unpack( receiveBufferPtr, faceUnpackList, 0);
  unpackedSize += elemManager.Unpack( receiveBufferPtr, elementAdjacencyReceiveList);

  std::cout<<"herro"<<std::endl;


}

void NeighborCommunicator::RebuildSyncLists( MeshLevel * const mesh,
                                             int const commID )
{
  NodeManager & nodeManager = *(mesh->getNodeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());
  ManagedGroup * const nodeNeighborData = nodeManager.GetGroup(nodeManager.groupKeys.neighborData)->
                                                      GetGroup( std::to_string( this->m_neighborRank ) );

  ManagedGroup * const faceNeighborData = faceManager.
                                            GetGroup( faceManager.groupKeys.neighborData )->
                                            GetGroup( std::to_string( this->m_neighborRank ) );


  localIndex_array & nodeGhostsToSend = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );
  localIndex_array & faceGhostsToSend = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToSend );

  localIndex_array const & nodeGhostsToReceive = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToReceive );
  localIndex_array const & faceGhostsToReceive = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToReceive );


  buffer_type & sendBuffer = SendBuffer(commID);
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int bufferSize = 0;
  bufferSize += CommBufferOps::Pack<false>( sendBufferPtr,
                                            nodeGhostsToReceive.data(),
                                            nodeGhostsToReceive.size(),
                                            nodeManager.m_localToGlobalMap );

  bufferSize += CommBufferOps::Pack<false>( sendBufferPtr,
                                            faceGhostsToReceive.data(),
                                            faceGhostsToReceive.size(),
                                            faceManager.m_localToGlobalMap );

  sendBuffer.resize(bufferSize);

  int packedSize = 0;
  packedSize += CommBufferOps::Pack<true>( sendBufferPtr,
                                            nodeGhostsToReceive.data(),
                                            nodeGhostsToReceive.size(),
                                            nodeManager.m_localToGlobalMap );

  packedSize += CommBufferOps::Pack<true>( sendBufferPtr,
                                            faceGhostsToReceive.data(),
                                            faceGhostsToReceive.size(),
                                            faceManager.m_localToGlobalMap );

  GEOS_ASSERT( bufferSize == packedSize, "Allocated Buffer Size is not equal to packed buffer size")

  this->MPI_iSendReceive( commID, MPI_COMM_WORLD );
  this->MPI_WaitAll(commID);

  buffer_type const & receiveBuffer = ReceiveBuffer(commID);
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  int unpackedSize = 0;
  unpackedSize += CommBufferOps::Unpack( receiveBufferPtr,
                                         nodeGhostsToSend,
                                         nodeManager.m_globalToLocalMap );
  unpackedSize += CommBufferOps::Unpack( receiveBufferPtr,
                                         faceGhostsToSend,
                                         faceManager.m_globalToLocalMap );


}


void NeighborCommunicator::PackBufferForSync( std::map<string, array<string> > const & fieldNames,
                                              MeshLevel * const mesh,
                                              int const commID )
{
  NodeManager & nodeManager = *(mesh->getNodeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());
  ManagedGroup * const nodeNeighborData = nodeManager.GetGroup(nodeManager.groupKeys.neighborData)->
                                                      GetGroup( std::to_string( this->m_neighborRank ) );

  ManagedGroup * const faceNeighborData = faceManager.
                                            GetGroup( faceManager.groupKeys.neighborData )->
                                            GetGroup( std::to_string( this->m_neighborRank ) );


  localIndex_array const & nodeGhostsToSend = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );
  localIndex_array const & faceGhostsToSend = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToSend );



  int bufferSize = 0;

  bufferSize += nodeManager.PackSize( fieldNames.at("node"), nodeGhostsToSend, 0 );
//  bufferSize += faceManager.PackSize( {}, nodeGhostsToSend, 1, 0 );
//  bufferSize += elemManager.PackSize( {}, elementAdjacencyList );
//
//
  buffer_type & sendBuffer = SendBuffer(commID);
  sendBuffer.resize(bufferSize);

  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = 0;
  packedSize += nodeManager.Pack( sendBufferPtr, fieldNames.at("node"), nodeGhostsToSend, 0 );
//  packedSize += faceManager.Pack( sendBufferPtr, {}, faceAdjacencyList, 1, 0 );
//  packedSize += elemManager.Pack( sendBufferPtr, {}, elementAdjacencyList );

  GEOS_ASSERT( bufferSize == packedSize, "Allocated Buffer Size is not equal to packed buffer size")
  this->MPI_iSendReceive( commID, MPI_COMM_WORLD );


//  MPI_WaitAll( commID );
}

void NeighborCommunicator::UnpackBufferForSync( MeshLevel * const mesh,
                                                int const commID )
{

  buffer_type const & receiveBuffer = ReceiveBuffer(commID);
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  NodeManager & nodeManager = *(mesh->getNodeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elemManager = *(mesh->getElemManager());

  ManagedGroup * const nodeNeighborData = nodeManager.GetGroup(nodeManager.groupKeys.neighborData)->
                                                      GetGroup( std::to_string( this->m_neighborRank ) );

  ManagedGroup * const faceNeighborData = faceManager.
                                            GetGroup( faceManager.groupKeys.neighborData )->
                                            GetGroup( std::to_string( this->m_neighborRank ) );

  localIndex_array & nodeGhostsToReceive = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToReceive );
  localIndex_array & faceGhostsToReceive = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToReceive );

  //
//
//  ElementRegionManager::ElementViewAccessor<localIndex_array> elementAdjacencyReceiveList;
//
//  elemManager.ConstructViewAccessor( ObjectManagerBase::viewKeyStruct::ghostsToReceiveString,
//                                     std::to_string( this->m_neighborRank ),
//                                     elementAdjacencyReceiveList );
  int unpackedSize = 0;

  unpackedSize += nodeManager.Unpack( receiveBufferPtr, nodeGhostsToReceive, 0);
//  unpackedSize += faceManager.Unpack( receiveBufferPtr, faceGhostsToReceive, 0);
//  unpackedSize += elemManager.Unpack( receiveBufferPtr, elementAdjacencyReceiveList);
//
//  std::cout<<"herro"<<std::endl;
//

}


} /* namespace geosx */
