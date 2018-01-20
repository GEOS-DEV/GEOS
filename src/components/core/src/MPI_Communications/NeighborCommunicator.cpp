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

NeighborCommunicator::NeighborCommunicator()
{
}

NeighborCommunicator::~NeighborCommunicator()
{
}

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
  MPI_Isend( sendBuffer,
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
  m_sendBufferSize[commID] = m_sendBuffer[commID].size();
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
                    m_sendBuffer[commID].size(),
                    mpiSendRequest,
                    m_receiveBuffer[commID].data(),
                    m_receiveBuffer[commID].size(),
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
      GetGroup(nodeManager->groupKeys.neighborData)->
      RegisterGroup( std::to_string( this->m_neighborRank ));

  ObjectManagerBase * const faceManager = mesh->getFaceManager();
  neighborGroups[numNeighborGroups++] = faceManager->
      GetGroup(faceManager->groupKeys.neighborData)->
      RegisterGroup( std::to_string( this->m_neighborRank ));

  ElementRegionManager * const elemManager = mesh->getElemManager();
  elemManager->forCellBlocks( [&]( ManagedGroup * const cellBlock ) -> void
  {
    neighborGroups[numNeighborGroups++] = cellBlock->
        GetGroup(faceManager->groupKeys.neighborData)->
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
  }



}

void NeighborCommunicator::FindGhosts( bool const contactActive,
                                       integer const depth,
                                       MeshLevel * const mesh )
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

  localIndex_array & nodeAdjacencyList = nodeNeighborData->getReference<localIndex_array>( nodeManager.viewKeys.ghostsToSend );
  localIndex_array  edgeAdjacencyList ;//= edgeNeighborData->getReference<localIndex_array>( edgeManager.viewKeys.ghostsToSend );;
  localIndex_array & faceAdjacencyList = faceNeighborData->getReference<localIndex_array>( faceManager.viewKeys.ghostsToSend );
  ElementRegionManager::ElementViewAccessor<localIndex_array> elementAdjacencyList;
  elemManager.ConstructViewAccessor( ObjectManagerBase::viewKeyStruct::ghostsToSendString,
                                     std::to_string( this->m_neighborRank ),
                                     elementAdjacencyList );

  mesh->GenerateAdjacencyLists( nodeNeighborData->getReference<localIndex_array>(nodeManager.viewKeys.matchedPartitionBoundaryObjects),
                                nodeAdjacencyList,
                                edgeAdjacencyList,
                                faceAdjacencyList,
                                elementAdjacencyList,
                                2 );


//
//  lSet & facesToSend;
//
//
//
//  //Currently, this MUST fill the objectsToSend associated with the following
//  //(see NeighborCommunication::SyncNames for the current list)
//  // 0: NodeManager -> PackNodes
//  // 1: EdgeManager -> PackEdges
//  // 2: FaceManager -> PackFaces
//  // 3: FiniteElementElementManager -> PackElements
//  // 4: DiscreteElementManager ->
//  // 5: DiscreteElementNodeManager ->
//  // 6: DiscreteElementFaceManager ->
//  // 7: EllipsoidalDiscreteElementManager ->
//
//  // so now we know which nodes are "shared" between the neighbors. what we do now is to collect all the objects
//  // attached to these nodes. These will be "ghost" objects on the neighbor.
//
//  lSet& allNodes = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementNodeManager];
//  lSet& allEdges = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementEdgeManager];
//  lSet& allFaces = tempNeighborData.objectLocalIndicesToSend[PhysicalDomainT::FiniteElementFaceManager];
//
//  // get list of elements connected to the matched nodes into "m_elementRegionsSendLocalIndices"
//  this->m_elementRegionsSendLocalIndices.clear();
//  m_domain->m_feElementManager.ConstructListOfIndexesFromMap( this->m_domain->m_feNodeManager.m_toElementsRelation,
//                                                              this->tempNeighborData.matchedIndices[PhysicalDomainT::FiniteElementNodeManager],
//                                                              this->m_elementRegionsSendLocalIndices,
//                                                              depth );
//
//  // now pack up the elements that are going over to the neighbor...also get the node that we are going to
//  // need to send.
//  m_domain
//          ->m_feElementManager.PackElements(
//      this->tempNeighborData.objectsToSend[PhysicalDomainT::FiniteElementElementManager],
//      allNodes,
//      allFaces,
//      this->m_elementRegionsSendLocalIndices,
//      this->m_domain->m_feNodeManager,
//      this->m_domain->m_feFaceManager,
//      true, true, true, true );
//
//  // add the "external" faces if the contact is on
//  if( contactActive )
//  {
//    const array<integer>& isExternalFace =
//        m_domain->m_feFaceManager.m_isExternal;
//
//    for( array<integer>::size_type a = 0 ;
//        a < m_domain->m_feFaceManager.m_numFaces ; ++a )
//    {
//      if( isExternalFace[a] == 1 )
//      {
//        bool allValidNodes = true;
//        for( localIndex_array::const_iterator
//        i = m_domain->m_feFaceManager.m_toNodesRelation[a].begin() ;
//            i != m_domain->m_feFaceManager.m_toNodesRelation[a].end() ; ++i )
//        {
//          const globalIndex gnode = m_domain->m_feNodeManager.m_localToGlobalMap[*i];
//          const int owningRank = GlobalIndexManager::OwningRank( gnode );
//          if( !( m_rankOfNeighborNeighbors.count( owningRank ) ) )
//          {
//            allValidNodes = false;
//          }
//        }
//        if( allValidNodes )
//        {
//          allFaces.insert( a );
//        }
//      }
//    }
//  }
//
//  for( lSet::const_iterator faceIndex = allFaces.begin() ;
//      faceIndex != allFaces.end() ; ++faceIndex )
//  {
//    for( localIndex_array::const_iterator
//    edgeIndex = m_domain->m_feFaceManager.m_toEdgesRelation[*faceIndex].begin() ;
//        edgeIndex != m_domain->m_feFaceManager.m_toEdgesRelation[*faceIndex].end()
//        ; ++edgeIndex )
//    {
//      const globalIndex gi =
//          m_domain->m_feEdgeManager.m_localToGlobalMap[*edgeIndex];
//      if( m_rankOfNeighborNeighbors.count( GlobalIndexManager::OwningRank( gi )
//                                                                           ) )
//      {
//        allEdges.insert( *edgeIndex );
//      }
//    }
//
//    for( localIndex_array::const_iterator
//    i = m_domain->m_feFaceManager.m_toNodesRelation[*faceIndex].begin() ;
//        i != m_domain->m_feFaceManager.m_toNodesRelation[*faceIndex].end() ;
//        ++i )
//    {
//      allNodes.insert( *i );
//    }
//  }
//
//  const PhysicalDomainT::ObjectDataStructureKeys keys[3] =
//      {
//        PhysicalDomainT::FiniteElementNodeManager,
//        PhysicalDomainT::FiniteElementEdgeManager,
//        PhysicalDomainT::FiniteElementFaceManager
//      };
//
//  for( int i = 0 ; i < 3 ; ++i )
//  {
//    const lSet& localSends =
//        tempNeighborData.objectLocalIndicesToSend[keys[i]];
//    globalIndex_array& globalSends =
//        tempNeighborData.objectGlobalIndicesToSend[keys[i]];
//
//    for( lSet::const_iterator a = localSends.begin() ; a != localSends.end() ; ++a
//        )
//    {
//      const globalIndex gIndex =
//          m_domain->GetObjectDataStructure( keys[i] ).m_localToGlobalMap[*a];
//      if( GlobalIndexManager::OwningRank( gIndex ) != m_rank &&
//          GlobalIndexManager::OwningRank( gIndex ) != m_neighborRank )
//      {
//        globalSends.push_back( gIndex );
//      }
//    }
//  }
}

} /* namespace geosx */
