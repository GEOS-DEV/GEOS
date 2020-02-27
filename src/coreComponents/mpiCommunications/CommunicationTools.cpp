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
 * @file CommunicationTools.cpp
 *
 */

#include "mpiCommunications/CommunicationTools.hpp"


#include "common/TimingMacros.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/ObjectManagerBase.hpp"

#include <algorithm>

namespace geosx
{

using namespace dataRepository;

CommunicationTools::CommunicationTools()
{
  // TODO Auto-generated constructor stub
}

CommunicationTools::~CommunicationTools()
{
  // TODO Auto-generated destructor stub
}


std::set< int > & CommunicationTools::getFreeCommIDs()
{
  static std::set< int > commIDs;
  static bool isInitialized = false;

  if( !isInitialized )
  {
    for( int a = 0 ; a < NeighborCommunicator::maxComm ; ++a )
    {
      commIDs.insert( a );
    }
    isInitialized = true;
  }

  return commIDs;
}


int CommunicationTools::reserveCommID()
{
  std::set< int > & commIDs = getFreeCommIDs();

  int rval = *( commIDs.begin() );
  commIDs.erase( rval );
  return rval;
}

void CommunicationTools::releaseCommID( int & ID )
{
  std::set< int > & commIDs = getFreeCommIDs();

  if( commIDs.count( ID ) > 0 )
  {
    GEOSX_ERROR( "Attempting to release commID that is already free" );
  }
  commIDs.insert( ID );
  ID = -1;
}

void CommunicationTools::AssignGlobalIndices( ObjectManagerBase & object,
                                              ObjectManagerBase const & compositionObject,
                                              array1d< NeighborCommunicator > & neighbors )
{
  GEOSX_MARK_FUNCTION;
  integer_array & ghostRank = object.getReference< integer_array >( object.m_ObjectManagerBaseViewKeys.ghostRank );
  ghostRank = -2;

  int const commSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  localIndex numberOfObjectsHere = object.size();
  localIndex_array numberOfObjects( commSize );
  localIndex_array glocalIndexOffset( commSize );
  MpiWrapper::allGather( numberOfObjectsHere, numberOfObjects );

  int const commRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );

  glocalIndexOffset[0] = 0;
  for( int rank = 1 ; rank < commSize ; ++rank )
  {
    glocalIndexOffset[rank] = glocalIndexOffset[rank - 1] + numberOfObjects[rank - 1];
  }

  // set the global indices as if they were all local to this process
  for( localIndex a = 0 ; a < object.size() ; ++a )
  {
    object.m_localToGlobalMap[a] = glocalIndexOffset[commRank] + a;
  }

  // get the relation to the composition object used that will be used to identify the main object. For example,
  // a face can be identified by its nodes.
  std::vector< std::vector< globalIndex > > objectToCompositionObject;
  object.ExtractMapFromObjectForAssignGlobalIndexNumbers( &compositionObject, objectToCompositionObject );

  // now arrange the data from objectToCompositionObject into a map "indexByFirstCompositionIndex", such that the key
  // is the lowest global index of the composition object that make up this object. The value of the map is a pair, with
  // the array being the remaining composition object global indices, and the second being the global index of the
  // object
  // itself.
  map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, localIndex > > > indexByFirstCompositionIndex;

  localIndex bufferSize = 0;
  for( std::size_t a = 0 ; a < objectToCompositionObject.size() ; ++a )
  {
    if( objectToCompositionObject[a].size() > 0 )
    {
      // set nodelist array
      std::vector< globalIndex > const & nodeList = objectToCompositionObject[a];

      // grab the first global index of the composition objects
      const globalIndex firstCompositionIndex = nodeList[0];

      // create a temporary to hold the pair
      std::pair< std::vector< globalIndex >, globalIndex > tempComp;

      // fill the array with the remaining composition object global indices
      tempComp.first.insert( tempComp.first.begin(), nodeList.begin() + 1, nodeList.end() );

      // set the second value of the pair to the localIndex of the object.
      tempComp.second = a;

      // push the tempComp onto the map.
      indexByFirstCompositionIndex[firstCompositionIndex].push_back( std::move( tempComp ) );
      bufferSize += 2 + nodeList.size();
    }
  }

  globalIndex_array objectToCompositionObjectSendBuffer;
  objectToCompositionObjectSendBuffer.reserve( bufferSize );

  // put the map into a buffer
  for( std::size_t a = 0 ; a < objectToCompositionObject.size() ; ++a )
  {
    if( objectToCompositionObject[a].size() > 0 )
    {
      std::vector< globalIndex > const & nodeList = objectToCompositionObject[a];
      objectToCompositionObjectSendBuffer.push_back( nodeList.size() );
      objectToCompositionObjectSendBuffer.push_back( object.m_localToGlobalMap[a] );
      for( std::size_t b = 0 ; b < nodeList.size() ; ++b )
      {
        objectToCompositionObjectSendBuffer.push_back( nodeList[b] );
      }
    }
  }


  MPI_iCommData commData;
  commData.resize( neighbors.size());
//  int commID = reserveCommID();
//
//  // send the composition buffers
//  {
//    int const sendSize = integer_conversion<int const>( objectToCompositionObjectSendBuffer.size() *
// sizeof(globalIndex));
//
//    for( localIndex in = 0 ; in < neighbors.size() ; ++in )
//    {
//      NeighborCommunicator & neighbor = neighbors[in];
//
//      neighbor.MPI_iSendReceive( reinterpret_cast<const char*>( objectToCompositionObjectSendBuffer.data() ),
//                                 sendSize,
//                                 commID,
//                                 MPI_COMM_GEOSX );
//    }
//    for( localIndex in = 0 ; in < neighbors.size() ; ++in )
//    {
//      neighbors[in].MPI_WaitAll( commID );
//    }
//  }

  array1d< int >  receiveBufferSizes( neighbors.size());
  array1d< globalIndex_array > receiveBuffers( neighbors.size());

  int const sendSize = integer_conversion< int const >( objectToCompositionObjectSendBuffer.size() );

  for( localIndex neighborIndex = 0 ; neighborIndex < neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    neighbor.MPI_iSendReceive( &sendSize,
                               1,
                               commData.mpiSizeSendBufferRequest[neighborIndex],
                               &(receiveBufferSizes[neighborIndex]),
                               1,
                               commData.mpiSizeRecvBufferRequest[neighborIndex],
                               commData.sizeCommID,
                               MPI_COMM_GEOSX );
  }


  for( localIndex count=0 ; count<neighbors.size() ; ++count )
  {
    int neighborIndex;
    MpiWrapper::Waitany( commData.size,
                         commData.mpiSizeRecvBufferRequest.data(),
                         &neighborIndex,
                         commData.mpiSizeRecvBufferStatus.data() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    receiveBuffers[neighborIndex].resize( receiveBufferSizes[neighborIndex] );
    neighbor.MPI_iSendReceive( objectToCompositionObjectSendBuffer.data(),
                               sendSize,
                               commData.mpiSendBufferRequest[neighborIndex],
                               receiveBuffers[neighborIndex].data(),
                               receiveBufferSizes[neighborIndex],
                               commData.mpiRecvBufferRequest[neighborIndex],
                               commData.commID,
                               MPI_COMM_GEOSX );

  }

  // unpack the data from neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager] to
  // the local arrays

  // object to receive the neighbor data
  // this baby is an Array (for each neighbor) of maps, with the key of lowest composition index, and a value
  // containing an array containing the std::pairs of the remaining composition indices, and the globalIndex of the
  // object.
  std::vector< map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, globalIndex > > > >
  neighborCompositionObjects( neighbors.size() );

  for( localIndex count=0 ; count<neighbors.size() ; ++count )
  {
    int neighborIndex;
    MpiWrapper::Waitany( commData.size,
                         commData.mpiRecvBufferRequest.data(),
                         &neighborIndex,
                         commData.mpiRecvBufferStatus.data() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    globalIndex const * recBuffer = receiveBuffers[neighborIndex].data();
    localIndex recBufferSize = receiveBufferSizes[neighborIndex];
    globalIndex const * endBuffer = recBuffer + recBufferSize;
    // iterate over data that was just received
    while( recBuffer < endBuffer )
    {
      // the first thing packed was the data size for a given object
      localIndex dataSize = integer_conversion< localIndex >( *recBuffer++ );

      // the second thing packed was the globalIndex of that object
      const globalIndex neighborGlobalIndex = *( recBuffer++ );

      // the global indices of the composition objects were next. they are ordered, so the lowest one is first.
      const globalIndex firstCompositionIndex = *( recBuffer++ );

      // the remaining composition object indices.
      std::vector< globalIndex > temp;
      for( localIndex b = 1 ; b < dataSize ; ++b )
      {
        temp.push_back( *( recBuffer++ ) );
      }

      // fill neighborCompositionObjects
      std::pair< std::vector< globalIndex >, globalIndex >
      tempComp( std::make_pair( std::move( temp ), std::move( neighborGlobalIndex ) ) );
//        tempComp( std::make_pair( temp, neighborGlobalIndex ) );

      neighborCompositionObjects[neighborIndex][firstCompositionIndex].push_back( tempComp );
    }
//    }
//
//
//  // now check to see if the global index is valid. We do this by checking the contents of neighborCompositionObjects
//  // with indexByFirstCompositionIndex
//  for( localIndex neighborIndex = 0 ; neighborIndex < neighbors.size() ; ++neighborIndex )
//  {
//    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    // Set iterators to the beginning of each indexByFirstCompositionIndex,
    // and neighborCompositionObjects[neighborNum].
    map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, localIndex > > >::const_iterator
      iter_local = indexByFirstCompositionIndex.begin();
    map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, globalIndex > > >::const_iterator
      iter_neighbor = neighborCompositionObjects[neighborIndex].begin();

    // now we continue the while loop as long as both of our iterators are in range.
    while( iter_local != indexByFirstCompositionIndex.end() &&
           iter_neighbor != neighborCompositionObjects[neighborIndex].end() )
    {
      // check to see if the map keys (first composition index) are the same.
      if( iter_local->first == iter_neighbor->first )
      {
        // first we loop over all local composition arrays (objects with the matched key)
        for( std::vector< std::pair< std::vector< globalIndex >, localIndex > >::const_iterator
             iter_local2 = iter_local->second.begin() ;
             iter_local2 != iter_local->second.end() ; ++iter_local2 )
        {
          // and loop over all of the neighbor composition arrays (objects with the matched key)
          for( std::vector< std::pair< std::vector< globalIndex >, globalIndex > >::const_iterator
               iter_neighbor2 = iter_neighbor->second.begin() ;
               iter_neighbor2 != iter_neighbor->second.end() ;
               ++iter_neighbor2 )
          {
            // now compare the composition arrays
            if( iter_local2->first.size() == iter_neighbor2->first.size() &&
                std::equal( iter_local2->first.begin(),
                            iter_local2->first.end(),
                            iter_neighbor2->first.begin() ) )
            {
              // they are equal, so we need to overwrite the global index for the object
              if( iter_neighbor2->second < object.m_localToGlobalMap[iter_local2->second] )
              {
                if( neighbor.NeighborRank() < commRank )
                {
                  object.m_localToGlobalMap[iter_local2->second] = iter_neighbor2->second;
                  ghostRank[iter_local2->second] = neighbor.NeighborRank();
                }
                else
                {
                  ghostRank[iter_local2->second] = -1;
                }

              }

              // we should break out of the iter_local2 loop since we aren't going to find another match.
              break;
            }
          }
        }
        ++iter_local;
        ++iter_neighbor;
      }
      else if( iter_local->first < iter_neighbor->first )
      {
        ++iter_local;
      }
      else if( iter_local->first > iter_neighbor->first )
      {
        ++iter_neighbor;
      }
    }
  }

  object.ConstructGlobalToLocalMap();

  globalIndex maxGlobalIndex = -1;
  for( localIndex a=0 ; a<object.m_localToGlobalMap.size() ; ++a )
  {
    maxGlobalIndex = std::max( maxGlobalIndex, object.m_localToGlobalMap[a] );
  }

  MpiWrapper::allReduce( &maxGlobalIndex,
                         &(object.m_maxGlobalIndex),
                         1,
                         MPI_MAX,
                         MPI_COMM_GEOSX );

}

void CommunicationTools::AssignNewGlobalIndices( ObjectManagerBase & object,
                                                 std::set< localIndex > const & indexList )
{
  int const thisRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const commSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  localIndex numberOfNewObjectsHere = indexList.size();
  localIndex_array numberOfNewObjects( commSize );
  localIndex_array glocalIndexOffset( commSize );
  MpiWrapper::allGather( numberOfNewObjectsHere, numberOfNewObjects );


  glocalIndexOffset[0] = 0;
  for( int rank = 1 ; rank < commSize ; ++rank )
  {
    glocalIndexOffset[rank] = glocalIndexOffset[rank - 1] + numberOfNewObjects[rank - 1];
  }

  localIndex nIndicesAssigned = 0;
  for( localIndex const newLocalIndex : indexList )
  {
    GEOSX_ERROR_IF( object.m_localToGlobalMap[newLocalIndex] != -1,
                    "Local object " << newLocalIndex << " should be new but already has a global index "
                                    << object.m_localToGlobalMap[newLocalIndex] );

    object.m_localToGlobalMap[newLocalIndex] = object.m_maxGlobalIndex + glocalIndexOffset[thisRank] + nIndicesAssigned + 1;
    object.m_globalToLocalMap[object.m_localToGlobalMap[newLocalIndex]] = newLocalIndex;

    nIndicesAssigned += 1;
  }

  globalIndex maxGlobalIndex = -1;
  for( localIndex a=0 ; a<object.m_localToGlobalMap.size() ; ++a )
  {
    maxGlobalIndex = std::max( maxGlobalIndex, object.m_localToGlobalMap[a] );
  }

  MpiWrapper::allReduce( &maxGlobalIndex,
                         &(object.m_maxGlobalIndex),
                         1,
                         MPI_MAX,
                         MPI_COMM_GEOSX );
}

void
CommunicationTools::
  AssignNewGlobalIndices( ElementRegionManager & elementManager,
                          std::map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems )
{
  int const thisRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const commSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

  localIndex numberOfNewObjectsHere = 0;

  for( auto const & iter : newElems )
  {
    std::set< localIndex > const & indexList = iter.second;
    numberOfNewObjectsHere += indexList.size();
  }

  localIndex_array numberOfNewObjects( commSize );
  localIndex_array glocalIndexOffset( commSize );
  MpiWrapper::allGather( numberOfNewObjectsHere, numberOfNewObjects );


  glocalIndexOffset[0] = 0;
  for( int rank = 1 ; rank < commSize ; ++rank )
  {
    glocalIndexOffset[rank] = glocalIndexOffset[rank - 1] + numberOfNewObjects[rank - 1];
  }

  localIndex nIndicesAssigned = 0;
  globalIndex maxGlobalIndex = -1;

  for( auto const & iter : newElems )
  {
    localIndex const er = iter.first.first;
    localIndex const esr = iter.first.second;
    std::set< localIndex > const & indexList = iter.second;

    ElementSubRegionBase * const subRegion = elementManager.GetRegion( er )->GetSubRegion( esr );

    for( localIndex const newLocalIndex : indexList )
    {
      GEOSX_ERROR_IF( subRegion->m_localToGlobalMap[newLocalIndex] != -1,
                      "Local object " << newLocalIndex << " should be new but already has a global index "
                                      << subRegion->m_localToGlobalMap[newLocalIndex] );

      subRegion->m_localToGlobalMap[newLocalIndex] = elementManager.m_maxGlobalIndex + glocalIndexOffset[thisRank] + nIndicesAssigned + 1;
      subRegion->m_globalToLocalMap[subRegion->m_localToGlobalMap[newLocalIndex]] = newLocalIndex;

      nIndicesAssigned += 1;
    }
    for( localIndex a=0 ; a<subRegion->m_localToGlobalMap.size() ; ++a )
    {
      maxGlobalIndex = std::max( maxGlobalIndex, subRegion->m_localToGlobalMap[a] );
    }
  }


  MpiWrapper::allReduce( &maxGlobalIndex,
                         &(elementManager.m_maxGlobalIndex),
                         1,
                         MPI_MAX,
                         MPI_COMM_GEOSX );
}

void
CommunicationTools::
  FindMatchedPartitionBoundaryObjects( ObjectManagerBase * const group,
                                       array1d< NeighborCommunicator > & allNeighbors )//,
//array1d< array1d<localIndex> > & matchedPartitionBoundaryObjects )
{
  GEOSX_MARK_FUNCTION;
  integer_array & domainBoundaryIndicator = group->getReference< integer_array >( group->m_ObjectManagerBaseViewKeys.domainBoundaryIndicator );

  array1d< globalIndex > globalPartitionBoundaryObjectsIndices;
  group->ConstructGlobalListOfBoundaryObjects( globalPartitionBoundaryObjectsIndices );


//  array1d<NeighborCommunicator> & allNeighbors = this->getReference< array1d<NeighborCommunicator> >(
// viewKeys.neighbors );

  // send the size of the partitionBoundaryObjects to neighbors
  {
    array1d< array1d< globalIndex > > neighborPartitionBoundaryObjects( allNeighbors.size() );
//    matchedPartitionBoundaryObjects.resize( allNeighbors.size() );

    int commID = reserveCommID();

    for( localIndex i=0 ; i<allNeighbors.size() ; ++i )
    {
      allNeighbors[i].MPI_iSendReceive( globalPartitionBoundaryObjectsIndices,
                                        neighborPartitionBoundaryObjects[i],
                                        commID, MPI_COMM_GEOSX );
    }

    for( localIndex i=0 ; i<allNeighbors.size() ; ++i )
    {
      NeighborCommunicator const & neighbor = allNeighbors[i];
      localIndex_array &
      matchedPartitionBoundaryObjects = group->GetGroup( group->m_ObjectManagerBaseGroupKeys.neighborData )->
                                          GetGroup( std::to_string( neighbor.NeighborRank()))->
                                          getReference< localIndex_array >( group->m_ObjectManagerBaseViewKeys.matchedPartitionBoundaryObjects );

      allNeighbors[i].MPI_WaitAll( commID );
      localIndex localCounter = 0;
      localIndex neighborCounter = 0;
      while( localCounter < globalPartitionBoundaryObjectsIndices.size() &&
             neighborCounter < neighborPartitionBoundaryObjects[i].size() )
      {
        if( globalPartitionBoundaryObjectsIndices[localCounter] == neighborPartitionBoundaryObjects[i][neighborCounter] )
        {
          localIndex const localMatchedIndex = group->m_globalToLocalMap.at( globalPartitionBoundaryObjectsIndices[localCounter] );
          matchedPartitionBoundaryObjects.push_back( localMatchedIndex );
          domainBoundaryIndicator[ localMatchedIndex ] = 2;
          ++localCounter;
          ++neighborCounter;
        }
        else if( globalPartitionBoundaryObjectsIndices[localCounter] > neighborPartitionBoundaryObjects[i][neighborCounter] )
        {
          ++neighborCounter;
        }
        else
        {
          ++localCounter;
        }
      }
    }
    releaseCommID( commID );
  }
}



void CommunicationTools::FindGhosts( MeshLevel * const meshLevel,
                                     array1d< NeighborCommunicator > & neighbors,
                                     bool use_nonblocking )
{
  GEOSX_MARK_FUNCTION;
  int commID = CommunicationTools::reserveCommID();

  GEOSX_MARK_BEGIN( "Neighbor wait loop" );
  {
    int neighbor_count = neighbors.size( );
    auto send = [&] ( int idx )
      {
        neighbors[idx].PrepareAndSendGhosts( false, 1, meshLevel, commID );
        return neighbors[idx].GetSizeRecvRequest( commID );
      };
    auto post_recv = [&] ( int idx )
      {
        neighbors[idx].PostRecv( commID );
        return neighbors[idx].GetRecvRequest( commID );
      };
    auto proc_recv = [&] ( int idx )
      {
        neighbors[idx].UnpackGhosts( meshLevel, commID );
        return MPI_REQUEST_NULL;
      };
    std::vector< std::function< MPI_Request ( int ) > > phases = { send, post_recv, proc_recv };
    if( use_nonblocking )
    {
      MpiWrapper::ActiveWaitSomeCompletePhase( neighbor_count, phases );
    }
    else
    {
      MpiWrapper::ActiveWaitOrderedCompletePhase( neighbor_count, phases );
    }
  }
  GEOSX_MARK_END( "Neighbor wait loop" );

  meshLevel->getNodeManager()->SetReceiveLists();
  meshLevel->getEdgeManager()->SetReceiveLists();
  meshLevel->getFaceManager()->SetReceiveLists();

  // at present removing this barrier allows a nondeterministic mpi error to happen on lassen
  //   it occurs less than 5% of the time and happens when a process enters the recv phase in
  //   the sync list exchange while another process still has not unpacked the ghosts received from
  //   the first process. Depending on the mpi implementation the sync send from the first process
  //   can be recv'd by the second process instead of the ghost send which has already been sent but
  //   not necessarily recieved.
  // Some restructuring to ensure this can't happen ( can also probably just change the send/recv tagging )
  //   can eliminate this. But at present runtimes are the same in either case, as time is mostly just
  //   shifted from the waitall in UnpackAndRebuildSyncLists since the processes are more 'in-sync' when
  //   hitting that point after introducing this barrier.
  MpiWrapper::Barrier( );

  {
    int neighbor_count = neighbors.size( );
    auto send = [&] ( int idx )
      {
        neighbors[idx].PrepareAndSendSyncLists( meshLevel, commID );
        return neighbors[idx].GetSizeRecvRequest( commID );
      };
    auto post_recv  = [&] ( int idx )
      {
        neighbors[idx].PostRecv( commID );
        return neighbors[idx].GetRecvRequest( commID );
      };
    auto proc_recv = [&] ( int idx )
      {
        neighbors[idx].UnpackAndRebuildSyncLists( meshLevel, commID );
        return MPI_REQUEST_NULL;
      };
    std::vector< std::function< MPI_Request ( int ) > > phases = { send, post_recv, proc_recv };
    if( use_nonblocking )
    {
      MpiWrapper::ActiveWaitSomeCompletePhase( neighbor_count, phases );
    }
    else
    {
      MpiWrapper::ActiveWaitOrderedCompletePhase( neighbor_count, phases );
    }
  }

  meshLevel->getNodeManager()->FixUpDownMaps( false );
  meshLevel->getEdgeManager()->FixUpDownMaps( false );
  meshLevel->getFaceManager()->FixUpDownMaps( false );
  for( localIndex er=0 ; er<meshLevel->getElemManager()->numRegions() ; ++er )
  {
    ElementRegionBase * const elemRegion = meshLevel->getElemManager()->GetRegion( er );
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      ElementSubRegionBase * const subRegion = elemRegion->GetSubRegion( esr );
      subRegion->FixUpDownMaps( false );
    }
  }
  meshLevel->getNodeManager()->CompressRelationMaps();
  meshLevel->getEdgeManager()->CompressRelationMaps();
  CommunicationTools::releaseCommID( commID );
}


void CommunicationTools::SynchronizePackSendRecvSizes( const std::map< string, string_array > & fieldNames,
                                                       MeshLevel * const mesh,
                                                       array1d< NeighborCommunicator > & neighbors,
                                                       MPI_iCommData & icomm,
                                                       bool on_device )
{
  GEOSX_MARK_FUNCTION;
  icomm.fieldNames.insert( fieldNames.begin(), fieldNames.end() );
  icomm.resize( neighbors.size() );

  for( localIndex neighborIndex=0 ; neighborIndex<neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    int const bufferSize = neighbor.PackCommSizeForSync( fieldNames, mesh, icomm.commID, on_device );

    neighbor.MPI_iSendReceiveBufferSizes( icomm.commID,
                                          icomm.mpiSizeSendBufferRequest[neighborIndex],
                                          icomm.mpiSizeRecvBufferRequest[neighborIndex],
                                          MPI_COMM_GEOSX );

    neighbor.resizeSendBuffer( icomm.commID, bufferSize );
  }
}


void CommunicationTools::SynchronizePackSendRecv( const std::map< string, string_array > & fieldNames,
                                                  MeshLevel * const mesh,
                                                  array1d< NeighborCommunicator > & neighbors,
                                                  MPI_iCommData & icomm,
                                                  bool on_device )
{
  GEOSX_MARK_FUNCTION;

  MPI_iCommData sizeComm;
  for( localIndex neighborIndex=0 ; neighborIndex<neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    neighbor.PackCommBufferForSync( fieldNames, mesh, icomm.commID, on_device );
  }

  for( localIndex count=0 ; count<neighbors.size() ; ++count )
  {
    int neighborIndex;
    MpiWrapper::Waitany( icomm.size,
                         icomm.mpiSizeRecvBufferRequest.data(),
                         &neighborIndex,
                         icomm.mpiSizeRecvBufferStatus.data() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.MPI_iSendReceiveBuffers( icomm.commID,
                                      icomm.mpiSendBufferRequest[neighborIndex],
                                      icomm.mpiRecvBufferRequest[neighborIndex],
                                      MPI_COMM_GEOSX );
  }


}

void CommunicationTools::SynchronizeUnpack( MeshLevel * const mesh,
                                            array1d< NeighborCommunicator > & neighbors,
                                            MPI_iCommData & icomm,
                                            bool on_device )
{
  GEOSX_MARK_FUNCTION;

#if 0
//  forall_in_range<RAJA::omp_parallel_for_exec>( 0, neighbors.size() ,
//                                                   [&] ( int neighborIndex )->void

#pragma omp parallel for
  for( size_t neighborIndex=0 ; neighborIndex<neighbors.size() ; ++neighborIndex )
  {
    std::cout<<"thread "<<omp_get_thread_num()<<" of "<<omp_get_num_threads()<<std::endl;
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    MPI_Wait( &( icomm.mpiRecvBufferRequest[neighborIndex] ),
              &( icomm.mpiRecvBufferStatus[neighborIndex] ) );

    neighbor.UnpackBufferForSync( icomm.fieldNames, mesh, icomm.commID );

  }
#else

  // unpack the buffers
  for( localIndex count=0 ; count<neighbors.size() ; ++count )
  {
    int neighborIndex;
    MpiWrapper::Waitany( icomm.size,
                         icomm.mpiRecvBufferRequest.data(),
                         &neighborIndex,
                         icomm.mpiRecvBufferStatus.data() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    neighbor.UnpackBufferForSync( icomm.fieldNames, mesh, icomm.commID, on_device );
  }

#endif
  MpiWrapper::Waitall( icomm.size,
                       icomm.mpiSizeSendBufferRequest.data(),
                       icomm.mpiSizeSendBufferStatus.data() );

  MpiWrapper::Waitall( icomm.size,
                       icomm.mpiSendBufferRequest.data(),
                       icomm.mpiSendBufferStatus.data() );

}

void CommunicationTools::SynchronizeFields( const std::map< string, string_array > & fieldNames,
                                            MeshLevel * const mesh,
                                            array1d< NeighborCommunicator > & neighbors,
                                            bool on_device )
{
  MPI_iCommData icomm;
  SynchronizePackSendRecvSizes( fieldNames, mesh, neighbors, icomm, on_device );
  SynchronizePackSendRecv( fieldNames, mesh, neighbors, icomm, on_device );
  SynchronizeUnpack( mesh, neighbors, icomm, on_device );
}


} /* namespace geosx */
