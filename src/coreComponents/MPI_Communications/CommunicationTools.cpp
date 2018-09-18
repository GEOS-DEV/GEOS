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
 * CommunicationTools.cpp
 *
 *  Created on: Jan 6, 2018
 *      Author: settgast
 */

#include "CommunicationTools.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "NeighborCommunicator.hpp"
#include "managers/DomainPartition.hpp"
#include "common/TimingMacros.hpp"

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

int CommunicationTools::MPI_Size( MPI_Comm const & comm )
{
  int size;
  MPI_Comm_size( comm, &size );
  return size;
}

int CommunicationTools::MPI_Rank( MPI_Comm const & comm )
{
  int rank;
  MPI_Comm_rank( comm, &rank );
  return rank;
}


std::set<int> & CommunicationTools::getFreeCommIDs()
{
  static std::set<int> commIDs;
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
  std::set<int> & commIDs = getFreeCommIDs();

  int rval = *( commIDs.begin() );
  commIDs.erase( rval );
  return rval;
}

void CommunicationTools::releaseCommID( int & ID )
{
  std::set<int> & commIDs = getFreeCommIDs();

  if( commIDs.count( ID ) > 0 )
  {
    GEOS_ERROR( "Attempting to release commID that is already free" );
  }
  commIDs.insert( ID );
  ID = -1;
}

void CommunicationTools::AssignGlobalIndices( ObjectManagerBase & object,
                                              ObjectManagerBase const & compositionObject,
                                              array1d<NeighborCommunicator> & neighbors )
{
  GEOSX_MARK_FUNCTION;

  integer_array & ghostRank = object.getReference<integer_array>( object.m_ObjectManagerBaseViewKeys.ghostRank );
  ghostRank = -2;

  int const commSize = MPI_Size( MPI_COMM_GEOSX );
  localIndex numberOfObjectsHere = object.size();
  localIndex_array numberOfObjects( commSize );
  localIndex_array glocalIndexOffset( commSize );
  MPI_Allgather( reinterpret_cast<char*>( &numberOfObjectsHere ),
                 sizeof(localIndex),
                 MPI_CHAR,
                 reinterpret_cast<char*>( numberOfObjects.data() ),
                 sizeof(localIndex),
                 MPI_CHAR,
                 MPI_COMM_GEOSX );

  int const commRank = MPI_Rank( MPI_COMM_GEOSX );

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
  array1d<globalIndex_array> objectToCompositionObject;
  object.ExtractMapFromObjectForAssignGlobalIndexNumbers( &compositionObject, objectToCompositionObject );

  // now arrange the data from objectToCompositionObject into a map "indexByFirstCompositionIndex", such that the key
  // is the lowest global index of the composition object that make up this object. The value of the map is a pair, with
  // the
  // array being the remaining composition object global indices, and the second being the global index of the object
  // itself.
  map<globalIndex, array1d<std::pair<globalIndex_array, localIndex> > > indexByFirstCompositionIndex;

//  for( array1d<globalIndex_array>::const_iterator a = objectToCompositionObject.begin() ;
//      a != objectToCompositionObject.end() ;
//      ++a )

  localIndex bufferSize = 0;
  for( localIndex a = 0 ; a < objectToCompositionObject.size() ; ++a )
  {
    // set nodelist array
    globalIndex_array const & nodeList = objectToCompositionObject[a];

    // grab the first global index of the composition objects
    const globalIndex firstCompositionIndex = nodeList[0];

    // create a temporary to hold the pair
    std::pair<globalIndex_array, globalIndex> tempComp;

    // fill the array with the remaining composition object global indices
    tempComp.first.insert( tempComp.first.begin(), nodeList.begin() + 1, nodeList.end() );

    // set the second value of the pair to the localIndex of the object.
    tempComp.second = a;

    // push the tempComp onto the map.
    indexByFirstCompositionIndex[firstCompositionIndex].push_back( tempComp );
    bufferSize += 2 + nodeList.size();
  }

  globalIndex_array objectToCompositionObjectSendBuffer;
  objectToCompositionObjectSendBuffer.reserve( bufferSize );

  // put the map into a buffer
  for( localIndex a = 0 ; a < objectToCompositionObject.size() ; ++a )
  {
    globalIndex_array const & nodeList = objectToCompositionObject[a];
    objectToCompositionObjectSendBuffer.push_back( nodeList.size() );
    objectToCompositionObjectSendBuffer.push_back( object.m_localToGlobalMap[a] );
    for( localIndex b = 0 ; b < nodeList.size() ; ++b )
    {
      objectToCompositionObjectSendBuffer.push_back( nodeList[b] );
    }
  }

  int commID = reserveCommID();

  // send the composition buffers
  {
    int const sendSize = integer_conversion<int const>( objectToCompositionObjectSendBuffer.size() * sizeof(globalIndex));

    for( localIndex in = 0 ; in < neighbors.size() ; ++in )
    {
      NeighborCommunicator & neighbor = neighbors[in];

      neighbor.MPI_iSendReceive( reinterpret_cast<const char*>( objectToCompositionObjectSendBuffer.data() ),
                                 sendSize,
                                 commID,
                                 MPI_COMM_GEOSX );
    }
    for( localIndex in = 0 ; in < neighbors.size() ; ++in )
    {
      neighbors[in].MPI_WaitAll( commID );
    }
  }

  // unpack the data from neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager] to
  // the local arrays

  // object to receive the neighbor data
  // this baby is and Array (for each neighbor) of maps, with the key of lowest composition index, and a value
  // containing
  // an array containing the std::pairs of the remaining composition indices, and the globalIndex of the object.
  array1d<map<globalIndex, array1d<std::pair<globalIndex_array, globalIndex> > > >
  neighborCompositionObjects( neighbors.size() );

  {
    for( localIndex neighborIndex = 0 ; neighborIndex < neighbors.size() ; ++neighborIndex )
    {
      NeighborCommunicator & neighbor = neighbors[neighborIndex];

      globalIndex const * recBuffer = reinterpret_cast<globalIndex const *>( neighbor.ReceiveBuffer( commID ).data() );
      localIndex recBufferSize = integer_conversion<localIndex>( neighbor.ReceiveBuffer( commID ).size() / sizeof(globalIndex));
      globalIndex const * endBuffer = recBuffer + recBufferSize;
      // iterate over data that was just received
      while( recBuffer < endBuffer )
      {
        // the first thing packed was the data size for a given object
        localIndex dataSize = integer_conversion<localIndex>( *recBuffer++ );

        // the second thing packed was the globalIndex of that object
        const globalIndex neighborGlobalIndex = *( recBuffer++ );

        // the global indices of the composition objects were next. they are ordered, so the lowest one is first.
        const globalIndex firstCompositionIndex = *( recBuffer++ );

        // the remaining composition object indices.
        globalIndex_array temp;
        for( localIndex b = 1 ; b < dataSize ; ++b )
        {
          temp.push_back( *( recBuffer++ ) );
        }

        // fill neighborCompositionObjects
        std::pair<globalIndex_array, globalIndex> tempComp( std::make_pair( temp, neighborGlobalIndex ) );
        neighborCompositionObjects[neighborIndex][firstCompositionIndex].push_back( tempComp );
      }
    }
  }
  releaseCommID( commID );

  // now check to see if the global index is valid. We do this by checking the contents of neighborCompositionObjects
  // with indexByFirstCompositionIndex
  for( localIndex neighborIndex = 0 ; neighborIndex < neighbors.size() ; ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    // it only matters if the neighbor rank is lower than this rank
//    if( neighbor.NeighborRank() < commRank )
    {
      // Set iterators to the beginning of each indexByFirstCompositionIndex,
      // and neighborCompositionObjects[neighborNum].
      map<globalIndex, array1d<std::pair<globalIndex_array, localIndex> > >::const_iterator
        iter_local = indexByFirstCompositionIndex.begin();
      map<globalIndex, array1d<std::pair<globalIndex_array, globalIndex> > >::const_iterator
        iter_neighbor = neighborCompositionObjects[neighborIndex].begin();

      // now we continue the while loop as long as both of our iterators are in range.
      while( iter_local != indexByFirstCompositionIndex.end() &&
             iter_neighbor != neighborCompositionObjects[neighborIndex].end() )
      {
        // check to see if the map keys (first composition index) are the same.
        if( iter_local->first == iter_neighbor->first )
        {
          // first we loop over all local composition arrays (objects with the matched key)
          for( array1d<std::pair<globalIndex_array, localIndex> >::const_iterator
               iter_local2 = iter_local->second.begin() ;
               iter_local2 != iter_local->second.end() ; ++iter_local2 )
          {
            // and loop over all of the neighbor composition arrays (objects with the matched key)
            for( array1d<std::pair<globalIndex_array, globalIndex> >::const_iterator
                 iter_neighbor2 = iter_neighbor->second.begin() ;
                 iter_neighbor2 != iter_neighbor->second.end() ;
                 ++iter_neighbor2 )
            {
              // now compare the composition arrays
              if( iter_local2->first.size() == iter_neighbor2->first.size() &&
                  std::equal( iter_local2->first.begin(), iter_local2->first.end(), iter_neighbor2->first.begin() ) )
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
  }

  object.ConstructGlobalToLocalMap();
}

void
CommunicationTools::
FindMatchedPartitionBoundaryObjects( ObjectManagerBase * const group,
                                     array1d<NeighborCommunicator> & allNeighbors )//,
//array1d< array1d<localIndex> > & matchedPartitionBoundaryObjects )
{
  GEOSX_MARK_FUNCTION;

  integer_array const & ghostRank = group->getReference<integer_array>( group->m_ObjectManagerBaseViewKeys.ghostRank );
  integer_array & domainBoundaryIndicator = group->getReference<integer_array>( group->m_ObjectManagerBaseViewKeys.domainBoundaryIndicator );
  globalIndex_array const & localToGlobal = group->getReference<globalIndex_array>( group->m_ObjectManagerBaseViewKeys.localToGlobalMap );

  array1d<globalIndex> globalPartitionBoundaryObjectsIndices;
  group->ConstructGlobalListOfBoundaryObjects( globalPartitionBoundaryObjectsIndices );


//  array1d<NeighborCommunicator> & allNeighbors = this->getReference< array1d<NeighborCommunicator> >(
// viewKeys.neighbors );

  // send the size of the partitionBoundaryObjects to neighbors
  {
    array1d< array1d<globalIndex> > neighborPartitionBoundaryObjects( allNeighbors.size() );
//    matchedPartitionBoundaryObjects.resize( allNeighbors.size() );

    int commID = reserveCommID();

    for( int i=0 ; i<allNeighbors.size() ; ++i )
    {
      allNeighbors[i].MPI_iSendReceive( globalPartitionBoundaryObjectsIndices,
                                        neighborPartitionBoundaryObjects[i],
                                        commID, MPI_COMM_GEOSX );
    }

    for( int i=0 ; i<allNeighbors.size() ; ++i )
    {
      NeighborCommunicator const & neighbor = allNeighbors[i];
      localIndex_array &
      matchedPartitionBoundaryObjects = group->GetGroup( group->m_ObjectManagerBaseGroupKeys.neighborData )->
                                        GetGroup( std::to_string( neighbor.NeighborRank()))->
                                        getReference<localIndex_array>( group->m_ObjectManagerBaseViewKeys.matchedPartitionBoundaryObjects );

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
                                     array1d<NeighborCommunicator> & neighbors )
{
  GEOSX_MARK_FUNCTION;
  int commID = CommunicationTools::reserveCommID();

  for( auto & neighbor : neighbors )
  {
    neighbor.FindAndPackGhosts( false, 2, meshLevel, commID );
  }

  for( auto & neighbor : neighbors )
  {
    neighbor.MPI_WaitAll( commID );
    neighbor.UnpackGhosts( meshLevel, commID );
  }
  meshLevel->getNodeManager()->SetReceiveLists();
  meshLevel->getEdgeManager()->SetReceiveLists();
  meshLevel->getFaceManager()->SetReceiveLists();

  for( auto & neighbor : neighbors )
  {
    neighbor.RebuildSyncLists( meshLevel, commID );
  }

  CommunicationTools::releaseCommID( commID );
}



void CommunicationTools::SynchronizeFields( const std::map<string, string_array >& fieldNames,
                                            MeshLevel * const mesh,
                                            array1d<NeighborCommunicator> & neighbors )
{

  int commID = CommunicationTools::reserveCommID();

  for( auto & neighbor : neighbors )
  {
    neighbor.PackBufferForSync( fieldNames, mesh, commID );
  }


  for( auto & neighbor : neighbors )
  {
    neighbor.MPI_WaitAll( commID );
    neighbor.UnpackBufferForSync( fieldNames, mesh, commID );
  }
  CommunicationTools::releaseCommID( commID );

//
//  // send and receive buffers
//  for( int neighborIndex=0 ; neighborIndex<m_neighbors.size() ; ++neighborIndex )
//  {
//    NeighborCommunication& neighbor = m_neighbors[neighborIndex];
//    neighbor.PackBuffer( fieldNames, commID  );
//    neighbor.SendReceiveBuffers( commID, mpiSendBufferRequest[neighborIndex], mpiRecvBufferRequest[neighborIndex] );
//
//  }
//
//
//  // unpack the buffers
//  for( int count=0 ; count<m_neighbors.size() ; ++count )
//  {
//    int neighborIndex;
//    MPI_Waitany( mpiRecvBufferRequest.size(), mpiRecvBufferRequest.data(), &neighborIndex, mpiRecvBufferStatus.data()
// );
//
//    NeighborCommunication& neighbor = this->m_neighbors[neighborIndex];
//    neighbor.UnpackBuffer( fieldNames );
//  }
//
//  MPI_Waitall( mpiSendBufferRequest.size(), mpiSendBufferRequest.data(), mpiSendBufferStatus.data() );
//

}


} /* namespace geosx */
