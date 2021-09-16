/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EmebeddedSurfacesParallelSynchronization.cpp
 */

#include "EmbeddedSurfacesParallelSynchronization.hpp"

#include "common/GeosxMacros.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ExtrinsicMeshData.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"


namespace geosx
{

using namespace dataRepository;

namespace
{
void packNewObjectsToGhosts( NeighborCommunicator * const neighbor,
                             int commID,
                             MeshLevel & mesh,
                             NewObjectLists & newObjects )
{
  EmbeddedSurfaceNodeManager & nodeManager = mesh.getEmbSurfNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  int neighborRank = neighbor->neighborRank();

  localIndex_array newNodesToSend;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > newSurfaceGhostsToSend;

  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > newElemsToSend;
  array1d< array1d< localIndex_array > > newElemsToSendData;

  arrayView1d< localIndex const > const & edgeGhostsToSend = edgeManager.getNeighborData( neighbor->neighborRank() ).ghostsToSend();
  arrayView1d< localIndex > const & parentIndex = nodeManager.getExtrinsicData< extrinsicMeshData::ParentEdgeIndex >();

  for( auto const ni : newObjects.newNodes )
  {
    forAll< serialPolicy >( edgeGhostsToSend.size(), [=, &newNodesToSend] ( localIndex const a )
    {
      if( edgeGhostsToSend[a] == parentIndex[ni] )
      {
        // a node is sent if the edge on which it was created had to be sent.
        newNodesToSend.emplace_back( ni );
      }
    } );
  }

  newElemsToSendData.resize( elemManager.numRegions() );
  newElemsToSend.resize( elemManager.numRegions() );
  //
  if( newObjects.newElements.size() > 0 )  // this map may be completely empty on some ranks.
  {
    elemManager.forElementSubRegionsComplete< EmbeddedSurfaceSubRegion >(
      [&]( localIndex const er, localIndex const esr, ElementRegionBase &, EmbeddedSurfaceSubRegion & subRegion )
    {

      FixedToManyElementRelation const & surfaceElementsToCells = subRegion.getToCellRelation();

      for( localIndex const & k : newObjects.newElements.at( {er, esr} ) )
      {
        localIndex const elemRegionIndex  = surfaceElementsToCells.m_toElementRegion[k][0];
        localIndex const elemSubRegionIndex = surfaceElementsToCells.m_toElementSubRegion[k][0];
        localIndex const elemIndex = surfaceElementsToCells.m_toElementIndex[k][0];

        CellElementRegion const & elemRegion = elemManager.getRegion< CellElementRegion >( elemRegionIndex );
        CellElementSubRegion const & elemSubRegion = elemRegion.getSubRegion< CellElementSubRegion >( elemSubRegionIndex );

        arrayView1d< localIndex const > const & elemGhostsToSend =
          elemSubRegion.getNeighborData( neighborRank ).ghostsToSend();

        forAll< serialPolicy >( elemGhostsToSend.size(), [=, &newSurfaceGhostsToSend] ( localIndex const a )
        {
          if( elemIndex == elemGhostsToSend[a] )
          {
            newSurfaceGhostsToSend[ {er, esr} ].insert( k );
            return;
          }
        } );
      }
    } );
  }

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newElemsToSendData[er].resize( elemRegion.numSubRegions() );
    newElemsToSend[er].resize( elemRegion.numSubRegions() );
    if( newSurfaceGhostsToSend.size() > 0 )  // this map may be completely empty on some ranks.
    {
      elemRegion.forElementSubRegionsIndex< EmbeddedSurfaceSubRegion >( [&]( localIndex const esr,
                                                                             EmbeddedSurfaceSubRegion & subRegion )
      {
        localIndex_array & surfaceElemGhostsToSend = subRegion.getNeighborData( neighborRank ).ghostsToSend();
        surfaceElemGhostsToSend.move( LvArray::MemorySpace::host );

        for( localIndex const & k : newSurfaceGhostsToSend.at( {er, esr} ) )
        {
          newElemsToSendData[er][esr].emplace_back( k );
          surfaceElemGhostsToSend.emplace_back( k );
        }
        newElemsToSend[er][esr] = newElemsToSendData[er][esr];
      } );
    }
  }



  parallelDeviceEvents sizeEvents;
  int bufferSize = 0;

  bufferSize += nodeManager.packGlobalMapsSize( newNodesToSend, 0 );
  bufferSize += elemManager.PackGlobalMapsSize( newElemsToSend );

  bufferSize += nodeManager.packUpDownMapsSize( newNodesToSend );
  bufferSize += elemManager.PackUpDownMapsSize( newElemsToSend );

  bufferSize += nodeManager.packSize( {}, newNodesToSend, 0, false, sizeEvents );
  bufferSize += elemManager.PackSize( {}, newElemsToSend );

  neighbor->resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor->sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  parallelDeviceEvents packEvents;
  int packedSize = 0;

  packedSize += nodeManager.packGlobalMaps( sendBufferPtr, newNodesToSend, 0 );
  packedSize += elemManager.PackGlobalMaps( sendBufferPtr, newElemsToSend );

  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, newNodesToSend );
  packedSize += elemManager.PackUpDownMaps( sendBufferPtr, newElemsToSend );

  packedSize += nodeManager.pack( sendBufferPtr, {}, newNodesToSend, 0, false, packEvents );
  packedSize += elemManager.Pack( sendBufferPtr, {}, newElemsToSend );


  GEOSX_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );


  GEOSX_UNUSED_VAR( commID );
//  neighbor->mpiISendReceive( commID, MPI_COMM_GEOSX );
}

void unpackNewToGhosts( NeighborCommunicator * const neighbor,
                        int commID,
                        MeshLevel & mesh )
{
  int unpackedSize = 0;

  EmbeddedSurfaceNodeManager & nodeManager = mesh.getEmbSurfNodeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  buffer_type const & receiveBuffer = neighbor->receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  localIndex_array newGhostNodes;
  ElementRegionManager::ElementReferenceAccessor< localIndex_array > newGhostElems;
  array1d< array1d< localIndex_array > > newGhostElemsData;
  newGhostElems.resize( elemManager.numRegions() );
  newGhostElemsData.resize( elemManager.numRegions() );


  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    newGhostElemsData[er].resize( elemRegion.numSubRegions() );
    newGhostElems[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      newGhostElems[er][esr].set( newGhostElemsData[er][esr] );
    }
  }

  parallelDeviceEvents events;

  unpackedSize += nodeManager.unpackGlobalMaps( receiveBufferPtr, newGhostNodes, 0 );
  unpackedSize += elemManager.UnpackGlobalMaps( receiveBufferPtr, newGhostElems );

  unpackedSize += nodeManager.unpackUpDownMaps( receiveBufferPtr, newGhostNodes, true, true );
  unpackedSize += elemManager.UnpackUpDownMaps( receiveBufferPtr, newGhostElems, true );

  unpackedSize += nodeManager.unpack( receiveBufferPtr, newGhostNodes, 0, false, events );
  unpackedSize += elemManager.Unpack( receiveBufferPtr, newGhostElems );

  waitAllDeviceEvents( events );

  elemManager.forElementSubRegionsComplete< ElementSubRegionBase >(
    [&]( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & subRegion )
  {
    localIndex_array & elemGhostsToReceive = subRegion.getNeighborData( neighbor->neighborRank() ).ghostsToReceive();

    if( newGhostElemsData[er][esr].size() > 0 )
    {
      forAll< serialPolicy >( newGhostElemsData[er][esr].size(), [=, &elemGhostsToReceive] ( localIndex const k )
      {
        localIndex const & newElemIndex = newGhostElemsData[er][esr][k];
        elemGhostsToReceive.emplace_back( newElemIndex );
      } );
    }
  } );
}

void packFracturedToGhosts( NeighborCommunicator * const neighbor,
                            int commID,
                            MeshLevel & mesh,
                            string const fractureRegionName )
{
  ElementRegionManager & elemManager = mesh.getElemManager();

  int neighborRank = neighbor->neighborRank();

  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > > elemsToSend;
  array1d< array1d< localIndex_array > > elemsToSendData;

  elemsToSendData.resize( elemManager.numRegions() );
  elemsToSend.resize( elemManager.numRegions() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    elemsToSendData[er].resize( elemRegion.numSubRegions() );
    elemsToSend[er].resize( elemRegion.numSubRegions() );

    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                       CellElementSubRegion & subRegion )
    {
      // we send all the ghosts
      arrayView1d< localIndex const >  const & elemsGhostsToSend = subRegion.getNeighborData( neighborRank ).ghostsToSend();
      elemsGhostsToSend.move( LvArray::MemorySpace::host );
      forAll< serialPolicy >( elemsGhostsToSend.size(), [=, &elemsToSendData] ( localIndex const k )
      {
        elemsToSendData[er][esr].emplace_back( elemsGhostsToSend[k] );
      } );
      elemsToSend[er][esr] = elemsToSendData[er][esr];
    } );
  }

  parallelDeviceEvents sizeEvents;
  int bufferSize = 0;

  bufferSize += elemManager.packFracturedElementsSize( elemsToSend, fractureRegionName );

  neighbor->resizeSendBuffer( commID, bufferSize );

  buffer_type & sendBuffer = neighbor->sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  parallelDeviceEvents packEvents;
  int packedSize = 0;

  packedSize += elemManager.packFracturedElements( sendBufferPtr, elemsToSend, fractureRegionName );

  GEOSX_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );
}

void unpackFracturedToGhosts( NeighborCommunicator * const neighbor,
                              int commID,
                              MeshLevel & mesh,
                              string const fractureRegionName )
{
  int unpackedSize = 0;

  ElementRegionManager & elemManager = mesh.getElemManager();

  buffer_type const & receiveBuffer = neighbor->receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  ElementRegionManager::ElementReferenceAccessor< localIndex_array > ghostElems;
  array1d< array1d< localIndex_array > > ghostElemsData;
  ghostElems.resize( elemManager.numRegions() );
  ghostElemsData.resize( elemManager.numRegions() );

  for( localIndex er=0; er<elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    ghostElemsData[er].resize( elemRegion.numSubRegions() );
    ghostElems[er].resize( elemRegion.numSubRegions() );
    for( localIndex esr=0; esr<elemRegion.numSubRegions(); ++esr )
    {
      ghostElems[er][esr].set( ghostElemsData[er][esr] );
    }
  }

  parallelDeviceEvents events;

  unpackedSize += elemManager.unpackFracturedElements( receiveBufferPtr, ghostElems, fractureRegionName );

  waitAllDeviceEvents( events );
}


}

void EmebeddedSurfacesParallelSynchronization::synchronizeNewSurfaces( MeshLevel & mesh,
                                                                       std::vector< NeighborCommunicator > & neighbors,
                                                                       NewObjectLists & newObjects,
                                                                       int const mpiCommOrder )
{
  //************************************************************************************************
  // We need to send over the new embedded surfaces and related objects for those whose parents are ghosts on neighbors.

  MPI_iCommData commData( CommunicationTools::getInstance().getCommID() );
  commData.resize( neighbors.size());
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewObjectsToGhosts( &neighbor,
                            commData.commID(),
                            mesh,
                            newObjects );

    neighbor.mpiISendReceiveBufferSizes( commData.commID(),
                                         commData.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOSX );

  }

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData.commID(),
                                     commData.mpiSendBufferRequest( neighborIndex ),
                                     commData.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOSX );
  }


  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;
    if( mpiCommOrder == 0 )
    {
      MpiWrapper::waitAny( commData.size(),
                           commData.mpiRecvBufferRequest(),
                           &neighborIndex,
                           commData.mpiRecvBufferStatus() );
    }
    else
    {
      MpiWrapper::wait( commData.mpiRecvBufferRequest() + count,
                        commData.mpiRecvBufferStatus() + count );
    }

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    unpackNewToGhosts( &neighbor, commData.commID(), mesh );
  }
}

void EmebeddedSurfacesParallelSynchronization::synchronizeFracturedElements( MeshLevel & mesh,
                                                                             std::vector< NeighborCommunicator > & neighbors,
                                                                             string const fractureRegionName )
{
  MPI_iCommData commData( CommunicationTools::getInstance().getCommID() );
  commData.resize( neighbors.size());
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packFracturedToGhosts( &neighbor,
                           commData.commID(),
                           mesh,
                           fractureRegionName );

    neighbor.mpiISendReceiveBufferSizes( commData.commID(),
                                         commData.mpiSendBufferSizeRequest( neighborIndex ),
                                         commData.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOSX );

  }

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData.commID(),
                                     commData.mpiSendBufferRequest( neighborIndex ),
                                     commData.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOSX );
  }


  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;

    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    unpackFracturedToGhosts( &neighbor, commData.commID(), mesh, fractureRegionName );
  }
}


} /* namespace geosx */
