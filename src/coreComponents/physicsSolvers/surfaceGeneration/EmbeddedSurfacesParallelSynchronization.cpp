/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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
#include "mpiCommunications/CommunicationTools.hpp"


namespace geosx
{

using namespace dataRepository;

EmebeddedSurfacesParallelSynchronization::EmebeddedSurfacesParallelSynchronization()
{}

EmebeddedSurfacesParallelSynchronization::~EmebeddedSurfacesParallelSynchronization()
{}

void EmebeddedSurfacesParallelSynchronization::synchronizeNewSurfaces( MeshLevel & mesh,
                                                                       std::vector< NeighborCommunicator > & neighbors )
{

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  //************************************************************************************************
  // We need to send over the new embedded surface objects whose parents are ghosts on neighbors.

  MPI_iCommData commData;
  commData.resize( neighbors.size());
  for( unsigned int neighborIndex=0; neighborIndex<neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    packNewModifiedObjectsToGhosts( &neighbor,
                                    commData.commID,
                                    mesh,
                                    modifiedObjects );

    neighbor.mpiISendReceiveBufferSizes( commData.commID,
                                         commData.mpiSizeSendBufferRequest[neighborIndex],
                                         commData.mpiSizeRecvBufferRequest[neighborIndex],
                                         MPI_COMM_GEOSX );

  }

  for( unsigned int count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitany( commData.size,
                         commData.mpiSizeRecvBufferRequest.data(),
                         &neighborIndex,
                         commData.mpiSizeRecvBufferStatus.data() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( commData.commID,
                                     commData.mpiSendBufferRequest[neighborIndex],
                                     commData.mpiRecvBufferRequest[neighborIndex],
                                     MPI_COMM_GEOSX );
  }


  for( unsigned int count=0; count<neighbors.size(); ++count )
  {

    int neighborIndex = count;
    if( mpiCommOrder == 0 )
    {
      MpiWrapper::waitany( commData.size,
                           commData.mpiRecvBufferRequest.data(),
                           &neighborIndex,
                           commData.mpiRecvBufferStatus.data() );
    }
    else
    {
      MpiWrapper::wait( commData.mpiRecvBufferRequest.data() + count,
                        commData.mpiRecvBufferStatus.data() + count );
    }

    NeighborCommunicator & neighbor = neighbors[neighborIndex];


    unpackNewModToGhosts( &neighbor, commData2.commID, mesh, receivedObjects );
  }

  modifiedObjects.insert( receivedObjects );

  nodeManager.fixUpDownMaps( false );
  edgeManager.fixUpDownMaps( false );
  faceManager.fixUpDownMaps( false );


  for( localIndex er = 0; er < elemManager.numRegions(); ++er )
  {
    ElementRegionBase & elemRegion = elemManager.getRegion( er );
    for( localIndex esr = 0; esr < elemRegion.numSubRegions(); ++esr )
    {
      elemRegion.getSubRegion( esr ).fixUpDownMaps( false );
    }
  }
}


} /* namespace geosx */
