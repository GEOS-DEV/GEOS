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
 * @file CommunicationTools.hpp
 */

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_

#include "CommID.hpp"

#include "common/MpiWrapper.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "mesh/FieldIdentifiers.hpp"

#include <set>

namespace geosx
{


class ObjectManagerBase;
class NodeManager;
class NeighborCommunicator;
class MeshLevel;
class ElementRegionManager;

class MPI_iCommData;



class CommunicationTools
{
public:
  CommunicationTools();
  ~CommunicationTools();

  static CommunicationTools & getInstance();

  void assignGlobalIndices( ObjectManagerBase & object,
                            NodeManager const & compositionObject,
                            std::vector< NeighborCommunicator > & neighbors );

  void assignNewGlobalIndices( ObjectManagerBase & object,
                               std::set< localIndex > const & indexList );

  void assignNewGlobalIndices( ElementRegionManager & elementManager,
                               std::map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems );

  void setupGhosts( MeshLevel & meshLevel,
                    std::vector< NeighborCommunicator > & neighbors,
                    bool use_nonblocking );

  CommID getCommID()
  { return CommID( m_freeCommIDs ); }

  void findMatchedPartitionBoundaryObjects( ObjectManagerBase & group,
                                            std::vector< NeighborCommunicator > & allNeighbors );

  void synchronizeFields( FieldIdentifiers const & fieldsToBeSync,
                          MeshLevel & mesh,
                          std::vector< NeighborCommunicator > & allNeighbors,
                          bool onDevice );

  void synchronizePackSendRecvSizes( FieldIdentifiers const & fieldsToBeSync,
                                     MeshLevel & mesh,
                                     std::vector< NeighborCommunicator > & neighbors,
                                     MPI_iCommData & icomm,
                                     bool onDevice );

  void synchronizePackSendRecv( FieldIdentifiers const & fieldsToBeSync,
                                MeshLevel & mesh,
                                std::vector< NeighborCommunicator > & allNeighbors,
                                MPI_iCommData & icomm,
                                bool onDevice );

  void asyncPack( FieldIdentifiers const & fieldsToBeSync,
                  MeshLevel & mesh,
                  std::vector< NeighborCommunicator > & neighbors,
                  MPI_iCommData & icomm,
                  bool onDevice,
                  parallelDeviceEvents & events );

  void asyncSendRecv( std::vector< NeighborCommunicator > & neighbors,
                      MPI_iCommData & icomm,
                      bool onDevice,
                      parallelDeviceEvents & events );

  void synchronizeUnpack( MeshLevel & mesh,
                          std::vector< NeighborCommunicator > & neighbors,
                          MPI_iCommData & icomm,
                          bool onDevice );

  bool asyncUnpack( MeshLevel & mesh,
                    std::vector< NeighborCommunicator > & neighbors,
                    MPI_iCommData & icomm,
                    bool onDevice,
                    parallelDeviceEvents & events,
                    MPI_Op op=MPI_REPLACE);

  void finalizeUnpack( MeshLevel & mesh,
                       std::vector< NeighborCommunicator > & neighbors,
                       MPI_iCommData & icomm,
                       bool onDevice,
                       parallelDeviceEvents & events,
                       MPI_Op op=MPI_REPLACE);

private:
  std::set< int > m_freeCommIDs;
  static CommunicationTools * m_instance;


};



} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_ */
