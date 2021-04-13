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
 * @file CommunicationTools.hpp
 */

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <set>

namespace geosx
{


class ObjectManagerBase;
class NodeManager;
class NeighborCommunicator;
class MeshLevel;
class ElementRegionManager;

class MPI_iCommData;

class CommID
{
public:
  CommID( std::set< int > & freeIDs ):
    m_freeIDs( freeIDs ),
    m_id( -1 )
  {
    GEOSX_ERROR_IF_EQ( freeIDs.size(), 0 );
    m_id = *freeIDs.begin();
    freeIDs.erase( freeIDs.begin() );
  }

  CommID( CommID && src ):
    m_freeIDs( src.m_freeIDs ),
    m_id( src.m_id )
  {
    src.m_id = -1;
  }

  ~CommID()
  {
    if( m_id < 0 )
    { return; }

    GEOSX_ERROR_IF( m_freeIDs.count( m_id ) > 0, "Attempting to release commID that is already free: " << m_id );

    m_freeIDs.insert( m_id );
    m_id = -1;
  }

  CommID( CommID const & ) = delete;
  CommID & operator=( CommID const & ) = delete;
  CommID & operator=( CommID && ) = delete;

  constexpr operator int()
  { return m_id; }

private:
  std::set< int > & m_freeIDs;
  int m_id = -1;
};


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

  void findGhosts( MeshLevel & meshLevel,
                   std::vector< NeighborCommunicator > & neighbors,
                   bool use_nonblocking );

  CommID getCommID()
  { return CommID( m_freeCommIDs ); }

  void findMatchedPartitionBoundaryObjects( ObjectManagerBase & group,
                                            std::vector< NeighborCommunicator > & allNeighbors );

  void synchronizeFields( const std::map< string, string_array > & fieldNames,
                          MeshLevel & mesh,
                          std::vector< NeighborCommunicator > & allNeighbors,
                          bool onDevice );

  void synchronizePackSendRecvSizes( const std::map< string, string_array > & fieldNames,
                                     MeshLevel & mesh,
                                     std::vector< NeighborCommunicator > & neighbors,
                                     MPI_iCommData & icomm,
                                     bool onDevice );

  void synchronizePackSendRecv( const std::map< string, string_array > & fieldNames,
                                MeshLevel & mesh,
                                std::vector< NeighborCommunicator > & allNeighbors,
                                MPI_iCommData & icomm,
                                bool onDevice );

  void asyncPack( const std::map< string, string_array > & fieldNames,
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
                    parallelDeviceEvents & events );

  void finalizeUnpack( MeshLevel & mesh,
                       std::vector< NeighborCommunicator > & neighbors,
                       MPI_iCommData & icomm,
                       bool onDevice,
                       parallelDeviceEvents & events );

private:
  std::set< int > m_freeCommIDs;
  static CommunicationTools * m_instance;


};


class MPI_iCommData
{
public:

  MPI_iCommData():
    size( 0 ),
    commID( CommunicationTools::getInstance().getCommID() ),
    sizeCommID( CommunicationTools::getInstance().getCommID() ),
    fieldNames(),
    mpiSendBufferRequest(),
    mpiRecvBufferRequest(),
    mpiSendBufferStatus(),
    mpiRecvBufferStatus()
  {}

  ~MPI_iCommData()
  {}



  void resize( localIndex numMessages )
  {
    mpiSendBufferRequest.resize( numMessages );
    mpiRecvBufferRequest.resize( numMessages );
    mpiSendBufferStatus.resize( numMessages );
    mpiRecvBufferStatus.resize( numMessages );
    mpiSizeSendBufferRequest.resize( numMessages );
    mpiSizeRecvBufferRequest.resize( numMessages );
    mpiSizeSendBufferStatus.resize( numMessages );
    mpiSizeRecvBufferStatus.resize( numMessages );
    size = static_cast< int >(numMessages);
  }

  int size;
  int commID;
  int sizeCommID;
  std::map< string, string_array > fieldNames;

  array1d< MPI_Request > mpiSendBufferRequest;
  array1d< MPI_Request > mpiRecvBufferRequest;
  array1d< MPI_Status >  mpiSendBufferStatus;
  array1d< MPI_Status >  mpiRecvBufferStatus;

  array1d< MPI_Request > mpiSizeSendBufferRequest;
  array1d< MPI_Request > mpiSizeRecvBufferRequest;
  array1d< MPI_Status >  mpiSizeSendBufferStatus;
  array1d< MPI_Status >  mpiSizeRecvBufferStatus;
};



} /* namespace geosx */

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_ */
