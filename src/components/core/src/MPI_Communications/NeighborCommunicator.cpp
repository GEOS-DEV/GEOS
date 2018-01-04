/*
 * NeighborCommunicator.cpp
 *
 *  Created on: Jan 2, 2018
 *      Author: settgast
 */

#include "NeighborCommunicator.hpp"

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
                                             MPI_Comm mpiComm  )
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

} /* namespace geosx */
