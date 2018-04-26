/*
 * NeighborCommunicator.hpp
 *
 *  Created on: Jan 2, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_
#define SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_

#include<set>
#include<mpi.h>
#include<vector>
#include "common/DataTypes.hpp"

namespace geosx
{
inline int CommTag( int const senderRank, int const receiverRank, int const comm )
{
  int m_size;
  MPI_Comm_size( MPI_COMM_WORLD, &m_size );
  return senderRank * m_size + receiverRank + m_size * m_size * comm;
}

class MeshLevel;


class NeighborCommunicator
{
public:

  NeighborCommunicator();
//  ~NeighborCommunicator();


  void MPI_iSendReceive( char const * const sendBuffer,
                         int const sendSize,
                         MPI_Request& sendRequest,
                         char * const receiveBuffer,
                         int const receiveSize,
                         MPI_Request& receiveRequest,
                         int const commID,
                         MPI_Comm mpiComm );


  void MPI_iSendReceiveBufferSizes( int const commID,
                                    MPI_Comm mpiComm  );

  void MPI_iSendReceiveBufferSizes( int const commID,
                                    MPI_Request& mpiSendRequest,
                                    MPI_Request& mpiRecvRequest,
                                    MPI_Comm mpiComm  );

  void MPI_iSendReceiveBuffers( int const commID,
                                MPI_Comm mpiComm  );

  void MPI_iSendReceiveBuffers( int const commID,
                                MPI_Request& mpiSendRequest,
                                MPI_Request& mpiRecvRequest,
                                MPI_Comm mpiComm  );

  void MPI_iSendReceive( int const commID,
                         MPI_Comm mpiComm  );

  void MPI_iSendReceive( int const commID,
                         MPI_Request& mpiSendRequest,
                         MPI_Request& mpiRecvRequest,
                         MPI_Comm mpiComm  );

  template<typename T>
  void MPI_iSendReceive( T const * const sendBuffer,
                         int const sendSize,
                         MPI_Request& sendReq,
                         T * const recvBuffer,
                         int const recvSize,
                         MPI_Request& recvReq,
                         int const commID,
                         MPI_Comm mpiComm  )
  {
    MPI_iSendReceive( reinterpret_cast<char const*>( sendBuffer ),
                      sendSize * sizeof(T),
                      sendReq,
                      reinterpret_cast<char*>( recvBuffer ),
                      recvSize * sizeof(T),
                      recvReq,
                      commID,
                      mpiComm );
  }


  template<typename T>
  void MPI_iSendReceive( array<T> const & sendBuffer,
                         MPI_Request& sendReq,
                         array<T> & recvBuffer,
                         MPI_Request& recvReq,
                         int const commID,
                         MPI_Comm mpiComm  );


  template<typename T>
  void MPI_iSendReceive( array<T> const & sendBuffer,
                         array<T> & recvBuffer,
                         int const commID,
                         MPI_Comm mpiComm )
  {
    MPI_iSendReceive( sendBuffer,
                      m_mpiSendBufferRequest[commID],
                      recvBuffer,
                      m_mpiRecvBufferRequest[commID],
                      commID,
                      mpiComm  );
  }

  template<typename T>
  void MPI_iSendReceiveWait( array<T> const & sendBuffer,
                             array<T> & recvBuffer,
                             int const commID,
                             MPI_Comm mpiComm  )
  {
    MPI_iSendReceive( sendBuffer, recvBuffer, commID, mpiComm );
    MPI_WaitAll(commID);
  }


  void MPI_iSendReceive( char const * const sendBuffer,
                         int const sendSize,
                         int const commID,
                         MPI_Comm mpiComm  );

  void MPI_WaitAll( int const commID,
                    MPI_Request& mpiSendRequest,
                    MPI_Status& mpiSendStatus,
                    MPI_Request& mpiRecvRequest,
                    MPI_Status& mpiReceiveStatus );

  void MPI_WaitAll( int const commID );

  void FindAndPackGhosts( bool const contactActive,
                          integer const depth,
                          MeshLevel * const meshLevel,
                          int const commID );

  void UnpackGhosts( MeshLevel * const meshLevel,
                     int const commID );

  void PackBufferForSync( MeshLevel * const meshLevel,
                          int const commID );

  int Rank();

  void SetNeighborRank( int const rank ) { m_neighborRank = rank; }
  int NeighborRank() const { return m_neighborRank; }

  void Clear();

  static int constexpr maxComm = 100;

  buffer_type const & ReceiveBuffer( int commID ) const
  {
    return m_receiveBuffer[commID];
  }
  buffer_type & ReceiveBuffer( int commID )
  {
    return m_receiveBuffer[commID];
  }

  buffer_type const & SendBuffer( int commID ) const
  {
    return m_sendBuffer[commID];
  }
  buffer_type & SendBuffer( int commID )
  {
    return m_sendBuffer[commID];
  }

  void AddNeighborGroupToMesh( MeshLevel * const mesh ) const;

private:
  int m_neighborRank;

  int m_sendBufferSize[maxComm];
  int m_receiveBufferSize[maxComm];

  buffer_type m_sendBuffer[maxComm];
  buffer_type m_receiveBuffer[maxComm];

  MPI_Request m_mpiSendBufferRequest[maxComm];
  MPI_Request m_mpiRecvBufferRequest[maxComm];
  MPI_Status m_mpiSendBufferStatus[maxComm];
  MPI_Status m_mpiRecvBufferStatus[maxComm];

};




template<typename T>
void NeighborCommunicator::MPI_iSendReceive( array<T> const & sendBuffer,
                                             MPI_Request& sendReq,
                                             array<T> & recvBuffer,
                                             MPI_Request& recvReq,
                                             int const commID,
                                             MPI_Comm mpiComm  )
{
  m_sendBufferSize[commID] = sendBuffer.size();

  MPI_iSendReceive( &m_sendBufferSize[commID],
                    1,
                    sendReq,
                    &m_receiveBufferSize[commID],
                    1,
                    recvReq,
                    commID,
                    mpiComm );

  MPI_Waitall( 1, &( recvReq ), &( m_mpiRecvBufferStatus[commID] ) );
  MPI_Waitall( 1, &( sendReq ), &( m_mpiSendBufferStatus[commID] ) );

  recvBuffer.resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( sendBuffer.data(),
                    m_sendBufferSize[commID],
                    sendReq,
                    recvBuffer.data(),
                    m_receiveBufferSize[commID],
                    recvReq,
                    commID,
                    mpiComm );
}


} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_ */
