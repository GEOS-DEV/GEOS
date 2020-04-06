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
/**
 * @file MpiWrapper.hpp
 */

#ifndef GEOSX_MPICOMMUNICATIONS_MPIWRAPPER_HPP_
#define GEOSX_MPICOMMUNICATIONS_MPIWRAPPER_HPP_

#include "common/DataTypes.hpp"

#if defined(GEOSX_USE_MPI)
  #include <mpi.h>
  #define MPI_PARAM( x ) x
#else
  #define MPI_PARAM( x )
typedef int MPI_Comm;

  #define MPI_COMM_NULL      ((MPI_Comm)0x04000000)
  #define MPI_COMM_WORLD ((MPI_Comm)0x44000000)

typedef int MPI_Datatype;
  #define MPI_CHAR           ((MPI_Datatype)0x4c000101)
  #define MPI_SIGNED_CHAR    ((MPI_Datatype)0x4c000118)
  #define MPI_UNSIGNED_CHAR  ((MPI_Datatype)0x4c000102)
  #define MPI_BYTE           ((MPI_Datatype)0x4c00010d)
  #define MPI_WCHAR          ((MPI_Datatype)0x4c00040e)
  #define MPI_SHORT          ((MPI_Datatype)0x4c000203)
  #define MPI_UNSIGNED_SHORT ((MPI_Datatype)0x4c000204)
  #define MPI_INT            ((MPI_Datatype)0x4c000405)
  #define MPI_UNSIGNED       ((MPI_Datatype)0x4c000406)
  #define MPI_LONG           ((MPI_Datatype)0x4c000807)
  #define MPI_UNSIGNED_LONG  ((MPI_Datatype)0x4c000808)
  #define MPI_FLOAT          ((MPI_Datatype)0x4c00040a)
  #define MPI_DOUBLE         ((MPI_Datatype)0x4c00080b)
  #define MPI_LONG_DOUBLE    ((MPI_Datatype)0x4c00100c)
  #define MPI_LONG_LONG_INT  ((MPI_Datatype)0x4c000809)
  #define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
  #define MPI_LONG_LONG      MPI_LONG_LONG_INT

typedef int MPI_Op;

  #define MPI_MAX     (MPI_Op)(0x58000001)
  #define MPI_MIN     (MPI_Op)(0x58000002)
  #define MPI_SUM     (MPI_Op)(0x58000003)
  #define MPI_PROD    (MPI_Op)(0x58000004)
  #define MPI_LAND    (MPI_Op)(0x58000005)
  #define MPI_BAND    (MPI_Op)(0x58000006)
  #define MPI_LOR     (MPI_Op)(0x58000007)
  #define MPI_BOR     (MPI_Op)(0x58000008)
  #define MPI_LXOR    (MPI_Op)(0x58000009)
  #define MPI_BXOR    (MPI_Op)(0x5800000a)
  #define MPI_MINLOC  (MPI_Op)(0x5800000b)
  #define MPI_MAXLOC  (MPI_Op)(0x5800000c)
  #define MPI_REPLACE (MPI_Op)(0x5800000d)
  #define MPI_NO_OP   (MPI_Op)(0x5800000e)

  #define MPI_SUCCESS          0      /* Successful return code */

typedef int MPI_Request;

typedef int MPI_Info;
  #define MPI_INFO_NULL (MPI_Info)(0x60000000)

struct MPI_Status
{
  int junk;
};

#endif


namespace geosx
{

/**
 * @struct MpiWrapper
 * This struct is a wrapper for all mpi.h functions that are used in GEOSX, and provides a collection of
 * convenience functions to make using the raw mpi functions simpler.
 *
 * The static wrapper functions around the mpi.h function are named by removing the "MPI_" from the beginning
 * of the native mpi function name. For instance the "Comm_rank()" function calls "MPI_Comm_rank()". Since all
 * wrapper functions are static, the should be referred to by their scoped name, for example "MpiWrapper::Comm_rank()".
 */
struct MpiWrapper
{
public:

  /**
   * @enum Reduction
   * Strongly typed enum class for calling collective functions using MPI_Op
   */
  enum class Reduction
  {
    Max,  //!< Max
    Min,  //!< Min
    Sum,  //!< Sum
    Prod //!< Prod
  };

  MpiWrapper() = delete;
  virtual ~MpiWrapper() = delete;

  /**
   * @name FUNCTION GROUP for the direct wrappers around naitive MPI functions
   * @param[in] sendbuf  Pointer to the memory to read the sent data from.
   * @param[out] recvbuf Pointer to the memory to write the received data in.
   * @param[in] count    The number of data entries that are being communicated.
   * @param[in] datatype The MPI_Datatype that is being communicated.
   * @param[in] op       The collective MPI_Op to apply for the function.
   * @param[in] comm     The MPI_Comm communicator that the function will act on.
   *
   * Please see standard MPI documentation for a detailed description of the parameters for each function that
   * is being wrapped
   */
  ///@{

  static void Barrier( MPI_Comm const & MPI_PARAM( comm )=MPI_COMM_GEOSX )
  {
  #ifdef GEOSX_USE_MPI
    MPI_Barrier( comm );
  #endif
  }

//  static int Bcast( void * buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm );

  static int Cart_coords( MPI_Comm comm, int rank, int maxdims, int coords[] );

  static int Cart_create( MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                          int reorder, MPI_Comm * comm_cart );

  static int Cart_rank( MPI_Comm comm, const int coords[] );

  static int Comm_free( MPI_Comm * comm );

  inline static int Comm_rank( MPI_Comm const & MPI_PARAM( comm )=MPI_COMM_GEOSX )
  {
    int rank = 0;
  #ifdef GEOSX_USE_MPI
    MPI_Comm_rank( comm, &rank );
  #endif
    return rank;
  }

  inline static int Comm_size( MPI_Comm const & MPI_PARAM( comm )=MPI_COMM_GEOSX )
  {
    int size = 1;
#ifdef GEOSX_USE_MPI
    MPI_Comm_size( comm, &size );
#endif
    return size;
  }

  static int Finalize( void );


  static int Init( int * argc, char * * * argv );

  static int Test( MPI_Request * request, int * flag, MPI_Status * status );

  static int Wait( MPI_Request * request, MPI_Status * status );

  static int Waitany( int count, MPI_Request array_of_requests[], int * indx, MPI_Status * status );

  static int Waitsome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] );

  static int Waitall( int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[] );

  static double Wtime( void );


  /**
   * Wait on MPI_Requests to complete on at a time and trigger a callback to
   *  process the completion.
   * @param[in] count The number of MPI_Requests being processed.
   * @param[inout] array_of_requests The MPI_Requests to actively wait on.
   * @param[in] func A callable object accepting an integer denoting the MPI_Request index
   *              which has completed.
   * @return MPI_SUCCESS or an MPI_ERROR returned by internal calls to MPI_WaitAny.
   */
  static int ActiveWaitAny( const int count, MPI_Request array_of_requests[], std::function< void ( int ) > func );

  /**
   * Wait on MPI_Requests to complete on or more at a time and trigger a callback to
   *  process the completion.
   * @param[in] count The number of MPI_Requests being processed.
   * @param[inout] array_of_requests The MPI_Requests to actively wait on.
   * @param[in] func A callable object accepting an integer denoting the MPI_Request index
   *              which has completed.
   * @return MPI_SUCCESS or an MPI_ERROR returned by internal calls to MPI_WaitSome.
   */
  static int ActiveWaitSome( const int count, MPI_Request array_of_requests[], std::function< void ( int ) > func );

  /**
   * Active non-blocking phased communication with multiple participants,
   *  each participant in each phase cannot depend on the previous phases
   *  being complete for any participant other than itself.
   * @param[in] participants The number of participants in each phase
   * @param[in] phases A vector of function objects taking int and returning MPI_Requests
   *               denoting the state of that participant in that phase.
   * @note The only restriction on phase[N](index) being called is that phase[N-1](index)
   *       has been called and the MPI_Request returned by that call has completed.
   * @note One can add a final recv phase by having that phase return MPI_REQUEST_NULL.
   * @return MPI_SUCCESS or and MPI_ERROR from internal calls to MPI_WaitAny.
   */
  static int ActiveWaitSomePartialPhase( const int participants,
                                         std::vector< std::function< MPI_Request ( int ) > > & phases );

  /**
   * Active non-blocking phased communication with multiple participants,
   *  each participant in each phase may depend on the previous phases
   *  being fully complete prior to entry into a subsequent phase.
   * @param[in] participants The number of participants in each phase
   * @param[in] phases A vector of function objects taking int and returning MPI_Requests
   *               denoting the state of that participant in that phase.
   * @note The restriction on phase[N](index) being called is that phase[N-1](0 - (particpants-1))
   *       have all been called and the MPI_Requests from those calls have all been completed.
   * @note One can add a final recv phase by having that phase return MPI_REQUEST_NULL.
   * @return MPI_SUCCESS or and MPI_ERROR from internal calls to MPI_WaitAny.
   */
  static int ActiveWaitSomeCompletePhase( const int participants,
                                          std::vector< std::function< MPI_Request ( int ) > > & phases );

  /**
   * Active blocking phased communication with multiple participants,
   *  each participant in each phase may depend on the previous phases
   *  being fully complete prior to entry into a subsequent phase.
   * @param[in] participants The number of participants in each phase
   * @param[in] phases A vector of function objects taking int and returning MPI_Requests
   *               denoting the state of that participant in that phase.
   * @note The restriction on phase[N](index) being called is that phase[N-1](0 - (particpants-1))
   *       and phase[N](0..index) have all been called and the MPI_Requests from those calls have all
   *       been completed.
   * @note One can add a final recv phase by having that phase return MPI_REQUEST_NULL.
   * @return MPI_SUCCESS or and MPI_ERROR from internal calls to MPI_WaitAny.
   */
  static int ActiveWaitOrderedCompletePhase( const int participants,
                                             std::vector< std::function< MPI_Request ( int ) > > & phases );
  ///@}

#if !defined(GEOSX_USE_MPI)
  static std::map< int, std::pair< int, void * > > & getTagToPointersMap()
  {
    static std::map< int, std::pair< int, void * > > tagToPointers;
    return tagToPointers;
  }
#endif

  /**
   * @brief Strongly typed wrapper around MPI_Allgather.
   * @tparam T_SEND The pointer type for \p sendbuf
   * @tparam T_RECV The pointer type for \p recvbuf
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[in] sendcount The number of values to send.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] recvcount The number of values to receive.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return The return value of the underlying call to MPI_Allgather().
   */
  template< typename T_SEND, typename T_RECV >
  static int Allgather( T_SEND const * sendbuf,
                        int sendcount,
                        T_RECV * recvbuf,
                        int recvcount,
                        MPI_Comm comm );

  /**
   * @brief Convenience function for MPI_Allgather.
   * @tparam T The type to send/recieve. This must have a valid conversion to MPI_Datatype in getMpiType();
   * @param[in] myValue The value to send.
   * @param[out] allValues The values recived from each rank.
   */
  template< typename T >
  static void allGather( T const myValue, array1d< T > & allValues );

  template< typename T >
  static int allGather( arrayView1d< T const > const & sendbuf,
                        array1d< T > & recvbuf );

  /**
   * @brief Strongly typed wrapper around MPI_Allreduce.
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] count The number of values to send/receive.
   * @param[in] op The MPI_Op to perform.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return The return value of the underlying call to MPI_Allreduce().
   */
  template< typename T >
  static int allReduce( T const * sendbuf, T * recvbuf, int count, MPI_Op op, MPI_Comm comm );


  template< typename T >
  static int scan( T const * sendbuf, T * recvbuf, int count, MPI_Op op, MPI_Comm comm);

  template< typename T >
  static int exscan( T const * sendbuf, T * recvbuf, int count, MPI_Op op, MPI_Comm comm);

  /**
   * @brief Strongly typed wrapper around MPI_Bcast.
   * @param[in/out] buffer The pointer to the send/recv buffer
   * @param[in] count The number of data types to send.
   * @param[in] root The rank sending the data.
   * @param[in] comm The MPI_Comm over which the MPI_Bcast operates.
   * @return The return value of the underlying call to MPI_Bcast().
   */
  template< typename T >
  static int bcast( T * buffer, int count, int root, MPI_Comm comm );


  /**
   * @brief Convenience function for MPI_Broadcast.
   * @tparam T The type to send/recieve. This must have a valid conversion to MPI_Type in getMpiType();
   * @param[in/out] myValue The value to send from the \p srcRank to the other ranks.
   * @param srcRank The rank that is sending the \p value.
   */
  template< typename T >
  static void Broadcast( T & value, int srcRank = 0 );

  /**
   * @brief Strongly typed wrapper around MPI_Gather().
   * @tparam TS The pointer type for \p sendbuf
   * @tparam TR The pointer type for \p recvbuf
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[in] sendcount The number of values to send.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] recvcount The number of values to receive.
   * @param[in] root The rank recieving the data.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return
   */
  template< typename TS, typename TR >
  static int gather( TS const * const sendbuf,
                     int sendcount,
                     TR * const recvbuf,
                     int recvcount,
                     int root,
                     MPI_Comm comm );

  /**
   * @brief Strongly typed wrapper around MPI_Gatherv.
   * @tparam TS The pointer type for \p sendbuf
   * @tparam TR The pointer type for \p recvbuf
   * @param[in] sendbuf The pointer to the sending buffer.
   * @param[in] sendcount The number of values to send.
   * @param[out] recvbuf The pointer to the receive buffer.
   * @param[in] recvcount The number of values to receive.
   * @param[in] displs integer array (of length group size). Entry i specifies the displacement relative to recvbuf at
   *                   which to place the incoming data from process i (significant only at root).
   * @param[in] root The rank recieving the data.
   * @param[in] comm The MPI_Comm over which the gather operates.
   * @return
   */
  template< typename TS, typename TR >
  static int gatherv( TS const * const sendbuf,
                      int sendcount,
                      TR * const recvbuf,
                      const int * recvcounts,
                      const int * displs,
                      int root,
                      MPI_Comm comm );


  /**
   * @brief Returns an MPI_Datatype from a c type.
   * @tparam T The type for which we want an MPI_Datatype
   * @return The MPI_Datatype associated wtih \param T
   */
  template< typename T >
  static MPI_Datatype getMpiType();

  static std::size_t getSizeofMpiType( MPI_Datatype const type );

  /**
   * @brief Returns an MPI_Op associated with our strongly typed Reduction enum.
   * @param[in] op The value of the Reduction enum to get an MPI_Op for.
   * @return The MPI_Op associated with \p op.
   */
  static MPI_Op getMpiOp( Reduction const op );


  /**
   * @brief Strongly typed wrapper around MPI_Irecv()
   * @param[out] buf The pointer to the buffer that contains the data to be received.
   * @param[in] count The number of elements in \p buf
   * @param[in] source The rank of the source process within \p comm.
   * @param[in] tag The message tag that is be used to distinguish different types of messages
   * @param[in] comm The handle to the MPI_Comm
   * @param[out] request Pointer to the MPI_Request associated with this request.
   * @return
   */
  template< typename T >
  static int iRecv( T * const buf,
                    int count,
                    int source,
                    int tag,
                    MPI_Comm comm,
                    MPI_Request * request );

  /**
   * @brief Strongly typed wrapper around MPI_Isend()
   * @param[in] buf The pointer to the buffer that contains the data to be sent.
   * @param[in] count The number of elements in \p buf.
   * @param[in] dest The rank of the destination process within \p comm.
   * @param[in] tag The message tag that is be used to distinguish different types of messages.
   * @param[in] comm The handle to the MPI_Comm.
   * @param[out] request Pointer to the MPI_Request associated with this request.
   * @return
   */
  template< typename T >
  static int iSend( T const * const buf,
                    int count,
                    int dest,
                    int tag,
                    MPI_Comm comm,
                    MPI_Request * request );

  /**
   * @brief Convenience function for a MPI_Reduce using a MPI_MIN operation.
   * @param value the value to send into the reduction.
   * @return The minimum of all \p value across the ranks.
   */
  template< typename T >
  static T Min( T const & value );

  /**
   * @brief Convenience function for a MPI_Reduce using a MPI_MAX operation.
   * @param[in] value the value to send into the reduction.
   * @return The maximum of all \p value across the ranks.
   */
  template< typename T >
  static T Max( T const & value );

  /**
   * @brief Compute exclusive prefix sum and full sum
   * @tparam T type of local (rank) value
   * @tparam U type of global (sum) value
   * @param[in] value the local value
   * @return a pair where first is the prefix sum, second is the full sum
   */
  template< typename U, typename T >
  static std::pair< U, U > PrefixSum( T const & value );

  /**
   * @brief Convenience function for the MPI_Reduce function.
   * @tparam The type of data to reduce. Must be a valid MPI_Datatype.
   * @param value The value to send to the reduction.
   * @param op The Reduction enum to perform.
   * @return The value of reduction across all ranks
   */
  template< typename T >
  static T Reduce( T const & value, Reduction const op );

  /**
   * @brief Convenience function for a MPI_Reduce using a MPI_SUM operation.
   * @param[in] value the value to send into the reduction.
   * @return The sum of all \p value across the ranks.
   */
  template< typename T >
  static T Sum( T const & value );
};

template<> inline MPI_Datatype MpiWrapper::getMpiType< char >()           { return MPI_CHAR; }
template<> inline MPI_Datatype MpiWrapper::getMpiType< signed char >()    { return MPI_SIGNED_CHAR; }
template<> inline MPI_Datatype MpiWrapper::getMpiType< float >()          { return MPI_FLOAT; }
template<> inline MPI_Datatype MpiWrapper::getMpiType< double >()         { return MPI_DOUBLE; }
template<> inline MPI_Datatype MpiWrapper::getMpiType< int >()            { return MPI_INT; }
template<> inline MPI_Datatype MpiWrapper::getMpiType< long int >()       { return MPI_LONG; }
template<> inline MPI_Datatype MpiWrapper::getMpiType< long long int >()  { return MPI_LONG_LONG; }

inline MPI_Op MpiWrapper::getMpiOp( Reduction const op )
{
  switch( op )
  {
    case Reduction::Sum:
    {
      return MPI_SUM;
    }
    case Reduction::Min:
    {
      return MPI_MIN;
    }
    case Reduction::Max:
    {
      return MPI_MAX;
    }
    case Reduction::Prod:
    {
      return MPI_PROD;
    }
    default:
      GEOSX_ERROR( "Unsupported reduction operation" );
      return MPI_NO_OP;
  }
}

template< typename T_SEND, typename T_RECV >
int MpiWrapper::Allgather( T_SEND const * const sendbuf,
                           int sendcount,
                           T_RECV * const recvbuf,
                           int recvcount,
                           MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  return MPI_Allgather( sendbuf, sendcount, getMpiType< T_SEND >(), recvbuf, recvcount, getMpiType< T_RECV >(), comm );
#else
  static_assert( std::is_same< T_SEND, T_RECV >::value,
                 "MpiWrapper::Allgather() for serial run requires send and receive buffers are of the same type" );
  GEOSX_ERROR_IF_NE_MSG( sendcount, recvcount, "sendcount is not equal to recvcount." );
  *recvbuf = *sendbuf;
  return 0;
#endif
}


template< typename T >
void MpiWrapper::allGather( T const myValue, array1d< T > & allValues )
{
#ifdef GEOSX_USE_MPI
  int const mpiSize = Comm_size( MPI_COMM_GEOSX );
  allValues.resize( mpiSize );

  MPI_Datatype const MPI_TYPE = getMpiType< T >();

  MPI_Allgather( &myValue, 1, MPI_TYPE, allValues.data(), 1, MPI_TYPE, MPI_COMM_GEOSX );

#else
  allValues.resize( 1 );
  allValues[0] = myValue;
#endif
}

template< typename T >
int MpiWrapper::allGather( arrayView1d< T const > const & sendValues,
                           array1d< T > & allValues )
{
  int const sendSize = integer_conversion< int >( sendValues.size() );
#ifdef GEOSX_USE_MPI
  int const mpiSize = Comm_size( MPI_COMM_GEOSX );
  allValues.resize( mpiSize * sendSize );
  return MPI_Allgather( sendValues.data(),
                        sendSize,
                        getMpiType< T >(),
                        allValues.data(),
                        sendSize,
                        getMpiType< T >(),
                        MPI_COMM_GEOSX );

#else
  allValues.resize( sendSize );
  for( localIndex a=0 ; a<sendSize ; ++a )
  {
    allValues[a] = sendValues[a];
  }
  return 0;
#endif
}

template< typename T >
int MpiWrapper::allReduce( T const * const sendbuf,
                           T * const recvbuf,
                           int count,
                           MPI_Op MPI_PARAM( op ),
                           MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  MPI_Datatype const MPI_TYPE = getMpiType< T >();
  return MPI_Allreduce( sendbuf, recvbuf, count, MPI_TYPE, op, comm );
#else
  memcpy( recvbuf, sendbuf, count*sizeof(T) );
  return 0;
#endif
}

template< typename T >
int MpiWrapper::scan( T const * const sendbuf,
                      T * const recvbuf,
                      int count,
                      MPI_Op MPI_PARAM( op ),
                      MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  MPI_Datatype const MPI_TYPE = getMpiType< T >();
  return MPI_Scan( sendbuf, recvbuf, count, MPI_TYPE, op, comm );
#else
  memcpy( recvbuf, sendbuf, count*sizeof(T) );
  return 0;
#endif
}

template< typename T >
int MpiWrapper::exscan( T const * const MPI_PARAM( sendbuf ),
                      T * const recvbuf,
                      int count,
                      MPI_Op MPI_PARAM( op ),
                      MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  MPI_Datatype const MPI_TYPE = getMpiType< T >();
  return MPI_Exscan( sendbuf, recvbuf, count, MPI_TYPE, op, comm );
#else
  memset( recvbuf, 0, count*sizeof(T) );
  return 0;
#endif
}

template< typename T >
int MpiWrapper::bcast( T * const MPI_PARAM( buffer ),
                       int MPI_PARAM( count ),
                       int MPI_PARAM( root ),
                       MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  return MPI_Bcast( buffer, count, getMpiType< T >(), root, comm );
#else
  return 0;
#endif

}

template< typename T >
void MpiWrapper::Broadcast( T & MPI_PARAM( value ), int MPI_PARAM( srcRank ) )
{
#ifdef GEOSX_USE_MPI
  MPI_Datatype const mpiType = getMpiType< T >();
  MPI_Bcast( &value, 1, mpiType, srcRank, MPI_COMM_GEOSX );
#endif
}

template<>
inline
void MpiWrapper::Broadcast< std::string >( std::string & MPI_PARAM( value ), int MPI_PARAM( srcRank ) )
{
#ifdef GEOSX_USE_MPI
  int size = value.size();
  Broadcast( size, srcRank );

  value.resize( size );

  MPI_Bcast( const_cast< char * >( value.data() ), size, getMpiType< char >(), srcRank, MPI_COMM_GEOSX );
#endif
}

template< typename TS, typename TR >
int MpiWrapper::gather( TS const * const sendbuf,
                        int sendcount,
                        TR * const recvbuf,
                        int recvcount,
                        int MPI_PARAM( root ),
                        MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  return MPI_Gather( sendbuf, sendcount, getMpiType< TS >(), recvbuf, recvcount, getMpiType< TR >(), root, comm );
#else
  static_assert( std::is_same< TS, TR >::value,
                 "MpiWrapper::gather() for serial run requires send and receive buffers are of the same type" );
  std::size_t const sendBufferSize = sendcount * sizeof(TS);
  std::size_t const recvBufferSize = recvcount * sizeof(TR);
  GEOSX_ERROR_IF_NE_MSG( sendBufferSize, recvBufferSize, "size of send buffer and receive buffer are not equal" );
  memcpy( recvbuf, sendbuf, sendBufferSize );
  return 0;
#endif
}

template< typename TS, typename TR >
int MpiWrapper::gatherv( TS const * const sendbuf,
                         int sendcount,
                         TR * const recvbuf,
                         const int * recvcounts,
                         const int * MPI_PARAM( displs ),
                         int MPI_PARAM( root ),
                         MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  return MPI_Gatherv( sendbuf, sendcount, getMpiType< TS >(), recvbuf, recvcounts, displs, getMpiType< TR >(), root, comm );
#else
  static_assert( std::is_same< TS, TR >::value,
                 "MpiWrapper::gather() for serial run requires send and receive buffers are of the same type" );
  std::size_t const sendBufferSize = sendcount * sizeof(TS);
  std::size_t const recvBufferSize = recvcounts[0] * sizeof(TR);
  GEOSX_ERROR_IF_NE_MSG( sendBufferSize, recvBufferSize, "size of send buffer and receive buffer are not equal" );
  memcpy( recvbuf, sendbuf, sendBufferSize );
  return 0;
#endif
}

template< typename T >
int MpiWrapper::iRecv( T * const buf,
                       int count,
                       int MPI_PARAM( source ),
                       int tag,
                       MPI_Comm MPI_PARAM( comm ),
                       MPI_Request * MPI_PARAM( request ) )
{
#ifdef GEOSX_USE_MPI
  return MPI_Irecv( buf, count, getMpiType< T >(), source, tag, comm, request );
#else
  std::map< int, std::pair< int, void * > > & pointerMap = getTagToPointersMap();
  std::map< int, std::pair< int, void * > >::iterator iPointer = pointerMap.find( tag );

  if( iPointer==pointerMap.end() )
  {
    pointerMap.insert( {tag, {1, buf}
                       } );
  }
  else
  {
    GEOSX_ERROR_IF( iPointer->second.first != 0,
                    "Tag does is assigned, but pointer was not set by iSend." );
    memcpy( buf, iPointer->second.second, count*sizeof(T) );
    pointerMap.erase( iPointer );
  }
  return 0;
#endif
}

template< typename T >
int MpiWrapper::iSend( T const * const buf,
                       int count,
                       int MPI_PARAM( dest ),
                       int tag,
                       MPI_Comm MPI_PARAM( comm ),
                       MPI_Request * MPI_PARAM( request ) )
{
#ifdef GEOSX_USE_MPI
  return MPI_Isend( buf, count, getMpiType< T >(), dest, tag, comm, request );
#else
  std::map< int, std::pair< int, void * > > & pointerMap = getTagToPointersMap();
  std::map< int, std::pair< int, void * > >::iterator iPointer = pointerMap.find( tag );

  if( iPointer==pointerMap.end() )
  {
    pointerMap.insert( {tag, {0, const_cast< T * >(buf)}
                       } );
  }
  else
  {
    GEOSX_ERROR_IF( iPointer->second.first != 1,
                    "Tag does is assigned, but pointer was not set by iRecv." );
    memcpy( iPointer->second.second, buf, count*sizeof(T) );
    pointerMap.erase( iPointer );
  }
  return 0;
#endif
}

template< typename U, typename T >
std::pair< U, U > MpiWrapper::PrefixSum( T const & value )
{
  array1d< T > gather;
  allGather( value, gather );

  std::pair< U, U > rval( 0, 0 );

  int const mpiSize = Comm_size( MPI_COMM_GEOSX );
  int const mpiRank = Comm_rank( MPI_COMM_GEOSX );

  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    if( p < mpiRank )
    {
      rval.first += static_cast< U >( gather[p] );
    }
    rval.second += static_cast< U >( gather[p] );
  }

  return rval;
}


template< typename T >
T MpiWrapper::Reduce( T const & value, Reduction const MPI_PARAM( op ) )
{
  T result = value;
#ifdef GEOSX_USE_MPI
  MPI_Allreduce( &value, &result, 1, getMpiType< T >(), getMpiOp( op ), MPI_COMM_GEOSX );
#else
#endif
  return result;
}

template< typename T >
T MpiWrapper::Sum( T const & value )
{
  return MpiWrapper::Reduce( value, Reduction::Sum );
}

template< typename T >
T MpiWrapper::Min( T const & value )
{
  return MpiWrapper::Reduce( value, Reduction::Min );
}

template< typename T >
T MpiWrapper::Max( T const & value )
{
  return MpiWrapper::Reduce( value, Reduction::Max );
}



} /* namespace geosx */

#endif /* GEOSX_MPICOMMUNICATIONS_MPIWRAPPER_HPP_ */
