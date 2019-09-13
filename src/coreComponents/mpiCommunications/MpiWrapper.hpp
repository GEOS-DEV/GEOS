/*
 * MpiWrapper.hpp
 *
 *  Created on: Sep 11, 2019
 *      Author: settgast
 */

#ifndef GEOSX_MPICOMMUNICATIONS_MPIWRAPPER_HPP_
#define GEOSX_MPICOMMUNICATIONS_MPIWRAPPER_HPP_

#include "common/DataTypes.hpp"

#if defined(GEOSX_USE_MPI)
  #include <mpi.h>
  #define MPI_PARAM(x) x
#else
  #define MPI_PARAM(x)
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

  struct MPI_Status
  {
    int junk;
  };



#endif


namespace geosx
{

class MpiWrapper
{
public:

  enum class Reduction { Sum, Min, Max };

  MpiWrapper();
  virtual ~MpiWrapper();

  template< typename T >
  static MPI_Datatype getMpiType();

  static MPI_Op getMpiOp( Reduction const op );

  static int Cart_coords( MPI_Comm comm, int rank, int maxdims, int coords[] );

  static int Cart_create( MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                          int reorder, MPI_Comm *comm_cart );


  static int Comm_free( MPI_Comm *comm );

  static int Init( int *argc, char ***argv );

  static int Finalize( void );


  static int Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                       MPI_Op op, MPI_Comm comm);

  static int Allgather( const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                        int recvcount, MPI_Datatype recvtype, MPI_Comm comm);

  static void Barrier( MPI_Comm ) {};

  static int Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

  static int MPI_Size( MPI_Comm const & comm );
  static int MPI_Rank( MPI_Comm const & comm );
  static int Cart_rank(MPI_Comm comm, const int coords[] );

  template<typename T>
  static void allGather( T const myValue, array1d<T> & allValues );

  static int Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);

  static int Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root,
                  MPI_Comm comm);


  static int Isend( const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                        MPI_Comm comm, MPI_Request *request);
  static int Irecv( void *buf, int count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, MPI_Request *request);

  /**
   * @brief Compute exclusive prefix sum and full sum
   * @tparam T type of local (rank) value
   * @tparam U type of global (sum) value
   * @param value the local value
   * @return a pair where first is the prefix sum, second is the full sum
   */
  template< typename U, typename T >
  static std::pair<U, U> PrefixSum( T const & value );

  template<typename T>
  static void Broadcast( T & value, int srcRank = 0 );

  template< typename T >
  static T Reduce( T const & value, Reduction const op);

  template<typename T>
  static T Sum( T const & value );

  template<typename T>
  static T Min( T const & value );

  template<typename T>
  static T Max( T const & value );

  static int Wait(MPI_Request *request, MPI_Status *status);
  static int Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status);
  static int Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);

  static double Wtime(void);

};


template<> inline MPI_Datatype MpiWrapper::getMpiType<double>()         { return MPI_DOUBLE; }
template<> inline MPI_Datatype MpiWrapper::getMpiType<int>()            { return MPI_INT; }
template<> inline MPI_Datatype MpiWrapper::getMpiType<long int>()       { return MPI_LONG; }
template<> inline MPI_Datatype MpiWrapper::getMpiType<long long int>()  { return MPI_LONG_LONG; }

inline MPI_Op MpiWrapper::getMpiOp( Reduction const op )
{
  switch( op )
  {
    case Reduction::Sum:
      return MPI_SUM;
    case Reduction::Min:
      return MPI_MIN;
    case Reduction::Max:
      return MPI_MAX;
    default:
      GEOS_ERROR( "Unsupported reduction op" );
      return MPI_SUM;
  }
}

template<typename T>
void MpiWrapper::allGather( T const myValue, array1d<T> & allValues )
{
#ifdef GEOSX_USE_MPI
  int const mpiSize = MPI_Size( MPI_COMM_GEOSX );
  allValues.resize( mpiSize );

  MPI_Datatype const MPI_TYPE = getMpiType<T>();

  MPI_Allgather( &myValue, 1, MPI_TYPE, allValues.data(), 1, MPI_TYPE, MPI_COMM_GEOSX );

#else
  allValues.resize(1);
  allValues[0] = myValue;
#endif
}

template< typename U, typename T >
std::pair<U, U> MpiWrapper::PrefixSum( T const & value )
{
  array1d<T> gather;
  allGather( value, gather );

  std::pair<U, U> rval( 0, 0 );

  int const mpiSize = MPI_Size( MPI_COMM_GEOSX );
  int const mpiRank = MPI_Rank( MPI_COMM_GEOSX );

  for( localIndex p = 0 ; p < mpiSize ; ++p )
  {
    if( p < mpiRank )
    {
      rval.first += static_cast<U>( gather[p] );
    }
    rval.second += static_cast<U>( gather[p] );
  }

  return rval;
}

template<typename T>
void MpiWrapper::Broadcast( T & MPI_PARAM(value), int MPI_PARAM(srcRank) )
{
#ifdef GEOSX_USE_MPI
  MPI_Datatype const mpiType = getMpiType<T>();
  MPI_Bcast( &value, 1, mpiType, srcRank, MPI_COMM_GEOSX );
#endif
}

template< typename T >
T MpiWrapper::Reduce( T const & value, Reduction const MPI_PARAM(op) )
{
  T result = value;
#ifdef GEOSX_USE_MPI
  MPI_Allreduce( &value, &result, 1, getMpiType<T>(), getMpiOp( op ), MPI_COMM_GEOSX );
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
