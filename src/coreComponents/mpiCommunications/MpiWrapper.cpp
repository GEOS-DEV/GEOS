/*
 * MpiWrapper.cpp
 *
 *  Created on: Sep 11, 2019
 *      Author: settgast
 */

#include "MpiWrapper.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif


namespace geosx
{


MpiWrapper::MpiWrapper()
{
  // TODO Auto-generated constructor stub

}

MpiWrapper::~MpiWrapper()
{
  // TODO Auto-generated destructor stub
}


int MpiWrapper::Allgather( const void * sendbuf,
                           int sendcout,
                           MPI_Datatype sendtype,
                           void * recvbuf,
                           int recvcount,
                           MPI_Datatype recvtype,
                           MPI_Comm comm )
{
#ifdef GEOSX_USE_MPI
  return MPI_Allgather( sendbuf, sendcount, sendtype, recvbuf,
                        recvcount, recvtype, comm);
#else
  return 0;
#endif
}


int MpiWrapper::Allreduce( const void * sendbuf,
                           void *recvbuf,
                           int count,
                           MPI_Datatype datatype,
                           MPI_Op op,
                           MPI_Comm comm )
{
#ifdef GEOSX_USE_MPI
  return MPI_Allreduce( sendbuf,
                        recvbuf,
                        count,
                        datatype,
                        op,
                        comm );
#else
  return 0;
#endif
}




int MpiWrapper::Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
#ifdef GEOSX_USE_MPI
  return MPI_Bcast( buffer, count, datatype, root, comm)
#else
  return 0;
#endif

}



int MpiWrapper::Cart_coords( MPI_Comm comm, int rank, int maxdims, int coords[] )
{
#ifdef GEOSX_USE_MPI
  MPI_Cart_coords( comm, rank, maxdims, coords );
#else
  return 0;
#endif
}

int MpiWrapper::Cart_create( MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                        int reorder, MPI_Comm *comm_cart )
{
#ifdef GEOSX_USE_MPI
  MPI_Cart_create( comm_old, ndims, dims, periods, reorder, comm_cart );
#else
  return 0;
#endif
}


int MpiWrapper::Comm_free( MPI_Comm *comm )
{
#ifdef GEOSX_USE_MPI
  MPI_Comm_free( comm );
#else
  return 0;
#endif
}

int MpiWrapper::Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                       int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
{
#ifdef GEOSX_USE_MPI
  return MPI_Gather( sendbuf, sendcount, sendtype, recvbuf,
                     recvcount,  recvtype,  root,  comm );
#else
  return 0;
#endif

}

int MpiWrapper::Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root,
                MPI_Comm comm)
{
#ifdef GEOSX_USE_MPI
  return MPI_Gatherv( sendbuf, sendcount, sendtype, recvbuf,
                      recvcounts, displs, recvtype, root,
                      comm);
#else
  return 0;
#endif

}


int MpiWrapper::Init( int *argc, char ***argv )
{
#ifdef GEOSX_USE_MPI
  return MPI_Init( argc, argv );
#else
  return 0;
#endif

}

int MpiWrapper::Isend( const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                       MPI_Comm comm, MPI_Request *request)
{
#ifdef GEOSX_USE_MPI
  return MPI_Isend( buf, count, datatype, dest, tag, comm, request)
#else
  return 0;
#endif

}

int MpiWrapper::Irecv( void *buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, MPI_Request *request)
{
#ifdef GEOSX_USE_MPI
  return MPI_Irecv( buf, count, datatype, source, tag, comm, request)
#else
  return 0;
#endif

}


int MpiWrapper::Finalize()
{
#ifdef GEOSX_USE_MPI
  return MPI_Finalize();
#else
  return 0;
#endif
}


int MpiWrapper::MPI_Size( MPI_Comm const & comm )
{
  int size = 1;
#ifdef GEOSX_USE_MPI
  MPI_Comm_size( comm, &size );
#endif
  return size;
}

int MpiWrapper::MPI_Rank( MPI_Comm const & comm )
{
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( comm, &rank );
#endif
  return rank;
}

int MpiWrapper::Cart_rank( MPI_Comm comm, const int coords[] )
{
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Cart_rank( comm, coords, &rank );
#endif
  return rank;
}

int MpiWrapper::Wait(MPI_Request *request, MPI_Status *status)
{
#ifdef GEOSX_USE_MPI
  return MPI_Wait(request, status);
#endif
  return 0;
}

int MpiWrapper::Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status)
{
#ifdef GEOSX_USE_MPI
  return MPI_Waitany( count,  array_of_requests, indx, status);
#endif
  return 0;
}

int MpiWrapper::Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[])
{
#ifdef GEOSX_USE_MPI
  return MPI_Waitall( count,  array_of_requests, array_of_statuses );
#endif
  return 0;
}

double MpiWrapper::Wtime(void)
{
#ifdef GEOSX_USE_MPI
  return MPI_Wtime( );
#else
  return 0;
#endif

}




} /* namespace geosx */

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif


