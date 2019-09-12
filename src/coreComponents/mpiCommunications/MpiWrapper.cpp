/*
 * MpiWrapper.cpp
 *
 *  Created on: Sep 11, 2019
 *      Author: settgast
 */

#include "MpiWrapper.hpp"

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



int MpiWrapper::Init( int *argc, char ***argv )
{
#ifdef GEOSX_USE_MPI
  return MPI_Init( argc, argv );
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

int MpiWrapper::Cart_Rank( MPI_Comm comm, const int coords[] )
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




} /* namespace geosx */
