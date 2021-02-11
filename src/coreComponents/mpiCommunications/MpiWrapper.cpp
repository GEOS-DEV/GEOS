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
 * @file MpiWrapper.cpp
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

//int MpiWrapper::Bcast( void * buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
//{
//#ifdef GEOSX_USE_MPI
//  return MPI_Bcast( buffer, count, datatype, root, comm );
//#else
//  return 0;
//#endif
//
//}



int MpiWrapper::cartCoords( MPI_Comm comm, int rank, int maxdims, int coords[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Cart_coords( comm, rank, maxdims, coords );
#else
  return 0;
#endif
}

int MpiWrapper::cartCreate( MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                            int reorder, MPI_Comm * comm_cart )
{
#ifdef GEOSX_USE_MPI
  return MPI_Cart_create( comm_old, ndims, dims, periods, reorder, comm_cart );
#else
  return 0;
#endif
}

int MpiWrapper::cartRank( MPI_Comm comm, const int coords[] )
{
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Cart_rank( comm, coords, &rank );
#endif
  return rank;
}

void MpiWrapper::commFree( MPI_Comm & comm )
{
#ifdef GEOSX_USE_MPI
  MPI_CHECK_ERROR( MPI_Comm_free( &comm ) );
#else
//  comm = MPI_COMM_NULL;
#endif
}

std::size_t MpiWrapper::getSizeofMpiType( MPI_Datatype const type )
{
  if( type == MPI_CHAR )
  {
    return sizeof(char);
  }
  else if( type == MPI_FLOAT )
  {
    return sizeof(float);
  }
  else if( type == MPI_DOUBLE )
  {
    return sizeof(double);
  }
  else if( type == MPI_INT )
  {
    return sizeof(int);
  }
  else if( type == MPI_LONG )
  {
    return sizeof(long int);
  }
  else if( type == MPI_LONG_LONG )
  {
    return sizeof(long long int);
  }
  else
  {
    GEOSX_ERROR( "No conversion implemented for MPI_Datatype "<<type );
  }
  return 0;
}


int MpiWrapper::init( int * argc, char * * * argv )
{
#ifdef GEOSX_USE_MPI
  return MPI_Init( argc, argv );
#else
  return 0;
#endif
}

void MpiWrapper::finalize()
{
#ifdef GEOSX_USE_MPI
  MPI_CHECK_ERROR( MPI_Finalize() );
#endif
}


MPI_Comm MpiWrapper::commDup( MPI_Comm const comm )
{
#ifdef GEOSX_USE_MPI
  MPI_Comm duplicate;
  MPI_CHECK_ERROR( MPI_Comm_dup( comm, &duplicate ) );
  return duplicate;
#else
  return comm;
#endif
}

MPI_Comm MpiWrapper::commSplit( MPI_Comm const comm, int color, int key )
{
#ifdef GEOSX_USE_MPI
  MPI_Comm scomm;
  MPI_CHECK_ERROR( MPI_Comm_split( comm, color, key, &scomm ) );
  return scomm;
#else
  return comm;
#endif
}

int MpiWrapper::test( MPI_Request * request, int * flag, MPI_Status * status )
{
#ifdef GEOSX_USE_MPI
  return MPI_Test( request, flag, status );
#endif
  *flag = 0;
  return 0;
}

int MpiWrapper::testany( int count,  MPI_Request array_of_requests[], int * idx, int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Testany( count, array_of_requests, idx, flag, array_of_statuses );
#endif
  *flag = 0;
  return 0;
}

int MpiWrapper::testall( int count,  MPI_Request array_of_requests[], int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Testall( count, array_of_requests, flag, array_of_statuses );
#endif
  *flag = 0;
  return 0;
}

int MpiWrapper::wait( MPI_Request * request, MPI_Status * status )
{
#ifdef GEOSX_USE_MPI
  return MPI_Wait( request, status );
#endif
  return 0;
}

int MpiWrapper::waitany( int count, MPI_Request array_of_requests[], int * indx, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Waitany( count, array_of_requests, indx, array_of_statuses );
#endif
  return 0;
}

int MpiWrapper::waitsome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Waitsome( count, array_of_requests, outcount, array_of_indices, array_of_statuses );
#endif
  // *outcount = 0;
  return 0;
}

int MpiWrapper::waitall( int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Waitall( count, array_of_requests, array_of_statuses );
#endif
  return 0;
}

double MpiWrapper::wtime( void )
{
#ifdef GEOSX_USE_MPI
  return MPI_Wtime( );
#else
  return 0;
#endif

}

int MpiWrapper::activeWaitAny( const int count, MPI_Request array_of_requests[], std::function< void ( int ) > func )
{
  int cmp = 0;
  while( cmp < count )
  {
    int idx = 0;
    MPI_Status stat;
    int err = waitany( count, array_of_requests, &idx, &stat );
    if( err != MPI_SUCCESS )
      return err;
    if( idx != MPI_UNDEFINED )   // only if all(requests == MPI_REQUEST_NULL)
    {
      func( idx );
    }
    cmp++;
  }
  return MPI_SUCCESS;
}

int MpiWrapper::activeWaitSome( const int count, MPI_Request array_of_requests[], std::function< void ( int ) > func )
{
  int cmp = 0;
  while( cmp < count )
  {
    int rcvd = 0;
    std::vector< int > indices( count, -1 );
    std::vector< MPI_Status > stats( count );
    int err = waitsome( count, array_of_requests, &rcvd, &indices[0], &stats[0] );
    if( err != MPI_SUCCESS )
      return err;
    if( rcvd > 0 )
    {
      for( int ii = 0; ii < rcvd; ++ii )
      {
        if( indices[ii] != MPI_UNDEFINED )
        {
          func( indices[ii] );
        }
      }
    }
    cmp += rcvd;
  }
  return MPI_SUCCESS;
}

int MpiWrapper::activeWaitSomePartialPhase( const int participants,
                                            std::vector< std::function< MPI_Request ( int ) > > const & phases )
{
  const int num_phases = sizeof(phases.size());
  std::vector< MPI_Request > phase_requests( participants * num_phases, MPI_REQUEST_NULL );
  for( int idx = 0; idx < participants; ++idx )
  {
    phase_requests[idx] = phases[0]( idx );
  }
  auto phase_invocation = [&] ( int idx )
  {
    int phase = (idx / participants) + 1;
    int phase_idx = idx % participants;
    phase_requests[idx + participants] = phases[phase]( phase_idx );
  };
  return activeWaitSome( participants * num_phases, &phase_requests[0], phase_invocation );
}

int MpiWrapper::activeWaitSomeCompletePhase( const int participants,
                                             std::vector< std::function< MPI_Request ( int ) > > const & phases )
{
  const int num_phases = phases.size();
  std::vector< MPI_Request > phase_requests( num_phases * participants, MPI_REQUEST_NULL );
  for( int idx = 0; idx < participants; ++idx )
  {
    phase_requests[idx] = phases[0]( idx );
  }
  int err = 0;
  for( int phase = 1; phase < num_phases; ++phase )
  {
    int prev_phase = phase - 1;
    auto phase_wrapper = [&] ( int idx ) { phase_requests[ ( phase * participants ) + idx ] = phases[phase]( idx ); };
    err = activeWaitSome( participants, &phase_requests[prev_phase * participants], phase_wrapper );
    if( err != MPI_SUCCESS )
      break;
  }
  return err;
}

int MpiWrapper::activeWaitOrderedCompletePhase( const int participants,
                                                std::vector< std::function< MPI_Request ( int ) > > const & phases )
{
  const int num_phases = phases.size();
  std::vector< MPI_Request > phase_requests( participants );
  for( int idx = 0; idx < participants; ++idx )
  {
    phase_requests[idx] = phases[0]( idx );
  }
  for( int phase = 1; phase < num_phases; ++phase )
  {
    for( int idx = 0; idx < participants; ++idx )
    {
      MPI_Status stat;
      wait( &phase_requests[idx], &stat );
      phase_requests[idx] = phases[phase]( idx );
    }
  }
  return MPI_SUCCESS;
}

} /* namespace geosx */

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
