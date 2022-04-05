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

void MpiWrapper::barrier( MPI_Comm const & MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_MPI
  MPI_Barrier( comm );
#endif
}

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

int MpiWrapper::commRank( MPI_Comm const & MPI_PARAM( comm ) )
{
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( comm, &rank );
#endif
  return rank;
}

int MpiWrapper::commSize( MPI_Comm const & MPI_PARAM( comm ) )
{
  int size = 1;
#ifdef GEOSX_USE_MPI
  MPI_Comm_size( comm, &size );
#endif
  return size;
}

bool MpiWrapper::commCompare( MPI_Comm const & comm1, MPI_Comm const & comm2 )
{
#ifdef GEOSX_USE_MPI
  int result;
  MPI_Comm_compare( comm1, comm2, &result );
  return result == MPI_IDENT || result == MPI_CONGRUENT;
#else
  return comm1 == comm2;
#endif
}

bool MpiWrapper::initialized()
{
#ifdef GEOSX_USE_MPI
  int ret = false;
  MPI_CHECK_ERROR( MPI_Initialized( &ret ) );
  return ret;
#else
  return false;
#endif
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

int MpiWrapper::testAny( int count, MPI_Request array_of_requests[], int * idx, int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Testany( count, array_of_requests, idx, flag, array_of_statuses );
#endif
  *flag = 0;
  return 0;
}

int MpiWrapper::testSome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Testsome( count, array_of_requests, outcount, array_of_indices, array_of_statuses );
#endif
  *outcount = 0;
  return 0;
}

int MpiWrapper::testAll( int count, MPI_Request array_of_requests[], int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Testall( count, array_of_requests, flag, array_of_statuses );
#endif
  *flag = 0;
  return 0;
}

int MpiWrapper::check( MPI_Request * request, int * flag, MPI_Status * status )
{
#ifdef GEOSX_USE_MPI
  return MPI_Request_get_status( *request, flag, status );
#endif
  *flag = 0;
  return 0;
}

int MpiWrapper::checkAny( int count, MPI_Request array_of_requests[], int * idx, int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  bool found = false;
  int flagCache = -1;
  int rval = MPI_SUCCESS;
  std::vector< int > rvals( count );
  for( int jdx = 0; jdx < count; ++jdx )
  {
    *flag = 0;
    rvals[ jdx ] = MPI_Request_get_status( array_of_requests[ jdx ], flag, &array_of_statuses[ jdx ] );
    if( *flag && !found )
    {
      *idx = jdx;
      flagCache = *flag;
    }
    if( rvals[ jdx ] != MPI_SUCCESS )
    {
      rval = rvals[ jdx ];
    }
  }
  if( found )
  {
    *flag = flagCache;
  }
  return rval;
#endif
  *flag = 0;
  return 0;
}

int MpiWrapper::checkAll( int count, MPI_Request array_of_requests[], int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  // assume all passing, any that don't pass set the flag to false
  *flag = 1;
  int rval = MPI_SUCCESS;
  std::vector< int > rvals( count );
  int iFlag = 0;
  for( int idx = 0; idx < count; ++idx )
  {
    rvals[ idx ] = MPI_Request_get_status( array_of_requests[ idx ], &iFlag, &array_of_statuses[ idx ] );
    if( !iFlag )
    {
      *flag = iFlag;
    }
    if( rvals[ idx ] != MPI_SUCCESS )
    {
      rval = rvals[ idx ];
    }
  }
  return rval;
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

int MpiWrapper::waitAny( int count, MPI_Request array_of_requests[], int * indx, MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Waitany( count, array_of_requests, indx, array_of_statuses );
#endif
  return 0;
}

int MpiWrapper::waitSome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] )
{
#ifdef GEOSX_USE_MPI
  return MPI_Waitsome( count, array_of_requests, outcount, array_of_indices, array_of_statuses );
#endif
  // *outcount = 0;
  return 0;
}

int MpiWrapper::waitAll( int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[] )
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
    int err = waitAny( count, array_of_requests, &idx, &stat );
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
    int err = waitSome( count, array_of_requests, &rcvd, &indices[0], &stats[0] );
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
