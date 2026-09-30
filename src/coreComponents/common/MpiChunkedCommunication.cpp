/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MpiChunkedCommunication.cpp
 */

#include "MpiChunkedCommunication.hpp"

#include <algorithm>
#include <climits>
#include <stdexcept>

namespace geos::mpi
{
namespace
{
void validateChunkSize( std::uint64_t chunkBytes )
{
  if( chunkBytes == 0 || chunkBytes > INT_MAX )
    throw std::invalid_argument( "MPI byte chunk size must be in [1, INT_MAX]" );
}
} // namespace

void sendBytes( void const * buffer, std::uint64_t bytes, int destination, int tag, MPI_Comm comm, std::uint64_t chunkBytes )
{
  validateChunkSize( chunkBytes );
#ifdef GEOS_USE_MPI
  auto ptr = static_cast< unsigned char const * >( buffer );
  while( bytes > 0 )
  {
    int const chunk = static_cast< int >( std::min( bytes, chunkBytes ) );
    MPI_Send( ptr, chunk, MPI_BYTE, destination, tag, comm );
    ptr += chunk;
    bytes -= chunk;
  }
#else
  GEOS_UNUSED_VAR( buffer, destination, tag, comm );
  if( bytes )
    throw std::logic_error( "MPI byte send requires an MPI build" );
#endif
}

void receiveBytes( void * buffer, std::uint64_t bytes, int source, int tag, MPI_Comm comm, std::uint64_t chunkBytes )
{
  validateChunkSize( chunkBytes );
#ifdef GEOS_USE_MPI
  auto ptr = static_cast< unsigned char * >( buffer );
  while( bytes > 0 )
  {
    int const chunk = static_cast< int >( std::min( bytes, chunkBytes ) );
    MPI_Recv( ptr, chunk, MPI_BYTE, source, tag, comm, MPI_STATUS_IGNORE );
    ptr += chunk;
    bytes -= chunk;
  }
#else
  GEOS_UNUSED_VAR( buffer, source, tag, comm );
  if( bytes )
    throw std::logic_error( "MPI byte receive requires an MPI build" );
#endif
}

void exchangeBytes( void const * sendBuffer, std::uint64_t sendLength, void * receiveBuffer, std::uint64_t receiveLength, int peer, int tag,
                    MPI_Comm comm, std::uint64_t chunkBytes )
{
  validateChunkSize( chunkBytes );
#ifdef GEOS_USE_MPI
  auto send = static_cast< unsigned char const * >( sendBuffer );
  auto receive = static_cast< unsigned char * >( receiveBuffer );
  while( sendLength || receiveLength )
  {
    int const sending = static_cast< int >( std::min( sendLength, chunkBytes ) );
    int const receiving = static_cast< int >( std::min( receiveLength, chunkBytes ) );
    MPI_Sendrecv( send, sending, MPI_BYTE, peer, tag, receive, receiving, MPI_BYTE, peer, tag, comm, MPI_STATUS_IGNORE );
    if( sending )
      send += sending;
    if( receiving )
      receive += receiving;
    sendLength -= sending;
    receiveLength -= receiving;
  }
#else
  GEOS_UNUSED_VAR( sendBuffer, receiveBuffer, peer, tag, comm );
  if( sendLength || receiveLength )
    throw std::logic_error( "MPI byte exchange requires an MPI build" );
#endif
}

void exchangeManyBytes( std::vector< ByteExchange > const & exchanges, MPI_Request * requests, int tag, MPI_Comm comm,
                        std::uint64_t chunkBytes )
{
  validateChunkSize( chunkBytes );
#ifdef GEOS_USE_MPI
  std::uint64_t offset = 0;
  bool pending = true;
  while( pending )
  {
    pending = false;
    for( std::size_t i = 0; i < exchanges.size(); ++i )
    {
      auto const & exchange = exchanges[i];
      requests[2 * i] = MPI_REQUEST_NULL;
      requests[2 * i + 1] = MPI_REQUEST_NULL;
      if( offset < exchange.receiveBytes )
      {
        int const count = static_cast< int >( std::min( exchange.receiveBytes - offset, chunkBytes ) );
        auto * buffer = static_cast< unsigned char * >( exchange.receiveBuffer ) + offset;
        MPI_Irecv( buffer, count, MPI_BYTE, exchange.peer, tag, comm, &requests[2 * i] );
        pending = pending || exchange.receiveBytes - offset > chunkBytes;
      }
    }
    for( std::size_t i = 0; i < exchanges.size(); ++i )
    {
      auto const & exchange = exchanges[i];
      if( offset < exchange.sendBytes )
      {
        int const count = static_cast< int >( std::min( exchange.sendBytes - offset, chunkBytes ) );
        auto const * buffer = static_cast< unsigned char const * >( exchange.sendBuffer ) + offset;
        MPI_Isend( buffer, count, MPI_BYTE, exchange.peer, tag, comm, &requests[2 * i + 1] );
        pending = pending || exchange.sendBytes - offset > chunkBytes;
      }
    }
    if( !exchanges.empty() )
      MPI_Waitall( static_cast< int >( 2 * exchanges.size() ), requests, MPI_STATUSES_IGNORE );
    // Increment only when another chunk exists, avoiding overflow at UINT64_MAX.
    if( pending )
      offset += chunkBytes;
  }
#else
  GEOS_UNUSED_VAR( requests, tag, comm );
  for( auto const & exchange : exchanges )
    if( exchange.sendBytes || exchange.receiveBytes )
      throw std::logic_error( "MPI byte exchange requires an MPI build" );
#endif
}
} // namespace geos::mpi
