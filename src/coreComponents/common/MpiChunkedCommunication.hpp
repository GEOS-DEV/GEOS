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
 * @file MpiChunkedCommunication.hpp
 */

#ifndef GEOS_COMMON_MPICHUNKEDCOMMUNICATION_HPP
#define GEOS_COMMON_MPICHUNKEDCOMMUNICATION_HPP

#include "MpiWrapper.hpp"

#include <cstdint>
#include <vector>

namespace geos::mpi
{

/** Blocking byte transport without narrowing a whole payload to MPI's int count.
 * Callers validate buffer lengths and allocations before entering communication.
 * Every peer must use the same chunk size. A smaller value is useful in tests.
 */
constexpr std::uint64_t defaultChunkBytes = UINT64_C( 1 ) << 30;

void sendBytes( void const * buffer, std::uint64_t bytes, int destination, int tag, MPI_Comm comm,
                std::uint64_t chunkBytes = defaultChunkBytes );
void receiveBytes( void * buffer, std::uint64_t bytes, int source, int tag, MPI_Comm comm, std::uint64_t chunkBytes = defaultChunkBytes );

/** Exchange with one peer; both sides know the two lengths before calling.
 * Zero-length chunks participate when the lengths differ. No collective occurs.
 */
void exchangeBytes( void const * sendBuffer, std::uint64_t sendBytes, void * receiveBuffer, std::uint64_t receiveBytes, int peer, int tag,
                    MPI_Comm comm, std::uint64_t chunkBytes = defaultChunkBytes );

struct ByteExchange
{
  void const * sendBuffer;
  std::uint64_t sendBytes;
  void * receiveBuffer;
  std::uint64_t receiveBytes;
  int peer;
};

/** Exchange chunks with all peers concurrently, using at most two requests per peer.
 * The caller must preflight unique peers, matching lengths, valid buffers/chunk size,
 * and exchanges.size() <= INT_MAX / 2, and allocate 2 * exchanges.size() requests.
 * This routine allocates no memory after communication begins. All peers use the
 * same tag/chunk size. Buffers remain valid until this blocking routine returns.
 */
void exchangeManyBytes( std::vector< ByteExchange > const & exchanges, MPI_Request * requests, int tag, MPI_Comm comm,
                        std::uint64_t chunkBytes = defaultChunkBytes );

} // namespace geos::mpi
#endif
