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

/**
 * @brief Blocking send of a byte buffer in chunks.
 * @param buffer The bytes to send.
 * @param bytes Number of bytes to send.
 * @param destination The receiving rank.
 * @param tag The MPI tag.
 * @param comm The MPI communicator.
 * @param chunkBytes Largest chunk size; must match the receiver's.
 */
void sendBytes( void const * buffer, std::uint64_t bytes, int destination, int tag, MPI_Comm comm,
                std::uint64_t chunkBytes = defaultChunkBytes );
/**
 * @brief Blocking receive of a byte buffer sent with sendBytes.
 * @param buffer Destination of the received bytes.
 * @param bytes Number of bytes to receive.
 * @param source The sending rank.
 * @param tag The MPI tag.
 * @param comm The MPI communicator.
 * @param chunkBytes Largest chunk size; must match the sender's.
 */
void receiveBytes( void * buffer, std::uint64_t bytes, int source, int tag, MPI_Comm comm, std::uint64_t chunkBytes = defaultChunkBytes );

/** Exchange with one peer; both sides know the two lengths before calling.
 * Zero-length chunks participate when the lengths differ. No collective occurs.
 * @param sendBuffer The bytes to send.
 * @param sendBytes Number of bytes to send.
 * @param receiveBuffer Destination of the received bytes.
 * @param receiveBytes Number of bytes to receive.
 * @param peer The other rank.
 * @param tag The MPI tag.
 * @param comm The MPI communicator.
 * @param chunkBytes Largest chunk size; must match the peer's.
 */
void exchangeBytes( void const * sendBuffer, std::uint64_t sendBytes, void * receiveBuffer, std::uint64_t receiveBytes, int peer, int tag,
                    MPI_Comm comm, std::uint64_t chunkBytes = defaultChunkBytes );

/// One peer of an exchangeManyBytes call.
struct ByteExchange
{
  void const * sendBuffer;    ///< The bytes to send.
  std::uint64_t sendBytes;    ///< Number of bytes to send.
  void * receiveBuffer;       ///< Destination of the received bytes.
  std::uint64_t receiveBytes; ///< Number of bytes to receive.
  int peer;                   ///< The other rank.
};

/** Exchange chunks with all peers concurrently, using at most two requests per peer.
 * The caller must preflight unique peers, matching lengths, valid buffers/chunk size,
 * and exchanges.size() <= INT_MAX / 2, and allocate 2 * exchanges.size() requests.
 * This routine allocates no memory after communication begins. All peers use the
 * same tag/chunk size. Buffers remain valid until this blocking routine returns.
 * @param exchanges One entry per peer.
 * @param requests Storage for 2 * exchanges.size() MPI requests.
 * @param tag The MPI tag.
 * @param comm The MPI communicator.
 * @param chunkBytes Largest chunk size; must match the peers'.
 */
void exchangeManyBytes( std::vector< ByteExchange > const & exchanges, MPI_Request * requests, int tag, MPI_Comm comm,
                        std::uint64_t chunkBytes = defaultChunkBytes );

/**
 * @brief Sparse personalized exchange of byte buffers.
 * @details Each rank sends one buffer to each peer listed in @p outgoing and
 * receives the buffers that other ranks address to it. Receivers learn their
 * senders with synchronous sends and a nonblocking barrier (the NBX algorithm),
 * so the cost depends on the number of peers of a rank, not on the size of the
 * communicator. A buffer addressed to the calling rank is returned unchanged.
 * Collective over @p comm. Each call communicates on its own duplicate of
 * @p comm, so consecutive calls cannot exchange each other's messages.
 * @param outgoing The buffer for each destination rank.
 * @param comm The MPI communicator.
 * @param tag The MPI tag of the handshake; payloads use @p tag + 1.
 * @param chunkBytes Largest chunk size of a payload message.
 * @return The buffer received from each source rank.
 */
stdMap< int, stdVector< char > > sparseExchange( stdMap< int, stdVector< char > > const & outgoing,
                                                 MPI_Comm comm, int tag = 7401,
                                                 std::uint64_t chunkBytes = defaultChunkBytes );

} // namespace geos::mpi
#endif
