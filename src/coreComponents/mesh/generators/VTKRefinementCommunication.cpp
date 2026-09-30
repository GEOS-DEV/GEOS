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
 * @file VTKRefinementCommunication.cpp
 */

#include "VTKRefinementCommunication.hpp"

#include "common/MpiChunkedCommunication.hpp"
#include "common/TimingMacros.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>

namespace geos::vtk::refinement
{
namespace
{
void putInteger( Bytes & bytes, std::uint64_t value )
{
  for( int i = 0; i < 8; ++i )
  {
    bytes.push_back( static_cast< unsigned char >( value >> ( 8 * i ) ) );
  }
}
void putDouble( Bytes & bytes, double value )
{
  static_assert( sizeof( double ) == sizeof( std::uint64_t ) );
  std::uint64_t bits;
  std::memcpy( &bits, &value, sizeof( bits ) );
  putInteger( bytes, bits );
}
void putKey( Bytes & bytes, EntityKey const & key )
{
  putInteger( bytes, key.meshNamespace );
  putInteger( bytes, static_cast< std::uint64_t >( key.kind ) );
  putInteger( bytes, key.corners.size() );
  for( vtkIdType id : key.corners )
  {
    putInteger( bytes, static_cast< std::uint64_t >( id ) );
  }
}
void validateKey( EntityKey const & key )
{
  std::size_t const n = key.corners.size();
  switch( key.kind )
  {
  case EntityKind::vertex:
  case EntityKind::cell:
    if( n != 1 )
    {
      throw std::invalid_argument( "Vertex/cell refinement key must contain one ID" );
    }
    break;
  case EntityKind::edge:
    if( n != 2 || key.corners[0] >= key.corners[1] )
    {
      throw std::invalid_argument( "Noncanonical refinement edge key" );
    }
    break;
  case EntityKind::face:
    if( n < 3 || n > 11 || canonicalCycle( key.corners ) != key.corners )
    {
      throw std::invalid_argument( "Noncanonical refinement face key" );
    }
    break;
  default:
    throw std::invalid_argument( "Invalid refinement entity kind" );
  }
  for( vtkIdType id : key.corners )
  {
    if( id < 0 )
    {
      throw std::invalid_argument( "Negative refinement support ID" );
    }
  }
}
class Reader
{
public:
  explicit Reader( Bytes const & bytes ) : m_bytes( bytes ) {}
  bool done() const { return m_cursor == m_bytes.size(); }
  std::uint64_t integer()
  {
    if( m_bytes.size() - m_cursor < 8 )
    {
      throw std::invalid_argument( "Truncated refinement integral record" );
    }
    std::uint64_t value = 0;
    for( int i = 0; i < 8; ++i )
    {
      value |= static_cast< std::uint64_t >( m_bytes[m_cursor++] ) << ( 8 * i );
    }
    return value;
  }
  vtkIdType id()
  {
    std::uint64_t const value = integer();
    if( value > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) )
    {
      throw std::overflow_error( "Refinement received ID exceeds vtkIdType" );
    }
    return static_cast< vtkIdType >( value );
  }
  double real()
  {
    auto const bits = integer();
    double value;
    std::memcpy( &value, &bits, sizeof( value ) );
    if( !std::isfinite( value ) )
    {
      throw std::invalid_argument( "Nonfinite refinement point record" );
    }
    return value;
  }
  EntityKey key()
  {
    EntityKey result;
    result.meshNamespace = integer();
    auto const kind = integer();
    if( kind > static_cast< std::uint64_t >( EntityKind::cell ) )
    {
      throw std::invalid_argument( "Invalid wire entity kind" );
    }
    result.kind = static_cast< EntityKind >( kind );
    auto const n = integer();
    if( n > 11 )
    {
      throw std::invalid_argument( "Refinement wire key exceeds supported arity" );
    }
    for( std::uint64_t i = 0; i < n; ++i )
    {
      result.corners.push_back( id() );
    }
    validateKey( result );
    return result;
  }
  Participants participants( int size )
  {
    auto const n = integer();
    if( n == 0 || n > static_cast< std::uint64_t >( size ) )
    {
      throw std::invalid_argument( "Invalid refinement participant count" );
    }
    Participants result;
    for( std::uint64_t i = 0; i < n; ++i )
    {
      auto const rank = integer();
      if( rank >= static_cast< std::uint64_t >( size ) )
      {
        throw std::invalid_argument( "Invalid refinement participant rank" );
      }
      result.push_back( static_cast< int >( rank ) );
    }
    if( !std::is_sorted( result.begin(), result.end() ) || std::adjacent_find( result.begin(), result.end() ) != result.end() )
    {
      throw std::invalid_argument( "Noncanonical refinement participant list" );
    }
    return result;
  }
  Bytes payload()
  {
    auto const n = integer();
    if( n > m_bytes.size() - m_cursor )
    {
      throw std::invalid_argument( "Truncated refinement field record" );
    }
    Bytes result( m_bytes.begin() + m_cursor, m_bytes.begin() + m_cursor + static_cast< std::size_t >( n ) );
    m_cursor += n;
    return result;
  }

private:
  Bytes const & m_bytes;
  std::size_t m_cursor{};
};
std::uint64_t checkedSum( std::uint64_t a, std::uint64_t b )
{
  if( b > std::numeric_limits< std::uint64_t >::max() - a )
  {
    throw std::overflow_error( "Refinement global count overflow" );
  }
  return a + b;
}
EntityKey faceCornerSet( EntityKey key )
{
  if( key.kind == EntityKind::face )
  {
    std::sort( key.corners.begin(), key.corners.end() );
  }
  return key;
}
} // namespace

bool ChildCellKey::operator<( ChildCellKey const & other ) const
{
  return std::tie( meshNamespace, parentId, templateCode, ordinal ) <
         std::tie( other.meshNamespace, other.parentId, other.templateCode, other.ordinal );
}

Communication::Communication( MPI_Comm comm, std::uint64_t chunkBytes )
    : m_comm( MpiWrapper::commDup( comm ) ), m_rank( MpiWrapper::commRank( m_comm ) ), m_size( MpiWrapper::commSize( m_comm ) ),
      m_chunkBytes( chunkBytes )
{
  try
  {
    auto const low = MpiWrapper::allReduce( chunkBytes, MpiWrapper::Reduction::Min, m_comm );
    auto const high = MpiWrapper::allReduce( chunkBytes, MpiWrapper::Reduction::Max, m_comm );
    checked( "communication configuration",
             [&]
             {
               if( low != high || chunkBytes == 0 || chunkBytes > INT_MAX )
               {
                 throw std::invalid_argument( "Inconsistent or invalid MPI refinement chunk size" );
               }
             } );
  }
  catch( ... )
  {
    MpiWrapper::commFree( m_comm );
    throw;
  }
}
Communication::~Communication()
{
  if( m_neighborComm != MPI_COMM_NULL )
  {
    MpiWrapper::commFree( m_neighborComm );
  }
  MpiWrapper::commFree( m_comm );
}

void Communication::checked( std::string const & phase, std::function< void() > const & work ) const
{
  char diagnostic[1024]{};
  int failedRank = m_size;
  try
  {
    work();
  }
  catch( std::exception const & error )
  {
    failedRank = m_rank;
    std::strncpy( diagnostic, error.what(), sizeof( diagnostic ) - 1 );
  }
  catch( ... )
  {
    failedRank = m_rank;
    std::strncpy( diagnostic, "Unknown local refinement validation error", sizeof( diagnostic ) - 1 );
  }
  failedRank = MpiWrapper::allReduce( failedRank, MpiWrapper::Reduction::Min, m_comm );
  if( failedRank < m_size )
  {
    MpiWrapper::bcast( diagnostic, static_cast< int >( sizeof( diagnostic ) ), failedRank, m_comm );
    throw std::runtime_error( "Uniform refinement " + phase + ", rank " + std::to_string( failedRank ) + ": " + diagnostic );
  }
}

Communication::Mail Communication::exchangePayloads( Mail const & outgoing, std::map< int, std::uint64_t > const & incomingLengths )
{
  Mail incoming;
  std::set< int > peers;
  std::vector< mpi::ByteExchange > exchanges;
  std::vector< MPI_Request > requests;
  std::uint64_t bytesSent = m_statistics.payloadBytesSent, chunksSent = m_statistics.payloadChunksSent;
  checked( "receive allocation",
           [&]
           {
             std::uint64_t total = 0;
             for( auto const & [peer, n] : incomingLengths )
             {
               total = checkedSum( total, n );
               if( n > Bytes{}.max_size() || total > std::numeric_limits< std::size_t >::max() )
               {
                 throw std::overflow_error( "Refinement receive buffers exceed addressable storage" );
               }
             }
             for( auto const & [peer, n] : incomingLengths )
             {
               if( n )
               {
                 incoming[peer].resize( static_cast< std::size_t >( n ) );
               }
             }
             auto const self = outgoing.find( m_rank );
             if( self != outgoing.end() )
             {
               incoming[m_rank] = self->second;
             }
             for( auto const & mail : outgoing )
             {
               if( mail.first != m_rank && !mail.second.empty() )
               {
                 peers.insert( mail.first );
                 bytesSent = checkedSum( bytesSent, mail.second.size() );
                 chunksSent = checkedSum( chunksSent, mail.second.size() / m_chunkBytes + ( mail.second.size() % m_chunkBytes != 0 ) );
               }
             }
             for( auto const & mail : incoming )
             {
               if( mail.first != m_rank && !mail.second.empty() )
               {
                 peers.insert( mail.first );
               }
             }
             if( peers.size() > static_cast< std::size_t >( INT_MAX / 2 ) )
             {
               throw std::overflow_error( "Too many refinement MPI requests" );
             }
             exchanges.reserve( peers.size() );
             requests.resize( 2 * peers.size() );
             for( int peer : peers )
             {
               auto const send = outgoing.find( peer );
               auto const receive = incoming.find( peer );
               exchanges.push_back( { send == outgoing.end() ? nullptr : send->second.data(),
                                      send == outgoing.end() ? 0 : send->second.size(),
                                      receive == incoming.end() ? nullptr : receive->second.data(),
                                      receive == incoming.end() ? 0 : receive->second.size(), peer } );
             }
           } );
  mpi::exchangeManyBytes( exchanges, requests.data(), 1, m_comm, m_chunkBytes );
  m_statistics.payloadBytesSent = bytesSent;
  m_statistics.payloadChunksSent = chunksSent;
  return incoming;
}

Communication::Mail Communication::exchangeDirectory( Mail const & outgoing )
{
  std::vector< std::uint64_t > sending, receiving;
  std::uint64_t countBytes = 0, directoryExchanges = 0;
  checked( "directory packing lengths",
           [&]
           {
             sending.resize( m_size );
             receiving.resize( m_size );
             countBytes = checkedSum( m_statistics.countBytesSent, 8 * static_cast< std::uint64_t >( m_size ) );
             directoryExchanges = checkedSum( m_statistics.directoryExchanges, 1 );
             for( auto const & [peer, bytes] : outgoing )
             {
               if( peer < 0 || peer >= m_size )
               {
                 throw std::invalid_argument( "Invalid refinement directory destination" );
               }
               sending[peer] = bytes.size();
             }
           } );
#ifdef GEOS_USE_MPI
  MPI_Alltoall( sending.data(), 1, MPI_UINT64_T, receiving.data(), 1, MPI_UINT64_T, m_comm );
#else
  receiving = sending;
#endif
  std::map< int, std::uint64_t > lengths;
  checked( "directory receive lengths",
           [&]
           {
             for( int peer = 0; peer < m_size; ++peer )
             {
               if( receiving[peer] )
               {
                 lengths[peer] = receiving[peer];
               }
             }
           } );
  m_statistics.directoryExchanges = directoryExchanges;
  m_statistics.countBytesSent = countBytes;
  return exchangePayloads( outgoing, lengths );
}

void Communication::initializeNeighbors( Participants neighbors )
{
  std::sort( neighbors.begin(), neighbors.end() );
  neighbors.erase( std::unique( neighbors.begin(), neighbors.end() ), neighbors.end() );
  neighbors.erase( std::remove( neighbors.begin(), neighbors.end(), m_rank ), neighbors.end() );
  checked( "neighbor graph",
           [&]
           {
             for( int r : neighbors )
             {
               if( r < 0 || r >= m_size )
               {
                 throw std::invalid_argument( "Invalid refinement neighbor" );
               }
             }
           } );
  if( m_neighborComm != MPI_COMM_NULL )
  {
    MpiWrapper::commFree( m_neighborComm );
  }
#ifdef GEOS_USE_MPI
  int const degree = static_cast< int >( neighbors.size() );
  // Directory replies supply the same participant list to every participant,
  // so both local adjacency lists are already known. Avoid rediscovering the
  // incoming graph through MPI_Dist_graph_create's distributed edge routing.
  int const emptyNeighbor = 0;
  int const * adjacency = neighbors.empty() ? &emptyNeighbor : neighbors.data();
  MPI_Dist_graph_create_adjacent( m_comm, degree, adjacency, MPI_UNWEIGHTED, degree, adjacency, MPI_UNWEIGHTED, MPI_INFO_NULL, 0,
                                  &m_neighborComm );
  int incoming = 0, outgoing = 0, weighted = 0;
  MPI_Dist_graph_neighbors_count( m_neighborComm, &incoming, &outgoing, &weighted );
  checked( "neighbor graph allocation",
           [&]
           {
             m_sources.resize( incoming );
             m_destinations.resize( outgoing );
           } );
  MPI_Dist_graph_neighbors( m_neighborComm, incoming, m_sources.data(), MPI_UNWEIGHTED, outgoing, m_destinations.data(), MPI_UNWEIGHTED );
  checked( "neighbor symmetry",
           [&]
           {
             auto sources = m_sources, destinations = m_destinations;
             std::sort( sources.begin(), sources.end() );
             std::sort( destinations.begin(), destinations.end() );
             if( sources != destinations || destinations != neighbors )
             {
               throw std::invalid_argument( "Asymmetric refinement participant graph" );
             }
           } );
#else
  m_sources = neighbors;
  m_destinations = neighbors;
#endif
  m_neighbors = std::move( neighbors );
  m_discovered = true;
}

Communication::Mail Communication::exchangeNeighbors( Mail const & outgoing )
{
  std::vector< std::uint64_t > sending, receiving;
  std::uint64_t countBytes = 0, neighborExchanges = 0;
  checked( "neighbor packing lengths",
           [&]
           {
             sending.resize( m_destinations.size() );
             receiving.resize( m_sources.size() );
             countBytes = checkedSum( m_statistics.countBytesSent, 8 * m_destinations.size() );
             neighborExchanges = checkedSum( m_statistics.neighborExchanges, 1 );
             for( auto const & [peer, bytes] : outgoing )
             {
               GEOS_UNUSED_VAR( bytes );
               if( peer != m_rank && !std::binary_search( m_neighbors.begin(), m_neighbors.end(), peer ) )
               {
                 throw std::invalid_argument( "Fine refinement point uses an undiscovered participant" );
               }
             }
             for( std::size_t i = 0; i < m_destinations.size(); ++i )
             {
               auto found = outgoing.find( m_destinations[i] );
               if( found != outgoing.end() )
               {
                 sending[i] = found->second.size();
               }
             }
           } );
#ifdef GEOS_USE_MPI
  // MPICH validates the pointers when count is nonzero, including ranks with
  // zero graph degree. Such ranks still participate in this collective.
  std::uint64_t emptySend = 0, emptyReceive = 0;
  MPI_Neighbor_alltoall( sending.empty() ? &emptySend : sending.data(), 1, MPI_UINT64_T,
                         receiving.empty() ? &emptyReceive : receiving.data(), 1, MPI_UINT64_T, m_neighborComm );
#endif
  std::map< int, std::uint64_t > lengths;
  checked( "neighbor receive lengths",
           [&]
           {
             for( std::size_t i = 0; i < m_sources.size(); ++i )
             {
               if( receiving[i] )
               {
                 lengths[m_sources[i]] = receiving[i];
               }
             }
           } );
  m_statistics.neighborExchanges = neighborExchanges;
  m_statistics.countBytesSent = countBytes;
  return exchangePayloads( outgoing, lengths );
}

void Communication::includeContactNeighbors( Participants neighbors )
{
  // The temporary vertex directory has already supplied symmetric adjacency.
  // Unioning two symmetric graphs preserves symmetry without another directory.
  checked( "contact neighbor union", [&] { neighbors.insert( neighbors.end(), m_neighbors.begin(), m_neighbors.end() ); } );
  initializeNeighbors( std::move( neighbors ) );
}

void Communication::validateVolumeFaces( std::vector< MainFace > const & faces, std::uint64_t mainNamespace )
{
  GEOS_MARK_SCOPE( "uniformRefinement/coarseVolumeFaceValidation" );
  using Counts = std::array< std::uint64_t, 2 >;
  std::unordered_map< EntityKey, Counts, EntityKeyHash > local;
  Mail outgoing;
  checked( "coarse owned volume faces",
           [&]
           {
             for( auto const & face : faces )
             {
               if( face.owners != Participants{ m_rank } )
               {
                 throw std::invalid_argument( "Nonlocal volume-face validation owner" );
               }
               EntityKey const key = entityKey( EntityKind::face, face.globalCorners, mainNamespace );
               validateKey( key );
               Connectivity directed = face.globalCorners;
               std::rotate( directed.begin(), std::min_element( directed.begin(), directed.end() ), directed.end() );
               auto & counts = local[key];
               auto & count = counts[directed == key.corners ? 0 : 1];
               count = checkedSum( count, 1 );
               if( counts[0] + counts[1] > 2 )
               {
                 throw std::invalid_argument( "Nonmanifold local volume face" );
               }
             }
             for( auto const & [key, counts] : local )
             {
               EntityKey routing = key;
               std::sort( routing.corners.begin(), routing.corners.end() );
               int const home = static_cast< int >( stableHash( routing ) % m_size );
               auto & bytes = outgoing[home];
               putKey( bytes, key );
               putInteger( bytes, counts[0] );
               putInteger( bytes, counts[1] );
             }
           } );
  auto const incoming = exchangeDirectory( outgoing );
  checked( "global coarse volume face incidence",
           [&]
           {
             struct Incidence
             {
               EntityKey cycle;
               Counts counts{};
             };
             std::unordered_map< EntityKey, Incidence, EntityKeyHash > incidence;
             for( auto const & [peer, bytes] : incoming )
             {
               GEOS_UNUSED_VAR( peer );
               Reader reader( bytes );
               while( !reader.done() )
               {
                 auto const key = reader.key();
                 if( key.kind != EntityKind::face || key.meshNamespace != mainNamespace )
                 {
                   throw std::invalid_argument( "Invalid volume face record" );
                 }
                 Counts const counts{ reader.integer(), reader.integer() };
                 if( counts[0] > 2 || counts[1] > 2 || counts[0] + counts[1] < 1 || counts[0] + counts[1] > 2 )
                 {
                   throw std::invalid_argument( "Invalid volume face incidence count" );
                 }
                 EntityKey support = key;
                 std::sort( support.corners.begin(), support.corners.end() );
                 auto const [found, inserted] = incidence.emplace( std::move( support ), Incidence{ key, {} } );
                 if( !inserted && !( found->second.cycle == key ) )
                 {
                   throw std::invalid_argument( "Crossed global volume face cycles" );
                 }
                 for( int orientation = 0; orientation < 2; ++orientation )
                 {
                   found->second.counts[orientation] = checkedSum( found->second.counts[orientation], counts[orientation] );
                 }
                 if( found->second.counts[0] + found->second.counts[1] > 2 )
                 {
                   throw std::invalid_argument( "Nonmanifold global volume face" );
                 }
               }
             }
             for( auto const & [key, record] : incidence )
             {
               GEOS_UNUSED_VAR( key );
               if( record.counts[0] + record.counts[1] == 2 && ( record.counts[0] != 1 || record.counts[1] != 1 ) )
               {
                 throw std::invalid_argument( "Incident volume faces have equal outward orientation" );
               }
             }
           } );
}

std::vector< std::vector< SurfaceSide > > Communication::discoverSurfaceSides( std::vector< MainFace > const & localFaces,
                                                                               std::vector< CoarseSurface > const & surfaces,
                                                                               std::uint64_t mainNamespace )
{
  GEOS_MARK_SCOPE( "uniformRefinement/coarseSurfaceAssociations" );
  std::vector< EntityKey > anchors;
  std::vector< Connectivity > queryCorners;
  std::vector< std::vector< SurfaceSide > > result;
  std::unique_ptr< SurfaceAssociations > associations;
  checked( "surface association input",
           [&]
           {
             associations = std::make_unique< SurfaceAssociations >( localFaces, mainNamespace );
             std::set< vtkIdType > vertices;
             for( auto const & face : localFaces )
             {
               if( face.owners != Participants{ m_rank } )
               {
                 throw std::invalid_argument( "Nonlocal coarse volume-face owner" );
               }
               vertices.insert( face.globalCorners.begin(), face.globalCorners.end() );
             }
             result.resize( surfaces.size() );
             queryCorners.resize( surfaces.size() );
             for( std::size_t q = 0; q < surfaces.size(); ++q )
             {
               auto const & buckets = surfaces[q].cornerBuckets;
               if( buckets.size() < 3 || buckets.size() > 11 )
               {
                 throw std::invalid_argument( "Unsupported coarse associated surface arity" );
               }
               for( auto const & bucket : buckets )
               {
                 if( bucket.empty() || !std::is_sorted( bucket.begin(), bucket.end() ) || bucket.front() < 0 ||
                     std::adjacent_find( bucket.begin(), bucket.end() ) != bucket.end() )
                 {
                   throw std::invalid_argument( "Surface association requires normalized nonempty buckets" );
                 }
               }
               auto & corners = queryCorners[q];
               for( std::size_t c = 0; c < buckets.size(); ++c )
               {
                 corners.push_back( static_cast< vtkIdType >( c ) );
               }
               result[q] = associations->match( corners, buckets );
               vertices.insert( buckets.front().begin(), buckets.front().end() );
             }
             for( vtkIdType id : vertices )
             {
               anchors.push_back( entityKey( EntityKind::vertex, { id }, mainNamespace ) );
             }
           } );
  // This separate coarse routing graph includes query-only ranks. It must never
  // become the main registry's participant map.
  Communication routing( m_comm, m_chunkBytes );
  Sharing const sharing = routing.discoverSharing( anchors );
  includeContactNeighbors( routing.neighbors() );
  Mail outgoing;
  checked( "surface association queries",
           [&]
           {
             for( std::size_t q = 0; q < surfaces.size(); ++q )
             {
               Participants peers;
               for( vtkIdType anchor : surfaces[q].cornerBuckets.front() )
               {
                 auto const & ranks = sharing.at( entityKey( EntityKind::vertex, { anchor }, mainNamespace ) );
                 peers.insert( peers.end(), ranks.begin(), ranks.end() );
               }
               std::sort( peers.begin(), peers.end() );
               peers.erase( std::unique( peers.begin(), peers.end() ), peers.end() );
               for( int peer : peers )
               {
                 if( peer == m_rank )
                 {
                   continue;
                 }
                 Bytes & bytes = outgoing[peer];
                 putInteger( bytes, q );
                 putInteger( bytes, surfaces[q].cornerBuckets.size() );
                 for( auto const & bucket : surfaces[q].cornerBuckets )
                 {
                   putInteger( bytes, bucket.size() );
                   for( vtkIdType id : bucket )
                   {
                     putInteger( bytes, id );
                   }
                 }
               }
             }
           } );
  auto const incoming = routing.exchangeNeighbors( outgoing );
  outgoing.clear();
  checked( "surface association replies",
           [&]
           {
             for( auto const & [peer, bytes] : incoming )
             {
               Reader reader( bytes );
               while( !reader.done() )
               {
                 auto const query = reader.integer(), n = reader.integer();
                 if( n < 3 || n > 11 )
                 {
                   throw std::invalid_argument( "Invalid surface query arity" );
                 }
                 std::vector< Connectivity > buckets( n );
                 Connectivity corners;
                 for( std::uint64_t c = 0; c < n; ++c )
                 {
                   corners.push_back( static_cast< vtkIdType >( c ) );
                   auto const count = reader.integer();
                   if( count == 0 || count > bytes.size() / 8 )
                   {
                     throw std::invalid_argument( "Invalid collocation query length" );
                   }
                   for( std::uint64_t i = 0; i < count; ++i )
                   {
                     buckets[c].push_back( reader.id() );
                   }
                 }
                 auto const sides = associations->match( corners, buckets );
                 Bytes & reply = outgoing[peer];
                 putInteger( reply, query );
                 putInteger( reply, sides.size() );
                 for( auto const & side : sides )
                 {
                   // Every local face must come from an owned volume on this rank.
                   if( side.owners != Participants{ m_rank } )
                   {
                     throw std::invalid_argument( "Nonlocal coarse volume-face owner" );
                   }
                   putInteger( reply, side.mainCornersBySurface.size() );
                   for( vtkIdType id : side.mainCornersBySurface )
                   {
                     putInteger( reply, id );
                   }
                 }
               }
             }
           } );
  auto const replies = routing.exchangeNeighbors( outgoing );
  checked( "surface association installation",
           [&]
           {
             for( auto const & [peer, bytes] : replies )
             {
               Reader reader( bytes );
               while( !reader.done() )
               {
                 auto const q = reader.integer(), count = reader.integer();
                 if( q >= result.size() || count > bytes.size() / 8 )
                 {
                   throw std::invalid_argument( "Invalid surface association reply" );
                 }
                 for( std::uint64_t s = 0; s < count; ++s )
                 {
                   auto const n = reader.integer();
                   if( n != surfaces[q].cornerBuckets.size() )
                   {
                     throw std::invalid_argument( "Surface side reply arity mismatch" );
                   }
                   Connectivity mapped;
                   for( std::uint64_t c = 0; c < n; ++c )
                   {
                     vtkIdType const id = reader.id();
                     auto const & bucket = surfaces[q].cornerBuckets[c];
                     if( !std::binary_search( bucket.begin(), bucket.end(), id ) )
                     {
                       throw std::invalid_argument( "Unrequested main side corner" );
                     }
                     mapped.push_back( id );
                   }
                   result[q].push_back( { entityKey( EntityKind::face, mapped, mainNamespace ), std::move( mapped ), { peer } } );
                 }
               }
             }
             for( auto & sides : result )
             {
               std::map< EntityKey, SurfaceSide > unique;
               for( auto & side : sides )
               {
                 auto const [found, inserted] = unique.emplace( side.mainFace, side );
                 if( !inserted )
                 {
                   if( found->second.mainCornersBySurface != side.mainCornersBySurface )
                   {
                     throw std::invalid_argument( "Inconsistent surface side correspondence" );
                   }
                   auto & owners = found->second.owners;
                   owners.insert( owners.end(), side.owners.begin(), side.owners.end() );
                   std::sort( owners.begin(), owners.end() );
                   owners.erase( std::unique( owners.begin(), owners.end() ), owners.end() );
                 }
               }
               if( unique.empty() )
               {
                 throw std::invalid_argument( "Surface cell has no actual incident volume face" );
               }
               sides.clear();
               for( auto & [key, side] : unique )
               {
                 GEOS_UNUSED_VAR( key );
                 sides.push_back( std::move( side ) );
               }
             }
             m_statistics.directoryExchanges = checkedSum( m_statistics.directoryExchanges, routing.m_statistics.directoryExchanges );
             m_statistics.neighborExchanges = checkedSum( m_statistics.neighborExchanges, routing.m_statistics.neighborExchanges );
             m_statistics.payloadChunksSent = checkedSum( m_statistics.payloadChunksSent, routing.m_statistics.payloadChunksSent );
             m_statistics.payloadBytesSent = checkedSum( m_statistics.payloadBytesSent, routing.m_statistics.payloadBytesSent );
             m_statistics.countBytesSent = checkedSum( m_statistics.countBytesSent, routing.m_statistics.countBytesSent );
           } );
  return result;
}

std::vector< vtkIdType > Communication::resolveSupportIds( std::uint64_t generation, std::vector< SupportLookup > const & requests,
                                                           std::unordered_map< EntityKey, vtkIdType, EntityKeyHash > const & localIds )
{
  GEOS_MARK_SCOPE( "uniformRefinement/surfaceSupportIds" );
  Mail outgoing;
  std::vector< vtkIdType > result;
  checked( "surface support requests",
           [&]
           {
             if( !m_discovered || generation == 0 )
             {
               throw std::invalid_argument( "Surface support lookup requires a positive generation and graph" );
             }
             result.assign( requests.size(), -1 );
             for( std::size_t q = 0; q < requests.size(); ++q )
             {
               auto const & request = requests[q];
               validateKey( request.key );
               if( request.key.kind == EntityKind::cell || request.faceOwner < 0 || request.faceOwner >= m_size )
               {
                 throw std::invalid_argument( "Invalid actual surface support request" );
               }
               if( request.faceOwner == m_rank )
               {
                 result[q] = localIds.at( request.key );
                 if( result[q] < 0 )
                 {
                   throw std::invalid_argument( "Negative local main support ID" );
                 }
               }
               else
               {
                 auto & bytes = outgoing[request.faceOwner];
                 putInteger( bytes, generation );
                 putInteger( bytes, q );
                 putKey( bytes, request.key );
               }
             }
           } );
  auto const incoming = exchangeNeighbors( outgoing );
  outgoing.clear();
  checked( "surface support replies",
           [&]
           {
             for( auto const & [peer, bytes] : incoming )
             {
               Reader reader( bytes );
               while( !reader.done() )
               {
                 if( reader.integer() != generation )
                 {
                   throw std::invalid_argument( "Surface support generation mismatch" );
                 }
                 auto const q = reader.integer();
                 auto const key = reader.key();
                 vtkIdType const id = localIds.at( key );
                 if( id < 0 )
                 {
                   throw std::invalid_argument( "Negative resolved main support ID" );
                 }
                 auto & reply = outgoing[peer];
                 putInteger( reply, generation );
                 putInteger( reply, q );
                 putKey( reply, key );
                 putInteger( reply, id );
               }
             }
           } );
  auto const replies = exchangeNeighbors( outgoing );
  checked( "surface support installation",
           [&]
           {
             for( auto const & [peer, bytes] : replies )
             {
               Reader reader( bytes );
               while( !reader.done() )
               {
                 if( reader.integer() != generation )
                 {
                   throw std::invalid_argument( "Surface support reply generation mismatch" );
                 }
                 auto const q = reader.integer();
                 auto const key = reader.key();
                 vtkIdType const id = reader.id();
                 if( q >= requests.size() || requests[q].faceOwner != peer || !( requests[q].key == key ) || result[q] != -1 )
                 {
                   throw std::invalid_argument( "Unexpected or duplicate surface support reply" );
                 }
                 result[q] = id;
               }
             }
             if( std::find( result.begin(), result.end(), -1 ) != result.end() )
             {
               throw std::invalid_argument( "Missing main surface support ID" );
             }
           } );
  return result;
}

Sharing Communication::discoverSharing( std::vector< EntityKey > const & entities )
{
  GEOS_MARK_SCOPE( "uniformRefinement/coarseDiscovery" );
  Mail outgoing;
  checked( "coarse entity keys",
           [&]
           {
             std::set< EntityKey > unique( entities.begin(), entities.end() );
             for( auto const & key : unique )
             {
               validateKey( key );
               int const home = static_cast< int >( stableHash( faceCornerSet( key ) ) % static_cast< std::uint64_t >( m_size ) );
               putKey( outgoing[home], key );
             }
           } );
  Mail const incoming = exchangeDirectory( outgoing );
  std::map< EntityKey, std::set< int > > directory;
  Mail responses;
  checked( "coarse entity directory",
           [&]
           {
             std::map< EntityKey, Connectivity > cycles;
             for( auto const & [peer, bytes] : incoming )
             {
               Reader reader( bytes );
               while( !reader.done() )
               {
                 auto const key = reader.key();
                 directory[key].insert( peer );
                 if( key.kind == EntityKind::face )
                 {
                   auto const [it, inserted] = cycles.emplace( faceCornerSet( key ), key.corners );
                   if( !inserted && it->second != key.corners )
                   {
                     throw std::invalid_argument( "Incompatible shared face edge cycles" );
                   }
                 }
               }
             }
             for( auto const & [key, ranks] : directory )
             {
               for( int r : ranks )
               {
                 auto & bytes = responses[r];
                 putKey( bytes, key );
                 putInteger( bytes, ranks.size() );
                 for( int participant : ranks )
                 {
                   putInteger( bytes, participant );
                 }
               }
             }
           } );
  Mail const replies = exchangeDirectory( responses );
  Sharing result;
  Participants neighbors;
  checked( "coarse directory replies",
           [&]
           {
             for( auto const & [peer, bytes] : replies )
             {
               GEOS_UNUSED_VAR( peer );
               Reader reader( bytes );
               while( !reader.done() )
               {
                 auto key = reader.key();
                 auto ranks = reader.participants( m_size );
                 if( !std::binary_search( ranks.begin(), ranks.end(), m_rank ) || !result.emplace( key, ranks ).second )
                 {
                   throw std::invalid_argument( "Invalid coarse refinement directory reply" );
                 }
                 neighbors.insert( neighbors.end(), ranks.begin(), ranks.end() );
               }
             }
             std::set< EntityKey > const expected( entities.begin(), entities.end() );
             if( result.size() != expected.size() )
             {
               throw std::invalid_argument( "Missing coarse refinement entity reply" );
             }
             for( auto const & key : expected )
             {
               if( !result.count( key ) )
               {
                 throw std::invalid_argument( "Unexpected coarse refinement directory reply" );
               }
             }
           } );
  initializeNeighbors( std::move( neighbors ) );
  return result;
}

IdRange Communication::allocateRange( std::uint64_t localCount, vtkIdType base ) const
{
  GEOS_MARK_SCOPE( "uniformRefinement/pointAndCellIds" );
  auto const minimum = MpiWrapper::allReduce( base, MpiWrapper::Reduction::Min, m_comm );
  auto const maximum = MpiWrapper::allReduce( base, MpiWrapper::Reduction::Max, m_comm );
  checked( "ID count validation",
           [&]
           {
             if( minimum != maximum || base < 0 )
             {
               throw std::invalid_argument( "Invalid collective refinement ID base" );
             }
           } );
  // Each 32-bit limb summed over at most INT_MAX ranks fits in uint64_t.
  // This detects even uint64_t overflow without O(P) allgather storage or
  // overflowing MPI_SUM before validation. No custom MPI reduction is needed.
  constexpr std::uint64_t mask = UINT64_C( 0xffffffff );
  std::array< std::uint64_t, 2 > const limbs{ localCount & mask, localCount >> 32 };
  std::array< std::uint64_t, 2 > sums{};
  MpiWrapper::allReduce( limbs, sums, MpiWrapper::Reduction::Sum, m_comm );
  std::uint64_t total = 0;
  checked( "ID range overflow",
           [&]
           {
             std::uint64_t const high = sums[1] + ( sums[0] >> 32 );
             if( high > mask )
             {
               throw std::overflow_error( "Refinement global count overflow" );
             }
             total = ( high << 32 ) | ( sums[0] & mask );
             if( total && total - 1 > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() - base ) )
             {
               throw std::overflow_error( "Refinement active ID range exceeds vtkIdType" );
             }
           } );
  std::uint64_t offset = 0;
  MpiWrapper::exscan( &localCount, &offset, 1, MPI_SUM, m_comm );
  if( m_rank == 0 )
  {
    offset = 0;
  }
  return { localCount ? static_cast< vtkIdType >( static_cast< std::uint64_t >( base ) + offset ) : 0, total };
}

std::map< EntityKey, PointRecord > Communication::resolvePoints( std::uint64_t generation, std::vector< PointCreation > const & points,
                                                                 vtkIdType localExistingMaximum )
{
  return resolvePointRecords( generation, points, localExistingMaximum, false );
}

std::map< EntityKey, PointRecord > Communication::reconcileExistingPoints( std::vector< PointCreation > const & points )
{
  return resolvePointRecords( 0, points, -1, true );
}

std::map< EntityKey, PointRecord > Communication::resolvePointRecords( std::uint64_t generation,
                                                                       std::vector< PointCreation > const & points,
                                                                       vtkIdType localExistingMaximum, bool existing )
{
  GEOS_MARK_SCOPE( "uniformRefinement/sharedEntityExchange" );
  auto const minimum = MpiWrapper::allReduce( generation, MpiWrapper::Reduction::Min, m_comm );
  auto const maximum = MpiWrapper::allReduce( generation, MpiWrapper::Reduction::Max, m_comm );
  std::map< EntityKey, PointCreation const * > expected;
  std::uint64_t allocated = 0;
  checked( "point creation requests",
           [&]
           {
             if( minimum != maximum || ( !existing && generation == 0 ) || localExistingMaximum < -1 )
             {
               throw std::invalid_argument( "Invalid refinement generation or point maximum" );
             }
             if( !m_discovered )
             {
               throw std::logic_error( "Refinement coarse sharing discovery must precede point resolution" );
             }
             for( auto const & point : points )
             {
               validateKey( point.key );
               if( existing ? point.key.kind != EntityKind::vertex : point.key.kind == EntityKind::vertex )
               {
                 throw std::invalid_argument( "Wrong entity kind for existing/new refinement point record" );
               }
               auto const & ranks = point.participants;
               if( ranks.empty() || !std::is_sorted( ranks.begin(), ranks.end() ) ||
                   std::adjacent_find( ranks.begin(), ranks.end() ) != ranks.end() || ranks.front() < 0 || ranks.back() >= m_size ||
                   !std::binary_search( ranks.begin(), ranks.end(), m_rank ) )
               {
                 throw std::invalid_argument( "Invalid new point participants" );
               }
               for( int r : ranks )
               {
                 if( r != m_rank && !std::binary_search( m_neighbors.begin(), m_neighbors.end(), r ) )
                 {
                   throw std::invalid_argument( "New point participant not discovered on coarse mesh" );
                 }
               }
               if( point.key.kind == EntityKind::cell && ranks.size() != 1 )
               {
                 throw std::invalid_argument( "Volume-cell interior point must have one participant" );
               }
               if( !std::isfinite( point.supportScale ) || point.supportScale <= 0 )
               {
                 throw std::invalid_argument( "Invalid shared point support extent" );
               }
               for( double x : point.position )
               {
                 if( !std::isfinite( x ) )
                 {
                   throw std::invalid_argument( "Nonfinite new point coordinates" );
                 }
               }
               if( !expected.emplace( point.key, &point ).second )
               {
                 throw std::invalid_argument( "Duplicate local point creation key" );
               }
               if( !existing && ranks.front() == m_rank )
               {
                 ++allocated;
               }
             }
           } );
  vtkIdType oldMaximum = -1;
  IdRange range{ 0, 0 };
  if( !existing )
  {
    oldMaximum = MpiWrapper::allReduce( localExistingMaximum, MpiWrapper::Reduction::Max, m_comm );
    std::uint64_t const anyNew = MpiWrapper::allReduce( allocated, MpiWrapper::Reduction::Max, m_comm );
    checked( "new point ID base",
             [&]
             {
               if( anyNew && oldMaximum == std::numeric_limits< vtkIdType >::max() )
               {
                 throw std::overflow_error( "Refinement max point ID + 1 overflow" );
               }
             } );
    range = allocateRange( allocated, anyNew ? oldMaximum + 1 : 0 );
  }
  Mail outgoing;
  std::map< EntityKey, PointRecord > records;
  checked( "authoritative point records",
           [&]
           {
             std::uint64_t ordinal = 0;
             for( auto const & [key, point] : expected )
             {
               if( point->participants.front() == m_rank )
               {
                 vtkIdType const id =
                     existing ? key.corners.front() : static_cast< vtkIdType >( static_cast< std::uint64_t >( range.first ) + ordinal++ );
                 PointRecord record{ id, point->position, point->fields };
                 records.emplace( key, record );
                 for( int r : point->participants )
                 {
                   if( r != m_rank )
                   {
                     auto & bytes = outgoing[r];
                     putInteger( bytes, generation );
                     putKey( bytes, key );
                     putInteger( bytes, point->participants.size() );
                     for( int participant : point->participants )
                     {
                       putInteger( bytes, participant );
                     }
                     putInteger( bytes, record.globalId );
                     for( double x : record.position )
                     {
                       putDouble( bytes, x );
                     }
                     putInteger( bytes, record.fields.size() );
                     bytes.insert( bytes.end(), record.fields.begin(), record.fields.end() );
                   }
                 }
               }
             }
           } );
  auto const incoming = exchangeNeighbors( outgoing );
  checked( "resolved point records",
           [&]
           {
             for( auto const & [peer, bytes] : incoming )
             {
               Reader reader( bytes );
               while( !reader.done() )
               {
                 if( reader.integer() != generation )
                 {
                   throw std::invalid_argument( "Wrong refinement generation in shared point record" );
                 }
                 EntityKey const key = reader.key();
                 auto const found = expected.find( key );
                 if( found == expected.end() || found->second->participants.front() != peer )
                 {
                   throw std::invalid_argument( "Unexpected shared point allocator/support" );
                 }
                 if( reader.participants( m_size ) != found->second->participants )
                 {
                   throw std::invalid_argument( "Inconsistent shared point participants" );
                 }
                 PointRecord record;
                 record.globalId = reader.id();
                 for( double & x : record.position )
                 {
                   x = reader.real();
                 }
                 record.fields = reader.payload();
                 if( ( existing ? record.globalId != key.corners.front() : record.globalId <= oldMaximum ) ||
                     !records.emplace( key, std::move( record ) ).second )
                 {
                   throw std::invalid_argument( "Repeated/invalid shared point ID record" );
                 }
               }
             }
             if( records.size() != expected.size() )
             {
               throw std::invalid_argument( "Missing shared point creation record" );
             }
             std::set< std::pair< std::uint64_t, vtkIdType > > ids;
             for( auto const & [key, record] : records )
             {
               if( !ids.emplace( key.meshNamespace, record.globalId ).second )
               {
                 throw std::invalid_argument( "Repeated local refinement point global ID" );
               }
               auto const & predicted = expected.at( key )->position;
               for( int d = 0; d < 3; ++d )
               {
                 if( std::abs( record.position[d] - predicted[d] ) >
                     64 * std::numeric_limits< double >::epsilon() * std::max( std::abs( record.position[d] ), std::abs( predicted[d] ) ) +
                         1e-12 * expected.at( key )->supportScale )
                 {
                   throw std::invalid_argument( "Shared point support coordinate mismatch" );
                 }
               }
             }
           } );
  return records;
}
std::map< ChildCellKey, CellRecord > Communication::resolveCells( std::uint64_t generation, std::vector< CellCreation > const & cells,
                                                                  vtkIdType base )
{
  return resolveCellRecords( generation, cells, base, false );
}

std::map< ChildCellKey, CellRecord > Communication::reconcileExistingCells( std::vector< CellCreation > const & cells )
{
  return resolveCellRecords( 0, cells, 0, true );
}

std::map< ChildCellKey, CellRecord > Communication::resolveCellRecords( std::uint64_t generation, std::vector< CellCreation > const & cells,
                                                                        vtkIdType base, bool existing )
{
  GEOS_MARK_SCOPE( "uniformRefinement/surfaceCellIds" );
  auto const minimum = MpiWrapper::allReduce( generation, MpiWrapper::Reduction::Min, m_comm );
  auto const maximum = MpiWrapper::allReduce( generation, MpiWrapper::Reduction::Max, m_comm );
  std::map< ChildCellKey, CellCreation const * > expected;
  std::uint64_t allocated = 0;
  checked( "surface cell requests",
           [&]
           {
             if( !m_discovered || ( !existing && generation == 0 ) || minimum != maximum )
             {
               throw std::invalid_argument( "Invalid surface refinement generation or missing coarse discovery" );
             }
             for( auto const & cell : cells )
             {
               auto const & key = cell.key;
               auto const & ranks = cell.participants;
               if( key.parentId < 0 || key.templateCode > 255 || key.ordinal >= 22 || ( existing && key.ordinal != 0 ) )
               {
                 throw std::invalid_argument( "Invalid surface child identity" );
               }
               if( ranks.empty() || !std::is_sorted( ranks.begin(), ranks.end() ) || ranks.front() < 0 || ranks.back() >= m_size ||
                   std::adjacent_find( ranks.begin(), ranks.end() ) != ranks.end() ||
                   !std::binary_search( ranks.begin(), ranks.end(), m_rank ) )
               {
                 throw std::invalid_argument( "Invalid surface child participants" );
               }
               for( int rank : ranks )
               {
                 if( rank != m_rank && !std::binary_search( m_neighbors.begin(), m_neighbors.end(), rank ) )
                 {
                   throw std::invalid_argument( "Surface child uses an undiscovered participant" );
                 }
               }
               if( !expected.emplace( key, &cell ).second )
               {
                 throw std::invalid_argument( "Duplicate local surface child identity" );
               }
               allocated += ranks.front() == m_rank;
             }
           } );
  auto const range = existing ? IdRange{ 0, 0 } : allocateRange( allocated, base );
  Mail outgoing;
  std::map< ChildCellKey, CellRecord > result;
  checked( "surface cell packing",
           [&]
           {
             std::uint64_t ordinal = 0;
             for( auto const & [key, cell] : expected )
             {
               if( cell->participants.front() == m_rank )
               {
                 CellRecord record{ existing ? key.parentId
                                             : static_cast< vtkIdType >( static_cast< std::uint64_t >( range.first ) + ordinal++ ),
                                    cell->fields };
                 result.emplace( key, record );
                 for( int rank : cell->participants )
                 {
                   if( rank != m_rank )
                   {
                     auto & bytes = outgoing[rank];
                     putInteger( bytes, generation );
                     putInteger( bytes, key.meshNamespace );
                     putInteger( bytes, key.parentId );
                     putInteger( bytes, key.templateCode );
                     putInteger( bytes, key.ordinal );
                     putInteger( bytes, cell->participants.size() );
                     for( int participant : cell->participants )
                     {
                       putInteger( bytes, participant );
                     }
                     putInteger( bytes, record.globalId );
                     putInteger( bytes, record.fields.size() );
                     bytes.insert( bytes.end(), record.fields.begin(), record.fields.end() );
                   }
                 }
               }
             }
           } );
  auto const incoming = exchangeNeighbors( outgoing );
  checked( "surface cell records",
           [&]
           {
             for( auto const & [peer, bytes] : incoming )
             {
               Reader reader( bytes );
               while( !reader.done() )
               {
                 if( reader.integer() != generation )
                 {
                   throw std::invalid_argument( "Wrong surface cell generation" );
                 }
                 ChildCellKey const key{ reader.integer(), reader.id(), reader.integer(), reader.integer() };
                 auto const found = expected.find( key );
                 if( found == expected.end() || found->second->participants.front() != peer )
                 {
                   throw std::invalid_argument( "Unexpected surface child allocator or identity" );
                 }
                 if( reader.participants( m_size ) != found->second->participants )
                 {
                   throw std::invalid_argument( "Inconsistent surface child participants" );
                 }
                 CellRecord record{ reader.id(), reader.payload() };
                 if( ( existing ? record.globalId != key.parentId : record.globalId < base ) ||
                     !result.emplace( key, std::move( record ) ).second )
                 {
                   throw std::invalid_argument( "Repeated or invalid surface child record" );
                 }
               }
             }
             if( result.size() != expected.size() )
             {
               throw std::invalid_argument( "Missing surface child record" );
             }
             std::set< std::pair< std::uint64_t, vtkIdType > > ids;
             for( auto const & [key, record] : result )
             {
               if( !ids.emplace( existing ? key.meshNamespace : 0, record.globalId ).second )
               {
                 throw std::invalid_argument( "Repeated allocated surface child ID" );
               }
             }
           } );
  return result;
}
} // namespace geos::vtk::refinement
