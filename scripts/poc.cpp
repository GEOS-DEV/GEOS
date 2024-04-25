/*
 * Compile: mpicxx -std=c++20 -I /tmp/json/include scan.cpp -o scan
 * Launch: mpirun -n 5 scan
 */

#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <span>
#include <string>
#include <vector>

std::map< std::set< int >, int > get_bucket_sizes( int rank )
{
  // for rank, intersection in enumerate(intersections):
  //         print(f" --> {rank}")
  //         for r, i in sorted(intersection.items()):
  //                 print(f"{{{{{', '.join(map(str, r))}}}, {len(i)}}},")

  switch( rank )
  {
    case 0:
      return {
        { { 0 },          799 },
        { { 0, 1 },       70 },
        { { 0, 1, 2 },    6 },
        { { 0, 1, 2, 3 }, 1 },
        { { 0, 1, 3 },    6 },
        { { 0, 2 },       118 },
        { { 0, 2, 3 },    6 },
        { { 0, 3 },       47 } };
    case 1:
      return {
        { { 0, 1 },       70 },
        { { 0, 1, 2 },    6 },
        { { 0, 1, 2, 3 }, 1 },
        { { 0, 1, 3 },    6 },
        { { 1 },          781 },
        { { 1, 2 },       45 },
        { { 1, 2, 3 },    5 },
        { { 1, 3 },       106 } };
    case 2:
      return {
        { { 0, 1, 2 },    6 },
        { { 0, 1, 2, 3 }, 1 },
        { { 0, 2 },       118 },
        { { 0, 2, 3 },    6 },
        { { 1, 2 },       45 },
        { { 1, 2, 3 },    5 },
        { { 2 },          769 },
        { { 2, 3 },       72 } };
    case 3:
      return {
        { { 0, 1, 2, 3 }, 1 },
        { { 0, 1, 3 },    6 },
        { { 0, 2, 3 },    6 },
        { { 0, 3 },       47 },
        { { 1, 2, 3 },    5 },
        { { 1, 3 },       106 },
        { { 2, 3 },       72 },
        { { 3 },          799 } };
    default:
      std::cout << "WRONG!!!" << std::endl;
      return {};
  };
}

/**
 * `rank` the current rank
 */
std::map< std::set< int >, int > update_bucket_offsets( std::map< std::set< int >, int > const & sizes,
                                                        std::map< std::set< int >, int > const & offsets,
                                                        int rank )
{
  std::map< std::set< int >, int > reduced;

  // Only consider the offsets that are still relevant (i.e. with high ranks)
  for( auto const & [ranks, offset]: offsets )
  {
    int const max_concerned = *std::max_element( ranks.begin(), ranks.cend() );
    if( max_concerned < rank )
    {
      continue;
    }
    reduced.emplace_hint( reduced.end(), ranks, offset );
  }

  // Add the offsets associated to the new buckets
  int next_offset = 0;
  for( auto const & [ranks, size]: sizes )
  {
    auto const it = reduced.find( ranks );
    if( it == reduced.end() )
    {
      reduced.emplace_hint( reduced.end(), ranks, next_offset );
      next_offset += size;
    }
    else
    {
      next_offset = it->second + size;  // Define the new offset from the last
    }
  }

  // Add an extra entry based for the following rank
  reduced.emplace_hint( reduced.end(), std::set< int >{ rank + 1 }, next_offset );

  return reduced;
}

template< class V >
std::map< std::set< int >, int > deserialize( V const & s )
{
  json const j = json::from_cbor( s, false );

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  return j.get< std::map< std::set< int >, int>>();
}

std::map< std::set< int >, int > deserialize( std::uint8_t const * data,
                                              int len )
{
  const std::span< const std::uint8_t > s( data, len );
  return deserialize( s );
}

std::vector< std::uint8_t > serialize( std::map< std::set< int >, int > const & in )
{
  return json::to_cbor( json( in ) );
}

/**
 * `in` contains the reduced result from the previous ranks
 * `inout` contains the input from the current rank and must overwritten with the new result.
 */
void f( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dptr )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // offsets provided by the previous rank(s)
  std::map< std::set< int >, int > offsets = deserialize( reinterpret_cast<std::uint8_t const *>(in), *len );

  // Sizes provided by the current rank, under the form of a pointer to the data.
  // No need to serialize, since we're on the same rank.
  std::uintptr_t addr;
  std::memcpy( &addr, inout, sizeof( std::uintptr_t ) );
  std::map< std::set< int >, int > const * sizes = reinterpret_cast<std::map< std::set< int >, int > *>(addr);

  std::map< std::set< int >, int > updated_offsets = update_bucket_offsets( *sizes, offsets, rank );

  // Serialize the updated offsets, so they get sent to the next rank.
  std::vector< std::uint8_t > const serialized = serialize( updated_offsets );
  std::memcpy( inout, serialized.data(), serialized.size() );
}

int main( int argc,
          char ** argv )
{
  int rank, size;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  constexpr int N = 2048;

  std::vector< std::uint8_t > sendBuff( N, 0 );
  std::vector< std::uint8_t > recvBuff( N, 0 );

  MPI_Op myOp;
  MPI_Op_create( f, false, &myOp );

  // For the rank 0, the `MPI_Scan` will not call the reduction operator.
  // So we need to reduce ourselves. Still we need to send this reduction to the following rank,
  // by copying to it to the send buffer.
  //
  // For the other ranks, the reduction operator will be called.
  // We'll then provide the data as a pointer to the instance on the current rank.
  // The reduction operator will then update the offsets and send them to the following rank.
  std::map< std::set< int >, int > const sizes = get_bucket_sizes( rank );
  std::map< std::set< int >, int > result;
  if( rank == 0 )
  {
    result = update_bucket_offsets( sizes, { { { 0, }, 0 } }, rank );
    std::vector< std::uint8_t > const bytes = serialize( result );
    std::memcpy( sendBuff.data(), bytes.data(), bytes.size() );
  }
  else
  {
    std::uintptr_t const addr = reinterpret_cast<std::uintptr_t>(&sizes);
    std::memcpy( sendBuff.data(), &addr, sizeof( std::uintptr_t ) );
  }

  MPI_Scan( sendBuff.data(), recvBuff.data(), N, MPI_BYTE, myOp, MPI_COMM_WORLD );

  if( rank != 0 )
  {
    result = deserialize( recvBuff );
  }

  std::cout << "on rank " << rank << " FINAL offsets -> " << json( result ) << std::endl;

  // Just be a little in order
  MPI_Barrier( MPI_COMM_WORLD );

  MPI_Finalize();

  return 0;
}
