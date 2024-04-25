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
  bool offroad = false;
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

std::map< std::set< int >, int > unserialize( void const * data,
                                              int len )
{

  std::uint8_t const * data0 = reinterpret_cast<std::uint8_t const *>(data);

  const std::span< const std::uint8_t > s( data0, len );
  json const j = json::from_cbor( s, false );

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  return j.get< std::map< std::set< int >, int>>();
}

std::vector< std::uint8_t > serialize( std::map< std::set< int >, int > const & in )
{
  return json::to_cbor( json( in ) );
}

/**
 * `inout` contains the input from the current rank and must overwritten with the new result.
 * `in` contains the reduced result from the previous ranks
 */
void f( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dptr )
{
  std::map< std::set< int >, int > offsets = unserialize( in, *len );  // offsets provided by the previous rank(s)
  std::map< std::set< int >, int > sizes = unserialize( inout, *len );  // sizes of the buckets provided by the current rank

  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  std::map< std::set< int >, int > updated_offsets = update_bucket_offsets( sizes, offsets, rank );


  std::vector< std::uint8_t > const serialized = serialize( updated_offsets );
  std::memcpy( inout, serialized.data(), *len );
}

int main( int argc,
          char ** argv )
{
  constexpr int N = 2048;

  int rank, size;
  std::uint8_t * local = (std::uint8_t *) ( malloc( N * sizeof( std::uint8_t ) ) );
  std::uint8_t * recv = (std::uint8_t *) ( malloc( N * sizeof( std::uint8_t ) ) );

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  MPI_Op myOp;
  MPI_Op_create( f, false, &myOp );

  // Input
  std::map< std::set< int >, int > sizes = get_bucket_sizes( rank );
  if( rank == 0 )
  {
    sizes = update_bucket_offsets( sizes, { { { 0, }, 0 } }, rank );
  }
  std::vector< std::uint8_t > const mm = serialize( sizes );
  std::memcpy( local, mm.data(), mm.size() );

  MPI_Scan( local, recv, N, MPI_BYTE, myOp, MPI_COMM_WORLD );

  std::map< std::set< int >, int > loc = unserialize( local, N );
  std::map< std::set< int >, int > rec = unserialize( recv, N );
  std::cout << "on rank " << rank << " FINAL local, recv -> " << json( loc ) << " | " << json( rec ) << std::endl;

  // Just be a little in order
  MPI_Barrier( MPI_COMM_WORLD );

  MPI_Finalize();

  free( local );
  free( recv );

  return 0;
}
