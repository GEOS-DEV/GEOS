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

#include "NewGhosting.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"

#include <NamedType/named_type.hpp>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDataSetSurfaceFilter.h>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <algorithm>
#include <utility>


namespace geos::ghosting
{

using NodeLocIdx = fluent::NamedType< localIndex, struct NodeLocIdxTag, fluent::Comparable, fluent::Printable >;
using NodeGlbIdx = fluent::NamedType< globalIndex, struct NodeGlbIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeLocIdx = fluent::NamedType< localIndex, struct EdgeLocIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable >;
using EdgeGlbIdx = fluent::NamedType< globalIndex, struct EdgeGlbIdxTag, fluent::Comparable, fluent::Printable >;
using FaceLocIdx = fluent::NamedType< localIndex, struct FaceLocIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable >;
using FaceGlbIdx = fluent::NamedType< globalIndex, struct FaceGlbIdxTag, fluent::Comparable, fluent::Printable >;
using CellLocIdx = fluent::NamedType< localIndex, struct CellLocIdxTag, fluent::Comparable, fluent::Printable >;
using CellGlbIdx = fluent::NamedType< globalIndex, struct CellGlbIdxTag, fluent::Comparable, fluent::Printable >;

using MpiRank = fluent::NamedType< int, struct MpiRankTag, fluent::Comparable, fluent::Printable, fluent::Addable >;


void to_json( json & j,
              const MpiRank & v )
{
  j = v.get();
}

void from_json( const json & j,
                MpiRank & v )
{
  v = MpiRank{ j.get< MpiRank::UnderlyingType >() };  // TODO use a `traits` instead
}

void to_json( json & j,
              const EdgeLocIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                EdgeLocIdx & v )
{
  v = EdgeLocIdx{ j.get< EdgeLocIdx::UnderlyingType >() };
}

void to_json( json & j,
              const FaceLocIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                FaceLocIdx & v )
{
  v = FaceLocIdx{ j.get< FaceLocIdx::UnderlyingType >() };
}

using Edge = std::tuple< NodeGlbIdx, NodeGlbIdx >;
using Face = std::vector< NodeGlbIdx >;

struct Exchange
{
  std::set< Edge > edges;
  std::set< Face > faces;
};

template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type *& buffer,
      Exchange const & exchange )
{
  localIndex sizeOfPackedChars = 0;
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, exchange.edges );
  sizeOfPackedChars += Pack< DO_PACKING >( buffer, exchange.faces );
  return sizeOfPackedChars;
}


std::set< vtkIdType > extractBoundaryCells( vtkSmartPointer< vtkDataSet > mesh )
{
  auto f = vtkDataSetSurfaceFilter::New();
  f->PassThroughCellIdsOn();
  f->PassThroughPointIdsOff();
  f->FastModeOff();

  string const originalCellsKey = "ORIGINAL_CELLS";
  f->SetOriginalCellIdsName( originalCellsKey.c_str() );
  auto boundaryMesh = vtkPolyData::New();
  f->UnstructuredGridExecute( mesh, boundaryMesh );
  vtkIdTypeArray const * originalCells = vtkIdTypeArray::FastDownCast( boundaryMesh->GetCellData()->GetArray( originalCellsKey.c_str() ) );

  std::set< vtkIdType > boundaryCellIdxs;
  for( auto i = 0; i < originalCells->GetNumberOfTuples(); ++i )
  {
    boundaryCellIdxs.insert( originalCells->GetValue( i ) );
  }

  return boundaryCellIdxs;
}

/**
 * @brief Order the nodes of the faces in a way that can be reproduced across the MPI ranks.
 * @param nodes The list of nodes as provided by the mesh.
 * @return A face with the nodes in the appropriate order
 * @details The nodes will be ordered in the following way.
 * First, we look for the lowest node index. It will become the first node.
 * Then we must pick the second node. We have to choices: just before or just after the first.
 * (we do not want to shuffle the nodes completely, we need to keep track of the order).
 * To do this, we select the nodes with the lowest index as well.
 * This also defines a direction in which we'll pick the other nodes.
 * For example, the face <tt>[2, 3, 1, 5, 9, 8]</tt> will become <tt>[1, 3, 2, 8, 9, 5]</tt>
 * because we'll start with @c 1 and then select the @c 3 over the @c 5.
 * Which defines the direction @c 2, @c 8, @c 9, @c 5.
 * @note This is the same pattern that we apply for edges.
 * Except that edges having only two nodes, it's not necessary to implement a dedicated function
 * and <tt>std::minmax</tt> is enough.
 */
Face reorderFaceNodes( std::vector< NodeGlbIdx > const & nodes )
{
  std::size_t const n = nodes.size();

  // Handles negative values of `i`.
  auto const modulo = [n]( integer const & i ) -> std::size_t
  {
    integer mod = i % n;
    if( mod < 0 )
    {
      mod += n;
    }
    return mod;
  };

  Face f;
  f.reserve( n );

  auto const it = std::min_element( nodes.cbegin(), nodes.cend() );
  std::size_t const minIdx = std::distance( nodes.cbegin(), it );
  int const increment = nodes[modulo( minIdx - 1 )] < nodes[modulo( minIdx + 1 )] ? -1 : 1;
  integer i = minIdx;
  for( std::size_t count = 0; count < n; ++count, i = i + increment )
  {
    f.emplace_back( nodes.at( modulo( i ) ) );
  }

  return f;
}

Exchange buildSpecificData( vtkSmartPointer< vtkDataSet > mesh,
                            std::set< vtkIdType > const & cellIds )
{
  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );

  std::vector< Edge > tmpEdges;
  std::vector< Face > tmpFaces;

  // Pre-allocation of the temporary vectors.
  {
    std::size_t numEdges = 0;
    std::size_t numFaces = 0;
    for( vtkIdType const & c: cellIds )
    {
      vtkCell * cell = mesh->GetCell( c );
      numEdges += cell->GetNumberOfEdges();
      numFaces += cell->GetNumberOfFaces();
    }
    tmpEdges.reserve( numEdges );
    tmpFaces.reserve( numFaces );
  }

  // Filling the temporary.
  for( auto c = 0; c < mesh->GetNumberOfCells(); ++c )
  {
    vtkCell * cell = mesh->GetCell( c );

    for( auto e = 0; e < cell->GetNumberOfEdges(); ++e )
    {
      vtkCell * edge = cell->GetEdge( e );
      vtkIdType const ln0 = edge->GetPointId( 0 );
      vtkIdType const ln1 = edge->GetPointId( 1 );
      vtkIdType const gn0 = globalPtIds->GetValue( ln0 );
      vtkIdType const gn1 = globalPtIds->GetValue( ln1 );
      tmpEdges.emplace_back( std::minmax( { NodeGlbIdx{ gn0 }, NodeGlbIdx{ gn1 } } ) );
    }

    for( auto f = 0; f < cell->GetNumberOfFaces(); ++f )
    {
      vtkCell * face = cell->GetFace( f );
      vtkIdList * pids = face->GetPointIds();
      std::vector< NodeGlbIdx > nodes( pids->GetNumberOfIds() );
      for( std::size_t i = 0; i < nodes.size(); ++i )
      {
        vtkIdType const lni = face->GetPointId( i );
        vtkIdType const gni = globalPtIds->GetValue( lni );
        nodes[i] = NodeGlbIdx{ gni };
      }
      tmpFaces.emplace_back( reorderFaceNodes( nodes ) );
    }
  }

  // Removing the duplicates by copying into a `std::set`.
  std::set< Edge > edges{ tmpEdges.cbegin(), tmpEdges.cend() };  // SortedArray requires the input to be sorted already.
  std::set< Face > faces{ tmpFaces.cbegin(), tmpFaces.cend() };

  return { std::move( edges ), std::move( faces ) };
}


array1d< globalIndex > convertExchange( Exchange const & exchange )
{
  std::size_t const edgeSize = 1 + 2 * std::size( exchange.edges );
  std::size_t faceSize = 1;
  for( Face const & face: exchange.faces )
  {
    faceSize += 1 + std::size( face );  // `+1` because we need to store the size so we know where to stop.
  }

  array1d< globalIndex > result;
  result.reserve( edgeSize + faceSize );

  result.emplace_back( std::size( exchange.edges ) );
  for( Edge const & edge: exchange.edges )
  {
    result.emplace_back( std::get< 0 >( edge ).get() );
    result.emplace_back( std::get< 1 >( edge ).get() );
  }
  result.emplace_back( std::size( exchange.faces ) );
  for( Face const & face: exchange.faces )
  {
    result.emplace_back( std::size( face ) );
    for( NodeGlbIdx const & n: face )
    {
      result.emplace_back( n.get() );
    }
  }

  return result;
}


Exchange convertExchange( array1d< globalIndex > const & input )
{
  Exchange exchange;

  globalIndex const numEdges = input[0];
  int i;
  for( i = 1; i < 2 * numEdges + 1; i += 2 )
  {
    NodeGlbIdx gn0{ input[i] }, gn1{ input[i + 1] };
    exchange.edges.insert( std::minmax( gn0, gn1 ) );
  }
  GEOS_ASSERT_EQ( std::size_t(numEdges), exchange.edges.size() );

  globalIndex const numFaces = input[i];
  for( ++i; i < input.size(); )
  {
    auto const s = input[i++];
    Face face( s );
    for( int j = 0; j < s; ++j, ++i )
    {
      face[j] = NodeGlbIdx{ input[i] };
    }
    exchange.faces.insert( face );
  }

  GEOS_ASSERT_EQ( std::size_t(numFaces), exchange.faces.size() );

  return exchange;
}


Exchange buildFullData( vtkSmartPointer< vtkDataSet > mesh )
{
  std::set< vtkIdType > cellIds;
  for( vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i )
  {
    cellIds.insert( cellIds.end(), i );
  }

  return buildSpecificData( mesh, cellIds );
}


Exchange buildExchangeData( vtkSmartPointer< vtkDataSet > mesh )
{
  return buildSpecificData( mesh, extractBoundaryCells( mesh ) );
}

struct Buckets
{
  std::map< std::set< MpiRank >, std::set< Edge > > edges;
  std::map< std::set< MpiRank >, std::set< Face > > faces;
};


/**
 * @brief
 * @param exchanged
 * @param neighbors excluding current rank
 */
Buckets buildIntersectionBuckets( std::map< MpiRank, Exchange > const & exchanged,
                                  MpiRank curRank,
                                  std::set< MpiRank > const & neighbors )
{
  std::map< Edge, std::set< MpiRank > > counts;  // TODO Use better intersection algorithms?
  // We "register" all the edges of the current rank: they are the only one we're interested in.
  for( Edge const & edge: exchanged.at( curRank ).edges )
  {
    counts.emplace_hint( counts.end(), edge, std::set< MpiRank >{ curRank } );
  }

  // We now loop on the neighbor edges.
  // If a neighbor has an edge in common with the current rank, they we store it.
  for( MpiRank const & neighborRank: neighbors )  // This does not include the current rank.
  {
    for( Edge const & edge: exchanged.at( neighborRank ).edges )
    {
      auto it = counts.find( edge );
      if( it != counts.cend() )  // TODO Extract `counts.cend()` out of the loop.
      {
        it->second.insert( neighborRank );
      }
    }
  }

  std::map< std::set< MpiRank >, std::set< Edge > > edgeBuckets;
  for( auto const & [edge, ranks]: counts )
  {
    if( ranks.find( curRank ) != ranks.cend() )
    {
      edgeBuckets[ranks].insert( edge );
    }
  }

  std::map< std::set< MpiRank >, std::set< Face > > faceBuckets;
  std::set< Face > curFaces = exchanged.at( curRank ).faces;
  for( MpiRank const & neighborRank: neighbors )  // This does not include the current rank.
  {
    for( Face const & face: exchanged.at( neighborRank ).faces )
    {
      auto it = curFaces.find( face );
      if( it != curFaces.cend() )
      {
        faceBuckets[{ curRank, neighborRank }].insert( *it );
        curFaces.erase( it );
      }
    }
  }
  faceBuckets[{ curRank }] = curFaces;

  // Checking if neighbors is too wide...  // TODO do we care?
  std::set< MpiRank > usefulNeighbors;
  for( auto const & [ranks, edges]: edgeBuckets )
  {
    if( not edges.empty() )
    {
      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
    }
  }
  std::vector< MpiRank > uselessNeighbors;
  std::set_difference( neighbors.cbegin(), neighbors.cend(), usefulNeighbors.cbegin(), usefulNeighbors.cend(), std::back_inserter( uselessNeighbors ) );
  // TODO... Remove the neighbors?

  return { edgeBuckets, faceBuckets };
}


// TODO Duplicated
std::map< MpiRank, Exchange > exchange( int commId,
                                        std::vector< NeighborCommunicator > & neighbors,
                                        Exchange const & data )
{
  MPI_iCommData commData( commId );
  integer const numNeighbors = LvArray::integerConversion< integer >( neighbors.size() );
  commData.resize( numNeighbors );


  array1d< globalIndex > const cv = convertExchange( data );

  for( integer i = 0; i < numNeighbors; ++i )
  {
    neighbors[i].mpiISendReceiveSizes( cv,
                                       commData.mpiSendBufferSizeRequest( i ),
                                       commData.mpiRecvBufferSizeRequest( i ),
                                       commId,
                                       MPI_COMM_GEOSX );
  }

  MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
  MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus() );

  array1d< array1d< globalIndex > > rawExchanged( neighbors.size() );

  for( integer i = 0; i < numNeighbors; ++i )
  {
    neighbors[i].mpiISendReceiveData( cv,
                                      commData.mpiSendBufferRequest( i ),
                                      rawExchanged[i],
                                      commData.mpiRecvBufferRequest( i ),
                                      commId,
                                      MPI_COMM_GEOSX );
  }
  MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );
  MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus() );

  std::map< MpiRank, Exchange > output;
  for( auto i = 0; i < numNeighbors; ++i )
  {
    output[MpiRank{ neighbors[i].neighborRank() }] = convertExchange( rawExchanged[i] );
  }
  return output;
}


std::size_t buildMaxBufferSize( Buckets const & buckets )
{
  auto const f = []( auto const & bucket ) -> std::size_t
  {
    std::size_t size = std::size( bucket );
    for( auto const & [ranks, _]: bucket )
    {
      size += std::size( ranks ) + 1 + 1;  // One `+1` for the size of the ranks set, the other one for the offset
    }
    return size;
  };

  return MpiWrapper::sum( f( buckets.edges ) + f( buckets.faces ) );  // Add the ScannedOffsets::{edge, face}Restart
}

struct BucketSizes
{
  using mapping = std::map< std::set< MpiRank >, localIndex >;
  mapping edges;
  mapping faces;
};

void to_json( json & j,
              const BucketSizes & v )
{
  j = json{ { "edges", v.edges }, { "faces", v.faces }  };
}

void from_json( const json & j,
                BucketSizes & v )
{
  v.edges = j.at( "edges" ).get< BucketSizes::mapping >();
  v.faces = j.at( "faces" ).get< BucketSizes::mapping >();
}

struct BucketOffsets
{
  std::map< std::set< MpiRank >, EdgeLocIdx > edges;
  std::map< std::set< MpiRank >, FaceLocIdx > faces;
};

void to_json( json & j,
              const BucketOffsets & v )
{
  j = json{ { "edges", v.edges }, { "faces", v.faces }  };
}

void from_json( const json & j,
                BucketOffsets & v )
{
  v.edges = j.at( "edges" ).get< std::map< std::set< MpiRank >, EdgeLocIdx > >();
  v.faces = j.at( "faces" ).get< std::map< std::set< MpiRank >, FaceLocIdx > >();
}

BucketSizes getBucketSize(Buckets const & buckets)
{
  BucketSizes output;

  for( auto const & [ranks, edges]: buckets.edges )
  {
    output.edges.emplace_hint( output.edges.end(), ranks, std::size( edges ) );
  }

  for( auto const & [ranks, faces]: buckets.faces )
  {
    output.faces.emplace_hint( output.faces.end(), ranks, std::size( faces ) );
  }

  return output;
}

/**
 * @brief
 * @tparam LOC_IDX The local index type of the geometrical quantity considered (typically @c EdgeLocIdx or @c FaceLocIdx).
 * @param sizes
 * @param offsets
 * @param curRank
 * @return
 */
template< typename LOC_IDX >
std::map< std::set< MpiRank >, LOC_IDX > updateBucketOffsets( std::map< std::set< MpiRank >, localIndex > const & sizes,
                                                              std::map< std::set< MpiRank >, LOC_IDX > const & offsets,
                                                              MpiRank curRank )
{
  std::map< std::set< MpiRank >, LOC_IDX > reducedOffsets;

  // Only consider the offsets that are still relevant (i.e. with high ranks)
  for( auto const & [ranks, offset]: offsets )
  {
    MpiRank const maxConcerned = *std::max_element( ranks.begin(), ranks.cend() );
    if( maxConcerned < curRank )
    {
      continue;
    }
    reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, offset );
  }

  // Add the offsets associated to the new buckets
  LOC_IDX nextOffset{ 0 };
  for( auto const & [ranks, size]: sizes )
  {
    auto const it = reducedOffsets.find( ranks );
    if( it == reducedOffsets.end() )
    {
      reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, nextOffset );
      nextOffset += LOC_IDX{ size };
    }
    else
    {
      nextOffset = it->second + LOC_IDX{ size };  // Define the new offset from the last
    }
  }

  // Add an extra entry based for the following rank
  reducedOffsets.emplace_hint( reducedOffsets.end(), std::set< MpiRank >{ curRank + MpiRank{ 1 } }, nextOffset );

  return reducedOffsets;
}

std::vector< std::uint8_t > serialize( BucketSizes const & sizes )
{
  return json::to_cbor( json( sizes ) );
}

std::vector< std::uint8_t > serialize( BucketOffsets const & offsets )
{
  return json::to_cbor( json( offsets ) );
}

/**
 * @brief
 * @tparam V Container of std::uint8_t
 * @param data
 * @return
 */
template< class V >
BucketOffsets deserialize( V const & data )
{
  return json::from_cbor( data, false ).template get< BucketOffsets >();
}

void f( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dptr )
{
  GEOS_ASSERT_EQ( *dptr, MPI_BYTE );

  MpiRank const curRank{ MpiWrapper::commRank() };

  // offsets provided by the previous rank(s)
  BucketOffsets const offsets = deserialize( Span< std::uint8_t >( (std::uint8_t *) in, *len ) );

  // Sizes provided by the current rank, under the form of a pointer to the data.
  // No need to serialize, since we're on the same rank.
  std::uintptr_t addr;
  std::memcpy( &addr, inout, sizeof( std::uintptr_t ) );
  BucketSizes const * sizes = reinterpret_cast<BucketSizes const *>(addr);

  BucketOffsets updatedOffsets;
  updatedOffsets.edges = updateBucketOffsets< EdgeLocIdx >( sizes->edges, offsets.edges, curRank );
  updatedOffsets.faces = updateBucketOffsets< FaceLocIdx >( sizes->faces, offsets.faces, curRank );

  // Serialize the updated offsets, so they get sent to the next rank.
  std::vector< std::uint8_t > const serialized = serialize( updatedOffsets );
  std::memcpy( inout, serialized.data(), serialized.size() );
}

void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< MpiRank > const & neighbors )
{
  // Now we exchange the data with our neighbors.
  MpiRank const curRank{ MpiWrapper::commRank() };

  std::vector< NeighborCommunicator > ncs;
  ncs.reserve( neighbors.size() );
  for( MpiRank const & rank: neighbors )
  {
    ncs.emplace_back( rank.get() );
  }

  CommID const commId = CommunicationTools::getInstance().getCommID();
  std::map< MpiRank, Exchange > exchanged = exchange( int( commId ), ncs, buildExchangeData( mesh ) );
  exchanged[MpiRank{ MpiWrapper::commRank() }] = buildFullData( mesh );

  Buckets const buckets = buildIntersectionBuckets( exchanged, curRank, neighbors );

  std::size_t const maxBufferSize = 100 * buildMaxBufferSize( buckets );

  std::vector< std::uint8_t > sendBuffer( maxBufferSize );
  std::vector< std::uint8_t > recvBuffer( maxBufferSize );

  BucketSizes const sizes = getBucketSize( buckets );
  BucketOffsets offsets;
  if( curRank == MpiRank{ 0 } )
  {
    // The `MPI_Scan` process will not call the reduction operator for rank 0.
    // So we need to reduce ourselves for ourselves.
    offsets.edges = updateBucketOffsets< EdgeLocIdx >( sizes.edges, { { { MpiRank{ 0 }, }, EdgeLocIdx{ 0 } } }, curRank );
    offsets.faces = updateBucketOffsets< FaceLocIdx >( sizes.faces, { { { MpiRank{ 0 }, }, FaceLocIdx{ 0 } } }, curRank );
    // Still we need to send this reduction to the following rank, by copying to it to the send buffer.
    std::vector< std::uint8_t > const bytes = serialize( offsets );
    std::memcpy( sendBuffer.data(), bytes.data(), bytes.size() );
  }
  else
  {
    // For the other ranks, the reduction operator will be called during the `Mpi_Scan` process.
    // So unlike for rank 0, we do not have to do it ourselves.
    // In order to provide the `sizes` to the reduction operator, since `sizes` will only be used on the current ranl,
    // we'll provide the information as a pointer to the instance.
    // The reduction operator will then compute the new offsets and send them to the following rank.
    std::uintptr_t const addr = reinterpret_cast<std::uintptr_t>(&sizes);
    std::memcpy( sendBuffer.data(), &addr, sizeof( std::uintptr_t ) );
  }

  MPI_Op op;
  MPI_Op_create( f, false, &op );

  MPI_Scan( sendBuffer.data(), recvBuffer.data(), maxBufferSize, MPI_BYTE, op, MPI_COMM_WORLD );

  if( curRank != MpiRank{ 0 } )
  {
    offsets = deserialize( recvBuffer );
  }

  std::cout << "offsets on rank " << curRank << " -> " << json( offsets ) << std::endl;
}

void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< int > const & neighbors )
{
  std::set< MpiRank > neighbors_;
  for( int const & rank: neighbors )
  {
    neighbors_.insert( MpiRank{ rank } );
  }

  return doTheNewGhosting( mesh, neighbors_ );
}

}  // end of namespace geos::ghosting
