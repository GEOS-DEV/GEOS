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
using EdgeLocIdx = fluent::NamedType< localIndex, struct EdgeLocIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeGlbIdx = fluent::NamedType< globalIndex, struct EdgeGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::PreIncrementable >;
using FaceLocIdx = fluent::NamedType< localIndex, struct FaceLocIdxTag, fluent::Comparable, fluent::Printable >;
using FaceGlbIdx = fluent::NamedType< globalIndex, struct FaceGlbIdxTag, fluent::Comparable, fluent::Printable, fluent::Addable, fluent::PreIncrementable >;
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

void from_json( const json & j,
                NodeGlbIdx & v )
{
  v = NodeGlbIdx{ j.get< NodeGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const NodeGlbIdx & v )
{
  j = v.get();
}

void to_json( json & j,
              const EdgeGlbIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                EdgeGlbIdx & v )
{
  v = EdgeGlbIdx{ j.get< EdgeGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const FaceGlbIdx & v )
{
  j = v.get();
}

void from_json( const json & j,
                FaceGlbIdx & v )
{
  v = FaceGlbIdx{ j.get< FaceGlbIdx::UnderlyingType >() };
}

void to_json( json & j,
              const CellGlbIdx & v )
{
  j = v.get();
}

using Edge = std::tuple< NodeGlbIdx, NodeGlbIdx >;
using Face = std::vector< NodeGlbIdx >;

struct Exchange
{
  std::set< NodeGlbIdx > nodes;
  std::set< Edge > edges;
  std::set< Face > faces;
};

void to_json( json & j,
              const Exchange & v )
{
  j = json{ { "nodes", v.nodes }, { "edges", v.edges }, { "faces", v.faces }  };
}

void from_json( const json & j,
                Exchange & v )
{
  v.nodes = j.at( "nodes" ).get< std::set< NodeGlbIdx > >();
  v.edges = j.at( "edges" ).get< std::set< Edge > >();
  v.faces = j.at( "faces" ).get< std::set< Face > >();
}

/**
 * @brief Extract the cells at the boundary of the mesh.
 * @param mesh The vtk mesh.
 * @return The vtk cell ids.
 */
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
  int const increment = nodes[modulo( minIdx - 1 )] < nodes[modulo( minIdx + 1 )] ? -1 : 1;  // TODO based on increment, I can say if the face is flipped or not.
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
  vtkIdTypeArray * gids = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );
  Span< vtkIdType > const globalPtIds( (vtkIdType *) gids->GetPointer( 0 ), gids->GetNumberOfTuples() );

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
      vtkIdType const gn0 = globalPtIds[ln0];
      vtkIdType const gn1 = globalPtIds[ln1];
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
        vtkIdType const gni = globalPtIds[lni];
        nodes[i] = NodeGlbIdx{ gni };
      }
      tmpFaces.emplace_back( reorderFaceNodes( nodes ) );
    }
  }

  std::set< NodeGlbIdx > nodes;
  std::transform( std::cbegin( globalPtIds ), std::cend( globalPtIds ),
                  std::inserter( nodes, std::end( nodes ) ), []( vtkIdType const & id )
                  { return NodeGlbIdx{ id }; } );

  // Removing the duplicates by copying into a `std::set`.
  std::set< Edge > edges{ tmpEdges.cbegin(), tmpEdges.cend() };  // SortedArray requires the input to be sorted already.
  std::set< Face > faces{ tmpFaces.cbegin(), tmpFaces.cend() };

  return { std::move( nodes ), std::move( edges ), std::move( faces ) };
}

array1d< std::uint8_t > convertExchange( Exchange const & exchange )
{
  std::vector< std::uint8_t > const tmp = json::to_cbor( json( exchange ) );
  array1d< std::uint8_t > result;
  result.reserve( tmp.size() );
  for( std::uint8_t const & t: tmp )
  {
    result.emplace_back( t );
  }
  return result;
}

Exchange convertExchange( array1d< std::uint8_t > const & input )
{
  std::vector< std::uint8_t > const tmp( std::cbegin( input ), std::cend( input ) );
  return json::from_cbor( tmp, false ).template get< Exchange >();
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
  std::map< std::set< MpiRank >, std::set< NodeGlbIdx > > nodes;
  std::map< std::set< MpiRank >, std::set< Edge > > edges;
  std::map< std::set< MpiRank >, std::set< Face > > faces;
};


/**
 * @brief Compute the intersection between for the ranks based on the information @p exchanged.
 * @param exchanged Geometrical information provided by the ranks (including the current rank).
 * @param curRank The current MPI rank.
 * @param neighbors excluding current rank
 * @return The intersection buckets.
 */
Buckets buildIntersectionBuckets( std::map< MpiRank, Exchange > const & exchanged,
                                  MpiRank curRank,
                                  std::set< MpiRank > const & neighbors )
{
  std::map< std::set< MpiRank >, std::set< NodeGlbIdx > > nodeBuckets;
  {  // Scope reduction // TODO DUPLICATED!
    std::map< NodeGlbIdx, std::set< MpiRank > > counts;  // TODO Use better intersection algorithms? + Duplicated
    // We "register" all the edges of the current rank: they are the only one we're interested in.
    for( NodeGlbIdx const & node: exchanged.at( curRank ).nodes )
    {
      counts.emplace_hint( counts.end(), node, std::set< MpiRank >{ curRank } );
    }

    // We now loop on the neighbor edges.
    // If a neighbor has an edge in common with the current rank, they we store it.
    for( MpiRank const & neighborRank: neighbors )  // This does not include the current rank.
    {
      for( NodeGlbIdx const & node: exchanged.at( neighborRank ).nodes )
      {
        auto it = counts.find( node );
        if( it != counts.cend() )  // TODO Extract `counts.cend()` out of the loop.
        {
          it->second.insert( neighborRank );
        }
      }
    }

    for( auto const & [node, ranks]: counts )
    {
      if( ranks.find( curRank ) != ranks.cend() )
      {
        nodeBuckets[ranks].insert( node );
      }
    }
  }

  std::map< std::set< MpiRank >, std::set< Edge > > edgeBuckets;
  {  // Scope reduction
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

    for( auto const & [edge, ranks]: counts )
    {
      if( ranks.find( curRank ) != ranks.cend() )
      {
        edgeBuckets[ranks].insert( edge );
      }
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
  for( auto const & [ranks, nodes]: nodeBuckets )
  {
    if( not nodes.empty() )
    {
      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
    }
  }
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

  return { nodeBuckets, edgeBuckets, faceBuckets };
}


// TODO Duplicated
std::map< MpiRank, Exchange > exchange( int commId,
                                        std::vector< NeighborCommunicator > & neighbors,
                                        Exchange const & data )
{
  MPI_iCommData commData( commId );
  integer const numNeighbors = LvArray::integerConversion< integer >( neighbors.size() );
  commData.resize( numNeighbors );


  array1d< std::uint8_t > const cv = convertExchange( data );

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

  array1d< array1d< std::uint8_t > > rawExchanged( neighbors.size() );

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


/**
 * @brief Estimate an upper bound of the serialized size of the size @p buckets.
 * @param buckets
 * @return Size in bytes.
 * @details @c MPI_Scan requires a fixed buffer size.
 * An a-priori estimation of the maximum size of the buffer leads to crazy amounts of memory
 * because of the combinatorial pattern of the rank intersections.
 * Therefore, we use an initial MPI communication to get a better estimation.
 * @node We don't care about the @c nodes because their global numbering is already done.
 */
std::size_t buildMaxBufferSize( Buckets const & buckets )  // TODO give the size buckets instead?
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
  std::map< std::set< MpiRank >, EdgeGlbIdx > edges;
  std::map< std::set< MpiRank >, FaceGlbIdx > faces;
};

void to_json( json & j,
              const BucketOffsets & v )
{
  j = json{ { "edges", v.edges }, { "faces", v.faces }  };
}

void from_json( const json & j,
                BucketOffsets & v )
{
  v.edges = j.at( "edges" ).get< std::map< std::set< MpiRank >, EdgeGlbIdx > >();
  v.faces = j.at( "faces" ).get< std::map< std::set< MpiRank >, FaceGlbIdx > >();
}

/**
 * @brief Compute the sizes of the intersection buckets.
 * @param buckets The intersection buckets.
 * @return The buckets of sizes.
 */
BucketSizes getBucketSize( Buckets const & buckets )
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
 * @tparam GLB_IDX The global index type of the geometrical quantity considered (typically @c EdgeGlbIdx or @c FaceGlbIdx).
 * @param sizes The sizes of the intersection bucket (from the current MPI rank).
 * @param offsets The offsets for each intersection bucket (from the MPI_Scan process, i.e. the previous ranks).
 * @param curRank Current MPI rank.
 * @return The offset buckets, updated with the sizes of the current MPI rank.
 */
template< typename GLB_IDX >
std::map< std::set< MpiRank >, GLB_IDX > updateBucketOffsets( std::map< std::set< MpiRank >, localIndex > const & sizes,
                                                              std::map< std::set< MpiRank >, GLB_IDX > const & offsets,
                                                              MpiRank curRank )
{
  std::map< std::set< MpiRank >, GLB_IDX > reducedOffsets;

  // We only keep the offsets that are still relevant to the current and higher ranks.
  // Note that `reducedOffsets` will be used by the _current_ rank too.
  // Therefore, we'll send information that may not be useful to higher ranks.
  // This is surely acceptable because the information is tiny.
  for( auto const & [ranks, offset]: offsets )
  {
    MpiRank const maxConcerned = *std::max_element( ranks.cbegin(), ranks.cend() );
    // TODO invert
//    if( maxConcerned >= curRank )
//    {
//      reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, offset );
//    }
    if( maxConcerned < curRank )
    {
      continue;
    }
    reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, offset );
  }

  // Add the offsets associated to the new size buckets from the current rank.
  GLB_IDX nextOffset{ 0 };
  for( auto const & [ranks, size]: sizes )
  {
    auto const it = reducedOffsets.find( ranks );
    if( it == reducedOffsets.end() )
    {
      reducedOffsets.emplace_hint( reducedOffsets.end(), ranks, nextOffset );
      nextOffset += GLB_IDX{ size };
    }
    else
    {
      nextOffset = it->second + GLB_IDX{ size };  // Define the new offset from the last
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

/**
 * @brief Custom MPI reduction function. Merges the sizes of the bucket from current rank into numbering offsets.
 * @param[in] in Contains the reduced result from the previous ranks. Needs to be unpacked into a @c BucketOffsets instance.
 * @param[inout] inout As an @e input, contains a pointer to the @c BucketSizes instance from the current rank.
 * As an @e output, will contained the (packed) instance of @c BucketOffsets after the @c reduction.
 * The @e output instance will be used by both current and next ranks.
 * @param[in] len Number of data.
 * @param[in] dataType Type of data. Mean to be @c MPI_BYTE for the current case.
 */
void f( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dataType )
{
  GEOS_ASSERT_EQ( *dataType, MPI_BYTE );

  MpiRank const curRank{ MpiWrapper::commRank() };

  // offsets provided by the previous rank(s)
  BucketOffsets const offsets = deserialize( Span< std::uint8_t >( (std::uint8_t *) in, *len ) );

  // Sizes provided by the current rank, under the form of a pointer to the data.
  // No need to serialize, since we're on the same rank.
  std::uintptr_t addr;
  std::memcpy( &addr, inout, sizeof( std::uintptr_t ) );
  BucketSizes const * sizes = reinterpret_cast<BucketSizes const *>(addr);

  BucketOffsets updatedOffsets;
  updatedOffsets.edges = updateBucketOffsets< EdgeGlbIdx >( sizes->edges, offsets.edges, curRank );
  updatedOffsets.faces = updateBucketOffsets< FaceGlbIdx >( sizes->faces, offsets.faces, curRank );

  // Serialize the updated offsets, so they get sent to the next rank.
  std::vector< std::uint8_t > const serialized = serialize( updatedOffsets );
  std::memcpy( inout, serialized.data(), serialized.size() );
}

struct MaxGlbIdcs
{
  NodeGlbIdx nodes;
  EdgeGlbIdx edges;
  FaceGlbIdx faces;
  CellGlbIdx cells;
};

void to_json( json & j,
              const MaxGlbIdcs & mo )
{
  j = json{ { "nodes", mo.nodes },
            { "edges", mo.edges },
            { "faces", mo.faces },
            { "cells", mo.cells } };
}

void g( void * in,
        void * inout,
        int * len,
        MPI_Datatype * dataType )
{
  GEOS_ASSERT_EQ( *len, 1 );

  MaxGlbIdcs const * i = reinterpret_cast<MaxGlbIdcs const *>(in);
  MaxGlbIdcs * io = reinterpret_cast<MaxGlbIdcs *>(inout);

  io->nodes = std::max( i->nodes, io->nodes );
  io->edges = std::max( i->edges, io->edges );
  io->faces = std::max( i->faces, io->faces );
  io->cells = std::max( i->cells, io->cells );
}

MaxGlbIdcs gatherOffset( vtkSmartPointer< vtkDataSet > mesh,
                         EdgeGlbIdx const & maxEdgeId,
                         FaceGlbIdx const & maxFaceId )
{
  MaxGlbIdcs offsets{ NodeGlbIdx{ 0 }, maxEdgeId, maxFaceId, CellGlbIdx{ 0 } };

  auto const extract = []( vtkDataArray * globalIds ) -> vtkIdType
  {
    vtkIdTypeArray * gids = vtkIdTypeArray::FastDownCast( globalIds );
    Span< vtkIdType > const s( (vtkIdType *) gids->GetPointer( 0 ), gids->GetNumberOfTuples() );
    return *std::max_element( s.begin(), s.end() );
  };

  offsets.nodes = NodeGlbIdx{ extract( mesh->GetPointData()->GetGlobalIds() ) };
  offsets.cells = CellGlbIdx{ extract( mesh->GetCellData()->GetGlobalIds() ) };

  // Otherwise, use `MPI_Type_create_struct`.
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, EdgeGlbIdx::UnderlyingType > );
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, FaceGlbIdx::UnderlyingType > );
  static_assert( std::is_same_v< NodeGlbIdx::UnderlyingType, CellGlbIdx::UnderlyingType > );

//  MPI_Datatype const underlying = internal::getMpiType< NodeGlbIdx::UnderlyingType >();
  MPI_Datatype t;
//  MPI_Type_contiguous( sizeof( MatrixOffsets ) / sizeof( underlying ), underlying, &t );
  MPI_Type_contiguous( 4, MPI_LONG_LONG_INT, &t );
  MPI_Type_commit( &t );

  MPI_Op op;
  MPI_Op_create( g, true, &op );

  MaxGlbIdcs result( offsets );

  MPI_Allreduce( &offsets, &result, 1, t, op, MPI_COMM_WORLD );

  return result;
}

struct MeshGraph  // TODO add the local <-> global mappings here?
{
  std::map< CellGlbIdx, std::set< FaceGlbIdx > > c2f;  // TODO What about the metadata (e.g. flip the face)
  std::map< FaceGlbIdx, std::set< EdgeGlbIdx > > f2e;  // TODO use Face here?
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n; // TODO use Edge here?
};

MeshGraph buildMeshGraph( vtkSmartPointer< vtkDataSet > mesh,  // TODO give a sub-mesh?
                          Buckets const & buckets,
                          BucketSizes const & sizes,
                          BucketOffsets const & offsets )
{
  MeshGraph result;

  for( auto const & [ranks, size]: sizes.edges )
  {
    EdgeGlbIdx i{ offsets.edges.at( ranks ) };
    for( Edge const & edge: buckets.edges.at( ranks ) )
    {
      result.e2n[i] = edge;
      ++i;
    }
  }

  // Simple inversion
  std::map< std::tuple< NodeGlbIdx, NodeGlbIdx >, EdgeGlbIdx > n2e;
  for( auto const & [e, n]: result.e2n )
  {
    n2e[n] = e;
  }

  std::map< std::vector< NodeGlbIdx >, FaceGlbIdx > n2f;
  for( auto const & [ranks, size]: sizes.faces )
  {
    FaceGlbIdx i{ offsets.faces.at( ranks ) };
    for( Face face: buckets.faces.at( ranks ) )  // Intentional copy for the future `emplace_back`.
    {
      n2f[face] = i;
      face.emplace_back( face.front() );
      for( std::size_t ii = 0; ii < face.size() - 1; ++ii )
      {
        NodeGlbIdx const & n0 = face[ii], & n1 = face[ii + 1];
        std::pair< NodeGlbIdx, NodeGlbIdx > const p0 = std::make_pair( n0, n1 );
        std::pair< NodeGlbIdx, NodeGlbIdx > const p1 = std::minmax( n0, n1 );
        result.f2e[i].insert( n2e.at( p1 ) );
        bool const flipped = p0 != p1;  // TODO store somewhere.
      }
      ++i;
    }
  }

  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );
  vtkIdTypeArray const * globalCellIds = vtkIdTypeArray::FastDownCast( mesh->GetCellData()->GetGlobalIds() );  // TODO do the mapping beforehand
  for( vtkIdType c = 0; c < mesh->GetNumberOfCells(); ++c )
  {
    vtkCell * cell = mesh->GetCell( c );
    CellGlbIdx const gci{ globalCellIds->GetValue( c ) };
    // TODO copy paste?
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
      std::vector< NodeGlbIdx > const reorderedNodes = reorderFaceNodes( nodes );
      result.c2f[gci].insert( n2f.at( reorderedNodes ) );
      // TODO... bool const flipped = ... compare nodes and reorderedNodes. Or ask `reorderFaceNodes` to tell
    }
  }

  return result;
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

//  std::cout << "exchanged on rank " << curRank << " -> " << json( exchanged[MpiRank{ MpiWrapper::commRank() }] ) << std::endl;

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
    offsets.edges = updateBucketOffsets< EdgeGlbIdx >( sizes.edges, { { { MpiRank{ 0 }, }, EdgeGlbIdx { 0 } } }, curRank );
    offsets.faces = updateBucketOffsets< FaceGlbIdx >( sizes.faces, { { { MpiRank{ 0 }, }, FaceGlbIdx { 0 } } }, curRank );
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

//  std::cout << "offsets on rank " << curRank << " -> " << json( offsets ) << std::endl;

  MpiRank const nextRank = curRank + MpiRank{ 1 };
  MaxGlbIdcs const matrixOffsets = gatherOffset( mesh, offsets.edges.at( { nextRank } ), offsets.faces.at( { nextRank } ) );
  std::cout << "matrixOffsets on rank " << curRank << " -> " << json( matrixOffsets ) << std::endl;

  MeshGraph const graph = buildMeshGraph( mesh, buckets, sizes, offsets );
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
