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

#include "Pods.hpp"

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"

#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_VectorOut.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_RowMatrixTransposer.h>
#include <Epetra_Vector.h>

#include "Indices.hpp"

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
  j = json{ { "nodes", v.nodes },
            { "edges", v.edges },
            { "faces", v.faces } };
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
  // TODO Better handle the boundary information, forgetting about the 3d cells and simply handling the outside shell.
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
        vtkIdType const nli = face->GetPointId( i );
        vtkIdType const ngi = globalPtIds[nli];
        nodes[i] = NodeGlbIdx{ ngi };
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

//void to_json( json & j,
//              const Buckets & v )
//{
//  j = json{ { "nodes", v.nodes }, { "edges", v.edges }, { "faces", v.faces }  };
//}

/**
 * @brief
 * @tparam T Typically NodeGloIdx or a container of NodeGlbIdx (for edges and faces)
 * @param exchanged
 * @param curRank
 * @param neighbors
 * @return
 */
template< typename T >
std::map< std::set< MpiRank >, std::set< T > >
buildIntersectionBuckets( std::map< MpiRank, std::set< T > const & > const & exchanged,
                          MpiRank curRank,
                          std::set< MpiRank > const & neighbors )
{
  std::map< T, std::set< MpiRank > > counts;  // TODO Use better intersection algorithms?
  // We "register" all the edges of the current rank: they are the only one we're interested in.
  for( T const & node: exchanged.at( curRank ) )
  {
    counts.emplace_hint( counts.end(), node, std::set< MpiRank >{ curRank } );
  }

  // We now loop on the neighbor edges.
  // If a neighbor has an edge in common with the current rank, they we store it.
  for( MpiRank const & neighborRank: neighbors )  // This does not include the current rank.
  {
    for( T const & node: exchanged.at( neighborRank ) )
    {
      auto it = counts.find( node );
      if( it != counts.cend() )  // TODO Extract `counts.cend()` out of the loop.
      {
        it->second.insert( neighborRank );
      }
    }
  }

  std::map< std::set< MpiRank >, std::set< T > > nodeBuckets;
  for( auto const & [node, ranks]: counts )
  {
    if( ranks.find( curRank ) != ranks.cend() )
    {
      nodeBuckets[ranks].insert( node );
    }
  }
  return nodeBuckets;
}


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
  std::map< MpiRank, std::set< NodeGlbIdx > const & > nodeInfo;
  for( auto const & [rank, exchange]: exchanged )
  {
    nodeInfo.emplace( rank, exchange.nodes );
  }
  std::map< std::set< MpiRank >, std::set< NodeGlbIdx > > const nodeBuckets = buildIntersectionBuckets( nodeInfo, curRank, neighbors );

  std::map< MpiRank, std::set< Edge > const & > edgeInfo;
  for( auto const & [rank, exchange]: exchanged )
  {
    edgeInfo.emplace( rank, exchange.edges );
  }
  std::map< std::set< MpiRank >, std::set< Edge > > const edgeBuckets = buildIntersectionBuckets( edgeInfo, curRank, neighbors );

  // For faces, the algorithm can be a tad simpler because faces can be shared by at most 2 ranks.
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

//  // Checking if neighbors is too wide...  // TODO do we care here?
//  std::set< MpiRank > usefulNeighbors;
//  for( auto const & [ranks, nodes]: nodeBuckets )
//  {
//    if( not nodes.empty() )
//    {
//      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
//    }
//  }
//  for( auto const & [ranks, edges]: edgeBuckets )
//  {
//    if( not edges.empty() )
//    {
//      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
//    }
//  }
//  std::vector< MpiRank > uselessNeighbors;
//  std::set_difference( neighbors.cbegin(), neighbors.cend(), usefulNeighbors.cbegin(), usefulNeighbors.cend(), std::back_inserter( uselessNeighbors ) );
//  // TODO... Remove the neighbors?

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
  j = json{ { "edges", v.edges },
            { "faces", v.faces } };
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
  j = json{ { "edges", v.edges },
            { "faces", v.faces } };
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
  reducedOffsets.emplace_hint( reducedOffsets.end(), std::set< MpiRank >{ curRank + 1_mpi }, nextOffset );

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
  std::map< FaceGlbIdx, std::set< EdgeGlbIdx > > f2e;
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n; // TODO use Edge here?
  std::set< NodeGlbIdx > nodes;
  std::set< FaceGlbIdx > otherFaces;  // Faces that are there but not owned.
  std::set< EdgeGlbIdx > otherEdges;  // Edges that are there but not owned.
  std::set< NodeGlbIdx > otherNodes;  // Nodes that are there but not owned.
  // TODO add the nodes here?
  // TODO add all types of connections here? How?
};

void to_json( json & j,
              const MeshGraph & v )  // For display
{
  j = json{ { "c2f", v.f2e },
            { "f2e", v.f2e },
            { "e2n", v.e2n } };
}


/**
 * @brief Builds the graph information for the owned elements only.
 * @param mesh
 * @param buckets
 * @param offsets
 * @param curRank
 * @return
 */
MeshGraph buildMeshGraph( vtkSmartPointer< vtkDataSet > mesh,  // TODO give a sub-mesh?
                          Buckets const & buckets,
                          BucketOffsets const & offsets,
                          MpiRank curRank )
{
  MeshGraph result;

  auto const isCurrentRankOwning = [&curRank]( std::set< MpiRank > const & ranks ) -> bool
  {
    return curRank == *std::min_element( std::cbegin( ranks ), std::cend( ranks ) );
  };

  for( auto const & [ranks, ns]: buckets.nodes )
  {
    std::set< NodeGlbIdx > & nodes = isCurrentRankOwning( ranks ) ? result.nodes : result.otherNodes;
    nodes.insert( std::cbegin( ns ), std::cend( ns ) );
  }

  // The `e2n` is a mapping for all the geometrical entities, not only the one owned like `result.e2n`.
  // TODO check that it is really useful.
  std::map< EdgeGlbIdx, std::tuple< NodeGlbIdx, NodeGlbIdx > > e2n;
  for( auto const & [ranks, edges]: buckets.edges )
  {
    bool const isOwning = isCurrentRankOwning( ranks );
    auto & m = isOwning ? result.e2n : e2n;
    EdgeGlbIdx i = offsets.edges.at( ranks );  // TODO hack
    for( Edge const & edge: edges )
    {
      m[i] = edge;
      if( not isOwning )
      {
        result.otherEdges.insert( i );  // TODO use the keys of e2n instead?
      }
      ++i;
    }
  }
  e2n.insert( std::cbegin( result.e2n ), std::cend( result.e2n ) );

  // Simple inversion
  std::map< std::tuple< NodeGlbIdx, NodeGlbIdx >, EdgeGlbIdx > n2e;
  for( auto const & [e, n]: e2n )
  {
    n2e[n] = e;  // TODO what about ownership?
  }

  for( auto const & [ranks, faces]: buckets.faces )
  {
    bool const isOwning = isCurrentRankOwning( ranks );

    if( isOwning )
    {
      FaceGlbIdx i = offsets.faces.at( ranks );
      for( Face face: faces )  // Intentional copy for the future `emplace_back`.
      {
        face.emplace_back( face.front() );  // Trick to build the edges.
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
    else
    {
      FaceGlbIdx const size = FaceGlbIdx{ intConv< FaceGlbIdx::UnderlyingType >( std::size( faces ) ) };
      for( FaceGlbIdx ii = offsets.faces.at( ranks ); ii < size; ++ii )
      {
        result.otherFaces.insert( ii );  // TODO insert iota
      }
    }
  }

  std::map< std::vector< NodeGlbIdx >, FaceGlbIdx > n2f;
  for( auto const & [ranks, faces]: buckets.faces )
  {
    FaceGlbIdx i = offsets.faces.at( ranks );
    for( Face const & face: faces )  // TODO hack
    {
      n2f[face] = i;
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
      std::vector< NodeGlbIdx > faceNodes( pids->GetNumberOfIds() );
      for( std::size_t i = 0; i < faceNodes.size(); ++i )
      {
        vtkIdType const nli = face->GetPointId( i );
        vtkIdType const ngi = globalPtIds->GetValue( nli );
        faceNodes[i] = NodeGlbIdx{ ngi };
      }
      std::vector< NodeGlbIdx > const reorderedFaceNodes = reorderFaceNodes( faceNodes );
      result.c2f[gci].insert( n2f.at( reorderedFaceNodes ) );
      // TODO... bool const flipped = ... compare nodes and reorderedNodes. Or ask `reorderFaceNodes` to tell
    }
  }

  return result;
}

std::unique_ptr< Epetra_CrsMatrix > makeTranspose( Epetra_CrsMatrix & input,
                                                   bool makeDataContiguous = true )
{
  Epetra_RowMatrixTransposer transposer( &input );  // This process does not modify the original `ghosted` matrix.
  Epetra_CrsMatrix * tr = nullptr;  // The transposer returns a pointer we must handle ourselves.
  transposer.CreateTranspose( makeDataContiguous, tr );

  std::unique_ptr< Epetra_CrsMatrix > ptr;
  ptr.reset( tr );
  return ptr;
}

/**
 * Contains the full result of the ghosting
 */
struct Ghost
{
  std::map< NodeGlbIdx, MpiRank > nodes;
  std::map< EdgeGlbIdx, MpiRank > edges;
  std::map< FaceGlbIdx, MpiRank > faces;
  std::map< CellGlbIdx, MpiRank > cells;
  std::set< MpiRank > neighbors;
};


Epetra_CrsMatrix multiply( int commSize,
                           Epetra_CrsMatrix const & indicator,
                           Epetra_CrsMatrix & upward )
{
  Epetra_Map const & ownedMap = upward.RowMap();
  Epetra_Map const & mpiMap = indicator.RangeMap();

  // Upward (n -> e -> f -> c)

  Epetra_CrsMatrix result_u0_0( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( upward, false, indicator, true, result_u0_0, false );
  result_u0_0.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-0.mat", result_u0_0 );

  Epetra_CrsMatrix result_u0_1( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( upward, false, result_u0_0, false, result_u0_1, false );
  result_u0_1.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-1.mat", result_u0_1 );

  Epetra_CrsMatrix result_u0_2( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( upward, false, result_u0_1, false, result_u0_2, false );
  result_u0_2.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-2.mat", result_u0_2 );

  // Downward (c -> f -> e -> n)
  auto tDownward = makeTranspose( upward );  // TODO check the algorithm to understand what's more relevant.

  Epetra_CrsMatrix result_d0_0( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( *tDownward, false, result_u0_2, false, result_d0_0, false );
  result_d0_0.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-4.mat", result_d0_0 );

  Epetra_CrsMatrix result_d0_1( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( *tDownward, false, result_d0_0, false, result_d0_1, false );
  result_d0_1.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-5.mat", result_d0_1 );

  Epetra_CrsMatrix result_d0_2( Epetra_DataAccess::Copy, ownedMap, commSize, false );
  EpetraExt::MatrixMatrix::Multiply( *tDownward, false, result_d0_1, false, result_d0_2, false );
  result_d0_2.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/result-6.mat", result_d0_2 );

  return result_d0_2;
}

class FindGeometricalType
{
public:
  enum Geom
  {
    NODE,
    EDGE,
    FACE,
    CELL
  };

  FindGeometricalType( std::size_t const edgeOffset,
                       std::size_t const faceOffset,
                       std::size_t const cellOffset )
    : m_edgeOffset( intConv< int >( edgeOffset ) ),
      m_faceOffset( intConv< int >( faceOffset ) ),
      m_cellOffset( intConv< int >( cellOffset ) )
  { }

  Geom operator()( int const & index ) const
  {
    if( index < m_edgeOffset )
    {
      return Geom::NODE;
    }
    else if( index < m_faceOffset )
    {
      return Geom::EDGE;
    }
    else if( index < m_cellOffset )
    {
      return Geom::FACE;
    }
    else
    {
      return Geom::CELL;
    }
  }

  template< typename GLOBAL_INDEX >
  [[nodiscard]] GLOBAL_INDEX as( int const & index ) const
  {
    if constexpr( std::is_same_v< GLOBAL_INDEX, NodeGlbIdx > )
    {
      return NodeGlbIdx{ intConv< NodeGlbIdx::UnderlyingType >( index ) };
    }
    else if constexpr( std::is_same_v< GLOBAL_INDEX, EdgeGlbIdx > )
    {
      return EdgeGlbIdx{ intConv< EdgeGlbIdx::UnderlyingType >( index - m_edgeOffset ) };
    }
    else if constexpr( std::is_same_v< GLOBAL_INDEX, FaceGlbIdx > )
    {
      return FaceGlbIdx{ intConv< FaceGlbIdx::UnderlyingType >( index - m_faceOffset ) };
    }
    else if constexpr( std::is_same_v< GLOBAL_INDEX, CellGlbIdx > )
    {
      return CellGlbIdx{ intConv< CellGlbIdx::UnderlyingType >( index - m_cellOffset ) };
    }
    else
    {
      GEOS_ERROR( "Internal Error." );
    }
  }

private:
  int const m_edgeOffset;
  int const m_faceOffset;
  int const m_cellOffset;
};


Ghost assembleAdjacencyMatrix( MeshGraph const & graph,
                               MaxGlbIdcs const & gis,
                               MpiRank curRank )
{
  std::size_t const edgeOffset = gis.nodes.get() + 1;
  std::size_t const faceOffset = edgeOffset + gis.edges.get() + 1;
  std::size_t const cellOffset = faceOffset + gis.faces.get() + 1;

  std::size_t const n = cellOffset + gis.cells.get() + 1;  // Total number of entries in the graph.

  std::size_t const numOwnedNodes = std::size( graph.nodes );
  std::size_t const numOwnedEdges = std::size( graph.e2n );
  std::size_t const numOwnedFaces = std::size( graph.f2e );
  std::size_t const numOwnedCells = std::size( graph.c2f );
  std::size_t const numOwned = numOwnedNodes + numOwnedEdges + numOwnedFaces + numOwnedCells;

  std::size_t const numOtherNodes = std::size( graph.otherNodes );
  std::size_t const numOtherEdges = std::size( graph.otherEdges );
  std::size_t const numOtherFaces = std::size( graph.otherFaces );
  std::size_t const numOther = numOtherNodes + numOtherEdges + numOtherFaces;

  std::vector< int > ownedGlbIdcs, numEntriesPerRow;  // TODO I couldn't use a vector of `std::size_t`
  std::vector< std::vector< int > > indices;
  ownedGlbIdcs.reserve( numOwned );
  numEntriesPerRow.reserve( numOwned );
  indices.reserve( numOwned );

  std::vector< int > otherGlbIdcs;  // TODO I couldn't use a vector of `std::size_t`
  otherGlbIdcs.reserve( numOther );
  for( NodeGlbIdx const & ngi: graph.otherNodes )
  {
    otherGlbIdcs.emplace_back( ngi.get() );
  }
  for( EdgeGlbIdx const & egi: graph.otherEdges )
  {
    otherGlbIdcs.emplace_back( egi.get() + edgeOffset );
  }
  for( FaceGlbIdx const & fgi: graph.otherFaces )
  {
    otherGlbIdcs.emplace_back( fgi.get() + faceOffset );
  }
  GEOS_ASSERT_EQ( numOther, std::size( otherGlbIdcs ) );

  for( NodeGlbIdx const & ngi: graph.nodes )
  {
    auto const i = ngi.get();
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( 1 );
    std::vector< int > const tmp( 1, ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }
  for( auto const & [egi, nodes]: graph.e2n )
  {
    auto const i = egi.get() + edgeOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::tuple_size_v< decltype( nodes ) > + 1 );  // `+1` comes from the identity
    std::vector< int > const tmp{ int( std::get< 0 >( nodes ).get() ), int( std::get< 1 >( nodes ).get() ), ownedGlbIdcs.back() };
    indices.emplace_back( tmp );
  }
  for( auto const & [fgi, edges]: graph.f2e )
  {
    auto const i = fgi.get() + faceOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::size( edges ) + 1 );  // `+1` comes from the identity
    std::vector< int > tmp;
    tmp.reserve( numEntriesPerRow.back() );
    for( EdgeGlbIdx const & egi: edges )
    {
      tmp.emplace_back( egi.get() + edgeOffset );
    }
    tmp.emplace_back( ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }
  for( auto const & [cgi, faces]: graph.c2f )
  {
    auto const i = cgi.get() + cellOffset;
    ownedGlbIdcs.emplace_back( i );
    numEntriesPerRow.emplace_back( std::size( faces ) + 1 );  // `+1` comes from the identity
    std::vector< int > tmp;
    tmp.reserve( numEntriesPerRow.back() );
    for( FaceGlbIdx const & fgi: faces )
    {
      tmp.emplace_back( fgi.get() + faceOffset );
    }
    tmp.emplace_back( ownedGlbIdcs.back() );
    indices.emplace_back( tmp );
  }

  GEOS_ASSERT_EQ( numOwned, std::size( ownedGlbIdcs ) );
  GEOS_ASSERT_EQ( numOwned, std::size( numEntriesPerRow ) );
  GEOS_ASSERT_EQ( numOwned, std::size( indices ) );
  for( std::size_t i = 0; i < numOwned; ++i )
  {
    GEOS_ASSERT_EQ( indices[i].size(), std::size_t( numEntriesPerRow[i] ) );
  }

  Epetra_MpiComm const & comm = Epetra_MpiComm( MPI_COMM_GEOSX );
  Epetra_Map const ownedMap( n, numOwned, ownedGlbIdcs.data(), 0, comm );

  Epetra_CrsMatrix upward( Epetra_DataAccess::Copy, ownedMap, numEntriesPerRow.data(), true );

  for( std::size_t i = 0; i < numOwned; ++i )
  {
    std::vector< int > const & rowIndices = indices[i];
    std::vector< double > const rowValues( std::size( rowIndices ), 1. );
    GEOS_ASSERT_EQ( std::size( rowIndices ), std::size_t( numEntriesPerRow[i] ) );
    upward.InsertGlobalValues( ownedGlbIdcs[i], std::size( rowIndices ), rowValues.data(), rowIndices.data() );
  }

//  upward.FillComplete();
  upward.FillComplete( ownedMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/adj.mat", upward );

  int const commSize( MpiWrapper::commSize() );

  // Now let's build the domain indicator matrix.
  // It's rectangular, one dimension being the number of MPI ranks, the other the number of nodes in the mesh graph.
  Epetra_Map const mpiMap( commSize, 0, comm );  // Let the current rank get the appropriate index in this map.
  Epetra_CrsMatrix indicator( Epetra_DataAccess::Copy, mpiMap, numOwned + numOther, true );

  std::vector< double > const ones( n, 1. );
  indicator.InsertGlobalValues( curRank.get(), numOwned, ones.data(), ownedGlbIdcs.data() );
  indicator.InsertGlobalValues( curRank.get(), numOther, ones.data(), otherGlbIdcs.data() );
  indicator.FillComplete( ownedMap, mpiMap );

  // TODO Could we use an Epetra_Vector as a diagonal matrix?
  std::vector< double > myRank( 1, curRank.get() );
  Epetra_CrsMatrix ownership( Epetra_DataAccess::Copy, ownedMap, 1, true );
  for( auto const & i: ownedGlbIdcs )
  {
    ownership.InsertGlobalValues( i, 1, myRank.data(), &i );
  }
  ownership.FillComplete( ownedMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ownership.mat", ownership );

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "indicator.NumGlobalCols() = " << indicator.NumGlobalCols() );
    GEOS_LOG_RANK( "indicator.NumGlobalRows() = " << indicator.NumGlobalRows() );
    GEOS_LOG_RANK( "ownership.NumGlobalCols() = " << ownership.NumGlobalCols() );
    GEOS_LOG_RANK( "ownership.NumGlobalRows() = " << ownership.NumGlobalRows() );
    GEOS_LOG_RANK( "ownership diag = " << std::boolalpha << ownership.LowerTriangular() and ownership.UpperTriangular() );
  }
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/indicator.mat", indicator );

  Epetra_CrsMatrix ghosted( multiply( commSize, indicator, upward ) );
  ghosted.PutScalar( 1. );  // This can be done after the `FillComplete`.
  ghosted.FillComplete( mpiMap, ownedMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ghosted.mat", ghosted );

  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghosted->NumGlobalCols() = " << ghosted.NumGlobalCols() );
    GEOS_LOG_RANK( "ghosted->NumGlobalRows() = " << ghosted.NumGlobalRows() );
  }
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/ghosted-after.mat", ghosted );

  Epetra_CrsMatrix ghostInfo( Epetra_DataAccess::Copy, mpiMap, 1, false );
  EpetraExt::MatrixMatrix::Multiply( ghosted, true, ownership, false, ghostInfo, false );
  ghostInfo.FillComplete( ownedMap, mpiMap );
  if( curRank == 0_mpi )
  {
    GEOS_LOG_RANK( "ghostInfo->NumGlobalCols() = " << ghostInfo.NumGlobalCols() );
    GEOS_LOG_RANK( "ghostInfo->NumGlobalRows() = " << ghostInfo.NumGlobalRows() );
  }

  Epetra_CrsMatrix identity( Epetra_DataAccess::Copy, ownedMap, 1, true );
  for( int const & i: ownedGlbIdcs )
  {
    identity.InsertGlobalValues( i, 1, ones.data(), &i );
  }
  identity.FillComplete();

  auto tDownward = makeTranspose( upward );  // TODO give it to multiply!
//  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/tDownward-before.mat", *tDownward );
  Epetra_Vector const zeros( ownedMap, true );
  tDownward->ReplaceDiagonalValues( zeros );  // TODO directly build zeros in the function call.
//  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/tDownward-after.mat", *tDownward );

  int extracted = 0;
  std::vector< double > extractedValues( n );
  std::vector< int > extractedIndices( n );
  ghostInfo.ExtractGlobalRowCopy( curRank.get(), int( n ), extracted, extractedValues.data(), extractedIndices.data() );
  extractedValues.resize( extracted );
  extractedIndices.resize( extracted );
  GEOS_LOG_RANK( "extracted = " << extracted );

  std::vector< int > missingIndices;
  missingIndices.reserve( n );

//  int extractedEdges = 0;
//  std::vector< double > extractedEdgesValues( n );
//  std::vector< int > extractedEdgesIndices( n );
  FindGeometricalType const getGeomType( edgeOffset, faceOffset, cellOffset );
  using Geom = FindGeometricalType::Geom;

  Ghost ghost;
  for( int i = 0; i < extracted; ++i )  // TODO Can we `zip` instead?
  {
    int const & index = extractedIndices[i];

    MpiRank const owner = MpiRank{ int( extractedValues[i] ) };
    if( owner != curRank )
    {
      ghost.neighbors.insert( owner );
    }

    switch( getGeomType( index ) )
    {
      case Geom::NODE:
      {
        NodeGlbIdx const ngi = getGeomType.as< NodeGlbIdx >( index );
        ghost.nodes.emplace( ngi, owner );
        break;
      }
      case Geom::EDGE:
      {
        EdgeGlbIdx const egi = getGeomType.as< EdgeGlbIdx >( index );
        ghost.edges.emplace( egi, owner );
        if( graph.e2n.find( egi ) == graph.e2n.cend() and graph.otherEdges.find( egi ) == graph.otherEdges.cend() )
        {
          missingIndices.emplace_back( index );
        }
        break;
      }
      case Geom::FACE:
      {
        FaceGlbIdx const fgi = getGeomType.as< FaceGlbIdx >( index );
        ghost.faces.emplace( fgi, owner );
        if( graph.f2e.find( fgi ) == graph.f2e.cend() and graph.otherFaces.find( fgi ) == graph.otherFaces.cend() )
        {
          missingIndices.emplace_back( index );
        }
        break;
      }
      case Geom::CELL:
      {
        CellGlbIdx const cgi = getGeomType.as< CellGlbIdx >( index );
        ghost.cells.emplace( cgi, owner );
        if( graph.c2f.find( cgi ) == graph.c2f.cend() )
        {
          missingIndices.emplace_back( index );
        }
        break;
      }
      default:
      {
        GEOS_ERROR( "Internal error" );
      }
    }
  }

  GEOS_LOG( "Gathering the missing mappings." );

  std::size_t const numMissingIndices = std::size( missingIndices );
  std::size_t const numGlobalMissingIndices = MpiWrapper::sum( numMissingIndices );

//  Epetra_Map const missingIndicesMap( intConv< long long int >( numGlobalMissingIndices ), intConv< int >( numMissingIndices ), 0, comm );
  Epetra_Map const missingIndicesMap( intConv< int >( numGlobalMissingIndices ), intConv< int >( numMissingIndices ), 0, comm );

//  std::size_t const offset = MpiWrapper::prefixSum< std::size_t >( numMissingIndices );
  int const offset = *missingIndicesMap.MyGlobalElements();
//  GEOS_LOG_RANK( "numMissingIndices, numGlobalMissingIndices, offset = " << numMissingIndices << ", " << numGlobalMissingIndices << ", " << offset );

  Epetra_CrsMatrix missing( Epetra_DataAccess::Copy, missingIndicesMap, 1, true );
  for( std::size_t i = 0; i < numMissingIndices; ++i )
  {
    missing.InsertGlobalValues( offset + i, 1, ones.data(), &missingIndices[i] );
  }
  missing.FillComplete( ownedMap, missingIndicesMap );

  // TODO Don't put ones in downard: reassemble with the ordering (e.g. edges point to a first and last node)
  Epetra_CrsMatrix missingMappings( Epetra_DataAccess::Copy, missingIndicesMap, 1, false );
  EpetraExt::MatrixMatrix::Multiply( missing, false, *tDownward, true, missingMappings, false );
  missingMappings.FillComplete( ownedMap, missingIndicesMap );
  EpetraExt::RowMatrixToMatrixMarketFile( "/tmp/matrices/missingMappings.mat", missingMappings );

  MeshGraph neighboringConnexions;

  int ext = 0;
  std::vector< double > extValues( n );
  std::vector< int > extIndices( n );
//  double * extValues = nullptr;
//  int * extIndices = nullptr;

  for( int i = 0; i < int( numMissingIndices ); ++i )
  {
    missingMappings.ExtractGlobalRowCopy( offset + i, int( n ), ext, extValues.data(), extIndices.data() );
//    missingMappings.ExtractGlobalRowView( intConv< int >( offset + i ), ext, extValues );
//    missingMappings.ExtractGlobalRowView( intConv< int >( offset + i ), ext, extValues, extIndices );
//    Span< double > const s0( extValues, ext );
//    Span< int > const s1( extIndices, ext );
//    GEOS_LOG( "ext 0 = " << std::size( s0 ) );
//    GEOS_LOG( "ext 1 = " << std::size( s1 ) );
//    GEOS_LOG_RANK( "ext = " << ext );
    int const index = missingIndices[i];
    switch( getGeomType( index ) )
    {
      case Geom::EDGE:
      {
        GEOS_ASSERT_EQ( ext, 2 );  // TODO check that val != 0?
        EdgeGlbIdx const egi = getGeomType.as< EdgeGlbIdx >( index );
        NodeGlbIdx const ngi0 = getGeomType.as< NodeGlbIdx >( extIndices[0] );
        NodeGlbIdx const ngi1 = getGeomType.as< NodeGlbIdx >( extIndices[1] );
        neighboringConnexions.e2n[egi] = std::minmax( { ngi0, ngi1 } );
        break;
      }
      case Geom::FACE:
      {
        GEOS_ASSERT_EQ( ext, 4 );  // TODO temporary check for the dev
        FaceGlbIdx const fgi = getGeomType.as< FaceGlbIdx >( index );
        std::set< EdgeGlbIdx > & tmp = neighboringConnexions.f2e[fgi];
        for( int ii = 0; ii < ext; ++ii )
        {
          tmp.insert( EdgeGlbIdx{ intConv< globalIndex >( extIndices[ii] - edgeOffset ) } );
        }
        break;
      }
      case Geom::CELL:
      {
        GEOS_ASSERT_EQ( ext, 6 );  // TODO temporary check for the dev
        CellGlbIdx const cgi = getGeomType.as< CellGlbIdx >( index );
        std::set< FaceGlbIdx > & tmp = neighboringConnexions.c2f[cgi];
        for( int ii = 0; ii < ext; ++ii )
        {
          tmp.insert( FaceGlbIdx{ intConv< globalIndex >( extIndices[ii] - faceOffset ) } );
        }
        break;
      }
      default:
      {
        GEOS_ERROR( "Internal error." );
      }
    }
  }

  if( curRank == 2_mpi )
  {
    GEOS_LOG_RANK( "neighboringConnexions = " << json( neighboringConnexions ) );
    std::map< MpiRank, std::set< CellGlbIdx > > tmp;
    for( auto const & [geom, rk]: ghost.cells )
    {
      tmp[rk].insert( geom );
    }
    for( auto const & [rk, geom]: tmp )
    {
      GEOS_LOG_RANK( "ghost geom = " << rk << ": " << json( geom ) );
    }
  }
  GEOS_LOG_RANK( "my final neighbors are " << json( ghost.neighbors ) );

  return ghost;


  // TODO to build the edges -> nodes map containing the all the ghost (i.e. the ghosted ones as well),
  // Let's create a vector (or matrix?) full of ones where we have edges and multiply using the adjacency matrix.
}

/**
 * @brief
 * @tparam GI A global index type
 * @param m
 * @return
 */
template< class GI >
std::tuple< std::vector< MpiRank >, std::vector< GI > >
buildGhostRankAndL2G( std::map< GI, MpiRank > const & m )
{
  std::size_t const size = std::size( m );

  std::vector< MpiRank > ghostRank;
  ghostRank.reserve( size );
  std::vector< GI > l2g;
  l2g.reserve( size );
  for( auto const & [t, rank]: m )
  {
    ghostRank.emplace_back( rank );
    l2g.emplace_back( t );
  }

  return std::make_tuple( ghostRank, l2g );
}

void buildPods( Ghost const & ghost )
{
  std::size_t const numNodes = std::size( ghost.nodes );
  std::size_t const numEdges = std::size( ghost.edges );
  std::size_t const numFaces = std::size( ghost.faces );

  auto [ghostRank, l2g] = buildGhostRankAndL2G( ghost.edges );

  NodeMgrImpl const nodeMgr( NodeLocIdx{ intConv< localIndex >( numNodes ) } );
  EdgeMgrImpl const edgeMgr( EdgeLocIdx{ intConv< localIndex >( numEdges ) }, std::move( ghostRank ), std::move( l2g ) );
  FaceMgrImpl const faceMgr( FaceLocIdx{ intConv< localIndex >( numFaces ) } );
}

std::unique_ptr< generators::MeshMappings > doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
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
  if( curRank == 0_mpi )
  {
    // The `MPI_Scan` process will not call the reduction operator for rank 0.
    // So we need to reduce ourselves for ourselves.
    offsets.edges = updateBucketOffsets< EdgeGlbIdx >( sizes.edges, { { { 0_mpi, }, 0_egi } }, curRank );
    offsets.faces = updateBucketOffsets< FaceGlbIdx >( sizes.faces, { { { 0_mpi, }, 0_fgi } }, curRank );
    // Still we need to send this reduction to the following rank, by copying to it to the send buffer.
    std::vector< std::uint8_t > const bytes = serialize( offsets );
    std::memcpy( sendBuffer.data(), bytes.data(), bytes.size() );
  }
  else
  {
    // For the other ranks, the reduction operator will be called during the `Mpi_Scan` process.
    // So unlike for rank 0, we do not have to do it ourselves.
    // In order to provide the `sizes` to the reduction operator,
    // since `sizes` will only be used on the current rank,
    // we'll provide the information as a pointer to the instance.
    // The reduction operator will then compute the new offsets and send them to the following rank.
    std::uintptr_t const addr = reinterpret_cast<std::uintptr_t>(&sizes);
    std::memcpy( sendBuffer.data(), &addr, sizeof( std::uintptr_t ) );
  }

  MPI_Op op;
  MPI_Op_create( f, false, &op );

  MPI_Scan( sendBuffer.data(), recvBuffer.data(), maxBufferSize, MPI_BYTE, op, MPI_COMM_WORLD );

  if( curRank != 0_mpi )
  {
    offsets = deserialize( recvBuffer );
  }

  std::cout << "offsets on rank " << curRank << " -> " << json( offsets ) << std::endl;

  MpiRank const nextRank = curRank + 1_mpi;
  MaxGlbIdcs const matrixOffsets = gatherOffset( mesh, offsets.edges.at( { nextRank } ) - 1_egi, offsets.faces.at( { nextRank } ) - 1_fgi );
  std::cout << "matrixOffsets on rank " << curRank << " -> " << json( matrixOffsets ) << std::endl;

  MeshGraph const graph = buildMeshGraph( mesh, buckets, offsets, curRank );  // TODO change into buildOwnedMeshGraph?
//  if( curRank == MpiRank{ 1 } )
//  {
//    std::cout << "My graph is " << json( graph ) << std::endl;
//  }

  Ghost const ghost = assembleAdjacencyMatrix( graph, matrixOffsets, curRank );

  buildPods( ghost );

  return {};
}

std::unique_ptr< generators::MeshMappings > doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                                                              std::set< int > const & neighbors )
{
  std::set< MpiRank > neighbors_;
  for( int const & rank: neighbors )
  {
    neighbors_.insert( MpiRank{ rank } );
  }
  GEOS_LOG_RANK( "my initial neighbors are " << json( neighbors_ ) );

  return doTheNewGhosting( mesh, neighbors_ );
}

}  // end of namespace geos::ghosting
