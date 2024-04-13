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

#include <algorithm>
#include <utility>


namespace geos::ghosting
{

using NodeLocIdx = fluent::NamedType< localIndex, struct NodeLocIdxTag, fluent::Comparable, fluent::Printable >;
using NodeGlbIdx = fluent::NamedType< globalIndex, struct NodeGlbIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeLocIdx = fluent::NamedType< localIndex, struct EdgeLocIdxTag, fluent::Comparable, fluent::Printable >;
using EdgeGlbIdx = fluent::NamedType< globalIndex, struct EdgeGlbIdxTag, fluent::Comparable, fluent::Printable >;
using FaceLocIdx = fluent::NamedType< localIndex, struct FaceLocIdxTag, fluent::Comparable, fluent::Printable >;
using FaceGlbIdx = fluent::NamedType< globalIndex, struct FaceGlbIdxTag, fluent::Comparable, fluent::Printable >;
using CellLocIdx = fluent::NamedType< localIndex, struct CellLocIdxTag, fluent::Comparable, fluent::Printable >;
using CellGlbIdx = fluent::NamedType< globalIndex, struct CellGlbIdxTag, fluent::Comparable, fluent::Printable >;

using MpiRank = fluent::NamedType< int, struct MpiRankTag, fluent::Comparable, fluent::Printable >;

//struct Edge
//{
//  std::pair< NodeGlbIdx, NodeGlbIdx > nodes;
//};

//struct Face
//{
//  std::vector< NodeGlbIdx > nodes;
//};

using Edge = std::pair< NodeGlbIdx, NodeGlbIdx >;
using Face = std::set< NodeGlbIdx >;

struct Exchange
{
  std::set< Edge > edges;
  std::set< Face > faces;
//  SortedArray< Edge > edges;
//  SortedArray< Face > faces;
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
      tmpEdges.emplace_back( std::minmax( { gn0, gn1 } ) );
    }

    for( auto f = 0; f < cell->GetNumberOfFaces(); ++f )
    {
      vtkCell * face = cell->GetFace( f );
      vtkIdList * pids = face->GetPointIds();
      Face ff;
      for( auto i = 0; i < pids->GetNumberOfIds(); ++i )
      {
        vtkIdType const lni = face->GetPointId( i );
        vtkIdType const gni = globalPtIds->GetValue( lni );
        ff.insert( NodeGlbIdx{ gni } );
      }
      tmpFaces.emplace_back( ff );
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
    result.emplace_back( edge.first.get() );
    result.emplace_back( edge.second.get() );
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
    std::set< NodeGlbIdx > face;
    for( int j = 0; j < s; ++j, ++i )
    {
      face.insert( NodeGlbIdx{ input[i] } );
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


struct ScannedOffsets
{
  std::map< std::set< MpiRank >, EdgeGlbIdx > edgeOffsets;
  std::map< std::set< MpiRank >, FaceGlbIdx > faceOffsets;
  EdgeGlbIdx edgeRestart;
  FaceGlbIdx faceRestart;
};


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

  std::size_t const maxBufferSize = buildMaxBufferSize( buckets );
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
