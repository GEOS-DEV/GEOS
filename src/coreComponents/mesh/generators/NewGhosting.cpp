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


Exchange buildSpecificData( vtkSmartPointer< vtkDataSet > mesh, std::set< vtkIdType > const & cellIds )
{
  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );

  std::vector< Edge > tmpEdges;
  std::vector< Face > tmpFaces;

  // Pre-allocation of the temporary vectors.
  {
    std::size_t numEdges = 0;
    std::size_t numFaces = 0;
    for( vtkIdType const & c : cellIds )
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
  array1d< globalIndex > result;
  result.reserve( 2 * exchange.edges.size() );
  for( Edge const & edge: exchange.edges )
  {
    result.emplace_back( edge.first.get() );
    result.emplace_back( edge.second.get() );
  }

  return result;
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

Exchange convertExchange( array1d< globalIndex > const & input )
{
  Exchange exchange;

  for( auto i = 0; i < input.size(); ++ ++i )
  {
    NodeGlbIdx gn0{ input[i] }, gn1{ input[i + 1] };
    exchange.edges.insert( std::minmax( gn0, gn1 ) );
  }

  return exchange;
}

/**
 * @brief
 * @param exchanged
 * @param neighbors excluding current rank
 */
std::map< std::set< MpiRank >, std::set< Edge > > findOverlappingEdges( std::map< MpiRank, Exchange > const & exchanged,
                                                                        MpiRank curRank,
                                                                        std::set< MpiRank > const & neighbors )
{
  GEOS_LOG_RANK( "Starting findOverlappingEdges" );
  GEOS_LOG_RANK( "exchanged.size() = " << exchanged.size() );

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

  std::map< std::set< MpiRank >, std::set< Edge > > ranksToIntersections;
  for( auto const & [edge, ranks]: counts )
  {
    if( ranks.find( curRank ) != ranks.cend() )
    {
      ranksToIntersections[ranks].insert( edge );
    }
  }

  // Checking if neighbors is too wide...  // TODO do we care?
  std::set< MpiRank > usefulNeighbors;
  for( auto const & [ranks, edges]: ranksToIntersections )
  {
    if( not edges.empty() )
    {
      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
    }
  }
  std::vector< MpiRank > uselessNeighbors;
  std::set_difference( neighbors.cbegin(), neighbors.cend(), usefulNeighbors.cbegin(), usefulNeighbors.cend(), std::back_inserter( uselessNeighbors ) );
  // TODO... Remove the neighbors?
  GEOS_LOG_RANK( "Ending findOverlappingEdges" );

  return ranksToIntersections;
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

using ScannedOffsets = std::map< std::set< MpiRank >, integer >;


void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< MpiRank > const & neighbors )
{
  // Now we exchange the data with our neighbors.
  MpiRank const curRank{ MpiWrapper::commRank() };

  std::vector< NeighborCommunicator > ncs;
  for( MpiRank const & rank: neighbors )
  {
    ncs.emplace_back( rank.get() );
  }

  CommID const commId = CommunicationTools::getInstance().getCommID();
  std::map< MpiRank, Exchange > exchanged = exchange( int( commId ), ncs, buildExchangeData( mesh ) );
  exchanged[MpiRank{ MpiWrapper::commRank() }] = buildFullData( mesh );

  std::map< std::set< MpiRank >, std::set< Edge > > const overlappingEdges = findOverlappingEdges( exchanged, curRank, neighbors );
}

void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< int > const & neighbors )
{
  std::set< MpiRank > neighbors_;
  for( int const & rank: neighbors )
  {
    neighbors_.insert( MpiRank{ rank } );
  }

  return doTheNewGhosting(mesh, neighbors_);
}

}  // end of namespace geos::ghosting
