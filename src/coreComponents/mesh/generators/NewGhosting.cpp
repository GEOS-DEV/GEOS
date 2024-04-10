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

#include <vtkPointData.h>

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
using Face = std::vector< NodeGlbIdx >;

struct Exchange
{
  std::set< Edge > edges;
//  std::set< Face > faces;
};

Exchange buildExchangeData( vtkSmartPointer< vtkDataSet > mesh )
{
  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );

  std::vector< Edge > tmp;

  // Pre-allocation of the temporary.
  {
    std::size_t numEdges = 0;
    for( auto c = 0; c < mesh->GetNumberOfCells(); ++c )
    {
      vtkCell * cell = mesh->GetCell( c );
      numEdges += cell->GetNumberOfEdges();
    }
    tmp.reserve( numEdges );
  }

  // Filling the temporary.
  for( auto c = 0; c < mesh->GetNumberOfCells(); ++c )
  {
    vtkCell * cell = mesh->GetCell( c );
    for( auto e = 0; e < cell->GetNumberOfEdges(); ++e )
    {
      vtkCell * edge = cell->GetEdge( e );
      vtkIdType const ln0 = edge->GetPointId( 0 ), ln1 = edge->GetPointId( 1 );
      vtkIdType const gn0 = globalPtIds->GetValue( ln0 ), gn1 = globalPtIds->GetValue( ln1 );
      tmp.emplace_back( std::minmax( { gn0, gn1 } ) );
    }
  }

  // Removing the duplicates by copying into a `std::set`.
  std::set< Edge > const edges{ tmp.cbegin(), tmp.cend() };

  return { edges };
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
 * @param neighborhood including current rank
 */
std::map< std::set< MpiRank >, std::set< Edge > > findOverlappingEdges( std::map< MpiRank, Exchange > const & exchanged,
                                                                        MpiRank curRank,
                                                                        std::set< MpiRank > const & neighborhood )
{
  GEOS_LOG_RANK( "Starting findOverlappingEdges" );
  GEOS_LOG_RANK( "exchanged.size() = " << exchanged.size() );

  std::map< Edge, std::set< MpiRank > > counts;  // TODO Use better intersection algorithms?
  for( MpiRank const & rank: neighborhood )
  {
    for( Edge const & edge: exchanged.at( rank ).edges )
    {
      counts[edge].insert( rank );
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

  // Checking if neighborhood is too wide...  // TODO do we care?
  std::set< MpiRank > usefulNeighbors;
  for( auto const & [ranks, edges]: ranksToIntersections )
  {
    if( not edges.empty() )
    {
      usefulNeighbors.insert( ranks.cbegin(), ranks.cend() );
    }
  }
  std::vector< MpiRank > uselessNeighbors;
  std::set_difference( neighborhood.cbegin(), neighborhood.cend(), usefulNeighbors.cbegin(), usefulNeighbors.cend(), std::back_inserter( uselessNeighbors ) );
  // TODO... Remove the neighbors?
  GEOS_LOG_RANK( "Ending findOverlappingEdges" );

  return ranksToIntersections;
}

// TODO Duplicated
std::map< MpiRank, Exchange > exchange( int commId,
                                        std::vector< NeighborCommunicator > & neighbors,
                                        Exchange && data )
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

  array1d< array1d< globalIndex > > tmpOutput( neighbors.size() );

  for( integer i = 0; i < numNeighbors; ++i )
  {
    neighbors[i].mpiISendReceiveData( cv,
                                      commData.mpiSendBufferRequest( i ),
                                      tmpOutput[i],
                                      commData.mpiRecvBufferRequest( i ),
                                      commId,
                                      MPI_COMM_GEOSX );
  }
  MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );
  MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus() );

  std::map< MpiRank, Exchange > output;
  for( auto i = 0; i < numNeighbors; ++i )
  {
    output[MpiRank{ neighbors[i].neighborRank() }] = convertExchange( tmpOutput[i] );
  }
  output[MpiRank{ MpiWrapper::commRank() }] = std::move( data );
  return output;
}


void doTheNewGhosting( vtkSmartPointer< vtkDataSet > mesh,
                       std::set< int > const & neighbors )
{
  // Now we exchange the data with our neighbors.
  MpiRank const curRank{ MpiWrapper::commRank() };

  std::vector< NeighborCommunicator > ncs;
  for( int const & rank: neighbors )
  {
    ncs.emplace_back( rank );
  }

  std::set< MpiRank > neighborhood;
  for( int const & rank: neighbors )
  {
    neighborhood.insert( MpiRank{ rank } );
  }
  neighborhood.insert( curRank );

  CommID const commId = CommunicationTools::getInstance().getCommID();
  std::map< MpiRank, Exchange > const exchanged = exchange( int( commId ), ncs, buildExchangeData( mesh ) );

  std::map< std::set< MpiRank >, std::set< Edge > > const overlappingEdges = findOverlappingEdges( exchanged, curRank, neighborhood );
}

}  // end of namespace geos::ghosting
