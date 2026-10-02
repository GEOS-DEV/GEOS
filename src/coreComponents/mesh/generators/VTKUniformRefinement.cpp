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

/** @file VTKUniformRefinement.cpp */
#include "VTKUniformRefinement.hpp"
#include "VTKRefinementSharing.hpp"
#include "VTKRefinementTemplates.hpp"
#include "VTKUtilities.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "common/logger/Logger.hpp"

#include <vtkArrayDispatch.h>
#include <vtkBitArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataArrayAccessor.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <tuple>
#include <unordered_set>

namespace geos::vtk
{
namespace
{
using namespace refinement;
constexpr char const * rootCellName = "_geosUniformRootCellId";
constexpr char const * parentCellName = "_geosUniformParentCellId";
constexpr char const * generationName = "_geosUniformGeneration";
constexpr char const * childName = "_geosUniformChildOrdinal";
constexpr char const * ownerName = "_geosUniformRootOwner";
constexpr char const * sourceTypeName = "_geosUniformSourceType";
constexpr char const * sourceAttributeName = "_geosUniformSourceAttribute";

// Keep cell context on the error path; successful refinement creates no
// diagnostic strings or type-erased callback allocations per parent.
template< typename WORK >
void withCellContext( std::string const & block, int level, vtkIdType parentId, int vtkType, WORK && work )
{
  try
  {
    work();
  }
  catch( std::exception const & error )
  {
    throw std::runtime_error( GEOS_FMT( "block '{}', level {}, parent global ID {}, VTK cell type {}: {}",
                                        block, level, parentId, vtkType, error.what() ) );
  }
}

std::uint64_t add( std::uint64_t a, std::uint64_t b )
{
  if( b > UINT64_MAX - a )
  {
    throw std::overflow_error( "Uniform refinement count addition overflow" );
  }
  return a + b;
}
std::uint64_t multiply( std::uint64_t a, std::uint64_t b )
{
  if( a && b > UINT64_MAX / a )
  {
    throw std::overflow_error( "Uniform refinement count multiplication overflow" );
  }
  return a * b;
}
void localCount( std::uint64_t count, char const * description = "local count" )
{
  auto const limit = std::min( static_cast< std::uint64_t >( std::numeric_limits< localIndex >::max() ),
                               static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) );
  if( count > limit || count > std::numeric_limits< std::size_t >::max() / sizeof( vtkIdType ) )
  {
    throw std::overflow_error( std::string( "Uniform refinement " ) + description + " exceeds GEOS/VTK storage" );
  }
}

struct ReadIds
{
  Connectivity & ids;
  template< typename Array > void operator()( Array * array ) const
  {
    vtkDataArrayAccessor< Array > access( array );
    using Value = typename decltype( access )::APIType;
    auto const limit = std::min( static_cast< std::uintmax_t >( std::numeric_limits< vtkIdType >::max() ),
                                 static_cast< std::uintmax_t >( std::numeric_limits< globalIndex >::max() ) );
    for( vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i )
    {
      Value const value = access.Get( i, 0 );
      if constexpr ( std::is_signed_v< Value > )
      {
        if( value < 0 )
        {
          throw std::invalid_argument( "Negative uniform refinement active ID" );
        }
      }
      if( static_cast< std::uintmax_t >( value ) > limit )
      {
        throw std::overflow_error( "Uniform refinement active ID exceeds storage" );
      }
      ids.push_back( static_cast< vtkIdType >( value ) );
    }
  }
};
Connectivity exactIds( vtkDataArray * array, vtkIdType count )
{
  if( !array && count == 0 )
  {
    return {};
  }
  if( array == nullptr || array->GetNumberOfComponents() != 1 || array->GetNumberOfTuples() != count )
  {
    throw std::invalid_argument( "Uniform refinement requires scalar active IDs with matching tuple counts" );
  }
  Connectivity ids;
  ids.reserve( count );
  if( !vtkArrayDispatch::DispatchByValueType< vtkArrayDispatch::Integrals >::Execute( array, ReadIds{ ids } ) )
  {
    throw std::invalid_argument( "Uniform refinement active IDs require integral storage" );
  }
  return ids;
}
struct ReadAttributes
{
  std::vector< int > & values;
  template< typename Array > void operator()( Array * array ) const
  {
    vtkDataArrayAccessor< Array > access( array );
    using Value = typename decltype( access )::APIType;
    for( vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i )
    {
      Value const value = access.Get( i, 0 );
      if constexpr ( std::is_integral_v< Value > )
      {
        if constexpr ( std::is_signed_v< Value > )
        {
          if( value < -1 )
          {
            throw std::invalid_argument( "Invalid uniform refinement region attribute" );
          }
        }
        if( value >= 0 && static_cast< std::uintmax_t >( value ) > static_cast< std::uintmax_t >( INT_MAX ) )
        {
          throw std::overflow_error( "Uniform refinement region attribute exceeds int" );
        }
      }
      else
      {
        if( !std::isfinite( value ) || value < -1 || static_cast< double >( value ) > static_cast< double >( INT_MAX ) )
        {
          throw std::invalid_argument( "Invalid uniform refinement region attribute" );
        }
        if( std::abs( static_cast< double >( value ) - static_cast< double >( static_cast< int >( value ) ) ) > 0 )
        {
          throw std::invalid_argument( "Uniform refinement region attribute must be integral" );
        }
      }
      values.push_back( static_cast< int >( value ) );
    }
  }
};

ElementType volumeType( Cell const & cell )
{
  if( cell.prismSides )
  {
    return static_cast< ElementType >( static_cast< int >( ElementType::Prism5 ) + cell.prismSides - 5 );
  }
  switch( cell.vtkType )
  {
    case VTK_TETRA:
      return ElementType::Tetrahedron;
    case VTK_PYRAMID:
      return ElementType::Pyramid;
    case VTK_WEDGE:
      return ElementType::Wedge;
    case VTK_HEXAHEDRON:
      return ElementType::Hexahedron;
    default:
      throw std::invalid_argument( "Unsupported normalized uniform refinement volume type" );
  }
}
bool isSurface( Cell const & cell ) { return cell.vtkType == VTK_TRIANGLE || cell.vtkType == VTK_QUAD || cell.vtkType == VTK_POLYGON; }
void countCell( CellCounts & counts, Cell const & cell )
{
  if( cell.prismSides )
  {
    ++counts.prisms[cell.prismSides];
  }
  else
  {
    switch( cell.vtkType )
    {
      case VTK_TETRA:
        ++counts.tetrahedra;
        break;
      case VTK_PYRAMID:
        ++counts.pyramids;
        break;
      case VTK_WEDGE:
        ++counts.wedges;
        break;
      case VTK_HEXAHEDRON:
        ++counts.hexahedra;
        break;
      default:
        throw std::invalid_argument( "Unsupported uniform refinement count type" );
    }
  }
}

void validateVolumeStorage( CellCounts const & counts )
{
  // GEOS arrays use localIndex for flattened offsets as well as object IDs.
  // Counting cells alone cannot bound connectivity allocation. Face-node
  // incidence conservatively bounds unique faces and mandatory cell/edge/face
  // arrays before fine allocation. Shared faces can make exact storage smaller.
  std::uint64_t incidence = add( multiply( counts.hexahedra, 24 ), multiply( counts.tetrahedra, 12 ) );
  incidence = add( incidence, add( multiply( counts.wedges, 18 ), multiply( counts.pyramids, 16 ) ) );
  for( int n = 5; n <= 11; ++n )
  {
    incidence = add( incidence, multiply( counts.prisms[n], 6 * n ) );
  }
  localCount( incidence, "conservative connectivity bound" );
}

struct Source
{
  ElementType type;
  int attribute;
};
struct State
{
  std::uint64_t ns;
  bool hasSurfaces = false;
  std::string name;
  std::vector< Coordinates > coordinates;
  Connectivity pointIds, cellIds;
  std::vector< Cell > cells;
  std::vector< bool > reverseSurface;
  std::vector< Source > sources;
  Connectivity roots, parents, generations, ordinals, rootOwners;
  vtkSmartPointer< vtkPointData > pointData;
  vtkSmartPointer< vtkCellData > cellData;
  std::vector< EntitySupport > interfaces;
  std::vector< Connectivity > buckets;
  std::vector< std::vector< SurfaceSide > > sides;
};

struct FieldStorage
{
  std::uint64_t tupleBytes{}, largestComponents{};
};

FieldStorage fieldStorage( vtkDataSetAttributes & data )
{
  FieldStorage result;
  for( int a = 0; a < data.GetNumberOfArrays(); ++a )
  {
    auto * array = data.GetAbstractArray( a );
    auto const components = static_cast< std::uint64_t >( array->GetNumberOfComponents() );
    result.largestComponents = std::max( result.largestComponents, components );
    if( auto * strings = vtkStringArray::SafeDownCast( array ) )
    {
      for( std::uint64_t c = 0; c < components; ++c )
      {
        std::uint64_t longest = 0;
        for( vtkIdType t = 0; t < strings->GetNumberOfTuples(); ++t )
          longest = std::max( longest, static_cast< std::uint64_t >( strings->GetValue( t * components + c ).size() ) );
        result.tupleBytes = add( result.tupleBytes, add( sizeof( vtkStdString ), add( longest, 1 ) ) );
      }
    }
    else
    {
      // One byte per bit is a conservative bound for packed vtkBitArray storage.
      auto const width = vtkBitArray::SafeDownCast( array ) ? 1 : array->GetDataTypeSize();
      if( width <= 0 )
        throw std::invalid_argument( "Unsupported refinement field storage forecast" );
      result.tupleBytes = add( result.tupleBytes, multiply( components, width ) );
    }
  }
  return result;
}

std::uint64_t volumeConnectivity( CellCounts const & counts )
{
  auto value = add( multiply( counts.hexahedra, 8 ), multiply( counts.tetrahedra, 4 ) );
  value = add( value, add( multiply( counts.wedges, 6 ), multiply( counts.pyramids, 5 ) ) );
  for( int n = 5; n <= 11; ++n )
    value = add( value, multiply( counts.prisms[n], 2 * n ) );
  return value;
}

std::uint64_t newVolumePointsUpperBound( CellCounts const & counts )
{
  auto value = add( multiply( counts.hexahedra, 19 ), multiply( counts.tetrahedra, 6 ) );
  value = add( value, add( multiply( counts.wedges, 12 ), multiply( counts.pyramids, 9 ) ) );
  for( int n = 5; n <= 11; ++n )
    value = add( value, multiply( counts.prisms[n], 4 * n + 3 ) );
  return value;
}

std::vector< RefinementResourceEstimate > resourceForecast( std::vector< State > const & states, int levels,
                                                            std::uint64_t retainedInputBytes, std::size_t neighbors,
                                                            std::uint64_t globalBucketWidth )
{
  std::vector< RefinementResourceEstimate > result( levels );
  for( auto const & state : states )
  {
    CellCounts volumes;
    std::uint64_t triangles = 0, quads = 0, polygonChildren = 0, polygonPoints = 0, coarseConnectivity = 0;
    for( auto const & cell : state.cells )
    {
      coarseConnectivity = add( coarseConnectivity, cell.points.size() );
      if( cell.vtkType == VTK_TRIANGLE )
        ++triangles;
      else if( cell.vtkType == VTK_QUAD )
        ++quads;
      else if( isSurface( cell ) )
      {
        polygonChildren = add( polygonChildren, cell.points.size() );
        polygonPoints = add( polygonPoints, add( cell.points.size(), 1 ) );
      }
      else
        countCell( volumes, cell );
    }
    // Final VTK arrays pad buckets to a communicator-wide capacity. A global
    // coarse maximum also bounds every fine bucket's actual volume-side IDs.
    std::uint64_t const bucketWidth = state.ns == 0 ? 0 : globalBucketWidth;
    FieldStorage const pointFields = fieldStorage( *state.pointData ), cellFields = fieldStorage( *state.cellData );
    auto points = static_cast< std::uint64_t >( state.coordinates.size() );
    double parentBytes = static_cast< double >( coarseConnectivity ) * sizeof( vtkIdType ) +
                         static_cast< double >( state.cells.size() ) * ( sizeof( Cell ) + 8 * sizeof( vtkIdType ) + cellFields.tupleBytes ) +
                         static_cast< double >( points ) * ( sizeof( Coordinates ) + sizeof( vtkIdType ) + pointFields.tupleBytes );
    for( int l = 0; l < levels; ++l )
    {
      points = add( points, newVolumePointsUpperBound( volumes ) );
      points = add( points, add( multiply( triangles, 3 ), add( multiply( quads, 5 ), polygonPoints ) ) );
      volumes = volumes.next();
      triangles = multiply( triangles, 4 );
      quads = add( multiply( quads, 4 ), polygonChildren );
      polygonChildren = polygonPoints = 0;
      auto const surfaces = add( triangles, quads );
      auto const cells = add( volumes.total(), surfaces );
      auto const connectivity = add( volumeConnectivity( volumes ), add( multiply( triangles, 3 ), multiply( quads, 4 ) ) );
      validateVolumeStorage( volumes );
      localCount( cells );
      localCount( multiply( points, 3 ), "conservative point-coordinate bound" );
      localCount( multiply( points, bucketWidth ), "conservative collocation bound" );
      for( auto const & [tuples, fields] : { std::pair{ points, pointFields }, std::pair{ cells, cellFields } } )
      {
        if( multiply( tuples, fields.largestComponents ) > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) )
          throw std::overflow_error( "Uniform refinement field tuple/component forecast exceeds vtkIdType" );
      }
      auto const fieldBytes = add( multiply( points, pointFields.tupleBytes ), multiply( cells, cellFields.tupleBytes ) );
      auto const collocationBytes = multiply( multiply( points, bucketWidth ), sizeof( vtkIdType ) );
      auto vtkBytes = add( fieldBytes, collocationBytes );
      vtkBytes = add( vtkBytes, multiply( points, sizeof( Coordinates ) + sizeof( vtkIdType ) ) );
      vtkBytes = add( vtkBytes, multiply( connectivity, sizeof( vtkIdType ) ) );
      vtkBytes = add( vtkBytes, add( multiply( cells, 9 * sizeof( vtkIdType ) + sizeof( unsigned char ) ), sizeof( vtkIdType ) ) );
      if( vtkBytes > std::numeric_limits< std::size_t >::max() )
        throw std::overflow_error( "Uniform refinement data forecast exceeds addressable storage" );
      // Mandatory owned-volume face/edge/node incidence; ghosts are a separate
      // worst case in which every coarse neighbor contributes its whole dataset.
      auto const incidence = add( multiply( volumes.hexahedra, 24 ), add( multiply( volumes.tetrahedra, 12 ),
                                                                          add( multiply( volumes.wedges, 18 ), multiply( volumes.pyramids, 16 ) ) ) );
      auto const geosBytes = add( multiply( incidence, 4 * sizeof( localIndex ) + 2 * sizeof( globalIndex ) ),
                                  multiply( points, 3 * sizeof( real64 ) + sizeof( globalIndex ) + sizeof( localIndex ) ) );
      double const childBytes = static_cast< double >( vtkBytes ) + static_cast< double >( cells ) *
                                ( sizeof( Cell ) + sizeof( Source ) + 7 * sizeof( vtkIdType ) );
      // Registry keys/supports, hash entries, child plans and descriptors. This
      // accounts for their payloads, not allocator-dependent capacities/overhead.
      double const scratch = static_cast< double >( points ) *
                             ( sizeof( PointRecipe ) + 2 * sizeof( EntityKey ) + 44 * sizeof( vtkIdType ) + 6 * sizeof( void * ) ) +
                             static_cast< double >( cells ) * ( 2 * sizeof( vtkIdType ) + sizeof( double ) + sizeof( ChildCellKey ) + sizeof( Participants ) );
      double const headerBytes = 256.0 + 8.0 * ( neighbors + 1 );
      double const exchange = 2.0 * neighbors *
                              ( static_cast< double >( points ) * ( headerBytes + pointFields.tupleBytes + bucketWidth * sizeof( vtkIdType ) ) +
                                static_cast< double >( surfaces ) * ( headerBytes + cellFields.tupleBytes ) );
      auto & estimate = result[l];
      estimate.volumeCells = add( estimate.volumeCells, volumes.total() );
      estimate.surfaceCellCopies = add( estimate.surfaceCellCopies, surfaces );
      estimate.pointCopiesUpperBound = add( estimate.pointCopiesUpperBound, points );
      estimate.connectivityEntries = add( estimate.connectivityEntries, connectivity );
      estimate.fieldBytesUpperBound = add( estimate.fieldBytesUpperBound, fieldBytes );
      estimate.vtkBytesUpperBound = add( estimate.vtkBytesUpperBound, vtkBytes );
      estimate.geosOwnedConnectivityBytes = add( estimate.geosOwnedConnectivityBytes, geosBytes );
      estimate.exchangeBytesUpperBound += exchange;
      estimate.modeledRefinerPeakBytes += parentBytes + childBytes + scratch + exchange + vtkBytes;
      parentBytes = childBytes;
    }
  }
  for( auto & estimate : result )
    estimate.modeledRefinerPeakBytes += retainedInputBytes;
  return result;
}

CommunicationStatistics communicationDelta( CommunicationStatistics const & current, CommunicationStatistics const & previous )
{
  return { current.directoryExchanges - previous.directoryExchanges,
           current.neighborExchanges - previous.neighborExchanges,
           current.payloadChunksSent - previous.payloadChunksSent,
           current.payloadBytesSent - previous.payloadBytesSent,
           current.countBytesSent - previous.countBytesSent };
}

template< std::size_t N >
std::array< bool, N > sumStatistics( std::array< std::uint64_t, N > const & local,
                                     std::array< std::uint64_t, N > & totals, MPI_Comm communicator )
{
  // Same two-limb reduction as ID allocation. It supports 64-bit localIndex
  // builds without letting MPI_SUM overflow before the result is inspected.
  constexpr std::uint64_t mask = UINT64_C( 0xffffffff );
  std::array< std::uint64_t, 2 * N > limbs{}, sums{};
  std::array< bool, N > overflow{};
  for( std::size_t i = 0; i < N; ++i )
  {
    limbs[2*i] = local[i] & mask;
    limbs[2*i+1] = local[i] >> 32;
  }
  MpiWrapper::allReduce( limbs, sums, MpiWrapper::Reduction::Sum, communicator );
  for( std::size_t i = 0; i < N; ++i )
  {
    auto const high = sums[2*i+1] + ( sums[2*i] >> 32 );
    overflow[i] = high > mask;
    totals[i] = overflow[i] ? UINT64_MAX : ( high << 32 ) | ( sums[2*i] & mask );
  }
  return overflow;
}

void reportForecast( std::vector< RefinementResourceEstimate > & estimates, bool report, Communication & comm, MPI_Comm communicator )
{
  GEOS_MARK_SCOPE_STR( "uniformRefinement/statistics" );
  // Nonempty volumes grow by at least six, so representable counts prove the
  // number of levels is smaller than the bit width; this is not a level cap.
  constexpr std::size_t width = std::numeric_limits< std::uint64_t >::digits;
  std::array< std::uint64_t, width > volumes{}, globalVolumes{};
  std::array< std::uint64_t, 7 * width > local{}, maximum{};
  std::array< double, 2 * width > model{}, maximumModel{};
  comm.checked( "resource report buffers", [&]
  {
    if( estimates.size() > width )
      throw std::logic_error( "Missing refinement growth preflight" );
    for( std::size_t l = 0; l < estimates.size(); ++l )
    {
      auto const & estimate = estimates[l];
      volumes[l] = estimate.volumeCells;
      std::size_t offset = 7 * l;
      for( auto value : { estimate.volumeCells, estimate.surfaceCellCopies, estimate.pointCopiesUpperBound,
                          estimate.connectivityEntries, estimate.fieldBytesUpperBound, estimate.vtkBytesUpperBound,
                          estimate.geosOwnedConnectivityBytes } )
        local[offset++] = value;
      model[2*l] = estimate.modeledRefinerPeakBytes;
      model[2*l+1] = estimate.exchangeBytesUpperBound;
    }
  } );
  // Global active-cell IDs must fit even when reporting is disabled.
  auto const overflow = sumStatistics( volumes, globalVolumes, communicator );
  if( report )
  {
    MpiWrapper::allReduce( local, maximum, MpiWrapper::Reduction::Max, communicator );
    MpiWrapper::allReduce( model, maximumModel, MpiWrapper::Reduction::Max, communicator );
  }
  comm.checked( "resource report", [&]
  {
    for( std::size_t l = 0; l < estimates.size(); ++l )
    {
      auto const limit = std::min( static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ),
                                   static_cast< std::uint64_t >( std::numeric_limits< globalIndex >::max() ) );
      if( overflow[l] || ( globalVolumes[l] && globalVolumes[l] - 1 > limit ) )
        throw std::overflow_error( "Uniform refinement forecast exceeds global cell ID storage" );
      if( report )
        estimates[l].geosGhostConnectivityBytesModel =
          static_cast< double >( maximum[7*l+6] ) * ( comm.size() - 1 );
      if( report && comm.rank() == 0 )
        GEOS_LOG( GEOS_FMT( "Uniform refinement forecast level {}: ownedCells={}, maxOwnedCells={}, maxPointCopiesBound={}, "
                            "maxFieldBytesBound={}, maxVtkBytesBound={}, maxRefinerPeakBytesModel={}, "
                            "maxGeosOwnedConnectivityBytesModel={}, maxExchangeBytesBound={}, maxGeosGhostConnectivityBytesModel={}",
                            l + 1, globalVolumes[l], maximum[7*l], maximum[7*l+2], maximum[7*l+4], maximum[7*l+5],
                            maximumModel[2*l], maximum[7*l+6], maximumModel[2*l+1],
                            static_cast< double >( maximum[7*l+6] ) * ( comm.size() - 1 ) ) );
    }
  } );
}

void reportLevel( RefinementLevelStatistics const & stats, int level, bool report, Communication & comm, MPI_Comm communicator )
{
  GEOS_MARK_SCOPE_STR( "uniformRefinement/statistics" );
  if( !report )
    return;
  std::array< std::uint64_t, 4 > const localCounts{ stats.ownedVolumeCells, stats.mainPointCopies, stats.ownedMainPoints, stats.sharedMainPointCopies };
  std::array< std::uint64_t, 4 > totals{};
  std::array< std::uint64_t, 6 > const localMax{ stats.ownedVolumeCells, stats.communication.directoryExchanges, stats.communication.neighborExchanges,
                                                 stats.communication.payloadChunksSent, stats.communication.payloadBytesSent, stats.communication.countBytesSent };
  std::array< std::uint64_t, 6 > maximum{};
  auto const overflow = sumStatistics( localCounts, totals, communicator );
  MpiWrapper::allReduce( localMax, maximum, MpiWrapper::Reduction::Max, communicator );
  comm.checked( "level statistics report", [&]
  {
    if( std::any_of( overflow.begin(), overflow.end(), []( bool value ) { return value; } ) )
    {
      if( comm.rank() == 0 )
        GEOS_LOG( "Uniform refinement global statistics exceed uint64_t; local statistics remain available" );
      return;
    }
    if( report && comm.rank() == 0 )
      GEOS_LOG( GEOS_FMT( "Uniform refinement level {}: ownedCells={}, maxOwnedCells={}, meanOwnedCells={}, uniqueMainPoints={}, "
                          "mainPointCopies={}, sharedMainPointCopies={}, maxRankDirectoryExchanges={}, maxRankNeighborExchanges={}, "
                          "maxRankPayloadMessages={}, maxRankPayloadBytes={}, maxRankCountBytes={}",
                          level, totals[0], maximum[0], static_cast< double >( totals[0] ) / comm.size(), totals[2], totals[1], totals[3],
                          maximum[1], maximum[2], maximum[3], maximum[4], maximum[5] ) );
  } );
}

double supportExtent( std::vector< Coordinates > const & coordinates, Connectivity const & support )
{
  Coordinates low = coordinates.at( support.front() ), high = low;
  for( vtkIdType p : support )
  {
    for( int d = 0; d < 3; ++d )
    {
      low[d] = std::min( low[d], coordinates.at( p )[d] );
      high[d] = std::max( high[d], coordinates.at( p )[d] );
    }
  }
  double extent = 0;
  for( int d = 0; d < 3; ++d )
  {
    extent = std::max( extent, high[d] - low[d] );
  }
  if( !std::isfinite( extent ) || extent <= 0 )
  {
    throw std::invalid_argument( "Uniform refinement support has no finite positive extent" );
  }
  return extent;
}

TransferPolicies policiesFor( State const & state, UniformRefinementOptions const & options )
{
  auto policies = state.ns == 0 ? options.fields : TransferPolicies{};
  if( state.ns != 0 )
  {
    auto const found = options.faceBlockFields.find( state.name );
    if( found != options.faceBlockFields.end() )
    {
      policies = found->second;
    }
  }
  policies.excludedPointArrays.insert( "vtkOriginalPointIds" );
  policies.excludedCellArrays.insert( "vtkOriginalCellIds" );
  return policies;
}

void noReservedArrays( vtkDataSetAttributes & data )
{
  for( int a = 0; a < data.GetNumberOfArrays(); ++a )
  {
    auto const * name = data.GetAbstractArray( a )->GetName();
    if( name && std::string( name ).starts_with( "_geosUniform" ) )
    {
      throw std::invalid_argument( "Input already contains uniform refinement state" );
    }
  }
}

State readState( vtkDataSet & input, std::uint64_t ns, std::string name, UniformRefinementOptions const & options, int rank )
{
  State state;
  state.ns = ns;
  state.name = std::move( name );
  localCount( multiply( input.GetNumberOfPoints(), 3 ) );
  localCount( input.GetNumberOfCells() );
  noReservedArrays( *input.GetPointData() );
  noReservedArrays( *input.GetCellData() );
  state.pointIds = exactIds( input.GetPointData()->GetGlobalIds(), input.GetNumberOfPoints() );
  state.cellIds = exactIds( input.GetCellData()->GetGlobalIds(), input.GetNumberOfCells() );
  std::set< vtkIdType > uniqueCells( state.cellIds.begin(), state.cellIds.end() );
  if( uniqueCells.size() != state.cellIds.size() )
  {
    throw std::invalid_argument( "Duplicate rank-local uniform refinement cell IDs" );
  }
  state.coordinates.resize( state.pointIds.size() );
  for( std::size_t p = 0; p < state.pointIds.size(); ++p )
  {
    input.GetPoint( p, state.coordinates[p].data() );
  }
  PointRegistry points( state.coordinates, state.pointIds, ns );
  std::vector< int > attributes;
  auto * attribute = input.GetCellData()->GetArray( options.regionAttribute.c_str() );
  if( ns == 0 && attribute )
  {
    if( attribute->GetNumberOfComponents() != 1 || attribute->GetNumberOfTuples() != input.GetNumberOfCells() )
    {
      throw std::invalid_argument( "Uniform refinement region attribute tuple mismatch" );
    }
    if( !vtkArrayDispatch::Dispatch::Execute( attribute, ReadAttributes{ attributes } ) )
    {
      throw std::invalid_argument( "Unsupported uniform refinement region storage" );
    }
  }
  else
  {
    attributes.assign( state.cellIds.size(), -1 );
  }
  auto * ghosts = input.GetCellData()->GetGhostArray();
  if( ghosts && ( ghosts->GetNumberOfComponents() != 1 || ghosts->GetNumberOfTuples() != input.GetNumberOfCells() ) )
  {
    throw std::invalid_argument( "Uniform refinement input ghost tuple mismatch" );
  }
  for( vtkIdType c = 0; c < input.GetNumberOfCells(); ++c )
  {
    if( ghosts && ghosts->GetValue( c ) )
    {
      withCellContext( state.name, 0, state.cellIds[c], input.GetCellType( c ), []
      {
        throw std::invalid_argument( "Uniform refinement input contains unnormalized cell ghosts" );
      } );
    }
    vtkCell & inputCell = *input.GetCell( c );
    withCellContext( state.name, 0, state.cellIds[c], inputCell.GetCellType(), [&]
    {
      Cell cell;
      bool reverse = false;
      if( inputCell.GetCellDimension() == 3 && ns == 0 )
      {
        cell = normalizeCoarseCell( inputCell, points );
        validateGeometry( cell, points );
        state.sources.push_back( { volumeType( cell ), attributes[c] } );
      }
      else if( inputCell.GetCellDimension() == 2 &&
               ( inputCell.GetCellType() == VTK_TRIANGLE || inputCell.GetCellType() == VTK_QUAD || inputCell.GetCellType() == VTK_POLYGON ) )
      {
        state.hasSurfaces = true;
        Connectivity original, global;
        std::map< vtkIdType, vtkIdType > local;
        for( vtkIdType p = 0; p < inputCell.GetNumberOfPoints(); ++p )
        {
          vtkIdType const id = inputCell.GetPointId( p );
          if( id < 0 || static_cast< std::size_t >( id ) >= state.pointIds.size() )
          {
            throw std::invalid_argument( "Invalid surface corner index" );
          }
          original.push_back( id );
          global.push_back( state.pointIds[id] );
          local.emplace( state.pointIds[id], id );
        }
        if( original.size() < 3 || original.size() > 11 )
        {
          throw std::invalid_argument( "Unsupported uniform surface polygon arity" );
        }
        auto const canonical = canonicalCycle( global );
        for( vtkIdType id : canonical )
        {
          cell.points.push_back( local.at( id ) );
        }
        auto const first = std::find( original.begin(), original.end(), cell.points[0] );
        auto const index = static_cast< std::size_t >( first - original.begin() );
        reverse = cell.points[1] != original[( index + 1 ) % original.size()];
        cell.vtkType = cell.points.size() == 3 ? VTK_TRIANGLE : cell.points.size() == 4 ? VTK_QUAD : VTK_POLYGON;
        state.sources.push_back( { ElementType::Polygon, attributes[c] } );
      }
      else
      {
        throw std::invalid_argument( "Unsupported lower-dimensional or auxiliary uniform refinement cell" );
      }
      state.cells.push_back( std::move( cell ) );
      state.reverseSurface.push_back( reverse );
    } );
  }
  auto const policies = policiesFor( state, options );
  vtkNew< vtkPointData > pointInput;
  pointInput->ShallowCopy( input.GetPointData() );
  if( options.requiredPointArrays )
  {
    std::set< std::string > required;
    if( ns == 0 )
      required = *options.requiredPointArrays;
    else
      for( auto const & [arrayName, policy] : policies.pointArrays )
      {
        GEOS_UNUSED_VAR( policy );
        required.insert( arrayName );
      }
    for( int a = pointInput->GetNumberOfArrays() - 1; a >= 0; --a )
    {
      auto * array = pointInput->GetAbstractArray( a );
      auto const * arrayName = array->GetName();
      if( !arrayName || !required.count( arrayName ) )
        pointInput->RemoveArray( a );
    }
  }
  state.pointData = transferPointData( *pointInput, points, policies );
  Connectivity identity( state.cells.size() );
  std::iota( identity.begin(), identity.end(), 0 );
  vtkNew< vtkCellData > cellInput;
  cellInput->ShallowCopy( input.GetCellData() );
  auto const & requiredCells = ns == 0 ? options.requiredCellArrays : options.requiredFaceBlockCellArrays;
  if( requiredCells )
  {
    for( int a = cellInput->GetNumberOfArrays() - 1; a >= 0; --a )
    {
      auto const * arrayName = cellInput->GetAbstractArray( a )->GetName();
      if( !arrayName || !requiredCells->count( arrayName ) )
        cellInput->RemoveArray( a );
    }
  }
  state.cellData =
    transferCellData( *cellInput, state.cells.size(), identity, std::vector< double >( identity.size(), 1 ), policies );
  state.roots = state.cellIds;
  state.rootOwners.assign( state.cells.size(), rank );
  if( state.hasSurfaces )
  {
    state.sides.resize( state.cells.size() );
  }
  if( ns != 0 )
  {
    auto * collocation = vtkIdTypeArray::SafeDownCast( input.GetPointData()->GetArray( "collocated_nodes" ) );
    if( input.GetNumberOfPoints() == 0 && !collocation )
    {
      return state;
    }
    if( !collocation || collocation->GetNumberOfTuples() != input.GetNumberOfPoints() || collocation->GetNumberOfComponents() < 1 )
    {
      throw std::invalid_argument( "Auxiliary refinement requires typed collocated_nodes with matching tuples" );
    }
    state.buckets.resize( state.pointIds.size() );
    for( vtkIdType p = 0; p < input.GetNumberOfPoints(); ++p )
    {
      auto & bucket = state.buckets[p];
      for( int b = 0; b < collocation->GetNumberOfComponents(); ++b )
      {
        vtkIdType const id = collocation->GetTypedComponent( p, b );
        if( id < -1 )
        {
          throw std::invalid_argument( "Invalid negative collocation padding" );
        }
        if( id >= 0 )
        {
          bucket.push_back( id );
        }
      }
      std::sort( bucket.begin(), bucket.end() );
      bucket.erase( std::unique( bucket.begin(), bucket.end() ), bucket.end() );
    }
  }
  return state;
}

std::vector< EntitySupport > surfaceSupports( State const & state, int rank )
{
  std::vector< EntitySupport > result;
  std::unordered_map< EntityKey, std::size_t, EntityKeyHash > indices;
  auto insert = [&]( EntityKind kind, Connectivity corners )
  {
    Connectivity global;
    for( vtkIdType p : corners )
    {
      global.push_back( state.pointIds.at( p ) );
    }
    auto key = entityKey( kind, std::move( global ), state.ns );
    if( indices.emplace( key, result.size() ).second )
    {
      result.push_back( { std::move( key ), std::move( corners ), { rank } } );
    }
  };
  for( auto const & cell : state.cells )
  {
    insert( EntityKind::face, cell.points );
    for( std::size_t p = 0; p < cell.points.size(); ++p )
    {
      insert( EntityKind::vertex, { cell.points[p] } );
      insert( EntityKind::edge, { cell.points[p], cell.points[( p + 1 ) % cell.points.size()] } );
    }
  }
  return result;
}

void coarseSharing( std::vector< State > & states, Communication & comm, MPI_Comm communicator )
{
  std::vector< EntityKey > keys;
  std::vector< std::vector< EntitySupport > > candidates;
  bool localUnusedPoints = false;
  comm.checked( "coupled coarse sharing candidates",
                [&]
  {
    candidates.resize( states.size() );
    auto const & main = states[0];
    std::vector< Cell > volumes;
    Connectivity ids;
    for( std::size_t c = 0; c < main.cells.size(); ++c )
    {
      if( !isSurface( main.cells[c] ) )
      {
        volumes.push_back( main.cells[c] );
        ids.push_back( main.cellIds[c] );
      }
    }
    auto boundary = coarseBoundary( volumes, ids, main.pointIds, comm.rank() );
    candidates[0] = std::move( boundary.entities );
    keys = std::move( boundary.volumeIds );
    std::vector< bool > used( main.coordinates.size(), false );
    for( auto const & cell : main.cells )
      for( vtkIdType p : cell.points )
        used[p] = true;
    localUnusedPoints = std::find( used.begin(), used.end(), false ) != used.end();
    for( std::size_t s = 1; s < states.size(); ++s )
    {
      candidates[s] = surfaceSupports( states[s], comm.rank() );
    }
    for( std::size_t s = 0; s < states.size(); ++s )
    {
      for( auto const & entity : candidates[s] )
      {
        keys.push_back( entity.key );
      }
      // ID-only cell records identify actual marker/auxiliary replicas.
      for( std::size_t c = 0; c < states[s].cells.size(); ++c )
      {
        if( isSurface( states[s].cells[c] ) )
        {
          keys.push_back( { states[s].ns, EntityKind::cell, { states[s].cellIds[c] } } );
        }
      }
    }
  } );
  if( MpiWrapper::max( localUnusedPoints ? 1 : 0, communicator ) != 0 )
  {
    // Unused copies can refer to another rank's interior vertex. Only this
    // exceptional input needs all original IDs; ordinary meshes route boundary
    // vertices alone. Preserve the existing unused-point authority contract.
    comm.checked( "unused original point directory",
                  [&]
    {
      auto const & main = states[0];
      std::vector< bool > routed( main.coordinates.size(), false );
      for( auto const & entity : candidates[0] )
        if( entity.key.kind == EntityKind::vertex )
          routed[entity.localCorners.front()] = true;
      for( std::size_t p = 0; p < routed.size(); ++p )
        if( !routed[p] )
        {
          EntityKey key{ main.ns, EntityKind::vertex, { main.pointIds[p] } };
          keys.push_back( key );
          candidates[0].push_back( { std::move( key ), { static_cast< vtkIdType >( p ) }, { comm.rank() } } );
        }
    } );
  }
  Sharing sharing = comm.discoverSharing( keys );
  comm.checked( "coupled coarse sharing installation",
                [&]
  {
    std::unordered_set< vtkIdType > volumeIds;
    for( std::size_t c = 0; c < states[0].cells.size(); ++c )
    {
      if( !isSurface( states[0].cells[c] ) )
      {
        volumeIds.insert( states[0].cellIds[c] );
      }
    }
    for( auto const & [key, ranks] : sharing )
    {
      if( key.kind == EntityKind::cell && key.meshNamespace == 0 && ranks.size() != 1 )
      {
        if( volumeIds.count( key.corners[0] ) )
        {
          throw std::invalid_argument( "Duplicate distributed volume ownership" );
        }
      }
    }
    for( std::size_t s = 0; s < states.size(); ++s )
    {
      for( auto & entity : candidates[s] )
      {
        if( sharing.at( entity.key ).size() > 1 )
        {
          entity.participants = sharing.at( entity.key );
          states[s].interfaces.push_back( std::move( entity ) );
        }
      }
    }
  } );
  // Coarse authoritative tuples are installed before geometry/field planning.
  std::vector< PointCreation > points;
  std::vector< CellCreation > cells;
  comm.checked( "coupled coarse tuples",
                [&]
  {
    for( auto const & state : states )
    {
      PointFieldLayout pointFields( *state.pointData );
      // Existing vertices have no edge/face recipe. Use their
      // incident cell extent, rather than an implicit unit mesh,
      // to compare independently read shared coordinates.
      std::vector< double > vertexExtents;
      if( !state.interfaces.empty() )
      {
        vertexExtents.resize( state.coordinates.size(), 0 );
        for( auto const & cell : state.cells )
        {
          double const extent = supportExtent( state.coordinates, cell.points );
          for( vtkIdType p : cell.points )
          {
            vertexExtents[p] = std::max( vertexExtents[p], extent );
          }
        }
      }
      for( auto const & support : state.interfaces )
      {
        if( support.key.kind == EntityKind::vertex )
        {
          points.push_back( { support.key, support.participants, state.coordinates[support.localCorners.front()],
                              pointFields.pack( support.localCorners.front() ),
                              vertexExtents[support.localCorners.front()] } );
        }
      }
      CellFieldLayout cellFields( *state.cellData );
      for( std::size_t c = 0; c < state.cells.size(); ++c )
      {
        if( isSurface( state.cells[c] ) )
        {
          cells.push_back( { { state.ns, state.cellIds[c], static_cast< std::uint64_t >( state.cells[c].vtkType ), 0 },
                             sharing.at( { state.ns, EntityKind::cell, { state.cellIds[c] } } ),
                             cellFields.pack( c ) } );
        }
      }
    }
  } );
  auto const pointRecords = comm.reconcileExistingPoints( points );
  auto const cellRecords = comm.reconcileExistingCells( cells );
  comm.checked( "coupled coarse tuple installation",
                [&]
  {
    for( auto & state : states )
    {
      PointFieldLayout pointFields( *state.pointData );
      for( auto const & support : state.interfaces )
      {
        if( support.key.kind == EntityKind::vertex )
        {
          auto const & record = pointRecords.at( support.key );
          state.coordinates[support.localCorners.front()] = record.position;
          pointFields.install( support.localCorners.front(), record.fields );
        }
      }
      CellFieldLayout cellFields( *state.cellData );
      for( std::size_t c = 0; c < state.cells.size(); ++c )
      {
        if( isSurface( state.cells[c] ) )
        {
          auto const key =
            ChildCellKey{ state.ns, state.cellIds[c], static_cast< std::uint64_t >( state.cells[c].vtkType ), 0 };
          cellFields.install( c, cellRecords.at( key ).fields );
          state.rootOwners[c] = sharing.at( { state.ns, EntityKind::cell, { state.cellIds[c] } } ).front();
        }
      }
    }
  } );
  // Store marker participants in side owners only as a temporary ID-replica map.
  comm.checked( "coarse surface replica participants",
                [&]
  {
    for( auto & state : states )
    {
      for( std::size_t c = 0; c < state.cells.size(); ++c )
      {
        if( isSurface( state.cells[c] ) )
        {
          state.sides[c].push_back( { { state.ns, EntityKind::cell, { state.cellIds[c] } },
                                      {},
                                      sharing.at( { state.ns, EntityKind::cell, { state.cellIds[c] } } ) } );
        }
      }
    }
  } );
}

double surfaceArea( Connectivity const & face, PointRegistry const & points )
{
  Coordinates center{};
  auto const & anchor = points.position( face.front() );
  for( vtkIdType p : face )
  {
    for( int d = 0; d < 3; ++d )
    {
      center[d] += ( points.position( p )[d] - anchor[d] ) / face.size();
    }
  }
  double area = 0;
  for( std::size_t c = 0; c < face.size(); ++c )
  {
    double a[3]{}, b[3]{}, cross[3]{};
    for( int d = 0; d < 3; ++d )
    {
      a[d] = points.position( face[c] )[d] - anchor[d] - center[d];
      b[d] = points.position( face[( c + 1 ) % face.size()] )[d] - anchor[d] - center[d];
    }
    LvArray::tensorOps::crossProduct( cross, a, b );
    area += std::sqrt( LvArray::tensorOps::l2NormSquared< 3 >( cross ) ) / 2;
  }
  if( !std::isfinite( area ) || area <= 0 )
  {
    throw std::invalid_argument( "Degenerate uniform refinement surface" );
  }
  return area;
}

struct Level
{
  State next;
  std::unique_ptr< PointRegistry > points;
  std::unique_ptr< InterfaceSharing > interfaces;
  Connectivity parentIndices;
  std::vector< double > fractions;
  std::vector< ChildCellKey > cellKeys;
  std::vector< Participants > cellParticipants;
  NodeSetMembers nodeSets;
  /// Parents that cannot be refined, and the first few reasons.
  std::uint64_t invalidCells = 0;
  std::vector< std::string > invalidReasons;
};

// Domain-boundary faces of the volume cells, in local point indices: faces with
// one local incident volume that are not shared with another rank.
std::vector< Connectivity > boundaryFaces( State const & state )
{
  std::unordered_set< EntityKey, EntityKeyHash > shared;
  for( auto const & entity : state.interfaces )
  {
    if( entity.key.kind == EntityKind::face )
    {
      shared.insert( entity.key );
    }
  }
  std::unordered_map< EntityKey, std::pair< Connectivity, int >, EntityKeyHash > faces;
  for( auto const & cell : state.cells )
  {
    if( isSurface( cell ) )
    {
      continue;
    }
    for( auto & corners : cellFaces( cell ) )
    {
      Connectivity ids;
      for( vtkIdType p : corners )
      {
        ids.push_back( state.pointIds[p] );
      }
      auto & face = faces[entityKey( EntityKind::face, std::move( ids ), state.ns )];
      face.first = std::move( corners );
      ++face.second;
    }
  }
  std::vector< Connectivity > result;
  for( auto & [key, face] : faces )
  {
    if( face.second == 1 && !shared.count( key ) )
    {
      result.push_back( std::move( face.first ) );
    }
  }
  return result;
}

void planChildren( State const & state, Level & level, int generation, int rank )
{
  GEOS_MARK_SCOPE_STR( "uniformRefinement/buildChildren" );
  auto & points = *level.points;
  for( std::size_t p = 0; p < state.cells.size(); ++p )
  {
    auto const & parent = state.cells[p];
    // Record every parent that cannot be refined, so that the error reports
    // all of them instead of the first one.
    try
    {
      withCellContext( state.name, generation, state.cellIds[p], parent.vtkType, [&]
      {
        std::vector< Cell > children;
        std::vector< double > measures;
        if( isSurface( parent ) )
        {
          for( auto & face : subdivideFace( parent.points, points ) )
          {
            measures.push_back( surfaceArea( face, points ) );
            children.push_back( { face.size() == 3 ? VTK_TRIANGLE : VTK_QUAD, std::move( face ), 0 } );
          }
        }
        else
        {
          children = subdivideCell( parent, state.cellIds[p], points ).children;
          for( auto const & child : children )
          {
            measures.push_back( signedMeasure( child, points ) );
          }
        }
        double const measure = std::accumulate( measures.begin(), measures.end(), 0. );
        for( std::size_t c = 0; c < children.size(); ++c )
        {
          level.fractions.push_back( measures[c] / measure );
          level.parentIndices.push_back( p );
          if( state.hasSurfaces )
          {
            level.cellKeys.push_back( { state.ns, state.cellIds[p], static_cast< std::uint64_t >( parent.vtkType ), c } );
            level.cellParticipants.push_back( isSurface( parent ) ? Participants{ rank } : Participants{} );
          }
          level.next.roots.push_back( state.roots[p] );
          level.next.parents.push_back( state.cellIds[p] );
          level.next.generations.push_back( generation );
          level.next.ordinals.push_back( c );
          level.next.rootOwners.push_back( state.rootOwners[p] );
          level.next.sources.push_back( state.sources[p] );
          level.next.reverseSurface.push_back( state.reverseSurface[p] );
          level.next.cells.push_back( std::move( children[c] ) );
        }
      } );
    }
    catch( std::exception const & error )
    {
      ++level.invalidCells;
      if( level.invalidReasons.size() < 5 )
      {
        level.invalidReasons.push_back( error.what() );
      }
    }
  }
}

Level planLevel( State const & state, UniformRefinementOptions const & options, int generation, int rank )
{
  GEOS_MARK_SCOPE_STR( "uniformRefinement/levelPlan" );
  Level level;
  level.next.ns = state.ns;
  level.next.hasSurfaces = state.hasSurfaces;
  level.next.name = state.name;
  level.points = std::make_unique< PointRegistry >( state.coordinates, state.pointIds, state.ns );
  planChildren( state, level, generation, rank );
  if( level.invalidCells )
  {
    // The caller reports these cells collectively; the plan is not used.
    return level;
  }
  auto & points = *level.points;
  localCount( level.next.cells.size() );
  localCount( multiply( points.points().size(), 3 ) );
  level.interfaces = std::make_unique< InterfaceSharing >( points, state.interfaces, rank );
  auto const policies = policiesFor( state, options );
  {
    GEOS_MARK_SCOPE_STR( "uniformRefinement/transferFields" );
    // Point data is transferred after the participants agree on node-set membership.
    if( std::any_of( policies.pointArrays.begin(), policies.pointArrays.end(),
                     []( auto const & policy ) { return policy.second == PointTransferPolicy::nodeSet; } ) )
    {
      if( state.ns != 0 )
      {
        throw std::invalid_argument( "Refinement node sets are only supported on the main mesh" );
      }
      level.nodeSets = boundaryNodeSets( *state.pointData, points, policies, boundaryFaces( state ) );
    }
    level.next.cellData = transferCellData( *state.cellData, state.cells.size(), level.parentIndices, level.fractions, policies );
  }
  level.next.cellIds.resize( level.next.cells.size(), -1 );
  if( state.hasSurfaces )
  {
    level.next.sides.resize( level.next.cells.size() );
  }
  return level;
}

vtkSmartPointer< vtkIdTypeArray > idArray( char const * name, Connectivity const & ids )
{
  auto array = vtkSmartPointer< vtkIdTypeArray >::New();
  array->SetName( name );
  array->SetNumberOfTuples( ids.size() );
  for( std::size_t p = 0; p < ids.size(); ++p )
  {
    array->SetValue( p, ids[p] );
  }
  return array;
}

std::set< std::pair< ElementType, int > > sourceSchema( State const & main, Communication & comm, MPI_Comm MPI_PARAM( communicator ) )
{
  std::vector< int > local, counts, offsets, gathered;
  comm.checked( "coarse source schema",
                [&]
  {
    std::set< std::pair< ElementType, int > > unique;
    for( auto const & source : main.sources )
    {
      if( getElementDim( source.type ) == 3 )
      {
        unique.emplace( source.type, source.attribute );
      }
    }
    if( unique.size() > static_cast< std::size_t >( INT_MAX / 2 ) )
    {
      throw std::overflow_error( "Too many coarse source blocks" );
    }
    for( auto const & [type, attribute] : unique )
    {
      local.push_back( static_cast< int >( type ) );
      local.push_back( attribute );
    }
    counts.resize( comm.size() );
    offsets.resize( comm.size() );
  } );
  int const count = static_cast< int >( local.size() );
#ifdef GEOS_USE_MPI
  MPI_Allgather( &count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator );
#else
  counts[0] = count;
#endif
  comm.checked( "source schema allocation",
                [&]
  {
    std::uint64_t total = 0;
    for( int rank = 0; rank < comm.size(); ++rank )
    {
      if( counts[rank] < 0 || counts[rank] % 2 )
      {
        throw std::invalid_argument( "Invalid source schema count" );
      }
      offsets[rank] = static_cast< int >( total );
      total = add( total, counts[rank] );
      if( total > INT_MAX )
      {
        throw std::overflow_error( "Coarse source metadata exceeds MPI count" );
      }
    }
    gathered.resize( total );
  } );
#ifdef GEOS_USE_MPI
  int empty = 0;
  MPI_Allgatherv( local.empty() ? &empty : local.data(), count, MPI_INT, gathered.empty() ? & empty : gathered.data(), counts.data(),
                  offsets.data(), MPI_INT, communicator );
#else
  gathered = local;
#endif
  std::set< std::pair< ElementType, int > > schema;
  comm.checked( "source schema installation",
                [&]
  {
    std::set< ElementType > types;
    std::set< int > attributes;
    for( std::size_t i = 0; i < gathered.size(); i += 2 )
    {
      auto const type = static_cast< ElementType >( gathered[i] );
      if( gathered[i] < static_cast< int >( ElementType::Tetrahedron ) ||
          gathered[i] > static_cast< int >( ElementType::Prism11 ) || gathered[i + 1] < -1 )
      {
        throw std::invalid_argument( "Invalid coarse source block metadata" );
      }
      types.insert( type );
      attributes.insert( gathered[i + 1] );
    }
    // Existing buildCellMap registers the type/attribute cross product, even
    // for empty combinations. Retain that original selector namespace.
    for( auto type : types )
    {
      for( int attribute : attributes )
      {
        schema.emplace( type, attribute );
      }
    }
  } );
  return schema;
}

std::string blockName( ElementType type, int attribute )
{
  return ( attribute == -1 ? std::string{} : std::to_string( attribute ) + "_" ) + getElementTypeName( type );
}

std::vector< RefinementBlockDescriptor > blockDescriptors( State const & main, std::set< std::pair< ElementType, int > > const & schema,
                                                           int levels )
{
  GEOS_MARK_SCOPE_STR( "uniformRefinement/finalBlockDescriptors" );
  std::vector< RefinementBlockDescriptor > blocks;
  std::map< std::tuple< ElementType, int, ElementType >, std::size_t > indices;
  std::set< std::string > names;
  for( auto const & [source, attribute] : schema )
  {
    CellCounts counts;
    Cell representative;
    switch( source )
    {
      case ElementType::Tetrahedron:
        representative.vtkType = VTK_TETRA;
        break;
      case ElementType::Pyramid:
        representative.vtkType = VTK_PYRAMID;
        break;
      case ElementType::Wedge:
        representative.vtkType = VTK_WEDGE;
        break;
      case ElementType::Hexahedron:
        representative.vtkType = VTK_HEXAHEDRON;
        break;
      default:
        representative.vtkType = VTK_POLYHEDRON;
        representative.prismSides = static_cast< int >( source ) - static_cast< int >( ElementType::Prism5 ) + 5;
    }
    countCell( counts, representative );
    for( int l = 0; l < levels; ++l )
    {
      counts = counts.next();
    }
    std::vector< ElementType > types;
    if( counts.tetrahedra )
    {
      types.push_back( ElementType::Tetrahedron );
    }
    if( counts.pyramids )
    {
      types.push_back( ElementType::Pyramid );
    }
    if( counts.wedges )
    {
      types.push_back( ElementType::Wedge );
    }
    if( counts.hexahedra )
    {
      types.push_back( ElementType::Hexahedron );
    }
    for( auto type : types )
    {
      std::string const origin = blockName( source, attribute );
      std::string const name = type == source ? origin : origin + "__refined_" + getElementTypeName( type );
      if( !names.insert( name ).second )
      {
        throw std::invalid_argument( "Refined cell-block name collision" );
      }
      indices.emplace( std::make_tuple( source, attribute, type ), blocks.size() );
      blocks.push_back( { origin, name, source, type, attribute, {} } );
    }
  }
  for( std::size_t c = 0; c < main.cells.size(); ++c )
  {
    if( !isSurface( main.cells[c] ) )
    {
      blocks.at( indices.at( { main.sources[c].type, main.sources[c].attribute, volumeType( main.cells[c] ) } ) ).cells.push_back( c );
    }
  }
  return blocks;
}

std::set< std::string > auxiliaryNames( AllMeshes & meshes, Communication & comm, MPI_Comm MPI_PARAM( communicator ) )
{
  // Only the small block schema is replicated, never geometry/connectivity.
  std::vector< char > local, gathered;
  std::vector< int > counts, offsets;
  comm.checked( "auxiliary block names",
                [&]
  {
    for( auto const & [name, mesh] : meshes.getFaceBlocks() )
    {
      GEOS_UNUSED_VAR( mesh );
      if( name.empty() || name.find( '\0' ) != std::string::npos )
      {
        throw std::invalid_argument( "Invalid auxiliary block name" );
      }
      if( add( local.size(), add( name.size(), 1 ) ) > INT_MAX )
      {
        throw std::overflow_error( "Auxiliary name metadata exceeds MPI count" );
      }
      local.insert( local.end(), name.begin(), name.end() );
      local.push_back( '\0' );
    }
    counts.resize( comm.size() );
    offsets.resize( comm.size() );
  } );
  int const count = static_cast< int >( local.size() );
#ifdef GEOS_USE_MPI
  MPI_Allgather( &count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator );
#else
  counts[0] = count;
#endif
  comm.checked( "auxiliary name allocation",
                [&]
  {
    std::uint64_t total = 0;
    for( int rank = 0; rank < comm.size(); ++rank )
    {
      if( counts[rank] < 0 )
      {
        throw std::invalid_argument( "Invalid auxiliary name count" );
      }
      offsets[rank] = static_cast< int >( total );
      total = add( total, counts[rank] );
      if( total > INT_MAX )
      {
        throw std::overflow_error( "Auxiliary name metadata exceeds MPI count" );
      }
    }
    gathered.resize( total );
  } );
#ifdef GEOS_USE_MPI
  char empty = 0;
  MPI_Allgatherv( local.empty() ? &empty : local.data(), count, MPI_CHAR, gathered.empty() ? & empty : gathered.data(), counts.data(),
                  offsets.data(), MPI_CHAR, communicator );
#else
  gathered = local;
#endif
  std::set< std::string > names;
  comm.checked( "auxiliary name installation",
                [&]
  {
    std::string name;
    for( char c : gathered )
    {
      if( c )
      {
        name.push_back( c );
      }
      else
      {
        if( name.empty() )
        {
          throw std::invalid_argument( "Empty auxiliary block schema name" );
        }
        names.insert( std::move( name ) );
        name.clear();
      }
    }
    if( !name.empty() )
    {
      throw std::invalid_argument( "Truncated auxiliary block schema" );
    }
  } );
  return names;
}

vtkSmartPointer< vtkUnstructuredGrid > buildGrid( State const & state, int bucketCapacity )
{
  GEOS_MARK_SCOPE_STR( "uniformRefinement/finalVtk" );
  auto grid = vtkSmartPointer< vtkUnstructuredGrid >::New();
  auto points = vtkSmartPointer< vtkPoints >::New();
  points->SetDataTypeToDouble();
  points->SetNumberOfPoints( state.coordinates.size() );
  for( std::size_t p = 0; p < state.coordinates.size(); ++p )
  {
    points->SetPoint( p, state.coordinates[p].data() );
  }
  grid->SetPoints( points );
  grid->Allocate( state.cells.size() );
  for( std::size_t c = 0; c < state.cells.size(); ++c )
  {
    auto corners = state.cells[c].points;
    if( state.reverseSurface[c] )
    {
      std::reverse( corners.begin() + 1, corners.end() );
    }
    grid->InsertNextCell( state.cells[c].vtkType, corners.size(), corners.data() );
  }
  grid->GetPointData()->ShallowCopy( state.pointData );
  grid->GetCellData()->ShallowCopy( state.cellData );
  grid->GetPointData()->SetGlobalIds( idArray( "GlobalPointIds", state.pointIds ) );
  grid->GetCellData()->SetGlobalIds( idArray( "GlobalCellIds", state.cellIds ) );
  grid->GetCellData()->AddArray( idArray( rootCellName, state.roots ) );
  grid->GetCellData()->AddArray( idArray( parentCellName, state.parents ) );
  grid->GetCellData()->AddArray( idArray( generationName, state.generations ) );
  grid->GetCellData()->AddArray( idArray( childName, state.ordinals ) );
  grid->GetCellData()->AddArray( idArray( ownerName, state.rootOwners ) );
  Connectivity sourceTypes, attributes;
  for( auto const & source : state.sources )
  {
    sourceTypes.push_back( static_cast< vtkIdType >( source.type ) );
    attributes.push_back( source.attribute );
  }
  grid->GetCellData()->AddArray( idArray( sourceTypeName, sourceTypes ) );
  grid->GetCellData()->AddArray( idArray( sourceAttributeName, attributes ) );
  if( state.ns != 0 )
  {
    if( bucketCapacity < 1 )
    {
      throw std::invalid_argument( "Missing global collocation capacity" );
    }
    localCount( multiply( state.buckets.size(), bucketCapacity ) );
    auto buckets = vtkSmartPointer< vtkIdTypeArray >::New();
    buckets->SetName( "collocated_nodes" );
    buckets->SetNumberOfComponents( bucketCapacity );
    buckets->SetNumberOfTuples( state.buckets.size() );
    buckets->FillValue( -1 );
    for( std::size_t p = 0; p < state.buckets.size(); ++p )
    {
      for( std::size_t b = 0; b < state.buckets[p].size(); ++b )
      {
        buckets->SetTypedComponent( p, b, state.buckets[p][b] );
      }
    }
    grid->GetPointData()->AddArray( buckets );
  }
  return grid;
}
} // namespace

void validateRefinedTransform( vtkDataSet & mesh, Coordinates const & translation, Coordinates const & scale, MPI_Comm communicator )
{
  Communication comm( communicator );
  comm.checked( "refined physical coordinate transform", [&]
  {
    int negativeAxes = 0;
    for( int d = 0; d < 3; ++d )
    {
      if( !std::isfinite( translation[d] ) || !std::isfinite( scale[d] ) || std::abs( scale[d] ) <= 0 )
      {
        throw std::invalid_argument( "Refined coordinate transform must be finite and nonsingular" );
      }
      negativeAxes += scale[d] < 0;
    }
    if( negativeAxes % 2 )
    {
      throw std::invalid_argument( "Refined coordinate transform reverses orientation" );
    }
    std::vector< Coordinates > physical( mesh.GetNumberOfPoints() );
    for( vtkIdType p = 0; p < mesh.GetNumberOfPoints(); ++p )
    {
      mesh.GetPoint( p, physical[p].data() );
      for( int d = 0; d < 3; ++d )
      {
        physical[p][d] = ( physical[p][d] + translation[d] ) * scale[d];
      }
    }
    PointRegistry points( std::move( physical ), exactIds( mesh.GetPointData()->GetGlobalIds(), mesh.GetNumberOfPoints() ) );
    for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); ++c )
    {
      vtkCell & cell = *mesh.GetCell( c );
      if( cell.GetCellDimension() == 3 )
      {
        validateGeometry( normalizeCell( cell ), points );
      }
    }
  } );
}

UniformRefinementResult refineUniformly( AllMeshes & meshes, int levels, UniformRefinementOptions const & options, MPI_Comm communicator )
{
  if( levels == 0 )
  {
    return {};
  }
  GEOS_MARK_SCOPE_STR( "uniformRefinement/total" );
  Communication comm( communicator, options.chunkBytes );
  UniformRefinementResult result;
  std::vector< State > states;
  auto const names = auxiliaryNames( meshes, comm, communicator );
  comm.checked( "coupled uniform refinement input",
                [&]
  {
    GEOS_MARK_SCOPE_STR( "uniformRefinement/preflight" );
    if( levels < 0 )
    {
      throw std::invalid_argument( "uniformRefinement must be nonnegative" );
    }
    if( !meshes.getMainMesh() )
    {
      throw std::invalid_argument( "Missing main uniform refinement mesh" );
    }
    states.push_back( readState( *meshes.getMainMesh(), 0, "main", options, comm.rank() ) );
    std::uint64_t ns = 1;
    for( auto const & name : names )
    {
      auto const found = meshes.getFaceBlocks().find( name );
      vtkSmartPointer< vtkDataSet > mesh;
      if( found == meshes.getFaceBlocks().end() )
      {
        mesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
      }
      else
      {
        mesh = found->second;
      }
      if( !mesh )
      {
        throw std::invalid_argument( "Missing auxiliary uniform refinement dataset" );
      }
      states.push_back( readState( *mesh, ns++, name, options, comm.rank() ) );
    }
  } );
  std::uint64_t localVolumes = 0;
  for( auto const & cell : states[0].cells )
    localVolumes += !isSurface( cell );
  std::uint64_t localBucketWidth = 0;
  for( auto const & state : states )
    for( auto const & bucket : state.buckets )
      localBucketWidth = std::max( localBucketWidth, static_cast< std::uint64_t >( bucket.size() ) );
  std::array< std::uint64_t, 6 > const configuration{ localVolumes, static_cast< std::uint64_t >( levels ),
                                                      static_cast< std::uint64_t >( INT_MAX - levels ), options.reportStatistics ? 1u : 0u, options.reportStatistics ? 0u : 1u, localBucketWidth };
  std::array< std::uint64_t, 6 > maximumConfiguration{};
  MpiWrapper::allReduce( configuration, maximumConfiguration, MpiWrapper::Reduction::Max, communicator );
  comm.checked( "uniform refinement growth preflight", [&]
  {
    GEOS_MARK_SCOPE_STR( "uniformRefinement/preflight" );
    if( maximumConfiguration[1] + maximumConfiguration[2] != INT_MAX ||
        ( maximumConfiguration[3] && maximumConfiguration[4] ) )
      throw std::invalid_argument( "Uniform refinement levels/statistics configuration differs across ranks" );
    if( maximumConfiguration[0] == 0 )
      throw std::invalid_argument( "Positive uniform refinement requires coarse volume cells" );
    // All supported volume templates grow by at least six. Check this common
    // lower bound first so empty ranks cannot loop through an enormous level
    // count while another rank fails its local recurrence.
    auto minimum = maximumConfiguration[0];
    for( int l = 0; l < levels; ++l )
    {
      minimum = multiply( minimum, 6 );
      localCount( minimum, "minimum volume growth" );
    }
    // Checked recurrences prove final volume and surface-cell representation
    // before subdivision. Exact point/connectivity checks follow each plan.
    for( auto const & state : states )
    {
      CellCounts counts;
      std::uint64_t triangles = 0, quads = 0;
      for( auto const & cell : state.cells )
      {
        if( isSurface( cell ) )
        {
          if( cell.points.size() == 3 )
          {
            triangles = add( triangles, 4 );
          }
          else
          {
            quads = add( quads, cell.points.size() );
          }
        }
        else
        {
          countCell( counts, cell );
        }
      }
      for( int l = 0; l < levels; ++l )
      {
        counts = counts.next();
        validateVolumeStorage( counts );
        if( l )
        {
          triangles = multiply( triangles, 4 );
          quads = multiply( quads, 4 );
        }
        localCount( add( counts.total(), add( triangles, quads ) ) );
      }
    }
  } );
  std::vector< MainFace > mainFaces;
  std::vector< CoarseSurface > surfaces;
  std::vector< std::pair< std::size_t, std::size_t > > surfaceIndices;
  comm.checked( "coupled coarse surface queries",
                [&]
  {
    for( auto const & cell : states[0].cells )
    {
      if( !isSurface( cell ) )
      {
        for( auto const & face : cellFaces( cell ) )
        {
          Connectivity ids;
          for( vtkIdType p : face )
          {
            ids.push_back( states[0].pointIds[p] );
          }
          mainFaces.push_back( { std::move( ids ), { comm.rank() } } );
        }
      }
    }
  } );
  comm.validateVolumeFaces( mainFaces );
  coarseSharing( states, comm, communicator );
  auto const sources = sourceSchema( states[0], comm, communicator );
  comm.checked( "coupled coarse surface buckets",
                [&]
  {
    for( std::size_t s = 0; s < states.size(); ++s )
    {
      for( std::size_t c = 0; c < states[s].cells.size(); ++c )
      {
        if( isSurface( states[s].cells[c] ) )
        {
          CoarseSurface query;
          for( vtkIdType p : states[s].cells[c].points )
          {
            query.cornerBuckets.push_back( s == 0 ? Connectivity{ states[s].pointIds[p] } : states[s].buckets[p] );
          }
          surfaces.push_back( std::move( query ) );
          surfaceIndices.emplace_back( s, c );
        }
      }
    }
  } );
  auto const globalSurfaceCount =
    MpiWrapper::allReduce( static_cast< std::uint64_t >( surfaces.size() ), MpiWrapper::Reduction::Sum, communicator );
  std::vector< std::vector< SurfaceSide > > associations;
  if( globalSurfaceCount )
  {
    associations = comm.discoverSurfaceSides( mainFaces, surfaces );
  }
  // Keep actual volume-side owners separate from surface replica participants.
  std::vector< std::vector< Participants > > surfaceParticipants;
  comm.checked( "coupled coarse surface installation",
                [&]
  {
    surfaceParticipants.resize( states.size() );
    for( std::size_t s = 0; s < states.size(); ++s )
    {
      if( states[s].hasSurfaces )
      {
        surfaceParticipants[s].resize( states[s].cells.size() );
      }
    }
    for( std::size_t q = 0; q < surfaceIndices.size(); ++q )
    {
      auto const [s, c] = surfaceIndices[q];
      surfaceParticipants[s][c] = states[s].sides[c].front().owners;
      states[s].sides[c] = std::move( associations[q] );
      if( s == 0 && !std::binary_search( states[s].sides[c].front().owners.begin(), states[s].sides[c].front().owners.end(),
                                         comm.rank() ) )
      {
        throw std::invalid_argument( "Main surface marker is not on a locally owned volume face" );
      }
    }
  } );
  // Coarse-only face and routing payloads are released before fine growth.
  std::vector< MainFace >{}.swap( mainFaces );
  std::vector< CoarseSurface >{}.swap( surfaces );
  std::vector< std::vector< SurfaceSide > >{}.swap( associations );

  RefinementLevelStatistics running;
  comm.checked( "resource forecast and coarse statistics", [&]
  {
    GEOS_MARK_SCOPE_STR( "uniformRefinement/resourceForecast" );
    std::uint64_t retainedBytes = multiply( meshes.getMainMesh()->GetActualMemorySize(), 1024 );
    for( auto const & [name, mesh] : meshes.getFaceBlocks() )
    {
      GEOS_UNUSED_VAR( name );
      if( mesh )
        retainedBytes = add( retainedBytes, multiply( mesh->GetActualMemorySize(), 1024 ) );
    }
    result.resources = resourceForecast( states, levels, retainedBytes, comm.neighbors().size(), maximumConfiguration[5] );
    result.levels.reserve( levels );
    running.ownedVolumeCells = localVolumes;
    running.mainPointCopies = running.ownedMainPoints = states[0].coordinates.size();
    for( auto const & entity : states[0].interfaces )
    {
      if( entity.key.kind == EntityKind::vertex )
      {
        ++running.sharedMainPointCopies;
        running.ownedMainPoints -= entity.participants.front() != comm.rank();
      }
    }
    result.coarseCommunication = running.communication = comm.statistics();
  } );
  reportForecast( result.resources, options.reportStatistics, comm, communicator );
  reportLevel( running, 0, options.reportStatistics, comm, communicator );
  auto previousCommunication = result.coarseCommunication;
  for( int generation = 1; generation <= levels; ++generation )
  {
    // Coarse reconciliation checked full layouts on every shared vertex and
    // surface replica. Transfer preserves those layouts and fine participants
    // inherit their coarse supports, so fine records need only typed values.
    std::string const levelScope = "uniformRefinement/level" + std::to_string( generation );
    GEOS_MARK_SCOPE_STR( levelScope.c_str() );
    std::vector< Level > plans;
    comm.checked( "coupled level planning",
                  [&]
    {
      plans.reserve( states.size() );
      for( auto const & state : states )
      {
        plans.push_back( planLevel( state, options, generation, comm.rank() ) );
      }
      for( std::size_t s = 0; s < plans.size(); ++s )
      {
        for( std::size_t c = 0; c < plans[s].next.cells.size(); ++c )
        {
          if( isSurface( plans[s].next.cells[c] ) )
          {
            plans[s].cellParticipants[c] = surfaceParticipants[s][plans[s].parentIndices[c]];
          }
        }
      }
    } );
    // Report every parent that cannot be refined, on every rank, before any
    // further communication uses the incomplete plans.
    std::uint64_t localInvalid = 0;
    for( auto const & plan : plans )
    {
      localInvalid += plan.invalidCells;
    }
    std::uint64_t const globalInvalid = MpiWrapper::sum( localInvalid, communicator );
    comm.checked( "cell geometry",
                  [&]
    {
      if( localInvalid == 0 )
      {
        return;
      }
      std::string message = GEOS_FMT( "{} cell(s) cannot be refined at level {} ({} on this rank). "
                                      "Uniform refinement stops instead of making an invalid cell. "
                                      "Improve these cells in the coarse mesh, or use uniformRefinement=\"{}\". "
                                      "First cells on this rank:",
                                      globalInvalid, generation, localInvalid, generation - 1 );
      for( auto const & plan : plans )
      {
        for( auto const & reason : plan.invalidReasons )
        {
          message += "\n  - " + reason;
        }
      }
      throw std::invalid_argument( message );
    } );
    // Ranks sharing a new point may see different boundary faces around it.
    // The point joins a node set when any of its participants places it there.
    std::vector< SharedFlag > nodeSetFlags;
    comm.checked( "shared node-set flags",
                  [&]
    {
      for( auto const & plan : plans )
      {
        for( vtkIdType p = plan.points->originalSize(); p < static_cast< vtkIdType >( plan.points->points().size() ); ++p )
        {
          std::uint64_t label = 0;
          for( auto const & [name, members] : plan.nodeSets )
          {
            GEOS_UNUSED_VAR( name );
            if( members[p] && plan.interfaces->participants( { p } ).size() > 1 )
            {
              nodeSetFlags.push_back( { plan.points->points()[p].key, plan.interfaces->participants( { p } ), label } );
            }
            ++label;
          }
        }
      }
    } );
    auto const sharedNodeSets = comm.unionSharedFlags( generation, nodeSetFlags );
    comm.checked( "coupled point fields",
                  [&]
    {
      GEOS_MARK_SCOPE_STR( "uniformRefinement/transferFields" );
      for( std::size_t s = 0; s < plans.size(); ++s )
      {
        auto & plan = plans[s];
        std::vector< std::vector< bool > * > byLabel;
        for( auto & [name, members] : plan.nodeSets )
        {
          GEOS_UNUSED_VAR( name );
          byLabel.push_back( &members );
        }
        for( auto const & [key, label] : sharedNodeSets )
        {
          if( key.meshNamespace != plan.next.ns )
          {
            continue;
          }
          vtkIdType const p = plan.points->findPoint( key );
          if( p < plan.points->originalSize() || label >= byLabel.size() )
          {
            throw std::invalid_argument( "Shared node-set flag for an unknown refinement point or set" );
          }
          ( *byLabel[label] )[p] = true;
        }
        plan.next.pointData = transferPointData( *states[s].pointData, *plan.points, policiesFor( states[s], options ), plan.nodeSets );
      }
    } );
    std::unordered_map< EntityKey, vtkIdType, EntityKeyHash > mainSupportIds;
    // Point ranges are per dataset namespace. Main IDs and each auxiliary point
    // space preserve their own old IDs; auxiliary cell ranges are separate below.
    for( std::size_t s = 0; s < plans.size(); ++s )
    {
      auto & plan = plans[s];
      std::vector< PointCreation > creations;
      comm.checked( "coupled point tuples",
                    [&]
      {
        PointFieldLayout fields( *plan.next.pointData );
        for( vtkIdType p = plan.points->originalSize(); p < static_cast< vtkIdType >( plan.points->points().size() ); ++p )
        {
          auto const & recipe = plan.points->points()[p];
          auto const & participants = plan.interfaces->participants( { p } );
          creations.push_back(
            { recipe.key, participants, recipe.position,
              participants.size() > 1 ? fields.pack( p, FieldTupleFormat::valuesOnly ) : Bytes{},
              supportExtent( states[s].coordinates, recipe.support ) } );
          if( s == 0 )
          {
            running.ownedMainPoints += participants.front() == comm.rank();
            running.sharedMainPointCopies += participants.size() > 1;
          }
        }
      } );
      vtkIdType const maximum = states[s].pointIds.empty() ? -1 : *std::max_element( states[s].pointIds.begin(), states[s].pointIds.end() );
      auto const records = comm.resolvePoints( generation, creations, maximum );
      comm.checked( "coupled point installation",
                    [&]
      {
        auto & next = plan.next;
        next.pointIds = states[s].pointIds;
        next.coordinates.reserve( plan.points->points().size() );
        PointFieldLayout fields( *next.pointData );
        for( vtkIdType p = 0; p < static_cast< vtkIdType >( plan.points->points().size() ); ++p )
        {
          auto const & recipe = plan.points->points()[p];
          if( p < plan.points->originalSize() )
          {
            next.coordinates.push_back( recipe.position );
          }
          else
          {
            auto const & record = records.at( recipe.key );
            next.pointIds.push_back( record.globalId );
            next.coordinates.push_back( record.position );
            if( plan.interfaces->participants( { p } ).size() > 1 )
            {
              fields.install( p, record.fields, FieldTupleFormat::valuesOnly );
            }
          }
          // Only auxiliary collocation requests use this map.
          if( s == 0 && plans.size() > 1 )
          {
            mainSupportIds.emplace( recipe.key, next.pointIds[p] );
          }
        }
        next.interfaces = plan.interfaces->fineSupports( next.pointIds );
      } );
    }
    std::uint64_t volumeCount = 0;
    comm.checked( "coupled volume counts",
                  [&]
    {
      for( auto const & cell : plans[0].next.cells )
      {
        volumeCount += !isSurface( cell );
      }
    } );
    auto const volumeRange = comm.allocateRange( volumeCount, 0 );
    running.ownedVolumeCells = volumeCount;
    running.mainPointCopies = plans[0].points->points().size();
    comm.checked( "coupled volume IDs",
                  [&]
    {
      vtkIdType id = volumeRange.first;
      for( std::size_t c = 0; c < plans[0].next.cells.size(); ++c )
      {
        if( !isSurface( plans[0].next.cells[c] ) )
        {
          plans[0].next.cellIds[c] = id++;
        }
      }
    } );
    vtkIdType surfaceBase = 0;
    comm.checked( "main surface prefix",
                  [&]
    {
      if( volumeRange.total > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) )
      {
        throw std::overflow_error( "Volume prefix leaves no representable surface base" );
      }
      surfaceBase = static_cast< vtkIdType >( volumeRange.total );
    } );
    vtkIdType auxiliaryBase = 0;
    std::vector< std::vector< Participants > > nextParticipants;
    comm.checked( "surface participant allocation", [&] { nextParticipants.resize( states.size() ); } );
    for( std::size_t s = 0; s < plans.size(); ++s )
    {
      auto & plan = plans[s];
      std::vector< CellCreation > creations;
      comm.checked( "coupled surface-cell tuples",
                    [&]
      {
        CellFieldLayout fields( *plan.next.cellData );
        if( plan.next.hasSurfaces )
        {
          nextParticipants[s].resize( plan.next.cells.size() );
        }
        for( std::size_t c = 0; c < plan.next.cells.size(); ++c )
        {
          if( isSurface( plan.next.cells[c] ) )
          {
            creations.push_back( { plan.cellKeys[c], plan.cellParticipants[c],
                                   plan.cellParticipants[c].size() > 1 ? fields.pack( c, FieldTupleFormat::valuesOnly ) : Bytes{} } );
            nextParticipants[s][c] = plan.cellParticipants[c];
          }
        }
      } );
      vtkIdType const base = s == 0 ? surfaceBase : auxiliaryBase;
      IdRange range{};
      auto const records = comm.resolveCells( generation, creations, base, &range );
      comm.checked( "coupled surface-cell installation",
                    [&]
      {
        CellFieldLayout fields( *plan.next.cellData );
        for( std::size_t c = 0; c < plan.next.cells.size(); ++c )
        {
          if( isSurface( plan.next.cells[c] ) )
          {
            auto const & record = records.at( plan.cellKeys[c] );
            plan.next.cellIds[c] = record.globalId;
            if( plan.cellParticipants[c].size() > 1 )
            {
              fields.install( c, record.fields, FieldTupleFormat::valuesOnly );
            }
          }
        }
        if( range.total > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) - base )
        {
          throw std::overflow_error( "Auxiliary cell prefix exceeds vtkIdType" );
        }
        if( s )
        {
          auxiliaryBase = base + static_cast< vtkIdType >( range.total );
        }
      } );
    }
    // Rebuild every auxiliary child side and bucket from actual parent faces.
    // No product of endpoint buckets and no coordinate identity is used.
    for( std::size_t s = 1; s < plans.size(); ++s )
    {
      auto & plan = plans[s];
      std::vector< SupportLookup > requests;
      using Binding = std::tuple< vtkIdType, std::size_t, vtkIdType >;
      std::map< Binding, std::size_t > bindings;
      comm.checked( "coupled auxiliary support planning",
                    [&]
      {
        GEOS_MARK_SCOPE_STR( "uniformRefinement/fractureAssociations" );
        for( std::size_t c = 0; c < plan.next.cells.size(); ++c )
        {
          vtkIdType const parent = plan.parentIndices[c];
          for( std::size_t side = 0; side < states[s].sides[parent].size(); ++side )
          {
            for( vtkIdType point : plan.next.cells[c].points )
            {
              auto const [found, inserted] = bindings.emplace( Binding{ parent, side, point }, requests.size() );
              GEOS_UNUSED_VAR( found );
              if( !inserted )
              {
                continue;
              }
              auto const & recipe = plan.points->points()[point];
              auto const & oldSide = states[s].sides[parent][side];
              requests.push_back( { SurfaceAssociations::supportKey( recipe.key.kind, recipe.support,
                                                                     states[s].cells[parent].points, oldSide ),
                                    oldSide.owners.front() } );
            }
          }
        }
      } );
      auto const ids = comm.resolveSupportIds( generation, requests, mainSupportIds );
      comm.checked( "coupled auxiliary associations",
                    [&]
      {
        GEOS_MARK_SCOPE_STR( "uniformRefinement/fractureAssociations" );
        plan.next.buckets = states[s].buckets;
        plan.next.buckets.resize( plan.next.pointIds.size() );
        for( auto const & [binding, request] : bindings )
        {
          auto const [parent, side, point] = binding;
          GEOS_UNUSED_VAR( parent, side );
          if( point >= plan.points->originalSize() )
          {
            plan.next.buckets[point].push_back( ids[request] );
          }
        }
        for( auto & bucket : plan.next.buckets )
        {
          std::sort( bucket.begin(), bucket.end() );
          bucket.erase( std::unique( bucket.begin(), bucket.end() ), bucket.end() );
        }
        for( std::size_t c = 0; c < plan.next.cells.size(); ++c )
        {
          vtkIdType const parent = plan.parentIndices[c];
          for( std::size_t side = 0; side < states[s].sides[parent].size(); ++side )
          {
            Connectivity corners;
            for( vtkIdType point : plan.next.cells[c].points )
            {
              corners.push_back( ids[bindings.at( { parent, side, point } )] );
            }
            plan.next.sides[c].push_back(
              { entityKey( EntityKind::face, corners ), std::move( corners ), states[s].sides[parent][side].owners } );
          }
        }
      } );
    }
    comm.checked( "coupled level output",
                  [&]
    {
      GEOS_MARK_SCOPE_STR( "uniformRefinement/validation" );
      for( auto const & plan : plans )
      {
        PointRegistry finalPoints( plan.next.coordinates, plan.next.pointIds, plan.next.ns );
        for( std::size_t c = 0; c < plan.next.cells.size(); ++c )
        {
          auto const & cell = plan.next.cells[c];
          withCellContext( plan.next.name, generation, plan.next.parents[c], cell.vtkType, [&]
          {
            if( !isSurface( cell ) )
              validateGeometry( cell, finalPoints );
          } );
        }
        if( std::find( plan.next.cellIds.begin(), plan.next.cellIds.end(), -1 ) != plan.next.cellIds.end() )
        {
          throw std::invalid_argument( "Unassigned refined cell ID" );
        }
      }
    } );
    // No old VTK dataset is committed on failure; parent states are released
    // after this collectively successful level, rather than retaining L meshes.
    comm.checked( "level statistics installation", [&]
    {
      running.communication = communicationDelta( comm.statistics(), previousCommunication );
      result.levels.push_back( running );
    } );
    previousCommunication = comm.statistics();
    reportLevel( running, generation, options.reportStatistics, comm, communicator );
    states.clear();
    for( auto & plan : plans )
    {
      states.push_back( std::move( plan.next ) );
    }
    surfaceParticipants = std::move( nextParticipants );
  }
  std::vector< vtkSmartPointer< vtkUnstructuredGrid > > output;
  vtkIdType localMainMaximum = -1, localAuxiliaryMaximum = -1;
  for( vtkIdType id : states[0].cellIds )
  {
    localMainMaximum = std::max( localMainMaximum, id );
  }
  for( std::size_t s = 1; s < states.size(); ++s )
  {
    for( vtkIdType id : states[s].cellIds )
    {
      localAuxiliaryMaximum = std::max( localAuxiliaryMaximum, id );
    }
  }
  auto const mainMaximum = MpiWrapper::allReduce( localMainMaximum, MpiWrapper::Reduction::Max, communicator );
  auto const auxiliaryMaximum = MpiWrapper::allReduce( localAuxiliaryMaximum, MpiWrapper::Reduction::Max, communicator );
  comm.checked( "final auxiliary import offset",
                [&]
  {
    if( auxiliaryMaximum >= 0 )
    {
      auto const limit = std::min( static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ),
                                   static_cast< std::uint64_t >( std::numeric_limits< globalIndex >::max() ) );
      if( mainMaximum < 0 || add( add( static_cast< std::uint64_t >( mainMaximum ), 1 ), auxiliaryMaximum ) > limit )
      {
        throw std::overflow_error( "Final fracture import offset exceeds global ID storage" );
      }
    }
  } );
  for( std::size_t s = 0; s < states.size(); ++s )
  {
    int capacity = 0;
    comm.checked( "collocation capacity",
                  [&]
    {
      for( auto const & bucket : states[s].buckets )
      {
        if( bucket.size() > static_cast< std::size_t >( INT_MAX ) )
        {
          throw std::overflow_error( "Collocation capacity exceeds int" );
        }
        capacity = std::max( capacity, static_cast< int >( bucket.size() ) );
      }
    } );
    capacity = MpiWrapper::allReduce( capacity, MpiWrapper::Reduction::Max, communicator );
    comm.checked( "coupled final VTK construction", [&] { output.push_back( buildGrid( states[s], std::max( capacity, 1 ) ) ); } );
  }
  stdMap< string, vtkSmartPointer< vtkDataSet > > blocks;
  comm.checked( "coupled final commit preparation",
                [&]
  {
    for( std::size_t s = 1; s < states.size(); ++s )
    {
      blocks.emplace( states[s].name, output[s] );
    }
    result.blocks = blockDescriptors( states[0], sources, levels );
    result.neighbors = comm.neighbors();
    result.communication = comm.statistics();
  } );
  meshes.getFaceBlocks().swap( blocks );
  meshes.setMainMesh( output[0] );
  return result;
}
} // namespace geos::vtk
