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

#include "mesh/generators/VTKMeshScattering.hpp"
#ifdef GEOS_USE_MPI
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#endif

#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"
#include "common/MpiWrapper.hpp"
#include "common/MpiChunkedCommunication.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/math.hpp"
#include "LvArray/src/system.hpp"

#include <vtkAbstractArray.h>
#include <vtkBitArray.h>
#include <vtkFieldData.h>
#include <vtkStringArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkExtractCells.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPartitionedDataSet.h>
#include <vtkPoints.h>
#include <vtkVersionMacros.h>
#ifdef GEOS_USE_MPI
#include <vtkRedistributeDataSetFilter.h>
#include <vtkMultiProcessController.h>
#if VTK_VERSION_NUMBER == VTK_VERSION_CHECK( 9, 7, 0 )
#include <vtkBoundingBox.h>
#include <vtkCellCenters.h>
#include <vtkDIYKdTreeUtilities.h>
#endif
#include <vtkMPIController.h>
#include <vtkMPI.h>
#endif
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>

namespace geos
{
namespace vtk
{

namespace
{

// ============================================================================
// Buffer serialization helpers
// ============================================================================

void appendBytes( stdVector< char > & buf, void const * data, int64_t n )
{
  int64_t const pos = static_cast< int64_t >( buf.size() );
  buf.resize( static_cast< size_t >( pos + n ) );
  std::memcpy( buf.data() + pos, data, static_cast< size_t >( n ) );
}

template< typename T >
void appendValue( stdVector< char > & buf, T val )
{
  appendBytes( buf, &val, sizeof( T ) );
}

template< typename T >
T readValue( char const * & ptr )
{
  T val;
  std::memcpy( &val, ptr, sizeof( T ) );
  ptr += sizeof( T );
  return val;
}

// ============================================================================
// MPI large-message helpers (handles buffers > 2 GB)
// ============================================================================

void mpiSendLarge( void const * buf, int64_t count, integer dest, integer tag, MPI_Comm comm )
{
  GEOS_ERROR_IF( count < 0, "Negative mesh-scatter byte count" );
  mpi::sendBytes( buf, static_cast< std::uint64_t >( count ), dest, tag, comm );
}

void mpiRecvLarge( void * buf, int64_t count, integer src, integer tag, MPI_Comm comm )
{
  GEOS_ERROR_IF( count < 0, "Negative mesh-scatter byte count" );
  mpi::receiveBytes( buf, static_cast< std::uint64_t >( count ), src, tag, comm );
}

// ============================================================================
// Pack / unpack vtkDataSetAttributes (cell data or point data)
// ============================================================================

constexpr int NUM_ATTR_TYPES = vtkDataSetAttributes::NUM_ATTRIBUTES;

void packString( stdVector< char > & buf, char const * value )
{
  int64_t const size = value ? static_cast< int64_t >( std::strlen( value ) ) : -1;
  appendValue( buf, size );
  if( size > 0 )
    appendBytes( buf, value, size );
}

void packDataArrays( stdVector< char > & buf, vtkFieldData * data )
{
  appendValue( buf, static_cast< int32_t >( data->GetNumberOfArrays() ) );
  for( int a = 0; a < data->GetNumberOfArrays(); ++a )
  {
    vtkAbstractArray * arr = data->GetAbstractArray( a );
    packString( buf, arr->GetName() );
    appendValue( buf, static_cast< int32_t >( arr->GetNumberOfComponents() ) );
    appendValue( buf, static_cast< int32_t >( arr->GetDataType() ) );
    appendValue( buf, static_cast< int64_t >( arr->GetNumberOfTuples() ) );
    for( int c = 0; c < arr->GetNumberOfComponents(); ++c )
      packString( buf, arr->GetComponentName( c ) );
    if( auto * strings = vtkStringArray::SafeDownCast( arr ) )
    {
      for( vtkIdType v = 0; v < arr->GetNumberOfValues(); ++v )
      {
        auto const & value = strings->GetValue( v );
        appendValue( buf, static_cast< int64_t >( value.size() ) );
        if( !value.empty() )
          appendBytes( buf, value.data(), value.size() );
      }
    }
    else if( auto * bits = vtkBitArray::SafeDownCast( arr ) )
    {
      for( vtkIdType v = 0; v < arr->GetNumberOfValues(); ++v )
        appendValue( buf, static_cast< unsigned char >( bits->GetValue( v ) ) );
    }
    else
    {
      auto * numeric = vtkDataArray::SafeDownCast( arr );
      GEOS_ERROR_IF( numeric == nullptr, "Unsupported abstract array in mesh scatter" );
      int64_t const count = arr->GetNumberOfValues();
      int64_t const width = arr->GetDataTypeSize();
      GEOS_ERROR_IF( width <= 0 || count > std::numeric_limits< int64_t >::max() / width,
                     "Mesh scatter array byte count overflow" );
      if( count > 0 )
        appendBytes( buf, numeric->GetVoidPointer( 0 ), count * width );
    }
  }
  // Indices preserve roles even for unnamed arrays or string pedigree IDs.
  auto * attrs = vtkDataSetAttributes::SafeDownCast( data );
  int roles[NUM_ATTR_TYPES];
  std::fill_n( roles, NUM_ATTR_TYPES, -1 );
  if( attrs )
    attrs->GetAttributeIndices( roles );
  for( int role : roles )
    appendValue( buf, static_cast< int32_t >( role ) );
}

std::pair< bool, string > unpackString( char const * & ptr )
{
  int64_t const size = readValue< int64_t >( ptr );
  if( size < 0 )
    return { false, {} };
  string value( ptr, size );
  ptr += size;
  return { true, std::move( value ) };
}

void unpackDataArrays( char const * & ptr, vtkFieldData * data )
{
  int32_t const nArrays = readValue< int32_t >( ptr );
  for( int a = 0; a < nArrays; ++a )
  {
    auto const name = unpackString( ptr );
    int32_t const nComp = readValue< int32_t >( ptr );
    int32_t const dataType = readValue< int32_t >( ptr );
    int64_t const nTuples = readValue< int64_t >( ptr );
    vtkSmartPointer< vtkAbstractArray > arr;
    arr.TakeReference( vtkAbstractArray::CreateArray( dataType ) );
    GEOS_ERROR_IF( arr == nullptr || nComp <= 0 || nTuples < 0,
                   "Invalid mesh scatter array metadata" );
    if( name.first )
      arr->SetName( name.second.c_str() );
    arr->SetNumberOfComponents( nComp );
    arr->SetNumberOfTuples( nTuples );
    for( int c = 0; c < nComp; ++c )
    {
      auto const component = unpackString( ptr );
      if( component.first )
        arr->SetComponentName( c, component.second.c_str() );
    }
    if( auto * strings = vtkStringArray::SafeDownCast( arr ) )
    {
      for( vtkIdType v = 0; v < arr->GetNumberOfValues(); ++v )
      {
        int64_t const size = readValue< int64_t >( ptr );
        GEOS_ERROR_IF( size < 0, "Invalid mesh scatter string length" );
        strings->SetValue( v, string( ptr, size ) );
        ptr += size;
      }
    }
    else if( auto * bits = vtkBitArray::SafeDownCast( arr ) )
    {
      for( vtkIdType v = 0; v < arr->GetNumberOfValues(); ++v )
      {
        auto const value = readValue< unsigned char >( ptr );
        GEOS_ERROR_IF( value > 1, "Invalid mesh scatter bit value" );
        bits->SetValue( v, value );
      }
    }
    else
    {
      auto * numeric = vtkDataArray::SafeDownCast( arr );
      GEOS_ERROR_IF( numeric == nullptr, "Unsupported abstract array in mesh scatter" );
      int64_t const count = arr->GetNumberOfValues();
      int64_t const width = arr->GetDataTypeSize();
      GEOS_ERROR_IF( width <= 0 || count > std::numeric_limits< int64_t >::max() / width,
                     "Mesh scatter array byte count overflow" );
      if( count > 0 )
        std::memcpy( numeric->GetVoidPointer( 0 ), ptr, count * width );
      ptr += count * width;
    }
    data->AddArray( arr );
  }
  auto * attrs = vtkDataSetAttributes::SafeDownCast( data );
  for( int t = 0; t < NUM_ATTR_TYPES; ++t )
  {
    int const role = readValue< int32_t >( ptr );
    if( attrs && role >= 0 )
      attrs->SetActiveAttribute( role, t );
  }
}

// ============================================================================
// Pack / unpack a vtkCellArray
//
// The offsets and connectivity of a vtkCellArray may be stored as either 32 or
// 64 bit integers depending on how the mesh was built; ConvertToDefaultStorage()
// normalizes them to vtkIdType so that the raw buffers can be copied as-is.
// A null array (e.g. the face arrays of a mesh without polyhedra) is encoded
// with a negative size.
// ============================================================================

void appendCellArray( stdVector< char > & buf, vtkCellArray * cells )
{
  if( cells == nullptr )
  {
    appendValue( buf, int64_t( -1 ) );
    return;
  }

  cells->ConvertToDefaultStorage();

  vtkDataArray * offsets = cells->GetOffsetsArray();
  vtkDataArray * conn = cells->GetConnectivityArray();

  int64_t const nOffsets = offsets->GetNumberOfValues();
  int64_t const connSize = conn->GetNumberOfValues();

  appendValue( buf, nOffsets );
  appendValue( buf, connSize );
  appendBytes( buf, offsets->GetVoidPointer( 0 ), nOffsets * sizeof( vtkIdType ) );
  appendBytes( buf, conn->GetVoidPointer( 0 ), connSize * sizeof( vtkIdType ) );
}

vtkSmartPointer< vtkCellArray > readCellArray( char const * & ptr )
{
  int64_t const nOffsets = readValue< int64_t >( ptr );
  if( nOffsets < 0 )
  {
    return nullptr;
  }

  int64_t const connSize = readValue< int64_t >( ptr );

  vtkNew< vtkIdTypeArray > offsets;
  offsets->SetNumberOfValues( nOffsets );
  std::memcpy( offsets->GetVoidPointer( 0 ), ptr, nOffsets * sizeof( vtkIdType ) );
  ptr += nOffsets * sizeof( vtkIdType );

  vtkNew< vtkIdTypeArray > conn;
  conn->SetNumberOfValues( connSize );
  std::memcpy( conn->GetVoidPointer( 0 ), ptr, connSize * sizeof( vtkIdType ) );
  ptr += connSize * sizeof( vtkIdType );

  auto cells = vtkSmartPointer< vtkCellArray >::New();
  cells->SetData( offsets, conn );
  return cells;
}

// ============================================================================
// Pack / unpack a full vtkUnstructuredGrid + assignment vector
// ============================================================================
void packGrid( vtkUnstructuredGrid * grid,
               stdVector< integer > const & assignment,
               stdVector< char > & buf )
{
  buf.clear();

  int64_t const nPoints = grid->GetNumberOfPoints();
  int64_t const nCells = grid->GetNumberOfCells();
  appendValue( buf, nPoints );
  appendValue( buf, nCells );

  // Points (always serialized as real64)
  if( nPoints > 0 )
  {
    vtkPoints * points = grid->GetPoints();
    if( auto * coords = vtkDoubleArray::SafeDownCast( points->GetData() ) )
    {
      appendBytes( buf, coords->GetPointer( 0 ), nPoints * 3 * sizeof( real64 ) );
    }
    else
    {
      for( int64_t i = 0; i < nPoints; ++i )
      {
        real64 p[3];
        points->GetPoint( i, p );
        appendBytes( buf, p, sizeof( p ) );
      }
    }
  }

  // Cell types, offsets, connectivity
  if( nCells > 0 )
  {
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK( 9, 6, 0 )
    vtkUnsignedCharArray * const types = vtkUnsignedCharArray::SafeDownCast( grid->GetCellTypes() );
#else
    vtkUnsignedCharArray * const types = grid->GetCellTypesArray();
#endif
    if( types != nullptr )
    {
      appendBytes( buf, types->GetPointer( 0 ), nCells * sizeof( unsigned char ) );
    }
    else
    {
      stdVector< unsigned char > cellTypes( nCells );
      for( int64_t i = 0; i < nCells; ++i )
      {
        cellTypes[i] = static_cast< unsigned char >( grid->GetCellType( i ) );
      }
      appendBytes( buf, cellTypes.data(), nCells * sizeof( unsigned char ) );
    }

    appendCellArray( buf, grid->GetCells() );

    // VTK_POLYHEDRON cells are not fully described by the connectivity array: their
    // face description lives in two separate arrays that must travel with the mesh.
    // Both are null for meshes without polyhedra.
    appendCellArray( buf, grid->GetPolyhedronFaceLocations() );
    appendCellArray( buf, grid->GetPolyhedronFaces() );
  }

  // Field data arrays
  packDataArrays( buf, grid->GetCellData() );
  packDataArrays( buf, grid->GetPointData() );
  packDataArrays( buf, grid->GetFieldData() );

  // Assignment vector
  int64_t const assignSize = static_cast< int64_t >( assignment.size() );
  appendValue( buf, assignSize );
  if( assignSize > 0 )
  {
    appendBytes( buf, assignment.data(), assignSize * sizeof( integer ) );
  }
}


std::pair< vtkSmartPointer< vtkUnstructuredGrid >, stdVector< integer > >
unpackGrid( stdVector< char > const & buf )
{
  char const * ptr = buf.data();

  int64_t const nPoints = readValue< int64_t >( ptr );
  int64_t const nCells = readValue< int64_t >( ptr );

  auto grid = vtkSmartPointer< vtkUnstructuredGrid >::New();

  // Points
  if( nPoints > 0 )
  {
    vtkNew< vtkPoints > points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints( nPoints );
    vtkDoubleArray * const coords = vtkDoubleArray::SafeDownCast( points->GetData() );
    std::memcpy( coords->GetPointer( 0 ), ptr, nPoints * 3 * sizeof( real64 ) );
    ptr += nPoints * 3 * sizeof( real64 );
    grid->SetPoints( points );
  }

  // Cells
  if( nCells > 0 )
  {
    vtkNew< vtkUnsignedCharArray > types;
    types->SetNumberOfValues( nCells );
    std::memcpy( types->GetVoidPointer( 0 ), ptr, nCells * sizeof( unsigned char ) );
    ptr += nCells * sizeof( unsigned char );

    vtkSmartPointer< vtkCellArray > cellArray = readCellArray( ptr );
    vtkSmartPointer< vtkCellArray > faceLocations = readCellArray( ptr );
    vtkSmartPointer< vtkCellArray > faces = readCellArray( ptr );

    if( faces == nullptr )
    {
      // SetCells() must not be used when polyhedra are present: it would reinterpret
      // the connectivity of those cells as a face stream and read out of bounds.
      unsigned char const * const typeBegin = types->GetPointer( 0 );
      GEOS_ERROR_IF( std::find( typeBegin, typeBegin + nCells, VTK_POLYHEDRON ) != typeBegin + nCells,
                     "Mesh scattering: polyhedral cells were received without their face description." );
      grid->SetCells( types, cellArray );
    }
    else
    {
      grid->SetPolyhedralCells( types, cellArray, faceLocations, faces );
    }
  }

  // Field data
  unpackDataArrays( ptr, grid->GetCellData() );
  unpackDataArrays( ptr, grid->GetPointData() );
  unpackDataArrays( ptr, grid->GetFieldData() );

  // Assignment
  int64_t const assignSize = readValue< int64_t >( ptr );
  stdVector< integer > assignment( assignSize );
  if( assignSize > 0 )
  {
    std::memcpy( assignment.data(), ptr, assignSize * sizeof( integer ) );
    ptr += assignSize * sizeof( integer );
  }

  return { grid, std::move( assignment ) };
}


// ============================================================================
// Split cells into two subsets: those with assignment < mid (kept locally)
// and those with assignment >= mid (sent to the partner rank).
// ============================================================================

struct SplitResult
{
  vtkSmartPointer< vtkUnstructuredGrid > loMesh;
  stdVector< integer > loAssignment;
  vtkSmartPointer< vtkUnstructuredGrid > hiMesh;
  stdVector< integer > hiAssignment;
};

SplitResult splitByMid( vtkUnstructuredGrid * mesh,
                        stdVector< integer > const & assignment,
                        integer mid )
{
  stdVector< vtkIdType > loCells, hiCells;
  stdVector< integer > loAssign, hiAssign;
  loCells.reserve( assignment.size() / 2 );
  hiCells.reserve( assignment.size() / 2 );
  loAssign.reserve( assignment.size() / 2 );
  hiAssign.reserve( assignment.size() / 2 );

  for( vtkIdType i = 0; i < static_cast< vtkIdType >( assignment.size() ); ++i )
  {
    if( assignment[i] < mid )
    {
      loCells.push_back( i );
      loAssign.push_back( assignment[i] );
    }
    else
    {
      hiCells.push_back( i );
      hiAssign.push_back( assignment[i] );
    }
  }

  auto extract = []( vtkUnstructuredGrid * inputMesh,
                     stdVector< vtkIdType > const & ids )
                 -> vtkSmartPointer< vtkUnstructuredGrid >
  {
    if( ids.empty() )
    {
      return vtkSmartPointer< vtkUnstructuredGrid >::New();
    }
    vtkNew< vtkExtractCells > extractor;
    extractor->SetInputData( inputMesh );
    extractor->SetCellIds( ids.data(), static_cast< vtkIdType >( ids.size() ) );
    extractor->Update();
    auto result = vtkSmartPointer< vtkUnstructuredGrid >::New();
    result->ShallowCopy( extractor->GetOutput() );
    return result;
  };

  // Extract the smaller subset first, then release the original mesh reference
  // before extracting the larger one to reduce peak memory.
  SplitResult result;
  if( hiCells.size() <= loCells.size() )
  {
    result.hiMesh = extract( mesh, hiCells );
    result.hiAssignment = std::move( hiAssign );
    result.loMesh = extract( mesh, loCells );
    result.loAssignment = std::move( loAssign );
  }
  else
  {
    result.loMesh = extract( mesh, loCells );
    result.loAssignment = std::move( loAssign );
    result.hiMesh = extract( mesh, hiCells );
    result.hiAssignment = std::move( hiAssign );
  }
  return result;
}


// ============================================================================
// Centroid computation
// ============================================================================

stdVector< stdArray< real64, 3 > >
computeCentroids( vtkDataSet & mesh )
{
  vtkIdType const n = mesh.GetNumberOfCells();
  stdVector< stdArray< real64, 3 > > centroids( n );
  vtkNew< vtkIdList > ptIds;

  for( vtkIdType c = 0; c < n; ++c )
  {
    mesh.GetCellPoints( c, ptIds );
    real64 cx = 0.0, cy = 0.0, cz = 0.0;
    vtkIdType const nPts = ptIds->GetNumberOfIds();
    for( vtkIdType i = 0; i < nPts; ++i )
    {
      real64 p[3];
      mesh.GetPoint( ptIds->GetId( i ), p );
      cx += p[0];
      cy += p[1];
      cz += p[2];
    }
    if( nPts > 0 )
    {
      real64 const inv = 1.0 / nPts;
      centroids[c] = { cx * inv, cy * inv, cz * inv };
    }
    else
    {
      centroids[c] = { 0.0, 0.0, 0.0 };
    }
  }
  return centroids;
}


// ============================================================================
// Cell-rank assignment: contiguous (index-based, no geometry)
// ============================================================================

stdVector< integer >
computeCellRanksContiguous( vtkIdType nCells, integer size )
{
  stdVector< integer > ranks( nCells );
  vtkIdType const perRank = nCells / size;
  vtkIdType const remainder = nCells % size;

  // Ranks [0, remainder) get (perRank+1) cells, the rest get perRank
  vtkIdType cell = 0;
  for( integer r = 0; r < size; ++r )
  {
    vtkIdType const count = perRank + ( r < remainder ? 1 : 0 );
    for( vtkIdType i = 0; i < count; ++i )
    {
      ranks[cell++] = r;
    }
  }
  return ranks;
}


// ============================================================================
// Cell-rank assignment: Cartesian grid partition
// ============================================================================

stdVector< integer >
computeCellRanksCartesian( vtkDataSet & mesh, integer nx, integer ny, integer nz, real64 const ( &bounds )[6] )
{
  vtkIdType const numCells = mesh.GetNumberOfCells();

  real64 const xMin = bounds[0], xMax = bounds[1];
  real64 const yMin = bounds[2], yMax = bounds[3];
  real64 const zMin = bounds[4], zMax = bounds[5];

  GEOS_ERROR_IF( nx > 1 && xMax <= xMin,
                 GEOS_FMT( "computeCellRanksCartesian: nx={} but mesh has zero extent in x ([{}, {}])", nx, xMin, xMax ) );
  GEOS_ERROR_IF( ny > 1 && yMax <= yMin,
                 GEOS_FMT( "computeCellRanksCartesian: ny={} but mesh has zero extent in y ([{}, {}])", ny, yMin, yMax ) );
  GEOS_ERROR_IF( nz > 1 && zMax <= zMin,
                 GEOS_FMT( "computeCellRanksCartesian: nz={} but mesh has zero extent in z ([{}, {}])", nz, zMin, zMax ) );

  real64 const dx = nx > 1 ? ( xMax - xMin ) / nx : 1.0;
  real64 const dy = ny > 1 ? ( yMax - yMin ) / ny : 1.0;
  real64 const dz = nz > 1 ? ( zMax - zMin ) / nz : 1.0;

  auto centroids = computeCentroids( mesh );
  stdVector< integer > ranks( numCells );

  for( vtkIdType c = 0; c < numCells; ++c )
  {
    integer const ix = std::clamp( static_cast< integer >( ( centroids[c][0] - xMin ) / dx ), 0, nx - 1 );
    integer const iy = std::clamp( static_cast< integer >( ( centroids[c][1] - yMin ) / dy ), 0, ny - 1 );
    integer const iz = std::clamp( static_cast< integer >( ( centroids[c][2] - zMin ) / dz ), 0, nz - 1 );
    // Rank ordering: ix + nx*(iy + ny*iz)  (X-fastest, Z-slowest).
    // This matches the convention GEOS uses for -x/-y/-z grid decomposition.
    ranks[c] = ix + nx * ( iy + ny * iz );
  }
  return ranks;
}


// ============================================================================
// Cell-rank assignment: Recursive Coordinate Bisection (nth_element)
// ============================================================================

stdVector< integer >
computeCellRanksRCB( vtkDataSet & mesh, integer size )
{
  auto centroids = computeCentroids( mesh );
  vtkIdType const n = mesh.GetNumberOfCells();

  stdVector< vtkIdType > indices( n );
  std::iota( indices.begin(), indices.end(), 0 );

  stdVector< integer > ranks( n );

  constexpr real64 realMax = std::numeric_limits< real64 >::max();

  std::function< void( vtkIdType, vtkIdType, integer, integer ) > bisect;
  bisect = [&]( vtkIdType begin, vtkIdType end, integer rankLo, integer rankHi )
  {
    if( rankHi - rankLo == 1 )
    {
      for( vtkIdType i = begin; i < end; ++i )
      {
        ranks[indices[i]] = rankLo;
      }
      return;
    }

    // Find the dimension with the largest spread
    real64 lo[3] = { realMax, realMax, realMax };
    real64 hi[3] = { -realMax, -realMax, -realMax };
    for( vtkIdType i = begin; i < end; ++i )
    {
      auto const & c = centroids[indices[i]];
      for( integer d = 0; d < 3; ++d )
      {
        lo[d] = LvArray::math::min( lo[d], c[d] );
        hi[d] = LvArray::math::max( hi[d], c[d] );
      }
    }
    integer bestDim = 0;
    if( hi[1] - lo[1] > hi[bestDim] - lo[bestDim] )
      bestDim = 1;
    if( hi[2] - lo[2] > hi[bestDim] - lo[bestDim] )
      bestDim = 2;

    // Split proportionally to balance cell counts between left/right rank groups
    integer const leftParts = ( rankHi - rankLo ) / 2;
    integer const totalParts = rankHi - rankLo;
    vtkIdType const splitAt = begin + ( ( end - begin ) * leftParts ) / totalParts;
    integer const rankMid = rankLo + leftParts;

    std::nth_element( indices.begin() + begin,
                      indices.begin() + splitAt,
                      indices.begin() + end,
                      [&]( vtkIdType a, vtkIdType b )
    { return centroids[a][bestDim] < centroids[b][bestDim]; } );

    bisect( begin, splitAt, rankLo, rankMid );
    bisect( splitAt, end, rankMid, rankHi );
  };

  bisect( 0, n, 0, size );
  return ranks;
}

// ============================================================================
// Distributed input: every rank assigns its own cells, then cells move once,
// directly to their destination. No rank collects the whole mesh, and every
// collective below exchanges data whose size is independent of the cell count.
// ============================================================================

/// Global index of the first local cell, with cells numbered in rank order.
vtkIdType firstGlobalCell( vtkIdType localCells, MPI_Comm comm )
{
  vtkIdType first = 0;
#ifdef GEOS_USE_MPI
  MpiWrapper::exscan( &localCells, &first, 1, MPI_SUM, comm );
  if( MpiWrapper::commRank( comm ) == 0 )
  {
    first = 0;
  }
#else
  GEOS_UNUSED_VAR( localCells, comm );
#endif
  return first;
}

/// Same blocks as computeCellRanksContiguous, applied to the rank-ordered global index.
stdVector< integer >
distributedRanksContiguous( vtkIdType localCells, vtkIdType totalCells, integer size, MPI_Comm comm )
{
  vtkIdType const first = firstGlobalCell( localCells, comm );
  vtkIdType const perRank = totalCells / size;
  vtkIdType const remainder = totalCells % size;
  vtkIdType const largeBlocks = remainder * ( perRank + 1 );
  stdVector< integer > ranks( localCells );
  for( vtkIdType i = 0; i < localCells; ++i )
  {
    vtkIdType const g = first + i;
    ranks[i] = static_cast< integer >( g < largeBlocks ? g / ( perRank + 1 ) : remainder + ( g - largeBlocks ) / perRank );
  }
  return ranks;
}

/// Same grid as the serial method: the global bounds come from one reduction.
stdVector< integer >
distributedRanksCartesian( vtkDataSet & mesh, integer nx, integer ny, integer nz, MPI_Comm comm )
{
  stdVector< real64 > low( 3, std::numeric_limits< real64 >::max() ), high( 3, -std::numeric_limits< real64 >::max() );
  if( mesh.GetNumberOfPoints() > 0 )
  {
    real64 local[6];
    mesh.GetBounds( local );
    for( integer d = 0; d < 3; ++d )
    {
      low[d] = local[2 * d];
      high[d] = local[2 * d + 1];
    }
  }
  stdVector< real64 > globalLow( 3 ), globalHigh( 3 );
  MpiWrapper::allReduce( low, globalLow, MpiWrapper::Reduction::Min, comm );
  MpiWrapper::allReduce( high, globalHigh, MpiWrapper::Reduction::Max, comm );
  real64 const bounds[6] = { globalLow[0], globalHigh[0], globalLow[1], globalHigh[1], globalLow[2], globalHigh[2] };
  return computeCellRanksCartesian( mesh, nx, ny, nz, bounds );
}

/// Order-preserving unsigned encoding of a double.
std::uint64_t orderedBits( real64 const value )
{
  std::uint64_t bits;
  std::memcpy( &bits, &value, sizeof( bits ) );
  return ( bits >> 63 ) ? ~bits : bits | ( UINT64_C( 1 ) << 63 );
}

/**
 * Distributed recursive coordinate bisection.
 *
 * The rank ranges of each level depend only on the MPI size, so all ranks hold
 * the same list of ranges. For each range, one reduction gives its cell count
 * and centroid box, which select the split dimension and the number of cells
 * that go left, with the same rule as the serial method. The split key of a
 * cell is (centroid coordinate, global cell index). The keys are unique, so the
 * split is exact even when many centroids share a coordinate. A radix selection
 * finds the split key: each round reduces a 16-bin histogram per range and
 * keeps only the local cells in the selected bin. A level needs at most 32
 * rounds, and usually stops early when the selected bin holds one key.
 */
stdVector< integer >
distributedRanksRCB( vtkDataSet & mesh, integer size, MPI_Comm comm )
{
  vtkIdType const n = mesh.GetNumberOfCells();
  auto const centroids = computeCentroids( mesh );
  vtkIdType const first = firstGlobalCell( n, comm );

  struct Range
  {
    integer lo, hi;
  };
  stdVector< Range > ranges{ { 0, size } };
  stdVector< integer > range( n, 0 );
  constexpr int digitBits = 4;
  constexpr int bins = 1 << digitBits;
  constexpr int keyBits = 128;

  while( std::any_of( ranges.begin(), ranges.end(), []( Range const & r ) { return r.hi - r.lo > 1; } ) )
  {
    std::size_t const m = ranges.size();
    stdVector< int64_t > counts( m, 0 ), globalCounts( m );
    stdVector< real64 > low( 3 * m, std::numeric_limits< real64 >::max() ), high( 3 * m, -std::numeric_limits< real64 >::max() );
    for( vtkIdType c = 0; c < n; ++c )
    {
      std::size_t const j = range[c];
      ++counts[j];
      for( integer d = 0; d < 3; ++d )
      {
        low[3 * j + d] = std::min( low[3 * j + d], centroids[c][d] );
        high[3 * j + d] = std::max( high[3 * j + d], centroids[c][d] );
      }
    }
    stdVector< real64 > globalLow( 3 * m ), globalHigh( 3 * m );
    MpiWrapper::allReduce( counts, globalCounts, MpiWrapper::Reduction::Sum, comm );
    MpiWrapper::allReduce( low, globalLow, MpiWrapper::Reduction::Min, comm );
    MpiWrapper::allReduce( high, globalHigh, MpiWrapper::Reduction::Max, comm );

    // Split dimension and left count of each range, as in the serial method.
    stdVector< integer > dims( m, 0 );
    stdVector< int64_t > remaining( m, 0 );
    stdVector< char > active( m, 0 );
    for( std::size_t j = 0; j < m; ++j )
    {
      integer const parts = ranges[j].hi - ranges[j].lo;
      if( parts < 2 || globalCounts[j] == 0 )
      {
        continue;
      }
      active[j] = 1;
      for( integer d = 1; d < 3; ++d )
      {
        if( globalHigh[3 * j + d] - globalLow[3 * j + d] > globalHigh[3 * j + dims[j]] - globalLow[3 * j + dims[j]] )
        {
          dims[j] = d;
        }
      }
      remaining[j] = globalCounts[j] * ( parts / 2 ) / parts;
    }

    // 128-bit split keys: centroid coordinate, then global cell index.
    auto keyWord = [&]( vtkIdType c, int word ) -> std::uint64_t
    {
      return word == 0 ? orderedBits( centroids[c][dims[range[c]]] ) : static_cast< std::uint64_t >( first + c );
    };
    auto digit = [&]( vtkIdType c, int round ) -> int
    {
      int const bit = keyBits - digitBits * ( round + 1 );
      return static_cast< int >( ( keyWord( c, bit < 64 ? 1 : 0 ) >> ( bit % 64 ) ) & ( bins - 1 ) );
    };

    // Radix selection of the remaining[j]-th smallest key of every range.
    stdVector< stdVector< vtkIdType > > candidates( m );
    for( vtkIdType c = 0; c < n; ++c )
    {
      if( active[range[c]] )
      {
        candidates[range[c]].push_back( c );
      }
    }
    stdVector< std::uint64_t > pivot( 2 * m, 0 );
    int rounds = 0;
    bool unique = false;
    while( rounds < keyBits / digitBits && !unique )
    {
      stdVector< int64_t > histogram( bins * m, 0 ), globalHistogram( bins * m );
      for( std::size_t j = 0; j < m; ++j )
      {
        for( vtkIdType c : candidates[j] )
        {
          ++histogram[bins * j + digit( c, rounds )];
        }
      }
      MpiWrapper::allReduce( histogram, globalHistogram, MpiWrapper::Reduction::Sum, comm );
      unique = true;
      for( std::size_t j = 0; j < m; ++j )
      {
        if( !active[j] )
        {
          continue;
        }
        int b = 0;
        while( remaining[j] >= globalHistogram[bins * j + b] )
        {
          remaining[j] -= globalHistogram[bins * j + b];
          ++b;
        }
        int const bit = keyBits - digitBits * ( rounds + 1 );
        pivot[2 * j + ( bit < 64 ? 1 : 0 )] |= static_cast< std::uint64_t >( b ) << ( bit % 64 );
        unique = unique && globalHistogram[bins * j + b] == 1;
        auto & list = candidates[j];
        list.erase( std::remove_if( list.begin(), list.end(), [&]( vtkIdType c ) { return digit( c, rounds ) != b; } ), list.end() );
      }
      ++rounds;
    }

    // Keys whose leading 4*rounds bits precede the pivot go left. The pivot's
    // prefix identifies a single key, or the whole key after 32 rounds.
    int const prefixBits = digitBits * rounds;
    auto prefix = [&]( std::uint64_t high, std::uint64_t low ) -> std::pair< std::uint64_t, std::uint64_t >
    {
      if( prefixBits <= 64 )
      {
        return { prefixBits == 0 ? 0 : high >> ( 64 - prefixBits ), 0 };
      }
      return { high, prefixBits == 128 ? low : low >> ( 128 - prefixBits ) };
    };
    stdVector< integer > next( m );
    stdVector< Range > nextRanges;
    for( std::size_t j = 0; j < m; ++j )
    {
      next[j] = static_cast< integer >( nextRanges.size() );
      integer const parts = ranges[j].hi - ranges[j].lo;
      if( parts < 2 )
      {
        nextRanges.push_back( ranges[j] );
      }
      else
      {
        integer const mid = ranges[j].lo + parts / 2;
        nextRanges.push_back( { ranges[j].lo, mid } );
        nextRanges.push_back( { mid, ranges[j].hi } );
      }
    }
    for( vtkIdType c = 0; c < n; ++c )
    {
      std::size_t const j = range[c];
      bool left = false;
      if( active[j] )
      {
        left = prefix( keyWord( c, 0 ), keyWord( c, 1 ) ) < prefix( pivot[2 * j], pivot[2 * j + 1] );
      }
      range[c] = next[j] + ( ranges[j].hi - ranges[j].lo > 1 && !left ? 1 : 0 );
    }
    ranges = std::move( nextRanges );
  }

  stdVector< integer > ranks( n );
  for( vtkIdType c = 0; c < n; ++c )
  {
    ranks[c] = ranges[range[c]].lo;
  }
  return ranks;
}

/// Split a local mesh into one piece per destination rank, halving the destination range at each step.
void splitByDestination( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                         stdVector< integer > assignment,
                         stdMap< integer, vtkSmartPointer< vtkUnstructuredGrid > > & pieces )
{
  if( assignment.empty() )
  {
    return;
  }
  auto const [low, high] = std::minmax_element( assignment.begin(), assignment.end() );
  if( *low == *high )
  {
    pieces.get_inserted( *low ) = mesh;
    return;
  }
  integer const mid = *low + ( *high - *low + 1 ) / 2;
  auto split = splitByMid( mesh, assignment, mid );
  mesh = nullptr;
  assignment.clear();
  splitByDestination( std::move( split.loMesh ), std::move( split.loAssignment ), pieces );
  splitByDestination( std::move( split.hiMesh ), std::move( split.hiAssignment ), pieces );
}

/**
 * Send each local cell directly to its assigned rank. Every rank may hold input
 * cells. Points shared by pieces from different ranks are merged by appendMeshParts.
 */
vtkSmartPointer< vtkUnstructuredGrid >
exchangeByRankAssignment( vtkUnstructuredGrid * mesh, stdVector< integer > assignment, MPI_Comm comm )
{
  integer const rank = MpiWrapper::commRank( comm );
  stdMap< integer, vtkSmartPointer< vtkUnstructuredGrid > > pieces;
  {
    auto local = vtkSmartPointer< vtkUnstructuredGrid >::New();
    local->ShallowCopy( mesh );
    splitByDestination( std::move( local ), std::move( assignment ), pieces );
  }
  stdVector< vtkSmartPointer< vtkUnstructuredGrid > > received;
  if( pieces.count( rank ) )
  {
    received.push_back( pieces.at( rank ) );
    pieces.erase( rank );
  }
  stdMap< int, stdVector< char > > outgoing;
  for( auto & [peer, piece] : pieces )
  {
    packGrid( piece, {}, outgoing.get_inserted( peer ) );
    piece = nullptr;
  }
  pieces.clear();
  auto incoming = mpi::sparseExchange( outgoing, comm );
  outgoing.clear();
  for( auto & [peer, buffer] : incoming )
  {
    GEOS_UNUSED_VAR( peer );
    received.push_back( unpackGrid( buffer ).first );
    buffer = {};
  }
  if( received.empty() )
  {
    return vtkSmartPointer< vtkUnstructuredGrid >::New();
  }
  if( received.size() == 1 )
  {
    return received.front();
  }
  stdVector< vtkUnstructuredGrid * > parts;
  for( auto const & piece : received )
  {
    parts.push_back( piece.GetPointer() );
  }
  return appendMeshParts( parts );
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

// ----------------------------------------------------------------------------
// Binary-tree scatter
//
// Only rank 0 holds mesh data and the assignment vector initially.
// At each level, the "root" of each active subrange extracts the upper
// half, packs it into a raw buffer, sends it to the midpoint rank, and
// keeps only the lower half.
// ----------------------------------------------------------------------------
vtkSmartPointer< vtkUnstructuredGrid >
scatterByRankAssignment( vtkUnstructuredGrid * inputMesh,
                         stdVector< integer > assignment,
                         MPI_Comm comm )
{
  integer const rank = MpiWrapper::commRank( comm );
  integer const size = MpiWrapper::commSize( comm );

  // Working mesh: only rank 0 starts with data
  vtkSmartPointer< vtkUnstructuredGrid > workingMesh;
  if( rank == 0 )
  {
    workingMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
    workingMesh->ShallowCopy( inputMesh );
  }

  integer lo = 0;
  integer hi = size;

  while( hi - lo > 1 )
  {
    integer const mid = lo + ( hi - lo ) / 2;

    if( rank == lo )
    {
      // Sender: always send bufSize to mid (even if mesh is empty) to avoid receiver deadlock.
      stdVector< char > buffer;
      if( workingMesh && workingMesh->GetNumberOfCells() > 0 )
      {
        // Split into lower (keep) and upper (send) halves.
        auto split = splitByMid( workingMesh, assignment, mid );
        workingMesh = nullptr;   // release original before packing

        packGrid( split.hiMesh, split.hiAssignment, buffer );
        split.hiMesh = nullptr;
        split.hiAssignment.clear();

        // Keep the lower half.
        workingMesh = std::move( split.loMesh );
        assignment = std::move( split.loAssignment );
      }

      int64_t bufSize = static_cast< int64_t >( buffer.size() );
      MpiWrapper::send( &bufSize, 1, mid, 0, comm );
      if( bufSize > 0 )
      {
        mpiSendLarge( buffer.data(), bufSize, mid, 1, comm );
      }

      hi = mid;
    }
    else if( rank == mid )
    {
      // Receiver: receive from rank lo

      int64_t bufSize = 0;
      MPI_Request request = MPI_REQUEST_NULL;
      MPI_Status status{};
      MpiWrapper::iRecv( &bufSize, 1, lo, 0, comm, &request );
      MpiWrapper::wait( &request, &status );

      if( bufSize > 0 )
      {
        stdVector< char > buffer( bufSize );
        mpiRecvLarge( buffer.data(), bufSize, lo, 1, comm );

        auto [mesh, assign] = unpackGrid( buffer );
        workingMesh = mesh;
        assignment = std::move( assign );
      }
      else
      {
        workingMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
        assignment.clear();
      }

      lo = mid;
    }
    else
    {
      // Inactive at this level, just narrow the range
      if( rank < mid )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }
  }

  if( !workingMesh )
  {
    workingMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
  }
  return workingMesh;
}

// ----------------------------------------------------------------------------
// Per-cell rank assignment dispatcher (rank 0 only)
// ----------------------------------------------------------------------------
stdVector< integer >
computeCellRanks( ScatterMethod method,
                  vtkDataSet & mesh,
                  arrayView1d< integer const > cartesianPartitions,
                  integer numRanks )
{
  GEOS_ERROR_IF( method == ScatterMethod::kdtree,
                 "computeCellRanks: kdtree method does not produce a per-cell rank assignment; "
                 "use scatterMesh() for kdtree." );

  vtkIdType const numCells = mesh.GetNumberOfCells();
  stdVector< integer > cellRanks;
  switch( method )
  {
    case ScatterMethod::contiguous:
      cellRanks = computeCellRanksContiguous( numCells, numRanks );
      break;
    case ScatterMethod::cartesian:
    {
      GEOS_ERROR_IF( cartesianPartitions.size() < 3,
                     "Cartesian method requires 3 partition values (nx, ny, nz)" );
      integer const nx = cartesianPartitions[0], ny = cartesianPartitions[1], nz = cartesianPartitions[2];
      GEOS_ERROR_IF( nx * ny * nz != numRanks,
                     GEOS_FMT( "partition grid {}x{}x{} = {} does not match MPI size {}",
                               nx, ny, nz, nx * ny * nz, numRanks ) );
      real64 bounds[6];
      mesh.GetBounds( bounds );
      cellRanks = computeCellRanksCartesian( mesh, nx, ny, nz, bounds );
      break;
    }
    case ScatterMethod::rcb:
      cellRanks = computeCellRanksRCB( mesh, numRanks );
      break;
    default:
      GEOS_ERROR( GEOS_FMT( "Unknown scatter method: {}", static_cast< integer >( method ) ) );
      break;
  }
  return cellRanks;
}


stdVector< integer >
computeCellRanksDistributed( ScatterMethod method,
                             vtkDataSet & mesh,
                             vtkIdType totalCells,
                             arrayView1d< integer const > cartesianPartitions,
                             MPI_Comm comm )
{
  integer const numRanks = MpiWrapper::commSize( comm );
  switch( method )
  {
    case ScatterMethod::contiguous:
      return distributedRanksContiguous( mesh.GetNumberOfCells(), totalCells, numRanks, comm );
    case ScatterMethod::cartesian:
    {
      GEOS_ERROR_IF( cartesianPartitions.size() < 3,
                     "Cartesian method requires 3 partition values (nx, ny, nz)" );
      integer const nx = cartesianPartitions[0], ny = cartesianPartitions[1], nz = cartesianPartitions[2];
      GEOS_ERROR_IF( nx * ny * nz != numRanks,
                     GEOS_FMT( "partition grid {}x{}x{} = {} does not match MPI size {}",
                               nx, ny, nz, nx * ny * nz, numRanks ) );
      return distributedRanksCartesian( mesh, nx, ny, nz, comm );
    }
    case ScatterMethod::rcb:
      return distributedRanksRCB( mesh, numRanks, comm );
    default:
      GEOS_ERROR( GEOS_FMT( "No distributed rank assignment for scatter method {}", static_cast< integer >( method ) ) );
  }
  return {};
}

vtkSmartPointer< vtkDataSet > scatterByBlock( vtkDataSet & mesh, MPI_Comm comm )
{
  GEOS_MARK_FUNCTION;
  int const rank = MpiWrapper::commRank( comm );
  int const size = MpiWrapper::commSize( comm );
  vtkIdType const localCells = mesh.GetNumberOfCells();
  vtkIdType const totalCells = MpiWrapper::sum( localCells, comm );
  if( size == 1 || totalCells == 0 )
  {
    auto result = vtkSmartPointer< vtkUnstructuredGrid >::New();
    result->ShallowCopy( &mesh );
    return result;
  }
#ifdef GEOS_USE_MPI
  vtkIdType firstCell = 0;
  MpiWrapper::exscan( &localCells, &firstCell, 1, MPI_SUM, comm );
  if( rank == 0 )
    firstCell = 0;
  vtkIdType const cellsPerRank = totalCells / size;
  vtkIdType const remainder = totalCells % size;
  vtkNew< vtkPartitionedDataSet > partitions;
  partitions->SetNumberOfPartitions( size );
  for( int r = 0; r < size; ++r )
  {
    vtkIdType const begin = r * cellsPerRank + std::min< vtkIdType >( r, remainder );
    vtkIdType const end = begin + cellsPerRank + ( r < remainder ? 1 : 0 );
    vtkIdType const localBegin = std::max( begin, firstCell ) - firstCell;
    vtkIdType const localEnd = std::min( end, firstCell + localCells ) - firstCell;
    vtkNew< vtkExtractCells > extractor;
    extractor->SetInputDataObject( &mesh );
    if( localEnd > localBegin )
      extractor->AddCellRange( localBegin, localEnd - 1 );
    LvArray::system::FloatingPointExceptionGuard guard;
    extractor->Update();
    partitions->SetPartition( r, extractor->GetOutput() );
  }
  auto result = vtk::redistribute( *partitions, comm );
  vtkIdType const after = MpiWrapper::sum( result->GetNumberOfCells(), comm );
  GEOS_ERROR_IF( after != totalCells,
                 GEOS_FMT( "Cell conservation failed during block fallback ({} -> {})", totalCells, after ) );
  return result;
#else
  GEOS_UNUSED_VAR( rank );
  return nullptr; // Multiple ranks are unavailable without MPI.
#endif
}

vtkSmartPointer< vtkDataSet >
scatterMesh( ScatterMethod method,
             vtkDataSet & mesh,
             arrayView1d< integer const > cartesianPartitions,
             MPI_Comm comm )
{
  GEOS_MARK_FUNCTION;

  integer const rank = MpiWrapper::commRank( comm );
  integer const size = MpiWrapper::commSize( comm );

  // Early exits
  vtkIdType const localCells = mesh.GetNumberOfCells();
  vtkIdType const totalCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Sum, comm );

  if( totalCells == 0 )
  {
    return vtkSmartPointer< vtkUnstructuredGrid >::New();
  }
  if( size == 1 )
  {
    auto copy = vtkSmartPointer< vtkUnstructuredGrid >::New();
    copy->ShallowCopy( &mesh );
    return copy;
  }

  // KdTree: legacy path using VTK's built-in redistribution
#ifdef GEOS_USE_MPI
  if( method == ScatterMethod::kdtree )
  {
    constexpr vtkIdType largeMeshCellCount = 10000000;
    GEOS_WARNING_IF( rank == 0 && totalCells > largeMeshCellCount,
                     GEOS_FMT( "Scattering {} cells with scatterMethod=kdtree. This method is significantly "
                               "slower than the alternatives on meshes of this size; "
                               "consider scatterMethod=rcb (geometry aware) instead.",
                               totalCells ) );
    vtkNew< vtkMPIController > controller;
    vtkMPICommunicatorOpaqueComm vtkComm( &comm );
    vtkNew< vtkMPICommunicator > communicator;
    communicator->InitializeExternal( &vtkComm );
    controller->SetCommunicator( communicator );
    // VTK's global controller is borrowed. Keep the previous object alive and
    // restore it before this block's local controller is destroyed.
    struct RestoreGlobalController
    {
      vtkSmartPointer< vtkMultiProcessController > previous;
      ~RestoreGlobalController() { vtkMultiProcessController::SetGlobalController( previous ); }
    } restoreController{ vtkMultiProcessController::GetGlobalController() };
    vtkMultiProcessController::SetGlobalController( controller );

    vtkNew< vtkRedistributeDataSetFilter > rdsf;
    rdsf->SetInputDataObject( &mesh );
    rdsf->SetNumberOfPartitions( size );
    rdsf->SetController( controller );
#if VTK_VERSION_NUMBER == VTK_VERSION_CHECK( 9, 7, 0 )
    // VTK 9.7 computes its cuts from dataset bounds but balances cell centers.
    // Some valid GEOS meshes have cell centers outside those bounds, which makes
    // VTK's DIY KdTree abort while building its local histogram. Extend VTK's
    // normal inflated domain to include those centers before generating cuts.
    vtkNew< vtkCellCenters > cellCenters;
    cellCenters->SetInputData( &mesh );
    cellCenters->Update();
    double cellCenterBounds[6];
    cellCenters->GetOutput()->GetBounds( cellCenterBounds );
    vtkBoundingBox localBounds;
    localBounds.AddBounds( mesh.GetBounds() );
    if( localBounds.IsValid() )
    {
      double constexpr boundingBoxLengthTolerance = 0.01;
      double constexpr boundingBoxInflationRatio = 0.01;
      double const xInflate = localBounds.GetLength( 0 ) < boundingBoxLengthTolerance
                              ? boundingBoxLengthTolerance
                              : boundingBoxInflationRatio * localBounds.GetLength( 0 );
      double const yInflate = localBounds.GetLength( 1 ) < boundingBoxLengthTolerance
                              ? boundingBoxLengthTolerance
                              : boundingBoxInflationRatio * localBounds.GetLength( 1 );
      double const zInflate = localBounds.GetLength( 2 ) < boundingBoxLengthTolerance
                              ? boundingBoxLengthTolerance
                              : boundingBoxInflationRatio * localBounds.GetLength( 2 );
      localBounds.Inflate( xInflate, yInflate, zInflate );
    }
    localBounds.AddBounds( cellCenterBounds );
    double correctedBounds[6];
    double const * correctedBoundsPtr = nullptr;
    if( localBounds.IsValid() )
    {
      localBounds.GetBounds( correctedBounds );
      correctedBoundsPtr = correctedBounds;
    }
    auto const cuts = vtkDIYKdTreeUtilities::GenerateCuts(
      &mesh, size, true, controller, correctedBoundsPtr );
    rdsf->SetUseExplicitCuts( true );
    rdsf->SetExplicitCuts( cuts );
#endif
    {
      // vtkRedistributeDataSetFilter uses VTK's XML writer internally to
      // serialize datasets exchanged by DIY. The writer may raise floating-point
      // exceptions while calculating progress for empty arrays. These exceptions
      // are harmless to VTK, but GEOS' enabled FPE traps turn them into SIGFPE.
      LvArray::system::FloatingPointExceptionGuard guard;
      rdsf->Update();
    }

    vtkSmartPointer< vtkDataSet > kdResult = vtkDataSet::SafeDownCast( rdsf->GetOutputDataObject( 0 ) );

    vtkIdType const localAfter = kdResult == nullptr ? 0 : kdResult->GetNumberOfCells();
    vtkIdType const totalAfter = MpiWrapper::allReduce( localAfter, MpiWrapper::Reduction::Sum, comm );

    if( totalAfter != totalCells )
    {
      GEOS_WARNING_IF( rank == 0,
                       GEOS_FMT( "VTK KdTree redistribution lost {} elements! Falling back to block redistribution.",
                                 totalCells - totalAfter ) );
      return scatterByBlock( mesh, comm );
    }
    return kdResult;
  }
#endif

  vtkUnstructuredGrid * inputGrid = vtkUnstructuredGrid::SafeDownCast( &mesh );
  GEOS_ERROR_IF( localCells > 0 && inputGrid == nullptr,
                 "input must be a vtkUnstructuredGrid" );
  vtkNew< vtkUnstructuredGrid > emptyGrid;
  if( inputGrid == nullptr )
  {
    inputGrid = emptyGrid.GetPointer();
  }

  vtkSmartPointer< vtkUnstructuredGrid > result;
  vtkIdType const rootCells = MpiWrapper::allReduce( rank == 0 ? localCells : vtkIdType{ 0 },
                                                     MpiWrapper::Reduction::Sum, comm );
  if( rootCells == totalCells )
  {
    // Serial input: rank 0 computes the assignment and ships cells through a binary tree.
    stdVector< integer > cellRanks;
    if( rank == 0 )
    {
      cellRanks = computeCellRanks( method, *inputGrid, cartesianPartitions, size );
    }
    result = scatterByRankAssignment( inputGrid, std::move( cellRanks ), comm );
  }
  else
  {
    // Distributed input (e.g. a .pvtu read in pieces): each rank assigns and sends its own cells.
    result = exchangeByRankAssignment( inputGrid,
                                       computeCellRanksDistributed( method, *inputGrid, totalCells, cartesianPartitions, comm ),
                                       comm );
  }

  // Validate cell conservation
  vtkIdType const localAfter = result->GetNumberOfCells();
  vtkIdType const totalAfter = MpiWrapper::allReduce( localAfter, MpiWrapper::Reduction::Sum, comm );

  GEOS_ERROR_IF( totalAfter != totalCells,
                 GEOS_FMT( "Cell conservation failed during scatter ({} -> {})", totalCells, totalAfter ) );

  return result;
}


} // namespace vtk
} // namespace geos
