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

#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/math.hpp"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkExtractCells.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkMultiProcessController.h>
#ifdef GEOS_USE_MPI
#include <vtkMPIController.h>
#include <vtkMPI.h>
#else
#include <vtkDummyController.h>
#endif
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>

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

constexpr int64_t MPI_CHUNK_SIZE = INT64_C( 1 ) << 30; // 1 GB

void mpiSendLarge( void const * buf, int64_t count, integer dest, integer tag, MPI_Comm comm )
{
  char const * ptr = static_cast< char const * >( buf );
  while( count > 0 )
  {
    integer const chunk = static_cast< integer >( LvArray::math::min( count, MPI_CHUNK_SIZE ) );
    MpiWrapper::send( ptr, chunk, dest, tag, comm );
    ptr += chunk;
    count -= chunk;
  }
}

void mpiRecvLarge( void * buf, int64_t count, integer src, integer tag, MPI_Comm comm )
{
  char * ptr = static_cast< char * >( buf );
  while( count > 0 )
  {
    integer const chunk = static_cast< integer >( LvArray::math::min( count, MPI_CHUNK_SIZE ) );
    MPI_Recv( ptr, chunk, MPI_BYTE, src, tag, comm, MPI_STATUS_IGNORE );
    ptr += chunk;
    count -= chunk;
  }
}

// ============================================================================
// Pack / unpack vtkDataSetAttributes (cell data or point data)
// ============================================================================

constexpr int NUM_ATTR_TYPES = vtkDataSetAttributes::NUM_ATTRIBUTES;

void packDataArrays( stdVector< char > & buf, vtkDataSetAttributes * attrs )
{
  // Count only vtkDataArray entries (skip string / abstract arrays)
  int32_t nArrays = 0;
  for( int a = 0; a < attrs->GetNumberOfArrays(); ++a )
  {
    if( attrs->GetArray( a ) != nullptr )
    {
      ++nArrays;
    }
  }
  appendValue( buf, nArrays );

  for( int a = 0; a < attrs->GetNumberOfArrays(); ++a )
  {
    vtkDataArray * arr = attrs->GetArray( a );
    if( arr == nullptr )
    {
      continue;
    }
    char const * name = arr->GetName() ? arr->GetName() : "";
    int32_t const nameLen = static_cast< int32_t >( std::strlen( name ) );
    appendValue( buf, nameLen );
    appendBytes( buf, name, nameLen );

    int32_t const nComp = arr->GetNumberOfComponents();
    int32_t const dataType = arr->GetDataType();
    int64_t const nTuples = arr->GetNumberOfTuples();
    int32_t const elemSize = arr->GetDataTypeSize();

    appendValue( buf, nComp );
    appendValue( buf, dataType );
    appendValue( buf, nTuples );
    appendValue( buf, elemSize );

    int64_t const dataBytes = nTuples * static_cast< int64_t >( nComp ) * static_cast< int64_t >( elemSize );
    appendBytes( buf, arr->GetVoidPointer( 0 ), dataBytes );
  }

  // Serialize active attribute designations (SCALARS=0 .. GLOBALIDS=5).
  // For each slot store the name length + name (empty string = none).
  for( int t = 0; t < NUM_ATTR_TYPES; ++t )
  {
    vtkDataArray * active = attrs->GetAttribute( t );
    char const * activeName = ( active && active->GetName() ) ? active->GetName() : "";
    int32_t const nameLen = static_cast< int32_t >( std::strlen( activeName ) );
    appendValue( buf, nameLen );
    if( nameLen > 0 )
    {
      appendBytes( buf, activeName, nameLen );
    }
  }
}

void unpackDataArrays( char const * & ptr, vtkDataSetAttributes * attrs )
{
  int32_t const nArrays = readValue< int32_t >( ptr );

  for( int a = 0; a < nArrays; ++a )
  {
    int32_t const nameLen = readValue< int32_t >( ptr );
    string name( ptr, nameLen );
    ptr += nameLen;

    int32_t const nComp = readValue< int32_t >( ptr );
    int32_t const dataType = readValue< int32_t >( ptr );
    int64_t const nTuples = readValue< int64_t >( ptr );
    int32_t const elemSize = readValue< int32_t >( ptr );

    vtkSmartPointer< vtkDataArray > arr( vtkDataArray::CreateDataArray( dataType ) );
    arr->SetName( name.c_str() );
    arr->SetNumberOfComponents( nComp );
    arr->SetNumberOfTuples( nTuples );

    int64_t const dataBytes = nTuples * static_cast< int64_t >( nComp ) * static_cast< int64_t >( elemSize );
    std::memcpy( arr->GetVoidPointer( 0 ), ptr, static_cast< size_t >( dataBytes ) );
    ptr += dataBytes;

    attrs->AddArray( arr );
  }

  // Restore active attribute designations
  for( int t = 0; t < NUM_ATTR_TYPES; ++t )
  {
    int32_t const nameLen = readValue< int32_t >( ptr );
    if( nameLen > 0 )
    {
      string activeName( ptr, nameLen );
      ptr += nameLen;
      attrs->SetActiveAttribute( activeName.c_str(), t );
    }
  }
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
    if( points->GetDataType() == VTK_DOUBLE )
    {
      appendBytes( buf, points->GetVoidPointer( 0 ), nPoints * 3 * sizeof( real64 ) );
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
    vtkUnsignedCharArray * types = grid->GetCellTypesArray();
    appendBytes( buf, types->GetVoidPointer( 0 ), nCells * sizeof( unsigned char ) );

    vtkCellArray * cells = grid->GetCells();
    vtkDataArray * offsets = cells->GetOffsetsArray();
    vtkDataArray * conn = cells->GetConnectivityArray();

    int64_t const nOffsets = offsets->GetNumberOfValues();
    int64_t const connSize = conn->GetNumberOfValues();

    appendValue( buf, connSize );
    appendBytes( buf, offsets->GetVoidPointer( 0 ), nOffsets * sizeof( vtkIdType ) );
    appendBytes( buf, conn->GetVoidPointer( 0 ), connSize * sizeof( vtkIdType ) );
  }

  // Field data arrays
  packDataArrays( buf, grid->GetCellData() );
  packDataArrays( buf, grid->GetPointData() );

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
    std::memcpy( points->GetVoidPointer( 0 ), ptr, nPoints * 3 * sizeof( real64 ) );
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

    int64_t const connSize = readValue< int64_t >( ptr );

    vtkNew< vtkIdTypeArray > offsets;
    offsets->SetNumberOfValues( nCells + 1 );
    std::memcpy( offsets->GetVoidPointer( 0 ), ptr, ( nCells + 1 ) * sizeof( vtkIdType ) );
    ptr += ( nCells + 1 ) * sizeof( vtkIdType );

    vtkNew< vtkIdTypeArray > conn;
    conn->SetNumberOfValues( connSize );
    std::memcpy( conn->GetVoidPointer( 0 ), ptr, connSize * sizeof( vtkIdType ) );
    ptr += connSize * sizeof( vtkIdType );

    vtkNew< vtkCellArray > cellArray;
    cellArray->SetData( offsets, conn );
    grid->SetCells( types, cellArray );
  }

  // Field data
  unpackDataArrays( ptr, grid->GetCellData() );
  unpackDataArrays( ptr, grid->GetPointData() );

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
computeCellRanksCartesian( vtkDataSet & mesh, integer nx, integer ny, integer nz )
{
  vtkIdType const numCells = mesh.GetNumberOfCells();

  real64 bounds[6];
  mesh.GetBounds( bounds );
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
// Binary-tree scatter
//
// Only rank 0 holds mesh data and the assignment vector initially.
// At each level, the "root" of each active subrange extracts the upper
// half, packs it into a raw buffer, sends it to the midpoint rank, and
// keeps only the lower half.
// ============================================================================

vtkSmartPointer< vtkUnstructuredGrid >
binaryTreeScatter( vtkUnstructuredGrid * inputMesh,
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
      MPI_Recv( &bufSize, 1, MPI_INT64_T, lo, 0, comm, MPI_STATUS_IGNORE );

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

} // anonymous namespace


// ============================================================================
// Public API
// ============================================================================

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

  GEOS_ERROR_IF( rank == 0 && localCells != totalCells,
                 GEOS_FMT( "Rank 0 must have the complete mesh. "
                           "Rank 0 has {} cells but total is {}", localCells, totalCells ) );

  // KdTree: legacy path using VTK's built-in redistribution
  if( method == ScatterMethod::kdtree )
  {
#ifdef GEOS_USE_MPI
    vtkNew< vtkMPIController > controller;
    vtkMPICommunicatorOpaqueComm vtkComm( &comm );
    vtkNew< vtkMPICommunicator > communicator;
    communicator->InitializeExternal( &vtkComm );
    controller->SetCommunicator( communicator );
#else
    vtkNew< vtkDummyController > controller;
#endif
    vtkMultiProcessController::SetGlobalController( controller );

    vtkNew< vtkRedistributeDataSetFilter > rdsf;
    rdsf->SetInputDataObject( &mesh );
    rdsf->SetNumberOfPartitions( size );
    rdsf->Update();

    vtkSmartPointer< vtkDataSet > kdResult = vtkDataSet::SafeDownCast( rdsf->GetOutputDataObject( 0 ) );

    vtkIdType const localAfter = kdResult->GetNumberOfCells();
    vtkIdType const totalAfter = MpiWrapper::allReduce( localAfter, MpiWrapper::Reduction::Sum, comm );

    if( totalAfter != totalCells )
    {
      GEOS_WARNING_IF( rank == 0,
                       GEOS_FMT( "{} cells lost during kdtree scatter ({} -> {}). "
                                 "Falling back to contiguous scatter.",
                                 totalCells - totalAfter, totalCells, totalAfter ) );
      // Fall through to the contiguous path below.
      method = ScatterMethod::contiguous;
    }
    else
    {
      return kdResult;
    }
  }

  // Compute cell to rank assignment (rank 0 only)
  stdVector< integer > cellRanks;

  if( rank == 0 )
  {
    switch( method )
    {
      case ScatterMethod::contiguous:
        cellRanks = computeCellRanksContiguous( totalCells, size );
        break;
      case ScatterMethod::cartesian:
      {
        GEOS_ERROR_IF( cartesianPartitions.size() < 3,
                       "Cartesian method requires 3 partition values (nx, ny, nz)" );
        integer const nx = cartesianPartitions[0], ny = cartesianPartitions[1], nz = cartesianPartitions[2];
        GEOS_ERROR_IF( nx * ny * nz != size,
                       GEOS_FMT( "partition grid {}x{}x{} = {} does not match MPI size {}",
                                 nx, ny, nz, nx * ny * nz, size ) );
        cellRanks = computeCellRanksCartesian( mesh, nx, ny, nz );
        break;
      }
      case ScatterMethod::rcb:
        cellRanks = computeCellRanksRCB( mesh, size );
        break;
      default:
        GEOS_ERROR( GEOS_FMT( "Unknown scatter method: {}", static_cast< integer >( method ) ) );
        break;
    }
  }

  // Scatter via binary tree
  auto * inputGrid = vtkUnstructuredGrid::SafeDownCast( &mesh );
  GEOS_ERROR_IF( rank == 0 && inputGrid == nullptr,
                 "input must be a vtkUnstructuredGrid" );

  // Non-rank-0 passes a dummy; binaryTreeScatter only uses rank 0's mesh.
  vtkNew< vtkUnstructuredGrid > dummyGrid;
  vtkUnstructuredGrid * gridPtr = ( rank == 0 ) ? inputGrid : dummyGrid.GetPointer();

  vtkSmartPointer< vtkUnstructuredGrid > result = binaryTreeScatter( gridPtr, std::move( cellRanks ), comm );

  // Validate cell conservation
  vtkIdType const localAfter = result->GetNumberOfCells();
  vtkIdType const totalAfter = MpiWrapper::allReduce( localAfter, MpiWrapper::Reduction::Sum, comm );

  GEOS_WARNING_IF( rank == 0 && totalAfter != totalCells,
                   GEOS_FMT( "{} cells lost during scatter ({} -> {})",
                             totalCells - totalAfter, totalCells, totalAfter ) );

  return result;
}


} // namespace vtk
} // namespace geos
