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

/**
 * @file testVTKRefinementCommunication.cpp
 */

#include "VTKRefinementTestMeshes.hpp"
#include "common/MpiChunkedCommunication.hpp"
#include "mesh/generators/VTKRefinementCommunication.hpp"
#include "mesh/generators/VTKRefinementFields.hpp"
#include "mesh/generators/VTKRefinementTemplates.hpp"
#include "mesh/generators/VTKRefinementSharing.hpp"

#include <gtest/gtest.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkTypeInt64Array.h>

#include <algorithm>
#include <climits>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>

using namespace geos;
using namespace geos::vtk::refinement;

namespace
{
EntityKey const vertex{ 0, EntityKind::vertex, { 71 } };
EntityKey const edge{ 0, EntityKind::edge, { 71, 81 } };
EntityKey const face{ 0, EntityKind::face, { 71, 81, 91, 101 } };
Participants allRanks( Communication const & comm )
{
  Participants ranks( comm.size() );
  std::iota( ranks.begin(), ranks.end(), 0 );
  return ranks;
}

struct LocalMesh
{
  std::vector< Coordinates > coordinates;
  Connectivity pointIds;
  std::vector< Cell > cells;
  Connectivity cellIds;
  std::vector< EntitySupport > supports;
};

LocalMesh cube( int x )
{
  LocalMesh mesh;
  mesh.coordinates = { { double( x ), 0, 0 }, { double( x + 1 ), 0, 0 }, { double( x + 1 ), 1, 0 }, { double( x ), 1, 0 },
                       { double( x ), 0, 1 }, { double( x + 1 ), 0, 1 }, { double( x + 1 ), 1, 1 }, { double( x ), 1, 1 } };
  for( auto const & point : mesh.coordinates )
  {
    mesh.pointIds.push_back( 1001 + 7 * ( 4 * static_cast< vtkIdType >( point[0] ) + 2 * static_cast< vtkIdType >( point[1] ) +
                                          static_cast< vtkIdType >( point[2] ) ) );
  }
  mesh.cells = { { VTK_HEXAHEDRON, { 0, 1, 2, 3, 4, 5, 6, 7 }, 0 } };
  mesh.cellIds = { 90001 + 13 * x };
  return mesh;
}

LocalMesh selectCells( LocalMesh const & complete, std::vector< int > const & selected )
{
  LocalMesh mesh;
  std::map< vtkIdType, vtkIdType > local;
  for( int index : selected )
  {
    Cell cell = complete.cells[index];
    for( vtkIdType & point : cell.points )
    {
      auto const insertion = local.emplace( point, mesh.coordinates.size() );
      if( insertion.second )
      {
        mesh.coordinates.push_back( complete.coordinates[point] );
        mesh.pointIds.push_back( complete.pointIds[point] );
      }
      point = insertion.first->second;
    }
    mesh.cells.push_back( std::move( cell ) );
    mesh.cellIds.push_back( complete.cellIds[index] );
  }
  return mesh;
}

void compareDescendants( LocalMesh const & mesh, LocalMesh const & reference )
{
  ASSERT_EQ( mesh.cells.size(), reference.cells.size() );
  for( std::size_t i = 0; i < mesh.cells.size(); ++i )
  {
    EXPECT_EQ( mesh.cells[i].vtkType, reference.cells[i].vtkType );
    ASSERT_EQ( mesh.cells[i].points.size(), reference.cells[i].points.size() );
    for( std::size_t j = 0; j < mesh.cells[i].points.size(); ++j )
    {
      for( int d = 0; d < 3; ++d )
      {
        EXPECT_NEAR( mesh.coordinates[mesh.cells[i].points[j]][d], reference.coordinates[reference.cells[i].points[j]][d], 1e-12 );
      }
    }
  }
}

std::vector< EntitySupport > incidence( LocalMesh const & mesh, int rank )
{
  return meshEntitySupports( mesh.cells, mesh.cellIds, mesh.pointIds, rank );
}

void refineDistributed( LocalMesh & mesh, Communication & comm, int generation )
{
  std::unique_ptr< PointRegistry > points;
  std::unique_ptr< SharingInheritance > inherited;
  std::unique_ptr< InterfaceSharing > interfaces;
  LocalMesh next;
  std::vector< PointCreation > creations;
  comm.checked( "connected mesh planning",
                [&]
                {
                  points = std::make_unique< PointRegistry >( mesh.coordinates, mesh.pointIds );
                  auto const boundary = coarseBoundary( mesh.cells, mesh.cellIds, mesh.pointIds, comm.rank() );
                  std::set< EntityKey > candidateKeys;
                  for( auto const & entity : boundary.entities )
                  {
                    candidateKeys.insert( entity.key );
                  }
                  for( auto const & entity : mesh.supports )
                  {
                    if( entity.participants.size() > 1 && entity.key.kind != EntityKind::cell && !candidateKeys.count( entity.key ) )
                    {
                      throw std::logic_error( "Boundary candidates omit a shared reference entity" );
                    }
                  }
                  for( std::size_t i = 0; i < mesh.cells.size(); ++i )
                  {
                    auto const split = subdivideCell( mesh.cells[i], mesh.cellIds[i], *points );
                    next.cells.insert( next.cells.end(), split.children.begin(), split.children.end() );
                  }
                  inherited = std::make_unique< SharingInheritance >( *points, mesh.supports, comm.rank() );
                  std::vector< EntitySupport > shared;
                  for( auto const & entity : mesh.supports )
                  {
                    if( entity.participants.size() > 1 && entity.key.kind != EntityKind::cell )
                    {
                      shared.push_back( entity );
                    }
                  }
                  interfaces = std::make_unique< InterfaceSharing >( *points, std::move( shared ), comm.rank() );
                  for( vtkIdType i = points->originalSize(); i < static_cast< vtkIdType >( points->points().size() ); ++i )
                  {
                    auto const & point = points->points()[i];
                    if( interfaces->participants( { i } ) != inherited->participants( { i } ) )
                    {
                      throw std::logic_error( "Interface-only point sharing differs from complete incidence" );
                    }
                    creations.push_back( { point.key, inherited->participants( { i } ), point.position, {} } );
                  }
                } );
  vtkIdType const maximum = mesh.pointIds.empty() ? -1 : *std::max_element( mesh.pointIds.begin(), mesh.pointIds.end() );
  auto const records = comm.resolvePoints( generation, creations, maximum );
  auto const range = comm.allocateRange( next.cells.size(), 0 );
  comm.checked( "connected mesh construction",
                [&]
                {
                  next.pointIds = mesh.pointIds;
                  for( vtkIdType i = 0; i < static_cast< vtkIdType >( points->points().size() ); ++i )
                  {
                    auto const & point = points->points()[i];
                    if( i < points->originalSize() )
                    {
                      next.coordinates.push_back( point.position );
                    }
                    else
                    {
                      auto const & record = records.at( point.key );
                      next.coordinates.push_back( record.position );
                      next.pointIds.push_back( record.globalId );
                    }
                  }
                  for( std::size_t i = 0; i < next.cells.size(); ++i )
                  {
                    next.cellIds.push_back( range.first + i );
                  }
                  next.supports = incidence( next, comm.rank() );
                  for( auto & entity : next.supports )
                  {
                    if( entity.key.kind != EntityKind::cell )
                    {
                      entity.participants = inherited->participants( entity.localCorners );
                      if( interfaces->participants( entity.localCorners ) != entity.participants )
                      {
                        throw std::logic_error( "Interface-only entity sharing differs from complete incidence" );
                      }
                    }
                  }
                  auto const sparse = interfaces->fineSupports( next.pointIds );
                  std::map< EntityKey, Participants > expected;
                  for( auto const & entity : next.supports )
                  {
                    if( entity.participants.size() > 1 )
                    {
                      expected.emplace( entity.key, entity.participants );
                    }
                  }
                  if( sparse.size() != expected.size() )
                  {
                    throw std::logic_error( "Interface-only fine incidence coverage mismatch" );
                  }
                  for( auto const & entity : sparse )
                  {
                    if( expected.at( entity.key ) != entity.participants )
                    {
                      throw std::logic_error( "Interface-only fine incidence participants mismatch" );
                    }
                  }
                } );
  mesh = std::move( next );
}

using Geometry = std::pair< int, std::vector< Coordinates > >;
std::set< Geometry > geometry( LocalMesh const & mesh )
{
  std::set< Geometry > result;
  for( auto const & cell : mesh.cells )
  {
    std::vector< Coordinates > corners;
    for( vtkIdType point : cell.points )
    {
      corners.push_back( mesh.coordinates[point] );
    }
    std::sort( corners.begin(), corners.end() );
    result.emplace( cell.vtkType, std::move( corners ) );
  }
  return result;
}

void refineSerialReference( LocalMesh & mesh )
{
  PointRegistry points( mesh.coordinates, mesh.pointIds );
  LocalMesh next;
  for( std::size_t i = 0; i < mesh.cells.size(); ++i )
  {
    auto const split = subdivideCell( mesh.cells[i], mesh.cellIds[i], points );
    next.cells.insert( next.cells.end(), split.children.begin(), split.children.end() );
  }
  next.pointIds = mesh.pointIds;
  for( vtkIdType i = 0; i < static_cast< vtkIdType >( points.points().size() ); ++i )
  {
    next.coordinates.push_back( points.position( i ) );
    if( i >= points.originalSize() )
    {
      next.pointIds.push_back( 1000000 + i );
    }
  }
  next.cellIds.resize( next.cells.size() );
  std::iota( next.cellIds.begin(), next.cellIds.end(), 0 );
  mesh = std::move( next );
}
} // namespace

TEST( VTKRefinementCommunication, BoundaryCandidatesCoverAllSharedCoarseEntities )
{
  Communication comm( MPI_COMM_GEOS );
  LocalMesh mesh;
  int const bx = 2 * ( comm.rank() % 2 ), by = 2 * ( ( comm.rank() / 2 ) % 2 ), bz = 2 * ( comm.rank() / 4 );
  auto point = []( int x, int y, int z ) { return x + 3 * ( y + 3 * z ); };
  for( int z = 0; z <= 2; ++z )
  {
    for( int y = 0; y <= 2; ++y )
    {
      for( int x = 0; x <= 2; ++x )
      {
        mesh.coordinates.push_back( { double( bx + x ), double( by + y ), double( bz + z ) } );
        vtkIdType const base = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741001 ) ) : 1001;
        mesh.pointIds.push_back( base + 3 * ( bx + x + 10 * ( by + y + 10 * ( bz + z ) ) ) );
      }
    }
  }
  for( int z = 0; z < 2; ++z )
  {
    for( int y = 0; y < 2; ++y )
    {
      for( int x = 0; x < 2; ++x )
      {
        mesh.cells.push_back( { VTK_HEXAHEDRON,
                                { point( x, y, z ), point( x + 1, y, z ), point( x + 1, y + 1, z ), point( x, y + 1, z ),
                                  point( x, y, z + 1 ), point( x + 1, y, z + 1 ), point( x + 1, y + 1, z + 1 ), point( x, y + 1, z + 1 ) },
                                0 } );
        mesh.cellIds.push_back( 11 + 13 * ( bx + x + 10 * ( by + y + 10 * ( bz + z ) ) ) );
      }
    }
  }
  mesh.supports = incidence( mesh, comm.rank() );
  std::vector< EntityKey > all;
  for( auto const & entity : mesh.supports )
  {
    all.push_back( entity.key );
  }
  auto const reference = comm.discoverSharing( all );
  auto const boundary = coarseBoundary( mesh.cells, mesh.cellIds, mesh.pointIds, comm.rank() );
  EXPECT_EQ( boundary.entities.size(), 98 );
  EXPECT_EQ( boundary.volumeIds.size(), 8 );
  auto keys = boundary.volumeIds;
  for( auto const & entity : boundary.entities )
  {
    keys.push_back( entity.key );
  }
  auto const candidates = comm.discoverSharing( keys );
  for( auto const & [key, ranks] : reference )
  {
    if( ranks.size() > 1 )
    {
      ASSERT_TRUE( candidates.count( key ) );
      EXPECT_EQ( candidates.at( key ), ranks );
    }
  }
  for( auto & entity : mesh.supports )
  {
    entity.participants = reference.at( entity.key );
  }
  for( int generation = 1; generation <= 2; ++generation )
  {
    EXPECT_NO_THROW( refineDistributed( mesh, comm, generation ) );
  }
}

TEST( VTKRefinementCommunication, ConnectedHexesInheritSharingAtTwoLevels )
{
  Communication comm( MPI_COMM_GEOS, 97 );
  LocalMesh mesh = cube( comm.rank() );
  LocalMesh reference = cube( comm.rank() );
  mesh.supports = incidence( mesh, comm.rank() );
  std::vector< EntityKey > keys;
  for( auto const & entity : mesh.supports )
  {
    keys.push_back( entity.key );
  }
  auto const sharing = comm.discoverSharing( keys );
  for( auto & entity : mesh.supports )
  {
    entity.participants = sharing.at( entity.key );
  }
  Connectivity const originalIds = mesh.pointIds;
  for( int generation = 1; generation <= 2; ++generation )
  {
    refineDistributed( mesh, comm, generation );
    refineSerialReference( reference );
    EXPECT_EQ( geometry( mesh ), geometry( reference ) );
    EXPECT_EQ( mesh.cells.size(), generation == 1 ? 8 : 64 );
    EXPECT_TRUE( std::equal( originalIds.begin(), originalIds.end(), mesh.pointIds.begin() ) );
    std::uint64_t ownedPoints = 0;
    for( auto const & entity : mesh.supports )
    {
      if( entity.key.kind == EntityKind::vertex && entity.participants.front() == comm.rank() )
      {
        ++ownedPoints;
      }
    }
    auto const unique = MpiWrapper::sum( ownedPoints, MPI_COMM_GEOS );
    EXPECT_EQ( unique, static_cast< std::uint64_t >( generation == 1 ? 9 * ( 2 * comm.size() + 1 ) : 25 * ( 4 * comm.size() + 1 ) ) );
    auto const cells = MpiWrapper::sum( static_cast< std::uint64_t >( mesh.cells.size() ), MPI_COMM_GEOS );
    EXPECT_EQ( cells, ( generation == 1 ? 8 : 64 ) * static_cast< std::uint64_t >( comm.size() ) );
    EXPECT_EQ( mesh.cellIds.front(), comm.rank() * static_cast< vtkIdType >( mesh.cells.size() ) );
    EXPECT_EQ( comm.statistics().directoryExchanges, 2 );
    EXPECT_EQ( comm.statistics().neighborExchanges, generation );
  }
}

TEST( VTKRefinementCommunication, RecursiveRefinementWithAnEmptyVolumeRank )
{
  Communication comm( MPI_COMM_GEOS, 101 );
  LocalMesh mesh;
  int const nonempty = std::max( 1, comm.size() - 1 );
  if( comm.rank() < nonempty )
  {
    mesh = cube( comm.rank() );
  }
  mesh.supports = incidence( mesh, comm.rank() );
  std::vector< EntityKey > keys;
  for( auto const & entity : mesh.supports )
  {
    keys.push_back( entity.key );
  }
  auto const sharing = comm.discoverSharing( keys );
  for( auto & entity : mesh.supports )
  {
    entity.participants = sharing.at( entity.key );
  }
  for( int generation = 1; generation <= 2; ++generation )
  {
    EXPECT_NO_THROW( refineDistributed( mesh, comm, generation ) );
    auto const count = MpiWrapper::sum( static_cast< std::uint64_t >( mesh.cells.size() ), MPI_COMM_GEOS );
    EXPECT_EQ( count, ( generation == 1 ? 8 : 64 ) * static_cast< std::uint64_t >( nonempty ) );
    if( comm.rank() >= nonempty )
    {
      EXPECT_TRUE( mesh.coordinates.empty() );
      EXPECT_TRUE( mesh.cells.empty() );
      EXPECT_TRUE( comm.neighbors().empty() );
    }
  }
}

TEST( VTKRefinementCommunication, AllMixedInterfacesConformAcrossRanksAtTwoLevels )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  int pairOrdinal = 0;
  for( int arity : { 3, 4 } )
  {
    using namespace testMeshes;
    std::vector< ReferenceCell > shapes =
        arity == 3
            ? std::vector< ReferenceCell >{ referenceCell( VTK_TETRA, 0 ), referenceCell( VTK_WEDGE, 0 ), referenceCell( VTK_PYRAMID, 1 ) }
            : std::vector< ReferenceCell >{ referenceCell( VTK_HEXAHEDRON, 0 ), referenceCell( VTK_WEDGE, 2 ),
                                            referenceCell( VTK_PYRAMID, 0 ) };
    if( arity == 4 )
    {
      for( int n = 5; n <= 11; ++n )
      {
        ReferenceCell prism;
        prism.face = 2;
        prism.cell = regularPrism( n, prism.xyz );
        shapes.push_back( std::move( prism ) );
      }
    }
    for( std::size_t a = 0; a < shapes.size(); ++a )
    {
      for( std::size_t b = a; b < shapes.size(); ++b, ++pairOrdinal )
      {
        SCOPED_TRACE( std::to_string( arity ) + ":" + std::to_string( a ) + "," + std::to_string( b ) );
        LocalMesh complete;
        complete.coordinates = arity == 3 ? std::vector< Coordinates >{ { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } }
                                          : std::vector< Coordinates >{ { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
        complete.cells.push_back( attach( shapes[a], 1, complete.coordinates ) );
        complete.cells.push_back( attach( shapes[b], -1, complete.coordinates ) );
        complete.cellIds = { 90001 + 26 * pairOrdinal, 90014 + 26 * pairOrdinal };
        for( std::size_t i = 0; i < complete.coordinates.size(); ++i )
        {
          complete.pointIds.push_back( 1001 + 7 * ( 1000 * pairOrdinal + i ) );
        }
        std::vector< int > selected;
        for( int side = 0; side < 2; ++side )
        {
          if( rank == ( pairOrdinal + side ) % size )
          {
            selected.push_back( side );
          }
        }
        LocalMesh mesh = selectCells( complete, selected );
        LocalMesh reference = mesh;
        Communication comm( MPI_COMM_GEOS, 257 );
        mesh.supports = incidence( mesh, rank );
        std::vector< EntityKey > keys;
        for( auto const & entity : mesh.supports )
        {
          keys.push_back( entity.key );
        }
        auto const sharing = comm.discoverSharing( keys );
        for( auto & entity : mesh.supports )
        {
          entity.participants = sharing.at( entity.key );
        }
        for( int generation = 1; generation <= 2; ++generation )
        {
          ASSERT_NO_THROW( refineDistributed( mesh, comm, generation ) );
          refineSerialReference( reference );
          compareDescendants( mesh, reference );
          int sharedFaces = 0;
          for( auto const & entity : mesh.supports )
          {
            if( entity.key.kind == EntityKind::face && entity.participants.size() > 1 )
            {
              ++sharedFaces;
            }
          }
          int const expected = size > 1 && !mesh.cells.empty() ? ( generation == 1 ? 4 : 16 ) : 0;
          EXPECT_EQ( sharedFaces, expected );
          auto descendants = [&]( Cell const & cell ) -> std::uint64_t
          {
            if( cell.prismSides )
            {
              return 2 * cell.prismSides * ( generation == 1 ? 1 : 8 );
            }
            return cell.vtkType == VTK_PYRAMID ? ( generation == 1 ? 10 : 92 ) : ( generation == 1 ? 8 : 64 );
          };
          EXPECT_EQ( MpiWrapper::sum( static_cast< std::uint64_t >( mesh.cells.size() ), MPI_COMM_GEOS ),
                     descendants( complete.cells[0] ) + descendants( complete.cells[1] ) );
          EXPECT_EQ( comm.statistics().directoryExchanges, 2 );
        }
      }
    }
  }
  EXPECT_EQ( pairOrdinal, 61 );
}

TEST( VTKRefinementCommunication, PolygonCapsShareIdsWithPermutedStartsAtTwoLevels )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  for( int n = 5; n <= 11; ++n )
  {
    SCOPED_TRACE( n );
    LocalMesh complete;
    Cell const lower = testMeshes::regularPrism( n, complete.coordinates );
    for( int i = 0; i < n; ++i )
    {
      auto point = complete.coordinates[n + i];
      point[0] += .3;
      point[1] -= .15;
      point[2] += 1;
      complete.coordinates.push_back( point );
    }
    Cell upper{ VTK_POLYHEDRON, {}, n };
    for( int layer = 0; layer < 2; ++layer )
    {
      for( int i = 0; i < n; ++i )
      {
        upper.points.push_back( ( layer + 1 ) * n + ( i + 3 ) % n );
      }
    }
    complete.cells = { lower, upper };
    complete.cellIds = { 9001, 9011 };
    vtkIdType const base = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741019 ) ) : 1001;
    for( std::size_t i = 0; i < complete.coordinates.size(); ++i )
    {
      complete.pointIds.push_back( base + 7 * i );
    }
    std::vector< int > selected;
    for( int side = 0; side < 2; ++side )
    {
      if( rank == ( n + side ) % size )
      {
        selected.push_back( side );
      }
    }
    LocalMesh mesh = selectCells( complete, selected ), reference = mesh;
    Communication comm( MPI_COMM_GEOS, 263 );
    mesh.supports = incidence( mesh, rank );
    std::vector< EntityKey > keys;
    for( auto const & entity : mesh.supports )
    {
      keys.push_back( entity.key );
    }
    auto const sharing = comm.discoverSharing( keys );
    for( auto & entity : mesh.supports )
    {
      entity.participants = sharing.at( entity.key );
    }
    for( int generation = 1; generation <= 2; ++generation )
    {
      ASSERT_NO_THROW( refineDistributed( mesh, comm, generation ) );
      refineSerialReference( reference );
      compareDescendants( mesh, reference );
      int sharedFaces = 0;
      std::uint64_t ownedPoints = 0;
      for( auto const & entity : mesh.supports )
      {
        if( entity.key.kind == EntityKind::face && entity.participants.size() > 1 )
        {
          ++sharedFaces;
        }
        if( entity.key.kind == EntityKind::vertex && entity.participants.front() == rank )
        {
          ++ownedPoints;
        }
      }
      EXPECT_EQ( sharedFaces, size > 1 && !mesh.cells.empty() ? n * ( generation == 1 ? 1 : 4 ) : 0 );
      // First cap: 2N+1 vertices, 3N edges, N quads. Its next pass adds
      // one point per edge and quad, giving 6N+1 vertices on each of nine planes.
      EXPECT_EQ( MpiWrapper::sum( ownedPoints, MPI_COMM_GEOS ), static_cast< std::uint64_t >( generation == 1 ? 10 * n + 5 : 54 * n + 9 ) );
      EXPECT_EQ( MpiWrapper::sum( static_cast< std::uint64_t >( mesh.cells.size() ), MPI_COMM_GEOS ),
                 static_cast< std::uint64_t >( generation == 1 ? 4 * n : 32 * n ) );
      EXPECT_EQ( comm.statistics().directoryExchanges, 2 );
    }
  }
}

TEST( VTKRefinementCommunication, DuplicateVolumeOwnersAreRejectedWithoutCellCenters )
{
  Communication comm( MPI_COMM_GEOS );
  if( comm.size() == 1 )
  {
    GTEST_SKIP() << "Requires duplicate ownership across ranks";
  }
  LocalMesh mesh;
  mesh.coordinates = testMeshes::referenceCell( VTK_TETRA, 0 ).xyz;
  for( auto & point : mesh.coordinates )
  {
    point[0] += 10 * comm.rank();
  }
  mesh.pointIds = { 71 + 4 * comm.rank(), 72 + 4 * comm.rank(), 73 + 4 * comm.rank(), 74 + 4 * comm.rank() };
  mesh.cells = { { VTK_TETRA, { 0, 1, 2, 3 }, 0 } };
  mesh.cellIds = { 8112 }; // Contradictory owners; tetrahedra create no cell-center point.
  mesh.supports = incidence( mesh, comm.rank() );
  std::vector< EntityKey > keys;
  for( auto const & entity : mesh.supports )
  {
    keys.push_back( entity.key );
  }
  auto const sharing = comm.discoverSharing( keys );
  for( auto & entity : mesh.supports )
  {
    entity.participants = sharing.at( entity.key );
  }
  EXPECT_THROW( comm.checked( "coarse volume ownership",
                              [&]
                              {
                                PointRegistry points( mesh.coordinates, mesh.pointIds );
                                SharingInheritance inherited( points, mesh.supports, comm.rank() );
                              } ),
                std::runtime_error );
}

TEST( VTKRefinementCommunication, ExistingVerticesKeepIdsAndValidateCoordinates )
{
  Communication comm( MPI_COMM_GEOS, 23 );
  vtkIdType const large = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741019 ) ) : 1001;
  EntityKey const mainVertex = entityKey( EntityKind::vertex, { large } );
  EntityKey const auxiliaryVertex = entityKey( EntityKind::vertex, { large }, 1 );
  auto const sharing = comm.discoverSharing( { mainVertex, auxiliaryVertex } );
  std::vector< PointCreation > originals{
      { mainVertex, sharing.at( mainVertex ), { 0, 0, 0 }, { static_cast< unsigned char >( comm.rank() ) } },
      { auxiliaryVertex, sharing.at( auxiliaryVertex ), { 1, 0, 0 }, { static_cast< unsigned char >( comm.rank() + 9 ) } } };
  auto const records = comm.reconcileExistingPoints( originals );
  EXPECT_EQ( records.at( mainVertex ).globalId, large );
  EXPECT_EQ( records.at( auxiliaryVertex ).globalId, large );
  EXPECT_EQ( records.at( mainVertex ).fields, ( Bytes{ 0 } ) );
  EXPECT_EQ( records.at( auxiliaryVertex ).fields, ( Bytes{ 9 } ) );
  EXPECT_EQ( comm.statistics().directoryExchanges, 2 );
  if( comm.size() > 1 )
  {
    if( comm.rank() == comm.size() - 1 )
    {
      originals[0].position[0] = 1e-3;
    }
    EXPECT_THROW( comm.reconcileExistingPoints( originals ), std::runtime_error );
    originals[0].position[0] = 0;
  }
  if( comm.rank() == comm.size() - 1 )
  {
    originals.push_back( originals.front() );
  }
  EXPECT_THROW( comm.reconcileExistingPoints( originals ), std::runtime_error );
  if( comm.rank() == comm.size() - 1 )
  {
    originals.pop_back();
  }
  EXPECT_NO_THROW( comm.reconcileExistingPoints( originals ) );
}

TEST( VTKRefinementCommunication, SharedTypedFieldsComeFromTheAllocator )
{
  Communication comm( MPI_COMM_GEOS, 19 );
  auto const sharing = comm.discoverSharing( { edge } );
  PointRegistry points( { { 0, 0, 0 }, { 1, 0, 0 } }, { 71, 81 } );
  vtkIdType const midpoint = points.edge( 0, 1 );
  vtkNew< vtkPointData > input;
  vtkNew< vtkDoubleArray > vector;
  vector->SetName( "velocity" );
  vector->SetNumberOfComponents( 3 );
  vector->SetNumberOfTuples( 2 );
  vector->SetComponentName( 0, "x" );
  vtkNew< vtkTypeInt64Array > label;
  label->SetName( "label" );
  label->SetNumberOfTuples( 2 );
  vtkTypeInt64 const large = INT64_C( 9007199254741019 );
  for( vtkIdType i = 0; i < 2; ++i )
  {
    label->SetValue( i, large + comm.rank() );
    for( int c = 0; c < 3; ++c )
    {
      vector->SetTypedComponent( i, c, 2 * i + c + comm.rank() );
    }
  }
  input->SetVectors( vector );
  input->AddArray( label );
  vtkSmartPointer< vtkPointData > transferred;
  std::unique_ptr< PointFieldLayout > fields;
  std::vector< PointCreation > creations;
  comm.checked( "typed point-field planning",
                [&]
                {
                  transferred = transferPointData( *input, points, {} );
                  fields = std::make_unique< PointFieldLayout >( *transferred );
                  creations = { { edge, sharing.at( edge ), points.position( midpoint ), fields->pack( midpoint ) } };
                } );
  auto const result = comm.resolvePoints( 1, creations, 81 );
  comm.checked( "typed point-field installation", [&] { fields->install( midpoint, result.at( edge ).fields ); } );
  EXPECT_EQ( vtkTypeInt64Array::SafeDownCast( transferred->GetArray( "label" ) )->GetValue( midpoint ), large );
  for( int c = 0; c < 3; ++c )
  {
    EXPECT_DOUBLE_EQ( transferred->GetVectors()->GetComponent( midpoint, c ), 1 + c );
  }
  if( comm.size() > 1 )
  {
    // A schema mismatch on only one rank fails at a collective boundary.
    if( comm.rank() == comm.size() - 1 )
    {
      transferred->GetVectors()->SetComponentName( 0, "different meaning" );
    }
    PointFieldLayout altered( *transferred );
    EXPECT_THROW( comm.checked( "one-rank field schema", [&] { altered.install( midpoint, result.at( edge ).fields ); } ),
                  std::runtime_error );
  }
  if( comm.rank() == comm.size() - 1 )
  {
    label->SetValue( 0, large + comm.rank() + 1 );
  }
  EXPECT_THROW( comm.checked( "one-rank categorical support", [&] { transferPointData( *input, points, {} ); } ), std::runtime_error );
}

TEST( VTKRefinementCommunication, DirectoryAndPointIdsUseFullKeys )
{
  Communication comm( MPI_COMM_GEOS,
                      11 ); // Force many small chunks without large allocations.
  auto const ranks = allRanks( comm );
  EntityKey const contactSide{ 1, EntityKind::edge, { 71, 81 } };
  auto const sharing = comm.discoverSharing( { vertex, edge, face, contactSide } );
  ASSERT_EQ( sharing.size(), 4 );
  EXPECT_EQ( sharing.at( edge ), ranks );
  EXPECT_EQ( sharing.at( contactSide ), ranks );
  vtkIdType const oldMaximum = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741091 ) ) : 1000;
  std::vector< PointCreation > points{
      { edge, ranks, { .5, 0, 0 }, { static_cast< unsigned char >( comm.rank() ) } },
      { face, ranks, { .5, .5, 0 }, {} },
      { contactSide, ranks, { .5, 0, 0 }, {} },
      { { 0, EntityKind::cell, { static_cast< vtkIdType >( 100 + comm.rank() ) } }, { comm.rank() }, { .5, .5, .5 }, {} } };
  auto const resolved = comm.resolvePoints( 1, points, oldMaximum );
  ASSERT_EQ( resolved.size(), 4 );
  EXPECT_NE( resolved.at( edge ).globalId, resolved.at( contactSide ).globalId );
  EXPECT_GT( resolved.at( edge ).globalId, oldMaximum );
  ASSERT_EQ( resolved.at( edge ).fields.size(), 1 );
  EXPECT_EQ( resolved.at( edge ).fields[0], 0 );
  auto const low = MpiWrapper::allReduce( resolved.at( edge ).globalId, MpiWrapper::Reduction::Min, MPI_COMM_GEOS );
  auto const high = MpiWrapper::allReduce( resolved.at( edge ).globalId, MpiWrapper::Reduction::Max, MPI_COMM_GEOS );
  EXPECT_EQ( low, high );
  auto const directoryExchanges = comm.statistics().directoryExchanges;
  vtkIdType maximum = oldMaximum;
  for( auto const & record : resolved )
  {
    maximum = std::max( maximum, record.second.globalId );
  }
  EntityKey const childEdge{ 0, EntityKind::edge, { 71, resolved.at( edge ).globalId } };
  auto const fine = comm.resolvePoints( 2, { { childEdge, ranks, { .25, 0, 0 }, {} } }, maximum );
  EXPECT_GT( fine.at( childEdge ).globalId, MpiWrapper::allReduce( maximum, MpiWrapper::Reduction::Max, MPI_COMM_GEOS ) );
  EXPECT_EQ( comm.statistics().directoryExchanges, directoryExchanges );
  EXPECT_EQ( comm.statistics().neighborExchanges, 2 );
  EXPECT_EQ( directoryExchanges, 2 );
  auto const repeat = comm.resolvePoints( 2, { { childEdge, ranks, { .25, 0, 0 }, {} } }, maximum );
  EXPECT_EQ( repeat.at( childEdge ).globalId, fine.at( childEdge ).globalId );
}

TEST( VTKRefinementCommunication, EdgeOnlyAndVertexOnlyParticipants )
{
  Communication comm( MPI_COMM_GEOS, 17 );
  std::vector< EntityKey > entities{ vertex };
  if( comm.rank() < 3 )
  {
    entities.push_back( edge );
  }
  if( comm.rank() < 2 )
  {
    entities.push_back( face );
  }
  auto const sharing = comm.discoverSharing( entities );
  EXPECT_EQ( sharing.at( vertex ), allRanks( comm ) );
  if( comm.rank() < 3 )
  {
    Participants expected( std::min( 3, comm.size() ) );
    std::iota( expected.begin(), expected.end(), 0 );
    EXPECT_EQ( sharing.at( edge ), expected );
  }
  std::vector< PointCreation > points;
  if( comm.rank() < 3 )
  {
    points.push_back( { edge, sharing.at( edge ), { .5, 0, 0 }, {} } );
  }
  if( comm.rank() < 2 )
  {
    points.push_back( { face, sharing.at( face ), { .5, .5, 0 }, {} } );
  }
  auto const result = comm.resolvePoints( 1, points, 100 );
  EXPECT_EQ( result.size(), points.size() );
  if( comm.rank() >= 3 )
  {
    EXPECT_TRUE( result.empty() );
  }
}

TEST( VTKRefinementCommunication, SharedEdgeAllocatorIsAnActualParticipant )
{
  Communication comm( MPI_COMM_GEOS, 13 );
  std::vector< EntityKey > entities{ vertex };
  if( comm.size() == 1 || comm.rank() > 0 )
  {
    entities.push_back( edge );
  }
  auto const sharing = comm.discoverSharing( entities );
  std::vector< PointCreation > points;
  int const allocator = comm.size() == 1 ? 0 : 1;
  if( comm.size() == 1 || comm.rank() > 0 )
  {
    points.push_back( { edge, sharing.at( edge ), { .5, 0, 0 }, { static_cast< unsigned char >( comm.rank() ) } } );
  }
  auto const result = comm.resolvePoints( 1, points, 100 );
  if( !points.empty() )
  {
    EXPECT_EQ( result.at( edge ).fields[0], allocator );
  }
  else
  {
    EXPECT_TRUE( result.empty() );
  }
}

TEST( VTKRefinementCommunication, ReplicatedSurfaceChildrenUseFullIdentitiesAndDisjointRanges )
{
  Communication comm( MPI_COMM_GEOS, 19 );
  vtkIdType const parent = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741019 ) ) : 10001;
  auto const sharing =
      comm.discoverSharing( { entityKey( EntityKind::cell, { parent }, 1 ), entityKey( EntityKind::cell, { parent }, 2 ) } );
  auto const ranks = sharing.at( entityKey( EntityKind::cell, { parent }, 1 ) );
  auto existing = comm.reconcileExistingCells( { { { 1, parent, 204, 0 }, ranks, { static_cast< unsigned char >( comm.rank() ) } },
                                                 { { 2, parent, 204, 0 }, ranks, { static_cast< unsigned char >( comm.rank() ) } } } );
  EXPECT_EQ( existing.size(), 2 );
  for( auto const & [key, record] : existing )
  {
    EXPECT_EQ( record.globalId, key.parentId );
    EXPECT_EQ( record.fields, ( Bytes{ 0 } ) );
  }
  vtkIdType base = 4001;
  for( std::uint64_t ns : { 1, 2 } )
  {
    std::vector< CellCreation > requests;
    for( int child = 3; child >= 0; --child )
    {
      requests.push_back( { { ns, parent, 204, static_cast< std::uint64_t >( child ) },
                            ranks,
                            { static_cast< unsigned char >( comm.rank() ), static_cast< unsigned char >( child ) } } );
    }
    auto const result = comm.resolveCells( 1, requests, base );
    ASSERT_EQ( result.size(), 4 );
    for( auto const & [key, record] : result )
    {
      EXPECT_EQ( record.globalId, base + static_cast< vtkIdType >( key.ordinal ) );
      EXPECT_EQ( record.fields[0], 0 );
      EXPECT_EQ( record.fields[1], key.ordinal );
    }
    base += 4;
  }
  EXPECT_EQ( comm.statistics().directoryExchanges, 2 );
  EXPECT_EQ( comm.statistics().neighborExchanges, 3 );
  std::vector< CellCreation > invalid{ { { 1, parent, 204, 0 }, ranks, {} } };
  if( comm.rank() == comm.size() - 1 )
  {
    invalid[0].key.ordinal = 22;
  }
  EXPECT_THROW( comm.resolveCells( 2, invalid, base ), std::runtime_error );
  invalid[0].key.ordinal = 0;
  EXPECT_NO_THROW( comm.resolveCells( 2, invalid, base ) );
  EXPECT_THROW( comm.resolveCells( 0, invalid, base ), std::runtime_error );
}

TEST( VTKRefinementCommunication, SparseIdRangesEmptyCountsAndExscanRankZero )
{
  Communication comm( MPI_COMM_GEOS );
  std::uint64_t const count = comm.rank() % 2 == 0 ? 3 : 0;
  auto const range = comm.allocateRange( count, 1001 );
  EXPECT_EQ( range.total, 3 * static_cast< std::uint64_t >( ( comm.size() + 1 ) / 2 ) );
  if( count )
  {
    EXPECT_EQ( range.first, 1001 + 3 * ( ( comm.rank() + 1 ) / 2 ) );
  }
  else
  {
    EXPECT_EQ( range.first, 0 );
  }
  auto const empty = comm.allocateRange( 0, std::numeric_limits< vtkIdType >::max() );
  EXPECT_EQ( empty.total, 0 );
  EXPECT_THROW( comm.allocateRange( comm.rank() == 0 ? 2 : 0, std::numeric_limits< vtkIdType >::max() ), std::runtime_error );
  EXPECT_THROW( comm.allocateRange( std::numeric_limits< std::uint64_t >::max(), 0 ), std::runtime_error );
  // Carry between the two count limbs and the largest representable active ID.
  std::uint64_t const wide = UINT64_C( 0x100000001 );
  if( sizeof( vtkIdType ) > 4 )
  {
    auto const carried = comm.allocateRange( wide, 0 );
    EXPECT_EQ( carried.total, wide * comm.size() );
    EXPECT_EQ( carried.first, static_cast< vtkIdType >( wide * comm.rank() ) );
  }
  else
  {
    EXPECT_THROW( comm.allocateRange( wide, 0 ), std::runtime_error );
  }
  auto const last = comm.allocateRange( comm.rank() == 0 ? 1 : 0, std::numeric_limits< vtkIdType >::max() );
  EXPECT_EQ( last.total, 1 );
  if( comm.rank() == 0 )
  {
    EXPECT_EQ( last.first, std::numeric_limits< vtkIdType >::max() );
  }
}

TEST( VTKRefinementCommunication, OneRankValidationFailureDoesNotHang )
{
  Communication comm( MPI_COMM_GEOS );
  EXPECT_THROW( comm.checked( "one-rank fixture",
                              [&]
                              {
                                if( comm.rank() == comm.size() - 1 )
                                {
                                  throw std::invalid_argument( "bad parent 987" );
                                }
                              } ),
                std::runtime_error );
  std::vector< EntityKey > entities{ vertex };
  if( comm.rank() == comm.size() - 1 )
  {
    entities.push_back( { 0, EntityKind::edge, { 81, 71 } } );
  }
  EXPECT_THROW( comm.discoverSharing( entities ), std::runtime_error );
  auto const sharing = comm.discoverSharing( { edge } );
  std::vector< PointCreation > points{ { edge, sharing.at( edge ), { .5, 0, 0 }, {} } };
  if( comm.rank() == comm.size() - 1 )
  {
    points[0].position[0] = std::numeric_limits< double >::quiet_NaN();
  }
  EXPECT_THROW( comm.resolvePoints( 1, points, 100 ), std::runtime_error );
  points[0].position = { .5, 0, 0 };
  EXPECT_NO_THROW( comm.resolvePoints( 1, points, 100 ) );
}

TEST( VTKRefinementCommunication, CrossedSharedFaceCyclesAreRejected )
{
  Communication comm( MPI_COMM_GEOS );
  EntityKey crossed = face;
  crossed.corners = canonicalCycle( { 71, 91, 81, 101 } );
  if( comm.size() == 1 )
  {
    EXPECT_THROW( comm.discoverSharing( { face, crossed } ), std::runtime_error );
  }
  else
  {
    EXPECT_THROW( comm.discoverSharing( { comm.rank() == 0 ? face : crossed } ), std::runtime_error );
  }
}

TEST( VTKRefinementCommunication, DistinctCommunicatorAndZeroDegreeRanks )
{
  int const worldRank = MpiWrapper::commRank( MPI_COMM_GEOS );
  MPI_Comm sub = MPI_COMM_GEOS;
#ifdef GEOS_USE_MPI
  MPI_Comm_split( MPI_COMM_GEOS, worldRank % 2, worldRank, &sub );
#endif
  {
    Communication comm( sub, 7 );
    auto const sharing = comm.discoverSharing( { vertex } );
    EXPECT_EQ( sharing.at( vertex ), allRanks( comm ) );
    EXPECT_EQ( comm.allocateRange( 1, 300 ).total, comm.size() );
    auto const empty = comm.discoverSharing( {} );
    EXPECT_TRUE( empty.empty() );
    EXPECT_TRUE( comm.neighbors().empty() );
    auto const result = comm.resolvePoints( 1, { { { 0, EntityKind::cell, { comm.rank() } }, { comm.rank() }, { 0, 0, 0 }, {} } }, -1 );
    EXPECT_EQ( result.size(), 1 );
  }
#ifdef GEOS_USE_MPI
  MpiWrapper::commFree( sub );
#endif
  GEOS_UNUSED_VAR( worldRank );
}

TEST( VTKRefinementCommunication, ChunkTransportUnequalAndZeroLengths )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  if( size < 2 )
  {
    GTEST_SKIP() << "Requires two ranks";
  }
  if( rank < 2 )
  {
    Bytes send( rank == 0 ? 31 : 5, static_cast< unsigned char >( rank + 1 ) ), receive( rank == 0 ? 5 : 31 );
    mpi::exchangeBytes( send.data(), send.size(), receive.data(), receive.size(), 1 - rank, 33, MPI_COMM_GEOS, 3 );
    EXPECT_TRUE( std::all_of( receive.begin(), receive.end(), [&]( unsigned char x ) { return x == 2 - rank; } ) );
    Bytes empty, nonempty( 19, 7 ), received( rank == 0 ? 19 : 0 );
    Bytes const & output = rank == 0 ? empty : nonempty;
    mpi::exchangeBytes( output.data(), output.size(), received.data(), received.size(), 1 - rank, 34, MPI_COMM_GEOS, 4 );
    if( rank == 0 )
    {
      EXPECT_EQ( received, nonempty );
    }
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );
}

TEST( VTKRefinementCommunication, ConcurrentChunkTransportUnequalPeersAndIsolatedRank )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  int const active = size > 2 ? size - 1 : size;
  std::vector< Bytes > send( size ), receive( size );
  std::vector< mpi::ByteExchange > exchanges;
  if( rank < active )
  {
    for( int peer = 0; peer < active; ++peer )
    {
      if( peer != rank )
      {
        // Deliberately unequal chunk counts, with one direction sometimes empty.
        auto length = []( int from, int to ) { return from % 2 ? 0 : ( from + 1 ) * ( to + 2 ) + 13; };
        send[peer].resize( length( rank, peer ), static_cast< unsigned char >( rank + 1 ) );
        receive[peer].resize( length( peer, rank ) );
        exchanges.push_back( { send[peer].data(), send[peer].size(), receive[peer].data(), receive[peer].size(), peer } );
      }
    }
  }
  std::vector< MPI_Request > requests( 2 * exchanges.size() );
  mpi::exchangeManyBytes( exchanges, requests.data(), 37, MPI_COMM_GEOS, 5 );
  for( int peer = 0; peer < size; ++peer )
  {
    EXPECT_TRUE( std::all_of( receive[peer].begin(), receive[peer].end(), [&]( unsigned char value ) { return value == peer + 1; } ) );
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );
}

TEST( VTKRefinementCommunication, RemoteSurfaceSidesAndFineSupportIdsUseCachedGraph )
{
  Communication comm( MPI_COMM_GEOS, 7 );
  comm.discoverSharing( {} );
  vtkIdType const base = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741001 ) ) : 1001;
  std::vector< MainFace > localFaces;
  std::vector< CoarseSurface > queries;
  std::unordered_map< EntityKey, vtkIdType, EntityKeyHash > localIds;
  std::vector< Connectivity > buckets( 4 );
  for( int s = 0; s < 3; ++s )
  {
    int const owner = comm.size() == 1 ? 0 : 1 + s % std::min( 2, comm.size() - 1 );
    Connectivity corners;
    for( int c = 0; c < 4; ++c )
    {
      vtkIdType const id = base + 10 * s + c;
      corners.push_back( id );
      buckets[c].push_back( id );
      if( owner == comm.rank() )
      {
        localIds[entityKey( EntityKind::vertex, { id } )] = id;
      }
    }
    if( owner == comm.rank() )
    {
      localFaces.push_back( { corners, { owner } } );
      for( int c = 0; c < 4; ++c )
      {
        localIds[entityKey( EntityKind::edge, { corners[c], corners[( c + 1 ) % 4] } )] = base + 100 + 10 * s + c;
      }
      localIds[entityKey( EntityKind::face, corners )] = base + 200 + s;
    }
  }
  if( comm.rank() == 0 )
  {
    queries.push_back( { buckets } );
  }
  auto const sides = comm.discoverSurfaceSides( localFaces, queries );
  auto const coarseDirectoryExchanges = comm.statistics().directoryExchanges;
  EXPECT_EQ( coarseDirectoryExchanges, 4 );
  if( comm.rank() == 0 )
  {
    ASSERT_EQ( sides.size(), 1 );
    ASSERT_EQ( sides[0].size(), 3 );
    // This rank owns no volume in the distributed case. It must not allocate
    // any main-side point or be added to the main point's participant list.
    for( auto const & side : sides[0] )
    {
      if( comm.size() > 1 )
      {
        EXPECT_NE( side.owners.front(), 0 );
      }
    }
  }
  else
  {
    EXPECT_TRUE( sides.empty() );
  }
  if( comm.rank() > 2 )
  {
    EXPECT_TRUE( comm.neighbors().empty() );
  }
  std::vector< SupportLookup > requests;
  std::vector< vtkIdType > expected;
  if( comm.rank() == 0 )
  {
    for( std::size_t s = 0; s < sides[0].size(); ++s )
    {
      auto const & side = sides[0][s];
      requests.push_back( { SurfaceAssociations::supportKey( EntityKind::vertex, { 0 }, { 0, 1, 2, 3 }, side ), side.owners.front() } );
      requests.push_back( { SurfaceAssociations::supportKey( EntityKind::edge, { 0, 1 }, { 0, 1, 2, 3 }, side ), side.owners.front() } );
      requests.push_back( { side.mainFace, side.owners.front() } );
      expected.insert( expected.end(), { base + 10 * static_cast< vtkIdType >( s ), base + 100 + 10 * static_cast< vtkIdType >( s ),
                                         base + 200 + static_cast< vtkIdType >( s ) } );
    }
  }
  EXPECT_EQ( comm.resolveSupportIds( 1, requests, localIds ), expected );
  EXPECT_EQ( comm.statistics().directoryExchanges, coarseDirectoryExchanges );
  for( auto & [key, id] : localIds )
  {
    if( key.kind != EntityKind::vertex )
    {
      ++id;
    }
  }
  for( std::size_t i = 0; i < expected.size(); ++i )
  {
    if( i % 3 != 0 )
    {
      ++expected[i];
    }
  }
  EXPECT_EQ( comm.resolveSupportIds( 2, requests, localIds ), expected );
  EXPECT_EQ( comm.statistics().directoryExchanges, coarseDirectoryExchanges );
  auto missing = localIds;
  if( comm.rank() == ( comm.size() == 1 ? 0 : 1 ) )
  {
    missing.erase( entityKey( EntityKind::edge, { base, base + 1 } ) );
  }
  EXPECT_THROW( comm.resolveSupportIds( 3, requests, missing ), std::runtime_error );
  EXPECT_EQ( comm.resolveSupportIds( 3, requests, localIds ), expected );
  if( comm.rank() == 0 )
  {
    queries[0].cornerBuckets[0] = { base + 3000 };
  }
  EXPECT_THROW( comm.discoverSurfaceSides( localFaces, queries ), std::runtime_error );
}

TEST( VTKRefinementCommunication, CoarseFullFaceValidationFindsHiddenNonmanifoldIncidence )
{
  Communication comm( MPI_COMM_GEOS, 5 );
  Connectivity const forward{ 71, 81, 91, 101 }, backward{ 71, 101, 91, 81 };
  std::vector< MainFace > faces;
  if( comm.rank() == 0 )
  {
    faces.push_back( { forward, { 0 } } );
  }
  if( comm.rank() == ( comm.size() == 1 ? 0 : 1 ) )
  {
    faces.push_back( { backward, { comm.rank() } } );
  }
  EXPECT_NO_THROW( comm.validateVolumeFaces( faces ) );
  EXPECT_TRUE( comm.neighbors().empty() );
  EXPECT_EQ( comm.statistics().directoryExchanges, 1 );
  // Rank zero's two incidences make the face locally internal, so boundary-only
  // discovery cannot see the third incidence owned by another rank.
  faces.clear();
  if( comm.rank() == 0 )
  {
    faces.push_back( { forward, { 0 } } );
    faces.push_back( { backward, { 0 } } );
  }
  if( comm.rank() == ( comm.size() == 1 ? 0 : 1 ) )
  {
    faces.push_back( { forward, { comm.rank() } } );
  }
  EXPECT_THROW( comm.validateVolumeFaces( faces ), std::runtime_error );
  faces.clear();
  if( comm.rank() == 0 )
  {
    faces.push_back( { forward, { 0 } } );
  }
  if( comm.rank() == ( comm.size() == 1 ? 0 : 1 ) )
  {
    faces.push_back( { forward, { comm.rank() } } );
  }
  EXPECT_THROW( comm.validateVolumeFaces( faces ), std::runtime_error );
  faces.clear();
  if( comm.rank() == 0 )
  {
    faces.push_back( { forward, { 0 } } );
  }
  if( comm.rank() == ( comm.size() == 1 ? 0 : 1 ) )
  {
    faces.push_back( { { 71, 91, 81, 101 }, { comm.rank() } } );
  }
  EXPECT_THROW( comm.validateVolumeFaces( faces ), std::runtime_error );
  EXPECT_NO_THROW( comm.validateVolumeFaces( {} ) );
}

int main( int argc, char ** argv )
{
  MpiWrapper::init( &argc, &argv );
  MPI_COMM_GEOS = MpiWrapper::commDup( MPI_COMM_WORLD );
  ::testing::InitGoogleTest( &argc, argv );
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) != 0 )
  {
    auto & listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release( listeners.default_result_printer() );
  }
  int const result = RUN_ALL_TESTS();
  MpiWrapper::finalize();
  return result;
}
