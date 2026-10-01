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

/** @file testVTKUniformRefinement.cpp */
#include "mesh/generators/VTKUniformRefinement.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include "VTKRefinementTestMeshes.hpp"

#include <gtest/gtest.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetReader.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedCharArray.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>

using namespace geos;
using namespace geos::vtk;
using namespace geos::vtk::refinement;

namespace
{
vtkIdType const sparseBase = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741001 ) ) : 1001;

vtkSmartPointer< vtkUnstructuredGrid > makeGrid( std::vector< Coordinates > const & xyz, Connectivity const & pointIds,
                                                 std::vector< Cell > const & cells, Connectivity const & cellIds )
{
  auto grid = vtkSmartPointer< vtkUnstructuredGrid >::New();
  vtkNew< vtkPoints > points;
  points->SetDataTypeToDouble();
  for( auto const & p : xyz )
  {
    points->InsertNextPoint( p.data() );
  }
  grid->SetPoints( points );
  for( auto const & cell : cells )
  {
    grid->InsertNextCell( cell.vtkType, cell.points.size(), cell.points.data() );
  }
  vtkNew< vtkIdTypeArray > pids, cids;
  pids->SetName( "inputPointIds" );
  cids->SetName( "inputCellIds" );
  for( vtkIdType id : pointIds )
  {
    pids->InsertNextValue( id );
  }
  for( vtkIdType id : cellIds )
  {
    cids->InsertNextValue( id );
  }
  grid->GetPointData()->SetGlobalIds( pids );
  grid->GetCellData()->SetGlobalIds( cids );
  vtkNew< vtkDoubleArray > affine;
  affine->SetName( "affine" );
  for( auto const & p : xyz )
  {
    affine->InsertNextValue( 2 * p[0] - 3 * p[1] + 4 * p[2] + 5 );
  }
  grid->GetPointData()->SetScalars( affine );
  vtkNew< vtkIntArray > attributes;
  attributes->SetName( "attribute" );
  for( auto const & cell : cells )
  {
    attributes->InsertNextValue( cell.vtkType == VTK_QUAD ? 77 : 1 );
  }
  grid->GetCellData()->AddArray( attributes );
  return grid;
}

vtkSmartPointer< vtkUnstructuredGrid > localCube( int x, bool marker )
{
  auto r = testMeshes::referenceCell( VTK_HEXAHEDRON, 0 );
  Connectivity ids;
  for( auto & p : r.xyz )
  {
    p[0] += x;
    ids.push_back( sparseBase +
                   7 * ( 4 * static_cast< vtkIdType >( p[0] ) + 2 * static_cast< vtkIdType >( p[1] ) + static_cast< vtkIdType >( p[2] ) ) );
  }
  std::vector< Cell > cells{ r.cell };
  Connectivity cellIds{ sparseBase + 1000 + 2 * x };
  if( marker )
  {
    cells.push_back( { VTK_QUAD, { 4, 5, 6, 7 }, 0 } );
    cellIds.push_back( sparseBase + 1001 + 2 * x );
  }
  return makeGrid( r.xyz, ids, cells, cellIds );
}

void checkAffine( vtkDataSet & grid )
{
  auto * data = vtkDoubleArray::SafeDownCast( grid.GetPointData()->GetScalars() );
  ASSERT_NE( data, nullptr );
  for( vtkIdType p = 0; p < grid.GetNumberOfPoints(); ++p )
  {
    Coordinates xyz{};
    grid.GetPoint( p, xyz.data() );
    EXPECT_NEAR( data->GetValue( p ), 2 * xyz[0] - 3 * xyz[1] + 4 * xyz[2] + 5, 1e-12 );
  }
}
} // namespace

TEST( VTKUniformRefinement, ZeroBypassesInvalidDataAndLeavesPointersAndArrays )
{
  vtkNew< vtkUnstructuredGrid > input;
  AllMeshes meshes( input, {} );
  auto const result = refineUniformly( meshes, 0, {}, MPI_COMM_NULL );
  EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
  EXPECT_EQ( input->GetPointData()->GetNumberOfArrays(), 0 );
  EXPECT_EQ( input->GetCellData()->GetNumberOfArrays(), 0 );
  EXPECT_TRUE( result.blocks.empty() );
  EXPECT_TRUE( result.neighbors.empty() );
  EXPECT_TRUE( result.resources.empty() );
  EXPECT_TRUE( result.levels.empty() );
}

TEST( VTKUniformRefinement, ConnectedHexesMarkersFieldsAndLineageThroughTwoLevels )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  auto input = localCube( rank, true );
  AllMeshes meshes( input, {} );
  auto const result = refineUniformly( meshes, 2, {}, MPI_COMM_GEOS );
  auto output = meshes.getMainMesh();
  EXPECT_NE( output.GetPointer(), input.GetPointer() );
  EXPECT_EQ( input->GetNumberOfCells(), 2 );
  EXPECT_EQ( output->GetNumberOfCells(), 80 );
  EXPECT_EQ( output->GetNumberOfPoints(), 125 );
  ASSERT_EQ( result.blocks.size(), 1 );
  EXPECT_EQ( result.blocks[0].sourceName, "1_hexahedra" );
  EXPECT_EQ( result.blocks[0].name, "1_hexahedra" );
  EXPECT_EQ( result.blocks[0].cells.size(), 64 );
  checkAffine( *output );
  auto * pids = vtkIdTypeArray::SafeDownCast( output->GetPointData()->GetGlobalIds() );
  auto * oldIds = vtkIdTypeArray::SafeDownCast( input->GetPointData()->GetGlobalIds() );
  for( int p = 0; p < 8; ++p )
  {
    EXPECT_EQ( pids->GetValue( p ), oldIds->GetValue( p ) );
  }
  auto * cids = vtkIdTypeArray::SafeDownCast( output->GetCellData()->GetGlobalIds() );
  auto * roots = vtkIdTypeArray::SafeDownCast( output->GetCellData()->GetArray( "_geosUniformRootCellId" ) );
  auto * generations = vtkIdTypeArray::SafeDownCast( output->GetCellData()->GetArray( "_geosUniformGeneration" ) );
  auto * owners = vtkIdTypeArray::SafeDownCast( output->GetCellData()->GetArray( "_geosUniformRootOwner" ) );
  auto * attributes = vtkIntArray::SafeDownCast( output->GetCellData()->GetArray( "attribute" ) );
  std::set< vtkIdType > unique;
  for( vtkIdType c = 0; c < output->GetNumberOfCells(); ++c )
  {
    bool const surface = output->GetCellType( c ) == VTK_QUAD;
    EXPECT_EQ( roots->GetValue( c ), sparseBase + 1000 + 2 * rank + surface );
    EXPECT_EQ( attributes->GetValue( c ), surface ? 77 : 1 );
    EXPECT_EQ( generations->GetValue( c ), 2 );
    EXPECT_EQ( owners->GetValue( c ), rank );
    EXPECT_EQ( cids->GetValue( c ) >= size * 64, surface );
    EXPECT_TRUE( unique.insert( cids->GetValue( c ) ).second );
  }
  EXPECT_EQ( result.communication.directoryExchanges, 5 );
  ASSERT_EQ( result.resources.size(), 2 );
  ASSERT_EQ( result.levels.size(), 2 );
  EXPECT_EQ( result.resources[0].volumeCells, 8 );
  EXPECT_EQ( result.resources[1].volumeCells, 64 );
  EXPECT_EQ( result.resources[1].surfaceCellCopies, 16 );
  EXPECT_GE( result.resources[1].pointCopiesUpperBound, static_cast< std::uint64_t >( output->GetNumberOfPoints() ) );
  EXPECT_GT( result.resources[1].fieldBytesUpperBound, result.resources[0].fieldBytesUpperBound );
  EXPECT_EQ( result.levels[0].mainPointCopies, 27 );
  EXPECT_EQ( result.levels[1].mainPointCopies, 125 );
  EXPECT_EQ( MpiWrapper::sum( result.levels[1].ownedMainPoints, MPI_COMM_GEOS ), 100 * size + 25 );
  std::uint64_t bytes = result.coarseCommunication.payloadBytesSent;
  for( auto const & level : result.levels )
  {
    EXPECT_EQ( level.communication.directoryExchanges, 0 );
    bytes += level.communication.payloadBytesSent;
  }
  EXPECT_EQ( bytes, result.communication.payloadBytesSent );
  if( size > 1 )
  {
    EXPECT_FALSE( result.neighbors.empty() );
  }
}

TEST( VTKUniformRefinement, ForecastsMatchGrowthAndCountSharedUnusedOriginalPoints )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  auto input = localCube( rank, false );
  input->GetPoints()->InsertNextPoint( -3, -3, -3 );
  vtkIdTypeArray::SafeDownCast( input->GetPointData()->GetGlobalIds() )->InsertNextValue( sparseBase - 7 );
  vtkDoubleArray::SafeDownCast( input->GetPointData()->GetScalars() )->InsertNextValue( -4 );
  AllMeshes meshes( input, {} );
  UniformRefinementOptions options;
  options.reportStatistics = true;
  auto const result = refineUniformly( meshes, 2, options, MPI_COMM_GEOS );
  auto output = meshes.getMainMesh();
  ASSERT_EQ( result.resources.size(), 2 );
  EXPECT_EQ( result.resources[0].volumeCells, 8 );
  EXPECT_EQ( result.resources[1].volumeCells, 64 );
  EXPECT_EQ( result.resources[0].pointCopiesUpperBound, 28 );
  EXPECT_EQ( result.resources[1].pointCopiesUpperBound, 180 );
  std::uint64_t connectivity = 0;
  for( vtkIdType c = 0; c < output->GetNumberOfCells(); ++c ) connectivity += output->GetCell( c )->GetNumberOfPoints();
  EXPECT_EQ( result.resources[1].connectivityEntries, connectivity );
  EXPECT_GT( result.resources[1].modeledRefinerPeakBytes, result.resources[1].vtkBytesUpperBound );
  EXPECT_TRUE( result.resources[1].geosGhostConnectivityBytesModel.has_value() );
  EXPECT_EQ( output->GetNumberOfPoints(), 126 );
  checkAffine( *output );
  EXPECT_EQ( MpiWrapper::sum( result.levels[1].ownedMainPoints, MPI_COMM_GEOS ), 100 * size + 26 );
  EXPECT_EQ( MpiWrapper::sum( result.levels[1].sharedMainPointCopies, MPI_COMM_GEOS ), size == 1 ? 0 : 50 * ( size - 1 ) + size );
}

TEST( VTKUniformRefinement, UnusedOriginalCopyMatchesAnInteriorVertexOnItsOtherRank )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  vtkIdType const centerId = sparseBase + 20000 + 7 * 13;
  vtkSmartPointer< vtkUnstructuredGrid > input;
  if( rank == 0 )
  {
    std::vector< Coordinates > xyz;
    Connectivity ids, cellIds;
    std::vector< Cell > cells;
    for( int z = 0; z < 3; ++z )
      for( int y = 0; y < 3; ++y )
        for( int x = 0; x < 3; ++x )
        {
          ids.push_back( sparseBase + 20000 + 7 * xyz.size() );
          xyz.push_back( { double( x - 2 ), double( y - 1 ), double( z - 1 ) } );
        }
    auto index = []( int x, int y, int z ) { return vtkIdType( 9 * z + 3 * y + x ); };
    for( int z = 0; z < 2; ++z )
      for( int y = 0; y < 2; ++y )
        for( int x = 0; x < 2; ++x )
        {
          cellIds.push_back( sparseBase + 30000 + cells.size() );
          cells.push_back( { VTK_HEXAHEDRON,
            { index( x, y, z ), index( x + 1, y, z ), index( x + 1, y + 1, z ), index( x, y + 1, z ),
              index( x, y, z + 1 ), index( x + 1, y, z + 1 ), index( x + 1, y + 1, z + 1 ), index( x, y + 1, z + 1 ) }, 0 } );
        }
    input = makeGrid( xyz, ids, cells, cellIds );
  }
  else
  {
    input = localCube( 3 * rank, false );
    input->GetPoints()->InsertNextPoint( -1, 0, 0 );
    vtkIdTypeArray::SafeDownCast( input->GetPointData()->GetGlobalIds() )->InsertNextValue( centerId );
    vtkDoubleArray::SafeDownCast( input->GetPointData()->GetScalars() )->InsertNextValue( 3 );
  }
  AllMeshes meshes( input, {} );
  auto const result = refineUniformly( meshes, 1, {}, MPI_COMM_GEOS );
  ASSERT_EQ( result.levels.size(), 1 );
  EXPECT_EQ( result.levels[0].mainPointCopies, rank == 0 ? 125 : 28 );
  EXPECT_EQ( result.levels[0].sharedMainPointCopies, size == 1 ? 0 : 1 );
  EXPECT_EQ( MpiWrapper::sum( result.levels[0].ownedMainPoints, MPI_COMM_GEOS ), 125 + 27 * ( size - 1 ) );
  EXPECT_EQ( result.levels[0].communication.directoryExchanges, 0 );
  checkAffine( *meshes.getMainMesh() );
}

TEST( VTKUniformRefinement, HugeLevelsWithEmptyRanksFailBeforeFineWork )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  auto input = rank == 0 ? localCube( 0, false ) : makeGrid( {}, {}, {}, {} );
  AllMeshes meshes( input, {} );
  EXPECT_THROW( refineUniformly( meshes, std::numeric_limits< int >::max(), {}, MPI_COMM_GEOS ), std::runtime_error );
  EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
  EXPECT_EQ( input->GetNumberOfCells(), rank == 0 ? 1 : 0 );
}

TEST( VTKUniformRefinement, DifferentRankLevelsAndReportingFailCollectively )
{
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) == 1 )
  {
    GTEST_SKIP() << "Requires multiple ranks";
  }
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  auto input = localCube( rank, false );
  AllMeshes meshes( input, {} );
  EXPECT_THROW( refineUniformly( meshes, rank == 0 ? 1 : 2, {}, MPI_COMM_GEOS ), std::runtime_error );
  UniformRefinementOptions options;
  options.reportStatistics = rank == 0;
  EXPECT_THROW( refineUniformly( meshes, 1, options, MPI_COMM_GEOS ), std::runtime_error );
  EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
}

TEST( VTKUniformRefinement, SourcePyramidsRemainSeparateFromOriginalTetrahedra )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  auto pyramid = testMeshes::referenceCell( VTK_PYRAMID, 0 );
  auto tet = testMeshes::referenceCell( VTK_TETRA, 0 );
  std::vector< Coordinates > xyz = pyramid.xyz;
  for( auto p : tet.xyz )
  {
    p[0] += 3;
    xyz.push_back( p );
  }
  for( auto & p : xyz )
  {
    p[0] += 10 * rank;
  }
  for( auto & p : tet.cell.points )
  {
    p += pyramid.xyz.size();
  }
  Connectivity pids;
  for( std::size_t p = 0; p < xyz.size(); ++p )
  {
    pids.push_back( sparseBase + 100 * rank + p );
  }
  auto input = makeGrid( xyz, pids, { pyramid.cell, tet.cell }, { sparseBase + 1000 + 2 * rank, sparseBase + 1001 + 2 * rank } );
  AllMeshes meshes( input, {} );
  auto const result = refineUniformly( meshes, 2, {}, MPI_COMM_GEOS );
  std::map< std::string, std::size_t > counts;
  for( auto const & block : result.blocks )
  {
    counts[block.name] = block.cells.size();
  }
  EXPECT_EQ( counts["1_pyramids"], 36 );
  EXPECT_EQ( counts["1_pyramids__refined_tetrahedra"], 56 );
  EXPECT_EQ( counts["1_tetrahedra"], 64 );
  EXPECT_EQ( meshes.getMainMesh()->GetNumberOfCells(), 156 );
}

TEST( VTKUniformRefinement, FractureRemoteSidesReplicasAndDisjointNamedBlocks )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  std::vector< Coordinates > xyz;
  Connectivity pids, cids;
  std::vector< Cell > cells;
  std::vector< Connectivity > buckets( 4 );
  for( int side = 0; side < 2; ++side )
  {
    auto cube = testMeshes::referenceCell( VTK_HEXAHEDRON, 0 );
    Connectivity const face = side == 0 ? Connectivity{ 1, 2, 6, 5 } : Connectivity{ 0, 3, 7, 4 };
    for( int c = 0; c < 4; ++c )
    {
      buckets[c].push_back( sparseBase + 100 * side + face[c] );
    }
    if( rank == ( size == 1 ? 0 : side ) )
    {
      vtkIdType const offset = xyz.size();
      for( std::size_t p = 0; p < cube.xyz.size(); ++p )
      {
        cube.xyz[p][0] += side;
        xyz.push_back( cube.xyz[p] );
        pids.push_back( sparseBase + 100 * side + p );
      }
      for( auto & p : cube.cell.points )
      {
        p += offset;
      }
      cells.push_back( cube.cell );
      cids.push_back( sparseBase + 1000 + side );
    }
  }
  auto main = makeGrid( xyz, pids, cells, cids );
  stdMap< string, vtkSmartPointer< vtkDataSet > > blocks;
  if( rank < std::min( size, 2 ) )
  {
    for( int block = 0; block < 2; ++block )
    {
      auto aux = makeGrid( { { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 }, { 1, 0, 1 } }, { 701, 721, 711, 731 },
                           { { VTK_QUAD, rank ? Connectivity{ 2, 1, 0, 3 } : Connectivity{ 0, 1, 2, 3 }, 0 } }, { sparseBase + 7000 } );
      vtkNew< vtkIdTypeArray > collocation;
      collocation->SetName( "collocated_nodes" );
      collocation->SetNumberOfComponents( 3 );
      collocation->SetNumberOfTuples( 4 );
      for( int p = 0; p < 4; ++p )
      {
        collocation->SetTypedComponent( p, 0, buckets[p][0] );
        collocation->SetTypedComponent( p, 1, buckets[p][1] );
        collocation->SetTypedComponent( p, 2, -1 );
      }
      aux->GetPointData()->AddArray( collocation );
      vtkNew< vtkStringArray > label;
      label->SetName( "label" );
      label->InsertNextValue( rank ? "replica" : "allocator" );
      aux->GetCellData()->AddArray( label );
      blocks.emplace( "fracture" + std::to_string( block ), aux );
    }
  }
  AllMeshes meshes( main, blocks );
  auto const result = refineUniformly( meshes, 2, {}, MPI_COMM_GEOS );
  ASSERT_EQ( meshes.getFaceBlocks().size(), 2 );
  std::set< vtkIdType > auxiliaryIds;
  std::map< vtkIdType, Coordinates > mainCoordinates;
  auto output = meshes.getMainMesh();
  auto * mainIds = vtkIdTypeArray::SafeDownCast( output->GetPointData()->GetGlobalIds() );
  for( vtkIdType p = 0; p < output->GetNumberOfPoints(); ++p )
  {
    Coordinates position{};
    output->GetPoint( p, position.data() );
    mainCoordinates.emplace( mainIds->GetValue( p ), position );
  }
  for( auto const & [name, aux] : meshes.getFaceBlocks() )
  {
    GEOS_UNUSED_VAR( name );
    EXPECT_EQ( aux->GetNumberOfCells(), rank < std::min( size, 2 ) ? 16 : 0 );
    EXPECT_EQ( aux->GetNumberOfPoints(), rank < std::min( size, 2 ) ? 25 : 0 );
    if( aux->GetNumberOfCells() == 0 )
    {
      continue;
    }
    checkAffine( *aux );
    auto * ids = vtkIdTypeArray::SafeDownCast( aux->GetCellData()->GetGlobalIds() );
    auto * label = vtkStringArray::SafeDownCast( aux->GetCellData()->GetAbstractArray( "label" ) );
    auto * collocation = vtkIdTypeArray::SafeDownCast( aux->GetPointData()->GetArray( "collocated_nodes" ) );
    ASSERT_EQ( collocation->GetNumberOfComponents(), 2 );
    for( vtkIdType c = 0; c < aux->GetNumberOfCells(); ++c )
    {
      EXPECT_TRUE( auxiliaryIds.insert( ids->GetValue( c ) ).second );
      EXPECT_EQ( label->GetValue( c ), "allocator" );
    }
    for( vtkIdType p = 0; p < aux->GetNumberOfPoints(); ++p )
    {
      Coordinates position{};
      aux->GetPoint( p, position.data() );
      for( int b = 0; b < 2; ++b )
      {
        auto const found = mainCoordinates.find( collocation->GetTypedComponent( p, b ) );
        if( found != mainCoordinates.end() )
        {
          EXPECT_EQ( found->second, position );
        }
      }
      if( size == 1 )
      {
        EXPECT_TRUE( mainCoordinates.count( collocation->GetTypedComponent( p, 0 ) ) );
        EXPECT_TRUE( mainCoordinates.count( collocation->GetTypedComponent( p, 1 ) ) );
      }
    }
  }
  EXPECT_EQ( result.communication.directoryExchanges, 5 );
}

TEST( VTKUniformRefinement, OneRankFailureIsCollectiveAndDoesNotReplaceInput )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  auto input = localCube( rank, false );
  AllMeshes meshes( input, {} );
  if( rank == 0 )
  {
    input->GetPointData()->SetGlobalIds( nullptr );
  }
  EXPECT_THROW( refineUniformly( meshes, 1, {}, MPI_COMM_GEOS ), std::runtime_error );
  EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
  EXPECT_EQ( input->GetNumberOfCells(), 1 );
  EXPECT_THROW( refineUniformly( meshes, -1, {}, MPI_COMM_GEOS ), std::runtime_error );
}

TEST( VTKUniformRefinement, InvalidCellDiagnosticsIncludeBlockLevelParentAndType )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  int const failingRank = MpiWrapper::commSize( MPI_COMM_GEOS ) - 1;
  for( int scenario = 0; scenario < 4; ++scenario )
  {
    bool const auxiliary = scenario == 2;
    bool const fineFailure = scenario == 1;
    auto input = localCube( 3 * rank, false );
    AllMeshes meshes( input, {} );
    vtkIdType const parentId = auxiliary ? sparseBase + 2000 : sparseBase + 1000 + 6 * failingRank;
    if( rank == failingRank )
    {
      if( auxiliary )
      {
        // Unsupported auxiliary lines fail before relational metadata is read.
        meshes.getFaceBlocks()["fault"] = makeGrid( { { 0, 0, 0 }, { 1, 0, 0 } }, { 3, 4 },
                                                    { { VTK_LINE, { 0, 1 }, 0 } }, { parentId } );
      }
      else if( scenario == 3 )
      {
        vtkNew< vtkUnsignedCharArray > ghosts;
        ghosts->SetName( vtkDataSetAttributes::GhostArrayName() );
        ghosts->InsertNextValue( vtkDataSetAttributes::DUPLICATECELL );
        input->GetCellData()->AddArray( ghosts );
      }
      else
      {
        for( vtkIdType p = 0; p < input->GetNumberOfPoints(); ++p )
        {
          Coordinates xyz{};
          input->GetPoint( p, xyz.data() );
          if( fineFailure )
          {
            // The width-two coarse cell is valid, but its x midpoint cannot
            // be represented at 1e16. Child validation must diagnose level one.
            xyz[0] = 1e16 + 2 * ( xyz[0] - 3 * rank );
          }
          else xyz[2] = 0;
          input->GetPoints()->SetPoint( p, xyz.data() );
        }
      }
    }
    try
    {
      ( void )refineUniformly( meshes, 2, {}, MPI_COMM_GEOS );
      FAIL() << "The rank-local invalid cell must fail on every rank";
    }
    catch( std::runtime_error const & error )
    {
      std::string const diagnostic = error.what();
      if( scenario == 3 )
      {
        EXPECT_NE( diagnostic.find( "unnormalized cell ghosts" ), std::string::npos );
      }
      EXPECT_NE( diagnostic.find( "rank " + std::to_string( failingRank ) ), std::string::npos );
      EXPECT_NE( diagnostic.find( auxiliary ? "block 'fault'" : "block 'main'" ), std::string::npos );
      EXPECT_NE( diagnostic.find( fineFailure ? "level 1" : "level 0" ), std::string::npos );
      EXPECT_NE( diagnostic.find( "parent global ID " + std::to_string( parentId ) ), std::string::npos );
      EXPECT_NE( diagnostic.find( "VTK cell type " + std::to_string( auxiliary ? VTK_LINE : VTK_HEXAHEDRON ) ),
                 std::string::npos );
    }
    EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
    EXPECT_EQ( input->GetNumberOfCells(), 1 );
  }
}

TEST( VTKUniformRefinement, RequiredPointArraysIgnoreUnimportedLabels )
{
  auto input = localCube( MpiWrapper::commRank( MPI_COMM_GEOS ), false );
  vtkNew< vtkStringArray > strings;
  strings->SetName( "unusedStrings" );
  vtkNew< vtkIntArray > labels;
  labels->SetName( "unusedLabels" );
  for( vtkIdType p = 0; p < input->GetNumberOfPoints(); ++p )
  {
    strings->InsertNextValue( std::to_string( p ) );
    labels->InsertNextValue( p );
  }
  input->GetPointData()->AddArray( strings );
  input->GetPointData()->AddArray( labels );
  UniformRefinementOptions options;
  options.requiredPointArrays = std::set< std::string >{ "affine" };
  AllMeshes meshes( input, {} );
  EXPECT_NO_THROW( refineUniformly( meshes, 2, options, MPI_COMM_GEOS ) );
  EXPECT_EQ( meshes.getMainMesh()->GetPointData()->GetAbstractArray( "unusedStrings" ), nullptr );
  EXPECT_EQ( meshes.getMainMesh()->GetPointData()->GetArray( "unusedLabels" ), nullptr );
  checkAffine( *meshes.getMainMesh() );
  options.requiredPointArrays->insert( "unusedLabels" );
  AllMeshes invalid( input, {} );
  EXPECT_THROW( refineUniformly( invalid, 1, options, MPI_COMM_GEOS ), std::runtime_error );
}

TEST( VTKUniformRefinement, RequiredCellArraysKeepVolumeAndSurfaceImportsWithoutUnusedCopies )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  std::vector< Coordinates > coordinates;
  Connectivity pointIds, cellIds;
  std::vector< Cell > cells;
  std::array< Connectivity, 2 > const sideFaces{ Connectivity{ 1, 2, 6, 5 }, Connectivity{ 0, 3, 7, 4 } };
  for( int side = 0; side < 2; ++side )
  {
    auto cube = testMeshes::referenceCell( VTK_HEXAHEDRON, 0 );
    vtkIdType const offset = coordinates.size();
    for( std::size_t p = 0; p < cube.xyz.size(); ++p )
    {
      cube.xyz[p][0] += side;
      cube.xyz[p][1] += 2 * rank;
      coordinates.push_back( cube.xyz[p] );
      pointIds.push_back( sparseBase + 100 * rank + 10 * side + p );
    }
    for( auto & point : cube.cell.points ) point += offset;
    cells.push_back( cube.cell );
    cellIds.push_back( sparseBase + 10000 + 2 * rank + side );
  }
  auto volume = makeGrid( coordinates, pointIds, cells, cellIds );
  std::vector< Coordinates > surfaceCoordinates;
  Connectivity surfacePointIds;
  vtkNew< vtkIdTypeArray > buckets;
  buckets->SetName( "collocated_nodes" );
  buckets->SetNumberOfComponents( 2 );
  for( int p = 0; p < 4; ++p )
  {
    surfaceCoordinates.push_back( coordinates[sideFaces[0][p]] );
    surfacePointIds.push_back( sparseBase + 20000 + 4 * rank + p );
    vtkIdType const bucket[2] = { pointIds[sideFaces[0][p]], pointIds[8 + sideFaces[1][p]] };
    buckets->InsertNextTypedTuple( bucket );
  }
  auto surface = makeGrid( surfaceCoordinates, surfacePointIds, { { VTK_QUAD, { 0, 1, 2, 3 }, 0 } },
                          { sparseBase + 30000 + rank } );
  surface->GetPointData()->AddArray( buckets );
  for( auto * grid : { volume.GetPointer(), surface.GetPointer() } )
  {
    vtkNew< vtkDoubleArray > selected, unused;
    selected->SetName( "selected" );
    selected->SetNumberOfTuples( grid->GetNumberOfCells() );
    selected->FillValue( 7 + rank );
    grid->GetCellData()->SetScalars( selected );
    unused->SetName( "unusedWideField" );
    unused->SetNumberOfComponents( 128 );
    unused->SetNumberOfTuples( grid->GetNumberOfCells() );
    unused->FillValue( 42 );
    grid->GetCellData()->AddArray( unused );
    vtkNew< vtkStringArray > strings;
    strings->SetName( "unusedLabels" );
    for( vtkIdType c = 0; c < grid->GetNumberOfCells(); ++c ) strings->InsertNextValue( "not imported" );
    grid->GetCellData()->AddArray( strings );
  }
  UniformRefinementOptions options;
  options.requiredPointArrays.emplace();
  options.requiredCellArrays = std::set< std::string >{ "selected", "attribute" };
  options.requiredFaceBlockCellArrays = std::set< std::string >{ "selected" };
  AllMeshes zero( volume, { { "fault", surface } } );
  EXPECT_NO_THROW( refineUniformly( zero, 0, options, MPI_COMM_GEOS ) );
  EXPECT_EQ( zero.getMainMesh().GetPointer(), volume.GetPointer() );
  EXPECT_EQ( zero.getFaceBlocks().at( "fault" ).GetPointer(), surface.GetPointer() );
  EXPECT_NE( volume->GetCellData()->GetArray( "unusedWideField" ), nullptr );
  EXPECT_NE( surface->GetCellData()->GetAbstractArray( "unusedLabels" ), nullptr );
  AllMeshes meshes( volume, { { "fault", surface } } );
  UniformRefinementResult result;
  ASSERT_NO_THROW( result = refineUniformly( meshes, 2, options, MPI_COMM_GEOS ) );
  ASSERT_EQ( result.resources.size(), 2 );
  // Main children retain one double plus an integer region label; surfaces
  // retain one double. The 128-component unused field contributes no growth.
  EXPECT_EQ( result.resources[0].fieldBytesUpperBound, 16 * 12 + 4 * 8 );
  EXPECT_EQ( result.resources[1].fieldBytesUpperBound, 128 * 12 + 16 * 8 );
  EXPECT_EQ( meshes.getMainMesh()->GetNumberOfCells(), 128 );
  EXPECT_EQ( meshes.getFaceBlocks().at( "fault" )->GetNumberOfCells(), 16 );
  EXPECT_NE( meshes.getMainMesh()->GetCellData()->GetArray( "attribute" ), nullptr );
  EXPECT_EQ( meshes.getFaceBlocks().at( "fault" )->GetCellData()->GetArray( "attribute" ), nullptr );
  for( auto const & grid : { meshes.getMainMesh(), meshes.getFaceBlocks().at( "fault" ) } )
  {
    EXPECT_EQ( grid->GetCellData()->GetArray( "unusedWideField" ), nullptr );
    EXPECT_EQ( grid->GetCellData()->GetAbstractArray( "unusedLabels" ), nullptr );
    auto * selected = vtkDoubleArray::SafeDownCast( grid->GetCellData()->GetScalars() );
    ASSERT_NE( selected, nullptr );
    EXPECT_EQ( selected->GetNumberOfTuples(), grid->GetNumberOfCells() );
    for( vtkIdType c = 0; c < grid->GetNumberOfCells(); ++c ) EXPECT_DOUBLE_EQ( selected->GetValue( c ), 7 + rank );
  }
}

TEST( VTKUniformRefinement, CoarseSchemaMismatchFailsBeforeValuesOnlyFineTransfer )
{
  int const size = MpiWrapper::commSize( MPI_COMM_GEOS );
  if( size == 1 )
  {
    GTEST_SKIP() << "Requires a shared coarse interface";
  }
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  auto input = localCube( rank, false );
  if( rank == size - 1 )
  {
    input->GetPointData()->GetScalars()->SetComponentName( 0, "different meaning" );
  }
  AllMeshes meshes( input, {} );
  EXPECT_THROW( refineUniformly( meshes, 2, {}, MPI_COMM_GEOS ), std::runtime_error );
  EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
  EXPECT_EQ( input->GetNumberOfCells(), 1 );
}

TEST( VTKUniformRefinement, EverySupportedEncodingThroughThreeLevels )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  for( char const * filename : { "supportedElements.vtk", "supportedElementsAsVTKPolyhedra.vtk" } )
  {
    vtkNew< vtkDataSetReader > reader;
    reader->SetFileName( ( std::string( VTK_REFINEMENT_FIXTURE_DIR ) + "/" + filename ).c_str() );
    reader->Update();
    auto input = vtkSmartPointer< vtkUnstructuredGrid >::New();
    input->DeepCopy( reader->GetOutput() );
    ASSERT_EQ( input->GetNumberOfCells(), 12 );
    vtkNew< vtkIdTypeArray > pointIds, cellIds;
    for( vtkIdType p = 0; p < input->GetNumberOfPoints(); ++p )
    {
      Coordinates xyz{};
      input->GetPoint( p, xyz.data() );
      xyz[0] += 100 * rank;
      input->GetPoints()->SetPoint( p, xyz.data() );
      pointIds->InsertNextValue( sparseBase + 10000 * rank + p );
    }
    for( vtkIdType c = 0; c < input->GetNumberOfCells(); ++c )
    {
      cellIds->InsertNextValue( sparseBase + 100000 + 100 * rank + c );
    }
    pointIds->SetName( "originalPointIds" );
    cellIds->SetName( "originalCellIds" );
    input->GetPointData()->SetGlobalIds( pointIds );
    input->GetCellData()->SetGlobalIds( cellIds );
    AllMeshes meshes( input, {} );
    auto const result = refineUniformly( meshes, 3, {}, MPI_COMM_GEOS );
    auto output = meshes.getMainMesh();
    EXPECT_EQ( output->GetNumberOfCells(), 10024 );
    std::size_t described = 0;
    for( auto const & block : result.blocks )
    {
      described += block.cells.size();
    }
    EXPECT_EQ( described, 10024 );
    EXPECT_EQ( result.communication.directoryExchanges, 3 );
    PointRegistry points(
        [&]
        {
          std::vector< Coordinates > xyz( output->GetNumberOfPoints() );
          for( vtkIdType p = 0; p < output->GetNumberOfPoints(); ++p )
          {
            output->GetPoint( p, xyz[p].data() );
          }
          return xyz;
        }(),
        [&]
        {
          Connectivity ids;
          auto * active = vtkIdTypeArray::SafeDownCast( output->GetPointData()->GetGlobalIds() );
          for( vtkIdType p = 0; p < output->GetNumberOfPoints(); ++p )
          {
            ids.push_back( active->GetValue( p ) );
          }
          return ids;
        }() );
    for( vtkIdType c = 0; c < output->GetNumberOfCells(); ++c )
    {
      EXPECT_NO_THROW( validateGeometry( normalizeCell( *output->GetCell( c ) ), points ) );
    }
  }
}

TEST( VTKUniformRefinement, ConservativeConnectivityBoundFailsBeforeFineAllocation )
{
  if( sizeof( localIndex ) != 4 )
  {
    GTEST_SKIP() << "This case targets the 32-bit local-index build";
  }
  auto input = localCube( MpiWrapper::commRank( MPI_COMM_GEOS ), false );
  AllMeshes meshes( input, {} );
  // One hex would have 134,217,728 children at level nine. Its cell count
  // fits int32; the conservative face-node incidence bound does not. Exact
  // unique connectivity can be smaller. No huge array is needed to exercise
  // this preflight bound and its collective failure path.
  EXPECT_THROW( refineUniformly( meshes, 9, {}, MPI_COMM_GEOS ), std::runtime_error );
  EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
  EXPECT_EQ( input->GetNumberOfCells(), 1 );
  EXPECT_EQ( input->GetNumberOfPoints(), 8 );
}

TEST( VTKUniformRefinement, PhysicalTransformIsValidatedWithoutChangingVtkCoordinates )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  auto input = localCube( rank, false );
  AllMeshes meshes( input, {} );
  refineUniformly( meshes, 1, {}, MPI_COMM_GEOS );
  auto output = meshes.getMainMesh();
  Coordinates original{};
  output->GetPoint( 0, original.data() );
  EXPECT_NO_THROW( validateRefinedTransform( *output, { 2000, -4, 9 }, { .1, 2, 3 }, MPI_COMM_GEOS ) );
  EXPECT_NO_THROW( validateRefinedTransform( *output, {}, { -1, -1, 1 }, MPI_COMM_GEOS ) );
  EXPECT_THROW( validateRefinedTransform( *output, {}, { -1, 1, 1 }, MPI_COMM_GEOS ), std::runtime_error );
  EXPECT_THROW( validateRefinedTransform( *output, { 1e30, 1e30, 1e30 }, { 1, 1, 1 }, MPI_COMM_GEOS ), std::runtime_error );
  Coordinates after{};
  output->GetPoint( 0, after.data() );
  EXPECT_EQ( original, after );
}

TEST( VTKUniformRefinement, SharedCoordinatesUseTheMeshExtent )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS ), size = MpiWrapper::commSize( MPI_COMM_GEOS );
  auto input = localCube( rank, false );
  for( vtkIdType p = 0; p < input->GetNumberOfPoints(); ++p )
  {
    Coordinates xyz{};
    input->GetPoint( p, xyz.data() );
    for( double & x : xyz )
    {
      x *= 1e-15;
    }
    input->GetPoints()->SetPoint( p, xyz.data() );
  }
  AllMeshes meshes( input, {} );
  EXPECT_NO_THROW( refineUniformly( meshes, 2, {}, MPI_COMM_GEOS ) );
  if( size > 1 )
  {
    // A unit-scale tolerance would silently accept a 0.1% mismatch on
    // this mesh. The rank-one cell remains valid but its shared vertices
    // must disagree collectively before replacing any input dataset.
    if( rank == 1 )
    {
      for( vtkIdType p = 0; p < input->GetNumberOfPoints(); ++p )
      {
        Coordinates xyz{};
        input->GetPoint( p, xyz.data() );
        xyz[1] += 1e-18;
        input->GetPoints()->SetPoint( p, xyz.data() );
      }
    }
    AllMeshes invalid( input, {} );
    EXPECT_THROW( refineUniformly( invalid, 1, {}, MPI_COMM_GEOS ), std::runtime_error );
    EXPECT_EQ( invalid.getMainMesh().GetPointer(), input.GetPointer() );
    EXPECT_EQ( input->GetNumberOfCells(), 1 );
  }
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
