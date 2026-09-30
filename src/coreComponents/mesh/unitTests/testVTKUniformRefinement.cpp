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

#include <algorithm>
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
  if( size > 1 )
  {
    EXPECT_FALSE( result.neighbors.empty() );
  }
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

TEST( VTKUniformRefinement, FlattenedConnectivityOverflowFailsBeforeFineAllocation )
{
  if( sizeof( localIndex ) != 4 )
  {
    GTEST_SKIP() << "This case targets the 32-bit local-index build";
  }
  auto input = localCube( MpiWrapper::commRank( MPI_COMM_GEOS ), false );
  AllMeshes meshes( input, {} );
  // One hex would have 134,217,728 children at level nine. Its cell count
  // fits int32; flattened face connectivity does not. No huge array is needed
  // to validate this preflight and its collective failure path.
  EXPECT_THROW( refineUniformly( meshes, 9, {}, MPI_COMM_GEOS ), std::runtime_error );
  EXPECT_EQ( meshes.getMainMesh().GetPointer(), input.GetPointer() );
  EXPECT_EQ( input->GetNumberOfCells(), 1 );
  EXPECT_EQ( input->GetNumberOfPoints(), 8 );
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
