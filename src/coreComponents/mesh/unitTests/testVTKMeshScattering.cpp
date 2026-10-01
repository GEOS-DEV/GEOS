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
 * @file testVTKMeshScattering.cpp
 * @brief Unit tests for VTKMeshScattering (MPI, 4 ranks).
 */

#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"
#include "mesh/generators/VTKMeshScattering.hpp"
#include "mesh/generators/VTKSuperCellPartitioning.hpp"
#include "mesh/generators/VTKUtilities.hpp"

#include <vtkBitArray.h>
#include <vtkStringArray.h>
#include <vtkTypeInt64Array.h>
#include <vtkDummyController.h>
#include <vtkMultiProcessController.h>
#include <vtkExtractCells.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSOADataArrayTemplate.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>

#include <gtest/gtest.h>
#include <set>

using namespace geos;
using namespace geos::vtk;

namespace
{

/**
 * @brief Build a regular 4x4x4 grid on rank 0.
 *
 * The mesh occupies [0, 4] x [0, 4] x [0, 4] with 64 cells,
 * 125 points, a cell data array ("CellId") and a point data array ("PointId").
 * Returns an empty grid on non-zero ranks.
 *
 * @param comm the communicator
 * @param asPolyhedra when true the cells are inserted as VTK_POLYHEDRON with an explicit
 *   6-face description instead of VTK_HEXAHEDRON. The geometry is identical; only the
 *   encoding differs. Polyhedra are the only cell type whose definition is not fully
 *   contained in the connectivity array, so they exercise code paths that all other
 *   cell types leave untested.
 */
vtkSmartPointer< vtkUnstructuredGrid > buildTestMesh( MPI_Comm comm, bool const asPolyhedra = false )
{
  integer const rank = MpiWrapper::commRank( comm );
  auto mesh = vtkSmartPointer< vtkUnstructuredGrid >::New();

  if( rank != 0 )
  {
    return mesh;
  }

  constexpr integer N = 4;  // cells per dimension
  constexpr integer NP = N + 1;

  // Points: (NP)^3 = 125 points
  vtkNew< vtkPoints > points;
  points->SetDataTypeToDouble();
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK( 9, 7, 0 )
  points->Reserve( NP * NP * NP );
#else
  points->Allocate( NP * NP * NP );
#endif

  for( integer k = 0; k <= N; ++k )
  {
    for( integer j = 0; j <= N; ++j )
    {
      for( integer i = 0; i <= N; ++i )
      {
        points->InsertNextPoint( static_cast< real64 >( i ),
                                 static_cast< real64 >( j ),
                                 static_cast< real64 >( k ) );
      }
    }
  }
  mesh->SetPoints( points );

  // Cells: N^3 = 64 cells
  mesh->Allocate( N * N * N );
  for( integer k = 0; k < N; ++k )
  {
    for( integer j = 0; j < N; ++j )
    {
      for( integer i = 0; i < N; ++i )
      {
        vtkIdType const base = i + j * NP + k * NP * NP;
        vtkIdType hex[8] = {
          base,
          base + 1,
          base + NP + 1,
          base + NP,
          base + NP * NP,
          base + NP * NP + 1,
          base + NP * NP + NP + 1,
          base + NP * NP + NP
        };

        if( asPolyhedra )
        {
          // The 6 quad faces of the hexahedron, in terms of the 8 cell points above.
          vtkNew< vtkCellArray > faces;
          auto const insertFace = [&]( integer const a, integer const b,
                                       integer const c, integer const d )
          {
            vtkIdType const face[4] = { hex[a], hex[b], hex[c], hex[d] };
            faces->InsertNextCell( 4, face );
          };
          insertFace( 0, 1, 2, 3 );  // bottom
          insertFace( 4, 5, 6, 7 );  // top
          insertFace( 0, 1, 5, 4 );
          insertFace( 1, 2, 6, 5 );
          insertFace( 2, 3, 7, 6 );
          insertFace( 3, 0, 4, 7 );

          mesh->InsertNextCell( VTK_POLYHEDRON, 8, hex, faces );
        }
        else
        {
          mesh->InsertNextCell( VTK_HEXAHEDRON, 8, hex );
        }
      }
    }
  }

  // Cell data: sequential cell IDs
  vtkNew< vtkIdTypeArray > cellIds;
  cellIds->SetName( "CellId" );
  cellIds->SetNumberOfTuples( N * N * N );
  for( vtkIdType i = 0; i < N * N * N; ++i )
  {
    cellIds->SetValue( i, i );
  }
  mesh->GetCellData()->AddArray( cellIds );

  // Point data: sequential point IDs
  vtkNew< vtkIdTypeArray > pointIds;
  pointIds->SetName( "PointId" );
  pointIds->SetNumberOfTuples( NP * NP * NP );
  for( vtkIdType i = 0; i < NP * NP * NP; ++i )
  {
    pointIds->SetValue( i, i );
  }
  mesh->GetPointData()->AddArray( pointIds );

  return mesh;
}

/**
 * @brief Helper: scatter with a given method and return the result as vtkUnstructuredGrid.
 */
vtkSmartPointer< vtkUnstructuredGrid >
scatter( ScatterMethod method,
         vtkUnstructuredGrid & mesh,
         arrayView1d< integer const > parts,
         MPI_Comm comm )
{
  vtkSmartPointer< vtkDataSet > result = scatterMesh( method, mesh, parts, comm );
  return vtkUnstructuredGrid::SafeDownCast( result );
}

} // anonymous namespace


class VTKMeshScatteringTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    comm = MPI_COMM_GEOS;
    rank = MpiWrapper::commRank( comm );
    size = MpiWrapper::commSize( comm );
    mesh = buildTestMesh( comm );

    parts.resize( 3 );
    parts[0] = 2; parts[1] = 2; parts[2] = 1;
  }

  MPI_Comm comm;
  integer rank;
  integer size;
  vtkSmartPointer< vtkUnstructuredGrid > mesh;
  array1d< integer > parts;

  static constexpr vtkIdType totalCells = 64;
  static constexpr vtkIdType totalPoints = 125;
};


/// All cells must survive the scatter (no loss, no duplication).
TEST_F( VTKMeshScatteringTest, CellConservation )
{
  stdVector< ScatterMethod > methods = {
    ScatterMethod::contiguous,
    ScatterMethod::cartesian,
    ScatterMethod::rcb,
    ScatterMethod::kdtree
  };

  for( auto method : methods )
  {
    auto result = scatter( method, *mesh, parts.toViewConst(), comm );
    ASSERT_NE( result, nullptr );

    vtkIdType const localCells = result->GetNumberOfCells();
    vtkIdType const globalCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Sum, comm );

    EXPECT_EQ( globalCells, totalCells )
      << "Cell conservation failed for method " << toString( method );
  }
}

TEST_F( VTKMeshScatteringTest, DoubleSoAPointsUseTypedTupleFallback )
{
  if( rank == 0 )
  {
    vtkNew< vtkSOADataArrayTemplate< double > > coordinates;
    coordinates->SetNumberOfComponents( 3 );
    coordinates->SetNumberOfTuples( mesh->GetNumberOfPoints() );
    for( vtkIdType p = 0; p < mesh->GetNumberOfPoints(); ++p )
    {
      double point[3];
      mesh->GetPoint( p, point );
      for( int d = 0; d < 3; ++d ) coordinates->SetTypedComponent( p, d, point[d] );
    }
    vtkNew< vtkPoints > points;
    points->SetData( coordinates );
    mesh->SetPoints( points );
  }
  for( auto method : { ScatterMethod::contiguous, ScatterMethod::cartesian, ScatterMethod::rcb } )
  {
    auto result = scatter( method, *mesh, parts.toViewConst(), comm );
    EXPECT_EQ( MpiWrapper::sum( result->GetNumberOfCells(), comm ), totalCells );
    for( vtkIdType p = 0; p < result->GetNumberOfPoints(); ++p )
    {
      double point[3];
      result->GetPoint( p, point );
      for( int d = 0; d < 3; ++d )
      {
        EXPECT_NEAR( point[d], std::round( point[d] ), 1e-12 );
        EXPECT_GE( point[d], 0 );
        EXPECT_LE( point[d], 4 );
      }
    }
  }
}

TEST_F( VTKMeshScatteringTest, BlockFallbackPreservesDistributedInput )
{
  // Feed the already distributed result back into the fallback; it must retain
  // every input rank's cells instead of shipping only the root's portion.
  auto root = scatterByBlock( *mesh, comm );
  EXPECT_EQ( MpiWrapper::sum( root->GetNumberOfCells(), comm ), totalCells );
  auto result = scatterByBlock( *root, comm );
  EXPECT_EQ( MpiWrapper::sum( result->GetNumberOfCells(), comm ), totalCells );
  EXPECT_GT( result->GetNumberOfCells(), 0 );
}


TEST_F( VTKMeshScatteringTest, StringsBitsFieldDataAndRolesSurviveCustomScatter )
{
  if( rank == 0 )
  {
    vtkNew< vtkStringArray > labels;
    labels->SetName( "labels" );
    labels->SetNumberOfComponents( 2 );
    labels->SetComponentName( 0, "material label" );
    labels->SetNumberOfTuples( totalCells );
    vtkNew< vtkBitArray > flags;
    flags->SetName( "flags" );
    flags->SetNumberOfComponents( 2 );
    flags->SetNumberOfTuples( totalCells );
    vtkNew< vtkTypeInt64Array > wide;
    wide->SetName( "wide" );
    wide->SetNumberOfTuples( totalCells );
    for( vtkIdType i = 0; i < totalCells; ++i )
    {
      labels->SetValue( 2*i, std::string( "a\0b", 3 ) + std::to_string( i ) );
      labels->SetValue( 2*i+1, "other" + std::to_string( i ) );
      flags->SetValue( 2*i, i % 2 );
      flags->SetValue( 2*i+1, i % 3 == 0 );
      wide->SetValue( i, INT64_C(9007199254741001) + i );
    }
    mesh->GetCellData()->AddArray( labels );
    vtkNew< vtkStringArray > pedigree;
    pedigree->SetName( "pedigree" );
    for( vtkIdType i = 0; i < totalCells; ++i ) pedigree->InsertNextValue( std::to_string( i ) );
    mesh->GetCellData()->SetPedigreeIds( pedigree );
    mesh->GetCellData()->SetScalars( flags );
    mesh->GetCellData()->AddArray( wide );
    vtkNew< vtkStringArray > pointLabels;
    pointLabels->SetName( "pointLabels" );
    for( vtkIdType i = 0; i < totalPoints; ++i ) pointLabels->InsertNextValue( std::to_string( i ) );
    mesh->GetPointData()->AddArray( pointLabels );
    vtkNew< vtkStringArray > metadata;
    metadata->SetName( "metadata" );
    metadata->InsertNextValue( std::string( "field\0data", 10 ) );
    mesh->GetFieldData()->AddArray( metadata );
  }
  // Exercise both the ordinary and fracture shipping paths with the same metadata.
  if( rank == 0 )
  {
    vtkNew< vtkIdTypeArray > atoms;
    atoms->SetName( "SuperCellId" );
    for( vtkIdType i = 0; i < totalCells; ++i ) atoms->InsertNextValue( i / 2 );
    mesh->GetCellData()->AddArray( atoms );
  }
  for( bool const fractures : { false, true } )
  for( auto method : { ScatterMethod::contiguous, ScatterMethod::cartesian, ScatterMethod::rcb } )
  {
    vtkSmartPointer< vtkDataSet > result;
    if( fractures ) result = redistributeBySuperCellBlocks( mesh, comm, method, parts.toViewConst() );
    else result = scatter( method, *mesh, parts.toViewConst(), comm );
    auto * labels = vtkStringArray::SafeDownCast( result->GetCellData()->GetAbstractArray( "labels" ) );
    auto * flags = vtkBitArray::SafeDownCast( result->GetCellData()->GetArray( "flags" ) );
    auto * wide = vtkTypeInt64Array::SafeDownCast( result->GetCellData()->GetArray( "wide" ) );
    auto * ids = vtkIdTypeArray::SafeDownCast( result->GetCellData()->GetArray( "CellId" ) );
    ASSERT_NE( labels, nullptr ); ASSERT_NE( flags, nullptr ); ASSERT_NE( wide, nullptr ); ASSERT_NE( ids, nullptr );
    auto * pedigree = vtkStringArray::SafeDownCast( result->GetCellData()->GetPedigreeIds() );
    ASSERT_NE( pedigree, nullptr );
    EXPECT_STREQ( pedigree->GetName(), "pedigree" );
    EXPECT_EQ( result->GetCellData()->GetScalars(), flags );
    EXPECT_STREQ( labels->GetComponentName( 0 ), "material label" );
    for( vtkIdType i = 0; i < result->GetNumberOfCells(); ++i )
    {
      auto const id = ids->GetValue( i );
      EXPECT_EQ( pedigree->GetValue( i ), std::to_string( id ) );
      EXPECT_EQ( labels->GetValue( 2*i ), std::string( "a\0b", 3 ) + std::to_string( id ) );
      EXPECT_EQ( labels->GetValue( 2*i+1 ), "other" + std::to_string( id ) );
      EXPECT_EQ( flags->GetValue( 2*i ), id % 2 );
      EXPECT_EQ( flags->GetValue( 2*i+1 ), id % 3 == 0 );
      EXPECT_EQ( wide->GetValue( i ), INT64_C(9007199254741001) + id );
    }
    auto * points = vtkStringArray::SafeDownCast( result->GetPointData()->GetAbstractArray( "pointLabels" ) );
    auto * pointIds = vtkIdTypeArray::SafeDownCast( result->GetPointData()->GetArray( "PointId" ) );
    ASSERT_NE( points, nullptr ); ASSERT_NE( pointIds, nullptr );
    for( vtkIdType i = 0; i < result->GetNumberOfPoints(); ++i ) EXPECT_EQ( points->GetValue( i ), std::to_string( pointIds->GetValue( i ) ) );
    auto * metadata = vtkStringArray::SafeDownCast( result->GetFieldData()->GetAbstractArray( "metadata" ) );
    ASSERT_NE( metadata, nullptr );
    EXPECT_EQ( metadata->GetValue( 0 ), std::string( "field\0data", 10 ) );
  }
}

TEST_F( VTKMeshScatteringTest, KdtreeAcceptsDistributedInputAndRestoresController )
{
  auto full = buildTestMesh( MPI_COMM_SELF );
  vtkNew< vtkIdList > selected;
  // Two input pieces across four ranks; the other ranks start empty.
  if( rank < 2 )
    for( vtkIdType i = 32*rank; i < 32*(rank+1); ++i ) selected->InsertNextId( i );
  vtkNew< vtkExtractCells > extract;
  extract->SetInputData( full );
  extract->SetCellList( selected );
  extract->Update();
  vtkSmartPointer< vtkMultiProcessController > previous = vtkMultiProcessController::GetGlobalController();
  vtkNew< vtkDummyController > expected;
  vtkMultiProcessController::SetGlobalController( expected );
  auto result = scatterMesh( ScatterMethod::kdtree, *extract->GetOutput(), parts.toViewConst(), comm );
  EXPECT_EQ( MpiWrapper::sum( result->GetNumberOfCells(), comm ), totalCells );
  EXPECT_EQ( vtkMultiProcessController::GetGlobalController(), expected.GetPointer() );
  EXPECT_EQ( vtkMultiProcessController::GetGlobalController()->GetNumberOfProcesses(), 1 );
  vtkMultiProcessController::SetGlobalController( previous );
}

TEST_F( VTKMeshScatteringTest, Int64GlobalIdsSurviveSurfaceNeighborDistribution )
{
  vtkIdType const base = INT64_C(9007199254741001);
  if( rank == 0 )
  {
    vtkIdType const boundary[4]{ 0, 1, 6, 5 };
    mesh->InsertNextCell( VTK_QUAD, 4, boundary );
    mesh->GetCellData()->Initialize();
    mesh->GetPointData()->Initialize();
    vtkNew< vtkTypeInt64Array > cells;
    cells->SetName( "GlobalCellIds" );
    for( vtkIdType i = 0; i <= totalCells; ++i ) cells->InsertNextValue( base + i );
    mesh->GetCellData()->SetGlobalIds( cells );
    vtkNew< vtkTypeInt64Array > points;
    points->SetName( "GlobalPointIds" );
    for( vtkIdType i = 0; i < totalPoints; ++i ) points->InsertNextValue( base + 1000 + i );
    mesh->GetPointData()->SetGlobalIds( points );
    // The XML reader promotes active IDs itself; this raw input exercises the
    // importer with Int64 arrays whose actual class is not vtkIdTypeArray.
    EXPECT_EQ( vtkIdTypeArray::SafeDownCast( mesh->GetCellData()->GetGlobalIds() ), nullptr );
    EXPECT_EQ( vtkIdTypeArray::SafeDownCast( mesh->GetPointData()->GetGlobalIds() ), nullptr );
  }
  stdMap< string, vtkSmartPointer< vtkDataSet > > fractures;
  auto distributed = redistributeMeshes( 0, mesh, fractures, comm, ScatterMethod::kdtree,
                                        parts.toViewConst(), PartitionMethod::parmetis, 0, 0, 1, "" );
  auto result = distributed.getMainMesh();
  auto * cells = vtkIdTypeArray::SafeDownCast( result->GetCellData()->GetGlobalIds() );
  auto * points = vtkIdTypeArray::SafeDownCast( result->GetPointData()->GetGlobalIds() );
  ASSERT_NE( cells, nullptr );
  ASSERT_NE( points, nullptr );
  array1d< globalIndex > local( result->GetNumberOfCells() );
  for( vtkIdType i = 0; i < result->GetNumberOfCells(); ++i ) local[i] = cells->GetValue( i );
  array1d< globalIndex > global;
  MpiWrapper::allGatherv( local.toViewConst(), global, comm );
  EXPECT_EQ( global.size(), totalCells + 1 );
  std::set< globalIndex > ids( global.begin(), global.end() );
  for( vtkIdType i = 0; i <= totalCells; ++i ) EXPECT_EQ( ids.count( base + i ), 1 );
  for( vtkIdType i = 0; i < result->GetNumberOfPoints(); ++i )
  {
    EXPECT_GE( points->GetValue( i ), base + 1000 );
    EXPECT_LT( points->GetValue( i ), base + 1000 + totalPoints );
  }
}

/// With 64 cells and 4 ranks, every rank must get exactly 16 cells.
TEST_F( VTKMeshScatteringTest, LoadBalance )
{
  stdVector< ScatterMethod > methods = {
    ScatterMethod::contiguous,
    ScatterMethod::cartesian,
    ScatterMethod::rcb
  };

  vtkIdType const expectedPerRank = totalCells / size;

  for( auto method : methods )
  {
    auto result = scatter( method, *mesh, parts.toViewConst(), comm );
    vtkIdType const localCells = result->GetNumberOfCells();

    EXPECT_EQ( localCells, expectedPerRank )
      << "Rank " << rank << " got " << localCells << " cells for method " << toString( method );
  }
}


/// Cell data and point data arrays must survive the scatter.
TEST_F( VTKMeshScatteringTest, DataArrayPreservation )
{
  auto result = scatter( ScatterMethod::rcb, *mesh, parts.toViewConst(), comm );

  // Cell data
  vtkDataArray * cellIdArr = result->GetCellData()->GetArray( "CellId" );
  ASSERT_NE( cellIdArr, nullptr ) << "CellId array lost during scatter";
  EXPECT_EQ( cellIdArr->GetNumberOfTuples(), result->GetNumberOfCells() );

  // Point data
  vtkDataArray * pointIdArr = result->GetPointData()->GetArray( "PointId" );
  ASSERT_NE( pointIdArr, nullptr ) << "PointId array lost during scatter";
  EXPECT_EQ( pointIdArr->GetNumberOfTuples(), result->GetNumberOfPoints() );
}


/// Cartesian partitioning with a 2x2x1 grid should place cells in the correct spatial quadrant.
TEST_F( VTKMeshScatteringTest, CartesianSpatialCorrectness )
{
  auto result = scatter( ScatterMethod::cartesian, *mesh, parts.toViewConst(), comm );

  // With a 2x2x1 partition on [0,4]^3, the split is at x=2 and y=2.
  // Rank = ix + 2*iy, where ix = (centroid.x < 2 ? 0 : 1), iy = (centroid.y < 2 ? 0 : 1).
  for( vtkIdType c = 0; c < result->GetNumberOfCells(); ++c )
  {
    real64 bounds[6];
    result->GetCell( c )->GetBounds( bounds );
    real64 const cx = ( bounds[0] + bounds[1] ) * 0.5;
    real64 const cy = ( bounds[2] + bounds[3] ) * 0.5;

    integer const expectedIx = ( cx < 2.0 ) ? 0 : 1;
    integer const expectedIy = ( cy < 2.0 ) ? 0 : 1;
    integer const expectedRank = expectedIx + 2 * expectedIy;

    EXPECT_EQ( expectedRank, rank )
      << "Cell centroid (" << cx << ", " << cy << ") on wrong rank";
  }
}


/// Scatter on a single rank (size==1 early exit) returns a deep copy.
TEST_F( VTKMeshScatteringTest, SingleRankNoOp )
{
  // Create a sub-communicator containing only rank 0.
  // On other ranks we just verify we don't crash.
  integer const color = ( rank == 0 ) ? 0 : MPI_UNDEFINED;
  MPI_Comm singleComm = MpiWrapper::commSplit( comm, color, rank );

  if( rank == 0 )
  {
    auto result = scatter( ScatterMethod::rcb, *mesh, parts.toViewConst(), singleComm );
    EXPECT_EQ( result->GetNumberOfCells(), totalCells );
    EXPECT_EQ( result->GetNumberOfPoints(), totalPoints );

    MpiWrapper::commFree( singleComm );
  }
}


/// An empty mesh on all ranks should return an empty grid without error.
TEST_F( VTKMeshScatteringTest, EmptyMesh )
{
  auto emptyMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
  auto result = scatter( ScatterMethod::rcb, *emptyMesh, parts.toViewConst(), comm );

  ASSERT_NE( result, nullptr );
  EXPECT_EQ( result->GetNumberOfCells(), 0 );
}


/// Super-cell atomicity: every super-cell must end up on a single rank, regardless of method.
/// Tags pairs of adjacent cells with the same SuperCellId, so the 64-cell grid becomes
/// 32 atoms of 2 cells each. After redistribution, each SuperCellId must be owned by exactly
/// one rank.
TEST_F( VTKMeshScatteringTest, SuperCellAtomicity )
{
  static constexpr vtkIdType cellsPerSuperCell = 2;
  static constexpr vtkIdType numSuperCells = totalCells / cellsPerSuperCell;

  if( rank == 0 )
  {
    vtkNew< vtkIdTypeArray > scIds;
    scIds->SetName( "SuperCellId" );
    scIds->SetNumberOfTuples( totalCells );
    for( vtkIdType c = 0; c < totalCells; ++c )
    {
      scIds->SetValue( c, c / cellsPerSuperCell );
    }
    mesh->GetCellData()->AddArray( scIds );
  }

  stdVector< ScatterMethod > const methods = {
    ScatterMethod::kdtree,
    ScatterMethod::contiguous,
    ScatterMethod::cartesian,
    ScatterMethod::rcb
  };

  for( auto method : methods )
  {
    vtkSmartPointer< vtkDataSet > result =
      redistributeBySuperCellBlocks( mesh, comm, method, parts.toViewConst() );
    ASSERT_NE( result, nullptr );

    // Cell conservation
    vtkIdType const localCells = result->GetNumberOfCells();
    vtkIdType const globalCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Sum, comm );
    EXPECT_EQ( globalCells, totalCells )
      << "Cell conservation failed for method " << toString( method );

    // Per-rank ownership: ownership[s] = 1 if this rank holds any cell of SuperCellId s.
    array1d< integer > ownership( numSuperCells );
    vtkIdTypeArray * scArr =
      vtkIdTypeArray::SafeDownCast( result->GetCellData()->GetArray( "SuperCellId" ) );
    ASSERT_NE( scArr, nullptr )
      << "SuperCellId array lost during redistribution for method " << toString( method );

    for( vtkIdType c = 0; c < localCells; ++c )
    {
      vtkIdType const s = scArr->GetValue( c );
      ASSERT_GE( s, 0 );
      ASSERT_LT( s, numSuperCells );
      ownership[s] = 1;
      if( method == ScatterMethod::kdtree )
      {
        // Legacy Morton blocks on this regular grid split y and z first.
        // RCB's x/y split gives different owners, so this pins the default.
        vtkIdType const y = (s / 2) % 4;
        vtkIdType const z = s / 8;
        EXPECT_EQ( rank, ( y >= 2 ? 1 : 0 ) + ( z >= 2 ? 2 : 0 ) );
      }
    }

    array1d< integer > totalOwners( numSuperCells );
    MpiWrapper::allReduce( ownership, totalOwners, MpiWrapper::Reduction::Sum, comm );

    for( vtkIdType s = 0; s < numSuperCells; ++s )
    {
      EXPECT_EQ( totalOwners[s], 1 )
        << "Super-cell " << s << " ended up on " << totalOwners[s]
        << " ranks (expected exactly 1) for method " << toString( method );
    }
  }
}


/// Empty rank repair must move whole atoms even when Cartesian or Morton bins are empty.
TEST_F( VTKMeshScatteringTest, FractureScatterFillsEmptyRanksWithWholeSuperCells )
{
  if( rank == 0 )
  {
    vtkNew< vtkIdTypeArray > atoms;
    atoms->SetName( "SuperCellId" );
    for( vtkIdType i = 0; i < totalCells; ++i ) atoms->InsertNextValue( i % 6 );
    mesh->GetCellData()->AddArray( atoms );
  }
  for( auto method : { ScatterMethod::kdtree, ScatterMethod::cartesian, ScatterMethod::rcb } )
  {
    auto result = redistributeBySuperCellBlocks( mesh, comm, method, parts.toViewConst() );
    EXPECT_GT( result->GetNumberOfCells(), 0 );
    EXPECT_EQ( MpiWrapper::sum( result->GetNumberOfCells(), comm ), totalCells );
    auto * ids = vtkIdTypeArray::SafeDownCast( result->GetCellData()->GetArray( "SuperCellId" ) );
    ASSERT_NE( ids, nullptr );
    for( vtkIdType atom = 0; atom < 6; ++atom )
    {
      bool present = false;
      for( vtkIdType i = 0; i < result->GetNumberOfCells(); ++i ) present |= ids->GetValue( i ) == atom;
      auto const minimum = MpiWrapper::min( present ? rank : size, comm );
      auto const maximum = MpiWrapper::max( present ? rank : -1, comm );
      EXPECT_EQ( minimum, maximum ) << "An atomic super-cell was split";
    }
  }
}

TEST_F( VTKMeshScatteringTest, PolyhedronFacePreservation )
{
  auto polyMesh = buildTestMesh( comm, true );

  stdVector< ScatterMethod > const methods = {
    ScatterMethod::contiguous,
    ScatterMethod::cartesian,
    ScatterMethod::rcb
  };

  for( auto method : methods )
  {
    auto result = scatter( method, *polyMesh, parts.toViewConst(), comm );
    ASSERT_NE( result, nullptr );

    vtkIdType const localCells = result->GetNumberOfCells();
    vtkIdType const globalCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Sum, comm );
    EXPECT_EQ( globalCells, totalCells )
      << "Polyhedral cells lost during scatter for method " << toString( method );

    for( vtkIdType c = 0; c < localCells; ++c )
    {
      ASSERT_EQ( result->GetCellType( c ), VTK_POLYHEDRON )
        << "Cell " << c << " changed type for method " << toString( method );

      vtkCell * cell = result->GetCell( c );
      EXPECT_EQ( cell->GetNumberOfPoints(), 8 )
        << "Cell " << c << " lost points for method " << toString( method );
      EXPECT_EQ( cell->GetNumberOfFaces(), 6 )
        << "Cell " << c << " lost its face description for method " << toString( method );
    }
  }
}

TEST_F( VTKMeshScatteringTest, PolyhedronFacePreservationWithSuperCells )
{
  static constexpr vtkIdType cellsPerSuperCell = 2;

  auto polyMesh = buildTestMesh( comm, true );

  if( rank == 0 )
  {
    vtkNew< vtkIdTypeArray > scIds;
    scIds->SetName( "SuperCellId" );
    scIds->SetNumberOfTuples( totalCells );
    for( vtkIdType c = 0; c < totalCells; ++c )
    {
      scIds->SetValue( c, c / cellsPerSuperCell );
    }
    polyMesh->GetCellData()->AddArray( scIds );
  }

  vtkSmartPointer< vtkDataSet > result =
    redistributeBySuperCellBlocks( polyMesh, comm, ScatterMethod::rcb, parts.toViewConst() );
  ASSERT_NE( result, nullptr );

  vtkIdType const localCells = result->GetNumberOfCells();
  vtkIdType const globalCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Sum, comm );
  EXPECT_EQ( globalCells, totalCells ) << "Polyhedral cells lost during super-cell redistribution";

  for( vtkIdType c = 0; c < localCells; ++c )
  {
    ASSERT_EQ( result->GetCellType( c ), VTK_POLYHEDRON ) << "Cell " << c << " changed type";
    EXPECT_EQ( result->GetCell( c )->GetNumberOfFaces(), 6 )
      << "Cell " << c << " lost its face description";
  }
}



int main( int argc, char * * argv )
{
  MpiWrapper::init( &argc, &argv );
  MPI_COMM_GEOS = MpiWrapper::commDup( MPI_COMM_WORLD );

  ::testing::InitGoogleTest( &argc, argv );

  integer const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  if( rank != 0 )
  {
    ::testing::TestEventListeners & listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release( listeners.default_result_printer() );
  }

  integer result = RUN_ALL_TESTS();

  MpiWrapper::finalize();
  return result;
}
