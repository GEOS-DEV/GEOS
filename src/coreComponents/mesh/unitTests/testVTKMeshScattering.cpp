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

#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::vtk;

namespace
{

/**
 * @brief Build a regular 4x4x4 hexahedral grid on rank 0.
 *
 * The mesh occupies [0, 4] x [0, 4] x [0, 4] with 64 hex cells,
 * 125 points, a cell data array ("CellId") and a point data array ("PointId").
 * Returns an empty grid on non-zero ranks.
 */
vtkSmartPointer< vtkUnstructuredGrid > buildTestMesh( MPI_Comm comm )
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
  points->Allocate( NP * NP * NP );

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

  // Hex cells: N^3 = 64 cells
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
        mesh->InsertNextCell( VTK_HEXAHEDRON, 8, hex );
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
