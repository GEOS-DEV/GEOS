/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testVTKHybridRedistribution.cpp
 */

#include "common/MpiWrapper.hpp"
#include "mesh/generators/VTKUtilities.hpp"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkMultiProcessController.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <gtest/gtest.h>

namespace geos
{
namespace vtk
{
namespace
{

vtkSmartPointer< vtkUnstructuredGrid > buildRootHexMesh( MPI_Comm const comm )
{
  vtkNew< vtkUnstructuredGrid > mesh;
  if( MpiWrapper::commRank( comm ) != 0 )
  {
    return mesh;
  }

  constexpr int N = 4;
  constexpr int NP = N + 1;
  vtkNew< vtkPoints > points;
  points->SetDataTypeToDouble();
  for( int k = 0; k <= N; ++k )
  {
    for( int j = 0; j <= N; ++j )
    {
      for( int i = 0; i <= N; ++i )
      {
        points->InsertNextPoint( i, j, k );
      }
    }
  }
  mesh->SetPoints( points );
  mesh->Allocate( N * N * N );
  for( int k = 0; k < N; ++k )
  {
    for( int j = 0; j < N; ++j )
    {
      for( int i = 0; i < N; ++i )
      {
        vtkIdType const base = i + j * NP + k * NP * NP;
        vtkIdType const hex[8] = {
          base, base + 1, base + NP + 1, base + NP,
          base + NP * NP, base + NP * NP + 1,
          base + NP * NP + NP + 1, base + NP * NP + NP
        };
        mesh->InsertNextCell( VTK_HEXAHEDRON, 8, hex );
      }
    }
  }

  vtkNew< vtkIdTypeArray > pointIds;
  pointIds->SetName( "GlobalPointIds" );
  pointIds->SetNumberOfTuples( mesh->GetNumberOfPoints() );
  for( vtkIdType point = 0; point < mesh->GetNumberOfPoints(); ++point )
  {
    pointIds->SetValue( point, point );
  }
  mesh->GetPointData()->SetGlobalIds( pointIds );

  vtkNew< vtkIdTypeArray > cellIds;
  cellIds->SetName( "GlobalCellIds" );
  cellIds->SetNumberOfTuples( mesh->GetNumberOfCells() );
  vtkNew< vtkDoubleArray > cellValue;
  cellValue->SetName( "cellValue" );
  cellValue->SetNumberOfTuples( mesh->GetNumberOfCells() );
  for( vtkIdType cell = 0; cell < mesh->GetNumberOfCells(); ++cell )
  {
    cellIds->SetValue( cell, cell );
    cellValue->SetValue( cell, 0.5 + cell );
  }
  mesh->GetCellData()->SetGlobalIds( cellIds );
  mesh->GetCellData()->AddArray( cellValue );
  return mesh;
}

TEST( VTKHybridRedistribution, DirectRootPreservesCellsFieldsAndExactNeighbors )
{
  MPI_Comm const comm = MPI_COMM_GEOS;
  ASSERT_EQ( MpiWrapper::commSize( comm ), 4 );
  vtkSmartPointer< vtkMultiProcessController > controller = getController();
  vtkMultiProcessController::SetGlobalController( controller );

  vtkSmartPointer< vtkUnstructuredGrid > mesh = buildRootHexMesh( comm );
  stdMap< string, vtkSmartPointer< vtkDataSet > > fractures;
  array1d< int > cartesianPartitions( 3 );
  cartesianPartitions[0] = 2;
  cartesianPartitions[1] = 2;
  cartesianPartitions[2] = 1;
  HybridPartitionOptions options;
  options.refinementPasses = 1;

  AllMeshes result = redistributeMeshes( 0,
                                         mesh,
                                         fractures,
                                         comm,
                                         ScatterMethod::rcb,
                                         cartesianPartitions.toViewConst(),
                                         PartitionMethod::parmetis,
                                         1,
                                         0,
                                         0,
                                         "",
                                         PartitionModel::hybrid,
                                         options );

  vtkSmartPointer< vtkDataSet > localMesh = result.getMainMesh();
  ASSERT_NE( localMesh, nullptr );
  EXPECT_GT( localMesh->GetNumberOfCells(), 0 );
  EXPECT_EQ( MpiWrapper::sum( localMesh->GetNumberOfCells(), comm ), 64 );
  EXPECT_TRUE( result.hasExactNeighborRanks() );
  EXPECT_FALSE( result.getExactNeighborRanks().empty() );

  vtkDataArray * const globalIds = localMesh->GetCellData()->GetGlobalIds();
  vtkDataArray * const values = localMesh->GetCellData()->GetArray( "cellValue" );
  ASSERT_NE( globalIds, nullptr );
  ASSERT_NE( values, nullptr );
  EXPECT_EQ( globalIds->GetNumberOfTuples(), localMesh->GetNumberOfCells() );
  EXPECT_EQ( values->GetNumberOfTuples(), localMesh->GetNumberOfCells() );

  int const rank = MpiWrapper::commRank( comm );
  int localMask = 0;
  for( int const neighbor : result.getExactNeighborRanks() )
  {
    EXPECT_NE( neighbor, rank );
    localMask |= 1 << neighbor;
  }
  array1d< int > masks;
  MpiWrapper::allGather( localMask, masks, comm );
  for( int neighbor = 0; neighbor < MpiWrapper::commSize( comm ); ++neighbor )
  {
    if( (localMask & (1 << neighbor)) != 0 )
    {
      EXPECT_NE( masks[neighbor] & (1 << rank), 0 );
    }
  }
}

} // namespace
} // namespace vtk
} // namespace geos

int main( int argc, char * * argv )
{
  geos::MpiWrapper::init( &argc, &argv );
  geos::MPI_COMM_GEOS = geos::MpiWrapper::commDup( MPI_COMM_WORLD );
  ::testing::InitGoogleTest( &argc, argv );

  if( geos::MpiWrapper::commRank( geos::MPI_COMM_GEOS ) != 0 )
  {
    ::testing::TestEventListeners & listeners = ::testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release( listeners.default_result_printer() );
  }

  int const result = RUN_ALL_TESTS();
  geos::MpiWrapper::finalize();
  return result;
}
