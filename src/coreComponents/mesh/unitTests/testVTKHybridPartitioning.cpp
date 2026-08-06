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
 * @file testVTKHybridPartitioning.cpp
 */

#include "mesh/generators/VTKHybridPartitioning.hpp"

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
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

vtkSmartPointer< vtkUnstructuredGrid > buildMixedMesh()
{
  vtkNew< vtkUnstructuredGrid > mesh;
  vtkNew< vtkPoints > points;
  points->SetDataTypeToDouble();
  for( vtkIdType point = 0; point < 13; ++point )
  {
    points->InsertNextPoint( static_cast< double >( point % 4 ),
                             static_cast< double >( (point / 4) % 2 ),
                             static_cast< double >( point / 8 ) );
  }
  mesh->SetPoints( points );
  mesh->Allocate( 4 );

  vtkIdType const hex[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
  vtkIdType const pyramid[5] = { 4, 5, 6, 7, 8 };
  vtkIdType const tetrahedron[4] = { 4, 5, 8, 9 };
  vtkIdType const wedge[6] = { 4, 5, 9, 10, 11, 12 };
  mesh->InsertNextCell( VTK_HEXAHEDRON, 8, hex );
  mesh->InsertNextCell( VTK_PYRAMID, 5, pyramid );
  mesh->InsertNextCell( VTK_TETRA, 4, tetrahedron );
  mesh->InsertNextCell( VTK_WEDGE, 6, wedge );

  vtkNew< vtkIdTypeArray > pointIds;
  pointIds->SetName( "GlobalPointIds" );
  pointIds->SetNumberOfTuples( mesh->GetNumberOfPoints() );
  for( vtkIdType point = 0; point < mesh->GetNumberOfPoints(); ++point )
  {
    pointIds->SetValue( point, 100 + point );
  }
  mesh->GetPointData()->SetGlobalIds( pointIds );

  vtkNew< vtkIdTypeArray > cellIds;
  cellIds->SetName( "GlobalCellIds" );
  cellIds->SetNumberOfTuples( mesh->GetNumberOfCells() );
  for( vtkIdType cell = 0; cell < mesh->GetNumberOfCells(); ++cell )
  {
    cellIds->SetValue( cell, 200 + cell );
  }
  mesh->GetCellData()->SetGlobalIds( cellIds );
  return mesh;
}

TEST( VTKHybridPartitioning, BuildsExactMixedTopologyAndSymmetricWeightedGraph )
{
  vtkSmartPointer< vtkUnstructuredGrid > mesh = buildMixedMesh();
  HybridPartitionOptions options;
  HybridPartitionTopology const topology = buildHybridPartitionTopology( *mesh, nullptr, options );

  EXPECT_EQ( topology.cellToVertex.size(), 4 );
  EXPECT_EQ( topology.vertexToCells.size(), 4 );
  EXPECT_EQ( topology.faceToVertices.size( 0 ), 3 );
  EXPECT_EQ( topology.vertexWeights.size( 0 ), 4 );
  EXPECT_EQ( topology.vertexWeights.size( 1 ), 3 );

  auto const graph = topology.graph.toViewConst();
  EXPECT_GT( graph.getOffsets()[graph.size()], 0 );
  EXPECT_EQ( topology.edgeWeights.size(), graph.getOffsets()[graph.size()] );
  for( pmet_idx_t vertex = 0; vertex < graph.size(); ++vertex )
  {
    for( pmet_idx_t const neighbor : graph[vertex] )
    {
      EXPECT_NE( vertex, neighbor );
      EXPECT_NE( std::find( graph[neighbor].begin(), graph[neighbor].end(), vertex ),
                 graph[neighbor].end() );
    }
  }
}

TEST( VTKHybridPartitioning, HonorsWeightFieldsAndAtomicSuperCells )
{
  vtkSmartPointer< vtkUnstructuredGrid > mesh = buildMixedMesh();

  vtkNew< vtkDoubleArray > fvmWeights;
  fvmWeights->SetName( "fvm_cost" );
  fvmWeights->SetNumberOfTuples( 4 );
  for( vtkIdType cell = 0; cell < 4; ++cell )
  {
    fvmWeights->SetValue( cell, cell + 1.0 );
  }
  mesh->GetCellData()->AddArray( fvmWeights );

  vtkNew< vtkIdTypeArray > superCellIds;
  superCellIds->SetName( "SuperCellId" );
  superCellIds->SetNumberOfTuples( 4 );
  superCellIds->SetValue( 0, 10 );
  superCellIds->SetValue( 1, 20 );
  superCellIds->SetValue( 2, 20 );
  superCellIds->SetValue( 3, 30 );
  mesh->GetCellData()->AddArray( superCellIds );

  HybridPartitionOptions options;
  options.fvmWeightField = "fvm_cost";
  options.imbalance.setValues< serialPolicy >( 1.0 );
  options.refinementPasses = 2;

  HybridPartitionTopology const topology = buildHybridPartitionTopology( *mesh, nullptr, options );
  EXPECT_EQ( topology.vertexToCells.size(), 3 );
  EXPECT_EQ( topology.cellToVertex[1], topology.cellToVertex[2] );
  EXPECT_GT( topology.vertexWeights( topology.cellToVertex[2], 0 ),
             topology.vertexWeights( topology.cellToVertex[0], 0 ) );

  HybridPartitionResult const result = partitionHybridMeshOnRoot( *mesh, nullptr, options, 2 );
  EXPECT_EQ( result.cellParts.size(), 4 );
  EXPECT_EQ( result.cellParts[1], result.cellParts[2] );
  EXPECT_LE( result.metrics.finalObjective, result.metrics.initialObjective + 1.0e-10 );
  EXPECT_EQ( result.rankNeighbors.size(), 2 );
  EXPECT_FALSE( result.rankNeighbors[0].empty() );
  EXPECT_FALSE( result.rankNeighbors[1].empty() );
}

TEST( VTKHybridPartitioning, ReportsUnsupportedTopologyForFallback )
{
  vtkSmartPointer< vtkUnstructuredGrid > mesh = buildMixedMesh();
  vtkIdType const line[2] = { 0, 1 };
  mesh->InsertNextCell( VTK_LINE, 2, line );

  string reason;
  EXPECT_FALSE( isHybridPartitioningSupported( *mesh, reason ) );
  EXPECT_NE( reason.find( "unsupported VTK type" ), string::npos );
}

} // namespace
} // namespace vtk
} // namespace geos
