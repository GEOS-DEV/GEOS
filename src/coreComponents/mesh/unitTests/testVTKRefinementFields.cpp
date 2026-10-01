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

/** @file testVTKRefinementFields.cpp */
#include "../generators/VTKRefinementFields.hpp"
#include "../generators/VTKRefinementTemplates.hpp"

#include <gtest/gtest.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkBitArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>
#include <vtkTypeInt64Array.h>
#include <vtkUnsignedCharArray.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

using namespace geos::vtk::refinement;

namespace
{
std::vector< Coordinates > const coordinates{ { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { .5, .5, 1 } };
Cell const pyramid{ VTK_PYRAMID, { 0, 1, 2, 3, 4 }, 0 };
vtkTypeInt64 const largeLabel = INT64_C( 9007199254741019 );

vtkSmartPointer< vtkPointData > pointData()
{
  auto data = vtkSmartPointer< vtkPointData >::New();
  vtkNew< vtkDoubleArray > vector;
  vector->SetName( "velocity" );
  vector->SetNumberOfComponents( 3 );
  vector->SetNumberOfTuples( 5 );
  vector->SetComponentName( 0, "x" );
  vector->SetComponentName( 1, "y" );
  vector->SetComponentName( 2, "z" );
  vtkNew< vtkFloatArray > tensor;
  tensor->SetName( "tensor" );
  tensor->SetNumberOfComponents( 9 );
  tensor->SetNumberOfTuples( 5 );
  vtkNew< vtkTypeInt64Array > labels;
  labels->SetName( "label" );
  labels->SetNumberOfComponents( 2 );
  labels->SetNumberOfTuples( 5 );
  vtkNew< vtkUnsignedCharArray > base;
  base->SetName( "base" );
  base->SetNumberOfTuples( 5 );
  vtkNew< vtkUnsignedCharArray > corner;
  corner->SetName( "corner" );
  corner->SetNumberOfTuples( 5 );
  vtkNew< vtkIdTypeArray > ids;
  ids->SetName( "originalIds" );
  ids->SetNumberOfTuples( 5 );
  vtkNew< vtkUnsignedCharArray > ghosts;
  ghosts->SetName( vtkDataSetAttributes::GhostArrayName() );
  ghosts->SetNumberOfTuples( 5 );
  for( vtkIdType i = 0; i < 5; ++i )
  {
    for( int c = 0; c < 3; ++c )
      vector->SetTypedComponent( i, c, c + coordinates[i][0] + 2 * coordinates[i][1] - coordinates[i][2] );
    for( int c = 0; c < 9; ++c )
      tensor->SetTypedComponent( i, c, static_cast< float >( c + coordinates[i][0] ) );
    labels->SetTypedComponent( i, 0, largeLabel );
    labels->SetTypedComponent( i, 1, -largeLabel );
    base->SetValue( i, i < 4 ? 1 : 0 );
    corner->SetValue( i, i == 0 ? 1 : 0 );
    ids->SetValue( i, 500 + i );
    ghosts->SetValue( i, 0 );
  }
  data->SetVectors( vector );
  data->SetTensors( tensor );
  data->AddArray( labels );
  data->AddArray( base );
  data->SetActiveScalars( "base" );
  data->AddArray( corner );
  data->SetGlobalIds( ids );
  data->AddArray( ghosts );
  return data;
}

TransferPolicies nodeSetPolicies()
{
  TransferPolicies policies;
  policies.pointArrays = { { "base", PointTransferPolicy::nodeSet }, { "corner", PointTransferPolicy::nodeSet } };
  return policies;
}
} // namespace

TEST( VTKRefinementFields, TypedAffineFieldsAndNodeSets )
{
  auto input = pointData();
  PointRegistry points( coordinates, { 500, 501, 502, 503, 504 } );
  subdivideCell( pyramid, 900, points );
  auto const center = points.cell( 900, pyramid.points );
  auto output = transferPointData( *input, points, nodeSetPolicies() );
  EXPECT_EQ( output->GetGlobalIds(), nullptr );
  EXPECT_EQ( output->GetArray( "originalIds" ), nullptr );
  EXPECT_EQ( output->GetArray( vtkDataSetAttributes::GhostArrayName() ), nullptr );
  EXPECT_STREQ( output->GetVectors()->GetName(), "velocity" );
  EXPECT_STREQ( output->GetScalars()->GetName(), "base" );
  EXPECT_STREQ( output->GetTensors()->GetName(), "tensor" );
  EXPECT_STREQ( output->GetVectors()->GetComponentName( 1 ), "y" );
  auto * labels = vtkTypeInt64Array::SafeDownCast( output->GetArray( "label" ) );
  ASSERT_NE( labels, nullptr );
  EXPECT_EQ( labels->GetNumberOfComponents(), 2 );
  for( vtkIdType i = 0; i < static_cast< vtkIdType >( points.points().size() ); ++i )
  {
    auto const & position = points.position( i );
    for( int c = 0; c < 3; ++c )
      EXPECT_NEAR( output->GetVectors()->GetComponent( i, c ), c + position[0] + 2 * position[1] - position[2], 1e-14 );
    for( int c = 0; c < 9; ++c )
      EXPECT_NEAR( output->GetTensors()->GetComponent( i, c ), c + position[0], 1e-6 );
    EXPECT_EQ( labels->GetTypedComponent( i, 0 ), largeLabel );
    EXPECT_EQ( labels->GetTypedComponent( i, 1 ), -largeLabel );
    if( i >= points.originalSize() )
    {
      EXPECT_DOUBLE_EQ( output->GetArray( "corner" )->GetComponent( i, 0 ), 0 );
    }
  }
  EXPECT_DOUBLE_EQ( output->GetArray( "base" )->GetComponent( points.face( { 0, 1, 2, 3 } ), 0 ), 1 );
  EXPECT_DOUBLE_EQ( output->GetArray( "base" )->GetComponent( points.edge( 0, 4 ), 0 ), 0 );
  EXPECT_DOUBLE_EQ( output->GetArray( "base" )->GetComponent( center, 0 ), 0 );
  // Existing input values remain byte-exact, including signed zero.
  vtkDoubleArray::SafeDownCast( input->GetVectors() )->SetTypedComponent( 0, 0, -0. );
  output = transferPointData( *input, points, nodeSetPolicies() );
  EXPECT_TRUE( std::signbit( output->GetVectors()->GetComponent( 0, 0 ) ) );
}

TEST( VTKRefinementFields, IntensiveAndMeasuredExtensivePyramidFields )
{
  vtkNew< vtkCellData > input;
  vtkNew< vtkDoubleArray > extensive;
  extensive->SetName( "inventory" );
  extensive->SetNumberOfComponents( 3 );
  extensive->SetNumberOfTuples( 1 );
  extensive->SetTypedComponent( 0, 0, 32 );
  extensive->SetTypedComponent( 0, 1, 64 );
  extensive->SetTypedComponent( 0, 2, 96 );
  vtkNew< vtkFloatArray > intensive;
  intensive->SetName( "state" );
  intensive->SetNumberOfComponents( 9 );
  intensive->SetNumberOfTuples( 1 );
  for( int c = 0; c < 9; ++c )
    intensive->SetTypedComponent( 0, c, static_cast< float >( c + 1 ) );
  vtkNew< vtkTypeInt64Array > label;
  label->SetName( "material" );
  label->InsertNextValue( largeLabel );
  vtkNew< vtkStringArray > text;
  text->SetName( "description" );
  text->InsertNextValue( "parent material" );
  vtkNew< vtkIdTypeArray > ids;
  ids->SetName( "oldCells" );
  ids->InsertNextValue( 12345 );
  input->SetVectors( extensive );
  input->SetTensors( intensive );
  input->AddArray( label );
  input->SetActiveScalars( "material" );
  input->AddArray( text );
  input->SetGlobalIds( ids );
  PointRegistry points( coordinates, { 500, 501, 502, 503, 504 } );
  auto const split = subdivideCell( pyramid, 900, points );
  std::vector< double > fractions;
  for( auto const & child : split.children )
    fractions.push_back( signedMeasure( child, points ) / signedMeasure( pyramid, points ) );
  TransferPolicies policies;
  policies.extensiveCellArrays = { "inventory" };
  auto output = transferCellData( *input, 1, Connectivity( split.children.size(), 0 ), fractions, policies );
  auto * copied = vtkTypeInt64Array::SafeDownCast( output->GetArray( "material" ) );
  EXPECT_STREQ( output->GetScalars()->GetName(), "material" );
  ASSERT_NE( copied, nullptr );
  double total[3]{};
  for( vtkIdType i = 0; i < static_cast< vtkIdType >( split.children.size() ); ++i )
  {
    EXPECT_EQ( copied->GetValue( i ), largeLabel );
    EXPECT_EQ( vtkStringArray::SafeDownCast( output->GetAbstractArray( "description" ) )->GetValue( i ), "parent material" );
    for( int c = 0; c < 3; ++c )
    {
      double const quantity = output->GetVectors()->GetComponent( i, c );
      total[c] += quantity;
      EXPECT_NEAR( quantity, 32 * ( c + 1 ) * ( split.children[i].vtkType == VTK_PYRAMID ? 1. / 8 : 1. / 16 ), 1e-14 );
    }
    for( int c = 0; c < 9; ++c )
      EXPECT_DOUBLE_EQ( output->GetTensors()->GetComponent( i, c ), c + 1 );
  }
  for( int c = 0; c < 3; ++c )
    EXPECT_NEAR( total[c], 32 * ( c + 1 ), 1e-12 );
  EXPECT_EQ( output->GetGlobalIds(), nullptr );
  EXPECT_EQ( output->GetArray( "oldCells" ), nullptr );
}

TEST( VTKRefinementFields, CellFractionsCheckEveryParentInArbitraryOrder )
{
  vtkNew< vtkCellData > input;
  vtkNew< vtkDoubleArray > quantity;
  quantity->SetName( "quantity" );
  quantity->InsertNextValue( 3 );
  quantity->InsertNextValue( 7 );
  input->AddArray( quantity );
  TransferPolicies policies;
  policies.extensiveCellArrays.insert( "quantity" );
  Connectivity const parents{ 1, 0, 1, 0 };
  std::vector< double > const fractions{ 0.75, 0.25, 0.25, 0.75 };
  auto output = transferCellData( *input, 2, parents, fractions, policies );
  std::vector< double > const expected{ 5.25, 0.75, 1.75, 2.25 };
  for( vtkIdType i = 0; i < 4; ++i )
  {
    EXPECT_DOUBLE_EQ( output->GetArray( "quantity" )->GetComponent( i, 0 ), expected[i] );
  }
  EXPECT_THROW( transferCellData( *input, 2, { 0 }, { 1 }, policies ), std::invalid_argument );
  EXPECT_THROW( transferCellData( *input, 2, parents, { 0.5, 0.25, 0.25, 0.75 }, policies ), std::invalid_argument );
}

TEST( VTKRefinementFields, CanonicalTuplesPreserveExactIntegersAndCheckSchema )
{
  auto input = pointData();
  PointRegistry points( coordinates, { 500, 501, 502, 503, 504 } );
  subdivideCell( pyramid, 900, points );
  auto data = transferPointData( *input, points, nodeSetPolicies() );
  PointFieldLayout original( *data );
  auto tuple = original.pack( 7 );
  vtkNew< vtkPointData > reordered;
  for( int i = data->GetNumberOfArrays() - 1; i >= 0; --i )
  {
    vtkSmartPointer< vtkAbstractArray > copy;
    copy.TakeReference( data->GetAbstractArray( i )->NewInstance() );
    copy->DeepCopy( data->GetAbstractArray( i ) );
    reordered->AddArray( copy );
  }
  reordered->SetActiveVectors( "velocity" );
  reordered->SetActiveScalars( "base" );
  reordered->SetActiveTensors( "tensor" );
  PointFieldLayout destination( *reordered );
  auto * label = vtkTypeInt64Array::SafeDownCast( reordered->GetArray( "label" ) );
  label->SetTypedComponent( 7, 0, 0 );
  label->SetTypedComponent( 7, 1, 0 );
  destination.install( 7, tuple );
  EXPECT_EQ( label->GetTypedComponent( 7, 0 ), largeLabel );
  EXPECT_EQ( label->GetTypedComponent( 7, 1 ), -largeLabel );
  EXPECT_EQ( destination.pack( 7 ), tuple );
  auto const values = original.pack( 7, FieldTupleFormat::valuesOnly );
  EXPECT_LT( values.size(), tuple.size() );
  label->SetTypedComponent( 7, 0, 0 );
  destination.install( 7, values, FieldTupleFormat::valuesOnly );
  EXPECT_EQ( destination.pack( 7 ), tuple );
  auto shortValues = values;
  shortValues.pop_back();
  EXPECT_THROW( destination.install( 7, shortValues, FieldTupleFormat::valuesOnly ), std::invalid_argument );
  EXPECT_EQ( destination.pack( 7 ), tuple );
  auto truncated = tuple;
  truncated.pop_back();
  EXPECT_THROW( destination.install( 7, truncated ), std::invalid_argument );
  tuple.push_back( 0 );
  EXPECT_THROW( destination.install( 7, tuple ), std::invalid_argument );
  EXPECT_THROW( original.pack( -1 ), std::invalid_argument );
  EXPECT_THROW( original.pack( points.points().size() ), std::invalid_argument );
  reordered->GetVectors()->SetComponentName( 0, "changed semantic" );
  PointFieldLayout changed( *reordered );
  EXPECT_THROW( changed.install( 7, original.pack( 7 ) ), std::invalid_argument );
}

TEST( VTKRefinementFields, InvalidOrUnimplementedPoliciesFail )
{
  auto input = pointData();
  PointRegistry points( coordinates, { 500, 501, 502, 503, 504 } );
  subdivideCell( pyramid, 900, points );
  auto policies = nodeSetPolicies();
  vtkTypeInt64Array::SafeDownCast( input->GetArray( "label" ) )->SetTypedComponent( 0, 0, largeLabel + 1 );
  EXPECT_THROW( transferPointData( *input, points, policies ), std::invalid_argument );
  vtkTypeInt64Array::SafeDownCast( input->GetArray( "label" ) )->SetTypedComponent( 0, 0, largeLabel );
  policies.pointArrays["label"] = PointTransferPolicy::continuous;
  EXPECT_THROW( transferPointData( *input, points, policies ), std::invalid_argument );
  policies = nodeSetPolicies();
  policies.pointArrays["missing"] = PointTransferPolicy::equal;
  EXPECT_THROW( transferPointData( *input, points, policies ), std::invalid_argument );
  policies = nodeSetPolicies();
  policies.pointArrays["velocity"] = PointTransferPolicy::nodeSet;
  EXPECT_THROW( transferPointData( *input, points, policies ), std::invalid_argument );
  policies = nodeSetPolicies();
  vtkDoubleArray::SafeDownCast( input->GetVectors() )->SetTypedComponent( 0, 0, std::numeric_limits< double >::infinity() );
  EXPECT_THROW( transferPointData( *input, points, policies ), std::invalid_argument );
  vtkNew< vtkCellData > cell;
  vtkNew< vtkTypeInt64Array > counts;
  counts->SetName( "integer" );
  counts->InsertNextValue( largeLabel );
  cell->AddArray( counts );
  policies.extensiveCellArrays = { "integer" };
  EXPECT_THROW( transferCellData( *cell, 1, { 0, 0 }, { .5, .5 }, policies ), std::invalid_argument );
  policies.extensiveCellArrays = { "missing" };
  EXPECT_THROW( transferCellData( *cell, 1, { 0, 0 }, { .5, .5 }, policies ), std::invalid_argument );
  policies.extensiveCellArrays.clear();
  EXPECT_THROW( transferCellData( *cell, 1, { 0, 0 }, { .4, .4 }, policies ), std::invalid_argument );
  EXPECT_THROW( transferCellData( *cell, 1, { 0, 1 }, { .5, .5 }, policies ), std::invalid_argument );
  EXPECT_THROW( transferCellData( *cell, 1, {}, {}, policies ), std::invalid_argument );
}

TEST( VTKRefinementFields, RelationalArraysRequireTheirSpecializedHandler )
{
  auto input = pointData();
  vtkNew< vtkIdTypeArray > associations;
  associations->SetName( "collocated_nodes" );
  associations->SetNumberOfComponents( 2 );
  associations->SetNumberOfTuples( 5 );
  for( vtkIdType i = 0; i < 5; ++i )
  {
    associations->SetTypedComponent( i, 0, 500 + i );
    associations->SetTypedComponent( i, 1, 1500 + i );
  }
  input->AddArray( associations );
  PointRegistry points( coordinates, { 500, 501, 502, 503, 504 } );
  subdivideCell( pyramid, 900, points );
  auto output = transferPointData( *input, points, nodeSetPolicies() );
  EXPECT_EQ( output->GetAbstractArray( "collocated_nodes" ), nullptr );
  auto policies = nodeSetPolicies();
  policies.pointArrays["collocated_nodes"] = PointTransferPolicy::equal;
  EXPECT_THROW( transferPointData( *input, points, policies ), std::invalid_argument );
}

TEST( VTKRefinementFields, SurfaceCellTupleCodecPreservesIntegersStringsBitsAndRoles )
{
  vtkNew< vtkCellData > source;
  vtkNew< vtkTypeInt64Array > category;
  category->SetName( "category" );
  category->InsertNextValue( largeLabel );
  source->AddArray( category );
  vtkNew< vtkDoubleArray > vector;
  vector->SetName( "vector" );
  vector->SetNumberOfComponents( 3 );
  vector->SetComponentName( 0, "x" );
  vector->SetNumberOfTuples( 1 );
  vector->SetTypedComponent( 0, 0, -0.0 );
  vector->SetTypedComponent( 0, 1, 1.25 );
  vector->SetTypedComponent( 0, 2, -3.5 );
  source->SetVectors( vector );
  vtkNew< vtkStringArray > labels;
  labels->SetName( "labels" );
  labels->SetNumberOfComponents( 2 );
  labels->SetNumberOfTuples( 1 );
  labels->SetValue( 0, std::string( "a\0b", 3 ) );
  labels->SetValue( 1, "fault" );
  source->AddArray( labels );
  vtkNew< vtkBitArray > flags;
  flags->SetName( "flags" );
  flags->SetNumberOfComponents( 2 );
  flags->SetNumberOfTuples( 1 );
  flags->SetValue( 0, 1 );
  flags->SetValue( 1, 0 );
  source->AddArray( flags );
  vtkNew< vtkCellData > target;
  target->DeepCopy( source );
  vtkTypeInt64Array::SafeDownCast( target->GetArray( "category" ) )->SetValue( 0, 2 );
  vtkStringArray::SafeDownCast( target->GetAbstractArray( "labels" ) )->SetValue( 0, "changed" );
  vtkBitArray::SafeDownCast( target->GetArray( "flags" ) )->SetValue( 0, 0 );
  CellFieldLayout from( *source ), to( *target );
  auto const tuple = from.pack( 0 );
  EXPECT_NO_THROW( to.install( 0, tuple ) );
  EXPECT_EQ( to.pack( 0 ), tuple );
  EXPECT_EQ( vtkTypeInt64Array::SafeDownCast( target->GetArray( "category" ) )->GetValue( 0 ), largeLabel );
  EXPECT_EQ( vtkStringArray::SafeDownCast( target->GetAbstractArray( "labels" ) )->GetValue( 0 ), std::string( "a\0b", 3 ) );
  EXPECT_TRUE( std::signbit( target->GetVectors()->GetComponent( 0, 0 ) ) );
  auto const values = from.pack( 0, FieldTupleFormat::valuesOnly );
  EXPECT_LT( values.size(), tuple.size() );
  vtkTypeInt64Array::SafeDownCast( target->GetArray( "category" ) )->SetValue( 0, 0 );
  vtkStringArray::SafeDownCast( target->GetAbstractArray( "labels" ) )->SetValue( 0, "changed again" );
  vtkBitArray::SafeDownCast( target->GetArray( "flags" ) )->SetValue( 0, 0 );
  to.install( 0, values, FieldTupleFormat::valuesOnly );
  EXPECT_EQ( to.pack( 0 ), tuple );
  auto shortValues = values;
  shortValues.pop_back();
  EXPECT_THROW( to.install( 0, shortValues, FieldTupleFormat::valuesOnly ), std::invalid_argument );
  EXPECT_EQ( to.pack( 0 ), tuple );
  auto malformed = tuple;
  malformed.pop_back();
  EXPECT_THROW( to.install( 0, malformed ), std::invalid_argument );
  EXPECT_EQ( to.pack( 0 ), tuple );
  malformed = tuple;
  malformed.push_back( 0 );
  EXPECT_THROW( to.install( 0, malformed ), std::invalid_argument );
  EXPECT_EQ( to.pack( 0 ), tuple );
  target->GetArray( "category" )->SetName( "different" );
  CellFieldLayout wrong( *target );
  EXPECT_THROW( wrong.install( 0, tuple ), std::invalid_argument );
}

int main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
