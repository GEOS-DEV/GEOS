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
 * @file testVTKRefinementTemplates.cpp
 */

#include "../generators/VTKRefinementTemplates.hpp"
#include "VTKRefinementTestMeshes.hpp"
#include "../generators/VTKRefinementSharing.hpp"
#include "../generators/VTKRefinementAssociations.hpp"

#include <gtest/gtest.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDataSetReader.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>

using namespace geos::vtk::refinement;

TEST( VTKRefinementTemplates, SurfaceAssociationsUseActualSidesAtJunctions )
{
  vtkIdType const base = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741001 ) ) : 1001;
  std::vector< MainFace > faces;
  std::vector< Connectivity > buckets( 4 );
  for( int side = 0; side < 3; ++side )
  {
    Connectivity corners;
    for( int corner = 0; corner < 4; ++corner )
    {
      vtkIdType const id = base + 10 * side + corner;
      corners.push_back( id );
      buckets[corner].push_back( id );
    }
    faces.push_back( { corners, { side } } );
  }
  // A replicated side adds its actual owner; equal coordinates are irrelevant.
  faces.push_back( { { base + 2, base + 1, base, base + 3 }, { 4 } } );
  // This face has an anchor match but is incomplete, so it creates no side.
  faces.push_back( { { base, base + 100, base + 101, base + 102 }, { 5 } } );
  SurfaceAssociations associations( faces, 7 );
  Connectivity const surface{ 2, 3, 0, 1 };
  auto sides = associations.match( surface, buckets );
  ASSERT_EQ( sides.size(), 3 );
  EXPECT_EQ( sides[0].owners, ( Participants{ 0, 4 } ) );
  for( int side = 0; side < 3; ++side )
  {
    EXPECT_EQ( sides[side].mainCornersBySurface,
               ( Connectivity{ base + 10 * side + 2, base + 10 * side + 3, base + 10 * side, base + 10 * side + 1 } ) );
    EXPECT_EQ( SurfaceAssociations::supportKey( EntityKind::edge, { 3, 2 }, surface, sides[side] ),
               entityKey( EntityKind::edge, { base + 10 * side + 2, base + 10 * side + 3 }, 7 ) );
    EXPECT_EQ( SurfaceAssociations::supportKey( EntityKind::face, surface, surface, sides[side] ), sides[side].mainFace );
    EXPECT_EQ( SurfaceAssociations::supportKey( EntityKind::face, { 0, 1, 2, 3 }, surface, sides[side] ), sides[side].mainFace );
    EXPECT_EQ( SurfaceAssociations::supportKey( EntityKind::vertex, { 0 }, surface, sides[side] ),
               entityKey( EntityKind::vertex, { base + 10 * side }, 7 ) );
    EXPECT_THROW( SurfaceAssociations::supportKey( EntityKind::edge, { 2, 0 }, surface, sides[side] ), std::invalid_argument );
    EXPECT_THROW( SurfaceAssociations::supportKey( EntityKind::face, { 2, 3, 0 }, surface, sides[side] ), std::invalid_argument );
    EXPECT_THROW( SurfaceAssociations::supportKey( EntityKind::cell, surface, surface, sides[side] ), std::invalid_argument );
  }
  EXPECT_EQ( associations.match( { 1, 0, 3, 2 }, buckets ).size(), 3 );
  EXPECT_THROW( associations.match( { 0, 2, 1, 3 }, buckets ), std::invalid_argument );
  auto ambiguous = buckets;
  ambiguous[1].push_back( base );
  EXPECT_THROW( associations.match( surface, ambiguous ), std::invalid_argument );
  auto missing = buckets;
  missing[0] = { base + 1000 };
  EXPECT_TRUE( associations.match( surface, missing ).empty() );
  EXPECT_THROW( associations.match( { 0, 0, 2, 3 }, buckets ), std::invalid_argument );
  EXPECT_THROW( associations.match( { 0, 1, 2, 4 }, buckets ), std::invalid_argument );
  auto crossed = faces;
  crossed.push_back( { { base, base + 2, base + 1, base + 3 }, { 0 } } );
  EXPECT_THROW( ( void )SurfaceAssociations( crossed ), std::invalid_argument );
}

namespace
{
PointRegistry registryFor( vtkDataSet & mesh )
{
  std::vector< Coordinates > xyz( mesh.GetNumberOfPoints() );
  Connectivity gids( xyz.size() );
  for( std::size_t i = 0; i < xyz.size(); ++i )
  {
    mesh.GetPoint( i, xyz[i].data() );
    // Exercise integral topology identity above the double precision limit.
    gids[i] =
      sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741001 ) + 7 * i ) : static_cast< vtkIdType >( 3 + 7 * i );
  }
  return PointRegistry( std::move( xyz ), std::move( gids ) );
}

void verifySubdivision( Cell const & parent, Subdivision const & split, PointRegistry const & registry )
{
  std::map< Connectivity, int > counts, boundary, orientations;
  for( auto const & face : split.faceChildren )
  {
    for( auto const & child : face )
    {
      ++boundary[canonicalCycle( child )];
    }
  }
  double sum = 0;
  for( auto const & child : split.children )
  {
    ASSERT_NO_THROW( validateGeometry( child, registry ) );
    EXPECT_GT( signedMeasure( child, registry ), 0 );
    sum += signedMeasure( child, registry );
    for( auto face : cellFaces( child ) )
    {
      Connectivity const key = canonicalCycle( face );
      ++counts[key];
      std::rotate( face.begin(), std::min_element( face.begin(), face.end() ), face.end() );
      orientations[key] += face == key ? 1 : -1;
    }
  }
  for( auto const & [face, count] : counts )
  {
    EXPECT_EQ( count, boundary.count( face ) ? 1 : 2 ) << "Missing boundary trace or unpaired internal face";
    if( !boundary.count( face ) )
    {
      EXPECT_EQ( orientations[face], 0 ) << "Internal outward orientations disagree";
    }
  }
  for( auto const & [face, count] : boundary )
  {
    EXPECT_EQ( count, 1 );
    EXPECT_EQ( counts[face], 1 );
  }
  double const measure = signedMeasure( parent, registry );
  EXPECT_NEAR( sum, measure, 1e-10 * measure );
}

vtkSmartPointer< vtkUnstructuredGrid > refine( vtkDataSet & mesh )
{
  PointRegistry registry = registryFor( mesh );
  std::vector< Cell > children;
  for( vtkIdType i = 0; i < mesh.GetNumberOfCells(); ++i )
  {
    Cell const parent = normalizeCell( *mesh.GetCell( i ) );
    Subdivision split;
    try
    {
      split = subdivideCell( parent, i, registry );
    }
    catch( std::exception const & error )
    {
      throw std::runtime_error( "Cell " + std::to_string( i ) + ", type " + std::to_string( parent.vtkType ) + ", signed measure " +
                                std::to_string( signedMeasure( parent, registry ) ) + ": " + error.what() );
    }
    verifySubdivision( parent, split, registry );
    children.insert( children.end(), split.children.begin(), split.children.end() );
  }
  auto output = vtkSmartPointer< vtkUnstructuredGrid >::New();
  vtkNew< vtkPoints > points;
  points->SetDataTypeToDouble();
  for( auto const & recipe : registry.points() )
  {
    points->InsertNextPoint( recipe.position.data() );
  }
  output->SetPoints( points );
  for( auto const & child : children )
  {
    output->InsertNextCell( child.vtkType, child.points.size(), child.points.data() );
  }
  return output;
}

using testMeshes::attach;
using testMeshes::referenceCell;
using testMeshes::ReferenceCell;
using testMeshes::regularPrism;

} // namespace

TEST( VTKRefinementTopology, FullKeysCyclesAndNamespaces )
{
  EXPECT_EQ( canonicalCycle( { 10, 20, 30, 40 } ), canonicalCycle( { 30, 20, 10, 40 } ) );
  EXPECT_NE( canonicalCycle( { 10, 20, 30, 40 } ), canonicalCycle( { 10, 30, 20, 40 } ) );
  EXPECT_THROW( canonicalCycle( { 1, 2, 1 } ), std::invalid_argument );
  EntityKey const a{ 0, EntityKind::edge, { 1, 2 } }, b{ 1, EntityKind::edge, { 1, 2 } };
  EXPECT_FALSE( a == b );
  EXPECT_NE( stableHash( a ), stableHash( b ) );
  struct ConstantHash
  {
    std::size_t operator()( EntityKey const & ) const { return 0; }
  };
  std::unordered_map< EntityKey, int, ConstantHash > collisions;
  collisions[a] = 1;
  collisions[b] = 2;
  EXPECT_EQ( collisions.size(), 2 );
  EXPECT_EQ( collisions[a], 1 );
  EXPECT_EQ( collisions[b], 2 );
  PointRegistry registry( { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { 0, 0, 0 } }, { 1, 2, 3, 4, 5 } );
  EXPECT_EQ( registry.edge( 0, 1 ), registry.edge( 1, 0 ) );
  EXPECT_NE( registry.edge( 0, 1 ), registry.edge( 4, 1 ) ); // Coincidence is not identity.
  EXPECT_EQ( registry.face( { 0, 1, 2, 3 } ), registry.face( { 2, 1, 0, 3 } ) );
  EXPECT_THROW( registry.face( { 0, 2, 1, 3 } ), std::invalid_argument );
  EXPECT_THROW( PointRegistry( { { 0, 0, 0 }, { 1, 0, 0 } }, { 1, 1 } ), std::invalid_argument );
}

TEST( VTKRefinementTemplates, ReportedShapesBothEncodingsThreeLevels )
{
  for( auto const * name : { "supportedElements.vtk", "supportedElementsAsVTKPolyhedra.vtk" } )
  {
    SCOPED_TRACE( name );
    vtkNew< vtkDataSetReader > reader;
    reader->SetFileName( ( std::string( VTK_REFINEMENT_FIXTURE_DIR ) + "/" + name ).c_str() );
    reader->Update();
    vtkSmartPointer< vtkDataSet > mesh = reader->GetOutput();
    ASSERT_EQ( mesh->GetNumberOfCells(), 12 );
    vtkIdType const expected[3] = { 154, 1244, 10024 };
    for( int level = 0; level < 3; ++level )
    {
      SCOPED_TRACE( level + 1 );
      ASSERT_NO_THROW( mesh = refine( *mesh ) );
      EXPECT_EQ( mesh->GetNumberOfCells(), expected[level] );
    }
  }
}

TEST( VTKRefinementTemplates, PolygonalPrismsAndAffineRecipes )
{
  for( int n = 5; n <= 11; ++n )
  {
    SCOPED_TRACE( n );
    std::vector< Coordinates > xyz;
    Cell const cell = regularPrism( n, xyz );
    Connectivity gids( xyz.size() );
    std::iota( gids.begin(), gids.end(), 0 );
    PointRegistry registry( xyz, gids );
    auto const split = subdivideCell( cell, 0, registry );
    EXPECT_EQ( split.children.size(), 2 * n );
    EXPECT_EQ( registry.points().size(), 6 * n + 3 );
    EXPECT_EQ( split.faceChildren[0].size(), n );
    EXPECT_EQ( split.faceChildren[1].size(), n );
    verifySubdivision( cell, split, registry );
    for( auto const & recipe : registry.points() )
    {
      double field = 0;
      for( auto const p : recipe.support )
      {
        field += ( 2 * xyz[p][0] - 3 * xyz[p][1] + 4 * xyz[p][2] + 7 ) / recipe.support.size();
      }
      EXPECT_NEAR( field, 2 * recipe.position[0] - 3 * recipe.position[1] + 4 * recipe.position[2] + 7, 1e-13 );
    }
  }
}

TEST( VTKRefinementTemplates, ConnectedMixedTriangleAndQuadInterfacesTwoLevels )
{
  for( int arity : { 3, 4 } )
  {
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
    for( std::size_t i = 0; i < shapes.size(); ++i )
    {
      for( std::size_t j = i; j < shapes.size(); ++j )
      {
        SCOPED_TRACE( std::to_string( arity ) + ":" + std::to_string( i ) + "," + std::to_string( j ) );
        std::vector< Coordinates > xyz = arity == 3 ? std::vector< Coordinates >{ { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } }
                                                    : std::vector< Coordinates >{ { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 } };
        Cell const a = attach( shapes[i], 1, xyz ), b = attach( shapes[j], -1, xyz );
        Connectivity gids( xyz.size() );
        std::iota( gids.begin(), gids.end(), 0 );
        PointRegistry registry( xyz, gids );
        auto const left = subdivideCell( a, 0, registry ), right = subdivideCell( b, 1, registry );
        verifySubdivision( a, left, registry );
        verifySubdivision( b, right, registry );
        std::set< Connectivity > l, r;
        for( auto const & f : left.faceChildren[shapes[i].face] )
        {
          l.insert( canonicalCycle( f ) );
        }
        for( auto const & f : right.faceChildren[shapes[j].face] )
        {
          r.insert( canonicalCycle( f ) );
        }
        EXPECT_EQ( l, r );
        EXPECT_EQ( l.size(), 4 );
        vtkNew< vtkUnstructuredGrid > grid;
        vtkNew< vtkPoints > points;
        points->SetDataTypeToDouble();
        for( auto const & p : registry.points() )
        {
          points->InsertNextPoint( p.position.data() );
        }
        grid->SetPoints( points );
        for( auto const * split : { &left, &right } )
        {
          for( auto const & child : split->children )
          {
            grid->InsertNextCell( child.vtkType, child.points.size(), child.points.data() );
          }
        }
        auto const fine = refine( *grid );
        std::map< Connectivity, int > interfaces;
        for( vtkIdType c = 0; c < fine->GetNumberOfCells(); ++c )
        {
          for( auto const & f : cellFaces( normalizeCell( *fine->GetCell( c ) ) ) )
          {
            bool plane = true;
            for( vtkIdType p : f )
            {
              plane = plane && std::abs( fine->GetPoint( p )[2] ) < 1e-12;
            }
            if( plane )
            {
              ++interfaces[canonicalCycle( f )];
            }
          }
        }
        EXPECT_EQ( interfaces.size(), 16 );
        for( auto const & f : interfaces )
        {
          EXPECT_EQ( f.second, 2 );
        }
      }
    }
  }
}

TEST( VTKRefinementTemplates, NativeAndPolyhedralGeometryEquivalent )
{
  vtkNew< vtkDataSetReader > native, poly;
  native->SetFileName( ( std::string( VTK_REFINEMENT_FIXTURE_DIR ) + "/supportedElements.vtk" ).c_str() );
  native->Update();
  poly->SetFileName( ( std::string( VTK_REFINEMENT_FIXTURE_DIR ) + "/supportedElementsAsVTKPolyhedra.vtk" ).c_str() );
  poly->Update();
  auto signature = []( vtkDataSet & mesh )
  {
    auto output = refine( mesh );
    std::vector< std::vector< std::int64_t > > cells;
    for( vtkIdType i = 0; i < output->GetNumberOfCells(); ++i )
    {
      vtkCell * cell = output->GetCell( i );
      std::vector< std::array< std::int64_t, 3 > > corners;
      for( vtkIdType j = 0; j < cell->GetNumberOfPoints(); ++j )
      {
        auto const p = output->GetPoint( cell->GetPointId( j ) );
        corners.push_back( { std::llround( p[0] * 1e5 ), std::llround( p[1] * 1e5 ), std::llround( p[2] * 1e5 ) } );
      }
      std::sort( corners.begin(), corners.end() );
      std::vector< std::int64_t > key{ cell->GetCellType() };
      for( auto const & p : corners )
      {
        key.insert( key.end(), p.begin(), p.end() );
      }
      cells.push_back( std::move( key ) );
    }
    std::sort( cells.begin(), cells.end() );
    return cells;
  };
  EXPECT_EQ( signature( *native->GetOutput() ), signature( *poly->GetOutput() ) );
}

TEST( VTKRefinementTemplates, PermutedPolyhedralFacesAndMalformedIncidence )
{
  std::vector< ReferenceCell > shapes{ referenceCell( VTK_TETRA, 0 ), referenceCell( VTK_PYRAMID, 0 ), referenceCell( VTK_WEDGE, 0 ),
                                       referenceCell( VTK_HEXAHEDRON, 0 ) };
  for( int n = 5; n <= 11; ++n )
  {
    ReferenceCell r;
    r.cell = regularPrism( n, r.xyz );
    shapes.push_back( std::move( r ) );
  }
  for( auto const & shape : shapes )
  {
    SCOPED_TRACE( shape.cell.vtkType );
    vtkNew< vtkPoints > xyz;
    xyz->SetDataTypeToDouble();
    for( auto const & p : shape.xyz )
    {
      xyz->InsertNextPoint( p.data() );
    }
    auto faces = cellFaces( shape.cell );
    std::reverse( faces.begin(), faces.end() );
    vtkNew< vtkCellArray > stream;
    for( auto face : faces )
    {
      std::reverse( face.begin(), face.end() );
      std::rotate( face.begin(), face.begin() + 1, face.end() );
      stream->InsertNextCell( face.size(), face.data() );
    }
    vtkNew< vtkUnstructuredGrid > poly;
    poly->SetPoints( xyz );
    poly->InsertNextCell( VTK_POLYHEDRON, shape.cell.points.size(), shape.cell.points.data(), stream );
    auto registry = registryFor( *poly );
    Cell const normalized = normalizeCell( *poly->GetCell( 0 ) );
    EXPECT_NO_THROW( subdivideCell( normalized, 0, registry ) );
    stream->InsertNextCell( faces[0].size(),
                            faces[0].data() ); // Duplicate face is never a canonical shape.
    vtkNew< vtkUnstructuredGrid > bad;
    bad->SetPoints( xyz );
    bad->InsertNextCell( VTK_POLYHEDRON, shape.cell.points.size(), shape.cell.points.data(), stream );
    EXPECT_THROW( normalizeCell( *bad->GetCell( 0 ) ), std::invalid_argument );
  }
}

TEST( VTKRefinementTemplates, CoarseFramesIgnoreLocalNumberingAndFaceStreamOrder )
{
  std::vector< ReferenceCell > shapes{ referenceCell( VTK_TETRA, 0 ), referenceCell( VTK_PYRAMID, 0 ), referenceCell( VTK_WEDGE, 0 ),
                                       referenceCell( VTK_HEXAHEDRON, 0 ) };
  for( int n = 5; n <= 11; ++n )
  {
    ReferenceCell r;
    r.cell = regularPrism( n, r.xyz );
    shapes.push_back( std::move( r ) );
  }
  for( auto const & shape : shapes )
  {
    SCOPED_TRACE( std::to_string( shape.cell.vtkType ) + ":" + std::to_string( shape.cell.prismSides ) );
    Connectivity expectedFrame;
    std::vector< std::vector< Coordinates > > expectedChildren[2];
    for( int permutation = 0; permutation < 3; ++permutation )
    {
      for( int encoding = 0; encoding < ( shape.cell.vtkType == VTK_HEXAHEDRON ? 3 : 2 ); ++encoding )
      {
        SCOPED_TRACE( std::to_string( permutation ) + ":" + std::to_string( encoding ) );
        Connectivity order( shape.xyz.size() ), mapped( order.size() );
        std::iota( order.begin(), order.end(), 0 );
        if( permutation == 1 )
        {
          std::reverse( order.begin(), order.end() );
        }
        if( permutation == 2 )
        {
          std::rotate( order.begin(), order.begin() + order.size() / 2, order.end() );
        }
        std::vector< Coordinates > coordinates;
        Connectivity ids;
        vtkIdType const base = sizeof( vtkIdType ) == 8 ? static_cast< vtkIdType >( UINT64_C( 9007199254741019 ) ) : 1001;
        for( std::size_t i = 0; i < order.size(); ++i )
        {
          mapped[order[i]] = i;
          coordinates.push_back( shape.xyz[order[i]] );
          ids.push_back( base + 11 * ( order.size() - 1 - order[i] ) );
        }
        Cell input = shape.cell;
        for( auto & point : input.points )
        {
          point = mapped[point];
        }
        vtkNew< vtkPoints > xyz;
        xyz->SetDataTypeToDouble();
        for( auto const & point : coordinates )
        {
          xyz->InsertNextPoint( point.data() );
        }
        vtkNew< vtkUnstructuredGrid > mesh;
        mesh->SetPoints( xyz );
        int const nativeType = input.prismSides == 5 ? VTK_PENTAGONAL_PRISM : input.prismSides == 6 ? VTK_HEXAGONAL_PRISM : input.vtkType;
        if( encoding == 1 || nativeType == VTK_POLYHEDRON )
        {
          auto faces = cellFaces( input );
          if( encoding == 1 )
          {
            std::reverse( faces.begin(), faces.end() );
          }
          vtkNew< vtkCellArray > stream;
          for( auto face : faces )
          {
            if( encoding == 1 )
            {
              std::reverse( face.begin(), face.end() );
              std::rotate( face.begin(), face.begin() + 1, face.end() );
            }
            stream->InsertNextCell( face.size(), face.data() );
          }
          mesh->InsertNextCell( VTK_POLYHEDRON, input.points.size(), input.points.data(), stream );
        }
        else
        {
          if( encoding == 2 )
          {
            std::swap( input.points[2], input.points[3] );
            std::swap( input.points[6], input.points[7] );
          }
          mesh->InsertNextCell( encoding == 2 ? VTK_VOXEL : nativeType, input.points.size(), input.points.data() );
        }
        PointRegistry registry( coordinates, ids );
        Cell const coarse = normalizeCoarseCell( *mesh->GetCell( 0 ), registry );
        Connectivity frame;
        for( vtkIdType point : coarse.points )
        {
          frame.push_back( ids[point] );
        }
        if( expectedFrame.empty() )
        {
          expectedFrame = frame;
        }
        EXPECT_EQ( frame, expectedFrame );
        auto const split = subdivideCell( coarse, 71, registry );
        verifySubdivision( coarse, split, registry );
        std::vector< Cell > cells = split.children;
        for( int generation = 0; generation < 2; ++generation )
        {
          std::vector< std::vector< Coordinates > > geometry;
          for( auto const & cell : cells )
          {
            std::vector< Coordinates > points;
            for( vtkIdType point : cell.points )
            {
              points.push_back( registry.position( point ) );
            }
            geometry.push_back( std::move( points ) );
          }
          if( expectedChildren[generation].empty() )
          {
            expectedChildren[generation] = geometry;
          }
          ASSERT_EQ( geometry.size(), expectedChildren[generation].size() );
          for( std::size_t c = 0; c < geometry.size(); ++c )
          {
            ASSERT_EQ( geometry[c].size(), expectedChildren[generation][c].size() );
            for( std::size_t p = 0; p < geometry[c].size(); ++p )
            {
              for( int d = 0; d < 3; ++d )
              {
                EXPECT_NEAR( geometry[c][p][d], expectedChildren[generation][c][p][d], 1e-12 );
              }
            }
          }
          if( generation == 0 )
          {
            std::vector< Coordinates > nextCoordinates;
            Connectivity nextIds = ids;
            for( vtkIdType i = 0; i < static_cast< vtkIdType >( registry.points().size() ); ++i )
            {
              nextCoordinates.push_back( registry.position( i ) );
              if( i >= registry.originalSize() )
              {
                nextIds.push_back( base + 11 * order.size() + i );
              }
            }
            registry = PointRegistry( std::move( nextCoordinates ), std::move( nextIds ) );
            std::vector< Cell > next;
            for( std::size_t c = 0; c < cells.size(); ++c )
            {
              // Keep inherited child connectivity; do not rebase it using fine IDs.
              auto const children = subdivideCell( cells[c], c, registry );
              next.insert( next.end(), children.children.begin(), children.children.end() );
            }
            cells = std::move( next );
          }
        }
      }
    }
  }
  vtkNew< vtkPoints > xyz;
  for( auto const & point : referenceCell( VTK_TETRA, 0 ).xyz )
  {
    xyz->InsertNextPoint( point.data() );
  }
  vtkNew< vtkUnstructuredGrid > inverted;
  inverted->SetPoints( xyz );
  Connectivity const wrong{ 0, 2, 1, 3 };
  inverted->InsertNextCell( VTK_TETRA, wrong.size(), wrong.data() );
  auto registry = registryFor( *inverted );
  EXPECT_THROW( subdivideCell( normalizeCoarseCell( *inverted->GetCell( 0 ), registry ), 71, registry ), std::invalid_argument );
}

TEST( VTKRefinementTemplates, PolygonalCapsShareTheSameTrace )
{
  for( int n = 5; n <= 11; ++n )
  {
    std::vector< Coordinates > xyz;
    Cell const lower = regularPrism( n, xyz );
    Cell upper{ VTK_POLYHEDRON, {}, n };
    for( int i = 0; i < n; ++i )
    {
      upper.points.push_back( n + ( i + 2 ) % n );
    }
    for( int i = 0; i < n; ++i )
    {
      auto p = xyz[upper.points[i]];
      p[0] += .3;
      p[1] -= .15;
      p[2] += 1;
      upper.points.push_back( xyz.size() );
      xyz.push_back( p );
    }
    Connectivity ids( xyz.size() );
    std::iota( ids.begin(), ids.end(), 0 );
    PointRegistry registry( xyz, ids );
    auto const a = subdivideCell( lower, 0, registry ), b = subdivideCell( upper, 1, registry );
    std::set< Connectivity > left, right;
    for( auto const & f : a.faceChildren[1] )
    {
      left.insert( canonicalCycle( f ) );
    }
    for( auto const & f : b.faceChildren[0] )
    {
      right.insert( canonicalCycle( f ) );
    }
    EXPECT_EQ( left, right );
    EXPECT_EQ( left.size(), n );
    verifySubdivision( lower, a, registry );
    verifySubdivision( upper, b, registry );
  }
}

TEST( VTKRefinementTemplates, TetrahedronDiagonalIsIndependentOfAllocatedIds )
{
  std::vector< Coordinates > xyz{ { 0, 0, 0 }, { 1, 0, 0 }, { .2, 1, 0 }, { .1, .3, .8 } };
  PointRegistry a( xyz, { 0, 1, 2, 3 } ), b( xyz, { 999, 31, 256, 90 } );
  Cell const parent{ VTK_TETRA, { 0, 1, 2, 3 }, 0 };
  auto const left = subdivideCell( parent, 9, a ), right = subdivideCell( parent, 987, b );
  ASSERT_EQ( left.children.size(), right.children.size() );
  for( std::size_t i = 0; i < left.children.size(); ++i )
  {
    EXPECT_EQ( left.children[i].points, right.children[i].points );
  }
}

TEST( VTKRefinementTemplates, PyramidUnequalFractions )
{
  PointRegistry registry( { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 }, { .5, .5, 1 } }, { 0, 1, 2, 3, 4 } );
  Cell const cell{ VTK_PYRAMID, { 0, 1, 2, 3, 4 }, 0 };
  auto const split = subdivideCell( cell, 0, registry );
  ASSERT_EQ( split.children.size(), 10 );
  EXPECT_EQ( registry.points().size(), 14 );
  for( auto const & child : split.children )
  {
    EXPECT_NEAR( signedMeasure( child, registry ) / signedMeasure( cell, registry ), child.vtkType == VTK_PYRAMID ? 1. / 8 : 1. / 16,
                 1e-14 );
  }
  verifySubdivision( cell, split, registry );
}

TEST( VTKRefinementTemplates, CountsAndOverflowWithoutAllocation )
{
  CellCounts counts;
  counts.tetrahedra = 1;
  counts.pyramids = 1;
  counts.wedges = 1;
  counts.hexahedra = 2;
  for( int n = 5; n <= 11; ++n )
  {
    counts.prisms[n] = 1;
  }
  EXPECT_EQ( counts.total(), 12 );
  for( auto expected : { 154, 1244, 10024 } )
  {
    counts = counts.next();
    EXPECT_EQ( counts.total(), expected );
  }
  counts.hexahedra = std::numeric_limits< std::uint64_t >::max() / 8 + 1;
  EXPECT_THROW( counts.next(), std::overflow_error );
}

TEST( VTKRefinementTemplates, ThinAndSkewedSupportedShapesThroughTwoLevels )
{
  std::vector< ReferenceCell > shapes{ referenceCell( VTK_TETRA, 0 ), referenceCell( VTK_PYRAMID, 0 ),
                                       referenceCell( VTK_WEDGE, 0 ), referenceCell( VTK_HEXAHEDRON, 0 ) };
  for( int n = 5; n <= 11; ++n )
  {
    ReferenceCell prism;
    prism.cell = regularPrism( n, prism.xyz );
    shapes.push_back( std::move( prism ) );
  }
  for( auto const & shape : shapes )
    for( double const thickness : { 1., 1e-4 } )
    {
      SCOPED_TRACE( std::to_string( shape.cell.vtkType ) + ":" + std::to_string( shape.cell.prismSides ) +
                    ":" + std::to_string( thickness ) );
      std::vector< Coordinates > coordinates;
      for( auto const & p : shape.xyz )
        coordinates.push_back( { 3 + .9 * p[0] + .15 * p[1] + .1 * p[2],
                                 -2 + .3 * p[0] + 1.1 * p[1] + .2 * p[2],
                                 7 + thickness * ( 1.2 * p[2] + .05 * p[0] - .1 * p[1] ) } );
      auto mesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
      vtkNew< vtkPoints > points;
      points->SetDataTypeToDouble();
      for( auto const & p : coordinates )
        points->InsertNextPoint( p.data() );
      mesh->SetPoints( points );
      if( shape.cell.prismSides )
      {
        vtkNew< vtkCellArray > faces;
        for( auto const & face : cellFaces( shape.cell ) )
          faces->InsertNextCell( face.size(), face.data() );
        mesh->InsertNextCell( VTK_POLYHEDRON, shape.cell.points.size(), shape.cell.points.data(), faces );
      }
      else
        mesh->InsertNextCell( shape.cell.vtkType, shape.cell.points.size(), shape.cell.points.data() );
      for( int level = 0; level < 2; ++level )
      {
        // Each generation gets a registry for all of its now-existing points,
        // matching the production controller's per-level registry lifecycle.
        mesh = refine( *mesh );
      }
    }
}

TEST( VTKRefinementTemplates, ScaleTranslationAndInvalidGeometry )
{
  for( double scale : { 1e-6, 1., 1e6 } )
  {
    std::vector< Coordinates > xyz{ { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
      { .1, .1, 1 }, { 1.2, 0, 1 }, { 1.1, 1.2, 1 }, { 0, 1, 1 } };
    for( auto & p : xyz )
    {
      for( auto & x : p )
      {
        x = scale * ( x + 100 );
      }
    }
    PointRegistry registry( xyz, { 0, 1, 2, 3, 4, 5, 6, 7 } );
    Cell cell{ VTK_HEXAHEDRON, { 0, 1, 2, 3, 4, 5, 6, 7 }, 0 };
    EXPECT_NO_THROW( subdivideCell( cell, 0, registry ) );
  }
  PointRegistry registry( { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 } }, { 0, 1, 2, 3 } );
  EXPECT_THROW( subdivideCell( { VTK_TETRA, { 0, 1, 2, 3 }, 0 }, 0, registry ), std::invalid_argument );
  PointRegistry flat( { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { .25, .25, 0 } }, { 0, 1, 2, 3 } );
  EXPECT_THROW( subdivideCell( { VTK_TETRA, { 0, 1, 2, 3 }, 0 }, 0, flat ), std::invalid_argument );
  EXPECT_THROW( PointRegistry( { { std::numeric_limits< double >::infinity(), 0, 0 } }, { 0 } ), std::invalid_argument );
  std::vector< Coordinates > xyz;
  Cell star = regularPrism( 5, xyz );
  star.points = { 0, 2, 4, 1, 3, 5, 7, 9, 6, 8 };
  PointRegistry capRegistry( xyz, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } );
  EXPECT_THROW( subdivideCell( star, 0, capRegistry ), std::invalid_argument );
}

TEST( VTKRefinementTemplates, FlatWarpedPyramidNamesTheInvertedChild )
{
  // A valid but flat pyramid with a warped base, taken from a field mesh. The
  // template's center child is inverted, and the error must say so.
  PointRegistry registry( { { 0, 0, 0 }, { -118.14174735, 0.18487936, -10.74911166 },
                            { -110.41312305, -145.77634666, -51.73707128 }, { 3.47049546, -123.45343273, -45.45207629 },
                            { -52.29993373, -67.27097194, -23.70069733 } },
                          { 0, 1, 2, 3, 4 } );
  try
  {
    subdivideCell( { VTK_PYRAMID, { 0, 1, 2, 3, 4 }, 0 }, 545413, registry );
    FAIL() << "The inverted center child must be rejected";
  }
  catch( std::invalid_argument const & error )
  {
    std::string const message = error.what();
    EXPECT_NE( message.find( "the parent pyramid is valid, but child 5 (a pyramid)" ), std::string::npos ) << message;
    EXPECT_NE( message.find( "too flat for its warped base" ), std::string::npos ) << message;
    EXPECT_NE( message.find( "parent 6.91e-04" ), std::string::npos ) << message;
  }
  // An inverted coarse cell is reported as such.
  PointRegistry inverted( { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, -1 } }, { 0, 1, 2, 3 } );
  try
  {
    subdivideCell( { VTK_TETRA, { 0, 1, 2, 3 }, 0 }, 0, inverted );
    FAIL() << "The inverted parent must be rejected";
  }
  catch( std::invalid_argument const & error )
  {
    EXPECT_NE( std::string( error.what() ).find( "the coarse tetrahedron is degenerate or inverted" ), std::string::npos ) << error.what();
  }
}

TEST( VTKRefinementTemplates, SharingFollowsExactSupportRatherThanEndpointOwners )
{
  auto const reference = referenceCell( VTK_HEXAHEDRON, 0 );
  Connectivity const ids{ 71, 81, 91, 101, 111, 121, 131, 141 };
  PointRegistry registry( reference.xyz, ids );
  auto entities = cellEntitySupports( reference.cell, 500, ids, 0 );
  EntityKey const edgeKey = entityKey( EntityKind::edge, { 71, 81 } );
  EntityKey const faceKey = entityKey( EntityKind::face, { 71, 81, 91, 101 } );
  for( auto & entity : entities )
  {
    if( entity.key.kind == EntityKind::vertex )
    {
      entity.participants = { 0, 1, 2 };
    }
    if( entity.key == edgeKey )
    {
      entity.participants = { 0, 2 };
    }
    if( entity.key == faceKey )
    {
      entity.participants = { 0, 1 };
    }
  }
  auto const split = subdivideCell( reference.cell, 500, registry );
  std::reverse( entities.begin(),
                entities.end() ); // Minimal support must not depend on incidence order.
  SharingInheritance inherited( registry, entities, 0 );
  auto const midpoint = registry.edge( 0, 1 ), adjacent = registry.edge( 1, 2 );
  auto const center = registry.cell( 500, reference.cell.points );
  auto const cap = registry.face( { 0, 1, 2, 3 } );
  EXPECT_EQ( inherited.participants( { 0 } ), ( Participants{ 0, 1, 2 } ) );
  EXPECT_EQ( inherited.participants( { midpoint } ), ( Participants{ 0, 2 } ) );
  EXPECT_EQ( inherited.participants( { 0, midpoint } ), ( Participants{ 0, 2 } ) );
  EXPECT_EQ( inherited.participants( { midpoint, adjacent } ), ( Participants{ 0, 1 } ) );
  EXPECT_EQ( inherited.participants( { cap } ), ( Participants{ 0, 1 } ) );
  EXPECT_EQ( inherited.participants( { center, cap } ), ( Participants{ 0 } ) );
  for( Cell const & child : split.children )
  {
    for( auto const & face : cellFaces( child ) )
    {
      EXPECT_NO_THROW( inherited.support( face ) );
    }
  }
  EXPECT_THROW( inherited.support( {} ), std::invalid_argument );
  EXPECT_THROW( inherited.support( { -1 } ), std::invalid_argument );
  auto missing = entities;
  missing.erase( std::remove_if( missing.begin(), missing.end(), [&]( EntitySupport const & e ) { return e.key == edgeKey; } ),
                 missing.end() );
  EXPECT_THROW( SharingInheritance( registry, missing, 0 ), std::invalid_argument );
  auto conflicting = entities;
  conflicting.push_back( entities.front() );
  conflicting.back().participants = { 0, 1 };
  EXPECT_THROW( SharingInheritance( registry, conflicting, 0 ), std::invalid_argument );
  auto replicated = entities;
  for( auto & entity : replicated )
  {
    if( entity.key.kind == EntityKind::cell )
    {
      entity.participants = { 0, 1 };
    }
  }
  EXPECT_THROW( SharingInheritance( registry, replicated, 0 ), std::invalid_argument );
}

TEST( VTKRefinementTemplates, InterfaceSharingRequiresClosedIncidenceAndPlannedPoints )
{
  auto const reference = referenceCell( VTK_HEXAHEDRON, 0 );
  Connectivity const ids{ 71, 81, 91, 101, 111, 121, 131, 141 };
  auto all = cellEntitySupports( reference.cell, 500, ids, 0 );
  std::vector< EntitySupport > shared;
  auto const capKey = entityKey( EntityKind::face, { 71, 81, 91, 101 } );
  for( auto entity : all )
  {
    if( entity.key.kind != EntityKind::cell &&
        std::all_of( entity.localCorners.begin(), entity.localCorners.end(), []( vtkIdType point ) { return point < 4; } ) )
    {
      entity.participants = { 0, 1 };
      shared.push_back( std::move( entity ) );
    }
  }
  ASSERT_EQ( shared.size(), 9 );
  auto missing = shared;
  missing.erase(
    std::remove_if( missing.begin(), missing.end(), []( EntitySupport const & entity ) { return entity.key.kind == EntityKind::edge; } ),
    missing.end() );
  PointRegistry registry( reference.xyz, ids );
  EXPECT_THROW( InterfaceSharing( registry, missing, 0 ), std::invalid_argument );
  InterfaceSharing inherited( registry, shared, 0 );
  EXPECT_THROW( inherited.fineSupports( ids ), std::invalid_argument );
  EXPECT_EQ( registry.points().size(), ids.size() );
  // Plan the complete volume before assigning IDs and installing its traces.
  subdivideCell( reference.cell, 500, registry );
  Connectivity fineIds = ids;
  for( vtkIdType i = registry.originalSize(); i < static_cast< vtkIdType >( registry.points().size() ); ++i )
  {
    fineIds.push_back( 1001 + i );
  }
  auto const fine = inherited.fineSupports( fineIds );
  EXPECT_EQ( fine.size(), 25 ); // 3x3 vertices, 12 edges and four quad faces.
  EXPECT_EQ( inherited.participants( { registry.cell( 500, reference.cell.points ) } ), ( Participants{ 0 } ) );
  EXPECT_EQ( inherited.participants( { registry.face( { 0, 1, 2, 3 } ) } ), ( Participants{ 0, 1 } ) );
  for( auto const & entity : fine )
  {
    EXPECT_EQ( entity.participants, ( Participants{ 0, 1 } ) );
  }
  auto inconsistent = shared;
  for( auto & entity : inconsistent )
  {
    if( entity.key == capKey )
    {
      entity.participants = { 0, 1, 2 };
    }
  }
  EXPECT_THROW( InterfaceSharing( registry, inconsistent, 0 ), std::invalid_argument );
  EXPECT_THROW( inherited.participants( {} ), std::invalid_argument );
  EXPECT_THROW( inherited.participants( { -1 } ), std::invalid_argument );
  auto duplicateCells = std::vector< Cell >{ reference.cell, reference.cell };
  EXPECT_THROW( coarseBoundary( duplicateCells, { 500, 501 }, ids, 0 ), std::invalid_argument );
  EXPECT_THROW( coarseBoundary( duplicateCells, { 500, 500 }, ids, 0 ), std::invalid_argument );
}

TEST( VTKRefinementTemplates, LocalIncidenceDeduplicatesSharedEntitiesAndRejectsDuplicateVolumes )
{
  auto const reference = referenceCell( VTK_HEXAHEDRON, 0 );
  Connectivity ids{ 71, 81, 91, 101, 111, 121, 131, 141 };
  PointRegistry registry( reference.xyz, ids );
  auto const split = subdivideCell( reference.cell, 500, registry );
  for( vtkIdType i = registry.originalSize(); i < static_cast< vtkIdType >( registry.points().size() ); ++i )
  {
    ids.push_back( 1001 + i );
  }
  Connectivity cells( split.children.size() );
  std::iota( cells.begin(), cells.end(), 4001 );
  auto const actual = meshEntitySupports( split.children, cells, ids, 3, 19 );
  std::map< EntityKey, Connectivity > expected;
  for( std::size_t i = 0; i < split.children.size(); ++i )
  {
    for( auto const & entity : cellEntitySupports( split.children[i], cells[i], ids, 3, 19 ) )
    {
      auto corners = entity.localCorners;
      std::sort( corners.begin(), corners.end() );
      expected.emplace( entity.key, std::move( corners ) );
    }
  }
  ASSERT_EQ( actual.size(), expected.size() );
  EXPECT_EQ( actual.size(), 125 ); // Complete incidence of a 2x2x2 hex grid.
  for( auto const & entity : actual )
  {
    auto corners = entity.localCorners;
    std::sort( corners.begin(), corners.end() );
    EXPECT_EQ( corners, expected.at( entity.key ) );
    EXPECT_EQ( entity.participants, ( Participants{ 3 } ) );
    EXPECT_EQ( entity.key.meshNamespace, 19 );
  }
  cells.back() = cells.front();
  EXPECT_THROW( meshEntitySupports( split.children, cells, ids, 3 ), std::invalid_argument );
  cells.pop_back();
  EXPECT_THROW( meshEntitySupports( split.children, cells, ids, 3 ), std::invalid_argument );
  EXPECT_TRUE( meshEntitySupports( {}, {}, {}, 3 ).empty() );
}

TEST( VTKRefinementTemplates, DatasetMeshWhenRequested )
{
  char const * path = std::getenv( "GEOS_REFINEMENT_TEST_MESH" );
  if( !path )
  {
    GTEST_SKIP() << "Set GEOS_REFINEMENT_TEST_MESH to exercise an external VTU mesh";
  }
  vtkNew< vtkXMLUnstructuredGridReader > reader;
  reader->SetFileName( path );
  reader->Update();
  ASSERT_GT( reader->GetOutput()->GetNumberOfCells(), 0 );
  vtkSmartPointer< vtkUnstructuredGrid > refined;
  ASSERT_NO_THROW( refined = refine( *reader->GetOutput() ) );
  RecordProperty( "inputCells", std::to_string( reader->GetOutput()->GetNumberOfCells() ) );
  RecordProperty( "refinedCells", std::to_string( refined->GetNumberOfCells() ) );
  CellCounts counts;
  for( vtkIdType i = 0; i < reader->GetOutput()->GetNumberOfCells(); ++i )
  {
    Cell const cell = normalizeCell( *reader->GetOutput()->GetCell( i ) );
    if( cell.prismSides )
    {
      ++counts.prisms[cell.prismSides];
    }
    else
    {
      switch( cell.vtkType )
      {
        case VTK_TETRA:
          ++counts.tetrahedra;
          break;
        case VTK_PYRAMID:
          ++counts.pyramids;
          break;
        case VTK_WEDGE:
          ++counts.wedges;
          break;
        case VTK_HEXAHEDRON:
          ++counts.hexahedra;
          break;
        default:
          FAIL() << "Unexpected normalized volume cell";
      }
    }
  }
  EXPECT_EQ( static_cast< std::uint64_t >( refined->GetNumberOfCells() ), counts.next().total() );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
