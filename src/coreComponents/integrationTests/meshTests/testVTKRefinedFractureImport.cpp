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

/** @file testVTKRefinedFractureImport.cpp */
#include "LvArray/src/system.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlockABC.hpp"
#include "mesh/generators/FaceBlockABC.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include "mesh/generators/VTKUniformRefinement.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include <gtest/gtest.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkExtractCells.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkTypeInt64Array.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

using namespace geos;
using namespace geos::dataRepository;

namespace
{
using Point = std::array< double, 3 >;
using Corners = std::vector< vtkIdType >;
globalIndex const idBase = sizeof( globalIndex ) == 8 && sizeof( vtkIdType ) == 8 ? INT64_C( 9007199254741001 ) : 10001;
CommandLineOptions commandLineOptions;

struct Fixture
{
  vtkSmartPointer< vtkUnstructuredGrid > main = vtkSmartPointer< vtkUnstructuredGrid >::New();
  std::map< string, vtkSmartPointer< vtkUnstructuredGrid > > blocks;

  Fixture()
  {
    vtkNew< vtkPoints > points;
    points->SetDataTypeToDouble();
    main->SetPoints( points );
  }

  void volume( int type, std::vector< Point > const & coordinates, std::vector< Corners > const & faces = {} )
  {
    Corners points;
    for( auto const & point : coordinates ) points.push_back( main->GetPoints()->InsertNextPoint( point.data() ) );
    if( faces.empty() ) main->InsertNextCell( type, points.size(), points.data() );
    else
    {
      vtkNew< vtkCellArray > polyFaces;
      for( auto const & face : faces )
      {
        Corners mapped;
        for( auto p : face ) mapped.push_back( points.at( p ) );
        polyFaces->InsertNextCell( mapped.size(), mapped.data() );
      }
      main->InsertNextCell( VTK_POLYHEDRON, points.size(), points.data(), polyFaces );
    }
  }

  void hex( double x, double y )
  {
    volume( VTK_HEXAHEDRON, { { x, y, 0 }, { x + 1, y, 0 }, { x + 1, y + 1, 0 }, { x, y + 1, 0 },
                             { x, y, 1 }, { x + 1, y, 1 }, { x + 1, y + 1, 1 }, { x, y + 1, 1 } } );
  }

  void ids( vtkUnstructuredGrid & grid )
  {
    vtkNew< vtkTypeInt64Array > points, cells;
    points->SetName( "GlobalPointIds" );
    cells->SetName( "GlobalCellIds" );
    for( vtkIdType p = 0; p < grid.GetNumberOfPoints(); ++p ) points->InsertNextValue( idBase + 17 * p );
    for( vtkIdType c = 0; c < grid.GetNumberOfCells(); ++c ) cells->InsertNextValue( idBase + 20000 + 13 * c );
    grid.GetPointData()->SetGlobalIds( points );
    grid.GetCellData()->SetGlobalIds( cells );
  }

  void surface( string const & name, int type, std::vector< std::vector< Point > > const & cells, double gap = 0 )
  {
    auto grid = vtkSmartPointer< vtkUnstructuredGrid >::New();
    vtkNew< vtkPoints > points;
    points->SetDataTypeToDouble();
    grid->SetPoints( points );
    std::map< Point, vtkIdType > indices;
    // This fixture intentionally identifies the specified junction vertices.
    // Main-volume points are always distinct, including coincident side copies.
    for( auto const & cell : cells )
    {
      Corners corners;
      for( auto const & point : cell )
      {
        auto const [it, inserted] = indices.emplace( point, indices.size() );
        if( inserted ) points->InsertNextPoint( point.data() );
        corners.push_back( it->second );
      }
      grid->InsertNextCell( type, corners.size(), corners.data() );
    }
    ids( *grid );
    std::vector< Corners > buckets( points->GetNumberOfPoints() );
    std::size_t width = 0;
    auto * mainIds = vtkTypeInt64Array::SafeDownCast( main->GetPointData()->GetGlobalIds() );
    for( auto const & [point, index] : indices )
    {
      for( vtkIdType p = 0; p < main->GetNumberOfPoints(); ++p )
      {
        Point position{};
        main->GetPoint( p, position.data() );
        bool const match = position == point || ( gap > 0 && std::equal_to< double >{}( position[0], point[0] ) &&
          std::equal_to< double >{}( position[1], point[1] ) && std::equal_to< double >{}( std::abs( position[2] - point[2] ), gap ) );
        if( match ) buckets[index].push_back( mainIds->GetValue( p ) );
      }
      if( buckets[index].size() < 2 ) throw std::runtime_error( "Fixture surface has fewer than two actual sides" );
      width = std::max( width, buckets[index].size() );
    }
    vtkNew< vtkIdTypeArray > collocation;
    collocation->SetName( "collocated_nodes" );
    collocation->SetNumberOfComponents( LvArray::integerConversion< int >( width + 1 ) );
    collocation->SetNumberOfTuples( buckets.size() );
    collocation->FillValue( -1 );
    for( std::size_t p = 0; p < buckets.size(); ++p )
      for( std::size_t b = 0; b < buckets[p].size(); ++b ) collocation->SetTypedComponent( p, b, buckets[p][b] );
    grid->GetPointData()->AddArray( collocation );
    blocks.emplace( name, grid );
  }

  void prepareMain()
  {
    // One contact component plus these disconnected atoms fills every coarse
    // rank without splitting the contact, as required by the normal scatter.
    for( int r = 1; r < MpiWrapper::commSize(); ++r ) hex( 10 + 2 * r, 5 );
    ids( *main );
  }
};

Fixture triangle( double gap )
{
  Fixture result;
  result.volume( VTK_TETRA, { { 0, 0, gap }, { 1, 0, gap }, { 0, 1, gap }, { 0, 0, 1 + gap } } );
  result.volume( VTK_TETRA, { { 0, 0, -gap }, { 0, 1, -gap }, { 1, 0, -gap }, { 0, 0, -1 - gap } } );
  result.prepareMain();
  result.surface( "triangle", VTK_TRIANGLE, { { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 } } }, gap );
  return result;
}

Fixture cap( int n )
{
  Fixture result;
  std::vector< Point > surface;
  for( int i = 0; i < n; ++i )
  {
    double const angle = 2 * std::acos( -1. ) * i / n;
    surface.push_back( { std::cos( angle ), std::sin( angle ), 0 } );
  }
  for( int side : { -1, 0 } )
  {
    std::vector< Point > coordinates;
    for( int z : { side, side + 1 } )
      for( auto point : surface )
      {
        point[2] = z;
        coordinates.push_back( point );
      }
    Corners bottom, top;
    for( int i = 0; i < n; ++i )
    {
      bottom.push_back( n - 1 - i );
      top.push_back( n + i );
    }
    std::vector< Corners > faces{ bottom, top };
    for( int i = 0; i < n; ++i ) faces.push_back( { i, ( i + 1 ) % n, n + ( i + 1 ) % n, n + i } );
    result.volume( VTK_POLYHEDRON, coordinates, faces );
  }
  result.prepareMain();
  result.surface( "cap", VTK_POLYGON, { surface } );
  return result;
}

Fixture junction()
{
  Fixture result;
  for( int x : { -1, 0 } )
    for( int y : { -1, 0 } ) result.hex( x, y );
  result.prepareMain();
  result.surface( "faultX", VTK_QUAD,
    { { { 0, -1, 0 }, { 0, 0, 0 }, { 0, 0, 1 }, { 0, -1, 1 } },
      { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 0, 1 } } } );
  result.surface( "faultY", VTK_QUAD,
    { { { -1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 }, { -1, 0, 1 } },
      { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 0, 1 }, { 0, 0, 1 } } } );
  return result;
}

Fixture strip()
{
  Fixture result;
  // Two disconnected contact atoms, joined by ordinary volume faces along y=0.
  // Each side reuses its own volume nodes; opposite contact sides remain distinct.
  for( int side = 0; side < 2; ++side )
  {
    vtkIdType const offset = result.main->GetNumberOfPoints();
    for( int z = 0; z < 2; ++z )
      for( int y = 0; y < 3; ++y )
        for( int x = 0; x < 2; ++x ) result.main->GetPoints()->InsertNextPoint( side + x - 1, y - 1, z );
    auto index = [offset]( int x, int y, int z ) { return offset + 6 * z + 2 * y + x; };
    for( int y = 0; y < 2; ++y )
    {
      Corners const corners{ index( 0, y, 0 ), index( 1, y, 0 ), index( 1, y + 1, 0 ), index( 0, y + 1, 0 ),
                             index( 0, y, 1 ), index( 1, y, 1 ), index( 1, y + 1, 1 ), index( 0, y + 1, 1 ) };
      result.main->InsertNextCell( VTK_HEXAHEDRON, corners.size(), corners.data() );
    }
  }
  // Exactly one indivisible atom per rank forces the two contact atoms apart.
  for( int r = 2; r < MpiWrapper::commSize(); ++r ) result.hex( 10 + 2 * r, 5 );
  result.ids( *result.main );
  result.surface( "strip", VTK_QUAD,
    { { { 0, -1, 0 }, { 0, 0, 0 }, { 0, 0, 1 }, { 0, -1, 1 } },
      { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 0, 1 } } } );
  return result;
}

class FixtureFiles
{
public:
  template< typename FACTORY > explicit FixtureFiles( FACTORY const & factory, bool parallelMain = false )
  {
    string error;
    if( MpiWrapper::commRank() == 0 )
    {
      try
      {
        LvArray::system::FloatingPointExceptionGuard guard;
        auto fixture = factory();
        auto const stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        m_folder = ( std::filesystem::temp_directory_path() / ( "tmp-geos-refined-fracture-" + std::to_string( stamp ) ) ).string();
        if( !std::filesystem::create_directory( m_folder ) ) throw std::runtime_error( "Fixture folder already exists" );
        if( parallelMain )
        {
          for( int piece = 0; piece < 2; ++piece )
          {
            vtkNew< vtkIdList > ids;
            for( vtkIdType c = piece; c < fixture.main->GetNumberOfCells(); c += 2 ) ids->InsertNextId( c );
            vtkNew< vtkExtractCells > extractor;
            extractor->SetInputData( fixture.main );
            extractor->SetCellList( ids );
            extractor->Update();
            vtkNew< vtkXMLUnstructuredGridWriter > writer;
            writer->SetFileName( ( m_folder + "/piece" + std::to_string( piece ) + ".vtu" ).c_str() );
            writer->SetInputData( extractor->GetOutput() );
            writer->SetDataModeToAscii();
            if( writer->Write() != 1 ) throw std::runtime_error( "Could not write parallel main piece" );
          }
          std::ofstream summary( m_folder + "/main.pvtu" );
          summary << R"xml(<?xml version="1.0"?>
<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian">
<PUnstructuredGrid GhostLevel="0">
<PPointData GlobalIds="GlobalPointIds"><PDataArray type="Int64" Name="GlobalPointIds"/></PPointData>
<PCellData GlobalIds="GlobalCellIds"><PDataArray type="Int64" Name="GlobalCellIds"/></PCellData>
<PPoints><PDataArray type="Float64" NumberOfComponents="3"/></PPoints>
<Piece Source="piece0.vtu"/><Piece Source="piece1.vtu"/>
</PUnstructuredGrid></VTKFile>)xml";
          if( !summary.good() ) throw std::runtime_error( "Could not write parallel main summary" );
          // vtkXMLMultiBlockDataReader does not read PVTU references. Keep the
          // XML import fixture serial; the direct importer below reads main.pvtu.
          vtkNew< vtkXMLUnstructuredGridWriter > mainWriter;
          mainWriter->SetFileName( ( m_folder + "/main.vtu" ).c_str() );
          mainWriter->SetInputData( fixture.main );
          mainWriter->SetDataModeToAscii();
          if( mainWriter->Write() != 1 ) throw std::runtime_error( "Could not write main XML fixture" );
          std::ofstream bundle( path() );
          bundle << R"xml(<?xml version="1.0"?><VTKFile type="vtkMultiBlockDataSet" version="1.0">
<vtkMultiBlockDataSet><DataSet index="0" name="main" file="main.vtu"/>)xml";
          int index = 1;
          for( auto const & [name, grid] : fixture.blocks )
          {
            vtkNew< vtkXMLUnstructuredGridWriter > writer;
            writer->SetFileName( ( m_folder + "/" + name + ".vtu" ).c_str() );
            writer->SetInputData( grid );
            writer->SetDataModeToAscii();
            if( writer->Write() != 1 ) throw std::runtime_error( "Could not write fracture block" );
            bundle << "<DataSet index=\"" << index++ << "\" name=\"" << name << "\" file=\"" << name << ".vtu\"/>";
          }
          bundle << "</vtkMultiBlockDataSet></VTKFile>";
          if( !bundle.good() ) throw std::runtime_error( "Could not write parallel mesh bundle" );
        }
        else
        {
          vtkNew< vtkMultiBlockDataSet > blocks;
          blocks->SetNumberOfBlocks( fixture.blocks.size() + 1 );
          blocks->SetBlock( 0, fixture.main );
          blocks->GetMetaData( 0U )->Set( blocks->NAME(), "main" );
          unsigned int index = 1;
          for( auto const & [name, grid] : fixture.blocks )
          {
            blocks->SetBlock( index, grid );
            blocks->GetMetaData( index++ )->Set( blocks->NAME(), name.c_str() );
          }
          vtkNew< vtkXMLMultiBlockDataWriter > writer;
          writer->SetFileName( path().c_str() );
          writer->SetInputData( blocks );
          writer->SetDataModeToAscii();
          if( writer->Write() != 1 ) throw std::runtime_error( "Could not write fixture" );
        }
      }
      catch( std::exception const & e ) { error = e.what(); }
    }
    MpiWrapper::broadcast( error );
    if( !error.empty() ) throw std::runtime_error( error );
    MpiWrapper::broadcast( m_folder );
  }

  ~FixtureFiles()
  {
    // This directory was created by this fixture and contains only its VTK files.
    if( MpiWrapper::commRank() == 0 && !m_folder.empty() )
    {
      std::error_code error;
      std::filesystem::remove_all( m_folder, error );
      EXPECT_FALSE( error ) << error.message();
    }
  }

  string path() const { return ( std::filesystem::path( m_folder ) / "mesh.vtm" ).string(); }
  string parallelMainPath() const { return m_folder + "/main.pvtu"; }

private:
  string m_folder;
};

template< typename VALIDATE >
void importFixture( FixtureFiles const & files, string const & faceBlocks, int level, VALIDATE const & validate,
                    string const & scatterMethod = "rcb" )
{
  GeosxState state( std::make_unique< CommandLineOptions >( commandLineOptions ) );
  auto const xml = GEOS_FMT( R"xml(<Mesh><VTKMesh name="mesh" file="{}" faceBlocks="{}"
    useGlobalIds="1" scatterMethod="{}" partitionRefinement="0" uniformRefinement="{}" /></Mesh>)xml",
    files.path(), faceBlocks, scatterMethod, level );
  xmlWrapper::xmlDocument document;
  document.loadString( xml );
  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );
  auto meshNode = document.getChild( "Mesh" );
  meshManager.processInputFileRecursive( document, meshNode );
  meshManager.postInputInitializationRecursive();
  DomainPartition domain( "domain", &root );
  meshManager.generateMeshes( domain );
  validate( domain.getMeshBody( "mesh" ).getGroup< CellBlockManagerABC >( keys::cellManager ) );
}

double validateLocalRelations( CellBlockManagerABC const & manager, FaceBlockABC const & block, int arity, bool isJunction, double gap )
{
  auto const positions = manager.getNodePositions();
  auto const ids = manager.getNodeLocalToGlobal();
  std::map< globalIndex, Point > nodes;
  for( localIndex p = 0; p < ids.size(); ++p ) nodes.emplace( ids[p], Point{ positions[p][0], positions[p][1], positions[p][2] } );
  auto const faces = manager.getFaceToNodes();
  auto const elemFaces = block.get2dElemToFaces();
  auto const elemCells = block.get2dElemToElems();
  auto const buckets = block.get2dElemsToCollocatedNodesBuckets();
  auto const elemEdges = block.get2dElemToEdges();
  auto const faceEdges = block.get2dFaceToEdge();
  auto const edgeElements = block.get2dFaceTo2dElems();
  for( localIndex f = 0; f < block.num2dFaces(); ++f )
  {
    EXPECT_GE( faceEdges[f], 0 );
    EXPECT_LT( faceEdges[f], manager.numEdges() );
    EXPECT_GE( edgeElements[f].size(), 1 );
    EXPECT_LE( edgeElements[f].size(), 2 );
  }
  double area = 0;
  for( localIndex e = 0; e < block.num2dElements(); ++e )
  {
    EXPECT_EQ( buckets[e].size(), arity );
    EXPECT_EQ( elemFaces[e].size(), 2 );
    EXPECT_EQ( elemCells.toCellIndex[e].size(), elemFaces[e].size() );
    EXPECT_EQ( elemCells.toBlockIndex[e].size(), elemFaces[e].size() );
    EXPECT_EQ( elemEdges[e].size(), arity );
    if( elemFaces[e].size() == 2 ) { EXPECT_NE( elemFaces[e][0], elemFaces[e][1] ); }
    for( localIndex edge : elemEdges[e] )
    {
      EXPECT_GE( edge, 0 );
      EXPECT_LT( edge, manager.numEdges() );
      bool mapped = false;
      for( localIndex f = 0; f < block.num2dFaces(); ++f )
        for( localIndex element : edgeElements[f] ) mapped = mapped || ( element == e && faceEdges[f] == edge );
      EXPECT_TRUE( mapped );
    }
    for( auto const & bucket : buckets[e] )
    {
      if( bucket.empty() ) { ADD_FAILURE() << "Empty collocation bucket"; continue; }
      auto const first = nodes.find( bucket[0] );
      if( first == nodes.end() ) { ADD_FAILURE() << "Missing contact point"; continue; }
      Point const point = first->second;
      bool const onJunction = isJunction && std::abs( point[0] ) < 1e-12 && std::abs( point[1] ) < 1e-12;
      EXPECT_EQ( bucket.size(), onJunction ? 4 : 2 );
      std::set< globalIndex > unique;
      std::set< double > sides;
      for( globalIndex id : bucket )
      {
        EXPECT_GE( id, idBase );
        EXPECT_TRUE( unique.insert( id ).second );
        auto const found = nodes.find( id );
        if( found == nodes.end() ) { ADD_FAILURE() << "Missing side point"; continue; }
        for( int d = 0; d < 3; ++d )
          if( gap > 0 && d == 2 ) EXPECT_NEAR( std::abs( found->second[d] ), gap, 1e-12 );
          else EXPECT_NEAR( found->second[d], point[d], 1e-12 );
        sides.insert( found->second[2] );
      }
      if( gap > 0 ) { EXPECT_EQ( sides.size(), 2 ); }
    }
    for( localIndex side = 0; side < elemFaces[e].size(); ++side )
    {
      localIndex const f = elemFaces[e][side];
      if( f < 0 || f >= faces.size() ) { ADD_FAILURE() << "Invalid incident face"; continue; }
      EXPECT_EQ( faces[f].size(), arity );
      auto const & volume = manager.getCellBlocks().getGroup< CellBlockABC >( elemCells.toBlockIndex[e][side] );
      localIndex const c = elemCells.toCellIndex[e][side];
      EXPECT_GE( c, 0 );
      EXPECT_LT( c, volume.numElements() );
      if( c < 0 || c >= volume.numElements() ) continue;
      auto const volumeFaces = volume.getElemToFaces();
      bool matched = false;
      for( localIndex vf = 0; vf < volumeFaces.size( 1 ); ++vf ) matched = matched || volumeFaces[c][vf] == f;
      EXPECT_TRUE( matched );
      std::set< globalIndex > faceIds;
      for( localIndex p : faces[f] ) faceIds.insert( ids[p] );
      for( auto const & bucket : buckets[e] )
      {
        std::size_t matches = 0;
        for( globalIndex id : bucket ) matches += faceIds.count( id );
        EXPECT_EQ( matches, 1 );
      }
      if( side == 0 )
      {
        Point normal{};
        for( localIndex p = 0; p < faces[f].size(); ++p )
        {
          Point const a = nodes.at( ids[faces[f][p]] ), b = nodes.at( ids[faces[f][( p + 1 ) % arity]] );
          for( int d = 0; d < 3; ++d ) normal[d] += a[( d + 1 ) % 3] * b[( d + 2 ) % 3] - a[( d + 2 ) % 3] * b[( d + 1 ) % 3];
        }
        area += .5 * std::sqrt( normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2] );
      }
    }
  }
  return area;
}

void validateOriginalPoints( CellBlockManagerABC const & manager, vtkUnstructuredGrid & original )
{
  auto const ids = manager.getNodeLocalToGlobal();
  auto const positions = manager.getNodePositions();
  array1d< globalIndex > localOriginalIds;
  globalIndex const maximum = idBase + 17 * ( original.GetNumberOfPoints() - 1 );
  for( localIndex p = 0; p < ids.size(); ++p )
  {
    if( ids[p] > maximum ) continue;
    EXPECT_GE( ids[p], idBase );
    EXPECT_EQ( ( ids[p] - idBase ) % 17, 0 );
    vtkIdType const index = ( ids[p] - idBase ) / 17;
    if( index < 0 || index >= original.GetNumberOfPoints() ) continue;
    Point expected{};
    original.GetPoints()->GetPoint( index, expected.data() );
    for( int d = 0; d < 3; ++d ) EXPECT_NEAR( positions[p][d], expected[d], 1e-12 );
    localOriginalIds.emplace_back( ids[p] );
  }
  array1d< globalIndex > allOriginalIds;
  MpiWrapper::allGatherv( localOriginalIds.toViewConst(), allOriginalIds );
  std::set< globalIndex > const unique( allOriginalIds.begin(), allOriginalIds.end() );
  EXPECT_EQ( unique.size(), std::size_t( original.GetNumberOfPoints() ) );
  for( vtkIdType p = 0; p < original.GetNumberOfPoints(); ++p ) EXPECT_TRUE( unique.count( idBase + 17 * p ) );
}

void validate( CellBlockManagerABC const & manager, std::vector< string > const & names, int volumeCount, int surfaceCount,
               int arity, int edgeCount, double expectedArea, bool isJunction = false, double gap = 0 )
{
  array1d< globalIndex > localVolumeIds;
  manager.getCellBlocks().forSubGroups< CellBlockABC >( [&]( CellBlockABC const & block )
  {
    localIndex const offset = localVolumeIds.size();
    localVolumeIds.resize( offset + block.numElements() );
    auto const ids = block.localToGlobalMap();
    for( localIndex e = 0; e < block.numElements(); ++e ) localVolumeIds[offset + e] = ids[e];
  } );
  array1d< globalIndex > globalVolumeIds;
  MpiWrapper::allGatherv( localVolumeIds.toViewConst(), globalVolumeIds );
  EXPECT_EQ( globalVolumeIds.size(), volumeCount );
  EXPECT_EQ( ( std::set< globalIndex >( globalVolumeIds.begin(), globalVolumeIds.end() ).size() ), std::size_t( volumeCount ) );
  for( globalIndex id : globalVolumeIds )
  {
    EXPECT_GE( id, 0 );
    EXPECT_LT( id, volumeCount );
  }
  EXPECT_EQ( manager.getFaceBlocks().numSubGroups(), names.size() );
  std::set< globalIndex > allSurfaceIds;
  for( auto const & name : names )
  {
    auto const & block = manager.getFaceBlocks().getGroup< FaceBlockABC >( name );
    EXPECT_EQ( MpiWrapper::sum( block.num2dFaces() ), edgeCount );
    auto const localIds = block.localToGlobalMap();
    array1d< globalIndex > globalIds;
    MpiWrapper::allGatherv( localIds.toViewConst(), globalIds );
    EXPECT_EQ( globalIds.size(), surfaceCount );
    for( globalIndex id : globalIds )
    {
      EXPECT_GE( id, volumeCount );
      EXPECT_LT( id, volumeCount + surfaceCount * names.size() );
      EXPECT_TRUE( allSurfaceIds.insert( id ).second );
    }
    double const area = validateLocalRelations( manager, block, arity, isJunction, gap );
    EXPECT_NEAR( MpiWrapper::sum( area ), expectedArea, 1e-11 );
  }
}

void initializeFixture( FixtureFiles const & files, std::vector< string > const & names, int level,
                        int volumeCount, int surfaceCount, int arity, double area, bool isJunction = false, double gap = 0,
                        bool requireGhosts = false )
{
  GeosxState state( std::make_unique< CommandLineOptions >( commandLineOptions ) );
  string faceBlocks = "{";
  string regions = R"xml(<CellElementRegion name="volume" cellBlocks="{*}" materialList="{}" />)xml";
  for( auto const & name : names )
  {
    if( faceBlocks.size() > 1 ) faceBlocks += ",";
    faceBlocks += name;
    regions += GEOS_FMT( R"xml(<SurfaceElementRegion name="{}" faceBlock="{}" defaultAperture="1e-4" materialList="{{}}" />)xml", name, name );
  }
  faceBlocks += "}";
  auto const xml = GEOS_FMT( R"xml(<Problem><Mesh><VTKMesh name="mesh" file="{}" faceBlocks="{}"
    useGlobalIds="1" scatterMethod="rcb" partitionRefinement="0" uniformRefinement="{}" /></Mesh>
    <ElementRegions>{}</ElementRegions></Problem>)xml", files.path(), faceBlocks, level, regions );
  auto & problem = state.getProblemManager();
  problem.parseInputString( xml );
  problem.problemSetup();
  auto const & mesh = problem.getDomainPartition().getMeshBody( "mesh" ).getBaseDiscretization();
  auto const & nodes = mesh.getNodeManager();
  auto const & faces = mesh.getFaceManager();
  auto const & elements = mesh.getElemManager();
  auto const nodePositions = nodes.referencePosition();
  auto const faceCenters = faces.faceCenter();
  auto const faceNormals = faces.faceNormal();
  auto verifyMaps = []( auto const & objects )
  {
    auto const ids = objects.localToGlobalMap();
    for( localIndex i = 0; i < objects.size(); ++i ) EXPECT_EQ( objects.globalToLocalMap( ids[i] ), i );
  };
  verifyMaps( nodes );
  verifyMaps( faces );
  verifyMaps( mesh.getEdgeManager() );
  localIndex ownedVolumes = 0;
  localIndex ghostVolumes = 0;
  elements.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cells )
  {
    verifyMaps( cells );
    for( localIndex c = 0; c < cells.size(); ++c )
    {
      EXPECT_GT( cells.getElementVolume()[c], 0 );
      ownedVolumes += cells.ghostRank()[c] < 0;
      ghostVolumes += cells.ghostRank()[c] >= 0;
    }
  } );
  EXPECT_EQ( MpiWrapper::sum( ownedVolumes ), volumeCount );
  auto const totalGhostVolumes = MpiWrapper::sum( ghostVolumes );
  if( requireGhosts && MpiWrapper::commSize() > 1 ) { EXPECT_GT( totalGhostVolumes, 0 ); }
  array1d< globalIndex > localSurfaceIds;
  for( auto const & name : names )
  {
    auto const & surfaces = elements.getRegion< SurfaceElementRegion >( name ).getUniqueSubRegion< FaceElementSubRegion >();
    verifyMaps( surfaces );
    localIndex owned = 0;
    double ownedArea = 0;
    auto const buckets = surfaces.get2dElemToCollocatedNodesBuckets();
    auto const & surfaceFaces = surfaces.faceList();
    auto const & relation = surfaces.getToCellRelation();
    auto const surfaceCenters = surfaces.getElementCenter();
    auto const & surfaceNodes = surfaces.nodeList();
    auto const & faceNodes = faces.nodeList();
    for( localIndex e = 0; e < surfaces.size(); ++e )
    {
      if( surfaces.ghostRank()[e] < 0 )
      {
        ++owned;
        ownedArea += surfaces.getElementArea()[e];
        localSurfaceIds.emplace_back( surfaces.localToGlobalMap()[e] );
      }
      EXPECT_EQ( buckets[e].size(), arity );
      for( auto const & bucket : buckets[e] )
      {
        if( bucket.empty() ) { ADD_FAILURE() << "Empty final collocation bucket"; continue; }
        auto const found = nodes.globalToLocalMap().find( bucket[0] );
        if( found == nodes.globalToLocalMap().end() ) { ADD_FAILURE() << "Missing final contact point"; continue; }
        auto const point = nodePositions[found->second];
        bool const atJunction = isJunction && std::abs( point[0] ) < 1e-12 && std::abs( point[1] ) < 1e-12;
        EXPECT_EQ( bucket.size(), atJunction ? 4 : 2 );
        for( globalIndex id : bucket )
        {
          auto const side = nodes.globalToLocalMap().find( id );
          if( side == nodes.globalToLocalMap().end() ) { ADD_FAILURE() << "Missing final side point"; continue; }
          for( int d = 0; d < 3; ++d )
            if( gap > 0 && d == 2 ) EXPECT_NEAR( std::abs( nodePositions[side->second][d] ), gap, 1e-12 );
            else EXPECT_NEAR( nodePositions[side->second][d], point[d], 1e-12 );
        }
      }
      for( int side = 0; side < 2; ++side )
      {
        localIndex const f = surfaceFaces[e][side];
        EXPECT_GE( f, 0 );
        EXPECT_LT( f, faces.size() );
        if( f < 0 || f >= faces.size() ) continue;
        EXPECT_EQ( faces.nodeList()[f].size(), arity );
        EXPECT_EQ( surfaceNodes[e].size(), 2 * arity );
        if( surfaceNodes[e].size() == 2 * arity )
        {
          if( level > 0 )
          {
            for( int p = 0; p < arity; ++p ) EXPECT_EQ( surfaceNodes[e][side * arity + p], faceNodes[f][p] );
          }
          else
          {
            // Develop reorders Kf1 without rebuilding the surface node halves.
            // Zero refinement preserves that ordering; membership still agrees.
            std::set< localIndex > surfaceSide, faceSide;
            for( int p = 0; p < arity; ++p )
            {
              surfaceSide.insert( surfaceNodes[e][side * arity + p] );
              faceSide.insert( faceNodes[f][p] );
            }
            EXPECT_EQ( surfaceSide.size(), std::size_t( arity ) );
            EXPECT_EQ( surfaceSide, faceSide );
          }
        }
        EXPECT_NEAR( faces.faceArea()[f], surfaces.getElementArea()[e], 1e-12 );
        localIndex const region = relation.m_toElementRegion[e][side];
        localIndex const subRegion = relation.m_toElementSubRegion[e][side];
        localIndex const cell = relation.m_toElementIndex[e][side];
        EXPECT_GE( region, 0 );
        EXPECT_GE( subRegion, 0 );
        EXPECT_GE( cell, 0 );
        if( region < 0 || subRegion < 0 || cell < 0 ) continue;
        auto const & incident = elements.getRegion( region ).getSubRegion< CellElementSubRegion >( subRegion );
        EXPECT_LT( cell, incident.size() );
        bool matched = false;
        if( cell < incident.size() )
        {
          for( int cf = 0; cf < incident.faceList().size( 1 ); ++cf ) matched = matched || incident.faceList()[cell][cf] == f;
          auto const cellCenters = incident.getElementCenter();
          double outward = 0;
          for( int d = 0; d < 3; ++d ) outward += faceNormals[f][d] * ( faceCenters[f][d] - cellCenters[cell][d] );
          EXPECT_GT( outward, 0 );
        }
        EXPECT_TRUE( matched );
      }
      if( surfaceFaces[e][0] >= 0 && surfaceFaces[e][0] < faces.size() &&
          surfaceFaces[e][1] >= 0 && surfaceFaces[e][1] < faces.size() )
      {
        double normalProduct = 0;
        for( int d = 0; d < 3; ++d )
        {
          // SurfaceElementSubRegion defines its center as the nodal barycenter,
          // which can differ from an area centroid on polygon-cap child quads.
          double center = 0;
          for( int side = 0; side < 2; ++side )
            for( localIndex p : faceNodes[surfaceFaces[e][side]] ) center += nodePositions[p][d];
          EXPECT_NEAR( surfaceCenters[e][d], center / ( 2 * arity ), 1e-12 );
          normalProduct += faceNormals[surfaceFaces[e][0]][d] * faceNormals[surfaceFaces[e][1]][d];
        }
        EXPECT_NEAR( normalProduct, -1, 1e-12 );
      }
    }
    EXPECT_EQ( MpiWrapper::sum( owned ), surfaceCount );
    EXPECT_NEAR( MpiWrapper::sum( ownedArea ), area, 1e-11 );
  }
  array1d< globalIndex > globalSurfaceIds;
  MpiWrapper::allGatherv( localSurfaceIds.toViewConst(), globalSurfaceIds );
  EXPECT_EQ( globalSurfaceIds.size(), surfaceCount * names.size() );
  EXPECT_EQ( ( std::set< globalIndex >( globalSurfaceIds.begin(), globalSurfaceIds.end() ).size() ), std::size_t( globalSurfaceIds.size() ) );
  for( globalIndex id : globalSurfaceIds )
  {
    if( level > 0 )
    {
      EXPECT_GE( id, volumeCount );
      EXPECT_LT( id, volumeCount + surfaceCount * names.size() );
    }
    else
    {
      // The coarse ID namespace is preserved at level zero.
      EXPECT_GT( id, idBase + 20000 + 13 * ( volumeCount - 1 ) );
    }
  }
}
} // namespace

TEST( VTKRefinedFractureImport, TriangleAndSeparatedSidesAtTwoLevels )
{
  for( double gap : { 0., .125 } )
  {
    auto const original = triangle( gap );
    FixtureFiles files( [&original] { return original; } );
    for( int level : { 1, 2 } )
    {
      SCOPED_TRACE( "triangle level " + std::to_string( level ) + " gap " + std::to_string( gap ) );
      int const growth = level == 1 ? 8 : 64;
      importFixture( files, "{triangle}", level, [&]( CellBlockManagerABC const & manager )
      {
        validateOriginalPoints( manager, *original.main );
        validate( manager, { "triangle" }, ( MpiWrapper::commSize() + 1 ) * growth, level == 1 ? 4 : 16, 3,
                  level == 1 ? 9 : 30, .5, false, gap );
      } );
    }
  }
}

TEST( VTKRefinedFractureImport, AllSevenPolygonCapAritiesAtTwoLevels )
{
  for( int n = 5; n <= 11; ++n )
  {
    auto const original = cap( n );
    FixtureFiles files( [&original] { return original; } );
    for( int level : { 1, 2 } )
    {
      SCOPED_TRACE( "cap " + std::to_string( n ) + " level " + std::to_string( level ) );
      int const volumeCount = 4 * n * ( level == 1 ? 1 : 8 ) + ( MpiWrapper::commSize() - 1 ) * ( level == 1 ? 8 : 64 );
      importFixture( files, "{cap}", level, [&]( CellBlockManagerABC const & manager )
      {
        validateOriginalPoints( manager, *original.main );
        validate( manager, { "cap" }, volumeCount, n * ( level == 1 ? 1 : 4 ), 4, n * ( level == 1 ? 3 : 10 ),
                  n * .5 * std::sin( 2 * std::acos( -1. ) / n ) );
      } );
    }
  }
}

TEST( VTKRefinedFractureImport, NamedBlocksAndFourSideJunctionAtTwoLevels )
{
  auto const original = junction();
  FixtureFiles files( [&original] { return original; } );
  for( int level : { 1, 2 } )
  {
    SCOPED_TRACE( "junction level " + std::to_string( level ) );
    importFixture( files, "{faultX, faultY}", level, [&]( CellBlockManagerABC const & manager )
    {
      validateOriginalPoints( manager, *original.main );
      validate( manager, { "faultX", "faultY" }, ( MpiWrapper::commSize() + 3 ) * ( level == 1 ? 8 : 64 ),
                level == 1 ? 8 : 32, 4, level == 1 ? 22 : 76, 2, true );
    } );
  }
}

TEST( VTKRefinedFractureImport, FullInitializationOfNamedJunctionRegions )
{
  FixtureFiles files( [] { return junction(); } );
  for( int level : { 1, 2 } )
  {
    SCOPED_TRACE( "full junction level " + std::to_string( level ) );
    initializeFixture( files, { "faultX", "faultY" }, level,
      ( MpiWrapper::commSize() + 3 ) * ( level == 1 ? 8 : 64 ), level == 1 ? 8 : 32, 4, 2, true );
  }
}

TEST( VTKRefinedFractureImport, FullInitializationOfTrianglesAndSeparatedSides )
{
  for( double gap : { 0., .125 } )
  {
    FixtureFiles files( [gap] { return triangle( gap ); } );
    for( int level : { 0, 1, 2 } )
    {
      SCOPED_TRACE( "full triangle level " + std::to_string( level ) + " gap " + std::to_string( gap ) );
      int const growth = level == 0 ? 1 : level == 1 ? 8 : 64;
      int const surfaceCount = level == 0 ? 1 : level == 1 ? 4 : 16;
      initializeFixture( files, { "triangle" }, level, ( MpiWrapper::commSize() + 1 ) * growth, surfaceCount, 3, .5, false, gap );
    }
  }
}

TEST( VTKRefinedFractureImport, FullInitializationOfAllSevenPolygonCaps )
{
  for( int n = 5; n <= 11; ++n )
  {
    FixtureFiles files( [n] { return cap( n ); } );
    for( int level : { 1, 2 } )
    {
      SCOPED_TRACE( "full cap " + std::to_string( n ) + " level " + std::to_string( level ) );
      int const volumeCount = 4 * n * ( level == 1 ? 1 : 8 ) + ( MpiWrapper::commSize() - 1 ) * ( level == 1 ? 8 : 64 );
      initializeFixture( files, { "cap" }, level, volumeCount, n * ( level == 1 ? 1 : 4 ), 4,
        n * .5 * std::sin( 2 * std::acos( -1. ) / n ) );
    }
  }
}

TEST( VTKRefinedFractureImport, FullInitializationOfConnectedFractureAcrossRankBoundaries )
{
  FixtureFiles files( [] { return strip(); } );
  for( int level : { 0, 1, 2 } )
  {
    SCOPED_TRACE( "connected strip level " + std::to_string( level ) );
    int const growth = level == 0 ? 1 : level == 1 ? 8 : 64;
    int const surfaces = level == 0 ? 2 : level == 1 ? 8 : 32;
    initializeFixture( files, { "strip" }, level, ( std::max( MpiWrapper::commSize(), 2 ) + 2 ) * growth,
                       surfaces, 4, 2, false, 0, true );
  }
}

TEST( VTKRefinedFractureImport, DeclaredSurfaceScalarAndVectorImportsSurviveRefinement )
{
  FixtureFiles files( []
  {
    auto fixture = triangle( 0 );
    auto const & surface = fixture.blocks.at( "triangle" );
    vtkNew< vtkFloatArray > aperture;
    aperture->SetName( "importAperture" );
    aperture->SetNumberOfTuples( surface->GetNumberOfCells() );
    aperture->FillValue( 2e-4f );
    surface->GetCellData()->SetScalars( aperture );
    vtkNew< vtkDoubleArray > tangent, unused;
    tangent->SetName( "importTangent" );
    tangent->SetNumberOfComponents( 3 );
    for( vtkIdType c = 0; c < surface->GetNumberOfCells(); ++c ) tangent->InsertNextTuple3( .6, .8, 0 );
    surface->GetCellData()->SetVectors( tangent );
    unused->SetName( "unusedWideField" );
    unused->SetNumberOfComponents( 128 );
    unused->SetNumberOfTuples( surface->GetNumberOfCells() );
    unused->FillValue( 42 );
    surface->GetCellData()->AddArray( unused );
    return fixture;
  } );
  for( int level : { 0, 1, 2 } )
  {
    SCOPED_TRACE( "surface field import level " + std::to_string( level ) );
    GeosxState state( std::make_unique< CommandLineOptions >( commandLineOptions ) );
    auto const xml = GEOS_FMT( R"xml(<Problem><Mesh><VTKMesh name="mesh" file="{}" faceBlocks="{{triangle}}"
      useGlobalIds="1" scatterMethod="rcb" partitionRefinement="0" uniformRefinement="{}"
      surfacicFieldsToImport="{{importAperture,importTangent}}"
      surfacicFieldsInGEOS="{{elementAperture,tangentVector1}}" /></Mesh>
      <ElementRegions><CellElementRegion name="volume" cellBlocks="{{*}}" materialList="{{}}" />
      <SurfaceElementRegion name="fault" faceBlock="triangle" defaultAperture="1e-4" materialList="{{}}" />
      </ElementRegions></Problem>)xml", files.path(), level );
    auto & problem = state.getProblemManager();
    problem.parseInputString( xml );
    problem.problemSetup();
    auto const & mesh = problem.getDomainPartition().getMeshBody( "mesh" ).getBaseDiscretization();
    auto const & surface = mesh.getElemManager().getRegion< SurfaceElementRegion >( "fault" )
      .getUniqueSubRegion< FaceElementSubRegion >();
    auto const aperture = surface.getElementAperture();
    auto const tangent = surface.getTangentVector1();
    localIndex owned = 0;
    for( localIndex c = 0; c < surface.size(); ++c )
    {
      if( surface.ghostRank()[c] >= 0 ) continue;
      ++owned;
      EXPECT_DOUBLE_EQ( aperture[c], static_cast< double >( 2e-4f ) );
      EXPECT_DOUBLE_EQ( tangent[c][0], .6 );
      EXPECT_DOUBLE_EQ( tangent[c][1], .8 );
      EXPECT_DOUBLE_EQ( tangent[c][2], 0 );
    }
    EXPECT_EQ( MpiWrapper::sum( owned ), level == 0 ? 1 : level == 1 ? 4 : 16 );
  }
}

TEST( VTKRefinedFractureImport, TwoPieceParallelMainWithFracturesKeepsEveryCell )
{
  int const ranks = MpiWrapper::commSize();
  if( ranks < 4 ) GTEST_SKIP() << "Two-piece input with fewer pieces than ranks requires at least four ranks";
  FixtureFiles files( [] { return triangle( 0 ); }, true );
  for( auto method : { geos::vtk::ScatterMethod::kdtree, geos::vtk::ScatterMethod::rcb } )
  {
    auto pieces = geos::vtk::loadAllMeshes( Path( files.parallelMainPath().c_str() ), "main", {} );
    vtkIdType const count = pieces.getMainMesh()->GetNumberOfCells();
    EXPECT_EQ( MpiWrapper::sum( count ), ranks + 1 );
    EXPECT_EQ( MpiWrapper::min( count ), 0 );
    vtkIdType const rootCount = MpiWrapper::sum( MpiWrapper::commRank() == 0 ? count : vtkIdType{ 0 } );
    EXPECT_LT( rootCount, ranks + 1 );
    auto bundle = geos::vtk::loadAllMeshes( Path( files.path().c_str() ), "main", { "triangle" } );
    array1d< integer > partitions( 3 );
    partitions[0] = ranks;
    partitions[1] = partitions[2] = 1;
    auto redistributed = geos::vtk::redistributeMeshes( 0, pieces.getMainMesh(), bundle.getFaceBlocks(), MPI_COMM_GEOS,
      method, partitions.toViewConst(), geos::vtk::PartitionMethod::parmetis, 0, 0, 1, "" );
    EXPECT_EQ( MpiWrapper::sum( redistributed.getMainMesh()->GetNumberOfCells() ), ranks + 1 );
    EXPECT_GT( redistributed.getMainMesh()->GetNumberOfCells(), 0 );
    geos::vtk::refineUniformly( redistributed, 2, {}, MPI_COMM_GEOS );
    EXPECT_EQ( MpiWrapper::sum( redistributed.getMainMesh()->GetNumberOfCells() ), ( ranks + 1 ) * 64 );
    EXPECT_EQ( MpiWrapper::sum( redistributed.getFaceBlocks().at( "triangle" )->GetNumberOfCells() ), 16 );
    for( int level : { 0, 1, 2 } )
    {
      int const growth = level == 0 ? 1 : level == 1 ? 8 : 64;
      importFixture( files, "{triangle}", level, [&]( CellBlockManagerABC const & manager )
      {
        localIndex cells = 0;
        manager.getCellBlocks().forSubGroups< CellBlockABC >( [&]( CellBlockABC const & block ) { cells += block.numElements(); } );
        EXPECT_EQ( MpiWrapper::sum( cells ), ( ranks + 1 ) * growth );
      }, method == geos::vtk::ScatterMethod::kdtree ? "kdtree" : "rcb" );
    }
  }
}

TEST( VTKRefinedFractureImport, GraphColoringIncludesIsolatedRanksAndCustomNeighbors )
{
  int const rank = MpiWrapper::commRank();
  int const ranks = MpiWrapper::commSize();
  SpatialPartition partition;
  partition.setPartitions( 1, 1, ranks );
  EXPECT_FALSE( partition.hasMetisNeighborList() );
  stdVector< int > neighbors;
  if( ranks > 1 && rank < 2 ) neighbors.push_back( 1 - rank );
  partition.setMetisNeighborList( neighbors );
  EXPECT_TRUE( partition.hasMetisNeighborList() );
  auto verifyColors = [&]( int neighbor )
  {
    int const color = neighbor < 0 ? partition.getColor() : partition.getColor(
      rank == 0 ? std::set< int >{ neighbor } : rank == neighbor ? std::set< int >{ 0 } : std::set< int >{} );
    array1d< int > colors( ranks );
    MpiWrapper::allgather( &color, 1, colors.data(), 1, MPI_COMM_GEOS );
    EXPECT_GE( color, 0 );
    if( ranks > 1 ) { EXPECT_NE( colors[0], colors[neighbor < 0 ? 1 : neighbor] ); }
  };
  verifyColors( -1 );
  if( ranks > 1 ) verifyColors( ranks - 1 );
  partition.setPartitions( 1, 1, ranks );
  EXPECT_FALSE( partition.hasMetisNeighborList() );
  EXPECT_TRUE( partition.getMetisNeighborList().empty() );
}

int main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
