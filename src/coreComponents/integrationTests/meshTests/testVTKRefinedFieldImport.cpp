/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/** @file testVTKRefinedFieldImport.cpp */
#include "common/MpiWrapper.hpp"
#include "LvArray/src/system.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/DomainPartition.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"

#include <gtest/gtest.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkTypeInt64Array.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <vector>

namespace geos
{
/** Test-only real solver registration: exercise the normal delayed import of
 * arbitrary regular scalar/vector/tensor fields alongside actual materials.
 */
class RefinementImportProbe : public SolidMechanicsLagrangianFEM
{
public:
  RefinementImportProbe( string const & name, dataRepository::Group * parent ): SolidMechanicsLagrangianFEM( name, parent ) {}
  static string catalogName() { return "RefinementImportProbe"; }
  string getCatalogName() const override { return catalogName(); }
  void registerDataOnMesh( dataRepository::Group & bodies ) override
  {
    SolidMechanicsLagrangianFEM::registerDataOnMesh( bodies );
    bodies.forSubGroups< MeshBody >( []( MeshBody & body )
    {
      body.forMeshLevels( []( MeshLevel & mesh )
      {
        mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( []( CellElementSubRegion & cells )
        {
          if( !cells.hasWrapper( "probeScalar" ) )
          {
            auto & scalar = cells.registerWrapper< array1d< real64 > >( "probeScalar" ).setDefaultValue( -1 ).reference();
            scalar.resize( cells.size() );
            scalar.setValues< serialPolicy >( -1 );
            auto & vector = cells.registerWrapper< array2d< real64 > >( "probeVector" ).reference();
            vector.resize( cells.size(), 3 );
            vector.setValues< serialPolicy >( -1 );
            auto & tensor = cells.registerWrapper< array2d< real64 > >( "probeTensor" ).reference();
            tensor.resize( cells.size(), 9 );
            tensor.setValues< serialPolicy >( -1 );
          }
        } );
      } );
    } );
  }
};
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, RefinementImportProbe, string const &, dataRepository::Group * const )
} // namespace geos

using namespace geos;
namespace
{
CommandLineOptions commandLineOptions;
globalIndex const idBase = sizeof( globalIndex ) == 8 && sizeof( vtkIdType ) == 8 ? INT64_C( 9007199254741001 ) : 10001;
using Point = std::array< double, 3 >;
std::array< string, 4 > const regions{ "hexRegion", "pyramidRegion", "tetRegion", "prismRegion" };
std::array< string, 4 > const materials{ "hexRock", "pyramidRock", "tetRock", "prismRock" };

real64 seed( globalIndex cell, int component ) { return 100 * ( cell + 1 ) + real64( component ) / 8; }
real64 stressSeed( globalIndex cell, int component ) { return 1e6 + 256 * cell + real64( component ) / 8; }
real64 bulkSeed( globalIndex cell ) { return 1.1e9 + 1e6 * cell; }

class FixtureFile
{
public:
  FixtureFile( bool polyhedra, bool singlePrecision )
  {
    string error;
    if( MpiWrapper::commRank() == 0 )
    {
      try
      {
        LvArray::system::FloatingPointExceptionGuard guard;
        vtkNew< vtkUnstructuredGrid > grid;
        vtkNew< vtkPoints > points;
        points->SetDataTypeToDouble();
        grid->SetPoints( points );
        std::map< Point, vtkIdType > sharedPoints;
        auto cell = [&]( int type, std::vector< Point > const & coordinates, std::vector< std::vector< int > > const & faces )
        {
          std::vector< vtkIdType > corners;
          for( auto const & point : coordinates )
          {
            auto found = sharedPoints.find( point );
            if( found == sharedPoints.end() ) found = sharedPoints.emplace( point, points->InsertNextPoint( point.data() ) ).first;
            corners.push_back( found->second );
          }
          if( polyhedra )
          {
            vtkNew< vtkCellArray > faceCells;
            for( auto const & face : faces )
            {
              std::vector< vtkIdType > mapped;
              for( int p : face ) mapped.push_back( corners.at( p ) );
              faceCells->InsertNextCell( mapped.size(), mapped.data() );
            }
            grid->InsertNextCell( VTK_POLYHEDRON, corners.size(), corners.data(), faceCells );
          }
          else grid->InsertNextCell( type, corners.size(), corners.data() );
        };
        int const copies = ( MpiWrapper::commSize() + 3 ) / 4;
        for( int copy = 0; copy < copies; ++copy )
        {
          double const shift = 10 * copy;
          auto translated = [shift]( std::vector< Point > coordinates )
          {
            for( auto & point : coordinates ) point[0] += shift;
            return coordinates;
          };
          cell( VTK_HEXAHEDRON, translated( { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
                                              { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } } ),
                { { 0, 3, 2, 1 }, { 4, 5, 6, 7 }, { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 }, { 3, 0, 4, 7 } } );
          cell( VTK_PYRAMID, translated( { { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 }, { 1, 0, 1 }, { 2, .5, .5 } } ),
                { { 0, 3, 2, 1 }, { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 } } );
          cell( VTK_TETRA, translated( { { 1, 0, 0 }, { 1, 1, 0 }, { 2, .5, .5 }, { 1.5, .5, -.5 } } ),
                { { 0, 2, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 0, 3 } } );
          cell( VTK_PENTAGONAL_PRISM,
                translated( { { 0, 0, 0 }, { 0, 1, 0 }, { -1, 1, 0 }, { -1.5, .5, 0 }, { -1, 0, 0 },
                              { 0, 0, 1 }, { 0, 1, 1 }, { -1, 1, 1 }, { -1.5, .5, 1 }, { -1, 0, 1 } } ),
                { { 4, 3, 2, 1, 0 }, { 5, 6, 7, 8, 9 }, { 0, 1, 6, 5 }, { 1, 2, 7, 6 },
                  { 2, 3, 8, 7 }, { 3, 4, 9, 8 }, { 4, 0, 5, 9 } } );
        }
        vtkNew< vtkTypeInt64Array > nodeIds;
        nodeIds->SetName( "pointIds" );
        for( vtkIdType p = 0; p < points->GetNumberOfPoints(); ++p ) nodeIds->InsertNextValue( idBase + 17 * p );
        grid->GetPointData()->SetGlobalIds( nodeIds );
        vtkNew< vtkTypeInt64Array > cellIds;
        cellIds->SetName( "cellIds" );
        for( vtkIdType c = 0; c < grid->GetNumberOfCells(); ++c ) cellIds->InsertNextValue( idBase + 20000 + 13 * c );
        grid->GetCellData()->SetGlobalIds( cellIds );
        auto addField = [&]( string const & name, int components, auto value )
        {
          vtkSmartPointer< vtkDataArray > field = singlePrecision ? vtkSmartPointer< vtkDataArray >( vtkSmartPointer< vtkFloatArray >::New() ) :
                                                                  vtkSmartPointer< vtkDataArray >( vtkSmartPointer< vtkDoubleArray >::New() );
          field->SetName( name.c_str() );
          field->SetNumberOfComponents( components );
          field->SetNumberOfTuples( grid->GetNumberOfCells() );
          for( int k = 0; k < components; ++k ) field->SetComponentName( k, ( "component" + std::to_string( k ) ).c_str() );
          for( vtkIdType c = 0; c < grid->GetNumberOfCells(); ++c )
            for( int k = 0; k < components; ++k ) field->SetComponent( c, k, value( c, k ) );
          grid->GetCellData()->AddArray( field );
        };
        addField( "attribute", 1, []( vtkIdType, int ) { return 7; } );
        addField( "scalar", 1, []( vtkIdType c, int k ) { return seed( c, k ); } );
        addField( "vector", 3, []( vtkIdType c, int k ) { return seed( c, k ); } );
        addField( "tensor", 9, []( vtkIdType c, int k ) { return seed( c, k ); } );
        grid->GetCellData()->SetActiveScalars( "scalar" );
        grid->GetCellData()->SetActiveVectors( "vector" );
        grid->GetCellData()->SetActiveTensors( "tensor" );
        for( auto const & material : materials )
        {
          addField( material + "Stress", 6, []( vtkIdType c, int k ) { return stressSeed( c, k ); } );
          // Use double for moduli so the seed is exactly represented.
          vtkNew< vtkDoubleArray > bulk;
          bulk->SetName( ( material + "Bulk" ).c_str() );
          for( vtkIdType c = 0; c < grid->GetNumberOfCells(); ++c ) bulk->InsertNextValue( bulkSeed( c ) );
          grid->GetCellData()->AddArray( bulk );
        }
        auto const stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        m_directory = std::filesystem::temp_directory_path() / ( "geos-refined-field-" + std::to_string( stamp ) );
        if( !std::filesystem::create_directory( m_directory ) ) throw std::runtime_error( "Fixture directory already exists" );
        vtkNew< vtkXMLUnstructuredGridWriter > writer;
        writer->SetFileName( path().c_str() );
        writer->SetInputData( grid );
        writer->SetDataModeToAscii();
        if( writer->Write() != 1 ) throw std::runtime_error( "Could not write fixture" );
      }
      catch( std::exception const & failure ) { error = failure.what(); }
    }
    MpiWrapper::broadcast( error );
    if( !error.empty() ) throw std::runtime_error( error );
    string directory = m_directory.string();
    MpiWrapper::broadcast( directory );
    m_directory = directory;
  }
  ~FixtureFile()
  {
    if( MpiWrapper::commRank() == 0 && !m_directory.empty() )
    {
      std::error_code error;
      std::filesystem::remove_all( m_directory, error );
      EXPECT_FALSE( error ) << error.message();
    }
  }
  string path() const { return ( m_directory / "mesh.vtu" ).string(); }
private:
  std::filesystem::path m_directory;
};

void initialize( FixtureFile const & fixture, int level )
{
  GeosxState state( std::make_unique< CommandLineOptions >( commandLineOptions ) );
  auto & problem = state.getProblemManager();
  string regionXml, materialXml;
  string fields = "{scalar,vector,tensor", targets = "{probeScalar,probeVector,probeTensor";
  std::array< string, 4 > const blocks{ "7_hexahedra", "7_pyramids", "7_tetrahedra", "7_pentagonalPrisms" };
  for( int kind = 0; kind < 4; ++kind )
  {
    regionXml += GEOS_FMT( R"xml(<CellElementRegion name="{}" cellBlocks="{{{}}}" materialList="{{{}}}"/>)xml",
                          regions[kind], blocks[kind], materials[kind] );
    materialXml += GEOS_FMT( R"xml(<ElasticIsotropic name="{}" defaultDensity="2700"
      defaultBulkModulus="5e9" defaultShearModulus="{}"/>)xml", materials[kind], ( kind + 1 ) * 1e9 );
    fields += "," + materials[kind] + "Stress," + materials[kind] + "Bulk";
    targets += "," + materials[kind] + "_stress," + materials[kind] + "_bulkModulus";
  }
  fields += "}";
  targets += "}";
  problem.parseInputString( GEOS_FMT( R"xml(<Problem>
    <Mesh><VTKMesh name="mesh" file="{}" useGlobalIds="1" scatterMethod="rcb" partitionRefinement="0"
      uniformRefinement="{}" fieldsToImport="{}" fieldNamesInGEOS="{}"/></Mesh>
    <Solvers gravityVector="{{0,0,0}}"><RefinementImportProbe name="probe" discretization="FE1"
      targetRegions="{{hexRegion,pyramidRegion,tetRegion,prismRegion}}"/></Solvers>
    <NumericalMethods><FiniteElements><FiniteElementSpace name="FE1" order="1"/></FiniteElements></NumericalMethods>
    <ElementRegions>{}</ElementRegions><Constitutive>{}</Constitutive>
    </Problem>)xml", fixture.path(), level, fields, targets, regionXml, materialXml ) );
  problem.problemSetup();
  auto & mesh = problem.getDomainPartition().getMeshBody( "mesh" ).getBaseDiscretization();
  globalIndex localGhosts = 0;
  int const copies = ( MpiWrapper::commSize() + 3 ) / 4;
  for( int kind = 0; kind < 4; ++kind )
  {
    globalIndex owned = 0;
    double volume = 0;
    auto const & region = mesh.getElemManager().getRegion( regions[kind] );
    EXPECT_EQ( region.getMaterialList(), ( string_array{ materials[kind] } ) );
    region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cells )
    {
      auto const ids = cells.localToGlobalMap();
      auto const ghosts = cells.ghostRank();
      auto const roots = level > 0 ? cells.getReference< array1d< globalIndex > >( "_geosUniformRootCellId" ).toViewConst() : ids;
      auto const scalar = cells.getReference< array1d< real64 > >( "probeScalar" ).toViewConst();
      auto const vector = cells.getReference< array2d< real64 > >( "probeVector" ).toViewConst();
      auto const tensor = cells.getReference< array2d< real64 > >( "probeTensor" ).toViewConst();
      auto const & material = cells.getConstitutiveModel< constitutive::ElasticIsotropic >( materials[kind] );
      auto const stress = material.getStress();
      auto const bulk = material.getBulkModulus();
      auto const shear = material.getShearModulus();
      for( localIndex c = 0; c < cells.size(); ++c )
      {
        EXPECT_EQ( cells.globalToLocalMap( ids[c] ), c );
        globalIndex const source = ( roots[c] - idBase - 20000 ) / 13;
        EXPECT_EQ( roots[c], idBase + 20000 + 13 * source );
        EXPECT_GE( source, 0 );
        EXPECT_LT( source, 4 * copies );
        EXPECT_EQ( source % 4, kind );
        EXPECT_DOUBLE_EQ( scalar[c], seed( source, 0 ) );
        for( int k = 0; k < 3; ++k ) EXPECT_DOUBLE_EQ( vector[c][k], seed( source, k ) );
        for( int k = 0; k < 9; ++k ) EXPECT_DOUBLE_EQ( tensor[c][k], seed( source, k ) );
        EXPECT_DOUBLE_EQ( bulk[c], bulkSeed( source ) );
        EXPECT_DOUBLE_EQ( shear[c], ( kind + 1 ) * 1e9 );
        EXPECT_GT( stress.size( 1 ), 0 );
        for( localIndex q = 0; q < stress.size( 1 ); ++q )
          for( int k = 0; k < 6; ++k ) EXPECT_DOUBLE_EQ( stress[c][q][k], stressSeed( source, k ) );
        if( ghosts[c] < 0 )
        {
          ++owned;
          volume += cells.getElementVolume()[c];
        }
        else ++localGhosts;
      }
    } );
    int const counts[3][4] = { { 1, 1, 1, 1 }, { 8, 10, 8, 10 }, { 64, 92, 64, 80 } };
    double const volumes[4] = { 1, 1. / 3, 1. / 8, 1.25 };
    EXPECT_EQ( MpiWrapper::sum( owned ), copies * counts[level][kind] );
    EXPECT_NEAR( MpiWrapper::sum( volume ), copies * volumes[kind], 1e-11 );
  }
  if( MpiWrapper::commSize() > 1 ) { EXPECT_GT( MpiWrapper::sum( localGhosts ), 0 ); }
}
} // namespace

TEST( VTKRefinedFieldImport, ScalarVectorTensorAndMaterialFieldsKeepSourceAlignment )
{
  for( bool polyhedra : { false, true } )
    for( bool singlePrecision : { false, true } )
    {
      FixtureFile fixture( polyhedra, singlePrecision );
      for( int level : { 0, 1, 2 } )
      {
        SCOPED_TRACE( "level " + std::to_string( level ) + " polyhedra " + std::to_string( polyhedra ) +
                      " float " + std::to_string( singlePrecision ) );
        initialize( fixture, level );
      }
    }
}

int main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
