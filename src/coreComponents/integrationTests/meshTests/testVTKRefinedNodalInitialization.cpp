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

/** @file testVTKRefinedNodalInitialization.cpp */
#include "common/MpiWrapper.hpp"
#include "LvArray/src/system.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

#include <gtest/gtest.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStringArray.h>
#include <vtkTypeInt64Array.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

using namespace geos;

namespace
{
CommandLineOptions commandLineOptions;
globalIndex const idBase = sizeof( globalIndex ) == 8 && sizeof( vtkIdType ) == 8 ? INT64_C( 9007199254741001 ) : 10001;
using TopologyKey = std::vector< globalIndex >;

class FixtureFile
{
public:
  FixtureFile( bool vertex, int copies )
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
        // Coordinates describe this fixture. The namespace deliberately keeps
        // unrelated coincident components distinct in the input topology.
        std::map< std::array< int, 4 >, vtkIdType > vertices;
        auto hex = [&]( int x, int y, int z, int space )
        {
          vtkIdType corners[8];
          constexpr int offsets[8][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 0, 1, 0 },
            { 0, 0, 1 }, { 1, 0, 1 }, { 1, 1, 1 }, { 0, 1, 1 } };
          for( int p = 0; p < 8; ++p )
          {
            std::array< int, 4 > const key{ x + offsets[p][0], y + offsets[p][1], z + offsets[p][2], space };
            auto found = vertices.find( key );
            if( found == vertices.end() )
            {
              vtkIdType const index = points->InsertNextPoint( key[0], key[1], key[2] );
              found = vertices.emplace( key, index ).first;
            }
            corners[p] = found->second;
          }
          grid->InsertNextCell( VTK_HEXAHEDRON, 8, corners );
        };
        for( int copy = 0; copy < copies; ++copy )
          for( int z = 0; z < ( vertex ? 2 : 1 ); ++z )
            for( int y = 0; y < 2; ++y )
              for( int x = 0; x < 2; ++x ) hex( x - 1, y - 1, vertex ? z - 1 : z, copy );
        int const groups = ( vertex ? 8 : 4 ) * copies;
        for( int cell = groups; cell < MpiWrapper::commSize(); ++cell ) hex( 10 + 2 * cell, 5, 0, copies );
        vtkNew< vtkTypeInt64Array > nodeIds;
        nodeIds->SetName( "pointIds" );
        for( vtkIdType p = 0; p < points->GetNumberOfPoints(); ++p ) nodeIds->InsertNextValue( idBase + 17 * p );
        grid->GetPointData()->SetGlobalIds( nodeIds );
        // Declared node sets must include new boundary nodes without extending
        // isolated corners. Unselected categorical/string metadata is irrelevant
        // to interpolation and must not prevent normal positive XML refinement.
        vtkNew< vtkUnsignedCharArray > left, bottom, corner;
        left->SetName( "leftMask" );
        bottom->SetName( "bottomMask" );
        corner->SetName( "cornerMask" );
        vtkNew< vtkIntArray > unusedInteger;
        unusedInteger->SetName( "unusedCategory" );
        vtkNew< vtkStringArray > unusedString;
        unusedString->SetName( "unusedLabel" );
        for( vtkIdType p = 0; p < points->GetNumberOfPoints(); ++p )
        {
          double position[3];
          points->GetPoint( p, position );
          bool const onLeft = std::abs( position[0] + 1 ) < 1e-12;
          bool const onBottom = std::abs( position[2] - ( vertex ? -1 : 0 ) ) < 1e-12;
          left->InsertNextValue( onLeft );
          bottom->InsertNextValue( onBottom );
          corner->InsertNextValue( onLeft && onBottom && std::abs( position[1] + 1 ) < 1e-12 );
          unusedInteger->InsertNextValue( static_cast< int >( p ) );
          unusedString->InsertNextValue( "point" + std::to_string( p ) );
        }
        for( vtkAbstractArray * array : { static_cast< vtkAbstractArray * >( left ),
                                          static_cast< vtkAbstractArray * >( bottom ),
                                          static_cast< vtkAbstractArray * >( corner ),
                                          static_cast< vtkAbstractArray * >( unusedInteger ),
                                          static_cast< vtkAbstractArray * >( unusedString ) } )
          grid->GetPointData()->AddArray( array );
        vtkNew< vtkTypeInt64Array > cellIds;
        cellIds->SetName( "cellIds" );
        for( vtkIdType c = 0; c < grid->GetNumberOfCells(); ++c ) cellIds->InsertNextValue( idBase + 20000 + 13 * c );
        grid->GetCellData()->SetGlobalIds( cellIds );
        auto const stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        m_directory = std::filesystem::temp_directory_path() / ( "geos-refined-nodal-" + std::to_string( stamp ) );
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

template< typename MANAGER >
void verifyMapsAndMaximum( MANAGER const & manager )
{
  globalIndex maximum = -1;
  auto const ids = manager.localToGlobalMap();
  for( localIndex i = 0; i < manager.size(); ++i )
  {
    EXPECT_EQ( manager.globalToLocalMap( ids[i] ), i );
    maximum = std::max( maximum, ids[i] );
  }
  EXPECT_EQ( manager.maxGlobalIndex(), MpiWrapper::max( maximum ) );
}

template< typename MANAGER >
void verifyEntityIds( MANAGER const & manager, NodeManager const & nodes, int arity )
{
  verifyMapsAndMaximum( manager );
  array1d< globalIndex > local;
  auto const nodeIds = nodes.localToGlobalMap();
  auto const entityIds = manager.localToGlobalMap();
  for( localIndex entity = 0; entity < manager.size(); ++entity )
  {
    TopologyKey key;
    auto const entityNodes = manager.nodeList()[entity];
    for( localIndex p = 0; p < entityNodes.size(); ++p )
      key.push_back( nodeIds[entityNodes[p]] );
    EXPECT_EQ( key.size(), arity );
    if( key.size() != std::size_t( arity ) )
      continue;
    std::sort( key.begin(), key.end() );
    local.emplace_back( entityIds[entity] );
    for( globalIndex point : key )
      local.emplace_back( point );
  }
  array1d< globalIndex > all;
  MpiWrapper::allGatherv( local.toViewConst(), all );
  std::map< globalIndex, TopologyKey > idToKey;
  std::map< TopologyKey, globalIndex > keyToId;
  for( localIndex offset = 0; offset < all.size(); offset += arity + 1 )
  {
    globalIndex const id = all[offset];
    TopologyKey key;
    for( int p = 0; p < arity; ++p )
      key.push_back( all[offset + 1 + p] );
    auto const byId = idToKey.emplace( id, key );
    EXPECT_EQ( byId.first->second, key );
    auto const byKey = keyToId.emplace( key, id );
    EXPECT_EQ( byKey.first->second, id );
  }
}

std::vector< std::set< TopologyKey > > gatherRankTopology( std::set< TopologyKey > const & keys, int arity )
{
  array1d< globalIndex > local;
  for( auto const & key : keys )
  {
    EXPECT_EQ( key.size(), arity );
    local.emplace_back( MpiWrapper::commRank() );
    for( auto point : key )
      local.emplace_back( point );
  }
  array1d< globalIndex > all;
  MpiWrapper::allGatherv( local.toViewConst(), all );
  std::vector< std::set< TopologyKey > > result( MpiWrapper::commSize() );
  for( localIndex offset = 0; offset < all.size(); offset += arity + 1 )
  {
    TopologyKey key;
    for( int p = 0; p < arity; ++p )
      key.push_back( all[offset + 1 + p] );
    result.at( all[offset] ).insert( key );
  }
  return result;
}

std::size_t sharedCount( std::set< TopologyKey > const & a, std::set< TopologyKey > const & b )
{
  std::vector< TopologyKey > common;
  std::set_intersection( a.begin(), a.end(), b.begin(), b.end(), std::back_inserter( common ) );
  return common.size();
}

void initialize( FixtureFile const & fixture, bool vertex, int copies, int level )
{
  GeosxState state( std::make_unique< CommandLineOptions >( commandLineOptions ) );
  auto & problem = state.getProblemManager();
  problem.parseInputString( GEOS_FMT(
                              R"xml(<Problem>
    <Mesh><VTKMesh name="mesh" file="{}" useGlobalIds="1" scatterMethod="rcb"
      partitionRefinement="0" uniformRefinement="{}" nodesetNames="{{leftMask,bottomMask,cornerMask}}" /></Mesh>
    <Solvers gravityVector="{{0,0,0}}"><SolidMechanicsLagrangianFEM name="solid"
      discretization="FE1" targetRegions="{{volume}}">
      <LinearSolverParameters solverType="gmres" preconditionerType="amg"/>
    </SolidMechanicsLagrangianFEM></Solvers>
    <NumericalMethods><FiniteElements><FiniteElementSpace name="FE1" order="1"/></FiniteElements></NumericalMethods>
    <ElementRegions><CellElementRegion name="volume" cellBlocks="{{*}}" materialList="{{rock}}"/></ElementRegions>
    <Constitutive><ElasticIsotropic name="rock" defaultDensity="2700"
      defaultBulkModulus="5e9" defaultShearModulus="4e9"/></Constitutive>
    </Problem>)xml", fixture.path(), level ) );
  problem.problemSetup();
  auto & domain = problem.getDomainPartition();
  auto & mesh = domain.getMeshBody( "mesh" ).getBaseDiscretization();
  auto & nodes = mesh.getNodeManager();
  auto const & faces = mesh.getFaceManager();
  auto const & edges = mesh.getEdgeManager();
  auto & solver = problem.getPhysicsSolverManager().getGroup< SolidMechanicsLagrangianFEM >( "solid" );
  solver.setupSystem( domain, solver.getDofManager(), solver.getLocalMatrix(), solver.getSystemRhs(), solver.getSystemSolution() );
  auto const & dofs = solver.getDofManager();
  auto const dofNumbers = nodes.getReference< array1d< globalIndex > >( dofs.getKey( "totalDisplacement" ) ).toViewConst();
  auto const nodeIds = nodes.localToGlobalMap();
  auto const ghosts = nodes.ghostRank();
  int const rank = MpiWrapper::commRank();
  int const ranks = MpiWrapper::commSize();
  int const subdivisions = 1 << level;
  int const groups = ( vertex ? 8 : 4 ) * copies;
  int const fillers = std::max( ranks - groups, 0 );
  globalIndex const expectedNodes = copies * ( 2 * subdivisions + 1 ) * ( 2 * subdivisions + 1 ) *
                                    ( vertex ? 2 * subdivisions + 1 : subdivisions + 1 ) +
                                    fillers * ( subdivisions + 1 ) * ( subdivisions + 1 ) * ( subdivisions + 1 );
  EXPECT_EQ( dofs.numGlobalDofs(), 3 * expectedNodes );
  verifyMapsAndMaximum( nodes );
  verifyEntityIds( faces, nodes, 4 );
  verifyEntityIds( edges, nodes, 2 );

  array1d< globalIndex > owned;
  array1d< globalIndex > ownedPositions;
  auto const positions = nodes.referencePosition();
  auto const & left = nodes.sets().getReference< SortedArray< localIndex > >( "leftMask" );
  auto const & bottom = nodes.sets().getReference< SortedArray< localIndex > >( "bottomMask" );
  auto const & corner = nodes.sets().getReference< SortedArray< localIndex > >( "cornerMask" );
  for( localIndex point = 0; point < nodes.size(); ++point )
  {
    bool const onLeft = std::abs( positions[point][0] + 1 ) < 1e-12;
    bool const onBottom = std::abs( positions[point][2] - ( vertex ? -1 : 0 ) ) < 1e-12;
    EXPECT_EQ( left.contains( point ), onLeft );
    EXPECT_EQ( bottom.contains( point ), onBottom );
    EXPECT_EQ( corner.contains( point ), onLeft && onBottom && std::abs( positions[point][1] + 1 ) < 1e-12 );
    EXPECT_GE( dofNumbers[point], 0 );
    EXPECT_EQ( dofNumbers[point] % 3, 0 );
    EXPECT_LT( dofNumbers[point] + 2, dofs.numGlobalDofs() );
    if( ghosts[point] < 0 )
    {
      owned.emplace_back( nodeIds[point] );
      owned.emplace_back( dofNumbers[point] );
      owned.emplace_back( rank );
      ownedPositions.emplace_back( nodeIds[point] );
      for( int d = 0; d < 3; ++d )
      {
        real64 const scaled = positions[point][d] * subdivisions;
        globalIndex const coordinate = std::llround( scaled );
        EXPECT_DOUBLE_EQ( scaled, real64( coordinate ) );
        ownedPositions.emplace_back( coordinate );
      }
    }
  }
  array1d< globalIndex > allOwned;
  MpiWrapper::allGatherv( owned.toViewConst(), allOwned );
  array1d< globalIndex > allPositions;
  MpiWrapper::allGatherv( ownedPositions.toViewConst(), allPositions );
  std::map< std::array< globalIndex, 3 >, std::set< globalIndex > > coordinateIds;
  for( localIndex i = 0; i < allPositions.size(); i += 4 )
    coordinateIds[{ allPositions[i + 1], allPositions[i + 2], allPositions[i + 3] }].insert( allPositions[i] );
  for( int z = vertex ? -subdivisions : 0; z <= subdivisions; ++z )
    for( int y = -subdivisions; y <= subdivisions; ++y )
      for( int x = -subdivisions; x <= subdivisions; ++x )
        EXPECT_EQ( ( coordinateIds[{ x, y, z }].size() ), copies );
  EXPECT_EQ( allOwned.size(), 3 * expectedNodes );
  std::map< globalIndex, std::pair< globalIndex, int > > authority;
  std::set< globalIndex > globalDofs;
  for( localIndex i = 0; i < allOwned.size(); i += 3 )
  {
    EXPECT_TRUE( authority.emplace( allOwned[i], std::make_pair( allOwned[i + 1], allOwned[i + 2] ) ).second );
    for( int component = 0; component < 3; ++component )
      EXPECT_TRUE( globalDofs.insert( allOwned[i + 1] + component ).second );
  }
  EXPECT_EQ( globalDofs.size(), 3 * expectedNodes );
  if( !globalDofs.empty() )
  {
    EXPECT_EQ( *globalDofs.begin(), 0 );
    EXPECT_EQ( *globalDofs.rbegin(), dofs.numGlobalDofs() - 1 );
  }
  for( localIndex point = 0; point < nodes.size(); ++point )
  {
    auto const found = authority.find( nodeIds[point] );
    if( found == authority.end() )
    {
      ADD_FAILURE() << "Missing nodal owner"; continue;
    }
    EXPECT_EQ( dofNumbers[point], found->second.first );
    EXPECT_EQ( ghosts[point] < 0 ? rank : ghosts[point], found->second.second );
  }

  std::set< TopologyKey > ownedNodes, ownedEdges, ownedFaces;
  localIndex ownedCells = 0;
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cells )
  {
    verifyMapsAndMaximum( cells );
    for( localIndex c = 0; c < cells.size(); ++c )
      if( cells.ghostRank()[c] < 0 )
      {
        ++ownedCells;
        EXPECT_GT( cells.getElementVolume()[c], 0 );
        // Cell-map rows can be strided in accelerator layouts. Index each
        // slice instead of requiring LvArray's contiguous iterator interface.
        auto const cellNodes = cells.nodeList()[c];
        for( localIndex p = 0; p < cellNodes.size(); ++p )
          ownedNodes.insert( { nodeIds[cellNodes[p]] } );
        auto const cellEdges = cells.edgeList()[c];
        for( localIndex e = 0; e < cellEdges.size(); ++e )
        {
          TopologyKey key;
          auto const edgeNodes = edges.nodeList()[cellEdges[e]];
          for( localIndex p = 0; p < edgeNodes.size(); ++p )
            key.push_back( nodeIds[edgeNodes[p]] );
          std::sort( key.begin(), key.end() );
          ownedEdges.insert( key );
        }
        auto const cellFaces = cells.faceList()[c];
        for( localIndex f = 0; f < cellFaces.size(); ++f )
        {
          TopologyKey key;
          auto const faceNodes = faces.nodeList()[cellFaces[f]];
          for( localIndex p = 0; p < faceNodes.size(); ++p )
            key.push_back( nodeIds[faceNodes[p]] );
          std::sort( key.begin(), key.end() );
          ownedFaces.insert( key );
        }
      }
  } );
  EXPECT_EQ( MpiWrapper::sum( ownedCells ), std::max( groups, ranks ) * subdivisions * subdivisions * subdivisions );
  auto const rankNodes = gatherRankTopology( ownedNodes, 1 );
  auto const rankEdges = gatherRankTopology( ownedEdges, 2 );
  auto const rankFaces = gatherRankTopology( ownedFaces, 4 );
  bool foundEdgeOnly = false, foundVertexOnly = false;
  for( int a = 0; a < ranks; ++a )
    for( int b = a + 1; b < ranks; ++b )
    {
      auto const sharedNodes = sharedCount( rankNodes[a], rankNodes[b] );
      auto const sharedEdges = sharedCount( rankEdges[a], rankEdges[b] );
      auto const sharedFaces = sharedCount( rankFaces[a], rankFaces[b] );
      foundEdgeOnly = foundEdgeOnly || ( sharedNodes > 1 && sharedEdges > 0 && sharedFaces == 0 );
      foundVertexOnly = foundVertexOnly || ( sharedNodes == 1 && sharedEdges == 0 && sharedFaces == 0 );
    }
  if( ranks >= groups )
  {
    EXPECT_TRUE( foundEdgeOnly );
    if( vertex )
    {
      EXPECT_TRUE( foundVertexOnly );
    }
  }

  // Deliberately corrupt ghost values, then use the normal solver communication
  // path. A position-only merge would also fail this ID-specific value oracle.
  auto const displacement = nodes.getField< fields::solidMechanics::totalDisplacement >().toView();
  auto probe = []( globalIndex id, int owner, int component )
  { return real64( id % 65537 + 100000 * owner + 3 * component ); };
  for( localIndex point = 0; point < nodes.size(); ++point )
    for( int component = 0; component < 3; ++component )
      displacement[point][component] = ghosts[point] < 0 ? probe( nodeIds[point], rank, component ) : -1000;
  FieldIdentifiers fieldsToSync;
  fieldsToSync.addFields( FieldLocation::Node, { "totalDisplacement" } );
  CommunicationTools::getInstance().synchronizeFields( fieldsToSync, mesh, domain.getNeighbors(), false );
  for( localIndex point = 0; point < nodes.size(); ++point )
    for( int component = 0; component < 3; ++component )
      EXPECT_DOUBLE_EQ( displacement[point][component], probe( nodeIds[point], ghosts[point] < 0 ? rank : ghosts[point], component ) );
}
} // namespace

TEST( VTKRefinedNodalInitialization, EdgeOnlyNeighborsAtThreeLevels )
{
  FixtureFile fixture( false, 1 );
  for( int level : { 0, 1, 2 } )
  {
    SCOPED_TRACE( "edge level " + std::to_string( level ) );
    initialize( fixture, false, 1, level );
  }
}

TEST( VTKRefinedNodalInitialization, VertexOnlyNeighborsAtThreeLevels )
{
  FixtureFile fixture( true, 1 );
  for( int level : { 0, 1, 2 } )
  {
    SCOPED_TRACE( "vertex level " + std::to_string( level ) );
    initialize( fixture, true, 1, level );
  }
}

TEST( VTKRefinedNodalInitialization, UnrelatedCoincidentComponentsRetainSeparateDofs )
{
  FixtureFile fixture( false, 2 );
  for( int level : { 0, 1, 2 } )
  {
    SCOPED_TRACE( "coincident components level " + std::to_string( level ) );
    initialize( fixture, false, 2, level );
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
