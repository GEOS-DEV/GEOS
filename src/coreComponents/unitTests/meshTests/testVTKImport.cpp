/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlockABC.hpp"

// special CMake-generated include
#include "tests/meshDirName.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::dataRepository;

template< class V >
void TestMeshImport( string const & meshFilePath, V const & validate )
{
  string const meshNode = "<Mesh ><VTKMesh name=\"mesh\" file=\"" + meshFilePath + "\" /></Mesh>";
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( meshNode.c_str(), meshNode.size() );
  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child( "Mesh" );

  conduit::Node node;
  Group root( "root", node );

  MeshManager meshManager( "mesh", &root );
  meshManager.processInputFileRecursive( xmlMeshNode );
  meshManager.postProcessInputRecursive();
  DomainPartition domain( "domain", &root );
  meshManager.generateMeshes( domain );

  // TODO Field import is not tested yet. Proper refactoring needs to be done first.

  validate( domain.getMeshBody( "mesh" ).getGroup< CellBlockManagerABC >( keys::cellManager ) );
}

TEST( VTKImport, cube )
{
  auto validate = []( CellBlockManagerABC const & cellBlockManager ) -> void
  {
    // `cube.vtk` is a cube made by 3 x 3 x 3 = 27 Hexahedron 3d-elements.
    // It contains 4 x 4 x 4 = 64 nodes.
    // On each face of the cube, you have 3 x 3 quad faces. Hence 9 x 6 = 54 quad 2d-elements.
    // On each edge of the cube, you have 3 Line elements. Hence 3 x 12 = 36 line 1d-elements.
    // On each vertex of the cube you have on Vertex element. Hence 8 Vertex 0d-cells.

    // The `cube.vtk` mesh contains an "attribute" field that is used to group cells together.
    // A region with id `-1` is considered as a non region.
    // For testing purpose, the "attribute" field was designed such that
    // - All 36 `Line` elements are in "region" 1 except the last two in regions -1 and 9.
    // - All 36 `Quad` elements are in "region" 2 except 4.
    //   Counting backwards from the end, quads number 0, 1, 3 and 4 with respectively regions 9 and -1, -1, -1.
    //   Those quads were selected such that they form a larger square, excluding the central node (number 55) from the region 2.
    //   This should appear in the test.
    // - All 36 `Hexahedron` elements are in "region" 3 except the last two in regions -1 and 9.
    // - All 36 `Vertex` elements are in "region" 4 except the last two in regions -1 and 9.

    // When run in parallel with two MPI ranks, the central hexahedra are on the splitting boundary.
    // The VTK default pattern is to assign the cell to one unique rank.
    // It happens to be the first one in our case.
    // This way, rank 0 (lower `x`) gets 18 hexaedra and 48 nodes,
    // while rank 1 (greater `x`) gets 9 hexahedra and 32 nodes.

    // This pattern could be influenced by settings the parameters of
    // vtkRedistributeDataSetFilter::SetBoundaryMode(...) to
    // ASSIGN_TO_ALL_INTERSECTING_REGIONS, ASSIGN_TO_ONE_REGION or SPLIT_BOUNDARY_CELLS.
    localIndex const expectedNumNodes = expected( 64, { 48, 32 } );
    ASSERT_EQ( cellBlockManager.numNodes(), expectedNumNodes );
    ASSERT_EQ( cellBlockManager.numEdges(), expected( 144, { 104, 64 } ) );
    ASSERT_EQ( cellBlockManager.numFaces(), expected( 108, { 75, 42 } ) );

    // The information in the tables is not filled yet. We can check the consistency of the sizes.
    ASSERT_EQ( cellBlockManager.getNodeToFaces().size(), expectedNumNodes );
    ASSERT_EQ( cellBlockManager.getNodeToElements().size(), expectedNumNodes );

    // We have all the 4 x 4  x 4 = 64 nodes in the "all" set.
    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), expectedNumNodes );

    // The "2" set are all the boundary nodes (64 - 8 inside nodes = 56),
    // minus an extra node that belongs to regions -1 and 9 only.
    SortedArray< localIndex > const & nodesRegion2 = cellBlockManager.getNodeSets().at( "2" );
    ASSERT_EQ( nodesRegion2.size(), expected( 55, { 39, 27 } ) );

    // Region "9" has only one quad, on the greater `x` direction.
    // This hex will belong to MPI rank 1.
    SortedArray< localIndex > const & nodesRegion9 = cellBlockManager.getNodeSets().at( "9" );
    ASSERT_EQ( nodesRegion9.size(), expected( 4, { 0, 4 } ) );

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    // 1 elements type on 3 regions ("-1", "3", "9") = 3 sub-groups
    std::array< std::pair< string, int >, 3 > const expectedCellBlocks =
    {
      {
        { "hexahedra", expected( 1, {  1, 0 } ) },
        { "3_hexahedra", expected( 25, { 17, 8 } ) },
        { "9_hexahedra", expected( 1, {  0, 1 } ) }
      }
    };
    ASSERT_EQ( cellBlockManager.getCellBlocks().numSubGroups(), expectedCellBlocks.size() );

    for( const auto & nameAndSize : expectedCellBlocks )
    {
      ASSERT_TRUE( cellBlockManager.getCellBlocks().hasGroup< CellBlockABC >( nameAndSize.first ) );
      CellBlockABC const * h = &cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( nameAndSize.first );
      localIndex const expectedSize = nameAndSize.second;

      // 8 nodes, 12 edges and 6 faces per hex.
      ASSERT_EQ( h->getElemToNodes().size( 1 ), 8 );
      ASSERT_EQ( h->getElemToEdges().size( 1 ), 12 );
      ASSERT_EQ( h->getElemToFaces().size( 1 ), 6 );

      ASSERT_EQ( h->size(), expectedSize );
      ASSERT_EQ( h->getElemToNodes().size( 0 ), expectedSize );
      ASSERT_EQ( h->getElemToEdges().size( 0 ), expectedSize );
      ASSERT_EQ( h->getElemToFaces().size( 0 ), expectedSize );
    }
  };

  string const cubeVTK = testMeshDir + "/cube.vtk";
  string const cubeVTU = testMeshDir + "/cube.vtu";
//  string const cubePVTU = testMeshDir + "/cube.pvtu";

  TestMeshImport( cubeVTK, validate );
  TestMeshImport( cubeVTU, validate );
//  TestMeshImport( cubePVTU, validate );
}

TEST( VTKImport, medley )
{
  SKIP_TEST_IN_PARALLEL( "Neither relevant nor implemented in parallel" );

  auto validate = []( CellBlockManagerABC const & cellBlockManager ) -> void
  {
    // `medley.vtk` is made of four elements.
    // - Element 0 is a pyramid, in region 0.
    // - Element 1 is an hexahedron, in region 1.
    // - Element 2 is a wedge, in region 2.
    // - Element 3 is a tetrahedron, in region 3.
    // All the elements belong to a region. Therefore, there is no "-1" region.
    // It contains 12 nodes, 24 edges, 17 faces.

    ASSERT_EQ( cellBlockManager.numNodes(), 12 );
    ASSERT_EQ( cellBlockManager.numEdges(), 24 );
    ASSERT_EQ( cellBlockManager.numFaces(), 17 );

    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), 12 );

    // 4 elements types x 4 regions = 16 sub-groups
    ASSERT_EQ( cellBlockManager.getCellBlocks().numSubGroups(), 16 );

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    CellBlockABC const & zone0 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "0_pyramids" );
    CellBlockABC const & zone1 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "1_hexahedra" );
    CellBlockABC const & zone2 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "2_wedges" );
    CellBlockABC const & zone3 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "3_tetrahedra" );

    std::vector< string > const elementNames{ "pyramids", "hexahedra", "wedges", "tetrahedra" };
    for( std::size_t prefix: { 0, 1, 2, 3 } )
    {
      for( std::size_t i = 0; i < 4; ++i )
      {
        string const name = std::to_string( prefix ) + "_" + elementNames[i];
        CellBlockABC const & zone = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( name );
        ASSERT_EQ( zone.size(), prefix == i ? 1 : 0 );
      }
    }

    // Pyramid
    ASSERT_EQ( zone0.getElemToNodes().size( 1 ), 5 );
    ASSERT_EQ( zone0.getElemToEdges().size( 1 ), 8 );
    ASSERT_EQ( zone0.getElemToFaces().size( 1 ), 5 );

    // Hexahedron
    ASSERT_EQ( zone1.getElemToNodes().size( 1 ), 8 );
    ASSERT_EQ( zone1.getElemToEdges().size( 1 ), 12 );
    ASSERT_EQ( zone1.getElemToFaces().size( 1 ), 6 );

    // Wedges
    ASSERT_EQ( zone2.getElemToNodes().size( 1 ), 6 );
    ASSERT_EQ( zone2.getElemToEdges().size( 1 ), 9 );
    ASSERT_EQ( zone2.getElemToFaces().size( 1 ), 5 );

    // Tetrahedron
    ASSERT_EQ( zone3.getElemToNodes().size( 1 ), 4 );
    ASSERT_EQ( zone3.getElemToEdges().size( 1 ), 6 );
    ASSERT_EQ( zone3.getElemToFaces().size( 1 ), 4 );

    for( auto const & z: { &zone0, &zone1, &zone2, &zone3 } )
    {
      ASSERT_EQ( z->size(), 1 );
      ASSERT_EQ( z->getElemToNodes().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToEdges().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToFaces().size( 0 ), 1 );
    }
  };

  string const medleyVTK = testMeshDir + "/medley.vtk";

  TestMeshImport( medleyVTK, validate );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
