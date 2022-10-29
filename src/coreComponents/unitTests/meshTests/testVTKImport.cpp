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
  string const meshNode = GEOSX_FMT( R"(<Mesh><VTKMesh name="mesh" file="{}" partitionRefinement="0"/></Mesh>)", meshFilePath );
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
    // For example, if it happens to be the first one:
    // - rank 0 (lower `x`) gets 18 hexaedra and 48 nodes,
    // - rank 1 (greater `x`) gets 9 hexahedra and 32 nodes.

    // This pattern could be influenced by settings the parameters of
    // vtkRedistributeDataSetFilter::SetBoundaryMode(...) to
    // ASSIGN_TO_ALL_INTERSECTING_REGIONS, ASSIGN_TO_ONE_REGION or SPLIT_BOUNDARY_CELLS.
    localIndex const expectedNumNodesRank1 = expected( 64, { 48, 32 } );
    localIndex const expectedNumNodesRank2 = expected( 64, { 32, 48 } );
    bool rankswap = cellBlockManager.numNodes() == expectedNumNodesRank2;
    localIndex const expectedNumNodes = rankswap ? expectedNumNodesRank2 : expectedNumNodesRank1;
    auto expectedSwap = [=] ( int seq, std::initializer_list< int > par )
    {
      std::vector< int > tmp( par );
      if( rankswap )
        return expected( seq, { tmp[1], tmp[0] } );
      else
        return expected( seq, par );
    };

    ASSERT_EQ( cellBlockManager.numNodes(), expectedNumNodes );
    ASSERT_EQ( cellBlockManager.numEdges(), expectedSwap( 144, { 104, 64 } ) );
    ASSERT_EQ( cellBlockManager.numFaces(), expectedSwap( 108, { 75, 42 } ) );

    // The information in the tables is not filled yet. We can check the consistency of the sizes.
    ASSERT_EQ( cellBlockManager.getNodeToFaces().size(), expectedNumNodes );
    ASSERT_EQ( cellBlockManager.getNodeToElements().toCellIndex.size(), expectedNumNodes );

    // We have all the 4 x 4  x 4 = 64 nodes in the "all" set.
    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), expectedNumNodes );

    if( cellBlockManager.getNodeSets().size()>1 )
    {
      // The "2" set are all the boundary nodes (64 - 8 inside nodes = 56),
      // minus an extra node that belongs to regions -1 and 9 only.
      SortedArray< localIndex > const & nodesRegion2 = cellBlockManager.getNodeSets().at( "2" );
      ASSERT_EQ( nodesRegion2.size(), expectedSwap( 55, { 39, 27 } ) );

      // Region "9" has only one quad, on the greater `x` direction.
      // This hex will belong to MPI rank 1.
      SortedArray< localIndex > const & nodesRegion9 = cellBlockManager.getNodeSets().at( "9" );
      ASSERT_EQ( nodesRegion9.size(), expectedSwap( 4, { 0, 4 } ) );

      // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
      // 1 elements type on 3 regions ("-1", "3", "9") = 3 sub-groups
      std::array< std::pair< string, int >, 3 > const expectedCellBlocks =
      {
        {
          { "hexahedra", expectedSwap( 1, {  1, 0 } ) },
          { "3_hexahedra", expectedSwap( 25, { 17, 8 } ) },
          { "9_hexahedra", expectedSwap( 1, {  0, 1 } ) }
        }
      };
      ASSERT_EQ( cellBlockManager.getCellBlocks().numSubGroups(), expectedCellBlocks.size() );

      for( const auto & nameAndSize : expectedCellBlocks )
      {
        ASSERT_TRUE( cellBlockManager.getCellBlocks().hasGroup< CellBlockABC >( nameAndSize.first ) );

        // here pb
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
    }
  };

  std::set< string > const meshFiles{ "cube.vtk",
                                      "cube_STRUCTURED_POINTS.vtk",
                                      "cube_RECTILINEAR_GRID.vtk",
                                      "cube_STRUCTURED_GRID.vtk",
                                      "cube_UNSTRUCTURED_GRID.vtk",
                                      "cube.vtu",
                                      //"cube.pvtu",
                                      "cube.vts",
                                      "cube.pvts",
                                      "cube.vtr",
                                      "cube.pvtr",
                                      "cube.vti",
                                      "cube.pvti" };
  for( string const & meshFile: meshFiles )
  {
    TestMeshImport( testMeshDir + "/" + meshFile, validate );
  }

}

TEST( VTKImport, medley )
{
  SKIP_TEST_IN_PARALLEL( "Neither relevant nor implemented in parallel" );

  auto validate = []( CellBlockManagerABC const & cellBlockManager ) -> void
  {
    // `medley.vtk` is made of eleven elements.
    // - Element 0 is a pyramid, in region 0.
    // - Element 1 is an hexahedron, in region 1.
    // - Element 2 is a wedge, in region 2.
    // - Element 3 is a tetrahedron, in region 3.
    // - Element 4 is a pentagonal prism, in region 4.
    // - Element 5 is an hexagonal prism, in region 5.
    // - Element 6 is an heptagonal prism, in region 6.
    // - Element 7 is an octagonal prism, in region 7.
    // - Element 8 is an nonagonal prism, in region 8.
    // - Element 9 is an decagonal prism, in region 9.
    // - Element 10 is an hendecagonal prism, in region 10.
    // All the elements belong to a region. Therefore, there is no "-1" region.
    // It contains 112 nodes, 180 edges, 84 faces.

    ASSERT_EQ( cellBlockManager.numNodes(), 112 );
    ASSERT_EQ( cellBlockManager.numEdges(), 180 );
    ASSERT_EQ( cellBlockManager.numFaces(), 84 );

    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), 112 );

    // 11 elements types x 11 regions = 121 sub-groups
    ASSERT_EQ( cellBlockManager.getCellBlocks().numSubGroups(), 121 );

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    CellBlockABC const & zone0 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "0_pyramids" );
    CellBlockABC const & zone1 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "1_hexahedra" );
    CellBlockABC const & zone2 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "2_wedges" );
    CellBlockABC const & zone3 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "3_tetrahedra" );
    CellBlockABC const & zone4 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "4_pentagonalPrisms" );
    CellBlockABC const & zone5 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "5_hexagonalPrisms" );
    CellBlockABC const & zone6 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "6_heptagonalPrisms" );
    CellBlockABC const & zone7 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "7_octagonalPrisms" );
    CellBlockABC const & zone8 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "8_nonagonalPrisms" );
    CellBlockABC const & zone9 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "9_decagonalPrisms" );
    CellBlockABC const & zone10 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "10_hendecagonalPrisms" );

    std::vector< string > const elementNames{ "pyramids",
                                              "hexahedra",
                                              "wedges",
                                              "tetrahedra",
                                              "pentagonalPrisms",
                                              "hexagonalPrisms",
                                              "heptagonalPrisms",
                                              "octagonalPrisms",
                                              "nonagonalPrisms",
                                              "decagonalPrisms",
                                              "hendecagonalPrisms" };
    for( std::size_t prefix: { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 } )
    {
      for( std::size_t i = 0; i < 11; ++i )
      {
        string const name = std::to_string( prefix ) + "_" + elementNames[i];
        CellBlockABC const & zone = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( name );
        ASSERT_EQ( zone.size(), prefix == i ? 1 : 0 );
      }
    }

    // Pyramid
    auto elementToNodes = zone0.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 5 );
    ASSERT_EQ( zone0.getElemToEdges().size( 1 ), 8 );
    ASSERT_EQ( zone0.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 3 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 0 );

    // Hexahedron
    elementToNodes = zone1.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 8 );
    ASSERT_EQ( zone1.getElemToEdges().size( 1 ), 12 );
    ASSERT_EQ( zone1.getElemToFaces().size( 1 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 3 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 7 );

    // Wedges
    elementToNodes = zone2.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 6 );
    ASSERT_EQ( zone2.getElemToEdges().size( 1 ), 9 );
    ASSERT_EQ( zone2.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 7 );

    // Tetrahedron
    elementToNodes = zone3.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 4 );
    ASSERT_EQ( zone3.getElemToEdges().size( 1 ), 6 );
    ASSERT_EQ( zone3.getElemToFaces().size( 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 11 );

    // Pentagonal prism
    elementToNodes = zone4.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 10 );
    ASSERT_EQ( zone4.getElemToEdges().size( 1 ), 15 );
    ASSERT_EQ( zone4.getElemToFaces().size( 1 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 13 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 14 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 3 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 15 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 16 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 17 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 7 );

    // Hexagonal prism
    elementToNodes = zone5.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 12 );
    ASSERT_EQ( zone5.getElemToEdges().size( 1 ), 18 );
    ASSERT_EQ( zone5.getElemToFaces().size( 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 18 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 19 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 20 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 21 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 22 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 23 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 24 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 25 );

    // Heptagonal prism
    elementToNodes = zone6.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 14 );
    ASSERT_EQ( zone6.getElemToEdges().size( 1 ), 21 );
    ASSERT_EQ( zone6.getElemToFaces().size( 1 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 19 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 26 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 27 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 28 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 29 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 30 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 20 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 23 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 31 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 32 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 33 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 34 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 35 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 24 );

    // Octagonal prism
    elementToNodes = zone7.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 16 );
    ASSERT_EQ( zone7.getElemToEdges().size( 1 ), 24 );
    ASSERT_EQ( zone7.getElemToFaces().size( 1 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 36 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 37 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 38 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 39 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 40 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 41 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 42 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 43 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 44 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 45 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 46 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 47 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 48 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 49 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 50 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 51 );

    // Nonagonal prism
    elementToNodes = zone8.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 18 );
    ASSERT_EQ( zone8.getElemToEdges().size( 1 ), 27 );
    ASSERT_EQ( zone8.getElemToFaces().size( 1 ), 11 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 52 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 53 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 54 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 55 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 56 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 57 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 58 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 59 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 60 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 61 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 62 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 63 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 64 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 65 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 66 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 67 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 68 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 69 );

    // Decagonal prism
    elementToNodes = zone9.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 20 );
    ASSERT_EQ( zone9.getElemToEdges().size( 1 ), 30 );
    ASSERT_EQ( zone9.getElemToFaces().size( 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 70 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 71 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 72 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 73 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 74 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 75 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 76 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 77 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 78 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 79 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 80 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 81 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 82 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 83 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 84 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 85 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 86 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 87 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 88 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 89 );

    // Hendecagonal prism
    elementToNodes = zone10.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 22 );
    ASSERT_EQ( zone10.getElemToEdges().size( 1 ), 33 );
    ASSERT_EQ( zone10.getElemToFaces().size( 1 ), 13 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 90 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 91 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 92 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 93 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 94 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 95 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 96 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 97 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 98 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 99 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 100 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 101 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 102 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 103 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 104 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 105 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 106 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 107 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 108 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 109 );
    EXPECT_EQ( elementToNodes( 0, 20 ), 110 );
    EXPECT_EQ( elementToNodes( 0, 21 ), 111 );

    for( auto const & z: { &zone0, &zone1, &zone2, &zone3, &zone4, &zone5, &zone6, &zone7, &zone8, &zone9, &zone10 } )
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

TEST( VTKImport, medley42 )
{
  SKIP_TEST_IN_PARALLEL( "Neither relevant nor implemented in parallel" );

  auto validate = []( CellBlockManagerABC const & cellBlockManager ) -> void
  {
    // `medley-42.vtk` is the same as `medley.vtk` with all eleven elements defined as polyhedron.
    // All the elements belong to a region. Therefore, there is no "-1" region.
    // It contains 112 nodes, 180 edges, 84 faces.

    ASSERT_EQ( cellBlockManager.numNodes(), 112 );
    ASSERT_EQ( cellBlockManager.numEdges(), 180 );
    ASSERT_EQ( cellBlockManager.numFaces(), 84 );

    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), 112 );

    // 11 elements types x 11 regions = 121 sub-groups
    ASSERT_EQ( cellBlockManager.getCellBlocks().numSubGroups(), 121 );

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    CellBlockABC const & zone0 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "0_pyramids" );
    CellBlockABC const & zone1 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "1_hexahedra" );
    CellBlockABC const & zone2 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "2_wedges" );
    CellBlockABC const & zone3 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "3_tetrahedra" );
    CellBlockABC const & zone4 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "4_pentagonalPrisms" );
    CellBlockABC const & zone5 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "5_hexagonalPrisms" );
    CellBlockABC const & zone6 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "6_heptagonalPrisms" );
    CellBlockABC const & zone7 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "7_octagonalPrisms" );
    CellBlockABC const & zone8 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "8_nonagonalPrisms" );
    CellBlockABC const & zone9 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "9_decagonalPrisms" );
    CellBlockABC const & zone10 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "10_hendecagonalPrisms" );

    std::vector< string > const elementNames{ "pyramids",
                                              "hexahedra",
                                              "wedges",
                                              "tetrahedra",
                                              "pentagonalPrisms",
                                              "hexagonalPrisms",
                                              "heptagonalPrisms",
                                              "octagonalPrisms",
                                              "nonagonalPrisms",
                                              "decagonalPrisms",
                                              "hendecagonalPrisms" };
    for( std::size_t prefix: { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 } )
    {
      for( std::size_t i = 0; i < 11; ++i )
      {
        string const name = std::to_string( prefix ) + "_" + elementNames[i];
        CellBlockABC const & zone = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( name );
        ASSERT_EQ( zone.size(), prefix == i ? 1 : 0 );
      }
    }

    // Pyramid
    auto elementToNodes = zone0.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 5 );
    ASSERT_EQ( zone0.getElemToEdges().size( 1 ), 8 );
    ASSERT_EQ( zone0.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 3 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 0 );

    // Hexahedron
    elementToNodes = zone1.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 8 );
    ASSERT_EQ( zone1.getElemToEdges().size( 1 ), 12 );
    ASSERT_EQ( zone1.getElemToFaces().size( 1 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 3 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 7 );

    // Wedges
    elementToNodes = zone2.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 6 );
    ASSERT_EQ( zone2.getElemToEdges().size( 1 ), 9 );
    ASSERT_EQ( zone2.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 7 );

    // Tetrahedron
    elementToNodes = zone3.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 4 );
    ASSERT_EQ( zone3.getElemToEdges().size( 1 ), 6 );
    ASSERT_EQ( zone3.getElemToFaces().size( 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 11 );

    // Pentagonal prism
    elementToNodes = zone4.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 10 );
    ASSERT_EQ( zone4.getElemToEdges().size( 1 ), 15 );
    ASSERT_EQ( zone4.getElemToFaces().size( 1 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 13 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 14 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 3 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 15 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 16 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 17 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 7 );

    // Hexagonal prism
    elementToNodes = zone5.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 12 );
    ASSERT_EQ( zone5.getElemToEdges().size( 1 ), 18 );
    ASSERT_EQ( zone5.getElemToFaces().size( 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 18 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 19 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 20 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 21 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 22 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 23 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 24 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 25 );

    // Heptagonal prism
    elementToNodes = zone6.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 14 );
    ASSERT_EQ( zone6.getElemToEdges().size( 1 ), 21 );
    ASSERT_EQ( zone6.getElemToFaces().size( 1 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 19 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 26 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 27 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 28 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 29 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 30 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 20 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 23 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 31 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 32 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 33 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 34 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 35 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 24 );

    // Octagonal prism
    elementToNodes = zone7.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 16 );
    ASSERT_EQ( zone7.getElemToEdges().size( 1 ), 24 );
    ASSERT_EQ( zone7.getElemToFaces().size( 1 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 36 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 37 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 38 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 39 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 40 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 41 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 42 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 43 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 44 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 45 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 46 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 47 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 48 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 49 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 50 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 51 );

    // Nonagonal prism
    elementToNodes = zone8.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 18 );
    ASSERT_EQ( zone8.getElemToEdges().size( 1 ), 27 );
    ASSERT_EQ( zone8.getElemToFaces().size( 1 ), 11 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 52 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 53 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 54 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 55 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 56 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 57 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 58 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 59 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 60 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 61 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 62 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 63 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 64 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 65 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 66 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 67 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 68 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 69 );

    // Decagonal prism
    elementToNodes = zone9.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 20 );
    ASSERT_EQ( zone9.getElemToEdges().size( 1 ), 30 );
    ASSERT_EQ( zone9.getElemToFaces().size( 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 70 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 71 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 72 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 73 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 74 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 75 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 76 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 77 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 78 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 79 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 80 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 81 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 82 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 83 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 84 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 85 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 86 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 87 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 88 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 89 );

    // Hendecagonal prism
    elementToNodes = zone10.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 22 );
    ASSERT_EQ( zone10.getElemToEdges().size( 1 ), 33 );
    ASSERT_EQ( zone10.getElemToFaces().size( 1 ), 13 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 90 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 91 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 92 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 93 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 94 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 95 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 96 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 97 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 98 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 99 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 100 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 101 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 102 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 103 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 104 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 105 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 106 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 107 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 108 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 109 );
    EXPECT_EQ( elementToNodes( 0, 20 ), 110 );
    EXPECT_EQ( elementToNodes( 0, 21 ), 111 );

    for( auto const & z: { &zone0, &zone1, &zone2, &zone3, &zone4, &zone5, &zone6, &zone7, &zone8, &zone9, &zone10 } )
    {
      ASSERT_EQ( z->size(), 1 );
      ASSERT_EQ( z->getElemToNodes().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToEdges().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToFaces().size( 0 ), 1 );
    }
  };

  string const medleyVTK = testMeshDir + "/medley-42.vtk";

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
