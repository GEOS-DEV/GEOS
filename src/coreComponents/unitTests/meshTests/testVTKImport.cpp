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
    // `supportedElements.vtk` is made of eleven elements.
    // - Element 0 is a tetrahedron, in region 0.
    // - Element 1 is a pyramid, in region 1.
    // - Element 2 is a wedge, in region 2.
    // - Element 3 is a hexahedron, in region 3.
    // - Element 4 is a pentagonal prism, in region 4.
    // - Element 5 is a hexagonal prism, in region 5.
    // - Element 6 is a heptagonal prism, in region 6.
    // - Element 7 is a octagonal prism, in region 7.
    // - Element 8 is a nonagonal prism, in region 8.
    // - Element 9 is a decagonal prism, in region 9.
    // - Element 10 is a hendecagonal prism, in region 10.
    // All the elements belong to a region. Therefore, there is no "-1" region.
    // It contains 135 nodes, 203 edges, 90 faces.

    ASSERT_EQ( cellBlockManager.numNodes(), 135 );
    ASSERT_EQ( cellBlockManager.numEdges(), 203 );
    ASSERT_EQ( cellBlockManager.numFaces(), 90 );

    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), 135 );

    // 11 elements types x 11 regions = 121 sub-groups
    ASSERT_EQ( cellBlockManager.getCellBlocks().numSubGroups(), 121 );

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    CellBlockABC const & zone0 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "0_tetrahedra" );
    CellBlockABC const & zone1 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "1_pyramids" );
    CellBlockABC const & zone2 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "2_wedges" );
    CellBlockABC const & zone3 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "3_hexahedra" );
    CellBlockABC const & zone4 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "4_pentagonalPrisms" );
    CellBlockABC const & zone5 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "5_hexagonalPrisms" );
    CellBlockABC const & zone6 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "6_heptagonalPrisms" );
    CellBlockABC const & zone7 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "7_octagonalPrisms" );
    CellBlockABC const & zone8 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "8_nonagonalPrisms" );
    CellBlockABC const & zone9 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "9_decagonalPrisms" );
    CellBlockABC const & zone10 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "10_hendecagonalPrisms" );

    std::vector< string > const elementNames{ "tetrahedra",
                                              "pyramids",
                                              "wedges",
                                              "hexahedra",
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

    // Tetrahedron
    auto elementToNodes = zone0.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 4 );
    ASSERT_EQ( zone0.getElemToEdges().size( 1 ), 6 );
    ASSERT_EQ( zone0.getElemToFaces().size( 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 0 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 3 );

    // Pyramid
    elementToNodes = zone1.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 5 );
    ASSERT_EQ( zone1.getElemToEdges().size( 1 ), 8 );
    ASSERT_EQ( zone1.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 8 );

    // Wedges
    elementToNodes = zone2.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 6 );
    ASSERT_EQ( zone2.getElemToEdges().size( 1 ), 9 );
    ASSERT_EQ( zone2.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 11 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 14 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 13 );

    // Hexahedron
    elementToNodes = zone3.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 8 );
    ASSERT_EQ( zone3.getElemToEdges().size( 1 ), 12 );
    ASSERT_EQ( zone3.getElemToFaces().size( 1 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 15 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 16 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 18 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 17 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 19 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 20 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 22 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 21 );

    // Pentagonal prism
    elementToNodes = zone4.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 10 );
    ASSERT_EQ( zone4.getElemToEdges().size( 1 ), 15 );
    ASSERT_EQ( zone4.getElemToFaces().size( 1 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 23 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 24 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 25 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 26 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 27 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 28 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 29 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 30 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 31 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 32 );

    // Hexagonal prism
    elementToNodes = zone5.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 12 );
    ASSERT_EQ( zone5.getElemToEdges().size( 1 ), 18 );
    ASSERT_EQ( zone5.getElemToFaces().size( 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 33 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 34 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 35 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 36 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 37 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 38 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 39 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 40 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 41 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 42 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 43 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 44 );

    // Heptagonal prism
    elementToNodes = zone6.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 14 );
    ASSERT_EQ( zone6.getElemToEdges().size( 1 ), 21 );
    ASSERT_EQ( zone6.getElemToFaces().size( 1 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 45 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 46 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 47 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 48 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 49 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 50 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 51 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 52 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 53 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 54 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 55 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 56 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 57 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 58 );

    // Octagonal prism
    elementToNodes = zone7.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 16 );
    ASSERT_EQ( zone7.getElemToEdges().size( 1 ), 24 );
    ASSERT_EQ( zone7.getElemToFaces().size( 1 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 59 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 60 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 61 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 62 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 63 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 64 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 65 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 66 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 67 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 68 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 69 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 70 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 71 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 72 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 73 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 74 );

    // Nonagonal prism
    elementToNodes = zone8.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 18 );
    ASSERT_EQ( zone8.getElemToEdges().size( 1 ), 27 );
    ASSERT_EQ( zone8.getElemToFaces().size( 1 ), 11 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 75 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 76 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 77 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 78 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 79 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 80 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 81 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 82 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 83 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 84 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 85 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 86 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 87 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 88 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 89 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 90 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 91 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 92 );

    // Decagonal prism
    elementToNodes = zone9.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 20 );
    ASSERT_EQ( zone9.getElemToEdges().size( 1 ), 30 );
    ASSERT_EQ( zone9.getElemToFaces().size( 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 93 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 94 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 95 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 96 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 97 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 98 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 99 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 100 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 101 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 102 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 103 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 104 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 105 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 106 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 107 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 108 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 109 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 110 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 111 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 112 );

    // Hendecagonal prism
    elementToNodes = zone10.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 22 );
    ASSERT_EQ( zone10.getElemToEdges().size( 1 ), 33 );
    ASSERT_EQ( zone10.getElemToFaces().size( 1 ), 13 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 113 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 114 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 115 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 116 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 117 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 118 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 119 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 120 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 121 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 122 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 123 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 124 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 125 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 126 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 127 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 128 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 129 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 130 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 131 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 132 );
    EXPECT_EQ( elementToNodes( 0, 20 ), 133 );
    EXPECT_EQ( elementToNodes( 0, 21 ), 134 );

    for( auto const & z: { &zone0, &zone1, &zone2, &zone3, &zone4, &zone5, &zone6, &zone7, &zone8, &zone9, &zone10 } )
    {
      ASSERT_EQ( z->size(), 1 );
      ASSERT_EQ( z->getElemToNodes().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToEdges().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToFaces().size( 0 ), 1 );
    }
  };


  string const medleyVTK = testMeshDir + "/supportedElements.vtk";

  TestMeshImport( medleyVTK, validate );
}

TEST( VTKImport, medley42 )
{
  SKIP_TEST_IN_PARALLEL( "Neither relevant nor implemented in parallel" );

  auto validate = []( CellBlockManagerABC const & cellBlockManager ) -> void
  {
    // `supportedElements42.vtk` is the same as `supportedElements.vtk` with all eleven elements defined as polyhedron.
    // - Element 0 is a tetrahedron, in region 0.
    // - Element 1 is a pyramid, in region 1.
    // - Element 2 is a wedge, in region 2.
    // - Element 3 is a hexahedron, in region 3.
    // - Element 4 is a pentagonal prism, in region 4.
    // - Element 5 is a hexagonal prism, in region 5.
    // - Element 6 is a heptagonal prism, in region 6.
    // - Element 7 is a octagonal prism, in region 7.
    // - Element 8 is a nonagonal prism, in region 8.
    // - Element 9 is a decagonal prism, in region 9.
    // - Element 10 is a hendecagonal prism, in region 10.
    // All the elements belong to a region. Therefore, there is no "-1" region.
    // It contains 135 nodes, 203 edges, 90 faces.

    ASSERT_EQ( cellBlockManager.numNodes(), 135 );
    ASSERT_EQ( cellBlockManager.numEdges(), 203 );
    ASSERT_EQ( cellBlockManager.numFaces(), 90 );

    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), 135 );

    // 11 elements types x 11 regions = 121 sub-groups
    ASSERT_EQ( cellBlockManager.getCellBlocks().numSubGroups(), 121 );

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    CellBlockABC const & zone0 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "0_tetrahedra" );
    CellBlockABC const & zone1 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "1_pyramids" );
    CellBlockABC const & zone2 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "2_wedges" );
    CellBlockABC const & zone3 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "3_hexahedra" );
    CellBlockABC const & zone4 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "4_pentagonalPrisms" );
    CellBlockABC const & zone5 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "5_hexagonalPrisms" );
    CellBlockABC const & zone6 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "6_heptagonalPrisms" );
    CellBlockABC const & zone7 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "7_octagonalPrisms" );
    CellBlockABC const & zone8 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "8_nonagonalPrisms" );
    CellBlockABC const & zone9 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "9_decagonalPrisms" );
    CellBlockABC const & zone10 = cellBlockManager.getCellBlocks().getGroup< CellBlockABC >( "10_hendecagonalPrisms" );

    std::vector< string > const elementNames{ "tetrahedra",
                                              "pyramids",
                                              "wedges",
                                              "hexahedra",
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

    // Tetrahedron
    auto elementToNodes = zone0.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 4 );
    ASSERT_EQ( zone0.getElemToEdges().size( 1 ), 6 );
    ASSERT_EQ( zone0.getElemToFaces().size( 1 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 0 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 1 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 2 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 3 );

    // Pyramid
    elementToNodes = zone1.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 5 );
    ASSERT_EQ( zone1.getElemToEdges().size( 1 ), 8 );
    ASSERT_EQ( zone1.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 4 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 8 );

    // Wedges
    elementToNodes = zone2.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 6 );
    ASSERT_EQ( zone2.getElemToEdges().size( 1 ), 9 );
    ASSERT_EQ( zone2.getElemToFaces().size( 1 ), 5 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 11 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 14 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 13 );

    // Hexahedron
    elementToNodes = zone3.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 8 );
    ASSERT_EQ( zone3.getElemToEdges().size( 1 ), 12 );
    ASSERT_EQ( zone3.getElemToFaces().size( 1 ), 6 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 15 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 16 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 18 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 17 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 19 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 20 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 22 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 21 );

    // Pentagonal prism
    elementToNodes = zone4.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 10 );
    ASSERT_EQ( zone4.getElemToEdges().size( 1 ), 15 );
    ASSERT_EQ( zone4.getElemToFaces().size( 1 ), 7 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 23 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 24 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 25 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 26 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 27 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 28 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 29 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 30 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 31 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 32 );

    // Hexagonal prism
    elementToNodes = zone5.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 12 );
    ASSERT_EQ( zone5.getElemToEdges().size( 1 ), 18 );
    ASSERT_EQ( zone5.getElemToFaces().size( 1 ), 8 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 33 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 34 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 35 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 36 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 37 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 38 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 39 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 40 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 41 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 42 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 43 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 44 );

    // Heptagonal prism
    elementToNodes = zone6.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 14 );
    ASSERT_EQ( zone6.getElemToEdges().size( 1 ), 21 );
    ASSERT_EQ( zone6.getElemToFaces().size( 1 ), 9 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 45 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 46 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 47 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 48 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 49 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 50 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 51 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 52 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 53 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 54 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 55 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 56 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 57 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 58 );

    // Octagonal prism
    elementToNodes = zone7.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 16 );
    ASSERT_EQ( zone7.getElemToEdges().size( 1 ), 24 );
    ASSERT_EQ( zone7.getElemToFaces().size( 1 ), 10 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 59 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 60 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 61 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 62 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 63 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 64 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 65 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 66 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 67 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 68 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 69 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 70 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 71 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 72 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 73 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 74 );

    // Nonagonal prism
    elementToNodes = zone8.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 18 );
    ASSERT_EQ( zone8.getElemToEdges().size( 1 ), 27 );
    ASSERT_EQ( zone8.getElemToFaces().size( 1 ), 11 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 75 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 76 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 77 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 78 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 79 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 80 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 81 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 82 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 83 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 84 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 85 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 86 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 87 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 88 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 89 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 90 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 91 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 92 );

    // Decagonal prism
    elementToNodes = zone9.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 20 );
    ASSERT_EQ( zone9.getElemToEdges().size( 1 ), 30 );
    ASSERT_EQ( zone9.getElemToFaces().size( 1 ), 12 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 93 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 94 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 95 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 96 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 97 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 98 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 99 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 100 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 101 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 102 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 103 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 104 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 105 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 106 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 107 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 108 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 109 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 110 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 111 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 112 );

    // Hendecagonal prism
    elementToNodes = zone10.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 22 );
    ASSERT_EQ( zone10.getElemToEdges().size( 1 ), 33 );
    ASSERT_EQ( zone10.getElemToFaces().size( 1 ), 13 );
    EXPECT_EQ( elementToNodes( 0, 0 ), 113 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 114 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 115 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 116 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 117 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 118 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 119 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 120 );
    EXPECT_EQ( elementToNodes( 0, 8 ), 121 );
    EXPECT_EQ( elementToNodes( 0, 9 ), 122 );
    EXPECT_EQ( elementToNodes( 0, 10 ), 123 );
    EXPECT_EQ( elementToNodes( 0, 11 ), 124 );
    EXPECT_EQ( elementToNodes( 0, 12 ), 125 );
    EXPECT_EQ( elementToNodes( 0, 13 ), 126 );
    EXPECT_EQ( elementToNodes( 0, 14 ), 127 );
    EXPECT_EQ( elementToNodes( 0, 15 ), 128 );
    EXPECT_EQ( elementToNodes( 0, 16 ), 129 );
    EXPECT_EQ( elementToNodes( 0, 17 ), 130 );
    EXPECT_EQ( elementToNodes( 0, 18 ), 131 );
    EXPECT_EQ( elementToNodes( 0, 19 ), 132 );
    EXPECT_EQ( elementToNodes( 0, 20 ), 133 );
    EXPECT_EQ( elementToNodes( 0, 21 ), 134 );

    for( auto const & z: { &zone0, &zone1, &zone2, &zone3, &zone4, &zone5, &zone6, &zone7, &zone8, &zone9, &zone10 } )
    {
      ASSERT_EQ( z->size(), 1 );
      ASSERT_EQ( z->getElemToNodes().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToEdges().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToFaces().size( 0 ), 1 );
    }
  };

  string const medleyVTK = testMeshDir + "/supportedElements42.vtk";

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
