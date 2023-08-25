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
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLMultiBlockDataWriter.h>

#include <gtest/gtest.h>
#include <conduit.hpp>

#include <filesystem>


using namespace geos;
using namespace geos::testing;
using namespace geos::dataRepository;


template< class V >
void TestMeshImport( string const & meshFilePath, V const & validate, string const fractureName="" )
{
  string const pattern = R"xml(
    <Mesh>
      <VTKMesh
        name="mesh"
        file="{}"
        partitionRefinement="0"
        useGlobalIds="0"
        {} />
    </Mesh>
  )xml";
  string const meshNode = GEOS_FMT( pattern, meshFilePath, fractureName.empty() ? "" : "faceBlocks=\"{" + fractureName + "}\"" );
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


class TestFractureImport : public ::testing::Test
{
protected:

  TestFractureImport() = default;

  inline static string const MULTI_BLOCK_NAME = "multi";

  std::filesystem::path m_vtkFile;

private:

  /// Folder where the vtk files will be written.
  std::filesystem::path m_vtkFolder;

  void SetUp() override
  {
    namespace fs = std::filesystem;

    fs::path const folder = fs::temp_directory_path();
    srand( (unsigned) time( nullptr ) );
    string const subFolder = "tmp-geos-vtk-" + std::to_string( rand() );
    m_vtkFolder = folder / subFolder;
    ASSERT_TRUE( fs::create_directory( m_vtkFolder ) );

    m_vtkFile = createFractureMesh( m_vtkFolder );
  }

  void TearDown() override
  {
    namespace fs = std::filesystem;

    // Carefully removing the files one by one and waiting the folders to be empty before removing them as well.
    // We do not want to remove important files!
    ASSERT_TRUE( fs::remove( m_vtkFolder / MULTI_BLOCK_NAME / ( MULTI_BLOCK_NAME + "_0.vtu" ) ) );
    ASSERT_TRUE( fs::remove( m_vtkFolder / MULTI_BLOCK_NAME / ( MULTI_BLOCK_NAME + "_1.vtu" ) ) );
    if( fs::is_empty( m_vtkFolder / MULTI_BLOCK_NAME ) )
    {
      ASSERT_TRUE( fs::remove( m_vtkFolder / MULTI_BLOCK_NAME ) );
    }
    ASSERT_TRUE( fs::remove( m_vtkFolder / ( MULTI_BLOCK_NAME + ".vtm" ) ) );
    if( fs::is_empty( m_vtkFolder ) )
    {
      ASSERT_TRUE( fs::remove( m_vtkFolder ) );
    }
  }

  static std::filesystem::path createFractureMesh( std::filesystem::path const & folder )
  {
    // The main mesh
    vtkNew< vtkUnstructuredGrid > main;
    {
      int constexpr numPoints = 16;
      double const pointsCoords[numPoints][3] = {
        { -1, 0, 0 },
        { -1, 1, 0 },
        { -1, 1, 1 },
        { -1, 0, 1 },
        { 0,  0, 0 },
        { 0,  1, 0 },
        { 0,  1, 1 },
        { 0,  0, 1 },
        { 0,  0, 0 },
        { 0,  1, 0 },
        { 0,  1, 1 },
        { 0,  0, 1 },
        { 1,  0, 0 },
        { 1,  1, 0 },
        { 1,  1, 1 },
        { 1,  0, 1 } };
      vtkNew< vtkPoints > points;
      points->Allocate( numPoints );
      for( double const * pointsCoord: pointsCoords )
      {
        points->InsertNextPoint( pointsCoord );
      }
      main->SetPoints( points );

      int constexpr numHexs = 2;
      vtkIdType const cubes[numHexs][8] = {
        { 0, 1, 2,  3,  4,  5,  6,  7 },
        { 8, 9, 10, 11, 12, 13, 14, 15 }
      };
      main->Allocate( numHexs );
      for( vtkIdType const * cube: cubes )
      {
        main->InsertNextCell( VTK_HEXAHEDRON, 8, cube );
      }

      vtkNew< vtkIdTypeArray > cellGlobalIds;
      cellGlobalIds->SetNumberOfComponents( 1 );
      cellGlobalIds->SetNumberOfTuples( numHexs );
      for( auto i = 0; i < numHexs; ++i )
      {
        cellGlobalIds->SetValue( i, i );
      }
      main->GetCellData()->SetGlobalIds( cellGlobalIds );

      vtkNew< vtkIdTypeArray > pointGlobalIds;
      pointGlobalIds->SetNumberOfComponents( 1 );
      pointGlobalIds->SetNumberOfTuples( numPoints );
      for( auto i = 0; i < numPoints; ++i )
      {
        pointGlobalIds->SetValue( i, i );
      }
      main->GetPointData()->SetGlobalIds( pointGlobalIds );
    }

    // The fracture mesh
    vtkNew< vtkUnstructuredGrid > fracture;
    {
      int constexpr numPoints = 4;
      double const pointsCoords[numPoints][3] = {
        { 0, 0, 0 },
        { 0, 1, 0 },
        { 0, 1, 1 },
        { 0, 0, 1 } };
      vtkNew< vtkPoints > points;
      points->Allocate( numPoints );
      for( double const * pointsCoord: pointsCoords )
      {
        points->InsertNextPoint( pointsCoord );
      }
      fracture->SetPoints( points );

      int constexpr numQuads = 1;
      vtkIdType const quad[numQuads][4] = { { 0, 1, 2, 3 } };
      fracture->Allocate( numQuads );
      for( vtkIdType const * q: quad )
      {
        fracture->InsertNextCell( VTK_QUAD, numPoints, q );
      }

      // Do not forget the collocated_nodes fields
      vtkNew< vtkIdTypeArray > collocatedNodes;
      collocatedNodes->SetName( "collocated_nodes" );
      collocatedNodes->SetNumberOfComponents( 2 );
      collocatedNodes->SetNumberOfTuples( numPoints );
      collocatedNodes->SetTuple2( 0, 4, 8 );
      collocatedNodes->SetTuple2( 1, 5, 9 );
      collocatedNodes->SetTuple2( 2, 6, 10 );
      collocatedNodes->SetTuple2( 3, 7, 11 );

      fracture->GetPointData()->AddArray( collocatedNodes );
    }

    vtkNew< vtkMultiBlockDataSet > multiBlock;
    multiBlock->SetNumberOfBlocks( 2 );
    multiBlock->SetBlock( 0, main );
    multiBlock->SetBlock( 1, fracture );
    multiBlock->GetMetaData( int( 0 ) )->Set( multiBlock->NAME(), "main" );
    multiBlock->GetMetaData( 1 )->Set( multiBlock->NAME(), "fracture" );

    vtkNew< vtkXMLMultiBlockDataWriter > writer;
    std::filesystem::path const vtkFile = folder / ( MULTI_BLOCK_NAME + ".vtm" );
    writer->SetFileName( vtkFile.c_str() );
    writer->SetInputData( multiBlock );
    writer->SetDataModeToAscii();
    writer->Write();

    return vtkFile;
  }
};


TEST_F( TestFractureImport, fracture )
{
  auto validate = []( CellBlockManagerABC const & cellBlockManager ) -> void
  {
    ASSERT_EQ( cellBlockManager.numNodes(), 16 );
    ASSERT_EQ( cellBlockManager.numEdges(), 24 );
    ASSERT_EQ( cellBlockManager.numFaces(), 12 );

    ASSERT_EQ( cellBlockManager.getFaceBlocks().numSubGroups(), 1 );
    FaceBlockABC const & faceBlock = cellBlockManager.getFaceBlocks().getGroup< FaceBlockABC >( 0 );
    ASSERT_EQ( faceBlock.num2dElements(), 1 );
    ASSERT_EQ( faceBlock.num2dFaces(), 4 );
    auto ecn = faceBlock.get2dElemsToCollocatedNodesBuckets();
    ASSERT_EQ( ecn[0].size(), 4 );
    for( int i = 0; i < 4; ++i )
    {
      auto bucket = ecn( 0, i );
      ASSERT_EQ( bucket.size(), 2 );
      std::set< globalIndex > result( bucket.begin(), bucket.end() );
      ASSERT_EQ( result, std::set< globalIndex >( { 4 + i, 8 + i } ) );
    }
  };

  TestMeshImport( m_vtkFile, validate, "fracture" );
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

TEST( VTKImport, supportedElements )
{
  SKIP_TEST_IN_PARALLEL( "Neither relevant nor implemented in parallel" );

  auto validate = []( CellBlockManagerABC const & cellBlockManager ) -> void
  {
    // `supportedElements.vtk` is made of twelve elements.
    // `supportedElementsAsVTKPolyhedra.vtk` is the same model with all elements defined as polyhedron.
    // - Element 0 is a tetrahedron, in region 0.
    // - Element 1 is a pyramid, in region 1.
    // - Element 2 is a wedge, in region 2.
    // - Element 3 is an hexahedron (converted from VTK_VOXEL), in region 3.
    // - Element 4 is an hexahedron, in region 3.
    // - Element 5 is a pentagonal prism, in region 4.
    // - Element 6 is an hexagonal prism, in region 5.
    // - Element 7 is an heptagonal prism, in region 6.
    // - Element 8 is an octagonal prism, in region 7.
    // - Element 9 is a nonagonal prism, in region 8.
    // - Element 10 is a decagonal prism, in region 9.
    // - Element 11 is an hendecagonal prism, in region 10.
    // All the elements belong to a region. Therefore, there is no "-1" region.
    // It contains 143 nodes, 215 edges, 96 faces.

    ASSERT_EQ( cellBlockManager.numNodes(), 143 );
    ASSERT_EQ( cellBlockManager.numEdges(), 215 );
    ASSERT_EQ( cellBlockManager.numFaces(), 96 );

    SortedArray< localIndex > const & allNodes = cellBlockManager.getNodeSets().at( "all" );
    ASSERT_EQ( allNodes.size(), 143 );

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
        // zone3 contains two hexahedra converted from VTK_VOXEL and VTK_HEXAHEDRON or corresponding VTK_POLYHEDRON
        if( prefix == 3 )
        {
          ASSERT_EQ( zone.size(), prefix == i ? 2 : 0 );
        }
        else
        {
          ASSERT_EQ( zone.size(), prefix == i ? 1 : 0 );
        }
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

    // Hexahedron (from VTK_VOXEL and VTK_HEXAHEDRON or corresponding VTK_POLYHEDRON)
    elementToNodes = zone3.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 8 );
    ASSERT_EQ( zone3.getElemToEdges().size( 1 ), 12 );
    ASSERT_EQ( zone3.getElemToFaces().size( 1 ), 6 );

    EXPECT_EQ( elementToNodes( 0, 0 ), 15 );
    EXPECT_EQ( elementToNodes( 0, 1 ), 16 );
    EXPECT_EQ( elementToNodes( 0, 2 ), 17 );
    EXPECT_EQ( elementToNodes( 0, 3 ), 18 );
    EXPECT_EQ( elementToNodes( 0, 4 ), 19 );
    EXPECT_EQ( elementToNodes( 0, 5 ), 20 );
    EXPECT_EQ( elementToNodes( 0, 6 ), 21 );
    EXPECT_EQ( elementToNodes( 0, 7 ), 22 );

    EXPECT_EQ( elementToNodes( 1, 0 ), 23 );
    EXPECT_EQ( elementToNodes( 1, 1 ), 24 );
    EXPECT_EQ( elementToNodes( 1, 2 ), 26 );
    EXPECT_EQ( elementToNodes( 1, 3 ), 25 );
    EXPECT_EQ( elementToNodes( 1, 4 ), 27 );
    EXPECT_EQ( elementToNodes( 1, 5 ), 28 );
    EXPECT_EQ( elementToNodes( 1, 6 ), 30 );
    EXPECT_EQ( elementToNodes( 1, 7 ), 29 );

    // Pentagonal prism
    elementToNodes = zone4.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 10 );
    ASSERT_EQ( zone4.getElemToEdges().size( 1 ), 15 );
    ASSERT_EQ( zone4.getElemToFaces().size( 1 ), 7 );
    for( int i=0; i < elementToNodes.size( 1 ); ++i )
    {
      EXPECT_EQ( elementToNodes( 0, i ), 31+i );
    }

    // Hexagonal prism
    elementToNodes = zone5.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 12 );
    ASSERT_EQ( zone5.getElemToEdges().size( 1 ), 18 );
    ASSERT_EQ( zone5.getElemToFaces().size( 1 ), 8 );
    for( int i=0; i < elementToNodes.size( 1 ); ++i )
    {
      EXPECT_EQ( elementToNodes( 0, i ), 41+i );
    }

    // Heptagonal prism
    elementToNodes = zone6.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 14 );
    ASSERT_EQ( zone6.getElemToEdges().size( 1 ), 21 );
    ASSERT_EQ( zone6.getElemToFaces().size( 1 ), 9 );
    for( int i=0; i < elementToNodes.size( 1 ); ++i )
    {
      EXPECT_EQ( elementToNodes( 0, i ), 53+i );
    }

    // Octagonal prism
    elementToNodes = zone7.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 16 );
    ASSERT_EQ( zone7.getElemToEdges().size( 1 ), 24 );
    ASSERT_EQ( zone7.getElemToFaces().size( 1 ), 10 );
    for( int i=0; i < elementToNodes.size( 1 ); ++i )
    {
      EXPECT_EQ( elementToNodes( 0, i ), 67+i );
    }

    // Nonagonal prism
    elementToNodes = zone8.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 18 );
    ASSERT_EQ( zone8.getElemToEdges().size( 1 ), 27 );
    ASSERT_EQ( zone8.getElemToFaces().size( 1 ), 11 );
    for( int i=0; i < elementToNodes.size( 1 ); ++i )
    {
      EXPECT_EQ( elementToNodes( 0, i ), 83+i );
    }

    // Decagonal prism
    elementToNodes = zone9.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 20 );
    ASSERT_EQ( zone9.getElemToEdges().size( 1 ), 30 );
    ASSERT_EQ( zone9.getElemToFaces().size( 1 ), 12 );
    for( int i=0; i < elementToNodes.size( 1 ); ++i )
    {
      EXPECT_EQ( elementToNodes( 0, i ), 101+i );
    }

    // Hendecagonal prism
    elementToNodes = zone10.getElemToNodes();
    ASSERT_EQ( elementToNodes.size( 1 ), 22 );
    ASSERT_EQ( zone10.getElemToEdges().size( 1 ), 33 );
    ASSERT_EQ( zone10.getElemToFaces().size( 1 ), 13 );
    for( int i=0; i < elementToNodes.size( 1 ); ++i )
    {
      EXPECT_EQ( elementToNodes( 0, i ), 121+i );
    }

    for( auto const & z: { &zone0, &zone1, &zone2, &zone4, &zone5, &zone6, &zone7, &zone8, &zone9, &zone10 } )
    {
      ASSERT_EQ( z->size(), 1 );
      ASSERT_EQ( z->getElemToNodes().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToEdges().size( 0 ), 1 );
      ASSERT_EQ( z->getElemToFaces().size( 0 ), 1 );
    }

    ASSERT_EQ( zone3.size(), 2 );
    ASSERT_EQ( zone3.getElemToNodes().size( 0 ), 2 );
    ASSERT_EQ( zone3.getElemToEdges().size( 0 ), 2 );
    ASSERT_EQ( zone3.getElemToFaces().size( 0 ), 2 );
  };

  string const medleyVTK = testMeshDir + "/supportedElements.vtk";
  TestMeshImport( medleyVTK, validate );

  string const medleyVTK42 = testMeshDir + "/supportedElementsAsVTKPolyhedra.vtk";
  TestMeshImport( medleyVTK42, validate );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
