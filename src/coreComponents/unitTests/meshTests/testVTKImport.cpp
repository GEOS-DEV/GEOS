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
#include "mainInterface/initialization.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "tests/meshDirName.hpp"
#include "mesh/MeshManager.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

// TODO unduplicate
#define SKIP_TEST_IF( COND, REASON ) \
  do \
  { \
    if( COND ) \
    { \
      GTEST_SKIP_( ": " REASON ); \
    } \
  } while( 0 )

#define SKIP_TEST_IN_SERIAL( REASON ) \
  do \
  { \
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize == 1, REASON ); \
  } while( 0 )

#define SKIP_TEST_IN_PARALLEL( REASON ) \
  do \
  { \
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize != 1, REASON ); \
  } while( 0 )


using namespace geosx;
using namespace geosx::dataRepository;

// TODO move in the functions?
static const std::string cubeVTK = testMeshDir + "/cube.vtk";
static const std::string cubeVTU = testMeshDir + + "/cube.vtu";
static const std::string cubePVTU = testMeshDir + + "/cube.pvtu";
static const std::string medleyVTK = testMeshDir + + "/medley.vtk";

template< class V >
void TestMeshImport( string const & meshFilePath, V const & validate )
{
  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );

  std::stringstream inputStreamMesh;
  inputStreamMesh <<
                  "<?xml version=\"1.0\" ?>" <<
                  "<Mesh xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
                  "  <VTKMeshGenerator name=\"Cube\" file=\"" << meshFilePath << "\" />"<<
                  "</Mesh>";
  const string inputStringMesh = inputStreamMesh.str();

  std::stringstream inputStreamRegion;
  inputStreamRegion <<
                    "<?xml version=\"1.0\" ?>" <<
                    "  <CellElementRegions xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
                    "  <CellElementRegion name=\"0\" cellBlocks=\"{hexahedron, 0_pyramids}\" materialList=\"{water, rock}\"/>" <<
                    "</ElementRegions>"; // FIXME What is this cellBlocks="{hexahedron}"?
  string inputStringRegion = inputStreamRegion.str();

  // Load the mesh
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringMesh.c_str(), inputStringMesh.size() );

  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child( "Mesh" );
  meshManager.processInputFileRecursive( xmlMeshNode );
  meshManager.postProcessInputRecursive();

  // Create the domain and generate the Mesh
  auto domain = std::make_unique< DomainPartition >( "domain", &root );
  meshManager.generateMeshes( *domain );
//  meshManager.importFields( *domain );

  MeshBody & meshBody = domain->getMeshBody( 0 );
  MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  // Create the ElementRegions
  xmlDocument.load_buffer( inputStringRegion.c_str(), inputStringRegion.size() );

  xmlWrapper::xmlNode xmlRegionNode = xmlDocument.child( "ElementRegions" );
  elemManager.processInputFileRecursive( xmlRegionNode );
  elemManager.postProcessInputRecursive();

  Group & cellBlockManager = domain->getGroup( keys::cellManager );

  elemManager.generateMesh( cellBlockManager );

  validate( meshLevel.getNodeManager(), cellBlockManager );
}

/**
 * @brief Returns the expected value depending on the MPI context.
 * @tparam T Type of the expected value.
 * @param[in] expectedSerial Expected value for serial case.
 * @param[in] expectedParallel Expected values for parallel MPI cases, for the current MPI rank.
 *        The length of the list should match the size of the MPI communicator @p comm.
 *        The @p i^th element of the list will be chosen for MPI rank @p i.
 * @param[in] comm The MPI_Comm communicator that the function will act on.
 * @return The expected value.
 *
 * @note This function is meant to be used to run the same test in serial or parallel environments.
 */
template< typename T >
T expected( T expectedSerial,
            std::initializer_list< T > expectedParallel,
            MPI_Comm const & comm=MPI_COMM_GEOSX )
{
  int const mpiSize = MpiWrapper::commSize( comm );
  if( mpiSize == 1 )
  {
    return expectedSerial;
  }
  else
  {
    GEOSX_ASSERT( expectedParallel.size() == std::size_t( mpiSize ) );
    std::vector< T > tmp( expectedParallel );
    return tmp[MpiWrapper::commRank( comm )];
  }
}


TEST( MyVTKImport, cube )
{
  auto validate = [](NodeManager const & nodeManager, Group const & cellBlockManager) -> void
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
    localIndex const expectedNumNodes = expected(64 , { 48, 32 });
    ASSERT_EQ( nodeManager.size(), expectedNumNodes );
    // The information in the tables is not filled yet. We can check the consistency of the sizes.
    ASSERT_EQ( nodeManager.faceList().size(), expectedNumNodes );
    ASSERT_EQ( nodeManager.elementList().size(), expectedNumNodes );

    // We have all the 4 x 4  x 4 = 64 nodes in the "all" set.
    WrapperBase const & allNodes = nodeManager.sets().getWrapperBase( "all" );
    ASSERT_EQ( allNodes.size(), expectedNumNodes );

    // The "-1" set contains 3 quads connected in L shape.
    // 2 of those quads touch in the center of the cube: 6 nodes for MPI rank 0.
    // Rank 1 only has one cube: 4 nodes.
    WrapperBase const & nodesNoRegion = nodeManager.sets().getWrapperBase( "-1" );
    ASSERT_EQ( nodesNoRegion.size(), expected( 8, { 6, 4 } ) );

    // The "2" set are all the boundary nodes (64 - 8 inside nodes = 56),
    // minus an extra node that belongs to regions -1 and 9 only.
    WrapperBase const & nodesRegion2 = nodeManager.sets().getWrapperBase( "2" );
    ASSERT_EQ( nodesRegion2.size(), expected( 55, { 39, 27 } ) );

    // Region "9" has only one quad, on the greater `x` direction.
    // This hex will belong to MPI rank 1.
    WrapperBase const & nodesRegion9 = nodeManager.sets().getWrapperBase( "9" );
    ASSERT_EQ( nodesRegion9.size(), expected( 4, { 0, 4 } ) );

    // 1 elements type on 3 regions ("-1", "3", "9") = 3 sub-groups
    ASSERT_EQ( cellBlockManager.getGroup( keys::cellBlocks ).numSubGroups(), 3 );
    GEOSX_LOG_RANK(cellBlockManager.getGroup( keys::cellBlocks ).getGroup(0).getName());
    GEOSX_LOG_RANK(cellBlockManager.getGroup( keys::cellBlocks ).getGroup(1).getName());
    GEOSX_LOG_RANK(cellBlockManager.getGroup( keys::cellBlocks ).getGroup(2).getName());

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    CellBlock const & hexs = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( "hexahedron" );
    CellBlock const & hexs3 = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( "3_hexahedron" );
    CellBlock const & hexs9 = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( "9_hexahedron" );

    for( CellBlock const * h: { &hexs, &hexs3, &hexs9 } )
    {
      // 8 nodes, 12 edges and 6 faces per hex.
      ASSERT_EQ( h->nodeList().size( 1 ), 8 );
      ASSERT_EQ( h->edgeList().size( 1 ), 12 );
      ASSERT_EQ( h->faceList().size( 1 ), 6 );
    }

    std::pair< CellBlock const *, localIndex > const p{ &hexs, expected( 1, { 1, 0 } ) },
      p3{ &hexs3, expected( 25, { 17, 8 } ) },
      p9{ &hexs9, expected( 1, { 0, 1 } ) };
    for( auto const & hs: { p, p3, p9 } )
    {
      CellBlock const * h = hs.first;
      localIndex const & expectedSize = hs.second;

      ASSERT_EQ( h->size(), expectedSize );
      ASSERT_EQ( h->nodeList().size( 0 ), expectedSize );
      ASSERT_EQ( h->edgeList().size( 0 ), expectedSize );
      ASSERT_EQ( h->faceList().size( 0 ), expectedSize );
    }

    // TODO importFields?
  };

  TestMeshImport( cubeVTK, validate );
  TestMeshImport( cubeVTU, validate );
//  TestMeshImport( cubePVTU, validate );
}

TEST( MyVTKImport, medley )
{
  SKIP_TEST_IN_PARALLEL("Neither relevant nor implemented in parallel");

  auto validate = [](NodeManager const & nodeManager, Group const & cellBlockManager) -> void
  {
    // `medley.vtk` is made of four elements.
    // - Element 0 is a pyramid, in region 0.
    // - Element 1 is an hexahedron, in region 1.
    // - Element 2 is a wedge, in region 2.
    // - Element 3 is a tetrahedron, in region 3.
    // All the elements belong to a region. Therefore, there is no "-1" region.
    // It contains 12 nodes.

   localIndex const numNodes = nodeManager.size();
    ASSERT_EQ( numNodes, 12 );

    WrapperBase const & allNodes = nodeManager.sets().getWrapperBase( "all" );
    ASSERT_EQ( allNodes.size(), 12 );

    // 4 elements types x 4 regions = 16 sub-groups
    ASSERT_EQ( cellBlockManager.getGroup( keys::cellBlocks ).numSubGroups(), 16 );

    // FIXME How to get the CellBlock as a function of the region, without knowing the naming pattern.
    CellBlock const & zone0 = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( "0_pyramids" );
    CellBlock const & zone1 = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( "1_hexahedron" );
    CellBlock const & zone2 = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( "2_wedges" );
    CellBlock const & zone3 = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( "3_tetrahedron" );

    for( string prefix: { "0", "1", "2", "3" } )
    {
      for( string elementType: { "pyramids", "hexahedron", "wedges", "tetrahedron" } )
      {
        string const name = prefix + "_" + elementType;
        CellBlock const & zone = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( name );
        ASSERT_EQ( zone.size(), name == "0_pyramids" or name == "1_hexahedron" or name == "2_wedges" or name == "3_tetrahedron" ? 1 : 0 );
      }
    }

    // Pyramid
    ASSERT_EQ( zone0.nodeList().size( 1 ), 5 );
    ASSERT_EQ( zone0.edgeList().size( 1 ), 8 );
    ASSERT_EQ( zone0.faceList().size( 1 ), 5 );

    // Hexahedron
    ASSERT_EQ( zone1.nodeList().size( 1 ), 8 );
    ASSERT_EQ( zone1.edgeList().size( 1 ), 12 );
    ASSERT_EQ( zone1.faceList().size( 1 ), 6 );

    // Wedges
    ASSERT_EQ( zone2.nodeList().size( 1 ), 6 );
    ASSERT_EQ( zone2.edgeList().size( 1 ), 9 );
    ASSERT_EQ( zone2.faceList().size( 1 ), 5 );

    // Tetrahedron
    ASSERT_EQ( zone3.nodeList().size( 1 ), 4 );
    ASSERT_EQ( zone3.edgeList().size( 1 ), 6 );
    ASSERT_EQ( zone3.faceList().size( 1 ), 4 );

    for( auto const & z: { &zone0, &zone1, &zone2, &zone3 } )
    {
      ASSERT_EQ( z->size(), 1 );
      ASSERT_EQ( z->nodeList().size( 0 ), 1 );
      ASSERT_EQ( z->edgeList().size( 0 ), 1 );
      ASSERT_EQ( z->faceList().size( 0 ), 1 );
    }
  };

  TestMeshImport( medleyVTK, validate );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
