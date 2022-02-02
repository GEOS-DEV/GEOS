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
//#include "mesh/generators/VTKMeshGenerator.hpp"
#include "tests/meshDirName.hpp"
#include "mesh/MeshManager.hpp"
//#include "mesh/mpiCommunications/CommunicationTools.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace geosx::dataRepository;

// TODO move in the functions?
static const std::string vtkFilePath = testMeshDir + "/cube.vtk";
static const std::string vtuFilePath = testMeshDir +  + "/cube.vtu";
static const std::string pvtuFilePath = testMeshDir +  + "/cube.pvtu";
static const std::string medleyFilePath = testMeshDir +  + "/medley.vtk";

void TestMeshImport( string const & meshFilePath )
{
  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );

  std::stringstream inputStreamMesh;
  inputStreamMesh <<
    "<?xml version=\"1.0\" ?>" <<
    "  <Mesh xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
    "  <VTKMeshGenerator name=\"Cube\" " <<
    "  file=\"" << meshFilePath.c_str()<< "\"/>"<<
    "</Mesh>";
  const string inputStringMesh = inputStreamMesh.str();

  std::stringstream inputStreamRegion;
  inputStreamRegion <<
    "<?xml version=\"1.0\" ?>" <<
    "  <CellElementRegions xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
    "  <CellElementRegion name=\"0\" cellBlocks=\"{hexahedron}\" materialList=\"{water, rock}\"/>" <<
    "</ElementRegions>";
  string inputStringRegion = inputStreamRegion.str();

  // Load the mesh
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringMesh.c_str(), inputStringMesh.size() );

  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child( "Mesh" );
  meshManager.processInputFileRecursive( xmlMeshNode );
  meshManager.postProcessInputRecursive();

  // Create the domain and generate the Mesh
  auto domain = std::unique_ptr< DomainPartition >( new DomainPartition( "domain", &root ) );
  meshManager.generateMeshes( *domain );

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

}

TEST( VTKImport, testVTK )
{
  TestMeshImport( vtkFilePath );
}

TEST( VTKImport, testVTU )
{
  TestMeshImport( vtuFilePath );
}

TEST( VTKImport, testPVTU )
{
  TestMeshImport( pvtuFilePath );
}

template< class V >
void MyTestMeshImport( string const & meshFilePath, V const & validate )
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
//  auto domain = std::unique_ptr< DomainPartition >( new DomainPartition( "domain", &root ) );
  DomainPartition * domain = new DomainPartition( "domain", &root );
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

  std::cout << "Hello" << std::endl;
}

TEST( MyVTKImport, serialCube )
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

//    NodeManager const & nodeManager = meshLevel.getNodeManager();
    localIndex const numNodes = nodeManager.size();
    ASSERT_EQ( numNodes, 64 );
    ASSERT_EQ( nodeManager.faceList().size(), 64 ); // The information are not filled yet.
    ASSERT_EQ( nodeManager.elementList().size(), 64 ); // The information are not filled yet.

    // We have all the 4 x 4  x 4 = 64 nodes in the "all" set.
    WrapperBase const & allNodes = nodeManager.sets().getWrapperBase( "all" );
    ASSERT_EQ( allNodes.size(), 64 );
    // The "-1" set contains 3 quads connected in L shape.
    WrapperBase const & noRegionNodes = nodeManager.sets().getWrapperBase( "-1" );
    ASSERT_EQ( noRegionNodes.size(), 8 );
    // The "2" set are all the boundary nodes (64 - 8 inside nodes = 56),
    // minus an extra node (node 55) that belongs to regions 2 and 9 only.
    WrapperBase const & nodesRegion2 = nodeManager.sets().getWrapperBase( "2" );
    ASSERT_EQ( nodesRegion2.size(), 55 );
    // Region "9" has only one quad.
    WrapperBase const & nodesRegion9 = nodeManager.sets().getWrapperBase( "9" );
    ASSERT_EQ( nodesRegion9.size(), 4 );

    // 1 elements type on 3 regions ("-1", "3", "9") = 3 sub-groups
    ASSERT_EQ( cellBlockManager.getGroup( keys::cellBlocks ).numSubGroups(), 3 );

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

    std::pair< CellBlock const *, localIndex > const p{ &hexs, 1 }, p3{ &hexs3, 25 }, p9{ &hexs9, 1 };
    for( auto const & hs: { p, p3, p9 } )
    {
      CellBlock const * h = hs.first;
      localIndex const & result = hs.second;

      ASSERT_EQ( h->size(), result );
      ASSERT_EQ( h->nodeList().size( 0 ), result );
      ASSERT_EQ( h->edgeList().size( 0 ), result );
      ASSERT_EQ( h->faceList().size( 0 ), result );
    }

    // TODO importFields?
  };

  MyTestMeshImport( vtkFilePath, validate );
  MyTestMeshImport( vtuFilePath, validate );
  MyTestMeshImport( vtkFilePath, validate );
}

TEST( MyVTKImport, medley )
{
//  CommunicationTools c;

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

  MyTestMeshImport( medleyFilePath, validate );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
