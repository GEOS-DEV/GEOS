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
#include "tests/meshFileNames.hpp"
#include "mesh/MeshManager.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace geosx::dataRepository;

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

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
