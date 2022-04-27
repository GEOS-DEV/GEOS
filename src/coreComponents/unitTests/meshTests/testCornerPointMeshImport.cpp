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
#include "dataRepository/xmlWrapper.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/generators/CellBlockManager.hpp"

// special CMake-generated include
#include "tests/meshDirName.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace geosx::dataRepository;

void testMeshImport( string const & inputStringMesh,
                     string const & inputStringRegion,
                     string const & meshBodyName,
                     string const & propertyToTest )
{
  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );

  // Load the mesh
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringMesh.c_str(), inputStringMesh.size() );

  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child( "Mesh" );
  meshManager.processInputFileRecursive( xmlMeshNode );
  meshManager.postProcessInputRecursive();

  // Create the domain and generate the Mesh
  DomainPartition domain( "domain", &root );
  meshManager.generateMeshes( domain );

  MeshBody & meshBody = domain.getMeshBody( meshBodyName );
  MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );
  NodeManager & nodeManager = meshLevel.getNodeManager();
  FaceManager const & faceManager = meshLevel.getFaceManager();
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  // Create the ElementRegions
  xmlDocument.load_buffer( inputStringRegion.c_str(), inputStringRegion.size() );

  xmlWrapper::xmlNode xmlRegionNode = xmlDocument.child( "ElementRegions" );
  elemManager.processInputFileRecursive( xmlRegionNode );
  elemManager.postProcessInputRecursive();

  CellBlockManager & cellBlockManager = meshBody.getGroup< CellBlockManager >( keys::cellManager );
  nodeManager.setGeometricalRelations( cellBlockManager, elemManager );

  elemManager.generateMesh( cellBlockManager );

  if( !propertyToTest.empty() )
  {
    elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( propertyToTest );
    } );
  }

  // Trigger import of the field
  meshManager.importFields( domain );

  // Check if the computed center match with the imported center
  if( !propertyToTest.empty() )
  {
    real64 const expectedSumPoreVolume = 1062756.27012;
    real64 sumPoreVolume = 0.0;
    real64 const eps = meshBody.getGlobalLengthScale() * 1e-8;
    elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
    {
      subRegion.calculateElementGeometricQuantities( nodeManager, faceManager );
      arrayView1d< real64 const > const elemVolume = subRegion.getElementVolume();
      arrayView1d< real64 const > const porosityProperty = subRegion.getReference< array1d< real64 > >( propertyToTest );

      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        sumPoreVolume += elemVolume[ei] * porosityProperty[ei];
      }
    } );
    EXPECT_LE( LvArray::math::abs( sumPoreVolume - expectedSumPoreVolume ), eps );
  }
}

TEST( CornerPointMeshImport, testECLIPSE )
{
  string const eclipseFilePath = testMeshDir + "/SPE10_small.GRDECL";

  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );

  std::stringstream inputStreamMesh;
  inputStreamMesh <<
    "<?xml version=\"1.0\" ?>" <<
    "  <Mesh xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
    " <CornerPointMesh" <<
    " name=\"mesh\"" <<
    " file=\"" << eclipseFilePath.c_str() << "\"" <<
    " fieldsToImport=\"{ PORO, PERM }\"" <<
    " fieldNamesInGEOSX=\"{ rockPorosity_referencePorosity, rockPerm_permeability }\"/>"
    "</Mesh>";
  const string inputStringMesh = inputStreamMesh.str();

  std::stringstream inputStreamRegion;
  inputStreamRegion <<
    "<?xml version=\"1.0\" ?>" <<
    "  <ElementRegions xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
    "  <CellElementRegion name=\"0\" cellBlocks=\"{ DEFAULT_HEX_0 }\" materialList=\"{ water, rock }\"/>" <<
    "</ElementRegions>";
  string inputStringRegion = inputStreamRegion.str();

  testMeshImport( inputStringMesh, inputStringRegion, "mesh", "rockPorosity_referencePorosity" );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
