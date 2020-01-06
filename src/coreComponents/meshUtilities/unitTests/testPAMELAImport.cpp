/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


// Source includes
#include "managers/initialization.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "tests/meshFileNames.hpp"
#include "meshUtilities/MeshManager.hpp"

// TPL includes
#include <gtest/gtest.h>


using namespace geosx;
using namespace geosx::dataRepository;

TEST( PAMELAImport, testXML )
{
  MeshManager meshManager("mesh", nullptr);

  // Load the mesh
  std::stringstream inputStreamMesh;
  inputStreamMesh <<
  "<Mesh>" <<
  "  <PAMELAMeshGenerator name=\"ToyModel\" " <<
  "  fieldsToImport=\"{barycenter}\""<<
  "  fieldNamesInGEOSX=\"{barycenter}\""<<
  "  file=\"" <<gmshFilePath.c_str()<< "\"/>"<<
  "</Mesh>";
  const string inputStringMesh = inputStreamMesh.str();
  std::cout << gmshFilePath.c_str() << std::endl;

  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringMesh.c_str(), inputStringMesh.size() );

  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child("Mesh");
  meshManager.ProcessInputFileRecursive( xmlMeshNode );
  meshManager.PostProcessInputRecursive();

  // Create the domain and generate the Mesh
  auto domain = std::unique_ptr< DomainPartition >( new DomainPartition( "domain", nullptr ) );
  meshManager.GenerateMeshes( domain.get() );

  Group * const meshBodies = domain->getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // Create the ElementRegions
  std::stringstream inputStreamRegion;
  inputStreamRegion <<
  "<ElementRegions>" <<
  "  <CellElementRegion name=\"0\" cellBlocks=\"{0_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <CellElementRegion name=\"1\" cellBlocks=\"{1_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <CellElementRegion name=\"2\" cellBlocks=\"{2_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <CellElementRegion name=\"3\" cellBlocks=\"{3_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <CellElementRegion name=\"4\" cellBlocks=\"{4_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <CellElementRegion name=\"5\" cellBlocks=\"{5_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <CellElementRegion name=\"6\" cellBlocks=\"{6_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <CellElementRegion name=\"7\" cellBlocks=\"{7_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "</ElementRegions>";
  const string inputStringRegion = inputStreamRegion.str();

  xmlDocument.load_buffer( inputStringRegion.c_str(), inputStringRegion.size() );

  xmlWrapper::xmlNode xmlRegionNode = xmlDocument.child("ElementRegions");
  elemManager->ProcessInputFileRecursive( xmlRegionNode );
  elemManager->PostProcessInputRecursive();

  Group const * const cellBlockManager = domain->GetGroup(keys::cellManager);

  // This method will call the CopyElementSubRegionFromCellBlocks that will trigger the property transfer.
  elemManager->GenerateMesh( cellBlockManager );


  // Check if the computed center match with the importer center
  auto centerProperty =  elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor> >( "barycenter" );
  elemManager->forElementRegions( [&]( ElementRegionBase * const elemRegion)->void
  {
    localIndex er = elemRegion->getIndexInParent();
    elemRegion->forElementSubRegionsIndex( [&]( localIndex const esr, auto * const elemSubRegion )
    {
      for( localIndex ei = 0; ei < elemSubRegion->size(); ei++ )
      {
        R1Tensor center = elemSubRegion->calculateElementCenter( ei, *nodeManager );
        R1Tensor centerFromProperty( centerProperty[er][esr][ei][0], 
                                     centerProperty[er][esr][ei][1],
                                     centerProperty[er][esr][ei][2] );
        center -= centerFromProperty;
        GEOSX_ERROR_IF( center.L2_Norm() > meshBody->getGlobalLengthScale() * 1e-8, "Property import of centers if wrong");
      }
    });
  });

}

int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
