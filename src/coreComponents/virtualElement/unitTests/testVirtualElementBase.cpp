/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "managers/ProblemManager.hpp"
#include "virtualElement/VirtualElementBase.hpp"
#include "virtualElement/ConformingVirtualElement_1.hpp"
#include "managers/DomainPartition.hpp"
#include "meshUtilities/MeshManager.hpp"

// TPL includes
#include "gtest/gtest.h"

using namespace geosx;

TEST( VirtualElementBase, compilation )
{
  string const inputStream=
    "<Problem>"
    "  <Mesh>"
    "    <InternalMesh"
    "      name=\"cube\""
    "      elementTypes=\"{C3D8}\""
    "      xCoords=\"{0.0, 1.0}\""
    "      yCoords=\"{0.0, 2.0}\""
    "      zCoords=\"{0.0, 3.0}\""
    "      nx=\"{10}\""
    "      ny=\"{11}\""
    "      nz=\"{12}\""
    "      cellBlockNames=\"{cb1}\""
    "    />"
    "  </Mesh>"
    "  <ElementRegions>"
    "    <CellElementRegion name=\"region1\" cellBlocks=\"{cb1}\""
    "                       materialList=\"{dummy_material}\" />"
    "  </ElementRegions>"
    "</Problem>";
  xmlWrapper::xmlDocument inputFile;
  xmlWrapper::xmlResult xmlResult = inputFile.load_buffer(inputStream.c_str(), inputStream.size());
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }
  xmlWrapper::xmlNode xmlProblemNode = inputFile.child( "Problem" );

  ProblemManager * problemManager = new ProblemManager( "Problem", nullptr );
  problemManager->InitializePythonInterpreter();
  problemManager->ProcessInputFileRecursive( xmlProblemNode );

  // Open mesh levels
  DomainPartition * domain  = problemManager->getDomainPartition();
  MeshManager * meshManager = problemManager->GetGroup< MeshManager >( problemManager->groupKeys.meshManager );
  meshManager->GenerateMeshLevels( domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * elementManager = mesh.getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager->getName().c_str() );
  elementManager->ProcessInputFileRecursive( topLevelNode );
  elementManager->PostProcessInputRecursive();
  std::cout << elementManager->getNumberOfElements() << std::endl << std::endl;
  std::cout << mesh.getNodeManager()->referencePosition() << std::endl << std::endl;
  problemManager->ProblemSetup();

  delete problemManager;
}
