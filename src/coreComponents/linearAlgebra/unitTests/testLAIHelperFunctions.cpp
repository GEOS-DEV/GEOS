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

/*
 * testLAIHelperFunctions.cpp
 *  Created on: Oct 29, 2019
 */

#include "gtest/gtest.h"

#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

using namespace geosx;

namespace
{
int global_argc;
char** global_argv;
}

class LAIHelperFunctionsTest : public ::testing::Test
{
public:
  /**
   * @brief Set the timer function.
   *
   * @param [out] time double actual time.
   */
  void setTimer( double &time )
  {
    time = MpiWrapper::Wtime();
  }

  /**
   * @brief Get elapsed time.
   *
   * @param [inout] time double
   * - on input: actual time from setTimer;
   * - on output: elapsed time.
   */
  void getElapsedTime( double &time )
  {
    time = MpiWrapper::Wtime() - time;
  }

protected:

  /**
   * @brief Provide a very simple xml input, with an internally generated mesh.
   */
  static void SetUpTestCase()
  {
    problemManager = new ProblemManager("Problem", nullptr);

    string const inputStream =
    "<Problem>"
    "  <Mesh>"
    "    <InternalMesh name=\"mesh1\""
    "                  elementTypes=\"{C3D8}\""
    "                  xCoords=\"{0, 1, 2, 3, 4}\""
    "                  yCoords=\"{0, 1}\""
    "                  zCoords=\"{0, 1}\""
    "                  nx=\"{4, 4, 4, 4}\""
    "                  ny=\"{4}\""
    "                  nz=\"{5}\""
    "                  cellBlockNames=\"{block1, block2, block3, block4}\"/>"
    "  </Mesh>"
    "  <ElementRegions>"
    "    <CellElementRegion name=\"region1\" cellBlocks=\"{block1}\" materialList=\"{dummy_material}\" />"
    "    <CellElementRegion name=\"region2\" cellBlocks=\"{block2}\" materialList=\"{dummy_material}\" />"
    "    <CellElementRegion name=\"region3\" cellBlocks=\"{block3}\" materialList=\"{dummy_material}\" />"
    "    <CellElementRegion name=\"region4\" cellBlocks=\"{block4}\" materialList=\"{dummy_material}\" />"
    "  </ElementRegions>"
    "</Problem>";

    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
    if (!xmlResult)
    {
      GEOS_LOG_RANK_0("XML parsed with errors!");
      GEOS_LOG_RANK_0("Error description: " << xmlResult.description());
      GEOS_LOG_RANK_0("Error offset: " << xmlResult.offset);
    }

    int mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
    dataRepository::Group * commandLine =
      problemManager->GetGroup<dataRepository::Group>( problemManager->groupKeys.commandLine );
    commandLine->registerWrapper<integer>( problemManager->viewKeys.xPartitionsOverride.Key() )->
      setApplyDefaultValue(mpiSize);

    xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( "Problem" );
    problemManager->InitializePythonInterpreter();
    problemManager->ProcessInputFileRecursive( xmlProblemNode );

    // Open mesh levels
    DomainPartition * domain  = problemManager->getDomainPartition();
    MeshManager * meshManager = problemManager->GetGroup<MeshManager>( problemManager->groupKeys.meshManager );
    meshManager->GenerateMeshLevels(domain);

    ElementRegionManager * elementManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager->getName().c_str() );
    elementManager->ProcessInputFileRecursive( topLevelNode );
    elementManager->PostProcessInputRecursive();

    problemManager->ProblemSetup();
  }

  /**
   * @brief Destructor.
   */
  static void TearDownTestCase()
  {
    delete problemManager;
    problemManager = nullptr;
  }

  static ProblemManager * problemManager;
};

ProblemManager * LAIHelperFunctionsTest::problemManager = nullptr; //!< the main problemManager.

TEST_F(LAIHelperFunctionsTest, Test_NodalVectorPermutation)
{
  DomainPartition * const domain = problemManager->getDomainPartition();

}


/**
 * @function main
 * @brief Main function to setup the GEOSX environment, read the xml file and run all cases.
 */
int main( int argc, char** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  // Global call will not work because CXXUtils has already been initialized in problemManager
  geosx::basicSetup( argc, argv );

  global_argc = argc;
  global_argv = new char*[static_cast<unsigned int>( global_argc )];
  for( int i = 0 ; i < argc ; ++i )
  {
    global_argv[i] = argv[i];
  }

  int const result = RUN_ALL_TESTS();

  delete[] global_argv;

  // Global call will not work because CXXUtils will be destructed by problemManager
  geosx::basicCleanup();

  return result;
}




