/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testSplitMesh.cpp
 * @brief Integration tests for loading meshes with predefined fractures via VTKMesh faceBlocks (serial).
 *
 * Verifies that a mesh read with pre-existing fracture surfaces embedded as faceBlocks contains
 * the correct number of surface elements, valid node coordinates, and the expected
 * Euler characteristic (χ = V-E+F-C).
 */

#include "testSplitMeshCommon.hpp"

CommandLineOptions g_commandLineOptions;

/// @brief Directory containing the test binary (and the copied mesh files).
///        Populated in main() from argv[0] so it is always correct regardless
///        of the build system (Makefiles, Ninja, Xcode) or how CTest sets CWD.
std::string g_testBinaryDir;

/**
 * @brief Parametrized test fixture for loading split (predefined-fracture) meshes.
 *
 * Parameters:
 *   (test-name, vtm-file, faceBlock-name, expected-surface-elements, expected-Euler-χ)
 */
class SplitMeshTest
  : public ::testing::TestWithParam<
    std::tuple< std::string,   // test case name
                std::string,   // .vtm mesh file (relative to testBinaryDir)
                std::string,   // faceBlock name (e.g. "Fault")
                localIndex,    // expected number of surface elements
                integer > >    // expected Euler characteristic
{
protected:
  void SetUp() override
  {
    if( !g_testBinaryDir.empty() )
    {
      testBinaryDir = g_testBinaryDir;
    }
    else
    {
      char const * envDir = std::getenv( "TEST_BINARY_DIR" );
      testBinaryDir = ( envDir != nullptr && envDir[0] != '\0' ) ? envDir : TEST_BINARY_DIR;
    }
  }

  std::string testBinaryDir;

  /// @brief Build the XML input string for this test case.
  std::string generateXmlInput( std::string const & meshFile,
                                std::string const & faceBlockName )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" checkEulerCharacteristic="1" useGlobalIds="1"
             file=")xml" << meshFile << R"xml("
             faceBlocks="{ )xml" << faceBlockName <<
      R"xml( }"/>
  </Mesh>
  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ emptyConstitutive }"/>
    <SurfaceElementRegion name=")xml" << faceBlockName << R"xml(" defaultAperture="1.0e-4"
                          faceBlock=")xml" << faceBlockName <<
      R"xml("
                          materialList="{ emptyConstitutive }"/>
  </ElementRegions>
  <Constitutive>
    <NullModel name="emptyConstitutive"/>
  </Constitutive>
</Problem>
)xml";
    return oss.str();
  }

  /// @brief Execute one parametrized test case end-to-end.
  void runTest( std::string const & testCaseName,
                std::string const & meshFileName,
                std::string const & faceBlockName,
                localIndex const expectedSurfaceElements,
                integer const expectedEuler )
  {
    std::string const xmlInput = generateXmlInput(
      testBinaryDir + "/" + meshFileName, faceBlockName );

    std::string const xmlFileName = testBinaryDir + "/test_split_mesh_" + testCaseName + ".xml";
    if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
    {
      std::ofstream xmlFile( xmlFileName );
      xmlFile << xmlInput;
    }
    MpiWrapper::barrier( MPI_COMM_GEOS );

    std::unique_ptr< CommandLineOptions > options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
    options->inputFileNames.clear();
    options->inputFileNames.push_back( xmlFileName );

    {
      GeosxState state( std::move( options ) );
      ASSERT_TRUE( state.initializeDataRepository() )
        << "Test " << testCaseName << ": Failed to initialize data repository for '"
        << meshFileName << "'";
      state.applyInitialConditions();

      ProblemManager & pm = state.getProblemManager();
      DomainPartition & domain = pm.getDomainPartition();
      MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

      NodeManager const & nodeManager = mesh.getNodeManager();
      EdgeManager const & edgeManager = mesh.getEdgeManager();
      FaceManager const & faceManager = mesh.getFaceManager();
      ElementRegionManager & elemManager = mesh.getElemManager();

      SplitMeshStats stats;
      stats.numNodes = nodeManager.size();
      stats.numEdges = edgeManager.size();
      stats.numFaces = faceManager.size();

      elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sr )
      {
        stats.numCells += sr.size();
      } );

      elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & sr )
      {
        stats.numSurfaceElements += sr.size();
      } );

      stats.eulerCharacteristic =
        computeEulerCharacteristic( nodeManager, edgeManager, faceManager, elemManager );

      GEOS_LOG_RANK_0( "========================================" );
      GEOS_LOG_RANK_0( "Test: " << testCaseName );
      GEOS_LOG_RANK_0( "  Nodes:            " << stats.numNodes );
      GEOS_LOG_RANK_0( "  Surface elements: " << stats.numSurfaceElements );
      GEOS_LOG_RANK_0( "  Euler χ:          " << stats.eulerCharacteristic );

      // Run events (VTK output)
      state.run();

      validateSplitMeshResults( testCaseName, meshFileName, stats,
                                expectedSurfaceElements, expectedEuler,
                                nodeManager, elemManager );

      GEOS_LOG_RANK_0( "========================================" );
    } // end scoped GeosxState

    if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
    {
      std::remove( xmlFileName.c_str() );
    }
  }
};

/// @brief Test body — dispatches to runTest().
TEST_P( SplitMeshTest, TopologyValidation )
{
  auto const & params = GetParam();
  runTest( std::get< 0 >( params ),
           std::get< 1 >( params ),
           std::get< 2 >( params ),
           std::get< 3 >( params ),
           std::get< 4 >( params ) );
}

// ---------------------------------------------------------------------------
// Test instantiation
//
// Column layout: test-name, vtm-file, faceBlock-name,
//                expected-surface-elements, expected-Euler-χ
// ---------------------------------------------------------------------------
// clang-format off
INSTANTIATE_TEST_SUITE_P(
  SplitMeshCases,
  SplitMeshTest,
  ::testing::Values(
    //                   test name              mesh file (vtm)                              faceBlock surf-elems         χ
    std::make_tuple( "DFN_5_fracs_hex", "DFN_5_fractures_hex_binarized.vtm", "Fault", localIndex( 196 ), integer( 2 ) ),
    std::make_tuple( "DFN_5_fracs_tet", "DFN_5_fractures_tet_binarized.vtm", "Fault", localIndex( 382 ), integer( 2 ) ),
    std::make_tuple( "full_span_hex_1", "fractured_full_span_mesh_hex_DFN_1.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "full_span_hex_2", "fractured_full_span_mesh_hex_DFN_2.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "full_span_hex_3", "fractured_full_span_mesh_hex_DFN_3.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "full_span_hex_12", "fractured_full_span_mesh_hex_DFN_12.vtm", "Fault", localIndex( 8 ), integer( 4 ) ),
    std::make_tuple( "full_span_hex_13", "fractured_full_span_mesh_hex_DFN_13.vtm", "Fault", localIndex( 8 ), integer( 4 ) ),
    std::make_tuple( "full_span_hex_23", "fractured_full_span_mesh_hex_DFN_23.vtm", "Fault", localIndex( 8 ), integer( 4 ) ),
    std::make_tuple( "full_span_hex_123", "fractured_full_span_mesh_hex_DFN_123.vtm", "Fault", localIndex( 12 ), integer( 8 ) ),
    std::make_tuple( "full_span_tet_1", "fractured_full_span_mesh_tet_DFN_1.vtm", "Fault", localIndex( 44 ), integer( 2 ) ),
    std::make_tuple( "full_span_tet_2", "fractured_full_span_mesh_tet_DFN_2.vtm", "Fault", localIndex( 44 ), integer( 2 ) ),
    std::make_tuple( "full_span_tet_3", "fractured_full_span_mesh_tet_DFN_3.vtm", "Fault", localIndex( 44 ), integer( 2 ) ),
    std::make_tuple( "full_span_tet_12", "fractured_full_span_mesh_tet_DFN_12.vtm", "Fault", localIndex( 88 ), integer( 4 ) ),
    std::make_tuple( "full_span_tet_13", "fractured_full_span_mesh_tet_DFN_13.vtm", "Fault", localIndex( 88 ), integer( 4 ) ),
    std::make_tuple( "full_span_tet_23", "fractured_full_span_mesh_tet_DFN_23.vtm", "Fault", localIndex( 88 ), integer( 4 ) ),
    std::make_tuple( "full_span_tet_123", "fractured_full_span_mesh_tet_DFN_123.vtm", "Fault", localIndex( 168 ), integer( 8 ) ),
    std::make_tuple( "mesh_hex_1", "fractured_mesh_hex_DFN_1.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "mesh_hex_2", "fractured_mesh_hex_DFN_2.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "mesh_hex_3", "fractured_mesh_hex_DFN_3.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "mesh_hex_12", "fractured_mesh_hex_DFN_12.vtm", "Fault", localIndex( 8 ), integer( 2 ) ),
    std::make_tuple( "mesh_hex_13", "fractured_mesh_hex_DFN_13.vtm", "Fault", localIndex( 8 ), integer( 2 ) ),
    std::make_tuple( "mesh_hex_23", "fractured_mesh_hex_DFN_23.vtm", "Fault", localIndex( 8 ), integer( 2 ) ),
    std::make_tuple( "mesh_hex_123", "fractured_mesh_hex_DFN_123.vtm", "Fault", localIndex( 12 ), integer( 2 ) ),
    std::make_tuple( "mesh_tet_1", "fractured_mesh_tet_DFN_1.vtm", "Fault", localIndex( 14 ), integer( 2 ) ),
    std::make_tuple( "mesh_tet_2", "fractured_mesh_tet_DFN_2.vtm", "Fault", localIndex( 14 ), integer( 2 ) ),
    std::make_tuple( "mesh_tet_3", "fractured_mesh_tet_DFN_3.vtm", "Fault", localIndex( 14 ), integer( 2 ) ),
    std::make_tuple( "mesh_tet_12", "fractured_mesh_tet_DFN_12.vtm", "Fault", localIndex( 32 ), integer( 2 ) ),
    std::make_tuple( "mesh_tet_13", "fractured_mesh_tet_DFN_13.vtm", "Fault", localIndex( 32 ), integer( 2 ) ),
    std::make_tuple( "mesh_tet_23", "fractured_mesh_tet_DFN_23.vtm", "Fault", localIndex( 32 ), integer( 2 ) ),
    std::make_tuple( "mesh_tet_123", "fractured_mesh_tet_DFN_123.vtm", "Fault", localIndex( 48 ), integer( 2 ) ),
    std::make_tuple( "wavy_full_span_hex_1", "fractured_wavy_full_span_mesh_hex_DFN_1.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "wavy_full_span_hex_2", "fractured_wavy_full_span_mesh_hex_DFN_2.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "wavy_full_span_hex_3", "fractured_wavy_full_span_mesh_hex_DFN_3.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "wavy_full_span_hex_12", "fractured_wavy_full_span_mesh_hex_DFN_12.vtm", "Fault", localIndex( 8 ), integer( 4 ) ),
    std::make_tuple( "wavy_full_span_hex_13", "fractured_wavy_full_span_mesh_hex_DFN_13.vtm", "Fault", localIndex( 8 ), integer( 4 ) ),
    std::make_tuple( "wavy_full_span_hex_23", "fractured_wavy_full_span_mesh_hex_DFN_23.vtm", "Fault", localIndex( 8 ), integer( 4 ) ),
    std::make_tuple( "wavy_full_span_hex_123", "fractured_wavy_full_span_mesh_hex_DFN_123.vtm", "Fault", localIndex( 12 ), integer( 8 ) ),
    std::make_tuple( "wavy_full_span_tet_1", "fractured_wavy_full_span_mesh_tet_DFN_1.vtm", "Fault", localIndex( 44 ), integer( 2 ) ),
    std::make_tuple( "wavy_full_span_tet_2", "fractured_wavy_full_span_mesh_tet_DFN_2.vtm", "Fault", localIndex( 44 ), integer( 2 ) ),
    std::make_tuple( "wavy_full_span_tet_3", "fractured_wavy_full_span_mesh_tet_DFN_3.vtm", "Fault", localIndex( 44 ), integer( 2 ) ),
    std::make_tuple( "wavy_full_span_tet_12", "fractured_wavy_full_span_mesh_tet_DFN_12.vtm", "Fault", localIndex( 88 ), integer( 4 ) ),
    std::make_tuple( "wavy_full_span_tet_13", "fractured_wavy_full_span_mesh_tet_DFN_13.vtm", "Fault", localIndex( 88 ), integer( 4 ) ),
    std::make_tuple( "wavy_full_span_tet_23", "fractured_wavy_full_span_mesh_tet_DFN_23.vtm", "Fault", localIndex( 88 ), integer( 4 ) ),
    std::make_tuple( "wavy_full_span_tet_123", "fractured_wavy_full_span_mesh_tet_DFN_123.vtm", "Fault", localIndex( 168 ), integer( 8 ) ),
    std::make_tuple( "wavy_hex_1", "fractured_wavy_mesh_hex_DFN_1.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "wavy_hex_2", "fractured_wavy_mesh_hex_DFN_2.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "wavy_hex_3", "fractured_wavy_mesh_hex_DFN_3.vtm", "Fault", localIndex( 4 ), integer( 2 ) ),
    std::make_tuple( "wavy_hex_12", "fractured_wavy_mesh_hex_DFN_12.vtm", "Fault", localIndex( 8 ), integer( 2 ) ),
    std::make_tuple( "wavy_hex_13", "fractured_wavy_mesh_hex_DFN_13.vtm", "Fault", localIndex( 8 ), integer( 2 ) ),
    std::make_tuple( "wavy_hex_23", "fractured_wavy_mesh_hex_DFN_23.vtm", "Fault", localIndex( 8 ), integer( 2 ) ),
    std::make_tuple( "wavy_hex_123", "fractured_wavy_mesh_hex_DFN_123.vtm", "Fault", localIndex( 12 ), integer( 2 ) ),
    std::make_tuple( "wavy_tet_1", "fractured_wavy_mesh_tet_DFN_1.vtm", "Fault", localIndex( 14 ), integer( 2 ) ),
    std::make_tuple( "wavy_tet_2", "fractured_wavy_mesh_tet_DFN_2.vtm", "Fault", localIndex( 14 ), integer( 2 ) ),
    std::make_tuple( "wavy_tet_3", "fractured_wavy_mesh_tet_DFN_3.vtm", "Fault", localIndex( 14 ), integer( 2 ) ),
    std::make_tuple( "wavy_tet_12", "fractured_wavy_mesh_tet_DFN_12.vtm", "Fault", localIndex( 32 ), integer( 2 ) ),
    std::make_tuple( "wavy_tet_13", "fractured_wavy_mesh_tet_DFN_13.vtm", "Fault", localIndex( 32 ), integer( 2 ) ),
    std::make_tuple( "wavy_tet_23", "fractured_wavy_mesh_tet_DFN_23.vtm", "Fault", localIndex( 32 ), integer( 2 ) ),
    std::make_tuple( "wavy_tet_123", "fractured_wavy_mesh_tet_DFN_123.vtm", "Fault", localIndex( 48 ), integer( 2 ) ),
    std::make_tuple( "t_shaped_wavy_hex", "t_shaped_wavy_mesh_hex_DFN_t1t2.vtm", "Fault", localIndex( 6 ), integer( 3 ) ),
    std::make_tuple( "t_shaped_wavy_tet", "t_shaped_wavy_mesh_tet_DFN_t1t2.vtm", "Fault", localIndex( 66 ), integer( 3 ) ),
    std::make_tuple( "y_shaped_wavy_hex", "y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtm", "Fault", localIndex( 6 ), integer( 3 ) ),
    std::make_tuple( "y_shaped_wavy_tet", "y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtm", "Fault", localIndex( 66 ), integer( 3 ) )
    )
  );
// clang-format on

int main( int argc, char * argv[] )
{
  // Derive binary directory from argv[0] for reliable mesh-file resolution.
  if( argc > 0 && argv[0] != nullptr )
  {
    std::vector< char > exePath( argv[0], argv[0] + std::strlen( argv[0] ) + 1 );
    g_testBinaryDir = ::dirname( exePath.data() );
  }

  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv, false );
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
