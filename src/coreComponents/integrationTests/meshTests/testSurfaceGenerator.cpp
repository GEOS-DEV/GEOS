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
 * @file testSurfaceGenerator.cpp
 * @brief Integration tests for SurfaceGenerator topology modification (serial, no MPI).
 *
 * Verifies node duplication, fracture element creation, and Euler characteristic (χ = V-E+F-C)
 * for boundary-cutting and internal fracture configurations on hex and tet meshes.
 */

#include "testSurfaceGeneratorCommon.hpp"

// Global flag to control debug printing (set to false to disable verbose output)
CommandLineOptions g_commandLineOptions;

/// @brief Directory containing the test binary (and the copied mesh files).
///        Populated in main() from argv[0] so it is always correct regardless
///        of the build system (Makefiles, Ninja, Xcode) or how CTest sets CWD.
std::string g_testBinaryDir;

/**
 * @brief Parametrized test fixture for SurfaceGenerator.
 *
 * Parameters: (test-name, mesh-file, node-set-names, χ_before, χ_after).
 */
class SurfaceGeneratorTest
  : public ::testing::TestWithParam< std::tuple< std::string, std::string, std::string, integer, integer > >
{
protected:
  void SetUp() override
  {
    // Priority: runtime argv[0]-derived dir > env var > compile-time macro.
    // g_testBinaryDir is always the actual directory of the executable, so it
    // is correct even when the Xcode generator bakes the wrong path into the
    // TEST_BINARY_DIR macro or when CTest sets a different CWD.
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
                                std::string const & nodeSetNames )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" file=")xml" << meshFile << R"xml(" nodesetNames=")xml" << nodeSetNames << R"xml(" partitionRefinement="0"/>
  </Mesh>
  
  <Solvers gravityVector="{0.0, 0.0, 0.0}">
    <SurfaceGenerator 
      name="SurfaceGen" 
      targetRegions="{ Region, Fracture }"
      fractureRegion="Fracture"
      initialRockToughness="1.0" 
      logLevel="0"/>
  </Solvers>
  
  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ emptyConstitutive }"/>
    <SurfaceElementRegion 
      name="Fracture" 
      defaultAperture="1.0e-4"
      faceBlock="faceElementSubRegion"
      materialList="{ emptyConstitutive }"/>
  </ElementRegions>
  
  <Constitutive>
    <NullModel name="emptyConstitutive"/>
  </Constitutive>
  
  <FieldSpecifications>
    <FieldSpecification 
      name="separableFace" 
      fieldName="isFaceSeparable" 
      initialCondition="1" 
      setNames=")xml" << nodeSetNames << R"xml(" 
      objectPath="faceManager" 
      scale="1"/>
    <FieldSpecification 
      name="frac" 
      initialCondition="1" 
      setNames=")xml" << nodeSetNames << R"xml(" 
      objectPath="faceManager" 
      fieldName="ruptureState" 
      scale="1"/>
  </FieldSpecifications>
  <Events maxTime="1.0e-10">
    <SoloEvent name="preFracture" target="/Solvers/SurfaceGen"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }

  /// @brief Execute one parametrized test case end-to-end.
  void runTest( std::string const & testCaseName,
                std::string const & meshFileName,
                std::string const & nodeSetNames,
                integer const expectedEulerBefore,
                integer const expectedEulerAfter )
  {
    // Prefix the mesh path with testBinaryDir so GEOS can find the .vtu files
    // regardless of the CWD set by CTest (which differs from TARGET_FILE_DIR).
    std::string const xmlInput = generateXmlInput(
      testBinaryDir + "/" + meshFileName,
      nodeSetNames );

    // Write XML next to the binary (testBinaryDir) so the path is always valid
    // whether the test is run directly or through CTest.
    // Only rank 0 writes; all ranks barrier before any rank tries to open it.
    std::string const xmlFileName = testBinaryDir + "/test_surface_gen_" + testCaseName + ".xml";
    if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
    {
      std::ofstream xmlFile( xmlFileName );
      xmlFile << xmlInput;
      xmlFile.close();
    }
    MpiWrapper::barrier( MPI_COMM_GEOS );

    // Setup GEOS
    std::unique_ptr< CommandLineOptions > options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
    options->inputFileNames.clear();
    options->inputFileNames.push_back( xmlFileName );

    // Scoped state to ensure GeosxState is fully destroyed before the next
    // test case constructs a new one (GeosxState is a singleton).
    {
      GeosxState state( std::move( options ) );
      ASSERT_TRUE( state.initializeDataRepository() )
        << "Test " << testCaseName << ": Failed to initialize data repository for mesh file '"
        << meshFileName << "' with node sets " << nodeSetNames;
      state.applyInitialConditions();

      // Get mesh before surface generation
      ProblemManager & pm = state.getProblemManager();
      DomainPartition & domain = pm.getDomainPartition();
      MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

      NodeManager & nodeManager = mesh.getNodeManager();
      EdgeManager & edgeManager = mesh.getEdgeManager();
      FaceManager & faceManager = mesh.getFaceManager();
      ElementRegionManager & elemManager = mesh.getElemManager();

    // Store initial statistics
    TopologyStats statsBefore;
    statsBefore.numNodes = nodeManager.size();
    statsBefore.numEdges = edgeManager.size();
    statsBefore.numFaces = faceManager.size();
    statsBefore.numCells = 0;
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
    {
      statsBefore.numCells += subRegion.size();
    } );

    // **NEW: Preprocess fracture topology to predict duplication BEFORE mesh split**
    std::vector< std::string > fractureNodeSets = parseNodeSetNames( nodeSetNames );
    ExpectedDuplication expected = preprocessFractureTopology( mesh, fractureNodeSets );
    
    GEOS_LOG_RANK_0( "========================================" );
    GEOS_LOG_RANK_0( "Test: " << testCaseName );
    GEOS_LOG_RANK_0( "Expected duplication:" );
    GEOS_LOG_RANK_0( "  Nodes to split: " << expected.numNodesToDuplicate
                     << " (new nodes created: " << expected.totalDuplicatedNodes << ")" );
    GEOS_LOG_RANK_0( "  Edges to split: " << expected.numEdgesToDuplicate
                     << " (new edges created: " << expected.totalDuplicatedEdges << ")" );
    GEOS_LOG_RANK_0( "  Faces to split: " << expected.numFacesToDuplicate
                     << " (new faces created: " << expected.totalDuplicatedFaces << ")" );

    // Compute Euler-Poincaré characteristic χ = V - E + F - C (should be 1 for connected solid).
    integer const eulerCharBeforeSplit = computeEulerCharacteristic( nodeManager, edgeManager, faceManager, elemManager );
    GEOS_LOG_RANK_0( "  Euler characteristic before split: " << eulerCharBeforeSplit );
    
    // Run the simulation (executes SurfaceGenerator and splits the mesh)
    state.run();

    // Get statistics after surface generation
    TopologyStats statsAfter;
    statsAfter.numNodes = nodeManager.size();
    statsAfter.numEdges = edgeManager.size();
    statsAfter.numFaces = faceManager.size();
    statsAfter.numCells = 0;
    statsAfter.numFractureElements = 0;
    
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
    {
      statsAfter.numCells += subRegion.size();
    } );
    
    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
    {
      statsAfter.numFractureElements += subRegion.size();
    } );

    statsAfter.numDuplicatedNodes = statsAfter.numNodes - statsBefore.numNodes;

    // **NEW: AFTER mesh split - Compare expeted vs computed values**
    GEOS_LOG_RANK_0( "Actual duplication:" );
    GEOS_LOG_RANK_0( "  New nodes created: " << statsAfter.numDuplicatedNodes );
    GEOS_LOG_RANK_0( "  Fracture elements created: " << statsAfter.numFractureElements );
    
    // Compute Euler-Poincaré characteristic χ = V - E + F - C after split
    integer const eulerCharAfterSplit = computeEulerCharacteristic( nodeManager, edgeManager, faceManager, elemManager );
    GEOS_LOG_RANK_0( "  Euler characteristic after split: " << eulerCharAfterSplit );
      
    // Store Euler characteristic
    statsAfter.eulerCharacteristic = eulerCharAfterSplit;

    // All validations are performed in one place
    validateSurfaceGeneratorResults( testCaseName,
                                     meshFileName,
                                     nodeSetNames,
                                     statsBefore,
                                     statsAfter,
                                     expected,
                                     expectedEulerBefore,
                                     eulerCharBeforeSplit,
                                     expectedEulerAfter,
                                     eulerCharAfterSplit,
                                     nodeManager );

    // Log statistics
    GEOS_LOG_RANK_0( "Summary:" );
    GEOS_LOG_RANK_0( "  Nodes before: " << statsBefore.numNodes << ", after: " << statsAfter.numNodes
                     << " (+" << statsAfter.numDuplicatedNodes << ")" );
    GEOS_LOG_RANK_0( "  Fracture elements: " << statsAfter.numFractureElements );
    GEOS_LOG_RANK_0( "========================================" );

    } // end scoped GeosxState
  }
};

/// @brief Test body — dispatches to runTest().
/// Only runs on exactly 1 MPI rank; skipped when the binary is launched under mpirun with >1 ranks.
TEST_P( SurfaceGeneratorTest, TopologyValidation )
{
  auto const & params = GetParam();
  std::string const & testCaseName = std::get< 0 >( params );
  std::string const & meshFileName = std::get< 1 >( params );
  std::string const & nodeSetNames = std::get< 2 >( params );
  integer const expectedEulerBefore = std::get< 3 >( params );
  integer const expectedEulerAfter = std::get< 4 >( params );

  runTest( testCaseName, meshFileName, nodeSetNames, expectedEulerBefore, expectedEulerAfter );
}

// ---------------------------------------------------------------------------
// NoBoundaryCut: single interior node at (0.5,0.5,0.5) is duplicated once per
// fracture plane sharing it, but Δχ = +1 regardless of additional planes → χ_after = 2.
// ---------------------------------------------------------------------------
// clang-format off
INSTANTIATE_TEST_SUITE_P(
  SurfaceGeneratorCases,
  SurfaceGeneratorTest,
  ::testing::Values(

    // =======================================================================
    // dfn_market meshes
    // -----------------------------------------------------------------------
    // Naming convention:
    //   fractured[_wavy][_full_span]_mesh_{hex|tet}_DFN_{suffix}.vtu
    //
    // Column layout:  test-name,  mesh-file,  node-sets,  χ_before, χ_after
    //
    // Node-set rules:
    //   • _DFN_1   → { f1_node_set }
    //   • _DFN_2   → { f2_node_set }
    //   • _DFN_3   → { f3_node_set }
    //   • _DFN_12  → { f1_node_set, f2_node_set }
    //   • _DFN_13  → { f1_node_set, f3_node_set }
    //   • _DFN_23  → { f2_node_set, f3_node_set }
    //   • _DFN_123 → { f1_node_set, f2_node_set, f3_node_set }

    // -----------------------------------------------------------------------
    // flat · no-boundary-cut · hex
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                             node sets                         χ_b  χ_a
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_1",   "fractured_mesh_hex_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_2",   "fractured_mesh_hex_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_3",   "fractured_mesh_hex_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_12",  "fractured_mesh_hex_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_13",  "fractured_mesh_hex_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_23",  "fractured_mesh_hex_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_123", "fractured_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // flat · no-boundary-cut · tet
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                             node sets                         χ_b  χ_a
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_1",   "fractured_mesh_tet_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_2",   "fractured_mesh_tet_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_3",   "fractured_mesh_tet_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_12",  "fractured_mesh_tet_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_13",  "fractured_mesh_tet_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_23",  "fractured_mesh_tet_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_123", "fractured_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // flat · full-span · hex
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                                      node sets                        χ_b  χ_a
    std::make_tuple( "Mkt_BndCut_hex_DFN_1",   "fractured_full_span_mesh_hex_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_2",   "fractured_full_span_mesh_hex_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_3",   "fractured_full_span_mesh_hex_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_12",  "fractured_full_span_mesh_hex_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_13",  "fractured_full_span_mesh_hex_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_23",  "fractured_full_span_mesh_hex_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_123", "fractured_full_span_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // flat · full-span · tet
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                                      node sets                        χ_b  χ_a
    std::make_tuple( "Mkt_BndCut_tet_DFN_1",   "fractured_full_span_mesh_tet_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_2",   "fractured_full_span_mesh_tet_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_3",   "fractured_full_span_mesh_tet_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_12",  "fractured_full_span_mesh_tet_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_13",  "fractured_full_span_mesh_tet_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_23",  "fractured_full_span_mesh_tet_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_123", "fractured_full_span_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // wavy · no-boundary-cut · hex
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                  node sets                        χ_b  χ_a
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_1",   "fractured_wavy_mesh_hex_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_2",   "fractured_wavy_mesh_hex_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_3",   "fractured_wavy_mesh_hex_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_12",  "fractured_wavy_mesh_hex_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_13",  "fractured_wavy_mesh_hex_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_23",  "fractured_wavy_mesh_hex_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_123", "fractured_wavy_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // wavy · no-boundary-cut · tet
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                  node sets                        χ_b  χ_a
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_1",   "fractured_wavy_mesh_tet_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_2",   "fractured_wavy_mesh_tet_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_3",   "fractured_wavy_mesh_tet_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_12",  "fractured_wavy_mesh_tet_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_13",  "fractured_wavy_mesh_tet_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_23",  "fractured_wavy_mesh_tet_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_123", "fractured_wavy_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // wavy · full-span · hex
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                           node sets                       χ_b  χ_a
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_1",   "fractured_wavy_full_span_mesh_hex_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_2",   "fractured_wavy_full_span_mesh_hex_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_3",   "fractured_wavy_full_span_mesh_hex_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_12",  "fractured_wavy_full_span_mesh_hex_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_13",  "fractured_wavy_full_span_mesh_hex_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_23",  "fractured_wavy_full_span_mesh_hex_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_123", "fractured_wavy_full_span_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // wavy · full-span · tet
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                           node sets                       χ_b  χ_a
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_1",   "fractured_wavy_full_span_mesh_tet_DFN_1.vtu",   "{ f1_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_2",   "fractured_wavy_full_span_mesh_tet_DFN_2.vtu",   "{ f2_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_3",   "fractured_wavy_full_span_mesh_tet_DFN_3.vtu",   "{ f3_node_set }",                           1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_12",  "fractured_wavy_full_span_mesh_tet_DFN_12.vtu",  "{ f1_node_set, f2_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_13",  "fractured_wavy_full_span_mesh_tet_DFN_13.vtu",  "{ f1_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_23",  "fractured_wavy_full_span_mesh_tet_DFN_23.vtu",  "{ f2_node_set, f3_node_set }",              1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_123", "fractured_wavy_full_span_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // T-shaped · boundary-cutting · hex & tet
    // -----------------------------------------------------------------------
    // F1 fully spans the domain (boundary-cutting); F2 is the stem terminating at F1.
    // The junction splits the solid into 3 parts → χ_after = 3.
    //                   test name                          mesh file                             node sets             χ_b  χ_a
    std::make_tuple( "Mkt_BndCut_t_shaped_hex_DFN_12", "t_shaped_wavy_mesh_hex_DFN_t1t2.vtu", "{ f1_node_set, f2_node_set }", 1, 3 ),
    std::make_tuple( "Mkt_BndCut_t_shaped_tet_DFN_12", "t_shaped_wavy_mesh_tet_DFN_t1t2.vtu", "{ f1_node_set, f2_node_set }", 1, 3 ),

    // -----------------------------------------------------------------------
    // Y-shaped · boundary-cutting · hex & tet
    // -----------------------------------------------------------------------
    // Three boundary-cutting fractures meeting along a common line (Y cross-section).
    // Splits the solid into 3 parts → χ_after = 3.
    //                   test name                          mesh file                                  node sets                          χ_b  χ_a
    std::make_tuple( "Mkt_BndCut_Y_shaped_hex_DFN_123", "y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 3 ),
    std::make_tuple( "Mkt_BndCut_Y_shaped_tet_DFN_123", "y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 3 )

  )
);

int main( int argc, char * argv[] )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv, false );
  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
