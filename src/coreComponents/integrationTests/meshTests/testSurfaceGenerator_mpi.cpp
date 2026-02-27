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
 * @file testSurfaceGenerator_mpi.cpp
 * @brief Integration tests for SurfaceGenerator topology modification (MPI parallel).
 *
 * Verifies node duplication, fracture element creation, and Euler characteristic (χ = V-E+F-C)
 * for boundary-cutting and internal fracture configurations on hex and tet meshes, running
 * with various partitioning configurations.
 */

#include "testSurfaceGeneratorCommon.hpp"

// Global flag to control debug printing (set to false to disable verbose output)
CommandLineOptions g_commandLineOptions;

/// @brief Directory containing the test binary (and the copied mesh files).
///        Populated in main() from argv[0] so it is always correct regardless
///        of the build system (Makefiles, Ninja, Xcode) or how CTest sets CWD.
std::string g_testBinaryDir;

/**
 * @brief MPI variant of the SurfaceGenerator integration test.
 *
 * Same topology verification as SurfaceGeneratorTest, but the domain is
 * distributed across 4 MPI ranks using a user-supplied partition tuple.
 *
 * Parametrized with
 *   std::tuple< std::string,           // test case name
 *               std::string,           // mesh file name
 *               std::string,           // node-set names
 *               integer,               // expected Euler characteristic before split
 *               integer,               // expected Euler characteristic after split
 *               localIndex,            // expected new nodes created (global, from serial ground truth)
 *               localIndex >           // expected fracture elements created (global, from serial ground truth)
 *
 * Because ::testing::Combine is used, the actual ParamType produced by GTest is
 *   std::tuple< MeshTuple, PartitionTuple >
 * where
 *   MeshTuple      = std::tuple< std::string, std::string, std::string, integer, integer, localIndex, localIndex >
 *   PartitionTuple = std::tuple< int, int, int >
 *
 * The expected node/fracture-element counts are hard-coded from the verified serial run.
 * Computing them from partitioned local topology (preprocessFractureTopology) is unreliable
 * because fracture-boundary classification depends on global topology, not the local partition.
 */
using SurfaceGenerator_mpiMeshTuple = std::tuple< std::string, std::string, std::string, integer, integer, localIndex, localIndex >;
using SurfaceGenerator_mpiPartitionTuple = std::tuple< int, int, int >;

class SurfaceGenerator_mpiTest
  : public ::testing::TestWithParam<
    std::tuple< SurfaceGenerator_mpiMeshTuple, SurfaceGenerator_mpiPartitionTuple > >
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

  /// @brief Reuse the same XML generator as the serial test.
  std::string generateXmlInput( std::string const & meshFile,
                                std::string const & nodeSetNames )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" file=")xml" << meshFile << R"xml(" nodesetNames=")xml" << nodeSetNames <<
      R"xml(" partitionMethod="parmetis"/>
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
    <SurfaceElementRegion name="Fracture" defaultAperture="1.0e-4" faceBlock="faceElementSubRegion" materialList="{ emptyConstitutive }"/>
  </ElementRegions>
  <Constitutive>
    <NullModel name="emptyConstitutive"/>
  </Constitutive>
  <FieldSpecifications>
    <FieldSpecification name="separableFace" fieldName="isFaceSeparable" initialCondition="1"
      setNames=")xml" << nodeSetNames << R"xml(" objectPath="faceManager" scale="1"/>
    <FieldSpecification name="frac" initialCondition="1"
      setNames=")xml" << nodeSetNames <<
      R"xml(" objectPath="faceManager"
      fieldName="ruptureState" scale="1"/>
  </FieldSpecifications>
  <Events maxTime="1.0e-10">
    <SoloEvent name="preFracture" target="/Solvers/SurfaceGen"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }
};

/**
 * @brief MPI variant of TopologyValidation.
 *
 * Replicates every assertion from the serial SurfaceGeneratorTest::runTest():
 *   A1  – node duplication count matches preprocessFractureTopology prediction
 *   A2  – fracture element count matches prediction
 *   A3  – no NaN node coordinates after split
 *   A4  – Euler characteristic before split equals expected value
 *   A5  – Euler characteristic after  split equals expected value
 *   A6  – no degenerate elements (duplicate node references)
 *
 * The Euler characteristic is computed via computeEulerCharacteristic
 * (manager-sizes production version from MeshGeneratorBase) both before
 * and after the split.
 */
TEST_P( SurfaceGenerator_mpiTest, TopologyValidation )
{
  auto const & params = GetParam();
  // ::testing::Combine produces std::tuple< MeshTuple, PartitionTuple >
  SurfaceGenerator_mpiMeshTuple const & meshParams = std::get< 0 >( params );
  SurfaceGenerator_mpiPartitionTuple const & partitions = std::get< 1 >( params );

  std::string const & testCaseName        = std::get< 0 >( meshParams );
  std::string const & meshFileName        = std::get< 1 >( meshParams );
  std::string const & nodeSetNames        = std::get< 2 >( meshParams );
  integer const expectedEulerBefore       = std::get< 3 >( meshParams );
  integer const expectedEulerAfter        = std::get< 4 >( meshParams );
  localIndex const expectedNewNodes       = std::get< 5 >( meshParams );
  localIndex const expectedFractureElems  = std::get< 6 >( meshParams );

  int const xPartitions = std::get< 0 >( partitions );
  int const yPartitions = std::get< 1 >( partitions );
  int const zPartitions = std::get< 2 >( partitions );

  // Check if we are running with the correct number of ranks for this partitioning.
  // This prevents the serial test runner (rank=1) from crashing when attempting
  // to run a 4-partition test.
  int const requiredRanks = xPartitions * yPartitions * zPartitions;
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) != requiredRanks )
  {
    GTEST_SKIP() << "Skipping MPI test expecting " << requiredRanks
                 << " ranks (running on " << MpiWrapper::commSize( MPI_COMM_GEOS ) << ")";
  }

  std::string const xmlInput = generateXmlInput( testBinaryDir + "/" + meshFileName, nodeSetNames );

  // Only rank 0 writes the XML; all ranks barrier before reading.
  std::string const xmlPath = testBinaryDir
                              + "/test_surface_gen_mpi_" + testCaseName
                              + "_" + std::to_string( xPartitions )
                              + "x" + std::to_string( yPartitions )
                              + "x" + std::to_string( zPartitions )
                              + ".xml";
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
  {
    std::ofstream ofs( xmlPath );
    ofs << xmlInput;
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );

  auto options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
  options->inputFileNames.push_back( xmlPath );
  options->problemName = "test_surface_gen_mpi_" + testCaseName;
  options->xPartitionsOverride = xPartitions;
  options->yPartitionsOverride = yPartitions;
  options->zPartitionsOverride = zPartitions;
  options->overridePartitionNumbers = true;

  // Scoped state to ensure GeosxState is fully destroyed before the next
  // test case constructs a new one (GeosxState is a singleton).
  {
    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() )
      << "Test " << testCaseName << ": Failed to initialize data repository for '"
      << meshFileName << "' with partitioning "
      << xPartitions << "x" << yPartitions << "x" << zPartitions;
    state.applyInitialConditions();

    // ------------------------------------------------------------------
    // Obtain managers from the base mesh level
    // ------------------------------------------------------------------
    ProblemManager & pm = state.getProblemManager();
    MeshLevel & mesh = pm.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization();
    NodeManager & nodeManager = mesh.getNodeManager();
    EdgeManager & edgeManager = mesh.getEdgeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // ------------------------------------------------------------------
    // Pre-split topology stats
    // ------------------------------------------------------------------
    TopologyStats statsBefore;
    statsBefore.numNodes = nodeManager.size();
    statsBefore.numEdges = edgeManager.size();
    statsBefore.numFaces = faceManager.size();
    statsBefore.numCells = 0;
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sr )
    {
      statsBefore.numCells += sr.size();
    } );

    // Snapshot owned-node count BEFORE the split.
    // getNumberOfLocalIndices() returns nodes with ghostRank <= -1 (owned by this rank).
    // This is the correct baseline for computing new-owned-node delta after the split.
    localIndex const ownedNodeCountBefore = nodeManager.getNumberOfLocalIndices();

    // ------------------------------------------------------------------
    // Build expected-duplication from hard-coded serial ground truth.
    //
    // preprocessFractureTopology is unreliable for partitioned meshes because
    // fracture-boundary classification depends on the *global* fracture topology
    // (which edges are on the fracture tip). On a local partition each rank only
    // sees a subset of the fracture; tip edges shared with a ghost face look
    // interior (count == 2) instead of boundary (count == 1), causing tip nodes
    // to be mis-classified as interior nodes and over-counting new nodes.
    //
    // The only reliable source for the global expected counts is the verified
    // serial run, so they are baked into the test tuple parameters.
    // ------------------------------------------------------------------
    ExpectedDuplication expected;
    expected.totalDuplicatedNodes = expectedNewNodes;
    expected.numFractureElements  = expectedFractureElems;

    GEOS_LOG_RANK_0( "========================================" );
    GEOS_LOG_RANK_0( "Test (MPI): " << testCaseName
                                    << " [" << xPartitions << "x" << yPartitions << "x" << zPartitions << "]" );
    GEOS_LOG_RANK_0( "Expected duplication (from serial ground truth):" );
    GEOS_LOG_RANK_0( "  New nodes:         " << expected.totalDuplicatedNodes );
    GEOS_LOG_RANK_0( "  Fracture elements: " << expected.numFractureElements );

    // ------------------------------------------------------------------
    // Euler characteristic BEFORE split
    // ------------------------------------------------------------------
    integer const eulerCharBeforeSplit =
      computeEulerCharacteristic( nodeManager, edgeManager, faceManager, elemManager );
    GEOS_LOG_RANK_0( "  Euler χ before split: " << eulerCharBeforeSplit );

    statsBefore.eulerCharacteristic = eulerCharBeforeSplit;

    // ------------------------------------------------------------------
    // Run SurfaceGenerator (the split)
    // ------------------------------------------------------------------
    state.run();

    // ------------------------------------------------------------------
    // Post-split topology stats
    // ------------------------------------------------------------------
    TopologyStats statsAfter;
    statsAfter.numNodes = nodeManager.size();
    statsAfter.numEdges = edgeManager.size();
    statsAfter.numFaces = faceManager.size();
    statsAfter.numCells = 0;
    statsAfter.numFractureElements = 0;

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & sr )
    {
      statsAfter.numCells += sr.size();
    } );
    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & sr )
    {
      // Use getNumberOfLocalIndices() instead of size() to count only owned
      // fracture elements (ghostRank <= -1). Summing sr.size() across ranks
      // would double-count ghost copies that appear on multiple ranks.
      statsAfter.numFractureElements += sr.getNumberOfLocalIndices();
    } );

    // ------------------------------------------------------------------
    // Count only OWNED new nodes (ghostRank <= -1) to avoid double-counting
    // ghost copies that appear on multiple ranks after the split.
    // ownedNodeCountBefore was snapshotted before state.run() so no
    // ghost-count assumptions are needed.
    // ------------------------------------------------------------------
    localIndex const ownedNodesAfter    = nodeManager.getNumberOfLocalIndices();
    localIndex const localOwnedNewNodes = ownedNodesAfter - ownedNodeCountBefore;

    statsAfter.numDuplicatedNodes = localOwnedNewNodes;

    integer const eulerCharAfterSplit =
      computeEulerCharacteristic( nodeManager, edgeManager, faceManager, elemManager );
    GEOS_LOG_RANK_0( "  Euler χ after  split:              " << eulerCharAfterSplit );

    statsAfter.eulerCharacteristic = eulerCharAfterSplit;

    // ------------------------------------------------------------------
    // Reduce actual counts to global totals.
    //
    // • localOwnedNewNodes  — each node is owned by exactly one rank, so
    //   MPI_SUM gives the true global new-node count with no double-counting.
    // • statsAfter.numFractureElements — each face element is owned by one rank; SUM is safe.
    // • expected values are global constants from the serial ground truth; no reduction needed.
    // ------------------------------------------------------------------
    localIndex const globalActualNewNodes =
      MpiWrapper::sum( localOwnedNewNodes, MPI_COMM_GEOS );
    localIndex const globalPredictedNewNodes = expected.totalDuplicatedNodes;
    localIndex const globalFractureElements =
      MpiWrapper::sum( statsAfter.numFractureElements, MPI_COMM_GEOS );
    localIndex const globalPredictedFractureElements = expected.numFractureElements;
    localIndex const globalNodesBefore =
      MpiWrapper::sum( ownedNodeCountBefore, MPI_COMM_GEOS );

    GEOS_LOG_RANK_0( "Actual duplication (global):" );
    GEOS_LOG_RANK_0( "  New nodes created:     " << globalActualNewNodes );
    GEOS_LOG_RANK_0( "  Fracture elements:     " << globalFractureElements );
    GEOS_LOG_RANK_0( "========================================" );

    // ------------------------------------------------------------------
    // All GTest assertions MUST fire on rank 0 only.
    // Running EXPECT_EQ on every rank produces N identical failure messages
    // and confuses the test runner into thinking there are N failures.
    // ------------------------------------------------------------------
    if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
    {
      // Build globally-reduced structs for the shared validation helper.
      TopologyStats globalStatsBefore = statsBefore;
      globalStatsBefore.numNodes = globalNodesBefore;

      TopologyStats globalStatsAfter = statsAfter;
      globalStatsAfter.numDuplicatedNodes  = globalActualNewNodes;
      globalStatsAfter.numFractureElements = globalFractureElements;
      globalStatsAfter.numNodes            = globalNodesBefore + globalActualNewNodes;

      ExpectedDuplication globalExpected = expected;
      globalExpected.totalDuplicatedNodes = globalPredictedNewNodes;
      globalExpected.numFractureElements  = globalPredictedFractureElements;

      validateSurfaceGeneratorResults( testCaseName,
                                       meshFileName,
                                       nodeSetNames,
                                       globalStatsBefore,
                                       globalStatsAfter,
                                       globalExpected,
                                       expectedEulerBefore,
                                       eulerCharBeforeSplit,
                                       expectedEulerAfter,
                                       eulerCharAfterSplit,
                                       nodeManager );
    }

  } // end scoped GeosxState

  // Cleanup XML after state is destroyed
  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
  {
    std::remove( xmlPath.c_str() );
  }
}

// ---------------------------------------------------------------------------
// MPI test instantiation
//
// All 6 unique partitions of 4 ranks (product x*y*z == 4):
//   (1,1,4)  (1,4,1)  (4,1,1)  (1,2,2)  (2,1,2)  (2,2,1)
//
// Meshes: the 4 new junction topologies that exercise the most demanding
// splitting paths (T-junction and Y-junction, hex and tet).
//
// Meshes: every multi-fracture (_DFN_123) configuration from the serial suite
// plus the T-junction and Y-junction topologies.
// Expected node/element counts are serial ground truth (all tests pass serially).
// ---------------------------------------------------------------------------
// clang-format off
INSTANTIATE_TEST_SUITE_P(
  SurfaceGenerator_mpiCases,
  SurfaceGenerator_mpiTest,
  ::testing::Combine(
    ::testing::Values(
      // -----------------------------------------------------------------------
      // flat · no-boundary-cut · hex · DFN_123     nodes:   7  elems:  12
      std::make_tuple( "Mkt_NoBndCut_hex_DFN_123",
                       "fractured_mesh_hex_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 2, localIndex( 7 ), localIndex( 12 ) ),
      // flat · no-boundary-cut · tet · DFN_123     nodes:  19  elems:  48
      std::make_tuple( "Mkt_NoBndCut_tet_DFN_123",
                       "fractured_mesh_tet_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 2, localIndex( 19 ), localIndex( 48 ) ),
      // -----------------------------------------------------------------------
      // flat · full-span · hex · DFN_123           nodes:  37  elems:  12
      std::make_tuple( "Mkt_BndCut_hex_DFN_123",
                       "fractured_full_span_mesh_hex_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 8, localIndex( 37 ), localIndex( 12 ) ),
      // flat · full-span · tet · DFN_123           nodes: 127  elems: 168
      std::make_tuple( "Mkt_BndCut_tet_DFN_123",
                       "fractured_full_span_mesh_tet_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 8, localIndex( 127 ), localIndex( 168 ) ),
      // -----------------------------------------------------------------------
      // wavy · no-boundary-cut · hex · DFN_123     nodes:   7  elems:  12
      std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_123",
                       "fractured_wavy_mesh_hex_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 2, localIndex( 7 ), localIndex( 12 ) ),
      // wavy · no-boundary-cut · tet · DFN_123     nodes:  19  elems:  48
      std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_123",
                       "fractured_wavy_mesh_tet_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 2, localIndex( 19 ), localIndex( 48 ) ),
      // -----------------------------------------------------------------------
      // wavy · full-span · hex · DFN_123           nodes:  37  elems:  12
      std::make_tuple( "Mkt_WavyBndCut_hex_DFN_123",
                       "fractured_wavy_full_span_mesh_hex_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 8, localIndex( 37 ), localIndex( 12 ) ),
      // wavy · full-span · tet · DFN_123           nodes: 127  elems: 168
      std::make_tuple( "Mkt_WavyBndCut_tet_DFN_123",
                       "fractured_wavy_full_span_mesh_tet_DFN_123.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 8, localIndex( 127 ), localIndex( 168 ) ),
      // -----------------------------------------------------------------------
      // T-shaped · boundary-cutting · hex          nodes:  15  elems:   6
      std::make_tuple( "Mkt_BndCut_t_shaped_hex_DFN_12",
                       "t_shaped_wavy_mesh_hex_DFN_t1t2.vtu",
                       "{ f1_node_set, f2_node_set }",
                       1, 3, localIndex( 15 ), localIndex( 6 ) ),
      // T-shaped · boundary-cutting · tet          nodes:  49  elems:  66
      std::make_tuple( "Mkt_BndCut_t_shaped_tet_DFN_12",
                       "t_shaped_wavy_mesh_tet_DFN_t1t2.vtu",
                       "{ f1_node_set, f2_node_set }",
                       1, 3, localIndex( 49 ), localIndex( 66 ) ),
      // -----------------------------------------------------------------------
      // Y-shaped · boundary-cutting · hex          nodes:  15  elems:   6
      std::make_tuple( "Mkt_BndCut_Y_shaped_hex_DFN_123",
                       "y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 3, localIndex( 15 ), localIndex( 6 ) ),
      // Y-shaped · boundary-cutting · tet          nodes:  49  elems:  66
      std::make_tuple( "Mkt_BndCut_Y_shaped_tet_DFN_123",
                       "y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set }",
                       1, 3, localIndex( 49 ), localIndex( 66 ) ),

      // Miscellaneous · no-boundary-cutting ·  hex          nodes:  180  elems:  196
      std::make_tuple( "Mkt_NoBndCut_5_fracs_hex_DFN",
                       "DFN_5_fractures_hex_binarized.vtu",
                       "{  f1_node_set, f2_node_set, f3_node_set, f4_node_set, f5_node_set }",
                       1, 2, localIndex( 180 ), localIndex( 196 ) ),

      // Miscellaneous · no-boundary-cutting · tet          nodes:  180  elems:  392
      // Note on χ_a = 12. χ_a should be 2 for a full topological construction of finite element spaces of order 2 or higher.
      // The simplest way to think about this is to split a single fracture discretized with two triangles when no high-order elements are used.
      // The Euler characteristic will then be incorrectly defined, since the internal edge will not be duplicated, interesting!.
      // This issue is ubiquitous in any workflow used for mesh splitting; whether internal or external to GEOS.
      // More importantly, this is a good example where the topological features of the mesh differ from the topological features of the function
      // spaces defined on the given mesh.
      std::make_tuple( "Mkt_NoBndCut_5_fracs_tet_DFN",
                       "DFN_5_fractures_tet_binarized.vtu",
                       "{ f1_node_set, f2_node_set, f3_node_set, f4_node_set, f5_node_set }",
                       1, 12, localIndex( 180 ), localIndex( 392 ) )

      ),
    ::testing::Values(
      std::make_tuple( 1, 1, 4 ),
      std::make_tuple( 1, 4, 1 ),
      std::make_tuple( 4, 1, 1 ),
      std::make_tuple( 1, 2, 2 ),
      std::make_tuple( 2, 1, 2 ),
      std::make_tuple( 2, 2, 1 )
      )
    )
  );


int main( int argc, char * argv[] )
{
  // Derive the directory containing this executable from argv[0].
  // dirname() modifies its argument, so work on a copy.
  // This is the most reliable way to find co-located mesh files regardless
  // of what CWD CTest or the shell has set.
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
