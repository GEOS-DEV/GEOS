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
 * @file testSplitMesh_mpi.cpp
 * @brief Integration tests for loading meshes with predefined fractures via VTKMesh faceBlocks (MPI parallel).
 *
 * Verifies that a mesh read with pre-existing fracture surfaces embedded as faceBlocks contains
 * the correct number of surface elements, valid node coordinates, and the expected
 * Euler characteristic (χ = V-E+F-C), running with various MPI partition configurations.
 */

#include "testSplitMeshCommon.hpp"

CommandLineOptions g_commandLineOptions;

/// @brief Directory containing the test binary (and the copied mesh files).
///        Populated in main() from argv[0] so it is always correct regardless
///        of the build system (Makefiles, Ninja, Xcode) or how CTest sets CWD.
std::string g_testBinaryDir;

/**
 * @brief MPI parametrized tuple types.
 *
 * MeshTuple:
 *   (test-name, vtm-file, faceBlock-name, expected-surface-elements, expected-Euler-χ)
 * PartitionTuple: (x, y, z)
 */
using SplitMesh_mpiMeshTuple      = std::tuple< std::string, std::string, std::string, localIndex, integer >;
using SplitMesh_mpiPartitionTuple = std::tuple< int, int, int >;

class SplitMesh_mpiTest
  : public ::testing::TestWithParam<
    std::tuple< SplitMesh_mpiMeshTuple, SplitMesh_mpiPartitionTuple > >
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
             faceBlocks="{ )xml" << faceBlockName << R"xml( }"/>
  </Mesh>
  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ emptyConstitutive }"/>
    <SurfaceElementRegion name=")xml" << faceBlockName << R"xml(" defaultAperture="1.0e-4"
                          faceBlock=")xml" << faceBlockName << R"xml("
                          materialList="{ emptyConstitutive }"/>
  </ElementRegions>
  <Constitutive>
    <NullModel name="emptyConstitutive"/>
  </Constitutive>
  <Outputs>
    <VTK name="vtkOutputM" outputRegionType="cell"/>
    <VTK name="vtkOutputF" outputRegionType="surface"/>
  </Outputs>
  <Events minTime="0.0" maxTime="0.1">
    <SoloEvent name="outputsM" target="/Outputs/vtkOutputM"/>
    <SoloEvent name="outputsF" target="/Outputs/vtkOutputF"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }
};

TEST_P( SplitMesh_mpiTest, TopologyValidation )
{
  auto const & params = GetParam();
  SplitMesh_mpiMeshTuple const & meshParams = std::get< 0 >( params );
  SplitMesh_mpiPartitionTuple const & partitions = std::get< 1 >( params );

  std::string const & testCaseName          = std::get< 0 >( meshParams );
  std::string const & meshFileName          = std::get< 1 >( meshParams );
  std::string const & faceBlockName         = std::get< 2 >( meshParams );
  localIndex const expectedSurfaceElements  = std::get< 3 >( meshParams );
  integer const expectedEuler               = std::get< 4 >( meshParams );

  int const xPartitions = std::get< 0 >( partitions );
  int const yPartitions = std::get< 1 >( partitions );
  int const zPartitions = std::get< 2 >( partitions );

  int const requiredRanks = xPartitions * yPartitions * zPartitions;
  if( MpiWrapper::commSize( MPI_COMM_GEOS ) != requiredRanks )
  {
    GTEST_SKIP() << "Skipping MPI test expecting " << requiredRanks
                 << " ranks (running on " << MpiWrapper::commSize( MPI_COMM_GEOS ) << ")";
  }

  std::string const xmlInput = generateXmlInput( testBinaryDir + "/" + meshFileName, faceBlockName );

  std::string const xmlPath = testBinaryDir
                              + "/test_split_mesh_mpi_" + testCaseName
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
  options->problemName          = "test_split_mesh_mpi_" + testCaseName;
  options->xPartitionsOverride  = xPartitions;
  options->yPartitionsOverride  = yPartitions;
  options->zPartitionsOverride  = zPartitions;
  options->overridePartitionNumbers = true;

  {
    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() )
      << "Test " << testCaseName << ": Failed to initialize data repository for '"
      << meshFileName << "' with partitioning "
      << xPartitions << "x" << yPartitions << "x" << zPartitions;
    state.applyInitialConditions();

    ProblemManager & pm = state.getProblemManager();
    MeshLevel & mesh = pm.getDomainPartition().getMeshBody( 0 ).getBaseDiscretization();

    NodeManager const & nodeManager = mesh.getNodeManager();
    EdgeManager const & edgeManager = mesh.getEdgeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // ---------- local stats ----------
    localIndex localSurfaceElements = 0;
    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & sr )
    {
      // Use getNumberOfLocalIndices() to count only owned elements (no double-counting).
      localSurfaceElements += sr.getNumberOfLocalIndices();
    } );

    integer const eulerChar =
      computeEulerCharacteristic( nodeManager, edgeManager, faceManager, elemManager );

    // ---------- global reduction ----------
    localIndex const globalSurfaceElements =
      MpiWrapper::sum( localSurfaceElements, MPI_COMM_GEOS );

    GEOS_LOG_RANK_0( "========================================" );
    GEOS_LOG_RANK_0( "Test (MPI): " << testCaseName
                                    << " [" << xPartitions << "x" << yPartitions << "x" << zPartitions << "]" );
    GEOS_LOG_RANK_0( "  Surface elements (global): " << globalSurfaceElements );
    GEOS_LOG_RANK_0( "  Euler χ:                   " << eulerChar );

    // Run events (VTK output)
    state.run();

    MpiWrapper::barrier( MPI_COMM_GEOS );

    // ---------- assertions on rank 0 only ----------
    if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
    {
      SplitMeshStats stats;
      stats.numSurfaceElements = globalSurfaceElements;
      stats.eulerCharacteristic = eulerChar;

      validateSplitMeshResults( testCaseName, meshFileName, stats,
                                expectedSurfaceElements, expectedEuler,
                                nodeManager, elemManager );
    }

    GEOS_LOG_RANK_0( "========================================" );

  } // end scoped GeosxState

  if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0 )
  {
    std::remove( xmlPath.c_str() );
  }
}

// ---------------------------------------------------------------------------
// MPI test instantiation
//
// Partitions of 4 ranks:  (1,1,4)  (1,4,1)  (4,1,1)  (1,2,2)  (2,1,2)  (2,2,1)
//
// Column layout: test-name, vtm-file, faceBlock-name,
//                expected-surface-elements, expected-Euler-χ
// ---------------------------------------------------------------------------
// clang-format off
INSTANTIATE_TEST_SUITE_P(
  SplitMesh_mpiCases,
  SplitMesh_mpiTest,
  ::testing::Combine(
    ::testing::Values(
      //                    test name              mesh file (vtm)                          faceBlock   surf-elems  χ
      std::make_tuple( "DFN_5_fracs_hex", "DFN_5_fractures_hex_binarized.vtm", "Fault",
                       localIndex( -1 ), integer( 1 ) )
      // NOTE: Set the expected surface element count to -1 as a placeholder until
      //       the actual value is known from a reference run.
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
// clang-format on

int main( int argc, char * argv[] )
{
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
