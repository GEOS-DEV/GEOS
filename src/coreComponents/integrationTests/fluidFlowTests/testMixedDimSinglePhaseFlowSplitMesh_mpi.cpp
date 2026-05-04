you/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testMixedDimSinglePhaseFlowSplitMesh_mpi.cpp
 *
 * MPI (multi-rank) variant of the mixed-dimensional single-phase flow test
 * for pre-split (VTM) meshes with fractures defined via faceBlocks.
 * Test fixture and TEST_P body live in testMixedDimSinglePhaseFlowSplitMeshFixture.hpp.
 */

#include "testMixedDimSinglePhaseFlowSplitMeshFixture.hpp"

CommandLineOptions g_commandLineOptions;

/**
 * @brief Parallel execution test cases (4 MPI ranks).
 *
 * The parameters are:
 * 1. Mesh file name (std::string): The VTM composite mesh file.
 * 2. bool: Whether to run the flow solver (true) or just check initialization (false).
 * 3. Partitioning (tuple<int, int, int>): Number of partitions in x, y, z directions.
 */
INSTANTIATE_TEST_SUITE_P(
  MixedDimSplitMeshPartitionedFlowCases,
  MixedDimSinglePhaseFlowSplitMeshTest,
  ::testing::Combine(
    ::testing::Values(
      // Flat hex split meshes
      "fractured_mesh_hex_DFN_123.vtm",

      // Wavy hex split meshes
      "fractured_wavy_mesh_hex_DFN_123.vtm",

      // Full span hex split meshes
      "fractured_full_span_mesh_hex_DFN_123.vtm",

      // Full span tet split meshes
      "fractured_full_span_mesh_tet_DFN_123.vtm",

      // T-shaped wavy split meshes
      "t_shaped_wavy_mesh_hex_DFN_t1t2.vtm",
      "t_shaped_wavy_mesh_tet_DFN_t1t2.vtm",

      // Y-shaped wavy split meshes
      "y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtm",
      "y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtm",

      // 5-fracture DFN split meshes
      "DFN_5_fractures_hex_binarized.vtm",
      "DFN_5_fractures_tet_binarized.vtm"
      ),
    ::testing::Bool(),
    ::testing::Values(
      std::make_tuple( 1, 1, 4 ),
      std::make_tuple( 1, 2, 2 )
      )
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
