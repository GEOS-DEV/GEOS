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
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testMixedDimHydrostaticEquilibriumSplitMesh_mpi.cpp
 *
 * MPI (multi-rank) variant of the mixed-dimensional hydrostatic equilibrium test
 * for pre-split (VTM) meshes with fractures defined via faceBlocks (WF2).
 * Test fixture and TEST_P body live in testMixedDimHydrostaticEquilibriumSplitMeshFixture.hpp.
 */

#include "testMixedDimHydrostaticEquilibriumSplitMeshFixture.hpp"

CommandLineOptions g_commandLineOptions;

/**
 * @brief Parallel execution test cases (4 MPI ranks).
 *
 * Mirrors the mesh set and partitioning schemes used in
 * testMixedDimHydrostaticEquilibrium_mpi.cpp (WF1) but uses the
 * corresponding .vtm composite files instead of .vtu files.
 *
 * 14 meshes × 2 partitioning schemes = 28 test cases.
 */
INSTANTIATE_TEST_SUITE_P(
  MixedDimHydrostaticEquilibriumSplitMeshCases,
  MixedDimHydrostaticEquilibriumSplitMeshTest,
  ::testing::Combine(
    ::testing::Values(
      // Flat tet meshes
      "fractured_mesh_tet_DFN_1.vtm",
      "fractured_mesh_tet_DFN_123.vtm",

      // Wavy tet meshes
      "fractured_wavy_mesh_tet_DFN_1.vtm",
      "fractured_wavy_mesh_tet_DFN_123.vtm",

      // Full span hex meshes
      "fractured_full_span_mesh_hex_DFN_1.vtm",
      "fractured_full_span_mesh_hex_DFN_123.vtm",

      // Full span tet meshes
      "fractured_full_span_mesh_tet_DFN_1.vtm",
      "fractured_full_span_mesh_tet_DFN_123.vtm",

      // T-shaped wavy meshes
      "t_shaped_wavy_mesh_hex_DFN_t1t2.vtm",
      "t_shaped_wavy_mesh_tet_DFN_t1t2.vtm",

      // Y-shaped wavy meshes
      "y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtm",
      "y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtm",

      // 5-fracture DFN market meshes
      "DFN_5_fractures_hex_binarized.vtm",
      "DFN_5_fractures_tet_binarized.vtm"
      ),
    ::testing::Values(
      std::make_tuple( 1, 4, 1 ),
      std::make_tuple( 2, 1, 2 )
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
