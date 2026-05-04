
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
 * @file testFEMConsistencySplitMesh_mpi.cpp
 *
 * MPI (multi-rank) variant of the FEM consistency test for split (pre-fractured) meshes.
 * Test fixture and TEST_P body live in testFEMConsistencySplitMeshFixture.hpp.
 */

#include "testFEMConsistencySplitMeshFixture.hpp"

CommandLineOptions g_commandLineOptions;

/**
 * @brief Parallel execution test cases (4 MPI ranks).
 *
 * The parameters are:
 * 1. Mesh file name (std::string): The VTM composite mesh file containing the geometry and fractures.
 * 2. s_xx (real64): Applied stress component XX.
 * 3. s_yy (real64): Applied stress component YY.
 * 4. s_zz (real64): Applied stress component ZZ.
 * 5. Partitioning (tuple<int, int, int>): Number of partitions in x, y, z directions.
 */
INSTANTIATE_TEST_SUITE_P(
  FractureStressCasesSplitMeshPartitioned,
  ConsistencySplitMeshTest,
  ::testing::Combine(
    ::testing::Values(
      // Flat hex split meshes - single and triple fracture only
      "fractured_mesh_hex_DFN_1.vtm",
      "fractured_mesh_hex_DFN_123.vtm",

      // Wavy hex split meshes - single and triple fracture only
      "fractured_wavy_mesh_hex_DFN_1.vtm",
      "fractured_wavy_mesh_hex_DFN_123.vtm",

      // Full span hex split meshes - single and triple fracture only
      "fractured_full_span_mesh_hex_DFN_1.vtm",
      "fractured_full_span_mesh_hex_DFN_123.vtm",

      // Full span tet split meshes - single and triple fracture only
      "fractured_full_span_mesh_tet_DFN_1.vtm",
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
    ::testing::Values( -1.0e6 ),     // s_xx
    ::testing::Values( -0.5e6 ),     // s_yy
    ::testing::Values( -2.0e6 ),     // s_zz
    ::testing::Values(
      std::make_tuple( 2, 2, 1 ),
      std::make_tuple( 4, 1, 1 )
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
