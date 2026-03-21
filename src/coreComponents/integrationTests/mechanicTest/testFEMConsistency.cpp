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
 * @file testFEMConsistency.cpp
 *
 * Serial (single-rank) variant of the FEM consistency test.
 * Test fixture and TEST_P body live in testFEMConsistencyFixture.hpp.
 */

#include "testFEMConsistencyFixture.hpp"

CommandLineOptions g_commandLineOptions;

/**
 * @brief Serial execution test cases (single rank, partition 1x1x1).
 */
INSTANTIATE_TEST_SUITE_P(
  FractureStressCasesSerial,
  ConsistencyTest,
  ::testing::Combine(
    ::testing::Values(
      // Flat hex meshes - single and triple fracture only
      "fractured_mesh_hex_DFN_1.vtu"//,
//      "fractured_mesh_hex_DFN_123.vtu"

      // Wavy hex meshes - single and triple fracture only
      // "fractured_wavy_mesh_hex_DFN_1.vtu",
      // "fractured_wavy_mesh_hex_DFN_123.vtu",

      // // Full span hex meshes - single and triple fracture only
      // "fractured_full_span_mesh_hex_DFN_1.vtu",
      // "fractured_full_span_mesh_hex_DFN_123.vtu",

      // // Full span tet meshes - single and triple fracture only
      // "fractured_full_span_mesh_tet_DFN_1.vtu",
      // "fractured_full_span_mesh_tet_DFN_123.vtu",

      // // T-shaped wavy meshes
      // "t_shaped_wavy_mesh_hex_DFN_t1t2.vtu",
      // "t_shaped_wavy_mesh_tet_DFN_t1t2.vtu",

      // // Y-shaped wavy meshes
      // "y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtu",
      // "y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtu",

      // // 5-fracture DFN market meshes
      // "DFN_5_fractures_hex_binarized.vtu",
      // "DFN_5_fractures_tet_binarized.vtu"
      ),
    ::testing::Values( -1.0e6 ),     // s_xx
    ::testing::Values( -0.5e6 ),     // s_yy
    ::testing::Values( -2.0e6 ),     // s_zz
    ::testing::Values(
      std::make_tuple( 1, 1, 1 )
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
