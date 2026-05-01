an/*
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
 * @file testFEMConsistencySplitMesh.cpp
 *
 * Serial (single-rank) variant of the FEM consistency test for split (pre-fractured) meshes.
 * Test fixture and TEST_P body live in testFEMConsistencySplitMeshFixture.hpp.
 */

#include "testFEMConsistencySplitMeshFixture.hpp"

CommandLineOptions g_commandLineOptions;

/**
 * @brief Serial execution test cases (single rank, partition 1x1x1).
 */
INSTANTIATE_TEST_SUITE_P(
  FractureStressCasesSplitMeshSerial,
  ConsistencySplitMeshTest,
  ::testing::Combine(
    ::testing::Values(
      // 5-fracture DFN split meshes
      "DFN_5_fractures_hex_binarized.vtm",
      "DFN_5_fractures_tet_binarized.vtm"
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
