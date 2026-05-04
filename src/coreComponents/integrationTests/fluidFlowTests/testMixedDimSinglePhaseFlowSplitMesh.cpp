/*
 d * ------------------------------------------------------------------------------------------------------------
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
 * @file testMixedDimSinglePhaseFlowSplitMesh.cpp
 *
 * Serial (single-rank) variant of the mixed-dimensional single-phase flow test
 * for pre-split (VTM) meshes with fractures defined via faceBlocks.
 * Test fixture and TEST_P body live in testMixedDimSinglePhaseFlowSplitMeshFixture.hpp.
 */

#include "testMixedDimSinglePhaseFlowSplitMeshFixture.hpp"

CommandLineOptions g_commandLineOptions;

/**
 * @brief Serial execution test cases (single rank, partition 1x1x1).
 */
INSTANTIATE_TEST_SUITE_P(
  MixedDimSplitMeshSerialFlowCases,
  MixedDimSinglePhaseFlowSplitMeshTest,
  ::testing::Combine(
    ::testing::Values(
      "fractured_mesh_hex_DFN_123.vtm"
      ),
    ::testing::Bool(),
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
