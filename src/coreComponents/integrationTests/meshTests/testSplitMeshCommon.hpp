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

#ifndef TESTSPLITMESHCOMMON_HPP
#define TESTSPLITMESHCOMMON_HPP

#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>               // std::strlen
#include <libgen.h>              // dirname (POSIX)

#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/generators/MeshGeneratorBase.hpp"
#include "common/MpiWrapper.hpp"

using namespace geos;

// Global flag to control debug printing
static constexpr bool ENABLE_SPLIT_MESH_DEBUG_PRINTS = false;

extern CommandLineOptions g_commandLineOptions;
extern std::string g_testBinaryDir;

/**
 * @brief Statistics collected from a split mesh.
 */
struct SplitMeshStats
{
  localIndex numNodes            = 0;
  localIndex numEdges            = 0;
  localIndex numFaces            = 0;
  localIndex numCells            = 0;
  localIndex numSurfaceElements  = 0;  ///< Elements in SurfaceElementRegion(s)
  integer eulerCharacteristic = 0;
};

/**
 * @brief Validate that the loaded split mesh satisfies basic sanity checks.
 *
 * Checks performed:
 *   B1 – Number of surface elements matches the expected count.
 *   B2 – All node coordinates are valid (no NaN).
 *   B3 – Euler characteristic matches expected value.
 *   B4 – No element has duplicate node references (degenerate Jacobian).
 */
void validateSplitMeshResults( std::string const & testCaseName,
                               std::string const & meshFileName,
                               SplitMeshStats const & stats,
                               localIndex const expectedSurfaceElements,
                               integer const expectedEuler,
                               NodeManager const & nodeManager,
                               ElementRegionManager const & elemManager )
{
  // B1: Surface element count
  bool const b1Pass = ( stats.numSurfaceElements == expectedSurfaceElements );
  GEOS_LOG_RANK_0( "Validating B1 " << ( b1Pass ? "(ok)" : "(notok)" )
                                    << ": Surface elements (Expected: " << expectedSurfaceElements
                                    << ", Actual: " << stats.numSurfaceElements << ")" );
  EXPECT_EQ( stats.numSurfaceElements, expectedSurfaceElements )
    << "Test " << testCaseName << ": Surface element count MISMATCH"
    << "\n  Expected: " << expectedSurfaceElements
    << "\n  Actual:   " << stats.numSurfaceElements
    << "\n  Mesh: "     << meshFileName;

  // B2: Node coordinate validity (no NaN)
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions =
    nodeManager.referencePosition();
  bool b2Pass = true;
  for( localIndex ni = 0; ni < nodeManager.size(); ++ni )
  {
    real64 const x = nodePositions[ni][0];
    real64 const y = nodePositions[ni][1];
    real64 const z = nodePositions[ni][2];
    if( std::isnan( x ) || std::isnan( y ) || std::isnan( z ) )
      b2Pass = false;
    EXPECT_FALSE( std::isnan( x ) || std::isnan( y ) || std::isnan( z ) )
      << "Test " << testCaseName << ": Node " << ni << " has invalid (NaN) coordinates"
      << "\n  Position: (" << x << ", " << y << ", " << z << ")"
      << "\n  Total nodes: " << nodeManager.size();
  }
  GEOS_LOG_RANK_0( "Validating B2 " << ( b2Pass ? "(ok)" : "(notok)" )
                                    << ": Node coordinate validity (nodes: " << nodeManager.size() << ")" );

  // B3: Euler characteristic
  bool const b3Pass = ( stats.eulerCharacteristic == expectedEuler );
  GEOS_LOG_RANK_0( "Validating B3 " << ( b3Pass ? "(ok)" : "(notok)" )
                                    << ": Euler χ (Expected: " << expectedEuler
                                    << ", Actual: " << stats.eulerCharacteristic << ")" );
  EXPECT_EQ( stats.eulerCharacteristic, expectedEuler )
    << "Test " << testCaseName << ": Euler characteristic MISMATCH"
    << "\n  Expected χ: " << expectedEuler
    << "\n  Actual χ:   " << stats.eulerCharacteristic
    << "\n  Mesh: "       << meshFileName;

  // B4: No degenerate elements (duplicate node references)
  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap = subRegion.nodeList();
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      localIndex const numNodesPerElem = elemToNodeMap.size( 1 );
      std::set< localIndex > uniqueNodes;
      for( localIndex i = 0; i < numNodesPerElem; ++i )
        uniqueNodes.insert( elemToNodeMap[ei][i] );
      EXPECT_EQ( (localIndex)uniqueNodes.size(), numNodesPerElem )
        << "Test " << testCaseName << ": Element " << ei
        << " in subRegion '" << subRegion.getName() << "' has duplicate node references";
    }
  } );
}

#endif // TESTSPLITMESHCOMMON_HPP
