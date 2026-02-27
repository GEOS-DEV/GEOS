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

#ifndef TESTSURFACEGENERATORCOMMON_HPP
#define TESTSURFACEGENERATORCOMMON_HPP

#include <gtest/gtest.h>
#include <fstream>
#include <tuple>
#include <sstream>
#include <set>
#include <map>
#include <queue>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>               // std::strlen
#include <libgen.h>              // dirname / basename (POSIX)

#include "mainInterface/GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "codingUtilities/Utilities.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
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

constexpr real64 COORDINATE_TOLERANCE = 1.0e-10;
constexpr real64 JACOBIAN_TOLERANCE = 1.0e-12;

// Global flag to control debug printing (set to false to disable verbose output)
static constexpr bool ENABLE_DEBUG_PRINTS = false;

extern CommandLineOptions g_commandLineOptions;
extern std::string g_testBinaryDir;

/**
 * @brief Represents a unique edge by its sorted node IDs
 */
struct Edge
{
  localIndex nodeA;
  localIndex nodeB;

  Edge( localIndex a, localIndex b )
    : nodeA( std::min( a, b )), nodeB( std::max( a, b ))
  {}

  bool operator<( Edge const & other ) const
  {
    return (nodeA < other.nodeA) || (nodeA == other.nodeA && nodeB < other.nodeB);
  }

  bool operator==( Edge const & other ) const
  {
    return nodeA == other.nodeA && nodeB == other.nodeB;
  }
};

/**
 * @brief Represents a unique facet by its sorted node IDs
 */
struct Facet
{
  std::vector< localIndex > nodes;

  Facet( std::initializer_list< localIndex > nodeList )
  {
    nodes.assign( nodeList );
    std::sort( nodes.begin(), nodes.end());
  }

  template< typename Iterator >
  Facet( Iterator begin, Iterator end )
  {
    nodes.assign( begin, end );
    std::sort( nodes.begin(), nodes.end());
  }

  Facet( arraySlice1d< localIndex const > nodeList )
  {
    nodes.assign( nodeList.begin(), nodeList.end());
    std::sort( nodes.begin(), nodes.end());
  }

  bool operator<( Facet const & other ) const
  {
    return nodes < other.nodes;
  }

  bool operator==( Facet const & other ) const
  {
    return nodes == other.nodes;
  }
};

/**
 * @brief Expected topology duplication statistics.
 *
 * Internal entities (not on the [0,1]³ boundary) are duplicated; boundary entities are not.
 */
struct ExpectedDuplication
{
  // Bulk mesh statistics (all cell elements)
  localIndex numNodesInBulk = 0;
  localIndex numEdgesInBulk = 0;
  localIndex numFacesInBulk = 0;
  localIndex numCellsInBulk = 0;

  // Fracture statistics
  localIndex numFractureNodes = 0;         // Total nodes in fracture nodesets
  localIndex numFractureEdges = 0;         // Total edges in fracture
  localIndex numFractureFaces = 0;         // Total faces in fracture
  localIndex numFractureElements = 0;      // Expected fracture elements (2D faces to create)

  // Duplication statistics (only internal entities)
  localIndex numNodesToDuplicate = 0;      // Internal fracture nodes (not on boundary)
  localIndex numEdgesToDuplicate = 0;      // Internal fracture edges (not on boundary)
  localIndex numFacesToDuplicate = 0;      // Fracture faces to duplicate

  localIndex totalDuplicatedNodes = 0;     // New nodes created
  localIndex totalDuplicatedEdges = 0;     // New edges created
  localIndex totalDuplicatedFaces = 0;     // New faces created

  std::set< localIndex > internalNodes;    // Nodes to duplicate (interior only)
  std::set< localIndex > boundaryNodes;    // Nodes on domain boundary (not duplicated)
  std::set< Edge > internalEdges;          // Edges to duplicate (interior only)
  std::set< Edge > boundaryEdges;          // Edges on domain boundary (not duplicated)

  localIndex numBoundaryNodesOnDomain = 0; // Count of fracture boundary nodes that touch domain boundary
};

/**
 * @brief Helper structure to store mesh topology statistics
 */
struct TopologyStats
{
  localIndex numNodes;
  localIndex numEdges;
  localIndex numFaces;
  localIndex numCells;
  localIndex numDuplicatedNodes;
  localIndex numFractureElements;
  integer eulerCharacteristic;
  integer numBodies;  // Number of separate domain fragments
};

/**
 * @brief Compute Euler-Poincaré characteristic χ = V - E + F - C from CellElementSubRegions only.
 *
 * χ = 1 for a single connected solid. Excludes FaceElementSubRegion (fracture surfaces).
 */
integer computeEulerCharacteristicTestHelper( NodeManager const & nodeManager,
                                              FaceManager const & faceManager,
                                              ElementRegionManager const & elemManager )
{
  (void) nodeManager;
  (void) faceManager;
  std::set< localIndex > uniqueNodes;
  std::set< Edge > uniqueEdges;
  std::set< Facet > uniqueFaces;
  localIndex numCells = 0;

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap = subRegion.nodeList();

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      ++numCells;
      localIndex const numNodesPerElem = elemToNodeMap.size( 1 );

      // Collect all nodes
      for( localIndex ni = 0; ni < numNodesPerElem; ++ni )
      {
        uniqueNodes.insert( elemToNodeMap[ei][ni] );
      }

      // Extract edges and faces based on element type
      if( numNodesPerElem == 4 )  // Tetrahedron
      {
        // Edges (6)
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][1] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][2] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][3] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][2] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][3] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][3] ) );

        // Faces (4)
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][2]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][3]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][2], elemToNodeMap[ei][3]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][1], elemToNodeMap[ei][2], elemToNodeMap[ei][3]} ) );
      }
      else if( numNodesPerElem == 8 )  // Hexahedron
      {
        // NOTE: GEOS uses ordering { 0, 1, 3, 2, 4, 5, 7, 6 } (see SiloFile.cpp:1338)
        // Bottom face: 0->1->3->2, Top face: 4->5->7->6

        // Edges (12): 4 bottom + 4 top + 4 vertical
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][1] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][3] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][3], elemToNodeMap[ei][2] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][0] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][4], elemToNodeMap[ei][5] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][5], elemToNodeMap[ei][7] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][7], elemToNodeMap[ei][6] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][6], elemToNodeMap[ei][4] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][4] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][5] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][3], elemToNodeMap[ei][7] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][6] ) );

        // Faces (6): bottom, top, and 4 sides
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][3], elemToNodeMap[ei][2]} ) );  // bottom
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][4], elemToNodeMap[ei][5], elemToNodeMap[ei][7], elemToNodeMap[ei][6]} ) );  // top
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][5], elemToNodeMap[ei][4]} ) );  // front
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][2], elemToNodeMap[ei][3], elemToNodeMap[ei][7], elemToNodeMap[ei][6]} ) );  // back
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][1], elemToNodeMap[ei][3], elemToNodeMap[ei][7], elemToNodeMap[ei][5]} ) );  // right
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][2], elemToNodeMap[ei][6], elemToNodeMap[ei][4]} ) );  // left
      }
      else if( numNodesPerElem == 6 )  // Prism
      {
        // Edges (9)
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][1] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][2] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][0] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][3], elemToNodeMap[ei][4] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][4], elemToNodeMap[ei][5] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][5], elemToNodeMap[ei][3] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][3] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][4] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][5] ) );

        // Faces (5)
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][2]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][3], elemToNodeMap[ei][4], elemToNodeMap[ei][5]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][4], elemToNodeMap[ei][3]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][1], elemToNodeMap[ei][2], elemToNodeMap[ei][5], elemToNodeMap[ei][4]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][2], elemToNodeMap[ei][0], elemToNodeMap[ei][3], elemToNodeMap[ei][5]} ) );
      }
      else if( numNodesPerElem == 5 )  // Pyramid
      {
        // Edges (8)
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][1] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][2] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][3] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][3], elemToNodeMap[ei][0] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][4] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][4] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][4] ) );
        uniqueEdges.insert( Edge( elemToNodeMap[ei][3], elemToNodeMap[ei][4] ) );
        // Faces (5)
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][2], elemToNodeMap[ei][3]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][4]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][1], elemToNodeMap[ei][2], elemToNodeMap[ei][4]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][2], elemToNodeMap[ei][3], elemToNodeMap[ei][4]} ) );
        uniqueFaces.insert( Facet( {elemToNodeMap[ei][3], elemToNodeMap[ei][0], elemToNodeMap[ei][4]} ) );
      }
    }
  } );

  localIndex const V = uniqueNodes.size();
  localIndex const E = uniqueEdges.size();
  localIndex const F = uniqueFaces.size();
  localIndex const C = numCells;

  // Compute Euler-Poincaré characteristic: χ = V - E + F - C
  // For a connected 3D solid: χ = 1
  integer const eulerCharacteristic = V - E + F - C;

  // Print diagnostic information
  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Euler characteristic computation:" );
    GEOS_LOG_RANK_0( "  V (nodes):  " << V );
    GEOS_LOG_RANK_0( "  E (edges):  " << E );
    GEOS_LOG_RANK_0( "  F (facets): " << F );
    GEOS_LOG_RANK_0( "  C (cells):  " << C );
    GEOS_LOG_RANK_0( "  χ = V - E + F - C = " << V << " - " << E << " + " << F << " - " << C << " = " << eulerCharacteristic );
  }

  return eulerCharacteristic;
}

/// @brief Assert all duplicated nodes have zero spatial distance from their parent.
void checkCoordinateCoincidence( NodeManager const & nodeManager,
                                 std::map< localIndex, localIndex > const & parentToChild )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions =
    nodeManager.referencePosition();

  for( auto const & [parentIdx, childIdx] : parentToChild )
  {
    real64 dx = nodePositions[childIdx][0] - nodePositions[parentIdx][0];
    real64 dy = nodePositions[childIdx][1] - nodePositions[parentIdx][1];
    real64 dz = nodePositions[childIdx][2] - nodePositions[parentIdx][2];
    real64 dist = std::sqrt( dx*dx + dy*dy + dz*dz );

    EXPECT_NEAR( dist, 0.0, COORDINATE_TOLERANCE )
      << "Duplicated node " << childIdx << " (parent " << parentIdx << ") has non-zero spatial distance: " << dist
      << "\n  Parent coordinates: (" << nodePositions[parentIdx][0] << ", " << nodePositions[parentIdx][1] << ", " << nodePositions[parentIdx][2] << ")"
      << "\n  Child coordinates:  (" << nodePositions[childIdx][0] << ", " << nodePositions[childIdx][1] << ", " << nodePositions[childIdx][2] << ")"
      << "\n  Expected tolerance: " << COORDINATE_TOLERANCE;
  }
}

/// @brief Assert no element has duplicate node references (degenerate Jacobian check).
void checkElementJacobians( ElementRegionManager const & elemManager )
{
  localIndex totalElements = 0;
  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    totalElements += subRegion.size();
  } );

  GEOS_LOG_RANK_0( "Checking element Jacobians for " << totalElements << " elements" );

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    // Note: Detailed Jacobian computation would require access to shape functions
    // Here we perform a simplified check that elements still reference valid nodes
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap =
      subRegion.nodeList();

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      localIndex const numNodesPerElem = elemToNodeMap.size( 1 );
      std::set< localIndex > uniqueNodes;
      for( localIndex i = 0; i < numNodesPerElem; ++i )
      {
        uniqueNodes.insert( elemToNodeMap[ei][i] );
      }

      // Element should have the correct number of unique nodes
      EXPECT_EQ( uniqueNodes.size(), numNodesPerElem )
        << "Element " << ei << " in subRegion '" << subRegion.getName()
        << "' has duplicate node references"
        << "\n  Expected unique nodes: " << numNodesPerElem
        << "\n  Actual unique nodes: " << uniqueNodes.size()
        << "\n  This indicates nodes on both sides of the fracture are incorrectly connected to the same element";
    }
  } );
}

/**
 * @brief Run all assertions for a single test case (node duplication, fracture elements,
 *        coordinate validity, and Euler characteristic).
 */
void validateSurfaceGeneratorResults( std::string const & testCaseName,
                                      std::string const & meshFileName,
                                      std::string const & nodeSetNames,
                                      TopologyStats const & statsBefore,
                                      TopologyStats const & statsAfter,
                                      ExpectedDuplication const & expected,
                                      integer const expectedEulerBefore,
                                      integer const eulerCharBeforeSplit,
                                      integer const expectedEulerAfter,
                                      integer const eulerCharAfterSplit,
                                      NodeManager const & nodeManager )
{
  // A1: Validate node duplication matches prediction
  // Note: totalDuplicatedNodes = total NEW nodes created (not counting original nodes)
  // For a node at N-fracture intersection: N new nodes are created
  bool const a1Pass = ( statsAfter.numDuplicatedNodes == expected.totalDuplicatedNodes );
  GEOS_LOG_RANK_0( "Validating A1 " << ( a1Pass ? "(ok)" : "(notok)" )
                                    << ": Node duplication prediction (Expected: " << expected.totalDuplicatedNodes
                                    << ", Actual: " << statsAfter.numDuplicatedNodes << ")" );
  EXPECT_EQ( statsAfter.numDuplicatedNodes, expected.totalDuplicatedNodes )
    << "Test " << testCaseName << ": Node duplication prediction MISMATCH"
    << "\n  Expected new nodes: " << expected.totalDuplicatedNodes
    << "\n  Actual new nodes:   " << statsAfter.numDuplicatedNodes
    << "\n  Nodes before split: " << statsBefore.numNodes
    << "\n  Nodes after split:  " << statsAfter.numNodes
    << "\n  Fracture node sets: " << nodeSetNames;

  // A2: Verify fracture elements were created and match the predicted count
  bool const a2Pass = ( statsAfter.numFractureElements == expected.numFractureElements );
  GEOS_LOG_RANK_0( "Validating A2 " << ( a2Pass ? "(ok)" : "(notok)" )
                                    << ": Fracture elements created (Expected: " << expected.numFractureElements
                                    << ", Actual: " << statsAfter.numFractureElements << ")" );
  EXPECT_EQ( statsAfter.numFractureElements, expected.numFractureElements )
    << "Test " << testCaseName << ": Fracture element count MISMATCH"
    << "\n  Expected fracture elements: " << expected.numFractureElements
    << "\n  Actual fracture elements:   " << statsAfter.numFractureElements
    << "\n  Mesh file: " << meshFileName
    << "\n  Fracture node sets: " << nodeSetNames
    << "\n  Nodes split: " << statsAfter.numDuplicatedNodes
    << "\n  POSSIBLE CAUSES:"
    << "\n    - Node sets are empty or not found in mesh"
    << "\n    - Faces are not marked as separable (ruptureState/isFaceSeparable)"
    << "\n    - SurfaceGenerator did not run or failed silently"
    << "\n    - Fracture face count prediction is wrong (check preprocessFractureTopology)";

  // A3: Validate all node coordinates are valid (not NaN)
  bool a3Pass = true;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions =
    nodeManager.referencePosition();
  for( localIndex ni = 0; ni < nodeManager.size(); ++ni )
  {
    real64 const x = nodePositions[ni][0];
    real64 const y = nodePositions[ni][1];
    real64 const z = nodePositions[ni][2];
    if( std::isnan( x ) || std::isnan( y ) || std::isnan( z ) )
      a3Pass = false;
    EXPECT_FALSE( std::isnan( x ) || std::isnan( y ) || std::isnan( z ) )
      << "Test " << testCaseName << ": Node " << ni << " has invalid (NaN) coordinates"
      << "\n  Position: (" << x << ", " << y << ", " << z << ")"
      << "\n  Total nodes in mesh: " << nodeManager.size()
      << "\n  This indicates memory corruption or uninitialized node positions";
  }
  GEOS_LOG_RANK_0( "Validating A3 " << ( a3Pass ? "(ok)" : "(notok)" )
                                    << ": Node coordinates validity (Total nodes: " << nodeManager.size() << ")" );

  // A4: Validate expected Euler characteristic before split
  bool const a4Pass = ( eulerCharBeforeSplit == expectedEulerBefore );
  GEOS_LOG_RANK_0( "Validating A4 " << ( a4Pass ? "(ok)" : "(notok)" )
                                    << ": Euler χ before split (Expected: " << expectedEulerBefore
                                    << ", Actual: " << eulerCharBeforeSplit << ")" );
  EXPECT_EQ( eulerCharBeforeSplit, expectedEulerBefore )
    << "Test " << testCaseName << ": Euler characteristic MISMATCH before split"
    << "\n  Expected χ: " << expectedEulerBefore
    << "\n  Actual χ:   " << eulerCharBeforeSplit
    << "\n  Mesh file:  " << meshFileName
    << "\n  NOTE: For a conformal mesh, expected χ = 1 (single connected solid)";

  // // A5: Validate expected Euler characteristic after split
  // bool const a5Pass = ( eulerCharAfterSplit == expectedEulerAfter );
  // GEOS_LOG_RANK_0( "Validating A5 " << ( a5Pass ? "(ok)" : "(notok)" )
  //                                   << ": Euler χ after split (Expected: " << expectedEulerAfter
  //                                   << ", Actual: " << eulerCharAfterSplit << ")" );
  EXPECT_EQ( eulerCharAfterSplit, expectedEulerAfter )
    << "Test " << testCaseName << ": Euler characteristic MISMATCH after split"
    << "\n  Expected χ: " << expectedEulerAfter
    << "\n  Actual χ:   " << eulerCharAfterSplit
    << "\n  Mesh file:  " << meshFileName
    << "\n  Nodes duplicated: " << statsAfter.numDuplicatedNodes
    << "\n  Fracture elements: " << statsAfter.numFractureElements;
}

/// @brief Parse "{ f1_node_set, f2_node_set, ... }" into individual set name strings.
std::vector< std::string > parseNodeSetNames( std::string const & nodeSetNames )
{
  std::vector< std::string > result;

  // Remove surrounding braces and whitespace
  size_t start = nodeSetNames.find( '{' );
  size_t end = nodeSetNames.find( '}' );

  if( start == std::string::npos || end == std::string::npos )
  {
    return result;
  }

  std::string content = nodeSetNames.substr( start + 1, end - start - 1 );

  // Split by comma
  std::istringstream ss( content );
  std::string token;

  while( std::getline( ss, token, ',' ) )
  {
    // Trim whitespace
    size_t first = token.find_first_not_of( " \t\n\r" );
    size_t last = token.find_last_not_of( " \t\n\r" );

    if( first != std::string::npos && last != std::string::npos )
    {
      result.push_back( token.substr( first, last - first + 1 ) );
    }
  }

  return result;
}

/**
 * @brief Compute expected node/edge/face duplication counts before the mesh split.
 *
 * Domain assumed to be [0,1]³. Entities on the boundary are not duplicated; interior
 * and intersection entities are. Returns an @c ExpectedDuplication for validation.
 */
ExpectedDuplication preprocessFractureTopology( MeshLevel & mesh,
                                                std::vector< std::string > const & fractureNodeSets )
{
  ExpectedDuplication expected;

  NodeManager & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions =
    nodeManager.referencePosition();

  // Ghost-rank views: entities with ghostRank > -1 are owned by another rank.
  // We must skip them so that MpiWrapper::sum across ranks gives correct global totals
  // without double-counting entities that appear on multiple ranks as ghosts.
  arrayView1d< integer const > const nodeGhostRank = nodeManager.ghostRank();
  arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();

  // Domain boundary tolerance for [0,1]³
  constexpr real64 boundaryTol = 1.0e-9;

  // Helper to check if a node is on the DOMAIN BOUNDARY (not to be duplicated)
  // A node is on the boundary if ANY coordinate is at 0 or 1
  auto isNodeOnBoundary = [&]( localIndex nodeIdx ) -> bool
  {
    real64 const x = nodePositions[nodeIdx][0];
    real64 const y = nodePositions[nodeIdx][1];
    real64 const z = nodePositions[nodeIdx][2];

    bool const onXBoundary = (std::abs( x ) < boundaryTol) || (std::abs( x - 1.0 ) < boundaryTol);
    bool const onYBoundary = (std::abs( y ) < boundaryTol) || (std::abs( y - 1.0 ) < boundaryTol);
    bool const onZBoundary = (std::abs( z ) < boundaryTol) || (std::abs( z - 1.0 ) < boundaryTol);

    return onXBoundary || onYBoundary || onZBoundary;
  };

  // Step 1: Compute bulk mesh statistics
  std::set< localIndex > bulkNodes;
  std::set< Edge > bulkEdges;
  std::set< Facet > bulkFaces;

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap = subRegion.nodeList();

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      ++expected.numCellsInBulk;
      localIndex const numNodesPerElem = elemToNodeMap.size( 1 );

      for( localIndex ni = 0; ni < numNodesPerElem; ++ni )
      {
        bulkNodes.insert( elemToNodeMap[ei][ni] );
      }

      if( numNodesPerElem == 8 )  // Hexahedron
      {
        // NOTE: GEOS uses ordering { 0, 1, 3, 2, 4, 5, 7, 6 } (see SiloFile.cpp:1338)
        // Bottom face: 0->1->3->2, Top face: 4->5->7->6

        // Create all 12 edges for this hexahedron
        Edge edges[12] = {
          Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][1] ),
          Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][3] ),
          Edge( elemToNodeMap[ei][3], elemToNodeMap[ei][2] ),
          Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][0] ),
          Edge( elemToNodeMap[ei][4], elemToNodeMap[ei][5] ),
          Edge( elemToNodeMap[ei][5], elemToNodeMap[ei][7] ),
          Edge( elemToNodeMap[ei][7], elemToNodeMap[ei][6] ),
          Edge( elemToNodeMap[ei][6], elemToNodeMap[ei][4] ),
          Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][4] ),
          Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][5] ),
          Edge( elemToNodeMap[ei][3], elemToNodeMap[ei][7] ),
          Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][6] )
        };

        // Insert each edge
        for( int i = 0; i < 12; ++i )
        {
          bulkEdges.insert( edges[i] );
        }

        // Create all 6 faces for this hexahedron
        Facet faces[6] = {
          Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][3], elemToNodeMap[ei][2]} ),  // bottom
          Facet( {elemToNodeMap[ei][4], elemToNodeMap[ei][5], elemToNodeMap[ei][7], elemToNodeMap[ei][6]} ),  // top
          Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][5], elemToNodeMap[ei][4]} ),  // front
          Facet( {elemToNodeMap[ei][2], elemToNodeMap[ei][3], elemToNodeMap[ei][7], elemToNodeMap[ei][6]} ),  // back
          Facet( {elemToNodeMap[ei][1], elemToNodeMap[ei][3], elemToNodeMap[ei][7], elemToNodeMap[ei][5]} ),  // right
          Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][2], elemToNodeMap[ei][6], elemToNodeMap[ei][4]} )   // left
        };

        // Insert each face
        for( int i = 0; i < 6; ++i )
        {
          bulkFaces.insert( faces[i] );
        }
      }
      else if( numNodesPerElem == 4 )  // Tetrahedron
      {
        // Edges (6)
        bulkEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][1] ) );
        bulkEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][2] ) );
        bulkEdges.insert( Edge( elemToNodeMap[ei][0], elemToNodeMap[ei][3] ) );
        bulkEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][2] ) );
        bulkEdges.insert( Edge( elemToNodeMap[ei][1], elemToNodeMap[ei][3] ) );
        bulkEdges.insert( Edge( elemToNodeMap[ei][2], elemToNodeMap[ei][3] ) );

        // Faces (4)
        bulkFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][2]} ) );
        bulkFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][1], elemToNodeMap[ei][3]} ) );
        bulkFaces.insert( Facet( {elemToNodeMap[ei][0], elemToNodeMap[ei][2], elemToNodeMap[ei][3]} ) );
        bulkFaces.insert( Facet( {elemToNodeMap[ei][1], elemToNodeMap[ei][2], elemToNodeMap[ei][3]} ) );
      }
    }
  } );

  expected.numNodesInBulk = bulkNodes.size();
  expected.numEdgesInBulk = bulkEdges.size();
  expected.numFacesInBulk = bulkFaces.size();

  // Step 2: Identify fracture entities (track per-fracture for intersections)
  std::set< localIndex > allFractureNodes;
  std::map< localIndex, std::set< integer > > nodeToFractureIds;  // Track which fractures each node belongs to
  std::set< Facet > fractureFaces;
  std::map< Facet, integer > facetToFractureId;  // Track which fracture each face belongs to

  integer fractureId = 0;
  for( std::string const & setName : fractureNodeSets )
  {
    if( !nodeManager.sets().hasWrapper( setName ) )
    {
      ++fractureId;
      continue;
    }

    std::set< localIndex > nodesInThisFracture;
    SortedArrayView< localIndex const > const nodeSet = nodeManager.getSet( setName ).toViewConst();
    for( localIndex const & nodeIdx : nodeSet )
    {
      // Skip ghost nodes: they are owned by another rank and would be double-counted
      // when MpiWrapper::sum accumulates predictions across all ranks.
      if( nodeGhostRank[nodeIdx] > -1 )
        continue;

      allFractureNodes.insert( nodeIdx );
      nodesInThisFracture.insert( nodeIdx );
      nodeToFractureIds[nodeIdx].insert( fractureId );
    }

    // Identify fracture faces for this specific fracture.
    // Skip ghost faces for the same MPI double-counting reason.
    for( localIndex fi = 0; fi < faceManager.size(); ++fi )
    {
      if( faceGhostRank[fi] > -1 )
        continue;

      localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( fi );
      bool allNodesInThisFracture = true;

      for( localIndex ni = 0; ni < numNodesInFace; ++ni )
      {
        if( nodesInThisFracture.find( faceToNodeMap[fi][ni] ) == nodesInThisFracture.end() )
        {
          allNodesInThisFracture = false;
          break;
        }
      }

      if( allNodesInThisFracture && numNodesInFace >= 3 )
      {
        std::vector< localIndex > faceNodes( numNodesInFace );
        for( localIndex ni = 0; ni < numNodesInFace; ++ni )
        {
          faceNodes[ni] = faceToNodeMap[fi][ni];
        }
        Facet facet( faceNodes.begin(), faceNodes.end() );
        fractureFaces.insert( facet );
        facetToFractureId[facet] = fractureId;
      }
    }

    ++fractureId;
  }

  expected.numFractureNodes = allFractureNodes.size();
  expected.numFractureFaces = fractureFaces.size();
  expected.numFractureElements = fractureFaces.size();

  // Step 3: Extract fracture edges from faces
  std::map< Edge, localIndex > edgeToFaceCount;
  std::set< Edge > fractureEdges;

  std::set< Facet > seenCombined;
  std::map< integer, std::set< Facet > > seenPerFracture;
  std::map< integer, std::map< Edge, localIndex > > perFractureEdgeCount;

  for( localIndex fi = 0; fi < faceManager.size(); ++fi )
  {
    // Skip ghost faces to avoid double-counting under MPI.
    if( faceGhostRank[fi] > -1 )
      continue;

    localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( fi );
    bool allNodesInFracture = true;

    for( localIndex ni = 0; ni < numNodesInFace; ++ni )
    {
      if( allFractureNodes.find( faceToNodeMap[fi][ni] ) == allFractureNodes.end() )
      {
        allNodesInFracture = false;
        break;
      }
    }

    if( allNodesInFracture && numNodesInFace >= 3 )
    {
      std::vector< localIndex > faceNodes( numNodesInFace );
      for( localIndex ni = 0; ni < numNodesInFace; ++ni )
        faceNodes[ni] = faceToNodeMap[fi][ni];
      Facet facet( faceNodes.begin(), faceNodes.end() );

      if( seenCombined.find( facet ) == seenCombined.end() )
      {
        seenCombined.insert( facet );
        for( localIndex ni = 0; ni < numNodesInFace; ++ni )
        {
          localIndex const nodeA = faceToNodeMap[fi][ni];
          localIndex const nodeB = faceToNodeMap[fi][(ni + 1) % numNodesInFace];
          Edge edge( nodeA, nodeB );
          fractureEdges.insert( edge );
          edgeToFaceCount[edge]++;
        }
      }

      auto const fracIt = facetToFractureId.find( facet );
      if( fracIt != facetToFractureId.end() )
      {
        integer const faceFractureId = fracIt->second;
        auto & seenForFrac = seenPerFracture[faceFractureId];
        if( seenForFrac.find( facet ) == seenForFrac.end() )
        {
          seenForFrac.insert( facet );
          for( localIndex ni = 0; ni < numNodesInFace; ++ni )
          {
            localIndex const nodeA = faceToNodeMap[fi][ni];
            localIndex const nodeB = faceToNodeMap[fi][(ni + 1) % numNodesInFace];
            Edge edge( nodeA, nodeB );
            perFractureEdgeCount[faceFractureId][edge]++;
          }
        }
      }
    }
  }

  std::map< integer, std::set< localIndex > > perFractureBoundaryNodes;
  for( auto const & [fid, edgeCount] : perFractureEdgeCount )
  {
    for( auto const & [edge, count] : edgeCount )
    {
      if( count == 1 )
      {
        perFractureBoundaryNodes[fid].insert( edge.nodeA );
        perFractureBoundaryNodes[fid].insert( edge.nodeB );
      }
    }
  }

  std::map< localIndex, localIndex > nodeCombBndEdgeCount;
  for( auto const & [edge, count] : edgeToFaceCount )
  {
    if( count == 1 )
    {
      nodeCombBndEdgeCount[edge.nodeA]++;
      nodeCombBndEdgeCount[edge.nodeB]++;
    }
  }

  expected.numFractureEdges = fractureEdges.size();

  // Step 4: Identify fracture boundary nodes
  std::set< localIndex > fractureBoundaryNodes;
  localIndex numBoundaryEdges = 0;
  localIndex numBoundaryEdgesOnDomain = 0;

  for( auto const & [edge, count] : edgeToFaceCount )
  {
    if( count == 1 )  // Edge is on fracture boundary
    {
      ++numBoundaryEdges;
      fractureBoundaryNodes.insert( edge.nodeA );
      fractureBoundaryNodes.insert( edge.nodeB );

      if( isNodeOnBoundary( edge.nodeA ) && isNodeOnBoundary( edge.nodeB ) )
      {
        ++numBoundaryEdgesOnDomain;
      }
    }
  }

  bool const isFullSpan = (numBoundaryEdges > 0) && (numBoundaryEdges == numBoundaryEdgesOnDomain);

  if( isFullSpan )
  {
    fractureBoundaryNodes.clear();
  }

  // Step 5: Apply duplication rules
  //
  // For boundary-cutting (isFullSpan=True) meshes the combined-boundary heuristic is correct:
  // fractureBoundaryNodes was already cleared, so every fracture node is "internal".
  //
  // For non-boundary-cutting (isFullSpan=False) meshes the combined-boundary approach is WRONG
  // when several fractures overlap: a node can sit on the perimeter of the *union* of all
  // fracture faces (landing in fractureBoundaryNodes) while still being interior to one or
  // more individual fractures.  GEOS splits such a node independently for each fracture it is
  // interior to, so we must use per-fracture boundary membership here.
  localIndex numBoundaryNodesTouchingDomain = 0;

  for( localIndex nodeIdx : allFractureNodes )
  {
    bool const isOnDomainBoundary = isNodeOnBoundary( nodeIdx );
    bool shouldDuplicate = false;

    if( isFullSpan )
    {
      // Boundary-cutting path: fractureBoundaryNodes is empty so every fracture node is
      // treated as "internal".  We still need the domain-boundary bookkeeping.
      bool const isOnFractureBoundary = (fractureBoundaryNodes.find( nodeIdx ) != fractureBoundaryNodes.end());
      if( !isOnFractureBoundary )
      {
        shouldDuplicate = true;
      }
      else if( isOnDomainBoundary )
      {
        shouldDuplicate = true;
        ++numBoundaryNodesTouchingDomain;
      }
    }
    else
    {
      // Non-boundary-cutting path: a node is duplicated if it is interior to at least one
      // individual fracture (i.e. NOT on the per-fracture rim of every fracture it belongs to).
      for( integer fid : nodeToFractureIds[nodeIdx] )
      {
        auto const it = perFractureBoundaryNodes.find( fid );
        bool const isOnThisFracBnd = (it != perFractureBoundaryNodes.end())
                                     && (it->second.count( nodeIdx ) > 0);
        if( !isOnThisFracBnd )
        {
          shouldDuplicate = true;
          break;
        }
      }
    }

    if( shouldDuplicate )
    {
      expected.internalNodes.insert( nodeIdx );
    }
    else
    {
      expected.boundaryNodes.insert( nodeIdx );
    }
  }

  expected.numBoundaryNodesOnDomain = numBoundaryNodesTouchingDomain;

  for( Edge const & edge : fractureEdges )
  {
    bool const nodeADuplicated = (expected.internalNodes.find( edge.nodeA ) != expected.internalNodes.end());
    bool const nodeBDuplicated = (expected.internalNodes.find( edge.nodeB ) != expected.internalNodes.end());

    if( nodeADuplicated || nodeBDuplicated )
    {
      expected.internalEdges.insert( edge );
    }
    else
    {
      expected.boundaryEdges.insert( edge );
    }
  }

  // Step 6: Compute per-node duplication counts.
  expected.numNodesToDuplicate = expected.internalNodes.size();
  expected.numEdgesToDuplicate = expected.internalEdges.size();
  expected.numFacesToDuplicate = expected.numFractureFaces;

  bool const isBoundaryCutting = isFullSpan;
  expected.totalDuplicatedNodes = 0;

  for( localIndex nodeIdx : expected.internalNodes )
  {
    auto const & fractureIds = nodeToFractureIds[nodeIdx];
    integer const numFractures = static_cast< integer >( fractureIds.size() );
    integer newNodesForThisNode;

    if( !isBoundaryCutting )
    {
      // Non-boundary-cutting: GEOS splits a node once per fracture it is *interior* to.
      // N fractures cutting independently through a node partition the neighbourhood into
      // 2^N sectors → 2^N - 1 new nodes are required.
      integer numFracturesInteriorTo = 0;
      for( integer fid : fractureIds )
      {
        auto const it = perFractureBoundaryNodes.find( fid );
        bool const isOnFracBnd = (it != perFractureBoundaryNodes.end())
                                 && (it->second.count( nodeIdx ) > 0);
        if( !isOnFracBnd )
          ++numFracturesInteriorTo;
      }
      newNodesForThisNode = std::max( integer{1}, (1 << numFracturesInteriorTo) - 1 );
    }
    else if( numFractures == 1 )
    {
      // Boundary-cutting, single fracture: 1 new node per internal node.
      newNodesForThisNode = 1;
    }
    else
    {
      // Boundary-cutting, multiple fractures: original intersection logic.
      integer numInteriorTo = 0;
      integer numBndOf     = 0;

      for( integer fid : fractureIds )
      {
        auto const it = perFractureBoundaryNodes.find( fid );
        bool const isOnFracBnd = (it != perFractureBoundaryNodes.end())
                                 && (it->second.count( nodeIdx ) > 0);
        if( isOnFracBnd )
          ++numBndOf;
        else
          ++numInteriorTo;
      }

      if( numBndOf == 0 )
      {
        newNodesForThisNode = (1 << numFractures) - 1;
      }
      else if( numInteriorTo == 0 )
      {
        auto const cntIt = nodeCombBndEdgeCount.find( nodeIdx );
        localIndex const nbnd = (cntIt != nodeCombBndEdgeCount.end()) ? cntIt->second : 0;
        newNodesForThisNode = (nbnd > 0)
          ? static_cast< integer >( std::max( localIndex{1}, nbnd - 1 ) )
          : static_cast< integer >( std::max( integer{1}, numFractures - 1 ) );
      }
      else
      {
        newNodesForThisNode = numFractures;
      }
    }

    expected.totalDuplicatedNodes += newNodesForThisNode;
  }

  expected.totalDuplicatedEdges = expected.numEdgesToDuplicate;
  expected.totalDuplicatedFaces = expected.numFacesToDuplicate;

  return expected;
}

#endif // TESTSURFACEGENERATORCOMMON_HPP
