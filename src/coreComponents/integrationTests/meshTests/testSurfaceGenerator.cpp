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
 * @file testSurfaceGenerator.cpp
 * @brief Integration tests for SurfaceGenerator topology modification (serial, no MPI).
 *
 * Verifies node duplication, fracture element creation, and Euler characteristic (χ = V-E+F-C)
 * for boundary-cutting and internal fracture configurations on hex and tet meshes.
 */

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

using namespace geos;

constexpr real64 COORDINATE_TOLERANCE = 1.0e-10;
constexpr real64 JACOBIAN_TOLERANCE = 1.0e-12;

// Global flag to control debug printing (set to false to disable verbose output)
static constexpr bool ENABLE_DEBUG_PRINTS = false;

CommandLineOptions g_commandLineOptions;

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
integer computeEulerCharacteristic( NodeManager const & nodeManager,
                                    FaceManager const & faceManager,
                                    ElementRegionManager const & elemManager )
{
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

/// @brief Compute χ = V - E + F - C from pre-computed topology stats.
integer computeEulerCharacteristic( TopologyStats const & stats )
{
  return stats.numNodes - stats.numEdges + stats.numFaces - stats.numCells;
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
 *        coordinate validity, Euler characteristic, and element Jacobians).
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
                                      NodeManager const & nodeManager,
                                      ElementRegionManager const & elemManager )
{
  // A1: Validate node duplication matches prediction
  // Note: totalDuplicatedNodes = total NEW nodes created (not counting original nodes)
  // For a node at N-fracture intersection: N new nodes are created
  GEOS_LOG_RANK_0( "Validating A1: Node duplication prediction (Expected: " << expected.totalDuplicatedNodes
                                                                            << ", Actual: " << statsAfter.numDuplicatedNodes << ")" );
  EXPECT_EQ( statsAfter.numDuplicatedNodes, expected.totalDuplicatedNodes )
    << "Test " << testCaseName << ": Node duplication prediction MISMATCH"
    << "\n  Expected new nodes: " << expected.totalDuplicatedNodes
    << "\n  Actual new nodes:   " << statsAfter.numDuplicatedNodes
    << "\n  Nodes before split: " << statsBefore.numNodes
    << "\n  Nodes after split:  " << statsAfter.numNodes
    << "\n  Fracture node sets: " << nodeSetNames;

  // A2: Verify fracture elements were created and match the predicted count
  GEOS_LOG_RANK_0( "Validating A2: Fracture elements created (Expected: " << expected.numFractureElements
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

  // A3: Verify nodes were duplicated
  GEOS_LOG_RANK_0( "Validating A3: Nodes were duplicated (Actual: " << statsAfter.numDuplicatedNodes << ")" );
  EXPECT_GT( statsAfter.numDuplicatedNodes, 0 )
    << "Test " << testCaseName << ": No nodes were duplicated"
    << "\n  Nodes before: " << statsBefore.numNodes
    << "\n  Nodes after:  " << statsAfter.numNodes
    << "\n  Predicted to split: " << expected.numNodesToDuplicate << " nodes"
    << "\n  Fracture elements created: " << statsAfter.numFractureElements
    << "\n  POSSIBLE CAUSES:"
    << "\n    - Fracture topology prediction found no internal nodes to split"
    << "\n    - All fracture nodes are on the boundary/tip (not duplicated)"
    << "\n    - SurfaceGenerator did not perform splitting";

  // A4: Validate all node coordinates are valid (not NaN)
  GEOS_LOG_RANK_0( "Validating A4: Node coordinates validity (Total nodes: " << nodeManager.size() << ")" );
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions =
    nodeManager.referencePosition();

  for( localIndex ni = 0; ni < nodeManager.size(); ++ni )
  {
    real64 const x = nodePositions[ni][0];
    real64 const y = nodePositions[ni][1];
    real64 const z = nodePositions[ni][2];

    EXPECT_FALSE( std::isnan( x ) || std::isnan( y ) || std::isnan( z ) )
      << "Test " << testCaseName << ": Node " << ni << " has invalid (NaN) coordinates"
      << "\n  Position: (" << x << ", " << y << ", " << z << ")"
      << "\n  Total nodes in mesh: " << nodeManager.size()
      << "\n  This indicates memory corruption or uninitialized node positions";
  }

  // A5: Validate expected Euler characteristic before split
  GEOS_LOG_RANK_0( "Validating A5: Euler χ before split (Expected: " << expectedEulerBefore
                                                                     << ", Actual: " << eulerCharBeforeSplit << ")" );
  EXPECT_EQ( eulerCharBeforeSplit, expectedEulerBefore )
    << "Test " << testCaseName << ": Euler characteristic MISMATCH before split"
    << "\n  Expected χ: " << expectedEulerBefore
    << "\n  Actual χ:   " << eulerCharBeforeSplit
    << "\n  Mesh file:  " << meshFileName
    << "\n  NOTE: For a conformal mesh, expected χ = 1 (single connected solid)";

  // A6: Validate expected Euler characteristic after split
  GEOS_LOG_RANK_0( "Validating A6: Euler χ after split (Expected: " << expectedEulerAfter
                                                                    << ", Actual: " << eulerCharAfterSplit << ")" );
  EXPECT_EQ( eulerCharAfterSplit, expectedEulerAfter )
    << "Test " << testCaseName << ": Euler characteristic MISMATCH after split"
    << "\n  Expected χ: " << expectedEulerAfter
    << "\n  Actual χ:   " << eulerCharAfterSplit
    << "\n  Mesh file:  " << meshFileName
    << "\n  Nodes duplicated: " << statsAfter.numDuplicatedNodes
    << "\n  Fracture elements: " << statsAfter.numFractureElements;

  // A7: Validate element Jacobians (no degenerate elements)
  GEOS_LOG_RANK_0( "Validating A7: Element Jacobians (checking for degenerate elements)" );
  checkElementJacobians( elemManager );
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
      allFractureNodes.insert( nodeIdx );
      nodesInThisFracture.insert( nodeIdx );
      nodeToFractureIds[nodeIdx].insert( fractureId );
    }

    if( ENABLE_DEBUG_PRINTS )
    {
      GEOS_LOG_RANK_0( "DEBUG: Fracture " << fractureId << " (" << setName << ") has "
                                          << nodesInThisFracture.size() << " nodes" );
    }

    // Identify fracture faces for this specific fracture
    for( localIndex fi = 0; fi < faceManager.size(); ++fi )
    {
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

  // Count intersection nodes (nodes shared by multiple fractures)
  if( ENABLE_DEBUG_PRINTS )
  {
    for( auto const & [nodeIdx, fractureIds] : nodeToFractureIds )
    {
      if( fractureIds.size() > 1 )
      {
        GEOS_LOG_RANK_0( "DEBUG: Node " << nodeIdx << " is at intersection of " << fractureIds.size() << " fractures" );
      }
    }
  }

  expected.numFractureNodes = allFractureNodes.size();
  expected.numFractureFaces = fractureFaces.size();
  expected.numFractureElements = fractureFaces.size();

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Found " << fractureFaces.size() << " fracture faces" );
  }

  // Step 3: Extract fracture edges from faces
  // For boundary detection, we need to know which edges belong to only 1 face
  std::map< Edge, localIndex > edgeToFaceCount;
  std::set< Edge > fractureEdges;

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Extracting edges from fracture faces" );
  }

  // Deduplicate faces before counting boundary edges to avoid classifying
  // T-junction intersection edges as interior (would inflate edgeToFaceCount).

  // Combined unique fracture faces (all fractures together)
  std::set< Facet > seenCombined;

  // Per-fracture unique face sets, keyed by fractureId
  std::map< integer, std::set< Facet > > seenPerFracture;

  // Per-fracture edge-to-face count (using unique faces only)
  std::map< integer, std::map< Edge, localIndex > > perFractureEdgeCount;

  for( localIndex fi = 0; fi < faceManager.size(); ++fi )
  {
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

      // Combined dedup: count each unique face only once for edgeToFaceCount
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

      // Per-fracture dedup: count each unique face only once per fracture
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

  // Build per-fracture boundary node sets.
  // A node is on the BOUNDARY of fracture F if it belongs to an edge that appears
  // exactly once in F's unique face set (i.e., the edge is on F's perimeter).
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

  // Build combined-perimeter edge count per node.
  // This is the number of combined-boundary edges (count==1 in the combined unique faces)
  // incident to each node. Used for the "boundary-of-all-fractures" node case.
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

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Found " << fractureEdges.size() << " fracture edges from " << fractureFaces.size() << " faces" );
  }

  // Step 4: Identify fracture boundary nodes
  // Fracture boundary edges are those that belong to only 1 face
  std::set< localIndex > fractureBoundaryNodes;
  std::set< localIndex > fractureInternalNodes;

  localIndex numBoundaryEdges = 0;
  localIndex numBoundaryEdgesOnDomain = 0;

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Analyzing fracture boundary edges:" );
  }

  for( auto const & [edge, count] : edgeToFaceCount )
  {
    if( count == 1 )  // Edge is on fracture boundary
    {
      ++numBoundaryEdges;
      fractureBoundaryNodes.insert( edge.nodeA );
      fractureBoundaryNodes.insert( edge.nodeB );

      // Check if this boundary edge is on the domain boundary
      bool const nodeAOnDomain = isNodeOnBoundary( edge.nodeA );
      bool const nodeBOnDomain = isNodeOnBoundary( edge.nodeB );
      bool const edgeOnDomain = nodeAOnDomain && nodeBOnDomain;

      if( ENABLE_DEBUG_PRINTS )
      {
        GEOS_LOG_RANK_0( "  Boundary edge (" << edge.nodeA << ", " << edge.nodeB << "): "
                                             << (edgeOnDomain ? "ON DOMAIN" : "INTERIOR") );
      }

      if( edgeOnDomain )
      {
        ++numBoundaryEdgesOnDomain;
      }
    }
  }

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Fracture boundary analysis:" );
    GEOS_LOG_RANK_0( "  Total boundary edges: " << numBoundaryEdges );
    GEOS_LOG_RANK_0( "  Boundary edges on domain: " << numBoundaryEdgesOnDomain );
    GEOS_LOG_RANK_0( "  Fracture boundary nodes found: " << fractureBoundaryNodes.size() );
  }

  // Heuristic: If ALL fracture boundary edges are on the domain boundary,
  // then this is a full-span fracture configuration with no interior tips
  bool const isFullSpan = (numBoundaryEdges > 0) && (numBoundaryEdges == numBoundaryEdgesOnDomain);

  if( isFullSpan )
  {
    if( ENABLE_DEBUG_PRINTS )
    {
      GEOS_LOG_RANK_0( "DEBUG: Detected FULL-SPAN fractures (all boundaries on domain)" );
      GEOS_LOG_RANK_0( "       -> ALL fracture nodes will be duplicated (no interior tips)" );
      GEOS_LOG_RANK_0( "       -> Clearing fractureBoundaryNodes set (was " << fractureBoundaryNodes.size() << " nodes)" );
    }
    // For full-span, clear the fracture boundary classification
    // All nodes will be treated as internal or intersection nodes
    fractureBoundaryNodes.clear();
  }
  else
  {
    if( ENABLE_DEBUG_PRINTS )
    {
      GEOS_LOG_RANK_0( "DEBUG: Detected PARTIAL fractures (interior tips exist)" );
      GEOS_LOG_RANK_0( "       -> " << fractureBoundaryNodes.size() << " nodes marked as fracture boundary" );
    }
  }

  // All other fracture nodes are internal
  for( localIndex nodeIdx : allFractureNodes )
  {
    if( fractureBoundaryNodes.find( nodeIdx ) == fractureBoundaryNodes.end() )
    {
      fractureInternalNodes.insert( nodeIdx );
    }
  }

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Found " << fractureInternalNodes.size() << " internal fracture nodes, "
                                     << fractureBoundaryNodes.size() << " fracture boundary nodes" );
  }

  // Step 5: Apply duplication rules
  // For full-span fractures, most nodes should be duplicated
  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Applying duplication rules to " << allFractureNodes.size() << " fracture nodes" );
  }

  localIndex numInternalFractureNodes = 0;
  localIndex numBoundaryNodesNotTouchingDomain = 0;
  localIndex numBoundaryNodesTouchingDomain = 0;
  localIndex numIntersectionNodes = 0;

  for( localIndex nodeIdx : allFractureNodes )
  {
    real64 const x = nodePositions[nodeIdx][0];
    real64 const y = nodePositions[nodeIdx][1];
    real64 const z = nodePositions[nodeIdx][2];

    integer const numFractures = nodeToFractureIds[nodeIdx].size();
    bool const isAtIntersection = (numFractures > 1);
    bool const isOnDomainBoundary = isNodeOnBoundary( nodeIdx );
    bool const isOnFractureBoundary = (fractureBoundaryNodes.find( nodeIdx ) != fractureBoundaryNodes.end());
    bool const isInternal = !isOnFractureBoundary;

    bool shouldDuplicate = false;
    std::string reason;

    // Simplified rules:
    // 1. Internal fracture nodes → ALWAYS duplicate
    // 2. Intersection nodes → ALWAYS duplicate (regardless of boundary status)
    // 3. Fracture boundary nodes touching domain → duplicate
    // 4. Fracture boundary nodes NOT touching domain → do NOT duplicate (fracture tips)

    if( isAtIntersection )
    {
      // Intersection nodes are ALWAYS duplicated
      shouldDuplicate = true;
      reason = "intersection node (" + std::to_string( numFractures ) + " fractures)";
      ++numIntersectionNodes;
    }
    else if( isInternal )
    {
      // Internal single-fracture nodes → duplicate
      shouldDuplicate = true;
      reason = "internal fracture node";
      ++numInternalFractureNodes;
    }
    else if( isOnFractureBoundary && isOnDomainBoundary )
    {
      // Fracture boundary touching domain → duplicate
      shouldDuplicate = true;
      reason = "fracture boundary touching domain";
      ++numBoundaryNodesTouchingDomain;
    }
    else if( isOnFractureBoundary && !isOnDomainBoundary )
    {
      // Fracture boundary NOT touching domain → do not duplicate (fracture tip)
      shouldDuplicate = false;
      reason = "fracture tip (boundary, no domain contact)";
      ++numBoundaryNodesNotTouchingDomain;
    }

    if( ENABLE_DEBUG_PRINTS )
    {
      GEOS_LOG_RANK_0( "  Node " << nodeIdx << " at (" << x << ", " << y << ", " << z << ") ["
                                 << numFractures << " fracture(s)]: "
                                 << reason << " -> " << (shouldDuplicate ? "DUPLICATE" : "NO DUPLICATE") );
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

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Classification summary:" );
    GEOS_LOG_RANK_0( "  Internal fracture nodes (duplicated): " << numInternalFractureNodes );
    GEOS_LOG_RANK_0( "  Intersection nodes (duplicated): " << numIntersectionNodes );
    GEOS_LOG_RANK_0( "  Boundary nodes touching domain (duplicated): " << numBoundaryNodesTouchingDomain );
    GEOS_LOG_RANK_0( "  Boundary nodes NOT touching domain (not duplicated): " << numBoundaryNodesNotTouchingDomain );
    GEOS_LOG_RANK_0( "  Total nodes to duplicate: " << expected.internalNodes.size() );
    GEOS_LOG_RANK_0( "  Total nodes NOT to duplicate: " << expected.boundaryNodes.size() );
  }

  // Store boundary-cutting indicator
  expected.numBoundaryNodesOnDomain = numBoundaryNodesTouchingDomain;

  // For edges, duplicate if at least one node is duplicated
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
  //
  // For boundary-cutting fractures, intersection nodes fall into three cases:
  //  A – interior to ALL fractures (true crossing):  new = 2^N - 1
  //  B – on boundary of ALL fractures (corner/edge): new = max(1, perimeter_edges - 1)
  //  C – mixed T-junction:                           new = N  (number of fractures)
  //
  // Single-fracture or non-full-span nodes always produce 1 new node.

  expected.numNodesToDuplicate = expected.internalNodes.size();
  expected.numEdgesToDuplicate = expected.internalEdges.size();
  expected.numFacesToDuplicate = expected.numFractureFaces;

  bool const isBoundaryCutting = isFullSpan;

  expected.totalDuplicatedNodes = 0;

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Computing duplication for each node (boundary-cutting: "
                     << (isBoundaryCutting ? "YES" : "NO") << "):" );
  }

  std::map< integer, localIndex > fractureCountDistribution;
  std::map< integer, localIndex > newNodesDistribution;

  for( localIndex nodeIdx : expected.internalNodes )
  {
    auto const & fractureIds = nodeToFractureIds[nodeIdx];
    integer const numFractures = static_cast< integer >( fractureIds.size() );
    integer newNodesForThisNode;

    if( !isBoundaryCutting || numFractures == 1 )
    {
      // Internal fracture or single-fracture: each node duplicated once
      newNodesForThisNode = 1;
    }
    else
    {
      // Classify this intersection node into Case A, B, or C
      integer numInteriorTo = 0;   // fractures for which this node is strictly interior
      integer numBndOf     = 0;    // fractures for which this node is on the boundary

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
        // Case A: interior to all fractures – true crossing
        newNodesForThisNode = (1 << numFractures) - 1;
      }
      else if( numInteriorTo == 0 )
      {
        // Case B: on boundary of all fractures – perimeter corner.
        // If the node sits on the combined fracture perimeter (nbnd > 0), the number of
        // new nodes equals the number of perimeter edges meeting at the node minus 1.
        // If the node is off the combined perimeter (nbnd == 0, e.g. interior Y-junction
        // nodes where all fractures terminate), use nf - 1.
        auto const cntIt = nodeCombBndEdgeCount.find( nodeIdx );
        localIndex const nbnd = (cntIt != nodeCombBndEdgeCount.end()) ? cntIt->second : 0;
        newNodesForThisNode = (nbnd > 0)
          ? static_cast< integer >( std::max( localIndex{1}, nbnd - 1 ) )
          : static_cast< integer >( std::max( integer{1}, numFractures - 1 ) );
      }
      else
      {
        // Case C: mixed – T-junction style
        newNodesForThisNode = numFractures;
      }
    }

    expected.totalDuplicatedNodes += newNodesForThisNode;
    fractureCountDistribution[numFractures]++;
    newNodesDistribution[newNodesForThisNode]++;

    if( ENABLE_DEBUG_PRINTS )
    {
      real64 const x = nodePositions[nodeIdx][0];
      real64 const y = nodePositions[nodeIdx][1];
      real64 const z = nodePositions[nodeIdx][2];
      GEOS_LOG_RANK_0( "  Node " << nodeIdx << " at (" << x << ", " << y << ", " << z << "): "
                                 << numFractures << " fracture(s) -> " << newNodesForThisNode << " new nodes" );
    }
  }

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Distribution summary:" );
    for( auto const & [numFractures, count] : fractureCountDistribution )
    {
      GEOS_LOG_RANK_0( "    " << count << " node(s) in " << numFractures << " fracture(s)" );
    }
    GEOS_LOG_RANK_0( "DEBUG: Total from summation: " << expected.totalDuplicatedNodes );
  }

  expected.totalDuplicatedEdges = expected.numEdgesToDuplicate;
  expected.totalDuplicatedFaces = expected.numFacesToDuplicate;

  if( ENABLE_DEBUG_PRINTS )
  {
    GEOS_LOG_RANK_0( "DEBUG: Final counts:" );
    GEOS_LOG_RANK_0( "  Unique nodes to duplicate: " << expected.numNodesToDuplicate );
    GEOS_LOG_RANK_0( "  Total new nodes created: " << expected.totalDuplicatedNodes );
    GEOS_LOG_RANK_0( "  Edges to duplicate: " << expected.totalDuplicatedEdges );
    GEOS_LOG_RANK_0( "  Faces to duplicate: " << expected.totalDuplicatedFaces );
  }

  return expected;
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
 * @brief Parametrized test fixture for SurfaceGenerator.
 *
 * Parameters: (test-name, mesh-file, node-set-names, χ_before, χ_after).
 */
class SurfaceGeneratorTest
  : public ::testing::TestWithParam< std::tuple< std::string, std::string, std::string, integer, integer > >
{
protected:
  void SetUp() override {}

  /// @brief Build the XML input string for this test case.
  std::string generateXmlInput( std::string const & meshFile,
                                std::string const & nodeSetNames )
  {
    std::ostringstream oss;
    oss << R"xml(<?xml version="1.0" ?>
<Problem>
  <Mesh>
    <VTKMesh name="mesh1" file=")xml" << meshFile << R"xml(" nodesetNames=")xml" << nodeSetNames <<
      R"xml(" />
  </Mesh>
  
  <Solvers gravityVector="{0.0, 0.0, 0.0}">
    <SurfaceGenerator 
      name="SurfaceGen" 
      targetRegions="{ Region }"
      fractureRegion="Fracture"
      initialRockToughness="1.0" 
      logLevel="0"/>
  </Solvers>
  
  <NumericalMethods>
    <FiniteElements>
      <FiniteElementSpace name="FE1" order="1"/>
    </FiniteElements>
  </NumericalMethods>
  
  <ElementRegions>
    <CellElementRegion name="Region" cellBlocks="{ * }" materialList="{ emptyConstitutive }"/>
    <SurfaceElementRegion 
      name="Fracture" 
      defaultAperture="1.0e-4" 
      materialList="{ emptyConstitutive }"/>
  </ElementRegions>
  
  <Constitutive>
    <NullModel name="emptyConstitutive"/>
  </Constitutive>
  
  <FieldSpecifications>
    <FieldSpecification 
      name="separableFace" 
      fieldName="isFaceSeparable" 
      initialCondition="1" 
      setNames=")xml" << nodeSetNames << R"xml(" 
      objectPath="faceManager" 
      scale="1"/>
    <FieldSpecification 
      name="frac" 
      initialCondition="1" 
      setNames=")xml" << nodeSetNames <<
      R"xml(" 
      objectPath="faceManager" 
      fieldName="ruptureState" 
      scale="1"/>
  </FieldSpecifications>
  <Events maxTime="1.0e-10">
    <SoloEvent name="preFracture" target="/Solvers/SurfaceGen"/>
  </Events>
</Problem>
)xml";
    return oss.str();
  }

  /// @brief Execute one parametrized test case end-to-end.
  void runTest( std::string const & testCaseName,
                std::string const & meshFileName,
                std::string const & nodeSetNames,
                integer const expectedEulerBefore,
                integer const expectedEulerAfter )
  {
    std::string const xmlInput = generateXmlInput(
      meshFileName,  // Mesh files are copied to same directory as executable
      nodeSetNames );

    // Write XML to temporary file
    std::string const xmlFileName = "test_surface_gen_" + testCaseName + ".xml";
    std::ofstream xmlFile( xmlFileName );
    xmlFile << xmlInput;
    xmlFile.close();

    // Setup GEOS
    std::unique_ptr< CommandLineOptions > options = std::make_unique< CommandLineOptions >( g_commandLineOptions );
    options->inputFileNames.clear();
    options->inputFileNames.push_back( xmlFileName );

    GeosxState state( std::move( options ) );
    ASSERT_TRUE( state.initializeDataRepository() )
      << "Test " << testCaseName << ": Failed to initialize data repository for mesh file '"
      << meshFileName << "' with node sets " << nodeSetNames;
    state.applyInitialConditions();

    // Get mesh before surface generation
    ProblemManager & pm = state.getProblemManager();
    DomainPartition & domain = pm.getDomainPartition();
    MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

    NodeManager & nodeManager = mesh.getNodeManager();
    EdgeManager & edgeManager = mesh.getEdgeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // Store initial statistics
    TopologyStats statsBefore;
    statsBefore.numNodes = nodeManager.size();
    statsBefore.numEdges = edgeManager.size();
    statsBefore.numFaces = faceManager.size();
    statsBefore.numCells = 0;
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
    {
      statsBefore.numCells += subRegion.size();
    } );

    // **NEW: Preprocess fracture topology to predict duplication BEFORE mesh split**
    std::vector< std::string > fractureNodeSets = parseNodeSetNames( nodeSetNames );
    ExpectedDuplication expected = preprocessFractureTopology( mesh, fractureNodeSets );

    GEOS_LOG_RANK_0( "========================================" );
    GEOS_LOG_RANK_0( "Test: " << testCaseName );
    GEOS_LOG_RANK_0( "Expected duplication:" );
    GEOS_LOG_RANK_0( "  Nodes to split: " << expected.numNodesToDuplicate
                                          << " (new nodes created: " << expected.totalDuplicatedNodes << ")" );
    GEOS_LOG_RANK_0( "  Edges to split: " << expected.numEdgesToDuplicate
                                          << " (new edges created: " << expected.totalDuplicatedEdges << ")" );
    GEOS_LOG_RANK_0( "  Faces to split: " << expected.numFacesToDuplicate
                                          << " (new faces created: " << expected.totalDuplicatedFaces << ")" );

    // Compute Euler-Poincaré characteristic χ = V - E + F - C (should be 1 for connected solid).
    //
    // Two independent implementations are compared here:
    //   1. Local (element-enumeration) version: re-derives V, E, F from the elem→node map.
    //      This is the reference used throughout the test.
    //   2. MeshGeneratorBase version: uses nodeManager.size(), edgeManager.size(),
    //      faceManager.size() and CellElementSubRegion sizes directly from the managers.
    //      This is the production implementation we want to exercise.
    // Both must agree.
    integer const eulerCharBeforeSplit = computeEulerCharacteristic( nodeManager, faceManager, elemManager );
    GEOS_LOG_RANK_0( "  Euler characteristic before split (local):   " << eulerCharBeforeSplit );

    integer const eulerCharBeforeSplitMgr = computeEulerCharacteristic( nodeManager, edgeManager, faceManager, elemManager );
    GEOS_LOG_RANK_0( "  Euler characteristic before split (managers): " << eulerCharBeforeSplitMgr );

    EXPECT_EQ( eulerCharBeforeSplit, eulerCharBeforeSplitMgr )
      << "Test " << testCaseName << ": computeEulerCharacteristic MISMATCH between implementations"
      << "\n  Local (elem-enumeration): " << eulerCharBeforeSplit
      << "\n  Manager-based:            " << eulerCharBeforeSplitMgr
      << "\n  Mesh file: " << meshFileName
      << "\n  The manager-based version (MeshGeneratorBase) returned a different result than the"
      << "\n  reference element-enumeration version. Check that EdgeManager and FaceManager counts"
      << "\n  include only bulk entities (no fracture / surface sub-region contributions).";

    // Run the simulation (executes SurfaceGenerator and splits the mesh)
    state.run();

    // Get statistics after surface generation
    TopologyStats statsAfter;
    statsAfter.numNodes = nodeManager.size();
    statsAfter.numEdges = edgeManager.size();
    statsAfter.numFaces = faceManager.size();
    statsAfter.numCells = 0;
    statsAfter.numFractureElements = 0;

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
    {
      statsAfter.numCells += subRegion.size();
    } );

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
    {
      statsAfter.numFractureElements += subRegion.size();
    } );

    statsAfter.numDuplicatedNodes = statsAfter.numNodes - statsBefore.numNodes;

    // **NEW: AFTER mesh split - Compare expeted vs computed values**
    GEOS_LOG_RANK_0( "Actual duplication:" );
    GEOS_LOG_RANK_0( "  New nodes created: " << statsAfter.numDuplicatedNodes );
    GEOS_LOG_RANK_0( "  Fracture elements created: " << statsAfter.numFractureElements );

    // Compute Euler-Poincaré characteristic χ = V - E + F - C after split
    integer const eulerCharAfterSplit = computeEulerCharacteristic( nodeManager, faceManager, elemManager );
    GEOS_LOG_RANK_0( "  Euler characteristic after split: " << eulerCharAfterSplit );

    // Store Euler characteristic
    statsAfter.eulerCharacteristic = eulerCharAfterSplit;
    statsAfter.numBodies = eulerCharAfterSplit;

    // All validations are performed in one place
    validateSurfaceGeneratorResults( testCaseName,
                                     meshFileName,
                                     nodeSetNames,
                                     statsBefore,
                                     statsAfter,
                                     expected,
                                     expectedEulerBefore,
                                     eulerCharBeforeSplit,
                                     expectedEulerAfter,
                                     eulerCharAfterSplit,
                                     nodeManager,
                                     elemManager );

    // Log statistics
    GEOS_LOG_RANK_0( "Summary:" );
    GEOS_LOG_RANK_0( "  Nodes before: " << statsBefore.numNodes << ", after: " << statsAfter.numNodes
                                        << " (+" << statsAfter.numDuplicatedNodes << ")" );
    GEOS_LOG_RANK_0( "  Fracture elements: " << statsAfter.numFractureElements );
    GEOS_LOG_RANK_0( "  Separate bodies: " << statsAfter.numBodies );
    GEOS_LOG_RANK_0( "========================================" );

    // Cleanup
    std::remove( xmlFileName.c_str() );
  }
};

/// @brief Test body — dispatches to runTest().
TEST_P( SurfaceGeneratorTest, TopologyValidation )
{
  auto const & params = GetParam();
  std::string const & testCaseName = std::get< 0 >( params );
  std::string const & meshFileName = std::get< 1 >( params );
  std::string const & nodeSetNames = std::get< 2 >( params );
  integer const expectedEulerBefore = std::get< 3 >( params );
  integer const expectedEulerAfter = std::get< 4 >( params );

  runTest( testCaseName, meshFileName, nodeSetNames, expectedEulerBefore, expectedEulerAfter );
}

// ---------------------------------------------------------------------------
// Note on expected Euler characteristic for NoBoundaryCut tests
// ---------------------------------------------------------------------------
// The mesh is a 4x4x4 hexahedral grid (64 cells, 125 nodes) in [0,1]^3.
// Each DFN fracture is a 2x2 patch of 4 quad faces fully embedded in the
// interior of the domain (no face touches the boundary).
//
// Before any split the solid has a single connected component:
//   χ = V - E + F - C = 1
//
// Each fracture contains exactly ONE interior node (the center of the 2x2
// patch).  Splitting that node "unzips" the internal connectivity:
//   ΔV = +1  (the center node is duplicated)
//   ΔE = +4  (the four in-plane edges radiating from the center are duplicated)
//   ΔF = +4  (the four quad faces sharing the center are duplicated)
//   ΔC =  0  (no new 3-D cells are created)
//
//   Δχ = ΔV - ΔE + ΔF - ΔC = 1 - 4 + 4 - 0 = +1
//   χ_after = 1 + 1 = 2
//
// The result χ = 2 holds for every NoBoundaryCut case (DFN1, DFN12,
// DFN13, DFN123) because all three fractures share the same single
// interior node at (0.5, 0.5, 0.5); activating additional fractures
// splits that node once more per fracture plane, but the incremental
// Δχ per additional split remains +1 only when a genuinely new interior
// node is processed — and in this mesh all subsequent fracture planes
// share the already-counted center node, so χ stays at 2.
// ---------------------------------------------------------------------------
// clang-format off
INSTANTIATE_TEST_SUITE_P(
  SurfaceGeneratorCases,
  SurfaceGeneratorTest,
  ::testing::Values(

    // =======================================================================
    // dfn_market meshes
    // -----------------------------------------------------------------------
    // Naming convention:
    //   fractured[_wavy][_full_span]_mesh_{hex|tet}_DFN_{suffix}.vtu
    //
    // Column layout:  test-name,  mesh-file,  node-sets,  χ_before, χ_after
    //
    // Node-set rules:
    //   • _DFN_1   → { f1_node_set }
    //   • _DFN_2   → { f2_node_set }
    //   • _DFN_3   → { f3_node_set }
    //   • _DFN_12  → { f1_node_set, f2_node_set }
    //   • _DFN_13  → { f1_node_set, f3_node_set }
    //   • _DFN_23  → { f2_node_set, f3_node_set }
    //   • _DFN_123 → { f1_node_set, f2_node_set, f3_node_set }

    // -----------------------------------------------------------------------
    // flat · no-boundary-cut · hex
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                             node sets                         χ_b  χ_a
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_1", "fractured_mesh_hex_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_2", "fractured_mesh_hex_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_3", "fractured_mesh_hex_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_12", "fractured_mesh_hex_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_13", "fractured_mesh_hex_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_23", "fractured_mesh_hex_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_hex_DFN_123", "fractured_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // flat · no-boundary-cut · tet
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                             node sets                         χ_b  χ_a
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_1", "fractured_mesh_tet_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_2", "fractured_mesh_tet_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_3", "fractured_mesh_tet_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_12", "fractured_mesh_tet_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_13", "fractured_mesh_tet_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_23", "fractured_mesh_tet_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_NoBndCut_tet_DFN_123", "fractured_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // flat · full-span · hex
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                                      node sets                        χ_b
    //  χ_a
    std::make_tuple( "Mkt_BndCut_hex_DFN_1", "fractured_full_span_mesh_hex_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_2", "fractured_full_span_mesh_hex_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_3", "fractured_full_span_mesh_hex_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_12", "fractured_full_span_mesh_hex_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_13", "fractured_full_span_mesh_hex_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_23", "fractured_full_span_mesh_hex_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_BndCut_hex_DFN_123", "fractured_full_span_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // flat · full-span · tet
    // -----------------------------------------------------------------------
    //                   test name                      mesh file                                      node sets                        χ_b
    //  χ_a
    std::make_tuple( "Mkt_BndCut_tet_DFN_1", "fractured_full_span_mesh_tet_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_2", "fractured_full_span_mesh_tet_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_3", "fractured_full_span_mesh_tet_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_12", "fractured_full_span_mesh_tet_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_13", "fractured_full_span_mesh_tet_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_23", "fractured_full_span_mesh_tet_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_BndCut_tet_DFN_123", "fractured_full_span_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // wavy · no-boundary-cut · hex
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                  node sets                        χ_b
    //  χ_a
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_1", "fractured_wavy_mesh_hex_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_2", "fractured_wavy_mesh_hex_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_3", "fractured_wavy_mesh_hex_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_12", "fractured_wavy_mesh_hex_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_13", "fractured_wavy_mesh_hex_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_23", "fractured_wavy_mesh_hex_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_hex_DFN_123", "fractured_wavy_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // wavy · no-boundary-cut · tet
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                  node sets                        χ_b
    //  χ_a
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_1", "fractured_wavy_mesh_tet_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_2", "fractured_wavy_mesh_tet_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_3", "fractured_wavy_mesh_tet_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_12", "fractured_wavy_mesh_tet_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_13", "fractured_wavy_mesh_tet_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_23", "fractured_wavy_mesh_tet_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyNoBndCut_tet_DFN_123", "fractured_wavy_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 2 ),

    // -----------------------------------------------------------------------
    // wavy · full-span · hex
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                           node sets
    //                       χ_b  χ_a
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_1", "fractured_wavy_full_span_mesh_hex_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_2", "fractured_wavy_full_span_mesh_hex_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_3", "fractured_wavy_full_span_mesh_hex_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_12", "fractured_wavy_full_span_mesh_hex_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_13", "fractured_wavy_full_span_mesh_hex_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_23", "fractured_wavy_full_span_mesh_hex_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_hex_DFN_123", "fractured_wavy_full_span_mesh_hex_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // wavy · full-span · tet
    // -----------------------------------------------------------------------
    //                   test name                           mesh file                                           node sets
    //                       χ_b  χ_a
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_1", "fractured_wavy_full_span_mesh_tet_DFN_1.vtu", "{ f1_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_2", "fractured_wavy_full_span_mesh_tet_DFN_2.vtu", "{ f2_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_3", "fractured_wavy_full_span_mesh_tet_DFN_3.vtu", "{ f3_node_set }", 1, 2 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_12", "fractured_wavy_full_span_mesh_tet_DFN_12.vtu", "{ f1_node_set, f2_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_13", "fractured_wavy_full_span_mesh_tet_DFN_13.vtu", "{ f1_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_23", "fractured_wavy_full_span_mesh_tet_DFN_23.vtu", "{ f2_node_set, f3_node_set }", 1, 4 ),
    std::make_tuple( "Mkt_WavyBndCut_tet_DFN_123", "fractured_wavy_full_span_mesh_tet_DFN_123.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 8 ),

    // -----------------------------------------------------------------------
    // T-shaped · boundary-cutting · hex & tet
    // -----------------------------------------------------------------------
    // F1 fully spans the domain (boundary-cutting); F2 is the stem terminating at F1.
    // The junction splits the solid into 3 parts → χ_after = 3.
    //                   test name                          mesh file                             node sets             χ_b  χ_a
    std::make_tuple( "Mkt_BndCut_t_shaped_hex_DFN_12", "t_shaped_wavy_mesh_hex_DFN_t1t2.vtu", "{ f1_node_set, f2_node_set }", 1, 3 ),
    std::make_tuple( "Mkt_BndCut_t_shaped_tet_DFN_12", "t_shaped_wavy_mesh_tet_DFN_t1t2.vtu", "{ f1_node_set, f2_node_set }", 1, 3 ),

    // -----------------------------------------------------------------------
    // Y-shaped · boundary-cutting · hex & tet
    // -----------------------------------------------------------------------
    // Three boundary-cutting fractures meeting along a common line (Y cross-section).
    // Splits the solid into 3 parts → χ_after = 3.
    //                   test name                          mesh file                                  node sets
    //                          χ_b  χ_a
    std::make_tuple( "Mkt_BndCut_Y_shaped_hex_DFN_123", "y_shaped_wavy_mesh_hex_DFN_y1y2y3.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 3 ),
    std::make_tuple( "Mkt_BndCut_Y_shaped_tet_DFN_123", "y_shaped_wavy_mesh_tet_DFN_y1y2y3.vtu", "{ f1_node_set, f2_node_set, f3_node_set }", 1, 3 )

    )
  );
// clang-format on

int main( int argc, char * argv[] )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv, false );

  int result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
