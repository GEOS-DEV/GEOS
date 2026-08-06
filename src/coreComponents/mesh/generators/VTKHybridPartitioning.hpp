/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKHybridPartitioning.hpp
 */

#ifndef GEOS_MESH_GENERATORS_VTKHYBRIDPARTITIONING_HPP_
#define GEOS_MESH_GENERATORS_VTKHYBRIDPARTITIONING_HPP_

#include "common/DataTypes.hpp"
#include "mesh/generators/METISInterface.hpp"

class vtkUnstructuredGrid;

namespace geos
{
namespace vtk
{

struct SuperCellInfo;

/**
 * @brief Configuration for root-local mixed FVM/FEM partitioning.
 */
struct HybridPartitionOptions
{
  HybridPartitionOptions()
    : imbalance( 3 )
  {
    imbalance.setValues< serialPolicy >( 0.05 );
  }

  /// Relative importance of exact cut faces in the objective and graph.
  real64 fvmCommunicationWeight = 1.0;

  /// Relative importance of shared-point replication in the objective and graph.
  real64 femCommunicationWeight = 1.0;

  /// Weak penalty for each distinct pair of neighboring ranks.
  real64 neighborPenalty = 0.0;

  /// Per-constraint relative imbalance tolerance: FVM, FEM, memory.
  array1d< real64 > imbalance;

  /// Deterministic METIS seed.
  pmet_idx_t seed = 2022;

  /// Optional scalar VTK cell arrays overriding topology-derived work estimates.
  string fvmWeightField;
  string femWeightField;
  string memoryWeightField;

  /// Explicit root graph memory limit in MiB; zero means unlimited.
  int64_t rootGraphMemoryLimitMB = 0;

  /// Number of exact objective-decreasing refinement passes.
  integer refinementPasses = 1;

  /// Backward-compatible additive cost for fracture-connected super-cells.
  integer fractureWeight = 0;

  /// Print detailed per-rank topology and load diagnostics when nonzero.
  integer diagnostics = 0;
};

/**
 * @brief Exact mixed-cell topology and its weighted sparse METIS graph.
 *
 * Rows in @c pointToVertices contain unique partition vertices incident on one
 * global point. Rows in @c vertexToPoints and @c vertexToFaces are the reverse
 * incidence needed by exact refinement. Super-cells are already contracted.
 */
struct HybridPartitionTopology
{
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > graph;
  array1d< pmet_idx_t > edgeWeights;
  array2d< pmet_idx_t > vertexWeights;

  array1d< localIndex > cellToVertex;
  ArrayOfArrays< localIndex, localIndex > vertexToCells;
  array1d< globalIndex > representativeGlobalCellIds;

  /// Every exact manifold interior face, including faces internal to a super-cell.
  array2d< localIndex > faceToVertices;
  array1d< real64 > exactFaceWeights;
  ArrayOfArrays< localIndex, localIndex > vertexToFaces;

  /// Exact global-point incidence after duplicate vertex removal.
  ArrayOfArrays< localIndex, localIndex > pointToVertices;
  array1d< globalIndex > pointGlobalIds;
  array1d< real64 > exactPointWeights;
  ArrayOfArrays< localIndex, localIndex > vertexToPoints;

  int64_t estimatedPeakBytes = 0;
};

/**
 * @brief Partition quality, balance, timing, and memory diagnostics.
 */
struct HybridPartitionMetrics
{
  int64_t numCells = 0;
  int64_t numPartitionVertices = 0;
  int64_t numGraphEdges = 0;
  int64_t numInteriorFaces = 0;
  int64_t numCutFaces = 0;
  real64 weightedFaceCut = 0.0;
  int64_t numSharedPoints = 0;
  int64_t pointReplication = 0;
  real64 weightedPointReplication = 0.0;
  int64_t maxRankNeighbors = 0;
  real64 averageRankNeighbors = 0.0;
  int64_t maxConnectedComponents = 0;
  array1d< real64 > imbalanceByConstraint;

  real64 initialObjective = 0.0;
  real64 finalObjective = 0.0;
  int64_t initialCutFaces = 0;
  int64_t initialPointReplication = 0;
  int64_t refinementMoves = 0;

  real64 graphBuildSeconds = 0.0;
  real64 initialPartitionSeconds = 0.0;
  real64 exactRefinementSeconds = 0.0;
  real64 redistributionSeconds = 0.0;
  int64_t estimatedRootGraphPeakBytes = 0;
};

/**
 * @brief Result of root-local hybrid partitioning.
 */
struct HybridPartitionResult
{
  /// Target rank for every input 3D cell.
  array1d< int64_t > cellParts;

  /// Exact first-neighbor lists for every target rank.
  stdVector< stdVector< int > > rankNeighbors;

  HybridPartitionMetrics metrics;
};

/**
 * @brief Check whether all cells use supported linear 3D VTK topologies.
 * @param mesh Root-local 3D unstructured grid.
 * @param reason Filled with a precise fallback reason on failure.
 * @return Whether hybrid topology construction is supported.
 */
bool isHybridPartitioningSupported( vtkUnstructuredGrid const & mesh,
                                    string & reason );

/**
 * @brief Conservatively estimate temporary and retained root graph memory.
 */
int64_t estimateHybridPartitionMemory( vtkUnstructuredGrid const & mesh );

/**
 * @brief Build exact faces/points, contracted vertices, weights, and symmetric CSR.
 */
HybridPartitionTopology buildHybridPartitionTopology(
  vtkUnstructuredGrid const & mesh,
  SuperCellInfo const * superCells,
  HybridPartitionOptions const & options );

/**
 * @brief Evaluate the exact configured hybrid objective for a vertex partition.
 */
real64 evaluateHybridObjective( HybridPartitionTopology const & topology,
                                arrayView1d< int64_t const > vertexParts,
                                int numParts,
                                HybridPartitionOptions const & options );

/**
 * @brief Partition a root-local 3D mesh and expand contracted vertices to cells.
 */
HybridPartitionResult partitionHybridMeshOnRoot(
  vtkUnstructuredGrid const & cells3D,
  SuperCellInfo const * superCells,
  HybridPartitionOptions const & options,
  int numParts );

} // namespace vtk
} // namespace geos

#endif // GEOS_MESH_GENERATORS_VTKHYBRIDPARTITIONING_HPP_
