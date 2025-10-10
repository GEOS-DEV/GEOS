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
 * @file NoOpEngine.hpp
 */

#ifndef GEOS_PARTITIONER_NOOPENGINE_HPP_
#define GEOS_PARTITIONER_NOOPENGINE_HPP_

#include "GraphPartitionEngine.hpp"

namespace geos
{

/**
 * @class NoOpEngine
 * @brief Graph partition engine that does not repartition
 *
 * NoOpEngine is a graph partitioning engine that:
 * - Does NOT repartition the graph
 * - Assigns each element to its current MPI rank
 * - Useful for testing or when mesh is already well-partitioned
 */
class NoOpEngine : public GraphPartitionEngine
{
public:

  /**
   * @brief Constructor
   * @param name The name of this engine instance
   * @param parent The parent group
   */
  explicit NoOpEngine( string const & name,
                       dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~NoOpEngine() override = default;

  /**
   * @brief Catalog name for factory registration
   * @return The catalog key
   */
  static string catalogName() { return "noop"; }

  /**
   * @brief Partition a graph (no-op: assigns each element to current rank)
   *
   * This method does NOT repartition. It simply assigns all local elements
   * to the current MPI rank, effectively preserving the existing distribution.
   *
   * @param graph Dual graph adjacency (ignored)
   * @param vertDist Vertex distribution across ranks
   * @param numParts Target number of partitions (ignored)
   * @param comm MPI communicator
   * @param numRefinements Number of refinement iterations (ignored)
   * @return Partition assignment: all local elements → current rank
   */
  array1d< pmet_idx_t > partition(
    ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
    arrayView1d< pmet_idx_t const > const & vertDist,
    pmet_idx_t const numParts,
    MPI_Comm const comm ) override;

  /**
   * @brief Build dual graph from mesh (fallback to ParMETIS)
   *
   * NoOpEngine does not implement its own mesh-to-dual conversion.
   * Instead, it falls back to ParMETIS implementation.
   *
   * @param elemToNodes Element-to-node connectivity
   * @param elemDist Element distribution across ranks
   * @param comm MPI communicator
   * @param minCommonNodes Minimum shared nodes for dual edge
   * @return Dual graph adjacency
   *
   * @throws If ParMETIS is not available (GEOS_USE_PARMETIS=OFF)
   */
  ArrayOfArrays< pmet_idx_t, pmet_idx_t >
  meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
              arrayView1d< pmet_idx_t const > const & elemDist,
              MPI_Comm const comm,
              int const minCommonNodes ) override;


  /**
   * @brief Get the number of refinement iterations
   */
  int getNumRefinements() const override { return 0; }
};

} // namespace geos

#endif // GEOS_PARTITIONER_NOOPENGINE_HPP_
