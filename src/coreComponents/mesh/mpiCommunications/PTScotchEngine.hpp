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
 * @file PTScotchEngine.hpp
 */

#ifndef GEOS_PARTITIONER_PTSCOTCHENGINE_HPP_
#define GEOS_PARTITIONER_PTSCOTCHENGINE_HPP_

#include "GraphPartitionEngine.hpp"

namespace geos
{

/**
 * @class PTScotchEngine
 * @brief Graph partitioning engine using the PT-Scotch library
 *
 * This engine wraps PT-Scotch (Parallel Scotch).
 */
class PTScotchEngine : public GraphPartitionEngine
{
public:

  /**
   * @brief Constructor
   * @param name The name of this engine instance
   * @param parent The parent group
   */
  explicit PTScotchEngine( string const & name,
                           dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~PTScotchEngine() override;

  /**
   * @brief Structure to hold XML input keys
   */
  struct viewKeyStruct
  {
    /// Partitioning strategy
    static constexpr char const * strategyString() { return "strategy"; }
  };

  /**
   * @brief Return the catalog name for factory registration
   * @return The string "PTScotch"
   */
  static string catalogName() { return "PTScotch"; }

  /**
   * @brief Partition a distributed graph using PT-Scotch
   *
   * Calls SCOTCH_dgraphPart to perform parallel graph partitioning.
   *
   * @param graph Edge connectivity (locally owned vertices)
   * @param vertDist Vertex distribution across ranks
   * @param numParts Target number of partitions
   * @param comm MPI communicator
   * @return Array mapping each local vertex to its target partition
   */
  array1d< pmet_idx_t > partition(
    ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
    arrayView1d< pmet_idx_t const > const & vertDist,
    pmet_idx_t const numParts,
    MPI_Comm comm ) override;


  /**
   * @brief Convert mesh to dual graph
   *
   * PT-Scotch doesn't have a native mesh-to-dual function, so we fall back
   * to ParMETIS implementation if available.
   */
  ArrayOfArrays< pmet_idx_t, pmet_idx_t >
  meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
              arrayView1d< pmet_idx_t const > const & elemDist,
              MPI_Comm comm,
              int const minCommonNodes ) override;

  /**
   * @brief Get the number of refinement iterations
   *
   * @note PT-Scotch doesn't expose refinements the same way as ParMETIS.
   *       This returns 0 to indicate it's not applicable.
   *
   * @return 0 (not applicable for PT-Scotch)
   */
  int getNumRefinements() const override { return 0; }

  /**
   * @brief Set the partitioning strategy
   * @param strategy Strategy string (e.g., "default", "quality", "speed")
   */
  void setStrategy( string const & strategy )
  {
    m_strategy = strategy;
  }

  /**
   * @brief Get the partitioning strategy
   * @return Strategy string
   */
  string const & getStrategy() const { return m_strategy; }

private:

  /// Partitioning strategy (default, quality, speed, etc.)
  string m_strategy;
};

} // namespace geos

#endif // GEOS_PARTITIONER_PTSCOTCHENGINE_HPP_
