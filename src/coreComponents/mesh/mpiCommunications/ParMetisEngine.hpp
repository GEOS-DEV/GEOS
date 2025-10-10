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
 * @file ParMetisEngine.hpp
 */

#ifndef GEOS_PARTITIONER_PARMETISENGINE_HPP_
#define GEOS_PARTITIONER_PARMETISENGINE_HPP_

#include "GraphPartitionEngine.hpp"

namespace geos
{

/**
 * @class ParMetisEngine
 * @brief Graph partitioning engine using the ParMETIS library
 */
class ParMetisEngine : public GraphPartitionEngine
{
public:

  /**
   * @brief Constructor
   * @param name The name of this engine instance
   * @param parent The parent group
   */
  explicit ParMetisEngine( string const & name,
                           dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~ParMetisEngine() override;

  /**
   * @brief Structure to hold XML input keys
   */
  struct viewKeyStruct
  {
    /// Number of refinement iterations
    static constexpr char const * numRefinementsString() { return "numRefinements"; }
  };

  /**
   * @brief Return the catalog name for factory registration
   * @return The string "ParMetis"
   */
  static string catalogName() { return "ParMetis"; }

  /**
   * @brief Partition a distributed graph using ParMETIS
   *
   * Calls ParMETIS_V3_PartKway to perform multilevel k-way partitioning.
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
   * @brief Convert element-node connectivity to element-element dual graph
   *
   * Given a mesh represented by element-node connectivity, construct the dual graph
   * where:
   * - Vertices = elements
   * - Edges = element pairs that share at least minCommonNodes nodes
   *
   * This is a utility function used by MeshPartitioner to build the graph from a mesh.
   *
   * @param elemToNodes Element-to-node connectivity map.
   *                    elemToNodes[i] = global node IDs of element i.
   * @param elemDist Element distribution across ranks.
   *                 elemDist[rank] = global index of first element on rank.
   * @param comm MPI communicator
   * @param minCommonNodes Minimum number of shared nodes to create an edge (typically 3 for 3D)
   * @return Dual graph where vertices=elements, edges=element adjacency
   *
   * @note This is a static utility, not part of the core partition() interface
   */
  /**
   * @brief Convert mesh to dual graph using ParMETIS (instance method)
   */
  ArrayOfArrays< pmet_idx_t, pmet_idx_t >
  meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
              arrayView1d< pmet_idx_t const > const & elemDist,
              MPI_Comm comm,
              int const minCommonNodes ) override;

  /**
   * @brief Static helper: ParMETIS mesh-to-dual implementation
   *
   * This can be called by other engines (e.g., PTScotch) as a fallback.
   */
  static ArrayOfArrays< pmet_idx_t, pmet_idx_t >
  parmetisMeshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
                      arrayView1d< pmet_idx_t const > const & elemDist,
                      MPI_Comm comm,
                      int const minCommonNodes );

  /**
   * @brief Get the number of refinement iterations
   * @return Number of refinement iterations
   */
  int getNumRefinements() const override { return m_numRefinements; }

  /**
   * @brief Set the number of refinement iterations
   * @param numRefinements Number of refinement iterations
   */
  void setNumRefinements( int const numRefinements )
  {
    GEOS_ERROR_IF( numRefinements < 0, "Number of refinements must be non-negative" );
    m_numRefinements = numRefinements;
  }

private:

  /// Number of refinement iterations
  int m_numRefinements;
};

} // namespace geos

#endif // GEOS_PARTITIONER_PARMETISENGINE_HPP_
