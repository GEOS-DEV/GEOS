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
 * @file GraphPartitionEngine.hpp
 */

#ifndef GEOS_PARTITIONER_GRAPHPARTITIONENGINE_HPP_
#define GEOS_PARTITIONER_GRAPHPARTITIONENGINE_HPP_

#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"
#include "LvArray/src/ArrayOfArrays.hpp"

namespace geos
{

#if defined(GEOS_USE_HIP) // still need int32 hypre for the current hip-capable build
/// Typedef to allow us to specify required parmetis/scotch integer type.
using pmet_idx_t = int32_t;
#else
/// Typedef to allow us to specify required parmetis/scotch integer type.
using pmet_idx_t = int64_t;
#endif

/**
 * @class GraphPartitionEngine
 * @brief Abstract interface for low-level graph partitioning algorithms
 *
 * This is a LOW-LEVEL engine for pure algorithms: graph -> partition IDs
 */
class GraphPartitionEngine : public dataRepository::Group
{
public:

  /**
   * @brief Constructor
   * @param name The name of this engine instance
   * @param parent The parent group
   */
  explicit GraphPartitionEngine( string const & name,
                                 dataRepository::Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~GraphPartitionEngine() override;

  /**
   * @brief Partition a distributed graph
   *
   * @param graph Edge connectivity (locally owned vertices).
   *              graph[i] contains the global indices of neighbors of local vertex i.
   * @param vertDist Vertex distribution across ranks.
   *                 vertDist[rank] = global index of first vertex on rank.
   *                 Length: commSize + 1
   * @param numParts Target number of partitions (typically = commSize)
   * @param comm MPI communicator
   * @return Array mapping each local vertex to its target partition ID [0, numParts)
   *
   * @pre graph.size() == vertDist[myRank+1] - vertDist[myRank]
   * @post return.size() == graph.size()
   * @post All return values in [0, numParts)
   */
  virtual array1d< pmet_idx_t > partition(
    ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
    arrayView1d< pmet_idx_t const > const & vertDist,
    pmet_idx_t numParts,
    MPI_Comm comm ) = 0;

  /**
   * @brief Convert element-node connectivity to element-element dual graph
   *
   * Given a mesh represented by element-node connectivity, construct the dual graph
   * where:
   * - Vertices = elements
   * - Edges = element pairs that share at least minCommonNodes nodes
   * Each engine must provide an implementation (typically using ParMETIS as fallback).
   *
   * @param elemToNodes Element-to-node connectivity map.
   *                    elemToNodes[i] = global node IDs of element i.
   * @param elemDist Element distribution across ranks.
   *                 elemDist[rank] = global index of first element on rank.
   * @param comm MPI communicator
   * @param minCommonNodes Minimum number of shared nodes to create an edge (typically 3 for 3D)
   * @return Dual graph where vertices=elements, edges=element adjacency
   */
  virtual ArrayOfArrays< pmet_idx_t, pmet_idx_t >
  meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
              arrayView1d< pmet_idx_t const > const & elemDist,
              MPI_Comm comm,
              int minCommonNodes ) = 0;

  /**
   * @brief Get the number of refinement iterations
   *
   * @return Number of refinement iterations
   */
  virtual int getNumRefinements() const = 0;

  /**
   * @brief Accessor for the singleton Catalog object
   *
   * @note This is an INTERNAL catalog, separate from the DomainPartitioner catalog.
   *
   * @return A reference to the Catalog object
   */
  using CatalogInterface = dataRepository::CatalogInterface< GraphPartitionEngine,
                                                             string const &,
                                                             dataRepository::Group * const >;

  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Get the catalog name for this partitioner type
   *
   */
  string getCatalogName() const
  {
    // Extract type name from the data context
    return getDataContext().toString();
  }

};

} // namespace geos

#endif // GEOS_PARTITIONER_GRAPHPARTITIONENGINE_HPP_
