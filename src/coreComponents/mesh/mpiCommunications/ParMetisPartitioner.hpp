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
 * @file ParMetisPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_PARMETISPARTITIONER_HPP_
#define GEOS_PARTITIONER_PARMETISPARTITIONER_HPP_

#include "GraphPartitioner.hpp"

namespace geos
{


/**
 * @class ParMetisPartitioner
 * @brief A partitioner that uses the ParMETIS library to perform graph-based partitioning.
 */
class ParMetisPartitioner : public GraphPartitioner
{
public:
  ParMetisPartitioner( string const & name,
                       Group * const parent );

  /**
   * @brief Structure to hold scoped key names for XML input.
   */
  struct viewKeyStruct
  {
    constexpr static char const * numRefinementsString() { return "numRefinements"; }
  };

  /**
   * @brief Return the name of the partitioner for the object Catalog factory.
   * @return The string "ParMetis".
   */
  static string catalogName() { return "ParMetis"; }

  /**
   * @brief Partitions a mesh using the ParMETIS library.
   * @param graph The input graph (edges of locally owned nodes).
   * @param vertDist The parallel distribution of vertices (vertex index offset on each rank).
   * @param numParts The target number of partitions.
   * @param comm The MPI communicator of processes to partition over.
   * @param numRefinements The number of partition refinement iterations.
   * @return An array of target partitions for each element in the local mesh.
   */
  array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                   arrayView1d< pmet_idx_t const > const & vertDist,
                                   pmet_idx_t const numParts,
                                   MPI_Comm comm,
                                   int const numRefinements ) override;

  /**
   * @brief Convert an element-node mesh map into a dual (element-element) graph.
   * @param elemToNodes The input mesh represented by its elem-node map.
   * @param elemDist The parallel distribution of elements.
   * @param comm The MPI communicator.
   * @param minCommonNodes Minimum number of shared nodes to create a graph edge.
   * @return A graph with an edge for every pair of elements that share at least minCommonNodes nodes.
   */
  static ArrayOfArrays< pmet_idx_t, pmet_idx_t >
  meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
              arrayView1d< pmet_idx_t const > const & elemDist,
              MPI_Comm comm,
              int const minCommonNodes );

  int getNumRefinements() const override { return m_numRefinements; }

private:
  int m_numRefinements;
};

} // namespace geos

#endif // GEOS_PARTITIONER_PARMETISPARTITIONER_HPP_
