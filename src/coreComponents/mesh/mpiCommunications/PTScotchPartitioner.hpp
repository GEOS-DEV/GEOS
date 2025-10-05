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
 * @file PTScotchPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_PTSCOTCHPARTITIONER_HPP_
#define GEOS_PARTITIONER_PTSCOTCHPARTITIONER_HPP_

#include "GraphPartitioner.hpp"

namespace geos
{

/**
 * @class PTScotchPartitioner
 * @brief A partitioner that uses the PT-Scotch library to perform graph-based partitioning.
 */
class PTScotchPartitioner : public GraphPartitioner
{
public:

  PTScotchPartitioner( string const & name,
                       Group * const parent );

  /**
   * @brief Return the name of the partitioner for the object Catalog factory.
   * @return The string "PTScotch".
   */
  static string catalogName() { return "PTScotch"; }

  /**
   * @brief Partitions a mesh using the PT-Scotch library.
   *
   * @param graph The input graph (edges of locally owned nodes).
   * @param vertDist The parallel distribution of vertices (ignored by PT-Scotch, but part of the interface).
   * @param numParts The target number of partitions.
   * @param comm The MPI communicator of processes to partition over.
   * @param numRefinements The number of partition refinement iterations (ignored by PT-Scotch).
   * @return An array of target partitions for each element in the local mesh.
   */
  array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                   arrayView1d< pmet_idx_t const > const & vertDist,
                                   pmet_idx_t const numParts,
                                   MPI_Comm comm,
                                   int const numRefinements ) override;

  int getNumRefinements() const override { return 0; }
};

} // namespace geos


#endif // GEOS_PARTITIONER_PTSCOTCHPARTITIONER_HPP_
