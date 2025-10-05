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
 * @file NoOpPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_NOOPPARTITIONER_HPP_
#define GEOS_PARTITIONER_NOOPPARTITIONER_HPP_

#include "GraphPartitioner.hpp"

namespace geos
{

/**
 * @class NoOpPartitioner
 * @brief A partitioner that performs no redistribution, assigning all local elements to their current rank.
 *
 * This partitioner is useful for serial runs, debugging, or scenarios where an initial
 * data distribution is guaranteed to be final. It fulfills the GraphPartitioner contract.
 */
class NoOpPartitioner : public GraphPartitioner
{
public:
  /**
   * @brief Constructor
   * @param name Name of the partitioner instance.
   * @param parent Parent group in the data repository.
   */
  NoOpPartitioner( string const & name, Group * const parent );

  /**
   * @brief Catalog name for factory registration.
   * @return The string "NoOp".
   */
  static string catalogName() { return "NoOp"; }

  /**
   * @brief The "no-op" partition method.
   *
   * This implementation does not redistribute any data. It returns a partition array
   * where every local element is assigned to the current MPI rank.
   *
   * @param graph Input graph (mostly ignored).
   * @param vertDist Vertex distribution, used to find the number of local elements.
   * @param numParts Number of partitions (ignored).
   * @param comm MPI communicator.
   * @param numRefinements Number of refinements (ignored).
   * @return An array where each entry is the ID of the current rank.
   */
  array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                   arrayView1d< pmet_idx_t const > const & vertDist,
                                   pmet_idx_t const numParts,
                                   MPI_Comm comm,
                                   int const numRefinements ) override;

  /**
   * @brief Gets the number of refinement steps configured for this partitioner.
   * @return The number of refinement steps.
   */
  int getNumRefinements() const override { return 0; }
};

} // namespace geos

#endif // GEOS_PARTITIONER_NOOPPARTITIONER_HPP_
