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
 * @file GraphPartitioner.hpp
 */

#ifndef GEOS_PARTITIONER_GRAPHPARTITIONER_HPP_
#define GEOS_PARTITIONER_GRAPHPARTITIONER_HPP_

#include "PartitionerBase.hpp"
#include "common/DataTypes.hpp"

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
 * @class GraphPartitioner
 * @brief An abstract base class for partitioners that operate on connectivity graphs (e.g., ParMETIS, PTScotch).
 */
class GraphPartitioner : public PartitionerBase
{
public:
  /// Inherit the constructor from PartitionerBase.
  using PartitionerBase::PartitionerBase;

  virtual ~GraphPartitioner() = default;

  /**
   * @brief Partition a mesh according to its dual graph.
   *
   * It takes a distributed graph representation of the mesh and returns an array indicating
   * the target partition for each local element.
   *
   * @param graph The input graph (edges of locally owned nodes).
   * @param vertDist The parallel distribution of vertices (vertex index offset on each rank).
   * @param numParts The target number of partitions.
   * @param comm The MPI communicator of processes to partition over.
   * @param numRefinements The number of partition refinement iterations.
   * @return An array of target partitions for each element in the local mesh.
   */
  virtual array1d< pmet_idx_t > partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                           arrayView1d< pmet_idx_t const > const & vertDist,
                                           pmet_idx_t const numParts,
                                           MPI_Comm comm,
                                           int const numRefinements ) = 0;

  // Provide concrete implementations for methods inherited from PartitionerBase
  void processCommandLineOverrides( unsigned int xparCL,
                                    unsigned int yparCL,
                                    unsigned int zparCL ) override;

  void setPartitionCounts( unsigned int xPartitions,
                           unsigned int yPartitions,
                           unsigned int zPartitions ) override;

  void setNeighborsRank( const std::vector< int > & neighborsRank ) override;

  /**
   * @brief Gets the number of refinement steps configured for this partitioner.
   * @return The number of refinement steps.
   */
  virtual int getNumRefinements() const = 0;

protected:
  /**
   * @brief Implementation of the graph coloring algorithm.
   */
  void color() override;
};

} // namespace geos

#endif // GEOS_PARTITIONER_GRAPHPARTITIONER_HPP_
