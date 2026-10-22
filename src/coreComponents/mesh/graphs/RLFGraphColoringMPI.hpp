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
 * @file RLFGraphColoringMPI.hpp
 */

#ifndef GEOS_GRAPH_RLFGRAPHCOLORINGMPI_HPP_
#define GEOS_GRAPH_RLFGRAPHCOLORINGMPI_HPP_

#include <vector>
#include <unordered_set>
#include <cstddef>
#include "common/DataTypes.hpp"
#include "GraphTools.hpp"
#include "GraphColoringBase.hpp"

namespace geos
{
namespace graph
{


/**
 * @class RLFGraphColoringMPI
 * @brief Distributed graph coloring using the Recursive Largest First (RLF) algorithm.
 *
 * This class implements a parallel version of the RLF algorithm for graph coloring
 * using MPI.
 */
class RLFGraphColoringMPI : public GraphColoringBase
{
public:

  /**
   * @brief Constructor.
   * @param comm MPI communicator (default: MPI_COMM_GEOS).
   */
  RLFGraphColoringMPI( MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Destructor.
   */
  ~RLFGraphColoringMPI();

  /**
   * @brief Returns the number of distinct colors used.
   * @param colors Vector of color assignments.
   * @return Number of unique colors.
   */
  size_t getNumberOfColors( const stdVector< int > & colors ) const;

  /**
   * @brief Returns the number of distinct colors used, assuming one node per rank.
   * @param color Color value.
   * @return Number of unique colors.
   */
  size_t getNumberOfColors( const int color ) const;

  /**
   * @brief Validates the coloring of a graph partition.
   * @param localAdjncy adjacency list.
   * @param localColor Color assigned to the local node.
   * @return True if the coloring is valid, false otherwise.
   */
  bool isColoringValid( const stdVector< size_t > & localAdjncy, const int localColor ) const;

  /**
   * @brief Colors a distributed graph.
   * @param localXadj Local adjacency list offsets.
   * @param localAdjncy Local adjacency list.
   * @return A vector of assigned colors.
   */
  stdVector< int > colorGraph( const stdVector< size_t > & localXadj, const stdVector< size_t > & localAdjncy ) override;

  /**
   * @brief Simplified coloring assuming one node per rank.
   * @param localAdjncy Local adjacency list.
   * @return Color of the node.
   */
  int colorGraph( const stdVector< size_t > & localAdjncy ) override;
};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_RLFGRAPHCOLORINGMPI_HPP_
