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
 * @file RLFGraphColoring.hpp
 */

#ifndef GEOS_GRAPH_RLFGRAPHCOLORING_HPP_
#define GEOS_GRAPH_RLFGRAPHCOLORING_HPP_

#include "GraphColoringBase.hpp"
#include "common/DataTypes.hpp"
#include <vector>

namespace geos
{
namespace graph
{


/**
 * @class RLFGraphColoring
 * @brief Graph coloring implementation using the Recursive Largest First (RLF) algorithm.
 *
 * This class implements the RLF algorithm to assign colors to graph nodes such that
 * no two adjacent nodes share the same color. It supports both full adjacency-based
 * coloring and simplified one-node-per-rank coloring.
 *
 * Reference:
 * Leighton, F. (1979). "A graph coloring algorithm for large scheduling problems".
 * Journal of Research of the National Bureau of Standards. 84 (6): 489–503.
 */
class RLFGraphColoring : public GraphColoringBase
{

public:

  /**
   * @brief Constructor.
   */
  RLFGraphColoring();

  /**
   * @brief Destructor.
   */
  ~RLFGraphColoring();

  /**
   * @brief Colors a graph.
   * @param xadj Adjacency list offsets.
   * @param adjncy Adjacency list.
   * @return A vector of assigned colors.
   */
  std::vector< int > colorGraph( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy ) override;

  /**
   * @brief Colors a graph assuming one node per rank.
   * @param localAdjncy Local adjacency list.
   * @return Color of the node.
   */
  int colorGraph( const std::vector< camp::idx_t > & localAdjncy ) override;

  /**
   * @brief Returns the number of distinct colors used.
   * @param colors Vector of color assignments.
   * @return Number of unique colors.
   */
  size_t getNumberOfColors( const std::vector< int > & colors ) const;

  /**
   * @brief Validates the coloring of a graph.
   * @param xadj Adjacency list offsets.
   * @param adjncy Adjacency list.
   * @param colors Vector of assigned colors.
   * @return True if coloring is valid, false otherwise.
   */
  bool isColoringValid( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy, const std::vector< int > & colors ) const;


private:
/**
 * @brief Colors a graph using the Recursive Largest First (RLF) algorithm.
 * @param xadj The adjacency list offsets for each node.
 * @param adjncy The adjacency list containing the neighbors of each node.
 * @return A vector where the index represents the node and the value represents the assigned color.
 */
  std::vector< int > RecursiveLargestFirstColoring( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy );
};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_RLFGRAPHCOLORING_HPP_
