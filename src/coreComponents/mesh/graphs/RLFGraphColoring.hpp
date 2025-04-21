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

#ifndef GEOS_GRAPH_RLFGraphColoring_HPP_
#define GEOS_GRAPH_RLFGraphColoring_HPP_

#include <vector>
#include "common/DataTypes.hpp"
#include "GraphColoringBase.hpp"

namespace geos
{
namespace graph
{

class RLFGraphColoring : public GraphColoringBase
{

public:
  RLFGraphColoring();
  ~RLFGraphColoring();

  std::vector< int > colorGraph( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy ) override;
  // Simplified version assuming one node per rank
  int colorGraph( const std::vector< camp::idx_t > & localAdjncy ) override;

  size_t getNumberOfColors( const std::vector< int > & colors ) const;
  bool isColoringValid( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy, const std::vector< int > & colors ) const;


private:
/**
 * @brief Colors a graph using the Recursive Largest First (RLF) algorithm.
 *
 * This function takes the adjacency list representation of a graph and colors the graph
 * such that no two adjacent nodes have the same color. The algorithm used is the Recursive
 * Largest First (RLF) algorithm.
 *
 * Reference:
 * Leighton, F. (1979). "A graph coloring algorithm for large scheduling problems".
 * Journal of Research of the National Bureau of Standards. 84 (6): 489–503.
 *
 * @param xadj The adjacency list offsets for each node.
 * @param adjncy The adjacency list containing the neighbors of each node.
 * @return A vector where the index represents the node and the value represents the assigned color.
 */
  std::vector< int > RecursiveLargestFirstColoring( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy );
};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_RLFGraphColoring_HPP_
