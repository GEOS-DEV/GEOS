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
 * @file GraphColoringBase.hpp
 */

#ifndef GEOS_GRAPH_GRAPHCOLORINGBASE_HPP_
#define GEOS_GRAPH_GRAPHCOLORINGBASE_HPP_

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"
#include <vector>


namespace geos
{
namespace graph
{

/**
 * @class GraphColoringBase
 * @brief Abstract base class for graph coloring strategies.
 *
 * Provides a common interface for graph coloring algorithms used in GEOS.
 * Derived classes must implement the coloring logic for adjacency graphs.
 */
class GraphColoringBase
{
public:

  /**
   * @brief Constructor.
   * @param comm MPI communicator (default: MPI_COMM_GEOS).
   */
  GraphColoringBase( MPI_Comm comm = MPI_COMM_GEOS ): m_comm( comm ) {}


  /**
   * @brief Pure virtual method to color a graph.
   * @param xadj Adjacency list offsets for each node.
   * @param adjncy Adjacency list containing neighbors of each node.
   * @return A vector of colors assigned to each node.
   */
  virtual std::vector< int > colorGraph( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy ) = 0;

  /**
   * @brief Pure virtual method to color a graph assuming one node per rank.
   * @param adjncy Adjacency list containing neighbors of each node.
   * @return Color of the node.
   */
  virtual int  colorGraph( const std::vector< camp::idx_t > & adjncy ) = 0;


  /**
   * @brief Virtual destructor.
   */
  virtual ~GraphColoringBase() {}



/**
 * @brief Checks the validity of the graph coloring.
 *
 * This function verifies that no two adjacent nodes in the graph have the same color.
 *
 * @param xadj The adjacency list offsets for each node.
 * @param adjncy The adjacency list containing the neighbors of each node.
 * @param coloring A vector where the index represents the node and the value represents the assigned color.*
 * @return True if the coloring is valid, false otherwise.
 */
  static bool isColoringValid( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy, const std::vector< int > & coloring );

  /**
   * @brief Checks the validity of the graph coloring assuming one node per rank.
   *
   * @param adjncy The adjacency list containing the neighbors of each node.
   * @param color Color assigned to the local node.
   * @param comm MPI communicator.
   *
   * @return True if the coloring is valid, false otherwise.
   */
  static bool isColoringValid( const std::vector< camp::idx_t > & adjncy, const int color, MPI_Comm comm );

/**
 * @brief Counts the number of distinct colors.
 *
 * This function takes a vector of integers representing colors
 * and returns the number of distinct (positive) colors present in the vector.
 *
 * @param colors A vector of integers representing colors.
 * @return The number of distinct colors in the vector.
 */
  static size_t getNumberOfColors( const std::vector< int > & colors );

  /**
   * @brief Counts the number of distinct colors  (parallel version).
   * @param colors Vector of color values.
   * @param comm MPI communicator.
   * @return Number of distinct colors.
   */
  static size_t getNumberOfColors( const std::vector< int > & colors, MPI_Comm comm );

  /**
   * @brief Counts the number of distinct colors assuming one node per rank (parallel version).
   * @param color Color of the node.
   * @param comm MPI communicator.
   * @return Number of distinct colors.
   */
  static size_t getNumberOfColors( const int color, MPI_Comm comm );


protected:
  MPI_Comm m_comm;  /**< MPI communicator used for parallel operations. */

};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_GRAPHCOLORINGBASE_HPP_
