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
 * @file GraphTools.hpp
 */

#ifndef GEOS_GRAPH_GRAPHTOOLS_HPP_
#define GEOS_GRAPH_GRAPHTOOLS_HPP_

#include "common/DataTypes.hpp"
#include <vector>
#include <unordered_set>


namespace geos
{
namespace graph
{

using camp::idx_t;

/**
 * @brief Validates the graph based on the given criteria.
 *
 * This function checks if the graph meets the following criteria:
 * 1. If a node is a neighbor to another node, the reverse must also be true.
 * 2. There are no duplicate nodes in the neighbor list.
 * 3. All neighbor indices are within the valid range of node indices.
 *
 * @param xadj The adjacency list offsets for each node.
 * @param adjncy The adjacency list containing the neighbors of each node.
 * @return True if the graph meets all criteria, false otherwise.
 */
bool isGraphValid( const std::vector< idx_t > & xadj, const std::vector< idx_t > & adjncy );

/**
 * @brief Calculates the degree of a node in the graph.
 *
 * This function returns the degree of the specified node, which is the number of edges connected to the node.
 *
 * @param node The index of the node whose degree is to be calculated.
 * @param xadj The adjacency list offsets for each node.
 * @return The degree of the node.
 */
size_t getGraphNodeDegree( idx_t node, const std::vector< idx_t > & xadj );

/**
 * @brief Retrieves the neighbors of a node in the graph.
 *
 * This function returns a set of neighbors for the specified node, based on the adjacency list representation of the graph.
 *
 * @param node The index of the node whose neighbors are to be retrieved.
 * @param xadj The adjacency list offsets for each node.
 * @param adjncy The adjacency list containing the neighbors of each node.
 * @return A set of indices representing the neighbors of the node.
 */
std::unordered_set< idx_t > getGraphNodeNeighbors( idx_t node, const std::vector< idx_t > & xadj, const std::vector< idx_t > & adjncy );



/**
 * @brief Generates a random graph.
 *
 * This function generates a random graph with the specified number of nodes and edges.
 *
 * @param num_nodes Number of nodes in the graph.
 * @param num_edges Number of edges in the graph.
 * @return A tuple containing xadj and adjncy.
 */
std::tuple< std::vector< idx_t >, std::vector< idx_t > > generateGraphRandom( size_t num_nodes, size_t num_edges );

/**
 * @brief Generates the adjacency list representation (xadj and adjncy) for a Cartesian domain decomposition in 3D.
 *
 * This function creates the adjacency list representation for a 3D Cartesian grid, where each node is connected
 * to its 6 neighboring nodes in the x, y, and z directions.
 *
 * @param nx Number of divisions along the x-axis.
 * @param ny Number of divisions along the y-axis.
 * @param nz Number of divisions along the z-axis.
 * @return A tuple containing xadj and adjncy.
 */
std::tuple< std::vector< idx_t >, std::vector< idx_t > > generateGraphCartPartitionning3D6( idx_t nx, idx_t ny, idx_t nz );

/**
 * @brief Generates the adjacency list representation (xadj and adjncy) for a Cartesian domain decomposition in 3D.
 *
 * This function creates the adjacency list representation for a 3D Cartesian grid, where each node is connected
 * to its 26 neighboring nodes in the x, y, and z directions, including edge and corner neighbors.
 *
 * @param nx Number of divisions along the x-axis.
 * @param ny Number of divisions along the y-axis.
 * @param nz Number of divisions along the z-axis.
 * @return A tuple containing xadj and adjncy.
 */
std::tuple< std::vector< idx_t >, std::vector< idx_t > > generateGraphCartPartitionning3D26( idx_t nx, idx_t ny, idx_t nz );


} // namespace geos
} // namespace graph

#endif // GEOS_GRAPH_GRAPHTOOLS_HPP_
