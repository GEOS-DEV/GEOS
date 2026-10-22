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
 * @file GraphToolsMPI.hpp
 */

#ifndef GEOS_GRAPH_GRAPHTOOLSMPI_HPP_
#define GEOS_GRAPH_GRAPHTOOLSMPI_HPP_

#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"
#include <vector>


namespace geos
{
namespace graph
{

/**
 * @brief Distributes graph data and returns local adjacency lists for each rank.
 *
 * This function scatters the graph data (xadj and adjncy) across multiple MPI ranks, assuming one node per rank.
 * It ensures that the number of ranks matches the size of xadj-1 and performs the necessary
 * checks for graph validity. Each rank receives a portion of the adjacency list.
 *
 * @param xadj The adjacency list offsets for each node.
 * @param adjncy The adjacency list containing the neighbors of each node.
 * @param comm The MPI communicator (default is MPI_COMM_GEOS).
 * @return A pair of vectors containing the local adjacency list offsets and neighbors for each rank.
 */
std::pair< stdVector< size_t >, stdVector< size_t > >
scatterGraphData( const stdVector< size_t > & xadj,
                  const stdVector< size_t > & adjncy,
                  MPI_Comm comm = MPI_COMM_GEOS );


/**
 * @brief Gathers distributed graph data from MPI ranks.
 *
 * This function collects local adjacency lists from all ranks and reconstructs the global xadj and adjncy arrays.
 *
 * @param localXadj Local adjacency list offsets.
 * @param localAdjncy Local adjacency list.
 * @param comm MPI communicator (default: MPI_COMM_GEOS).
 * @return A pair of vectors: global xadj and adjncy.
 */

std::pair< stdVector< size_t >, stdVector< size_t > >
gatherGraphData( const stdVector< size_t > & localXadj,
                 const stdVector< size_t > & localAdjncy,
                 MPI_Comm comm= MPI_COMM_GEOS );


/**
 * @brief Creates the xadj array from the local adjacency array.
 *
 * This function takes the local adjacency array and the MPI communicator,
 * and creates the xadj array which represents the adjacency list offsets for each node.
 *
 * @param localAdjncy The local adjacency array.
 * @param comm The MPI communicator.
 * @return The xadj array.
 */
stdVector< size_t > createXadjFromAdjncy( const stdVector< size_t > & localAdjncy, MPI_Comm comm );


/**
 * @brief Generates global vertex IDs for distributed graph nodes.
 *
 * This function assigns a unique global ID to each vertex in the distributed graph.
 *
 * @param localXadj Local adjacency list offsets.
 * @param comm MPI communicator.
 * @return A vector of global vertex IDs.
 */
stdVector< int > createVertexGlobalID( const stdVector< size_t > & localXadj, MPI_Comm comm );

} // namespace geos
} // namespace graph

#endif // GEOS_GRAPH_GRAPHTOOLSMPI_HPP_
