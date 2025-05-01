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

using camp::idx_t;

/**
 * @brief Distributes graph data and returns local adjacency lists for each rank.
 *
 * This function scatters the graph data (xadj and adjncy) across multiple MPI ranks, assuming one node per rank.
 * It ensures that the number of ranks matches the size of xadj-1 and performs the necessary
 * checks for graph validity. Each rank receives a portion of the adjacency list.
 *
 * @param xadj The adjacency list offsets for each node.
 * @param adjncy The adjacency list containing the neighbors of each node.
 * @param comm The MPI communicator (default is MPI_COMM_WORLD).
 * @return A pair of vectors containing the local adjacency list offsets and neighbors for each rank.
 */
std::pair< std::vector< camp::idx_t >, std::vector< camp::idx_t > >
scatterGraphData( const std::vector< camp::idx_t > & xadj,
                  const std::vector< camp::idx_t > & adjncy,
                  MPI_Comm comm = MPI_COMM_GEOS );


std::pair< std::vector< camp::idx_t >, std::vector< camp::idx_t > >
gatherGraphData( const std::vector< camp::idx_t > & localXadj,
                 const std::vector< camp::idx_t > & localAdjncy,
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
std::vector< camp::idx_t > createXadjFromAdjncy( const std::vector< camp::idx_t > & localAdjncy, MPI_Comm comm );


std::vector< int > createVertexGlobalID( const std::vector< camp::idx_t > & localXadj, MPI_Comm comm );


} // namespace geos
} // namespace graph

#endif // GEOS_GRAPH_GRAPHTOOLSMPI_HPP_
