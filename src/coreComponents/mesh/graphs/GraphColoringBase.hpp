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

//#include "common/MpiWrapper.hpp"
#include "common/DataTypes.hpp"
#include <vector>
#include <mpi.h>


namespace geos
{
namespace graph
{

class GraphColoringBase
{
public:
  GraphColoringBase( MPI_Comm comm = MPI_COMM_WORLD ): m_comm( comm ) {}

  virtual std::vector< int > colorGraph( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy ) = 0;
  virtual int  colorGraph( const std::vector< camp::idx_t > & adjncy ) = 0;

  virtual ~GraphColoringBase() {}   // Virtual destructor



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
 * @brief Counts the number of distinct colors in a vector.
 *
 * This function takes a vector of integers representing colors
 * and returns the number of distinct (positive) colors present in the vector.
 *
 * @param colors A vector of integers representing colors.
 * @return The number of distinct colors in the vector.
 */
  static size_t getNumberOfColors( const std::vector< int > & colors );



  static bool isColoringValid( const std::vector< camp::idx_t > & adjncy, const int color, MPI_Comm comm );
  static size_t getNumberOfColors( const std::vector< int > & colors, MPI_Comm comm );
  static size_t getNumberOfColors( const int color, MPI_Comm comm );


protected:
  MPI_Comm m_comm;   // MPI communicator

};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_GRAPHCOLORINGBASE_HPP_
