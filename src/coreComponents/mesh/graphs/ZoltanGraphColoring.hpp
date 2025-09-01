/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License: LGPL-2.1-only
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
 * @file ZoltanGraphColoring.hpp
 */

#ifndef GEOS_GRAPH_ZOLTANGRAPHCOLORING_HPP_
#define GEOS_GRAPH_ZOLTANGRAPHCOLORING_HPP_

#include "GraphColoringBase.hpp"
#include "common/MpiWrapper.hpp"
#include <zoltan.h>
#include <zoltan_cpp.h>
#include <vector>
#include <string>

namespace geos
{
namespace graph
{


/**
 * @class ZoltanGraphColoring
 * @brief Graph coloring implementation using the Zoltan library.
 *
 * This class uses Zoltan's distributed graph coloring capabilities to assign
 * colors to graph vertices in parallel. It supports both full adjacency-based
 * coloring and simplified one-node-per-rank coloring.
 */
class ZoltanGraphColoring : public GraphColoringBase
{
public:

  /**
   * @brief Constructor.
   * @param comm MPI communicator (default: MPI_COMM_GEOS).
   */
  ZoltanGraphColoring( MPI_Comm comm = MPI_COMM_GEOS );

  /**
   * @brief Destructor.
   */
  ~ZoltanGraphColoring();

/**
 * @brief Returns the number of distinct colors used.
 * @param colors Vector of color assignments.
 * @return Number of unique colors.
 */
  size_t getNumberOfColors( const std::vector< int > & colors ) const;

/**
 * @brief Returns the number of colors assuming one node per rank.
 * @param color Color value.
 * @return Number of colors.
 */
  size_t getNumberOfColors( const int color ) const;

/**
 * @brief Validates the coloring of a graph assuming one node per rank.
 * @param adjncy Adjacency list.
 * @param color Color of the node.
 * @return True if the coloring is valid, false otherwise.
 */
  bool isColoringValid( const std::vector< camp::idx_t > & adjncy, const int color ) const;

  /**
   * @brief Colors a graph.
   * @param xadj Adjacency list offsets.
   * @param adjncy Adjacency list.
   * @return A vector of assigned colors.
   */
  std::vector< int > colorGraph( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy ) override;

  /**
   * @brief Simplified coloring assuming one node per rank.
   * @param adjncy Local adjacency list.
   * @return Number of colors used.
   */
  int colorGraph( const std::vector< camp::idx_t > & adjncy ) override;


private:
  Zoltan * m_zz; /**< Pointer to Zoltan instance. */

  struct ZoltanGraph
  {
    int m_numVertices;
    std::vector< int > m_vertexGID;
    std::vector< int > m_xadj;
    std::vector< int > m_adjncy;
    int m_rank;
  };

  /**
   * @brief Zoltan callback to get number of vertices.
   */
  static int getNumberOfVertices( void * data, int * ierr );

  /**
   * @brief Zoltan callback to get vertex list.
   */
  static void getVertexList( void * data, int sizeGID, int sizeLID,
                             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                             int wgtDim, float * objWgts, int * ierr );

  /**
   * @brief Zoltan callback to get number of edges per vertex.
   */
  static void getNumEdgesList( void * data, int sizeGID, int sizeLID,
                               int numObj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int * numEdges, int * ierr );

  /**
   * @brief Zoltan callback to get edge list.
   */
  static void getEdgeList( void * data, int sizeGID, int sizeLID,
                           int numObj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int * numEdges, ZOLTAN_ID_PTR adjncy, int * nborProc,
                           int wgtDim, float * ewgts, int * ierr );
};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_ZOLTANGRAPHCOLORING_HPP_
