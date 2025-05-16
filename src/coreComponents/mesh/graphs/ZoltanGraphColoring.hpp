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

class ZoltanGraphColoring : public GraphColoringBase
{
public:
  ZoltanGraphColoring( MPI_Comm comm = MPI_COMM_GEOS );
  ~ZoltanGraphColoring();

  size_t getNumberOfColors( const std::vector< int > & colors ) const;
  size_t getNumberOfColors( const int color ) const;
  bool isColoringValid( const std::vector< camp::idx_t > & adjncy, const int color ) const;


  std::vector< int > colorGraph( const std::vector< camp::idx_t > & xadj, const std::vector< camp::idx_t > & adjncy ) override;

  // Simplified version assuming one node per rank
  int colorGraph( const std::vector< camp::idx_t > & adjncy ) override;


private:
  Zoltan * m_zz;

  struct ZoltanGraph
  {
    int m_numVertices;
    std::vector< int > m_vertexGID;
    std::vector< int > m_xadj;
    std::vector< int > m_adjncy;
    int m_rank;
  };

  static int getNumberOfVertices( void * data, int * ierr );
  static void getVertexList( void * data, int sizeGID, int sizeLID,
                             ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                             int wgtDim, float * objWgts, int * ierr );
  static void getNumEdgesList( void * data, int sizeGID, int sizeLID,
                               int numObj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                               int * numEdges, int * ierr );
  static void getEdgeList( void * data, int sizeGID, int sizeLID,
                           int numObj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                           int * numEdges, ZOLTAN_ID_PTR adjncy, int * nborProc,
                           int wgtDim, float * ewgts, int * ierr );
};

} // namespace graph
} // namespace geos

#endif // GEOS_GRAPH_ZOLTANGRAPHCOLORING_HPP_
