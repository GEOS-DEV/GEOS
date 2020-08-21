/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GraphEdge.hpp
 */

#ifndef GEOSX_MESH_GRAPHEDGE_HPP_
#define GEOSX_MESH_GRAPHEDGE_HPP_
#include "mesh/GraphVertex.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"


namespace geosx
{

/**
 * @class GraphEdge
 *
 * Represents connection between two points for the GraphBase.
 */
class GraphEdge
{
public:
  /**
   * @brief Constructor for GraphEdge object
   * @param [in] index of the edge
   * @param [in] neighbour1 one of the two GraphVertexs forming the connection
   * @param [in] neighbour2 the other GraphVertex
   * 
   */   
  GraphEdge( int const index, GraphVertex* neighbour1, GraphVertex* neighbour2, real64 transm);

  /// Destructor
  virtual ~GraphEdge();

  localIndex getIndice() const { return ind; }

  localIndex getIndice() { return ind; }

  GraphVertex* getN1() const { return n1; }

  GraphVertex* getN2() const { return n2; }

  
private:
  localIndex ind;
  GraphVertex* n1;
  GraphVertex* n2;
  real64 transmissibility;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_GRAPHEDGE_HPP_ */
