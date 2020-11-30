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

  GraphEdge()=delete;

  /**
   * @brief Constructor for GraphEdge object
   * @param [in] index of the edge
   * @param [in] neighbour1 one of the two GraphVertex forming the connection
   * @param [in] neighbour2 the other GraphVertex
   * 
   */   
  GraphEdge( int const index, std::shared_ptr<GraphVertex> neighbour1, std::shared_ptr<GraphVertex> neighbour2, real64 transm);

  /// Destructor
  virtual ~GraphEdge();

  localIndex getEdgeIndex() const { return m_EdgeIndex; }

  localIndex getEdgeIndex() { return m_EdgeIndex; }

  std::shared_ptr<GraphVertex> getVertex1() const { return m_vertex1; }

  std::shared_ptr<GraphVertex> getVertex1() { return m_vertex1; }

  std::shared_ptr<GraphVertex> getVertex2() const { return m_vertex2; }

  std::shared_ptr<GraphVertex> getVertex2() { return m_vertex2; }

  real64 getTransmissibility() const { return m_transmissibility; }

  real64 getTransmissibility() { return m_transmissibility; }

  
private:
  localIndex m_EdgeIndex;
  std::shared_ptr<GraphVertex> m_vertex1;
  std::shared_ptr<GraphVertex> m_vertex2;
  real64 m_transmissibility;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_GRAPHEDGE_HPP_ */
