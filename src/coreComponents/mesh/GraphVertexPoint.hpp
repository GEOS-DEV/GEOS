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
 * @file GraphVertexPoint.hpp
 */

#ifndef GEOSX_MESH_GRAPHVERTEXPOINT_HPP_
#define GEOSX_MESH_GRAPHVERTEXPOINT_HPP_
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mesh/GraphVertex.hpp"
namespace geosx
{

/**
 * @class GraphVertexPoint
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class GraphVertexPoint : public GraphVertex
{
public:

  /**
  * Constructor for GraphVertex object
  * @param [in] index of the vertex
  *
  */ 
  GraphVertexPoint( const int regionInd, const int subRegionInd, const int vertexInd);

  localIndex getCorrespondingId() const { return m_correspondingId; }
  
  void setCorrespondingId (int correspondingId) { m_correspondingId = correspondingId; }
    
private:
  localIndex m_correspondingId;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_GRAPHVERTEXPOINT_HPP_ */
