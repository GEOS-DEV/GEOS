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
 * @file GraphVertex.hpp
 */

#ifndef GEOSX_MESH_GRAPHVERTEX_HPP_
#define GEOSX_MESH_GRAPHVERTEX_HPP_
#include "meshUtilities/ComputationalGeometry.hpp"
namespace geosx
{

/**
 * @class GraphVertex
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class GraphVertex
{
public:
  /**
  * Constructor for GraphVertex object
  * @param [in] index of the vertex
  *
  */ 
  GraphVertex( const int regionInd, const int subRegionInd, const int vertexInd);

  /// Destructor
  virtual ~GraphVertex();

  globalIndex getGlobalVertexIndex() const { return m_globalVertexIndex; }

  globalIndex getGlobalVertexIndex() { return m_globalVertexIndex; }

  localIndex getLocalVertexIndex() const { return m_localVertexIndex; }

  localIndex getLocalVertexIndex() { return m_localVertexIndex; }

  localIndex getRegionIndex() const { return m_regionIndex; }

  localIndex getRegionIndex() { return m_regionIndex; }

  localIndex getSubRegionIndex() const { return m_subRegionIndex; }

  localIndex getSubRegionIndex() { return m_subRegionIndex; }

  localIndex getGhostIndex() const { return m_ghostIndex; }

  localIndex getGhostIndex() { return m_ghostIndex; }

  localIndex getCorrespondingId() const {return 0;}

  void setGhostIndex (int ghostIndex) { m_ghostIndex = ghostIndex; }

  void setLocalIndex (int elementIndex) { m_localVertexIndex = elementIndex; }
  
  void setCorrespondingId (int correspondingId) {std::cout<<"Could not write " << correspondingId << "\n";}

private:
  localIndex m_regionIndex;
  localIndex m_subRegionIndex;
  localIndex m_localVertexIndex;
  globalIndex m_globalVertexIndex;
  localIndex m_ghostIndex;
  };

} /* namespace geosx */

#endif /* GEOSX_MESH_GRAPHVERTEX_HPP_ */
