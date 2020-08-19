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
  GraphVertex( const int);

  /// Destructor
  virtual ~GraphVertex();

   localIndex getIndice() const { return ind; }

  
private:
  localIndex ind;
  };

} /* namespace geosx */

#endif /* GEOSX_MESH_GRAPHVERTEX_HPP_ */
