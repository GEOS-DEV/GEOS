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
 * @file Edge.hpp
 */

#ifndef GEOSX_MESH_EDGE_HPP_
#define GEOSX_MESH_EDGE_HPP_
#include "mesh/Vertice.hpp"

namespace geosx
{

/**
 * @class Edge
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class Edge
{
public:

  Edge( const int, Vertice*, Vertice*);

  /// Destructor
  virtual ~Edge();

  int getIndice() const { return ind; }


  
private:
  int ind;
  Vertice* n1;
  Vertice* n2;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_EDGE_HPP_ */
