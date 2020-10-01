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
 * @file GraphEdge.cpp
 */
#include <string>
#include <iostream>
#include "common/Path.hpp"
#include "mesh/GraphEdge.hpp"
#include "mesh/GraphVertex.hpp"

namespace geosx
{


GraphEdge::GraphEdge( int const index, std::shared_ptr<GraphVertex> neighbour1, std::shared_ptr<GraphVertex> neighbour2, real64 transm):
  m_EdgeIndex(index),m_vertex1(neighbour1),m_vertex2(neighbour2),m_transmissibility(transm)
  {    
    //std::cout<<"Constructing edge "<<m_EdgeIndex<<" between vertices "<<m_vertex1->getIndice()<<" and "<<m_vertex2->getIndice()<<"\n";
  }


GraphEdge::~GraphEdge()
  {
    //std::cout<<"Destructing edge "<<m_EdgeIndex<<" between vertices "<<m_vertex1->getGlobalVertexIndex()<<" and "<<m_vertex2->getGlobalVertexIndex()<<"\n";
  }

} /* namespace geosx */
