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
 * @file GraphEdge2.cpp
 */
#include <string>
#include <iostream>
#include "common/Path.hpp"
#include "mesh/GraphEdge2.hpp"
#include "mesh/GraphVertex.hpp"
#include "mesh/GraphVertexFace.hpp"


namespace geosx
{


GraphEdge2::GraphEdge2( int const index, std::shared_ptr<GraphVertex> neighbour1, std::shared_ptr<GraphVertexFace> neighbour2, real64 transm):
  m_EdgeIndex(index),m_vertex1(neighbour1),m_vertex2(neighbour2),m_transmissibility(transm)
  {    
    //std::cout<<"Constructing edge "<<m_EdgeIndex<<" between vertices "<<m_vertex1->getIndice()<<" and "<<m_vertex2->getIndice()<<"\n";
  }


GraphEdge2::~GraphEdge2()
  {
    //std::cout<<"Destructing edge "<<m_EdgeIndex<<" between vertices "<<m_vertex1->getIndice()<<" and "<<m_vertex2->getIndice()<<"\n";
  }

} /* namespace geosx */
