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


GraphEdge::GraphEdge( int const index, GraphVertex* neighbour1, GraphVertex* neighbour2, real64 transm):
  ind(index),n1(neighbour1),n2(neighbour2),transmissibility(transm)
  {    
    std::cout<<"Constructing edge "<<ind<<" between vertices "<<n1->getIndice()<<" and "<<n2->getIndice()<<"\n";
  }


GraphEdge::~GraphEdge()
  {
    //std::cout<<"Destructing edge "<<ind<<" between vertices "<<n1->getIndice()<<" and "<<n2->getIndice()<<"\n";
  }

} /* namespace geosx */
