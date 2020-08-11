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
 * @file Edge.cpp
 */
#include <string>
#include "common/Path.hpp"
#include "Edge.hpp"

namespace geosx
{



Edge::Edge( const int index, Vertice  neighbour1, Vertice neighbour2):
  ind(index),n1(neighbour1),n2(neighbour2)
  {    
  }


Edge::~Edge()
{}

} /* namespace geosx */
