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
 * @file GraphVertexFace.cpp
 */
#include <string>
#include <iostream>
#include "common/Path.hpp"
#include "GraphVertexFace.hpp"
#include "GraphVertex.hpp"
#include "mesh/GraphEdge.hpp"

namespace geosx
{

GraphVertexFace::GraphVertexFace( const int regionInd, const int subRegionInd, const int vertexInd):GraphVertex(regionInd, subRegionInd, vertexInd),m_correspondingId(-1)
  {
  }


} /* namespace geosx */
