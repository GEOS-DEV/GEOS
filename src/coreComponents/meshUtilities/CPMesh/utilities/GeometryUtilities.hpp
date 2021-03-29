/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GeometryUtilities.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CPMESH_UTILITIES_GEOMETRYUTILITIES_HPP_
#define GEOSX_MESHUTILITIES_CPMESH_UTILITIES_GEOMETRYUTILITIES_HPP_

namespace geosx
{

namespace CPMesh
{

namespace GeometryUtilities
{

struct Vertex
{
public:
  Vertex( real64 const & x, real64 const & y, real64 const & z )
    : m_localVertex( -1 ), m_x( x ), m_y( y ), m_z( z ) {}
  localIndex m_localVertex;
  real64 m_x;
  real64 m_y;
  real64 m_z;
};

struct CompareVertices
{
  bool operator()( const Vertex & v1, const Vertex & v2 ) const
  {
    if( !isZero( v1.m_x - v2.m_x ) )
    {
      return v1.m_x > v2.m_x;
    }
    else if( !isZero( v1.m_y - v2.m_y ) )
    {
      return v1.m_y > v2.m_y;
    }
    else
    {
      return v1.m_z > v2.m_z;
    }
  }
};

} // end namespace GeometryUtilities

} // end namespace CPMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CPMESH_GEOMETRYUTILITIES_HPP_
