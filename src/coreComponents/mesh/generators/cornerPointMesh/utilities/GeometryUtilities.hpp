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

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTMESH_UTILITIES_GEOMETRYUTILITIES_HPP_
#define GEOSX_MESHUTILITIES_CORNERPOINTMESH_UTILITIES_GEOMETRYUTILITIES_HPP_

namespace geosx
{

namespace cornerPointMesh
{

namespace geometryUtilities
{

/**
 * @struct Vertex
 * @brief a helper struct used to filter out duplicate nodes
 */
struct Vertex
{
public:

  /**
   * @brief Constructor
   * @param[in] x X coordinate of the vertex
   * @param[in] y Y coordinate of the vertex
   * @param[in] z Z coordinate of the vertex
   */
  Vertex( real64 const & x, real64 const & y, real64 const & z )
    : m_localIndex( -1 ), m_globalIndex( -1 ), m_x( x ), m_y( y ), m_z( z ) {}

  /// Local index of the vertex
  localIndex m_localIndex;

  /// Global index of the vertex
  globalIndex m_globalIndex;

  /// X coordinate of the vertex
  real64 m_x;

  /// Y coordinate of the vertex
  real64 m_y;

  /// Z coordinate of the vertex
  real64 m_z;
};


/**
 * @struct CompareVertices
 * @brief another helper struct used to filter out duplicate nodes
 */
struct CompareVertices
{

  /**
   * @brief Comparaison operator
   * @param[in] v1 first vertex
   * @param[in] v2 second vertex
   * @return true if v1 > v2, and false otherwise
   */
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

struct Edge
{

  Edge( const Vertex point0, const Vertex point1 )
  : m_point0( point0 ), m_point1( point1) {}

  /// Both ends of edge
  Vertex m_point0;
  Vertex m_point1;

  /**
   * @brief Comparison operator for edges. Every edge has two vertexes
   * @param[in] edge1, first array of two vertexes
   * @param[in] edge2, second array of two vertexes
   * @return true if edge1 and edge2 are possible to intersect, and false otherwise
   */
  bool compare( const Edge newEdge ) const
  {
    // this only compares if current edge is higher than newEdge, which returns true if it is and false otherwise
    real64 minZEdge1 = (m_point0.m_z <= m_point1.m_z) ? m_point0.m_z : m_point1.m_z ;
    real64 maxZEdge2 = (newEdge.m_point0.m_z >= newEdge.m_point1.m_z) ? newEdge.m_point0.m_z : newEdge.m_point1.m_z ;

    if ( minZEdge1 > maxZEdge2)
      return true;
    else
      return false;
  }
};

} // end namespace geometryUtilities

} // end namespace cornerPointMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CORNERPOINTMESH_GEOMETRYUTILITIES_HPP_
