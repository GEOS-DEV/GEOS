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
///@cond DO_NOT_DOCUMENT
#ifndef GEOSX_MESHUTILITIES_CORNERPOINTMESH_UTILITIES_GEOMETRYUTILITIES_HPP_
#define GEOSX_MESHUTILITIES_CORNERPOINTMESH_UTILITIES_GEOMETRYUTILITIES_HPP_

#include "codingUtilities/Utilities.hpp"

namespace geosx
{

namespace cornerPointMesh
{

namespace geometryUtilities
{

inline real64 roundWithPrecision( real64 const & num, localIndex const & digit )
{
  // return 0 if rounded number with given decimal precision is smaller than 1
  return ( round( num * pow( 10, digit ))/ pow( 10, digit ));
}

static real64 SMALL = 1e-15;

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
  Vertex( real64 const & x,
          real64 const & y,
          real64 const & z )
    : m_localIndex( -1 ),
    m_globalIndex( -1 ),
    m_x( x ),
    m_y( y ),
    m_z( z ),
    m_status( true )
  {}

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

  /// activeness
  bool m_status;

  Vertex( )
  {
    m_x = SMALL;
    m_y = SMALL;
    m_z = SMALL;
    m_status = false;
    m_localIndex = -1;
    m_globalIndex = -1;
  }

  void setIndex( localIndex const index )
  {
    m_localIndex = index;
  }
  void makeInactive()
  {
    m_status = false;
  }

  bool isActive() const
  {
    return m_status;
  }

  Vertex findLowerPoint( Vertex const & newPoint ) const
  {
    return (m_z <= newPoint.m_z) ? *this : newPoint;
  }

  Vertex findHigherPoint( Vertex const & newPoint ) const
  {
    return (m_z >= newPoint.m_z) ? *this : newPoint;
  }

  bool operator==( Vertex const & newPoint ) const
  {
    return ( isZero( m_x - newPoint.m_x ) &&
             isZero( m_y - newPoint.m_y ) &&
             isZero( m_z - newPoint.m_z ) );
  }

  bool operator<( Vertex const & newPoint ) const
  {
    if( !isZero( m_x - newPoint.m_x ) )
    {
      return m_x > newPoint.m_x;
    }
    else if( !isZero( m_y - newPoint.m_y ) )
    {
      return m_y > newPoint.m_y;
    }
    else
    {
      return m_z > newPoint.m_z;
    }
  }

  Vertex operator-( Vertex const & newPoint ) const
  {
    return Vertex( m_x - newPoint.m_x, m_y - newPoint.m_y, m_z - newPoint.m_z );
  }

  Vertex operator+( Vertex const & newPoint ) const
  {
    return Vertex( m_x + newPoint.m_x, m_y + newPoint.m_y, m_z + newPoint.m_z );
  }

  real64 mag() const
  {
    return sqrt( pow( m_x, 2 ) + pow( m_y, 2 ) + pow( m_z, 2 ));
  }

  Vertex crossProduct( Vertex const & p2 ) const
  {
    return Vertex( m_y*p2.m_z - m_z*p2.m_y, m_z*p2.m_x - m_x*p2.m_z, m_x*p2.m_y - m_y*p2.m_x );
  }

  real64 dotProduct( Vertex const & p2 ) const
  {
    return (m_x*p2.m_x + m_y*p2.m_y + m_z*p2.m_z);
  }

  Vertex scale( const real64 & f ) const
  {
    return Vertex( m_x * f, m_y * f, m_z * f );
  }

  real64 calcDist( Vertex const & p1 ) const
  {
    return sqrt( pow( m_x - p1.m_x, 2 ) + pow( m_y - p1.m_y, 2 ) + pow( m_z - p1.m_z, 2 ));
  }
};

struct Line
{

  Line( Vertex const & point0,
        Vertex const & point1 )
    : m_point0( point0 ),
    m_point1( point1 )
  {}

  /// Both ends of Line
  Vertex m_point0;
  Vertex m_point1;

  /**
   * @brief Comparison operator for edges. Every edge has two vertexes
   * @param[in] Line1, first array of two vertexes
   * @param[in] Line2, second array of two vertexes
   * @return true if edge1 and edge2 are possible to intersect, and false otherwise
   */
  bool compare( Line const & newLine ) const
  {
    // this only compares if current edge is higher than newEdge, which returns true if it is and false otherwise
    real64 const minZEdge1 = (m_point0.m_z <= m_point1.m_z) ? m_point0.m_z : m_point1.m_z;
    real64 const maxZEdge2 = (newLine.m_point0.m_z >= newLine.m_point1.m_z) ? newLine.m_point0.m_z : newLine.m_point1.m_z;

    return ( minZEdge1 > maxZEdge2 );
  }

  Vertex normalize( Vertex const & p1, Vertex const & p2 ) const
  {
    Vertex const p = p1 - p2;
    real64 const m = p.mag();
    return ( fabs( m - 0 ) < SMALL ) ? Vertex( 0.0, 0.0, 0.0 ) : Vertex( p.m_x / m, p.m_y/m, p.m_z /m );
  }

  bool findIntersectionPoint( Line const & newLine, Vertex & interPoint ) const
  {
    // two points have to be two lines at the same time. When two points are presented, the first point is returned.
    // This is an arbitrary choice.
    // TODO: use the isoparametric mapping to find the intersection points
    bool isOnLine1( false );
    bool isOnLine2( false );
    bool isNewPoint( false );

    Vertex const p1 = m_point0;
    Vertex const p2 = m_point1;
    Vertex const p3 = newLine.m_point0;
    Vertex const p4 = newLine.m_point1;

    Vertex const A = p1 - p3;
    Vertex const B = p2 - p1;
    Vertex const C = p4 - p3;

    // uv1 and uv2 are two unit vectors
    Vertex const uv1 = this->normalize( p1, p2 );
    Vertex const uv2 = this->normalize( p3, p4 );
    // Check for parallel lines

    Vertex const cp12 = uv1.crossProduct( uv2 );
    real64 const cp12Mag = cp12.mag();
    //real64 tol(1e-5);

    if( !isZero( roundWithPrecision( cp12Mag, 6 ) - 0.0 ))
    {
      real64 const ma = (A.dotProduct( C ) * C.dotProduct( B ) - A.dotProduct( B ) * C.dotProduct( C ))
                        / (B.dotProduct( B ) * C.dotProduct( C ) - C.dotProduct( B ) * C.dotProduct( B ));
      real64 const mb = ((ma * C.dotProduct( B )) + A.dotProduct( C ))/ C.dotProduct( C );

      // Calculate the point on line1 that is the closet point to line2
      Vertex Pa = p1 + B.scale( ma );
      // Calculate the point on line2 that is the closet point to line1
      Vertex const Pb = p3 + C.scale( mb );

      real64 const intersDist = Pa.calcDist( Pb );
      if( roundWithPrecision( ma, 3 ) >= 0.0 && roundWithPrecision( ma, 3 ) <= 1.0 )
      {
        isOnLine1 = true;
        if( isZero( roundWithPrecision( ma, 2 ) - 0.0 ))
        {
          isNewPoint = false;
          interPoint = p1;
        }
        else if( isZero( roundWithPrecision( ma, 2 )- 1.0 ))
        {
          isNewPoint = false;
          interPoint = p2;
        }
        else
        {
          isNewPoint = true;
        }
      }
      else
      {
        // not on line
        isOnLine1 = false;
      }

      isOnLine2 = (roundWithPrecision( mb, 3 ) >= 0.0 && roundWithPrecision( mb, 3 ) <= 1.0) ? true : false;

      if( roundWithPrecision( intersDist, 4 ) > 0.0 )
      {
        return false;
      }
      else
      {
        // arbitrarly return pa
        if( isOnLine1 && isOnLine2 )
        {
          if( isNewPoint )
          {
            Pa.setIndex( -1 );
            interPoint = Pa;
          }
        }
        else
        {
          return false;
        }
      }
    }
    else
    {
      return false;
    }

    return true;
  }
};

struct Face
{
  /*
   * A1 ---- A2
   * |       |
   * |       |
   * A3 ---- A4
   */

  Face( Vertex const & point0,
        Vertex const & point1,
        Vertex const & point2,
        Vertex const & point3 )
    : m_point0( point0 ),
    m_point1( point1 ),
    m_point2( point2 ),
    m_point3( point3 )
  {}

  Vertex m_point0;
  Vertex m_point1;
  Vertex m_point2;
  Vertex m_point3;

  void setIndexForFacePoints( localIndex const (&faceVertices)[ 4 ] )
  {
    m_point0.setIndex( faceVertices[0] );
    m_point1.setIndex( faceVertices[1] );
    m_point2.setIndex( faceVertices[2] );
    m_point3.setIndex( faceVertices[3] );
  }

  std::vector< Vertex > findIntersectionPoints( Face const & faceB, bool & isValid ) const
  {
    //Find intersection points between two faces. face0 and face1 are not coplanar, but both of them share the same pillars.
    //Here we set face0 as the targeted face, which contains A1 through A4, and face1 as the matching face, which contains
    //B1 through B4. The position of points belonging to face0 is as follow:
    //z
    //|__ x
    //       A1 -------- A2
    //       |           |
    //       |           |
    //       A3 -------- A4
    //
    //       3 -------- 2
    //       |          |
    //       |          |
    //       0 -------- 1
    //if two faces do not intersect, return null list. Otherwise, return a list of point coordinates, which include points from input
    //and new intersected points if any. For calculating those points and returning point array in a proper sequence, we follow the
    //following procedure:
    //(1) find the 4 points lying on the pillars are the lowest of the top lines
    //   and the highest of the bottom lines:
    //   p1 = argmin(A1, B1)
    //   p3 = argmin(A2, B2)
    //   p5 = argmax(A3, B3)
    //   p7 = argmax(A4, B4)
    //
    //(2) find the four points in between, which are the four possible intersections of the
    //   lines A12, B12, A34 and B34.
    //   p2 = A12 X B12
    //   p6 = A34 X B34
    //
    //   if p1 < p5, we need to find p8, which could be A12 X B34 or B12 X A34
    //   if p3 < p7, we need to find p4, which could be either A34 X B12 or A12 X B34
    //
    //(3) find intersection between line segments
    //
    //(4) return a list of coordinates in a designated order, clockwise, which points towards face0
    //p_array takes p0, p1, p2, p3, p4, p5, p6, p7
    //return sequence should be p0, p1, p2, p3, p6, p5, p4, p7

    // TODO: for now, we return a list of calculated points without giving their indices.
    // Among these points, some are new while some are not.
    // Such an information will be needed when updating the coordinates of all vertices.
    // For now, we do such a calculation after the new point is computed.

    // return a list points

    isValid = true;
    localIndex const numVertices = 8;
    std::vector< Vertex > newPoints;
    std::vector< Vertex > faceVertices( numVertices, Vertex() );

    // get temporary points: p1, p3, p5, p7
    faceVertices[0] = m_point3.findLowerPoint( faceB.m_point3 );
    faceVertices[2] = m_point2.findLowerPoint( faceB.m_point2 );

    faceVertices[4] = m_point0.findHigherPoint( faceB.m_point0 );
    faceVertices[6] = m_point1.findHigherPoint( faceB.m_point1 );

    Vertex p2( 0, 0, 0 );
    Vertex p6( 0, 0, 0 );
    // try to find p2
    Line const lineA12 = Line( m_point3, m_point2 );
    Line const lineB12 = Line( faceB.m_point3, faceB.m_point2 );
    bool const foundP2 = lineA12.findIntersectionPoint( lineB12, p2 );
    if( foundP2 )
    {
      faceVertices[1] = p2;
    }
    else
    {
      faceVertices[1].makeInactive();
    }
    // not a point or a line segement

    // try to find p6
    Line const lineA34 = Line( m_point0, m_point1 );
    Line const lineB34 = Line( faceB.m_point0, faceB.m_point1 );
    bool const foundP6 = lineA34.findIntersectionPoint( lineB34, p6 );
    if( foundP6 )
    {
      faceVertices[5] = p6;
    }
    else
    {
      faceVertices[5].makeInactive();
    }

    if( faceVertices[0].m_z < faceVertices[4].m_z )
    {
      // try to find p8
      Vertex p8( 0, 0, 0 );
      bool const foundP8 =
        ( isZero( faceVertices[0].m_z - faceB.m_point3.m_z ) )
  ? lineB12.findIntersectionPoint( lineA34, p8 )
  : lineB34.findIntersectionPoint( lineA12, p8 );
      if( foundP8 )
      {
        faceVertices[7] = p8;
      }
      else
      {
        isValid = false;
        GEOSX_THROW( "P8 cannot be found!!", InputError );
      }

      faceVertices[0].makeInactive();
      faceVertices[4].makeInactive();
    }

    if( faceVertices[2].m_z < faceVertices[6].m_z )
    {
      // try to find p4
      Vertex p4( 0, 0, 0 );
      bool const foundP4 =
        ( isZero( faceVertices[2].m_z - faceB.m_point2.m_z ) )
  ? lineA34.findIntersectionPoint( lineB12, p4 )
  : lineA12.findIntersectionPoint( lineB34, p4 );

      if( foundP4 )
      {
        faceVertices[3] = p4;
      }
      else
      {
        isValid = false;
        GEOSX_THROW( "P4 cannot be found!", InputError );
      }

      faceVertices[2].makeInactive();
      faceVertices[6].makeInactive();
    }

    // swap 6th and 4th element
    Vertex const tmpVertex = faceVertices[4];
    faceVertices[4] = faceVertices[6];
    faceVertices[6] = tmpVertex;

    for( localIndex i = 0; i < numVertices; ++i )
    {
      if( faceVertices[i].isActive() )
      {
        // valid vertex
        // intersection points found

        if( std::find( newPoints.begin(), newPoints.end(), faceVertices[i] ) ==  newPoints.end() )
        {
          newPoints.push_back( faceVertices[i] );
        }
      }
    }
    return newPoints;
  }
};

} // end namespace geometryUtilities

} // end namespace cornerPointMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CORNERPOINTMESH_GEOMETRYUTILITIES_HPP_
///@endcond DO_NOT_DOCUMENT
