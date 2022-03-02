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
#include "codingUtilities/Utilities.hpp"

namespace geosx
{

namespace cornerPointMesh
{

namespace geometryUtilities
{

inline real64 roundWithPrecision(const real64 num , const localIndex& digit )
{
  // return 0 if rounded number with given decimal precision is smaller than 1
  return ( round(num * pow(10, digit))/ pow(10, digit));
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

  Vertex( real64 const & x, real64 const & y, real64 const & z )
    : m_localIndex( -1 ),
      m_globalIndex( -1 ),
      m_x( x ), m_y( y ), m_z( z ), m_status(true) {}
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

  void setIndex(const localIndex& index)
  {
    m_localIndex = index;
  }
  void setStatus(const bool flag)
  {
    m_status = flag;
  }

  bool getStatus() const
  {
    return m_status;
  }

  Vertex findLowerPoint( const Vertex& newPoint ) const
  {
    return (m_z <= newPoint.m_z) ? *this : newPoint;
  }

  Vertex findHigherPoint( const Vertex& newPoint ) const
  {
    return (m_z >= newPoint.m_z) ? *this : newPoint;
  }

  bool operator==(const Vertex& newPoint) const
  {
    if ( isZero(m_x - newPoint.m_x) && isZero( m_y - newPoint.m_y) &&
        isZero( m_z - newPoint.m_z))
      return true;
    else
      return false;
  }

//  bool operator<( const Vertex& newPoint) const
//  {
//    // TODO: isZero( v1.m_x - v2.m_x ) use isZero
//    if ( !( isZero(m_x - newPoint.m_x) && isZero( m_y - newPoint.m_y) &&
//            isZero( m_z - newPoint.m_z)) )
//      return true;
//    else
//      return false;
//  }

  bool operator<( const Vertex& newPoint) const
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

  Vertex operator-(const Vertex& newPoint) const
  {
    return Vertex( m_x - newPoint.m_x, m_y - newPoint.m_y, m_z - newPoint.m_z );
  }

  Vertex operator+(const Vertex& newPoint) const
  {
    return Vertex( m_x + newPoint.m_x, m_y + newPoint.m_y, m_z + newPoint.m_z );
  }

  void print() const
  {
    std::cout << "(" << m_x << "," << m_y << "," << m_z << ")";
  }

  void printIndex() const
  {
    // Debugging purpose
    std::cout << m_localIndex;
  }
  real64 mag() const
  {
    return sqrt( pow( m_x, 2) + pow( m_y, 2) + pow(m_z, 2));
  }

  Vertex crossProduct(const Vertex& p2) const
  {
    return Vertex( m_y*p2.m_z - m_z*p2.m_y, m_z*p2.m_x - m_x*p2.m_z, m_x*p2.m_y - m_y*p2.m_x);
  }

  real64 dotProduct(const Vertex& p2) const
  {
    return (m_x*p2.m_x + m_y*p2.m_y + m_z*p2.m_z);
  }

  Vertex ptFactor(const real64& f) const
  {
    return Vertex(m_x * f, m_y * f, m_z * f );
  }

  real64 calcDist(const Vertex& p1)
  {
    return sqrt( pow(m_x - p1.m_x, 2) + pow(m_y - p1.m_y, 2) + pow(m_z - p1.m_z, 2));
  }
};

///**
// * @struct CompareVertices
// * @brief another helper struct used to filter out duplicate nodes
// */
//struct CompareVertices
//{
//
//  /**
//   * @brief Comparaison operator
//   * @param[in] v1 first vertex
//   * @param[in] v2 second vertex
//   * @return true if v1 > v2, and false otherwise
//   */
//  bool operator()( const Vertex & v1, const Vertex & v2 ) const
//  {
//    if( !isZero( v1.m_x - v2.m_x ) )
//    {
//      return v1.m_x > v2.m_x;
//    }
//    else if( !isZero( v1.m_y - v2.m_y ) )
//    {
//      return v1.m_y > v2.m_y;
//    }
//    else
//    {
//      return v1.m_z > v2.m_z;
//    }
//  }
//};

struct Line
{

  Line( const Vertex& point0, const Vertex& point1 )
  : m_point0( point0 ), m_point1( point1) {}

  /// Both ends of Line
  Vertex m_point0;
  Vertex m_point1;

  /**
   * @brief Comparison operator for edges. Every edge has two vertexes
   * @param[in] Line1, first array of two vertexes
   * @param[in] Line2, second array of two vertexes
   * @return true if edge1 and edge2 are possible to intersect, and false otherwise
   */
  bool compare( const Line& newLine ) const
  {
    // this only compares if current edge is higher than newEdge, which returns true if it is and false otherwise
    real64 minZEdge1 = (m_point0.m_z <= m_point1.m_z) ? m_point0.m_z : m_point1.m_z ;
    real64 maxZEdge2 = (newLine.m_point0.m_z >= newLine.m_point1.m_z) ? newLine.m_point0.m_z : newLine.m_point1.m_z;

    if ( minZEdge1 > maxZEdge2)
      return true;
    else
      return false;
  }

  Vertex normalize( Vertex& p1, Vertex& p2) const
  {
    Vertex p = p1 - p2;
    real64 m = p.mag();
    return ( fabs( m - 0) < SMALL ) ? Vertex(0.0, 0.0, 0.0) : Vertex(p.m_x / m, p.m_y/m, p.m_z /m);
  }

  bool findIntersectionPoint( const Line& newLine, Vertex& interPoint) const
  {
    // two points have to be two lines at the same time. When two points are presented, the first point is returned. This is an arbitrary choice.
    // TODO: use the isoparametric mapping to find the intersection points
    bool isOnLine1(false);
    bool isOnLine2(false);
    bool isNewPoint(false);

    Vertex p1 = m_point0;
    Vertex p2 = m_point1;
    Vertex p3 = newLine.m_point0;
    Vertex p4 = newLine.m_point1;

    Vertex A = p1 - p3;
    Vertex B = p2 - p1;
    Vertex C = p4 - p3;

    // uv1 and uv2 are two unit vectors
    Vertex uv1 = this->normalize( p1, p2) ;
    Vertex uv2 = this->normalize( p3, p4);
    // Check for parallel lines

    Vertex cp12 = uv1.crossProduct(uv2);
    real64 cp12Mag = cp12.mag();
    //real64 tol(1e-5);

    if ( !isZero( roundWithPrecision( cp12Mag, 6) - 0.0))
    {
      real64 ma = (A.dotProduct(C) * C.dotProduct(B) - A.dotProduct(B) * C.dotProduct(C))/
                  (B.dotProduct(B) * C.dotProduct(C) - C.dotProduct(B) * C.dotProduct(B));
      real64 mb = ((ma * C.dotProduct(B)) + A.dotProduct(C))/ C.dotProduct(C);

      // Calculate the point on line1 that is the closet point to line2
      Vertex Pa = p1 + B.ptFactor(ma);
      // Calculate the point on line2 that is the closet point to line1
      Vertex Pb = p3 + C.ptFactor(mb);

      real64 intersDist = Pa.calcDist( Pb );
      if (roundWithPrecision(ma, 3) >= 0.0 && roundWithPrecision(ma, 3) <= 1.0)
      {
        isOnLine1 = true;
        if ( isZero(roundWithPrecision(ma, 2) - 0.0))
        {
          isNewPoint = false;
          interPoint = p1;
        }
        else if ( isZero(roundWithPrecision(ma, 2)- 1.0))
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

      isOnLine2 = (roundWithPrecision(mb, 3) >= 0.0 && roundWithPrecision(mb, 3) <= 1.0) ? true : false;

      if ( roundWithPrecision(intersDist, 4) > 0.0)
      {
        return false;
      }
      else
      {
        // arbitrarly return pa
        if (isOnLine1 && isOnLine2)
        {
          if (isNewPoint)
          {
            Pa.setIndex(-1);
            interPoint = Pa;
          }
        }
        else
          return false;
      }
    }
    else
      return false;

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

  Face( const Vertex& point0, const Vertex& point1,
        const Vertex& point2, const Vertex& point3)
  : m_point0( point0 ), m_point1( point1),
    m_point2( point2 ), m_point3( point3){}

  Vertex m_point0;
  Vertex m_point1;
  Vertex m_point2;
  Vertex m_point3;

  void printPoints()
  {
    std::cout << "Point0: ";
    m_point0.print();

    std::cout << " Point1: ";
    m_point1.print();

    std::cout << " Point2: ";
    m_point2.print();

    std::cout << " Point3: ";
    m_point3.print();
  }

  void setIndexForFacePoints(localIndex const (&faceVertices)[ 4 ])
  {
    m_point0.setIndex(faceVertices[0]);
    m_point1.setIndex(faceVertices[1]);
    m_point2.setIndex(faceVertices[2]);
    m_point3.setIndex(faceVertices[3]);
  }

  std::vector<Vertex> findIntersectionPoints(const Face& faceB, bool& isValid ) const
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

    // TODO: for now, we return a list of calculated points without giving their indexes. Among these points, some are new while some are not.
    // such a information will be needed when updating the coordinates of all vertices. For now, we do such a calculation after the new point is
    // computed.

    // return a list points

    isValid = true;
    localIndex vertexNum = 8;
    // std::vector<Vertex> faceVertexes( vertexNum, Vertex()); // be careful about initializing this array
    std::vector<Vertex> newPointVector;
    std::vector<Vertex> faceVertexes(vertexNum, Vertex());
//    std::vector<int> pointIndex(vertexNum, -1);

    // get temporary points: p1, p3, p5, p7
    faceVertexes[0] = m_point3.findLowerPoint(faceB.m_point3);
    faceVertexes[2] = m_point2.findLowerPoint(faceB.m_point2);

    faceVertexes[4] = m_point0.findHigherPoint(faceB.m_point0);
    faceVertexes[6] = m_point1.findHigherPoint(faceB.m_point1);

    Vertex p2(0, 0, 0);
    Vertex p6(0, 0, 0);
    // try to find p2
    Line lineA12 = Line( m_point3, m_point2);
    Line lineB12 = Line( faceB.m_point3, faceB.m_point2);
    bool foundP2 = lineA12.findIntersectionPoint(lineB12, p2);
    if (foundP2)
      faceVertexes[1] = p2;
    else
      faceVertexes[1].setStatus(false);
      // not a point or a line segement

    // try to find p6
    Line lineA34 = Line( m_point0, m_point1);
    Line lineB34 = Line( faceB.m_point0, faceB.m_point1);
    bool foundP6 = lineA34.findIntersectionPoint(lineB34, p6);
    if (foundP6)
      faceVertexes[5] = p6;
    else
      faceVertexes[5].setStatus(false);

    if (faceVertexes[0].m_z < faceVertexes[4].m_z)
    {
      // try to find p8
      Vertex p8(0, 0, 0);
      bool foundP8 = ( isZero(faceVertexes[0].m_z - faceB.m_point3.m_z ) ) ? lineB12.findIntersectionPoint(lineA34, p8) :
                      lineB34.findIntersectionPoint(lineA12, p8);
      if (foundP8)
        faceVertexes[7] = p8;
      else
      {
        isValid = false;
        GEOSX_THROW( "P8 cannot be found!!", InputError );
      }

      faceVertexes[0].setStatus(false);
      faceVertexes[4].setStatus(false);

    }

    if (faceVertexes[2].m_z < faceVertexes[6].m_z)
    {
      // try to find p4
      Vertex p4(0, 0, 0);
      bool foundP4 = ( isZero( faceVertexes[2].m_z - faceB.m_point2.m_z)) ? lineA34.findIntersectionPoint(lineB12, p4) :
                       lineA12.findIntersectionPoint(lineB34, p4);

      if (foundP4)
        faceVertexes[3] = p4;
      else
      {
        isValid = false;
        GEOSX_THROW( "P4 cannot be found!", InputError );
      }

      faceVertexes[2].setStatus(false);
      faceVertexes[6].setStatus(false);

    }
    // swap 6th and 4th element
    auto tmpVertex = faceVertexes[4];
    faceVertexes[4] = faceVertexes[6];
    faceVertexes[6] = tmpVertex;

    for (localIndex i = 0; i < vertexNum; ++ i)
    {
      if ( faceVertexes[i].getStatus() )
      {
        // valida vertex
        // intersection points found

        if (std::find(newPointVector.begin(), newPointVector.end(), faceVertexes[i]) ==  newPointVector.end())
          newPointVector.push_back(faceVertexes[i]);
      }
    }
    return newPointVector;
  }
};

} // end namespace geometryUtilities

} // end namespace cornerPointMesh

} // end namespace geosx

#endif //GEOSX_MESHUTILITIES_CORNERPOINTMESH_GEOMETRYUTILITIES_HPP_
