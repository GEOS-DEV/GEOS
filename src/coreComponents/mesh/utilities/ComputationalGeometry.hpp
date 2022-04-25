/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ComputationalGeometry.hpp
 */

#ifndef GEOSX_MESH_UTILITIES_COMPUTATIONALGEOMETRY_HPP_
#define GEOSX_MESH_UTILITIES_COMPUTATIONALGEOMETRY_HPP_

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "finiteElement/elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "finiteElement/elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp"
#include "finiteElement/elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp"
#include "finiteElement/elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp"
#include "LvArray/src/output.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
namespace computationalGeometry
{

/// Machine epsilon for double-precision calculations
constexpr real64 machinePrecision = std::numeric_limits< real64 >::epsilon();

/**
 * @brief Calculate the intersection between a line and a plane.
 * @tparam LINEDIR_TYPE the type of @p lineDir
 * @tparam POINT_TYPE the type of @p linePoint
 * @tparam NORMAL_TYPE the type of @p planeNormal
 * @tparam ORIGIN_TYPE the type of @p planeOrigin
 * @tparam INTPOINT_TYPE the type of @p instersectionPoint
 * @param[in] lineDir vector defining direction of the line
 * @param[in] linePoint one point of the line
 * @param[in] planeNormal normal to plane
 * @param[in] planeOrigin plane origin
 * @param[out] intersectionPoint the intersection point
 */
template< typename LINEDIR_TYPE,
          typename POINT_TYPE,
          typename NORMAL_TYPE,
          typename ORIGIN_TYPE,
          typename INTPOINT_TYPE >
void LinePlaneIntersection( LINEDIR_TYPE const & lineDir,
                            POINT_TYPE const & linePoint,
                            NORMAL_TYPE const & planeNormal,
                            ORIGIN_TYPE const & planeOrigin,
                            INTPOINT_TYPE & intersectionPoint )
{
  /* Find intersection line plane
   * line equation: p - (d*lineDir + linePoing) = 0;
   * plane equation: ( p - planeOrigin) * planeNormal = 0;
   * d = (planeOrigin - linePoint) * planeNormal / (lineDir * planeNormal )
   * pInt = d*lineDir+linePoint;
   */
  real64 dummy[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( planeOrigin );
  LvArray::tensorOps::subtract< 3 >( dummy, linePoint );
  real64 const d = LvArray::tensorOps::AiBi< 3 >( dummy, planeNormal ) /
                   LvArray::tensorOps::AiBi< 3 >( lineDir, planeNormal );

  LvArray::tensorOps::copy< 3 >( intersectionPoint, linePoint );
  LvArray::tensorOps::scaledAdd< 3 >( intersectionPoint, lineDir, d );
}

/**
 * @brief Reorder a set of points counter-clockwise.
 * @tparam NORMAL_TYPE the type of @p normal
 * @param[in] points coordinates of the points
 * @param[in] normal vector normal to the plane
 * @return an std::vector containing the original indices of the reordered points.
 */
template< typename NORMAL_TYPE >
array1d< int >  orderPointsCCW( arrayView2d< real64 > const & points,
                                NORMAL_TYPE const & normal )
{
  localIndex const numPoints = points.size( 0 );

  array2d< real64 > orderedPoints( numPoints, 3 );

  array1d< int > indices( numPoints );
  array1d< real64 > angle( numPoints );

  // compute centroid of the set of points
  real64 centroid[3];
  LvArray::tensorOps::fill< 3 >( centroid, 0 );
  for( localIndex a = 0; a < numPoints; ++a )
  {
    LvArray::tensorOps::add< 3 >( centroid, points[ a ] );
    indices[ a ] = a;
  }

  LvArray::tensorOps::scale< 3 >( centroid, 1.0 / numPoints );

  real64 v0[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( centroid );
  LvArray::tensorOps::subtract< 3 >( v0, points[ 0 ] );
  LvArray::tensorOps::normalize< 3 >( v0 );

  // compute angles
  angle[ 0 ] = 0;
  for( localIndex a = 1; a < numPoints; ++a )
  {
    real64 v[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( centroid );
    LvArray::tensorOps::subtract< 3 >( v, points[ a ] );
    real64 const dot = LvArray::tensorOps::AiBi< 3 >( v, v0 );

    real64 crossProduct[ 3 ];
    LvArray::tensorOps::crossProduct( crossProduct, v, v0 );
    real64 const det = LvArray::tensorOps::AiBi< 3 >( normal, crossProduct );

    angle[ a ] = std::atan2( det, dot );
  }

  // sort the indices
  std::sort( indices.begin(), indices.end(), [&]( int i, int j ) { return angle[ i ] < angle[ j ]; } );

  // copy the points in the reorderedPoints array.
  for( localIndex a=0; a < numPoints; a++ )
  {
    // fill in with ordered
    LvArray::tensorOps::copy< 3 >( orderedPoints[ a ], points[ indices[ a ] ] );
  }

  for( localIndex a = 0; a < numPoints; a++ )
  {
    LvArray::tensorOps::copy< 3 >( points[a], orderedPoints[a] );
  }

  return indices;
}

/**
 * @brief Calculate the area of a polygon given the set of points in ccw order defining it.
 * @tparam NORMAL_TYPE the type of @p normal
 * @param[in] points coordinates of the points
 * @param[in] normal vector normal to the plane
 * @return the area of the polygon
 */
template< typename NORMAL_TYPE >
real64 ComputeSurfaceArea( arrayView2d< real64 const > const & points,
                           NORMAL_TYPE const && normal )
{
  real64 surfaceArea = 0.0;

  array2d< real64 > orderedPoints( points.size( 0 ), 3 );

  for( localIndex a = 0; a < points.size( 0 ); a++ )
  {
    LvArray::tensorOps::copy< 3 >( orderedPoints[a], points[a] );
  }

  orderPointsCCW( orderedPoints, normal );

  for( localIndex a = 0; a < points.size( 0 ) - 2; ++a )
  {
    real64 v1[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( orderedPoints[ a + 1 ] );
    real64 v2[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( orderedPoints[ a + 2 ] );

    LvArray::tensorOps::subtract< 3 >( v1, orderedPoints[ 0 ] );
    LvArray::tensorOps::subtract< 3 >( v2, orderedPoints[ 0 ] );

    real64 triangleNormal[ 3 ];
    LvArray::tensorOps::crossProduct( triangleNormal, v1, v2 );
    surfaceArea += LvArray::tensorOps::l2Norm< 3 >( triangleNormal );
  }

  return surfaceArea * 0.5;
}

/**
 * @brief Calculate the centroid of a convex 3D polygon as well as the normal and the rotation matrix.
 * @tparam CENTER_TYPE The type of @p center.
 * @tparam NORMAL_TYPE The type of @p normal.
 * @param[in] pointsIndices list of index references for the points array in
 *   order (CW or CCW) about the polygon loop
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal normal to the face
 * @param[in] areaTolerance tolerance used in the geometric computations
 * @return area of the convex 3D polygon
 * @details if area < - areaTolerance, this function will throw an error,
 *          and if (- areaTolerance <= area <= areaTolerance), the area is set to zero
 */
template< typename CENTER_TYPE, typename NORMAL_TYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 centroid_3DPolygon( arraySlice1d< localIndex const > const pointsIndices,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & points,
                           CENTER_TYPE && center,
                           NORMAL_TYPE && normal,
                           real64 const areaTolerance = 0.0 )
{
  real64 area = 0.0;
  LvArray::tensorOps::fill< 3 >( center, 0 );
  LvArray::tensorOps::fill< 3 >( normal, 0 );

  GEOSX_ERROR_IF_LT( pointsIndices.size(), 2 );
  for( localIndex a=0; a<(pointsIndices.size()-2); ++a )
  {
    real64 v1[ 3 ], v2[ 3 ], vc[ 3 ];

    LvArray::tensorOps::copy< 3 >( v1, points[ pointsIndices[ a + 1 ] ] );
    LvArray::tensorOps::copy< 3 >( v2, points[ pointsIndices[ a + 2 ] ] );

    LvArray::tensorOps::copy< 3 >( vc, points[ pointsIndices[ 0 ] ] );
    LvArray::tensorOps::add< 3 >( vc, v1 );
    LvArray::tensorOps::add< 3 >( vc, v2 );

    LvArray::tensorOps::subtract< 3 >( v1, points[ pointsIndices[ 0 ] ] );
    LvArray::tensorOps::subtract< 3 >( v2, points[ pointsIndices[ 0 ] ] );

    real64 triangleNormal[ 3 ];
    LvArray::tensorOps::crossProduct( triangleNormal, v1, v2 );
    real64 const triangleArea = LvArray::tensorOps::l2Norm< 3 >( triangleNormal );

    LvArray::tensorOps::add< 3 >( normal, triangleNormal );

    area += triangleArea;
    LvArray::tensorOps::scaledAdd< 3 >( center, vc, triangleArea );
  }
  if( area > areaTolerance )
  {
    LvArray::tensorOps::scale< 3 >( center, 1.0 / ( area * 3.0 ) );
    LvArray::tensorOps::normalize< 3 >( normal );
    area *= 0.5;
  }
  else if( area < -areaTolerance )
  {
    for( localIndex a=0; a<pointsIndices.size(); ++a )
    {
      GEOSX_LOG_RANK( "Points: " << points[ pointsIndices[ a ] ] << " " << pointsIndices[ a ] );
    }
    GEOSX_ERROR( "Negative area found : " << area );
  }
  else
  {
    return 0.0;
  }

  return area;
}

/**
 * @brief Change the orientation of the input vector to be consistent in a global sense.
 * @tparam NORMAL_TYPE type of @p normal
 * @param[inout] normal normal to the face
 */
template< typename NORMAL_TYPE >
GEOSX_HOST_DEVICE
void FixNormalOrientation_3D( NORMAL_TYPE && normal )
{
  real64 const orientationTolerance = 10 * machinePrecision;

  // Orient local normal in global sense.
  // First check: align with z direction
  if( normal[ 2 ] <= -orientationTolerance )
  {
    LvArray::tensorOps::scale< 3 >( normal, -1.0 );
  }
  else if( std::fabs( normal[ 2 ] ) < orientationTolerance )
  {
    // If needed, second check: align with y direction
    if( normal[ 1 ] <= -orientationTolerance )
    {
      LvArray::tensorOps::scale< 3 >( normal, -1.0 );
    }
    else if( fabs( normal[ 1 ] ) < orientationTolerance )
    {
      // If needed, third check: align with x direction
      if( normal[ 0 ] <= -orientationTolerance )
      {
        LvArray::tensorOps::scale< 3 >( normal, -1.0 );
      }
    }
  }
}

/**
 * @brief Calculate the rotation matrix for a face in the 3D space
 * @tparam NORMAL_TYPE type of @p normal
 * @tparam MATRIX_TYPE type of @p rotationMatrix
 * @param[in] normal normal to the face
 * @param[out] rotationMatrix rotation matrix for the face
 */
template< typename NORMAL_TYPE, typename MATRIX_TYPE >
GEOSX_HOST_DEVICE
void RotationMatrix_3D( NORMAL_TYPE const & normal,
                        MATRIX_TYPE && rotationMatrix )
{
  real64 m1[ 3 ] = { normal[ 2 ], 0.0, -normal[ 0 ] };
  real64 m2[ 3 ] = { 0.0, normal[ 2 ], -normal[ 1 ] };
  real64 const norm_m1 = LvArray::tensorOps::l2Norm< 3 >( m1 );
  real64 const norm_m2 = LvArray::tensorOps::l2Norm< 3 >( m2 );

  // If present, looks for a vector with 0 norm
  // Fix the uncertain case of norm_m1 very close to norm_m2
  if( norm_m1+1.e+2*machinePrecision > norm_m2 )
  {
    LvArray::tensorOps::crossProduct( m2, normal, m1 );
    LvArray::tensorOps::normalize< 3 >( m2 );
    LvArray::tensorOps::normalize< 3 >( m1 );
  }
  else
  {
    LvArray::tensorOps::crossProduct( m1, normal, m2 );
    LvArray::tensorOps::scale< 3 >( m1, -1 );
    LvArray::tensorOps::normalize< 3 >( m1 );
    LvArray::tensorOps::normalize< 3 >( m2 );
  }

  // Save everything in the standard form (3x3 rotation matrix)
  rotationMatrix[ 0 ][ 0 ] = normal[ 0 ];
  rotationMatrix[ 1 ][ 0 ] = normal[ 1 ];
  rotationMatrix[ 2 ][ 0 ] = normal[ 2 ];
  rotationMatrix[ 0 ][ 1 ] = m1[ 0 ];
  rotationMatrix[ 1 ][ 1 ] = m1[ 1 ];
  rotationMatrix[ 2 ][ 1 ] = m1[ 2 ];
  rotationMatrix[ 0 ][ 2 ] = m2[ 0 ];
  rotationMatrix[ 1 ][ 2 ] = m2[ 1 ];
  rotationMatrix[ 2 ][ 2 ] = m2[ 2 ];

  GEOSX_ERROR_IF( fabs( LvArray::tensorOps::determinant< 3 >( rotationMatrix ) - 1.0 ) > 1.e+1 * machinePrecision,
                  "Rotation matrix with determinant different from +1.0" );
}

/**
 * @brief Return the sign of a given value as an integer.
 * @tparam T type of value
 * @param val the value in question
 * @return -1, 0 or 1 depending on whether the value is negative, zero or positive
 */
template< typename T >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
int sign( T const val )
{
  return (T( 0 ) < val) - (val < T( 0 ));
}


/**
 * @brief Check if a point is inside a convex polyhedron (3D polygon)
 * @tparam POINT_TYPE type of @p point
 * @param[in] nodeCoordinates a global array of nodal coordinates
 * @param[in] faceIndices global indices of the faces of the cell
 * @param[in] facesToNodes map from face to nodes
 * @param[in] elemCenter coordinates of the element center
 * @param[in] point coordinates of the query point
 * @param[in] areaTolerance same as in centroid_3DPolygon
 * @return whether the point is inside
 *
 * @note For faces with n>3 nodes that are non-planar, average normal is used
 */
template< typename POINT_TYPE >
GEOSX_HOST_DEVICE
bool isPointInsidePolyhedron( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodeCoordinates,
                              arraySlice1d< localIndex const > const & faceIndices,
                              ArrayOfArraysView< localIndex const > const & facesToNodes,
                              POINT_TYPE const & elemCenter,
                              POINT_TYPE const & point,
                              real64 const areaTolerance = 0.0 )
{
  localIndex const numFaces = faceIndices.size();
  R1Tensor faceCenter, faceNormal, cellToFaceVec;
  int prevSign = 0;

  for( localIndex kf = 0; kf < numFaces; ++kf )
  {
    // compute the face normal at this face
    localIndex const faceIndex = faceIndices[kf];
    centroid_3DPolygon( facesToNodes[faceIndex], nodeCoordinates, faceCenter, faceNormal, areaTolerance );

    // make sure that the normal is outward pointing
    LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
    LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );
    if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    // compute the vector face center to query point
    LvArray::tensorOps::subtract< 3 >( faceCenter, point );
    int const s = sign( LvArray::tensorOps::AiBi< 3 >( faceNormal, faceCenter ) );

    // all dot products should be non-negative (we enforce outward normals)
    if( prevSign * s < 0 )
    {
      return false;
    }
    prevSign = s;
  }
  return true;
}


/**
 * @brief Compute the dimensions of the bounding box containing the element
 *   defined here by the coordinates of its vertices.
 * @tparam VEC_TYPE type of @p boxDims
 * @param[in] elemIndex index of the element in pointIndices.
 * @param[in] pointIndices the indices of the vertices in pointCoordinates.
 * @param[in] pointCoordinates the vertices coordinates.
 * @param[out] boxDims The dimensions of the bounding box.
 */
template< typename VEC_TYPE >
GEOSX_HOST_DEVICE
void getBoundingBox( localIndex const elemIndex,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const & pointIndices,
                     arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates,
                     VEC_TYPE && boxDims )
{
  // This holds the min coordinates of the set in each direction
  R1Tensor minCoords = { LvArray::NumericLimits< real64 >::max,
                         LvArray::NumericLimits< real64 >::max,
                         LvArray::NumericLimits< real64 >::max };

  // boxDims is used to hold the max coordinates.
  LvArray::tensorOps::fill< 3 >( boxDims, LvArray::NumericLimits< real64 >::min );

  // loop over all the vertices of the element to get the min and max coords
  for( localIndex a = 0; a < pointIndices.size( 1 ); ++a )
  {
    localIndex const id = pointIndices( elemIndex, a );
    for( localIndex d = 0; d < 3; ++d )
    {
      minCoords[ d ] = fmin( minCoords[ d ], pointCoordinates( id, d ) );
      boxDims[ d ] = fmax( boxDims[ d ], pointCoordinates( id, d ) );
    }
  }

  LvArray::tensorOps::subtract< 3 >( boxDims, minCoords );
}

/**
 * @brief Compute the volume of an element (tetrahedron, pyramid, wedge, hexahedron)
 * @tparam FE_TYPE the type of finite element space
 * @param[in] X vertices of the element
 * @return the volume of the element
 */
template< typename FE_TYPE >
GEOSX_HOST_DEVICE inline
real64 elementVolume( real64 const (&X)[FE_TYPE::numNodes][3] )
{
  real64 result{};
  for( localIndex q=0; q<FE_TYPE::numQuadraturePoints; ++q )
  {
    result = result + FE_TYPE::transformedQuadratureWeight( q, X );
  }
  return result;
}

///**
// * @brief Compute the volume of an hexahedron
// * @param[in] X vertices of the hexahedron
// * @return the volume of the hexahedron
// */
//GEOSX_HOST_DEVICE
//inline
//real64 hexVolume( real64 const (&X)[8][3] )
//{
//  return elementVolume< finiteElement::H1_Hexahedron_Lagrange1_GaussLegendre2 >( X );
//}
//
///**
// * @brief Compute the volume of an tetrahedron
// * @param[in] X vertices of the tetrahedron
// * @return the volume of the tetrahedron
// */
//GEOSX_HOST_DEVICE
//inline
//real64 tetVolume( real64 const (&X)[4][3] )
//{
//  return finiteElement::H1_Tetrahedron_Lagrange1_Gauss1::transformedQuadratureWeight( 0, X );
//}
//
///**
// * @brief Compute the volume of a wedge
// * @param[in] X vertices of the wedge
// * @return the volume of the wedge
// */
//GEOSX_HOST_DEVICE
//inline
//real64 wedgeVolume( real64 const (&X)[6][3] )
//{
//  constexpr integer numQuadraturePoints = 6;
//  real64 result =  finiteElement::H1_Wedge_Lagrange1_Gauss6::transformedQuadratureWeight( 0, X );
//  for( localIndex q=1; q<numQuadraturePoints; ++q )
//  {
//    result = result + finiteElement::H1_Wedge_Lagrange1_Gauss6::transformedQuadratureWeight( q, X );
//  }
//  return result;
//}
//
///**
// * @brief Compute the volume of a pyramid
// * @param[in] X vertices of the pyramid
// * @return the volume of the pyramid
// */
//GEOSX_HOST_DEVICE
//inline
//real64 pyramidVolume( real64 const (&X)[5][3] )
//{
//  constexpr integer numQuadraturePoints = 5;
//  real64 result =  finiteElement::H1_Pyramid_Lagrange1_Gauss5::transformedQuadratureWeight( 0, X );
//  for( localIndex q=1; q<numQuadraturePoints; ++q )
//  {
//    result = result + finiteElement::H1_Pyramid_Lagrange1_Gauss5::transformedQuadratureWeight( q, X );
//  }
//  return result;
//}

} /* namespace computationalGeometry */
} /* namespace geosx */

#endif /* GEOSX_MESH_UTILITIES_COMPUTATIONALGEOMETRY_HPP_ */
