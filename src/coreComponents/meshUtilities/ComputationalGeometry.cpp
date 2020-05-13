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
 * @file ComputationalGeometry.cpp
 */

#include "ComputationalGeometry.hpp"

namespace geosx
{


namespace computationalGeometry
{

/**
 * Calculates the intersection between a line and a plane
 * @param[in] vector defining direction of the line
 * @param[in] 1 point of the line
 * @param[in] normal to plane
 * @param[in] plane origin
 * @return area of the convex 3D polygon
 */
R1Tensor LinePlaneIntersection( R1Tensor lineDir,
                                R1Tensor linePoint,
                                R1Tensor planeNormal,
                                R1Tensor planeOrigin )
{
  /* Find intersection line plane
   * line equation: p - (d*lineDir + linePoing) = 0;
   * plane equation: ( p - planeOrigin) * planeNormal = 0;
   * d = (planeOrigin - linePoint) * planeNormal / (lineDir * planeNormal )
   * pInt = d*lineDir+linePoint;
   */
  R1Tensor dummy;
  real64 d;
  // Intersection
  R1Tensor pInt;

  dummy = planeOrigin;
  dummy -= linePoint;
  d  = Dot( dummy, planeNormal );
  d /= Dot( lineDir, planeNormal );


  pInt    = linePoint;
  pInt[0] += d * lineDir[0];
  pInt[1] += d * lineDir[1];
  pInt[2] += d * lineDir[2];

  return pInt;
}
/**
 * Calculates the area of a polygon given the set of points defining it
 * @param[in] coordinates of the points
 * @param[in] number of points
 * @param[in] unit normal vector to the surface
 * @return area
 */
real64 ComputeSurfaceArea( array1d< R1Tensor > const & points,
                           localIndex const numPoints,
                           R1Tensor const & normal )
{
  // reorder points counterclockwise

  array1d< R1Tensor > pointsReordered = orderPointsCCW( points, numPoints, normal );

  real64 surfaceArea = 0.0;
  R1Tensor v1, v2;
  const R1Tensor & x0 = pointsReordered[0];

  for( localIndex a=0; a<(numPoints-2); ++a )
  {
    v1  = pointsReordered[a+1];
    v2  = pointsReordered[a+2];

    v1 -= x0;
    v2 -= x0;

    R1Tensor triangleNormal;
    triangleNormal.Cross( v1, v2 );
    const real64 triangleArea = triangleNormal.Normalize();

    surfaceArea += triangleArea;
  }
  surfaceArea *= 0.5;
  return surfaceArea;
}

/**
 * Given a set of points on a plane it orders them counterclockwise
 * @param[in] coordinates of the points
 * @param[in] number of points
 * @param[in] unit normal vector to the surface
 * @return reordered points
 */
array1d< R1Tensor > orderPointsCCW( array1d< R1Tensor > const & points,
                                    localIndex const numPoints,
                                    R1Tensor const & normal )
{
  array1d< R1Tensor > orderedPoints( numPoints );
  R1Tensor p0 = points[0];
  R1Tensor centroid = p0;

  std::vector< int > indices( numPoints );
  indices[0] = 0;
  real64 dot, det;
  std::vector< real64 > angle( numPoints );

  // compute centroid of the set of points

  for( localIndex a=1; a < numPoints; a++ )
  {
    centroid += points[a];
    indices[a] = a;
  }
  centroid /= numPoints;

  R1Tensor v0, v;
  v0  = centroid;
  v0 -= points[0];
  v0.Normalize();

  // compute angles
  angle[0] = 0;
  //std::cout << std::endl;
  for( localIndex a=1; a < numPoints; a++ )
  {
    v        = centroid;
    v       -= points[a];
    dot      = Dot( v, v0 );
    det      = Dot( normal, Cross( v, v0 ));
    angle[a] = std::atan2( det, dot );
    // std::cout << angle[a] << " - ";
  }

  // sort the indices
  std::sort( indices.begin(), indices.end(), [&]( int i, int j ){return angle[i]<angle[j];} );
  // std::cout << std::endl;

  // copy the points in the reorderedPoints array.
  for( localIndex a=0; a < numPoints; a++ )
  {
    // fill in with ordered
    // std::cout << indices[a] << " - ";
    orderedPoints[a] = points[indices[a]];
  }
  //std::cout << std::endl;

  return orderedPoints;
}


/**
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( localIndex const * const pointsIndices,
                           localIndex const numPoints,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           real64 areaTolerance )
{
  R1Tensor v1, v2, vc;
  real64 area = 0.0;
  center = 0.0;
  normal=0.;

  if( numPoints > 2 )
  {
    R1Tensor const x0 = points[pointsIndices[0]];
    for( localIndex a=0; a<(numPoints-2); ++a )
    {
      v1  = points[pointsIndices[a+1]];
      v2  = points[pointsIndices[a+2]];

      vc  = x0;
      vc += v1;
      vc += v2;

      v1 -= x0;
      v2 -= x0;

      R1Tensor triangleNormal;
      triangleNormal.Cross( v1, v2 );
      const real64 triangleArea = triangleNormal.Normalize();
      triangleNormal *= triangleArea;
      normal += triangleNormal;
      area += triangleArea;
      vc *= triangleArea;
      center += vc;
    }
    if( area > areaTolerance )
    {
      center /= (area * 3.0);
      normal.Normalize();
      area *= 0.5;
    }
    else if( area < -areaTolerance )
    {
      for( localIndex a=0; a<numPoints; ++a )
      {
        GEOSX_LOG_RANK( "Points: " << points[pointsIndices[a]]( 0 ) << " "
                                   << points[pointsIndices[a]]( 1 ) << " "
                                   << points[pointsIndices[a]]( 2 ) << " "
                                   << pointsIndices[a] );
      }

      GEOSX_ERROR( "Negative area found : " << area );
    }
    else
    {
      return 0.;
    }
  }
  else if( numPoints == 1 )
  {
    center = points[pointsIndices[0]];
  }
  else if( numPoints == 2 )
  {
    center  = points[pointsIndices[0]];

    //For 2D elements, a face is actually an edge with two nodes. We treat the
    // length of this edge as the surface area and use it in the calculation of
    // tractions.
    R1Tensor x1_x0;
    x1_x0 = points[pointsIndices[1]];
    center += x1_x0;
    center *= 0.5;

    x1_x0 -= points[pointsIndices[0]];
    area = Dot( x1_x0, x1_x0 );
    area = sqrt( area );
  }
  return area;
}

real64 Centroid_3DPolygon( arrayView1d< localIndex const > const & pointsIndices,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           real64 areaTolerance )
{ return Centroid_3DPolygon( pointsIndices.data(), pointsIndices.size(), points, center, normal, areaTolerance ); }

/**
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] pointReferences 3D reference point list
 * @param[in] pointDisplacements 3D displacement list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon( localIndex const * const pointsIndices,
                           localIndex const numPoints,
                           arrayView1d< R1Tensor const > const & pointReferences,
                           arrayView1d< R1Tensor const > const & pointDisplacements,
                           R1Tensor & center,
                           R1Tensor & normal )
{
  R1Tensor v1, v2, vc;
  real64 area = 0.0;
  center = 0.0;

  if( numPoints==3 )
  {
    const localIndex a0 = pointsIndices[0];
    const localIndex a1 = pointsIndices[1];
    const localIndex a2 = pointsIndices[2];

    v1  = pointReferences[a1];
    v1 += pointDisplacements[a1];
    v2  = pointReferences[a2];
    v2 += pointDisplacements[a2];

    vc  = pointReferences[a0];
    vc += pointDisplacements[a0];
    vc += v1;
    vc += v2;

    v1 -= pointReferences[a0];
    v1 -= pointDisplacements[a0];
    v2 -= pointReferences[a0];
    v2 -= pointDisplacements[a0];

    normal.Cross( v1, v2 );
    const real64 triangleArea = normal.Normalize();
    area += triangleArea;
    area *= 0.5;
    vc /= 3.0;
    center += vc;
  }
  else if( numPoints==4 )
  {
    v1  = pointReferences[pointsIndices[3]];
    v1 += pointDisplacements[pointsIndices[3]];
    R1Tensor x3_x1( v1 );
    center += v1;

    v1  = pointReferences[pointsIndices[1]];
    v1 += pointDisplacements[pointsIndices[1]];
    x3_x1 -= v1;
    center += v1;

    v1  = pointReferences[pointsIndices[2]];
    v1 += pointDisplacements[pointsIndices[2]];
    R1Tensor x2_x0( v1 );
    center += v1;

    v1  = pointReferences[pointsIndices[0]];
    v1 += pointDisplacements[pointsIndices[0]];
    x2_x0 -= v1;
    center += v1;

    normal.Cross( x2_x0, x3_x1 );

    area = 0.5 * normal.Normalize();
    center *= 0.25;
  }
  else if( numPoints>4 )
  {
    const localIndex a0 = pointsIndices[0];
    for( localIndex a=0; a<(numPoints-2); ++a )
    {
      const localIndex a1 = pointsIndices[a+1];
      const localIndex a2 = pointsIndices[a+2];

      v1  = pointReferences[a1];
      v1 += pointDisplacements[a1];
      v2  = pointReferences[a2];
      v2 += pointDisplacements[a2];

      vc  = pointReferences[a0];
      vc += pointDisplacements[a0];
      vc += v1;
      vc += v2;

      v1 -= pointReferences[a0];
      v1 -= pointDisplacements[a0];
      v2 -= pointReferences[a0];
      v2 -= pointDisplacements[a0];

      normal.Cross( v1, v2 );
      const real64 triangleArea = normal.Normalize();
      area += triangleArea;
      vc *= triangleArea;
      center += vc;
    }
    if( area > 0.0 )
    {
      center /= (area * 3.0);
      area *= 0.5;
    }
    else
    {
      GEOSX_ERROR( "Zero area calculated!" );
    }
  }
  else if( numPoints==1 )
  {
    center = pointReferences[0];
    center += pointDisplacements[0];
  }
  else if( numPoints==2 )
  {
    center  = pointReferences[pointsIndices[0]];
    center += pointDisplacements[pointsIndices[0]];

    //For 2D elements, a face is actually an edge with two nodes. We treat the
    // length of this edge as the surface area and use it in the calculation of
    // tractions.
    R1Tensor x1_x0;
    x1_x0 = pointReferences[pointsIndices[1]];
    x1_x0 += pointDisplacements[pointsIndices[1]];
    center += x1_x0;
    center *= 0.5;

    x1_x0 -= pointReferences[pointsIndices[0]];
    x1_x0 -= pointDisplacements[pointsIndices[0]];
    area = Dot( x1_x0, x1_x0 );
    area = sqrt( area );

    x1_x0[2] = 0.0;
    x1_x0.Normalize();

    normal[0] = -x1_x0[1];
    normal[1] = x1_x0[0];
    normal[2] = 0.0;
  }

  return area;
}

template< typename T >
int sgn( T val )
{
  return (T( 0 ) < val) - (val < T( 0 ));
}

bool IsPointInsidePolyhedron( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodeCoordinates,
                              array1d< array1d< localIndex > > const & faceNodeIndicies,
                              R1Tensor const & point,
                              real64 const areaTolerance )
{
  localIndex const numFaces = faceNodeIndicies.size( 0 );
  R1Tensor faceCenter, faceNormal;
  int sign = 0;

  for( localIndex kf = 0; kf < numFaces; ++kf )
  {
    Centroid_3DPolygon( faceNodeIndicies[kf], nodeCoordinates, faceCenter, faceNormal, areaTolerance );
    faceCenter -= point;
    int const s = sgn( Dot( faceNormal, faceCenter ) );

    // all dot products should be non-negative (for outward normals) or non-positive (for inward normals)
    if( sign * s < 0 )
    {
      return false;
    }
    sign = s;
  }

  return true;
}

real64 Centroid_3DPolygon( arrayView1d< localIndex const > const & pointsIndices,
                           arrayView1d< R1Tensor const > const & pointReferences,
                           arrayView1d< R1Tensor const > const & pointDisplacements,
                           R1Tensor & center,
                           R1Tensor & normal )
{ return Centroid_3DPolygon( pointsIndices.data(), pointsIndices.size(), pointReferences, pointDisplacements, center, normal ); }

//GEOSX_HOST_DEVICE
//real64 HexVolume( R1Tensor const * const X )
//{
//  R1Tensor X7_X1( X[7] );
//  X7_X1 -= X[1];
//
//  R1Tensor X6_X0( X[6] );
//  X6_X0 -= X[0];
//
//  R1Tensor X7_X2( X[7] );
//  X7_X2 -= X[2];
//
//  R1Tensor X3_X0( X[3] );
//  X3_X0 -= X[0];
//
//  R1Tensor X5_X0( X[5] );
//  X5_X0 -= X[0];
//
//  R1Tensor X7_X4( X[7] );
//  X7_X4 -= X[4];
//
//  R1Tensor X7_X1plusX6_X0( X7_X1 );
//  X7_X1plusX6_X0 += X6_X0;
//
//  R1Tensor X7_X2plusX5_X0( X7_X2 );
//  X7_X2plusX5_X0 += X5_X0;
//
//  R1Tensor X7_X4plusX3_X0( X7_X4 );
//  X7_X4plusX3_X0 += X3_X0;
//
//  return 1.0/12.0 * ( Dot( X7_X1plusX6_X0, Cross( X7_X2, X3_X0 ) ) +
//                      Dot( X6_X0, Cross( X7_X2plusX5_X0, X7_X4 ) ) +
//                      Dot( X7_X1, Cross( X5_X0, X7_X4plusX3_X0 ) ) );
//}

real64 TetVolume( R1Tensor const * const X )
{
  R1Tensor X1_X0( X[1] );
  X1_X0 -= X[0];
  R1Tensor X2_X0( X[2] );
  X2_X0 -= X[0];
  R1Tensor X3_X0( X[3] );
  X3_X0 -= X[0];
  return std::fabs( Dot( X1_X0, Cross( X2_X0, X3_X0 )) / 6.0 );
}

real64 WedgeVolume( R1Tensor const * const X )
{
  R1Tensor tet1[4];
  tet1[0] = X[0];
  tet1[1] = X[1];
  tet1[2] = X[2];
  tet1[3] = X[4];
  R1Tensor tet2[4];
  tet2[0] = X[0];
  tet2[1] = X[2];
  tet2[2] = X[4];
  tet2[3] = X[5];
  R1Tensor tet3[4];
  tet3[0] = X[0];
  tet3[1] = X[3];
  tet3[2] = X[4];
  tet3[3] = X[5];
  return TetVolume( tet1 ) + TetVolume( tet2 ) + TetVolume( tet3 );
}


real64 PyramidVolume( R1Tensor const * const X )
{
  R1Tensor tet1[4];
  tet1[0] = X[0];
  tet1[1] = X[1];
  tet1[2] = X[2];
  tet1[3] = X[4];
  R1Tensor tet2[4];
  tet2[0] = X[0];
  tet2[1] = X[2];
  tet2[2] = X[3];
  tet2[3] = X[4];
  return TetVolume( tet1 ) + TetVolume( tet2 );
}



template< typename NODEMAP >
R1Tensor GetBoundingBox( localIndex elemIndex,
                         NODEMAP const & pointIndices,
                         arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates )
{
  localIndex constexpr dim = 3;

  // these arrays will store the min and max coordinates of the elem in each direction
  R1Tensor minCoords( 1e99 );
  R1Tensor maxCoords( -1e99 );

  // loop over all the vertices of the element to get the min and max coords
  for( localIndex a = 0; a < pointIndices.size( 1 ); ++a )
  {
    localIndex const id = pointIndices( elemIndex, a );
    R1Tensor const coords = pointCoordinates[id];

    for( localIndex d = 0; d < dim; ++d )
    {
      if( coords[d] < minCoords[d] )
      {
        minCoords[d] = coords[d];
      }
      if( coords[d] > maxCoords[d] )
      {
        maxCoords[d] = coords[d];
      }
    }
  }

  // compute the dimensions of the bounding box
  R1Tensor box( 0 );
  for( localIndex d = 0; d < dim; ++d )
  {
    box[d] = maxCoords[d] - minCoords[d];
  }

  return box;
}

template R1Tensor GetBoundingBox( localIndex elemIndex,
                                  InterObjectRelation< array2d< localIndex, RAJA::PERM_IJ > > const & pointIndices,
                                  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates );
template R1Tensor GetBoundingBox( localIndex elemIndex,
                                  InterObjectRelation< array2d< localIndex, RAJA::PERM_JI > > const & pointIndices,
                                  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates );


}



} /* namespace geosx */
