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
                           arrayView1d<R1Tensor const> const & points,
                           R1Tensor & center,
                           R1Tensor & normal,
                           real64 areaTolerance )
{
  R1Tensor v1,v2,vc;
  real64 area = 0.0;
  center = 0.0;
  normal=0.;

  if( numPoints > 2 )
  {
    const R1Tensor& x0 = points[pointsIndices[0]];
    for( localIndex a=0 ; a<(numPoints-2) ; ++a )
    {
      v1  = points[pointsIndices[a+1]];
      v2  = points[pointsIndices[a+2]];

      vc  = x0;
      vc += v1;
      vc += v2;

      v1 -= x0;
      v2 -= x0;

      R1Tensor triangleNormal;
      triangleNormal.Cross( v1,v2 );
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
      for( localIndex a=0 ; a<numPoints ; ++a )
      {
        GEOSX_LOG_RANK("Points: " << points[pointsIndices[a]](0) << " "
                      << points[pointsIndices[a]](1) << " "
                      << points[pointsIndices[a]](2) << " "
                      << pointsIndices[a]);
      }

      GEOSX_ERROR("Negative area found : " << area );
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
    area = Dot(x1_x0, x1_x0);
    area = sqrt(area);
  }
  return area;
}

real64 Centroid_3DPolygon( arrayView1d<localIndex const> const & pointsIndices,
                           arrayView1d<R1Tensor const> const & points,
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
                           arrayView1d<R1Tensor const> const & pointReferences,
                           arrayView1d<R1Tensor const> const & pointDisplacements,
                           R1Tensor & center,
                           R1Tensor & normal )
{
  R1Tensor v1,v2,vc;
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

    normal.Cross(v1,v2);
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
    R1Tensor x3_x1(v1);
    center += v1;

    v1  = pointReferences[pointsIndices[1]];
    v1 += pointDisplacements[pointsIndices[1]];
    x3_x1 -= v1;
    center += v1;

    v1  = pointReferences[pointsIndices[2]];
    v1 += pointDisplacements[pointsIndices[2]];
    R1Tensor x2_x0(v1);
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
    for( localIndex a=0 ; a<(numPoints-2) ; ++a )
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

      normal.Cross(v1,v2);
      const real64 triangleArea = normal.Normalize();
      area += triangleArea;
      vc *= triangleArea;
      center += vc;
    }
    if(area > 0.0)
    {
      center /= (area * 3.0);
      area *= 0.5;
    }
    else
    {
      GEOSX_ERROR("Zero area calculated!");
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
    area = Dot(x1_x0, x1_x0);
    area = sqrt(area);

    x1_x0[2] = 0.0;
    x1_x0.Normalize();

    normal[0] = -x1_x0[1];
    normal[1] = x1_x0[0];
    normal[2] = 0.0;
  }

  return area;
}

template<typename T>
int sgn( T val )
{
  return (T(0) < val) - (val < T(0));
}

bool IsPointInsidePolyhedron( arrayView1d<R1Tensor const> const & nodeCoordinates,
                              array1d<array1d<localIndex>> const & faceNodeIndicies,
                              R1Tensor const & point,
                              real64 const areaTolerance )
{
  localIndex const numFaces = faceNodeIndicies.size(0);
  R1Tensor faceCenter, faceNormal;
  int sign = 0;

  for (localIndex kf = 0; kf < numFaces; ++kf)
  {
    Centroid_3DPolygon( faceNodeIndicies[kf], nodeCoordinates, faceCenter, faceNormal, areaTolerance );
    faceCenter -= point;
    int const s = sgn( Dot( faceNormal, faceCenter ) );

    // all dot products should be non-negative (for outward normals) or non-positive (for inward normals)
    if (sign * s < 0)
    {
      return false;
    }
    sign = s;
  }

  return true;
}

real64 Centroid_3DPolygon( arrayView1d<localIndex const> const & pointsIndices,
                           arrayView1d<R1Tensor const> const & pointReferences,
                           arrayView1d<R1Tensor const> const & pointDisplacements,
                           R1Tensor & center,
                           R1Tensor & normal )
{ return Centroid_3DPolygon( pointsIndices.data(), pointsIndices.size(), pointReferences, pointDisplacements, center, normal ); }

real64 HexVolume( R1Tensor const * const X )
{
  R1Tensor X7_X1( X[7] );
  X7_X1 -= X[1];

  R1Tensor X6_X0( X[6] );
  X6_X0 -= X[0];

  R1Tensor X7_X2( X[7] );
  X7_X2 -= X[2];

  R1Tensor X3_X0( X[3] );
  X3_X0 -= X[0];

  R1Tensor X5_X0( X[5] );
  X5_X0 -= X[0];

  R1Tensor X7_X4( X[7] );
  X7_X4 -= X[4];

  R1Tensor X7_X1plusX6_X0( X7_X1 );
  X7_X1plusX6_X0 += X6_X0;

  R1Tensor X7_X2plusX5_X0( X7_X2 );
  X7_X2plusX5_X0 += X5_X0;

  R1Tensor X7_X4plusX3_X0( X7_X4 );
  X7_X4plusX3_X0 += X3_X0;

  return 1.0/12.0 * ( Dot( X7_X1plusX6_X0, Cross( X7_X2, X3_X0 ) ) +
                      Dot( X6_X0, Cross( X7_X2plusX5_X0, X7_X4 ) ) +
                      Dot( X7_X1, Cross( X5_X0, X7_X4plusX3_X0 ) ) );
}

real64 TetVolume( R1Tensor const * const X ) {
    R1Tensor X1_X0( X[1] );
    X1_X0 -= X[0];
    R1Tensor X2_X0( X[2] );
    X2_X0 -= X[0];
    R1Tensor X3_X0( X[3] );
    X3_X0 -= X[0];
    return std::fabs(Dot(X1_X0, Cross(X2_X0, X3_X0)) / 6.0);
}

real64 WedgeVolume( R1Tensor const * const X ) {
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
    return TetVolume(tet1) + TetVolume(tet2) + TetVolume(tet3);
}


real64 PyramidVolume( R1Tensor const * const X ) {
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
    return TetVolume(tet1) + TetVolume(tet2);
}
 
}
} /* namespace geosx */
