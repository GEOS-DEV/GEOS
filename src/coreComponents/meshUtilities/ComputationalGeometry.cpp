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
#include "LvArray/src/tensorOps.hpp"

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

static real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();

void FixNormalOrientation_3D( arraySlice1d< real64 > const normal )
{
  real64 const orientationTolerance = 1.e+1*machinePrecision;

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
    else if( std::fabs( normal[ 1 ] ) < orientationTolerance )
    {
      // If needed, third check: align with x direction
      if( normal[ 0 ] <= -orientationTolerance )
      {
        LvArray::tensorOps::scale< 3 >( normal, -1.0 );
      }
    }
  }
}

void RotationMatrix_3D( arraySlice1d< real64 const > const normal,
                        arraySlice2d< real64 > const rotationMatrix )
{
  real64 m1[ 3 ] = { normal[ 2 ], 0.0, -normal[ 0 ] };
  real64 m2[ 3 ] = { 0.0, normal[ 2 ], -normal[ 1 ] };
  real64 const norm_m1 = LvArray::tensorOps::normalize< 3 >( m1 );
  real64 const norm_m2 = LvArray::tensorOps::normalize< 3 >( m2 );

  // If present, looks for a vector with 0 norm
  // Fix the uncertain case of norm_m1 very close to norm_m2
  if( norm_m1+1.e+2*machinePrecision > norm_m2 )
  {
    LvArray::tensorOps::crossProduct( m2, normal, m1 );
    LvArray::tensorOps::normalize< 3 >( m2 );
  }
  else
  {
    LvArray::tensorOps::crossProduct( m1, normal, m2 );
    LvArray::tensorOps::scale< 3 >( m1, -1 );
    LvArray::tensorOps::normalize< 3 >( m1 );
  }

  // Save everything in the standard form (3x3 rotation matrix)
  rotationMatrix( 0, 0 ) = normal[ 0 ];
  rotationMatrix( 1, 0 ) = normal[ 1 ];
  rotationMatrix( 2, 0 ) = normal[ 2 ];
  rotationMatrix( 0, 1 ) = m1[ 0 ];
  rotationMatrix( 1, 1 ) = m1[ 1 ];
  rotationMatrix( 2, 1 ) = m1[ 2 ];
  rotationMatrix( 0, 2 ) = m2[ 0 ];
  rotationMatrix( 1, 2 ) = m2[ 1 ];
  rotationMatrix( 2, 2 ) = m2[ 2 ];

  GEOSX_ERROR_IF( std::fabs( LvArray::tensorOps::determinant< 3 >( rotationMatrix ) - 1.0 ) > 1.e+1*machinePrecision,
                  "Rotation matrix with determinant different from +1.0" );

  return;
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
  real64 faceCenter[ 3 ], faceNormal[ 3 ];
  int sign = 0;

  for( localIndex kf = 0; kf < numFaces; ++kf )
  {
    Centroid_3DPolygon( faceNodeIndicies[kf], nodeCoordinates, faceCenter, faceNormal, areaTolerance );

    faceCenter[ 0 ] -= point[ 0 ];
    faceCenter[ 1 ] -= point[ 1 ];
    faceCenter[ 2 ] -= point[ 2 ];
    int const s = sgn( faceNormal[ 0 ] * faceCenter[ 0 ] + faceNormal[ 1 ] * faceCenter[ 1 ] + faceNormal[ 2 ] * faceCenter[ 2 ] );

    // TODO: Replace the above with
    // LvArray::tensorOps::subtract< 3 >( faceCenter, point );
    // int const s = sgn( LvArray::tensorOps::innerProduct< 3 >( faceNormal, faceCenter ) );

    // all dot products should be non-negative (for outward normals) or non-positive (for inward normals)
    if( sign * s < 0 )
    {
      return false;
    }
    sign = s;
  }

  return true;
}

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
