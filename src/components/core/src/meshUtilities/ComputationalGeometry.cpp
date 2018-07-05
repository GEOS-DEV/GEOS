/*
 * ComputationalGeometry.cpp
 *
 *  Created on: Jun 26, 2018
 *      Author: settgast
 */

#include "ComputationalGeometry.hpp"

namespace geosx
{
namespace computationalGeometry
{

/**
 * @author settgast
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon(const localIndex_array& pointsIndices,
                         const array<R1Tensor>& points,
                         R1Tensor& center,
                         R1Tensor& normal )
{
  R1Tensor v1,v2,vc;
  const localIndex_array::size_type n = pointsIndices.size();
  real64 area = 0.0;
  center = 0.0;

  if( n>2 )
  {
    const R1Tensor& x0 = points[pointsIndices[0]];
    for( localIndex_array::size_type a=0 ; a<(n-2) ; ++a )
    {
      v1  = points[pointsIndices[a+1]];
      v2  = points[pointsIndices[a+2]];

      vc  = x0;
      vc += v1;
      vc += v2;

      v1 -= x0;
      v2 -= x0;

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
      center = 0.0;
      for( localIndex_array::size_type a=0 ; a<n ; ++a )
      {
        center += points[a];
      }
      std::cout << "Randy's bug: area = " << area << std::endl;
      for( localIndex_array::size_type a=0 ; a<n ; ++a )
        std::cout << points[pointsIndices[a]](0) << " "
                  << points[pointsIndices[a]](1) << " "
                  << points[pointsIndices[a]](2) << " "
                  << pointsIndices[a] << std::endl;
      GEOS_ERROR("");
      center /= n;
    }
  }
  else if( n==1 )
  {
    center = points[pointsIndices[0]];
  }
  else if( n==2 )
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

/**
 * @author settgast
 * Calculates the centroid of a convex 3D polygon as well as the normal
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] pointReferences 3D reference point list
 * @param[in] pointDisplacements 3D displacement list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
real64 Centroid_3DPolygon(const localIndex_array& pointsIndices,
                         const array<R1Tensor>& pointReferences,
                         const array<R1Tensor>& pointDisplacements,
                         R1Tensor& center,
                         R1Tensor& normal )
{
  R1Tensor v1,v2,vc;
  const localIndex_array::size_type n = pointsIndices.size();
  real64 area = 0.0;
  center = 0.0;

  if( n==3 )
  {
    const localIndex_array::size_type a0 = pointsIndices[0];
    const localIndex_array::size_type a1 = pointsIndices[1];
    const localIndex_array::size_type a2 = pointsIndices[2];

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
  else if( n==4 )
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
  else if( n>4 )
  {
    const localIndex_array::size_type a0 = pointsIndices[0];
    for( localIndex_array::size_type a=0 ; a<(n-2) ; ++a )
    {
      const localIndex_array::size_type a1 = pointsIndices[a+1];
      const localIndex_array::size_type a2 = pointsIndices[a+2];

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
      GEOS_ERROR("GeometryUtilities.cpp::Centroid_3DPolygon(): zero area calculated!!\n");
    }
  }
  else if( n==1 )
  {
    center = pointReferences[0];
    center += pointDisplacements[0];
  }
  else if( n==2 )
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

real64 HexVolume( array<R1Tensor> const & X )
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


}
} /* namespace geosx */
