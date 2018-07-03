/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * GeometryUtilities.cpp
 *
 *  Created on: Nov 19, 2011
 *      Author: settgast1
 */

#include "GeometryUtilities.h"

namespace GeometryUtilities
{

/**
 * @author Scott Johnson
 * @brief Calculates the polygonal area for a convex polygon with ordered points
 *
 * @param[in] pointsLocalProjection Ordered points of the convex polygon in
 * local 2D coordinates
 * @return Area of the polygon
 */
realT Area_2DPolygon(const array<R1TensorT<2> >& pointsLocalProjection)
{
  realT area = 0.;
  if(pointsLocalProjection.size() < 3)
    return area;
  //now, we have the list of vertices in order in pointsLocalProjection,
  //so get the area
  R1TensorT<2> uvi, uvip;
  for(int i = 1 ; i < (static_cast<int>(pointsLocalProjection.size())-1) ; ++i)
  {
    int ip = i+1;
    //get the cross product of uv(i+1)-uv(0) and uv(i) - uv(0)
    //uv(i)-uv(0)
    uvi = pointsLocalProjection[i];
    uvi -= pointsLocalProjection[0];

    uvip = pointsLocalProjection[ip];
    uvip -= pointsLocalProjection[0];

    //add magnitude of cross
    const realT tmp = uvip[0]*uvi[1] - uvi[0]*uvip[1];
    area += fabs(tmp);
  }

  //area is half of the stored value (1/2 |cross|)
  area *= 0.5;
  return area;
}


/**
 * @brief Finds the volume and centroid of the tetrahedron
 * @author Scott Johnson
 *
 * @param[in] x0 Point 0 on the tetrahedron
 * @param[in] x1 Point 1 on the tetrahedron
 * @param[in] x2 Point 2 on the tetrahedron
 * @param[in] x3 Point 3 on the tetrahedron
 * @param[out] x Centroid of the tetrahedron
 * @return Volume of the tetrahedron
 */
realT CentroidAndVolume_3DTetrahedron(const R1Tensor& x0, const R1Tensor& x1, const R1Tensor& x2, const R1Tensor& x3, R1Tensor& x)
{
  R1Tensor dx1(x1), dx2(x2), dx3(x3);
  dx1 -= x0;
  dx2 -= x0;
  dx3 -= x0;
  x.Cross(dx2, dx3);
  const realT volume = fabs(Dot(dx1, x)) / 6.0;
  x = dx1;
  x += dx2;
  x += dx3;
  x *= 0.25;
  x += x0;
  return volume;
}


/**
 * @author Scott Johnson
 * Calculates the centroid of a convex 2D polygon
 * @param[in] points 2D point list in order (CW or CCW) about the polygon loop
 * @param[out] center 2D center of the given ordered polygon point list
 * @return area of the convex 2D polygon
 */
realT Centroid_2DPolygon( const array<R1TensorT<2> >& points, R1TensorT<2>& center )
{
  R1TensorT<2> v0,v1,vc;
  const array<R1TensorT<2> >::size_type n = points.size();
  realT area = 0.0;
  center = 0.0;

  if (n > 2)
  {
    for (array<R1TensorT<2> >::size_type a = 0 ; a < (n - 2) ; ++a)
    {
      //get the directional vectors
      v0 = points[a + 1];
      v1 = points[a + 2];

      //get the centroid of the triangular patch
      vc = points[0];
      vc += v0;
      vc += v1;
      vc /= 3.0;

      v0 -= points[0];
      v1 -= points[0];

      //get the area of the quadrilateral; this will be halved later
      const realT area_a = v0(0) * v1(1) - v1(0) * v0(1);

      //weight the centroid appropriately
      vc *= area_a;

      //add to the sum of areas
      area += area_a;

      //add the contribution to the centroid
      center += vc;
    }

    if (area > 0.0)
    {
      //normalize the centroid
      center /= area;

      //halve the area (currently, sum of quad areas)
      area *= 0.5;
    }
    else
    {
      center = 0.0;
      for (lArray1d::size_type a = 0 ; a < n ; ++a)
      {
        center += points[a];
      }
      std::cout << "Randy's bug: area = " << area << std::endl;
      for (lArray1d::size_type a = 0 ; a < n ; ++a)
        std::cout << points[a](0) << " " << points[a](1) << std::endl;
      throw GPException("GeometryUtilities.cpp::122\n");
      center /= n;
    }
  }
  else if (n == 1)
  {
    center = points[0];
  }
  else if (n == 2)
  {
    center = points[0];
    center += points[1];
    center *= 0.5;
  }
  return area;
}

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
realT Centroid_3DPolygon(const lArray1d& pointsIndices,
                         const array<R1Tensor>& points,
                         R1Tensor& center,
                         R1Tensor& normal )
{
  R1Tensor v1,v2,vc;
  const lArray1d::size_type n = pointsIndices.size();
  realT area = 0.0;
  center = 0.0;

  if( n>2 )
  {
    const R1Tensor& x0 = points[pointsIndices[0]];
    for( lArray1d::size_type a=0 ; a<(n-2) ; ++a )
    {
      v1  = points[pointsIndices[a+1]];
      v2  = points[pointsIndices[a+2]];

      vc  = x0;
      vc += v1;
      vc += v2;

      v1 -= x0;
      v2 -= x0;

      normal.Cross(v1,v2);
      const realT triangleArea = normal.Normalize();
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
      for( lArray1d::size_type a=0 ; a<n ; ++a )
      {
        center += points[a];
      }
      std::cout << "Randy's bug: area = " << area << std::endl;
      for( lArray1d::size_type a=0 ; a<n ; ++a )
        std::cout << points[pointsIndices[a]](0) << " "
                  << points[pointsIndices[a]](1) << " "
                  << points[pointsIndices[a]](2) << " "
                  << pointsIndices[a] << std::endl;
      throw GPException("GeometryUtilities.cpp::205\n");
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
realT Centroid_3DPolygon(const lArray1d& pointsIndices,
                         const array<R1Tensor>& pointReferences,
                         const array<R1Tensor>& pointDisplacements,
                         R1Tensor& center,
                         R1Tensor& normal )
{
  R1Tensor v1,v2,vc;
  const lArray1d::size_type n = pointsIndices.size();
  realT area = 0.0;
  center = 0.0;

  if( n==3 )
  {
    const lArray1d::size_type a0 = pointsIndices[0];
    const lArray1d::size_type a1 = pointsIndices[1];
    const lArray1d::size_type a2 = pointsIndices[2];

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
    const realT triangleArea = normal.Normalize();
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
    const lArray1d::size_type a0 = pointsIndices[0];
    for( lArray1d::size_type a=0 ; a<(n-2) ; ++a )
    {
      const lArray1d::size_type a1 = pointsIndices[a+1];
      const lArray1d::size_type a2 = pointsIndices[a+2];

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
      const realT triangleArea = normal.Normalize();
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
      throw GPException("GeometryUtilities.cpp::Centroid_3DPolygon(): zero area calculated!!\n");
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

/**
 * @author annavarapusr1
 * Calculates the normal of a convex 3D polygon in the initial configuration
 * @param[in] pointIndices list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] pointReferences 3D reference point list
 * @param[out] normal Normal to the face
 * @return area of the convex 3D polygon
 */
void FaceNormal3DPolygonReferenceConfig(const lArray1d& pointsIndices,
                                        const array<R1Tensor>& pointReferences,
                                        R1Tensor& normal )
{
  R1Tensor v1,v2;
  const lArray1d::size_type n = pointsIndices.size();
  realT area = 0.0;

  if( n==3 )
  {
    const lArray1d::size_type a0 = pointsIndices[0];
    const lArray1d::size_type a1 = pointsIndices[1];
    const lArray1d::size_type a2 = pointsIndices[2];

    v1  = pointReferences[a1];
    v2  = pointReferences[a2];

    v1 -= pointReferences[a0];
    v2 -= pointReferences[a0];

    normal.Cross(v1,v2);
    const realT triangleArea = normal.Normalize();
    area += triangleArea;
    area *= 0.5;

  }
  else if( n==4 )
  {
    v1  = pointReferences[pointsIndices[3]];
    R1Tensor x3_x1(v1);

    v1  = pointReferences[pointsIndices[1]];
    x3_x1 -= v1;

    v1  = pointReferences[pointsIndices[2]];
    R1Tensor x2_x0(v1);

    v1  = pointReferences[pointsIndices[0]];
    x2_x0 -= v1;

    normal.Cross( x2_x0, x3_x1 );

    area = 0.5 * normal.Normalize();
  }
  else if( n>4 )
  {
    const lArray1d::size_type a0 = pointsIndices[0];
    for( lArray1d::size_type a=0 ; a<(n-2) ; ++a )
    {
      const lArray1d::size_type a1 = pointsIndices[a+1];
      const lArray1d::size_type a2 = pointsIndices[a+2];

      v1  = pointReferences[a1];
      v2  = pointReferences[a2];

      v1 -= pointReferences[a0];
      v2 -= pointReferences[a0];

      normal.Cross(v1,v2);
      const realT triangleArea = normal.Normalize();
      area += triangleArea;
    }
    if(area > 0.0)
    {
      area *= 0.5;
    }
    else
    {
      throw GPException("GeometryUtilities.cpp::Centroid_3DPolygon(): zero area calculated!!\n");
    }
  }
  else if( n==2 )
  {

    //For 2D elements, a face is actually an edge with two nodes. We treat the
    // length of this edge as the surface area and use it in the calculation of
    // tractions.
    R1Tensor x1_x0;
    x1_x0 = pointReferences[pointsIndices[1]];

    x1_x0 -= pointReferences[pointsIndices[0]];

    area = Dot(x1_x0, x1_x0);
    area = sqrt(area);

    x1_x0[2] = 0.0;
    x1_x0.Normalize();

    normal[0] = -x1_x0[1];
    normal[1] = x1_x0[0];
    normal[2] = 0.0;
  }
}

/**
 * Calculates the centroid of a convex 2D polygon as well as the normal
 * @param[in] pointsEntries list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] pointsCommonPlanes 3D point list
 * @param[out] center 2D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 2D polygon
 */
realT Centroid_3DPolygon( const lArray1d& pointsEntries,
                          const array<R1Tensor>& pointsCommonPlanes,
                          R1Tensor& center )
{
  R1Tensor normal;
  return Centroid_3DPolygon( pointsEntries,pointsCommonPlanes,center,normal );
}

/**
 * Calculates the centroid of a convex 2D polygon as well as the normal
 * @param[in] pointsEntries list of index references for the points array in
 * order (CW or CCW) about the polygon loop
 * @param[in] pointsCommonPlanes 3D point list
 * @param[out] center 2D center of the given ordered polygon point list
 * @param[out] normal Normal to the face
 * @return area of the convex 2D polygon
 */
realT Centroid_3DPolygon( const lArray1d& pointsEntries,
                          const array<R1Tensor>& pointReferences,
                          const array<R1Tensor>& pointDisplacements,
                          R1Tensor& center )
{
  R1Tensor normal;
  return Centroid_3DPolygon( pointsEntries,pointReferences,pointDisplacements,center,normal );
}



void FindProjectionInParentSpace ( const R1Tensor& xPoint,
                                   const R1Tensor& normal,
                                   const array<R1Tensor>& xNodes,
                                   R1Tensor& S,
                                   realT N[4],
                                   const bool initialGuess,
                                   const realT tol)
{
  R1Tensor x1_x0(xNodes[1]);
  x1_x0 -= xNodes[0];
  x1_x0 *= 0.25;

  R1Tensor x2_x3(xNodes[2]);
  x2_x3 -= xNodes[3];
  x2_x3 *= 0.25;

  R1Tensor x3_x0(xNodes[3]);
  x3_x0 -= xNodes[0];
  x3_x0 *= 0.25;

  R1Tensor x2_x1(xNodes[2]);
  x2_x1 -= xNodes[1];
  x2_x1 *= 0.25;

  R1Tensor x0px3(xNodes[0]);
  x0px3 += xNodes[3];
  x0px3 *= 0.25;

  R1Tensor x1px2(xNodes[1]);
  x1px2 += xNodes[2];
  x1px2 *= 0.25;

  R1Tensor xcf(x0px3);
  xcf += x1px2;

  realT Gtol = 0.0;

  for( unsigned int a=0 ; a<4 ; ++a )
  {
    const realT maxX = xNodes(a).MaxVal();
    if( maxX > Gtol )
      Gtol = maxX;
  }
  Gtol = (Gtol*tol)*(Gtol*tol);

  R1Tensor G;
  R2Tensor J;
  R2Tensor invJ;
  R1Tensor dS;
  R1Tensor ndS;

  const R1Tensor b(xPoint);


  if( !initialGuess )
  {
    // initial guess xi,eta,t = (0,0,0)

    // calculate residual
    G = xcf;
    G -= b;

    // calculate jacobian
    J.FillColumn(  0, x1_x0 );
    J.AddToColumn( 0, x2_x3 );
    J.FillColumn(  1, x3_x0 );
    J.AddToColumn( 1, x2_x1 );
    J.FillColumn(  2, normal );

    // solve for increment in solution
    invJ.Inverse(J);
    dS.AijBj(invJ,G);

    // ADD the increment. This is a subtract because the dS should have been
    // negative.
    S -= dS;
  }


  // now set up the iteration
  R1Tensor temp1;
  unsigned int count=0;
//    bool converged1 = false;
//    bool converged2 = false;
  const unsigned int maxiter = 5;
  for( count=0 ; count<maxiter ; ++count )
  {
    // extract the solution variables from the vector for easier use
    const realT& xi  = S[0];
    const realT& eta = S[1];
    const realT& t   = S[2];

    // construct the residual
    temp1  = x3_x0;
    temp1 *= eta;
    temp1 += x0px3;
    temp1 *= 1.0 - xi;
    G = temp1;

    temp1  = x2_x1;
    temp1 *= eta;
    temp1 += x1px2;
    temp1 *= 1.0 + xi;
    G += temp1;

    G -= b;
    temp1 = normal;
    temp1 *= t;
    G += temp1;


    if( Dot(G,G) < Gtol )
    {
//        converged1 = true;
//        if( converged2 )
      break;
    }

    // construct the jacobian

    temp1  = x2_x3;
    temp1 -= x1_x0;
    temp1 *= eta;
    temp1 += x1_x0;
    temp1 += x2_x3;
    J.FillColumn( 0, temp1 );

    temp1  = x2_x1;
    temp1 -= x3_x0;
    temp1 *= xi;
    temp1 += x3_x0;
    temp1 += x2_x1;
    J.FillColumn( 1, temp1 );


    J.FillColumn(  2, normal );

    // solve for increment in solution
    invJ.Inverse(J);
    dS.AijBj(invJ,G);

    // ADD the increment. This is a subtract because the dS should have been
    // negative.
    S -= dS;

    if( fabs(dS[0])<tol && fabs(dS[1])<tol && fabs(dS[2])<tol )
    {
//        converged2 = true;
//        if( converged1 )
      break;
    }

  }


  if( count == maxiter )
  {
    std::cout.precision(15);
    std::cout<<"  Parent coordinate inversion did not converge after iteration "<<count<<std::endl;
    //    std::cout<<"  Force = "<<force<<std::endl;
    //    std::cout<<"  xap = "<<xap<<std::endl;
    std::cout<<"    "<<S<<" | "<<dS<<" | "<<G<<" | "<<Dot(G,G)<<" | "<<Gtol<<std::endl;
  }
//    else
//    {
//    std::cout.precision(15);
//    std::cout<<"  FaceID "<<faceIndex<<": tolerances reached on iteration
// "<<count<<std::endl;
//    std::cout<<"  Force = "<<force<<std::endl;
//    std::cout<<"  xap = "<<xap<<std::endl;

  //     std::cout<<count<<"    "<<S<<" | "<<dS<<" | "<<G<<" | "<<Dot(G,G)<<" |
  // "<<Gtol<<std::endl;
//    }

  const realT xi  = S[0];
  const realT eta = S[1];
  const realT xi_x_eta = xi * eta;
  N[0] = 0.25 * (1.0 - xi - eta + xi_x_eta );
  N[1] = 0.25 * (1.0 + xi - eta - xi_x_eta );
  N[2] = 0.25 * (1.0 + xi + eta + xi_x_eta );
  N[3] = 0.25 * (1.0 - xi + eta - xi_x_eta );



}



/**
 * @brief Project a point onto the given line segment
 * @author Scott Johnson
 *
 * @param[in] point Point to project
 * @param[in] x0 Terminus 0 on the segment
 * @param[in] x1 Terminus 1 on the segment
 * @param[out] ndist Unsigned distance of original point normal to the segment
 * @param[out] udist Signed distance along x0-x1 from x0
 * @param[out] segmentLength Length of line segment x0-x1
 * @param[out] pointproj Projected point
 */
void ProjectPointToLineSegment( const R1Tensor& x0,
                                const R1Tensor& x1,
                                const R1Tensor& point,
                                realT& ndist,
                                realT& udist,
                                realT& segmentLength,
                                R1Tensor& pointproj)
{
  R1Tensor dx(point);
  dx -= x0;

  //get u-vector
  R1Tensor u(x1);
  u -= x0;
  segmentLength = u.Normalize();

  //get udist, which is the distance from x0 to the projection point
  udist = Dot(dx, u);

  //get pointproj
  pointproj = u;
  pointproj *= udist;
  pointproj += x0;

  //get ndist
  dx = point;
  dx -= pointproj;
  ndist = dx.L2_Norm();
}

realT ProjectPointToPlane( const R1Tensor& point,
                           const R1Tensor& xplane,
                           const R1Tensor& nplane,
                           R1Tensor& pointproj)
{
  //get the projection of the point onto the plane
  //along the plane's normal
  pointproj = point;
  pointproj -= xplane;
  R1Tensor dxn(nplane);
  const realT dn = Dot(nplane, pointproj);
  dxn *= dn;
  pointproj -= dxn;
  pointproj += xplane;
  return dn;
}

realT OrthogonalVectorComponent(  const R1Tensor& vtarget,
                                  const R1Tensor& vorthogonal,
                                  R1Tensor& vresidual)
{
  R1Tensor vtmp(vorthogonal);
  realT proj = Dot(vorthogonal, vtarget);
  vtmp *= proj;
  vresidual = vtarget;
  vresidual -= vtmp;
  return proj;
}

/**
 * @brief Project a point along a given vector into the given plane
 * @author Scott Johnson
 *
 * @param[in] point Point to project
 * @param[in] xplane Point in the plane
 * @param[in] nplane Unit normal to the plane
 * @param[in] vproj Unit vector to project the point along
 * @param[out] pointproj Projected point
 */
void ProjectPointToPlaneAlongUnitVector( const R1Tensor& point,
                                         const R1Tensor& xplane,
                                         const R1Tensor& nplane,
                                         const R1Tensor& vproj,
                                         R1Tensor& pointproj)
{
  //get the projection of the point onto the plane
  //along the plane's normal
  const realT dn = ProjectPointToPlane(point, xplane, nplane, pointproj);

  //handle the case of the point being in the plane already
  if( isZero(dn) )
  {
    return;
  }

  //get the projection of vproj onto the normal
  R1Tensor vn(nplane);
  const realT dnv = Dot(vproj, nplane);
  vn *= dnv;

  //make sure the vproj is not perpendicular to the normal
  if( isEqual(dnv,0) )
    throw GPException("Cannot find the projection of a point onto a line when the projection vector is parallel to the surface!");

  //get the projection of the vproj in the plane
  R1Tensor u = vproj;
  u -= vn;

  //scale it by (the normal distance of the point / normal component of the
  // vproj)
  //be careful with signs here:
  //if dn > 0 and dnv > 0 -> translate along -u
  //if dn > 0 and dnv < 0 -> translate along +u
  //if dn < 0 and dnv > 0 -> translate along +u
  //if dn < 0 and dnv < 0 -> translate along -u
  //i.e., -dn/dnv
  u *= -dn / dnv;

  //add it to the point in the plane
  pointproj += u;
}

/**
 * @date September 15, 2011
 * @author Scott Johnson
 * @brief Line intersection finds the intersection (if it exists) between two 2D
 * line segments
 *
 * Parallel segments as well as degenerate segments (i.e., collocated termini)
 * are handled in the logic
 * The point of intersection for the skew case is determined using dot products
 * of the normal component
 *
 * @param[in] uv1 The u-v coordinates for the first node of the segment on 1
 * @param[in] uv1p The u-v coordinates for the second node of the segment on 1
 * @param[in] uv2 The u-v coordinates for the first node of the segment on 2
 * @param[in] uv2p The u-v coordinates for the second node of the segment on 2
 * @param[out] intersection The u-v coordinates of the intersection (if it
 * exists)
 * @param[out] intersection2 The u-v coordinates of the second intersection (if
 * it exists ... only for parallel segment overlap)
 * @return Contact state: 0 if skew intersection, 1 if point to line
 * intersection, 2 if parallel segment overlap, 3 if not intersected
 */
localIndex LineIntersection( const R1TensorT<2>& uv1,
                             const R1TensorT<2>& uv1p,
                             const R1TensorT<2>& uv2,
                             const R1TensorT<2>& uv2p,
                             R1TensorT<2>& intersection,
                             R1TensorT<2>& intersection2,
                             const realT small)
{
  //segment 1 relative to its first node
  R1TensorT<2> uvl1(uv1p);
  uvl1 -= uv1;
  const realT dd_uvl1 = Dot(uvl1, uvl1);

  R1TensorT<2> uv21(uv2);
  uv21 -= uv1;
  R1TensorT<2> uv2p1(uv2p);
  uv2p1 -= uv1;

  //this is a point-to-segment contact where segment 1 is a point
  if( isZero(dd_uvl1) )
  {
    intersection = uv1;
    uv21.Normalize();
    uv2p1.Normalize();
    const realT tt = Dot(uv21, uv2p1);
    return ( isEqual(tt,0) || isEqual(tt,-1) ) ? 1 : 3;   // point or no contact
  }

  //normal to segment 1
  R1TensorT<2> normal1;
  {
    normal1(0) = -uvl1(1);
    normal1(1) = uvl1(0);
    normal1.Normalize();
  }

  //get the relative positions of segment 2's nodes to
  //segment 1's first node projected along the normal to 1
  realT dnuv21 = 0., dnuv2p1 = 0., dd = 0.0;
  {
    dnuv21 = Dot(normal1, uv21);
    dnuv2p1 = Dot(normal1, uv2p1);

    //if the nodes of 2 projected along 1's normal
    //do not lie on opposite sides of 1, then there
    //is no intersection
    if((dnuv21 * dnuv2p1) > 0)
      return 3;   // no contact

    //get the position of the intersection along line 2
    //as a percentage of the segment length
    dd = dnuv2p1 - dnuv21;

    //if the normal projection of segment 2 has no length, the segments are
    // parallel
    //now, we need to check cases and see whether the segments overlap
    if( isEqual(dd, 0) )
    {
      //if the segments are parallel but not on the same line, there is no
      // contact
      if( !isEqual(dnuv2p1,0) )
        return 3;

      //otherwise, the segments are parallel and on the same line,
      //but do the segments overlap at least a point on the line?
      const realT tp = Dot(uv2p1, uvl1) / dd_uvl1;
      const realT t =  Dot(uv21, uvl1)  / dd_uvl1;
      realT pa = 0, pb = 1;
      if(t > tp)
      {
        if(t<1)
          pb = t;
        if(tp>0)
          pa = tp;
      }
      else
      {
        if(tp<1)
          pb = tp;
        if(t>0)
          pa = t;
      }
      if(pb<0)
        return 3;          // no contact
      if(pa>1)
        return 3;          // no contact

      intersection = uvl1;
      intersection *= pa;
      intersection += uv1;

      intersection2 = uvl1;
      intersection2 *= pb;
      intersection2 += uv1;

      return 2;   // segment contact
    }  //END: dd == 0
  }

  //calculate intersection
  dd = fabs(dnuv21/dd);
  intersection = uv2p;
  intersection -= uv2;
  intersection *= dd;
  intersection += uv2;

  //project intersection onto line segment 1
  {
    R1TensorT<2> tmp(intersection);
    tmp -= uv1;
    dd = Dot(uvl1, tmp);
    dd /= dd_uvl1;

    //dd now holds the squared distance between node 1's first node and the
    // intersection
    //divided by the squared length of segment 1, so now check to make sure the
    //(normalized) distance (squared) is within tolerance ... i.e., the
    // intersection point
    //is not allowed to be too close to either end of segment 1
    return (dd > small && dd < (1. -  small)) ? 0 : 3;   // skew or no contact
  }
}


/**
 * @brief Determine whether the 2D point is in the 2D face defined by ordered
 * points
 * @author Scott Johnson
 *
 * The point in face is the cumulative OR for point in triangles composing the
 * face
 *
 * @param[in] uv List of 2D points on the face ordered CW
 * @param[in] uv0 2D point to test
 * @return flag true if inside
 */
bool PointInFace(const array<R1TensorT<2> >& uv, const R1TensorT<2>& uv0, const realT tol)
{
  array<R1TensorT<2> >::size_type sz = uv.size();
  if(sz < 3)
    throw GPException("Cannot define a face with less than 3 vertices");
  else
    --sz;
  for(array<R1TensorT<2> >::size_type ii = 1 ; ii < sz ; ++ii)
    if(PointInTriangle(uv[0], uv[ii], uv[ii+1], uv0, tol))
      return true;
  return false;
}

bool PointInTetrahedron(const R1Tensor& x,
                        const R1Tensor& x0,
                        const R1Tensor& x1,
                        const R1Tensor& x2,
                        const R1Tensor& x3)
{
  R1Tensor center(x0);
  center += x1;
  center += x2;
  center += x3;
  center *= 0.25;

  R1Tensor normal(0);

  R1Tensor xc0(x0), x01(x1), x02(x2), x0p(x);
  xc0 -= center;
  if(isZero(Dot(xc0,xc0)))
    throw GPException("PointInTetrahedron: Attempt to evaluate a degenerate tetrahedron");
  x01 -= x0;
  x02 -= x0;
  x0p -= x0;

  //0,1,2
  normal.Cross(x01, x02);
  if(Dot(normal,xc0) < 0)
    normal *= -1;
  if(Dot(normal,x0p) > 0)
    return false;

  R1Tensor x03(x3);
  x03 -= x0;
  //0,1,3
  normal.Cross(x01, x03);
  if(Dot(normal,xc0) < 0)
    normal *= -1;
  if(Dot(normal,x0p) > 0)
    return false;

  //0,2,3
  normal.Cross(x02, x03);
  if(Dot(normal,xc0) < 0)
    normal *= -1;
  if(Dot(normal,x0p) > 0)
    return false;

  //1,2,3
  R1Tensor x10(x0);
  x10 -= x1;
  x02 += x10;
  x03 += x10;
  x0p += x10;
  xc0 -= x10;
  normal.Cross(x02, x03);
  if(Dot(normal,xc0) < 0)
    normal *= -1;
  if(Dot(normal,x0p) > 0)
    return false;

  //all cases ok
  return true;
}

bool PointInPolyhedron(const R1Tensor& point,
                       const localIndex elementIndex,
                       const lArray2d& toFaces,
                       const lArray2d& toNodes,
                       const array<lArray1d>& faceToNodes,
                       const array<R1Tensor>& referencePosition,
                       const array<R1Tensor>& displacement,
                       const bool planarFaces)
{
  //Get element center
  R1Tensor xel(0.0);
  {
    R1Tensor position(0.0);
    R1Tensor xmin(std::numeric_limits<realT>::max());
    R1Tensor xmax(-std::numeric_limits<realT>::max());
    for (localIndex i = 0 ; i < toNodes.Dimension(1) ; i++)
    {
      position = referencePosition[toNodes(elementIndex, i)];
      position += displacement[toNodes(elementIndex, i)];
      xel += position;
      xmin.SetMin(position);
      xmax.SetMax(position);
    }
    xel *= 1.0 / toNodes.Dimension(1);

    //
    for (localIndex i = 0 ; i < nsdof ; i++)
    {
      if (point[i] > xmax[i] || point[i] < xmin[i])
        return false;
    }
  }

  //determine whether to reject the point
  R1Tensor center(0), normal(0);
  if(planarFaces)
  {
    for (localIndex i = 0 ; i < toFaces.Dimension(1) ; i++)
    {
      const localIndex faceIndex = toFaces(elementIndex, i);

      //make sure the normal is facing out
      const lArray1d& nodeList = faceToNodes[faceIndex];
      GeometryUtilities::Centroid_3DPolygon(nodeList, referencePosition, displacement, center,
                                            normal);
      R1Tensor xfc(center);
      xfc -= xel;
      if (Dot(normal, xfc) < 0)
        normal *= -1.0;

      //make sure the point does not lie outside of the planar approximation of
      // the face
      R1Tensor d(point);
      d -= center;
      if (Dot(normal, d) > 0)
        return false;
    }
    return true;
  }
  else
  {
    for (localIndex i = 0 ; i < toFaces.Dimension(1) ; i++)
    {
      const localIndex faceIndex = toFaces(elementIndex, i);

      const lArray1d& nodeList = faceToNodes[faceIndex];
      GeometryUtilities::Centroid_3DPolygon(nodeList, referencePosition, displacement, center,
                                            normal);

      //go through each edge, and use the face center and element center to
      // construct tets
      for(localIndex iNd = 0 ; iNd < nodeList.size() ; iNd++)
      {
        const localIndex ib = nodeList[iNd];
        const localIndex ia = iNd > 0 ? nodeList[iNd-1] : nodeList.back();

        R1Tensor xb(referencePosition[ib]);
        xb += displacement[ib];
        R1Tensor xa(referencePosition[ia]);
        xa += displacement[ia];

        if(PointInTetrahedron(point, xb, xa, center, xel))
          return true;
      }
    }

    //if the point does not lie in any of the constituent tets, then it is
    // outside
    return false;
  }
}

/**
 * @brief Determine whether the 2D point is in the 2D triangle
 * @author Scott Johnson
 *
 * This is via a barycentric technique: e.g.,
 * http://www.blackpawn.com/texts/pointinpoly/default.html
 *
 * @param[in] uva First point on the triangle
 * @param[in] uvb Second point on the triangle
 * @param[in] uvc Third point on the triangle
 * @param[in] uv0 2D point to test
 * @return flag true if inside
 */
bool PointInTriangle(
  const R1TensorT<2>& uva, const R1TensorT<2>& uvb,
  const R1TensorT<2>& uvc, const R1TensorT<2>& uv0,
  const realT tol)
{
  // Compute vectors
  R1TensorT<2> v0(uvc);
  v0 -= uva;
  R1TensorT<2> v1(uvb);
  v1 -= uva;
  R1TensorT<2> v2(uv0);
  v2 -= uva;

  // Compute dot products
  const realT dot00 = Dot(v0, v0);
  const realT dot01 = Dot(v0, v1);
  const realT dot02 = Dot(v0, v2);
  const realT dot11 = Dot(v1, v1);
  const realT dot12 = Dot(v1, v2);

  // Compute barycentric coordinates
  realT invDenom = (dot00 * dot11 - dot01 * dot01);
  if( isZero(invDenom) )
    return false;
  invDenom = 1.0 / invDenom;
  const realT u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  const realT v = (dot00 * dot12 - dot01 * dot02) * invDenom;

  // Check if point is in triangle
  return (u > -tol) && (v > -tol) && (u + v < (1+tol));
}

/**
 * @author Scott Johnson
 * @brief Finds the polygon of intersection for two convex polygons of the same
 * handedness for the special case where the polygons are coincident
 * @param[in] a1 Ordered 2D points associated with polygon 1
 * @param[in] a2 Ordered 2D points associated with polygon 2
 * @param[out] aa Ordered 2D points associated with the polygon of intersection
 * @return Number of points in the polygon of intersection
 */
array<R1TensorT<2> >::size_type CoincidentIntersection_2DPolygons( const array<R1TensorT<2> >& a1,
                                                                   const array<R1TensorT<2> >& a2,
                                                                   array<R1TensorT<2> >& aa)
{
  array<R1TensorT<2> > uv;
  uv.resize(3,static_cast< R1TensorT<2> >(0.0) );

  bool allIn1 = true, allIn2 = false;

  //FIRST, CHECK WHETHER ALL POINTS OF 2 ARE IN 1
  for (array<R1TensorT<2> >::size_type ii = 0 ; ii < a2.size() ; ++ii)
  {
    bool found = false;
    for (int i = 2 ; i < static_cast<int>(a1.size()) ; ++i)
    {
      if (PointInTriangle(a1[i-2], a1[i-1], a1[i], a2[ii]))
      {
        found = true;
        break;
      }
    }
    if (!found)
    {
      allIn1 = false;
      break;
    }
  }
  if(!allIn1)
  {
    allIn2 = true;
    for (array<R1TensorT<2> >::size_type ii = 0 ; ii < a1.size() ; ++ii)
    {
      bool found = false;
      for (int i = 2 ; i < static_cast<int>(a2.size()) ; ++i)
      {
        if (PointInTriangle(a2[i-2], a2[i-1], a2[i], a1[ii]))
        {
          found = true;
          break;
        }
      }
      if (!found)
      {
        allIn2 = false;
        break;
      }
    }
  }

  //determine whether 1 is in 2 or vice versa
  aa.clear();
  if (allIn1)
  {
    //all the points of polygon 2 are contained in 1
    //so populate the list with the points in 2
    for (array<R1TensorT<2> >::size_type i = 0 ; i < a2.size() ; ++i)
      aa.push_back(a2[i]);
  }
  else if(allIn2)
  {
    //all the points of polygon 1 are contained in 2
    //so populate the list with the points in 1
    for (array<R1TensorT<2> >::size_type i = 0 ; i < a1.size() ; ++i)
      aa.push_back(a1[i]);
  }
  return aa.size();
}



/**
 * @author Scott Johnson
 * @brief Determines the intersection points between two convex polygons of the
 * same handedness
 *
 * Intersections are determined using an O(N^2) technique comparing each segment
 * to each segment
 * Parallel segments as well as degenerate segments (i.e., collocated termini)
 * are handled in the logic
 * of LineIntersection
 *
 * @param[in] a1 Points on polygon 1
 * @param[in] a2 Points on polygon 2
 * @param[out] intersections Points where polygons 1 and 2 intersect
 */
void Intersections_2DPolygons(const array<R1TensorT<2> >& a1,
                              const array<R1TensorT<2> >& a2,
                              array<R1TensorT<2> >& intersections)
{
  //note: this is an O(N^2) algorithm ... shouldn't matter since typical faces
  // have
  //limited numbers of segments
  //intersections.clear();

  //identify intersections
  for (array<R1TensorT<2> >::size_type i1 = 1 ; i1 <= a1.size() ; ++i1)
  {
    array<R1TensorT<2> >::size_type ii1 = i1 == a1.size() ? 0 : i1;
    for (array<R1TensorT<2> >::size_type i2 = 1 ; i2 <= a2.size() ; ++i2)
    {
      array<R1TensorT<2> >::size_type ii2 = i2 == a2.size() ? 0 : i2;

      //find intersection of segments 1 and 2, otherwise go to the next
      // iteration
      R1TensorT<2> intersection, intersection2;
      localIndex contactType = LineIntersection(a1[i1-1],a1[ii1],a2[i2-1],a2[ii2],
                                                intersection, intersection2, 0.0);   // Call
                                                                                     // a
                                                                                     // different
                                                                                     // LineIntersection
      if(contactType < 3)   // no contact
      {
        intersections.push_back(intersection);
        if(contactType == 2)
          intersections.push_back(intersection2);
      }
    }  //foreach i2
  }  //foreach i1

  //return intersections.size();
}

/**
 * @brief Order the polygon points in counter-clockwise order, and remove all
 * redundant points
 * @author Scott Johnson
 *
 * Uses a generalized Graham Scan.
 * The given polygon points are ordered in CCW direction by obtaining the dot
 * product
 * with the pivot vector (in the non-pathological case, this is the unit vector
 * along 0-1,
 * but the logic admits the non-ideal case, too); for descending dot product,
 * the angular
 * displacement from the pivot increases. Because this is a convex loop, the
 * logic should hold.
 *
 * Ordering is secondarily dependent on the distance from the pivot node, such
 * that, for any set
 * of nodes with the same angular displacement, that with the largest distance
 * is kept, while the
 * others are removed. This produces a minimal description of the convex
 * polygon.
 *
 * @param[in,out] polygonPoints List of unordered points coming in, reordered
 * out
 */
void OrderPoints_2DPolygon(array<R1TensorT<2> >& polygonPoints, const realT tol  )
{
  if(polygonPoints.size() < 3)
    return;

  R1TensorT<2> vpivot = static_cast< R1TensorT<2> >(0.0);

  //(0) fill array of point structs
  array<PointStruct> arr;
  arr.push_back(vpivot);
  for(array<R1TensorT<2> >::size_type i = 0 ; i < polygonPoints.size() ; ++i)
  {
    arr.push_back(PointStruct(polygonPoints[i]));
    arr[0].m_point += polygonPoints[i];
  }
  arr[0].m_point /= polygonPoints.size();
  polygonPoints.clear();

  //(1) make sure the lowest coordinate point is first
  array<PointStruct>::iterator iter = arr.begin();

  //(2) get the pivot vector 0->1
  //remove any colocated points at the beginning of the array
  //and set the pivot vector as the unit vector between the first two
  //(non-colocated) points
  {
    //make sure index 0 stays at index 0
    ++iter;  //iter is now set at arr[1]

    while(iter != arr.end())
    {
      vpivot = (*iter).m_point;
      vpivot -= arr[0].m_point;
      if( isZero(Dot(vpivot, vpivot)) )
      {
        iter = arr.erase(iter);
      }
      else
      {
        ++iter;
        break;
      }
    }
    if(iter == arr.end())
      return;
    arr[1].m_value = vpivot.Normalize();
    arr[1].m_angle = 1.0;
  }   //iter is now set at arr[2]

  //(3) if the number of polygon points has been reduced below the threshold,
  // return
  if(arr.size() < 4)
    return;

  //(4) get angle values about the pivot
  {
    //get the additional vectors and find the cross products and dot products
    // (sine and cosine)
    R1TensorT<2> value;
    while(iter != arr.end())
    {
      value = (*iter).m_point;
      value -= arr[0].m_point;
      (*iter).m_value = value.Normalize();
      if( isEqual((*iter).m_value,0) )
      {
        iter = arr.erase(iter);
        continue;
      }

      //get the angle metric relative to vpivot
      realT sineValue = PointStruct::Cross(vpivot, value);
      realT cosineValue = Dot(vpivot, value);

      //the "angle" metric is in the range [1,-3) in counter-clockwise order
      // (note: at cosine=+/-1 : sine can be +/-0)
      (*iter).m_angle = sineValue >= 0 ? cosineValue : -2 - cosineValue;  //1,0,-1,-2,(-3
                                                                          // or
                                                                          // 1)

      ++iter;
    }
    if(arr.size() < 4)
      return;
  }

  //(4) sort points: counter-clockwise will go most positive to most negative,
  // starting with pivot as the first point
  //    sorting is primary on dot product and secondary on distance
  std::sort(arr.begin()+1,arr.end(), PointStruct::Compare);

  //(5) polygonPoints is now ordered, so go through and remove those with
  // duplicate angle values; note that the first
  //    (i.e., that with the largest distance value of the set of equal dot
  // product values) will be kept while the others
  //    are removed from the collection.
  /*  {
      iter = arr.begin();//polygonPoints[0]
   ++iter;//arr[1] -> corresponding to "last"
      realT last = (*iter).angle;
   ++iter;//arr[2] -> corresponding to "curr"
      while(iter != arr.end())
      {
        if(fabs((*iter).angle-last) < tol || fabs((*iter).value) < tol) {
          iter = arr.erase(iter);
        } else {
          last = (*iter).angle;
   ++iter;
        }
      }
     }*/



  //(6) refill polygonPoints with arr points
  if(arr.size() < 4)
    return;

  for(array<PointStruct>::size_type i = 1 ; i < arr.size() ; ++i)
    polygonPoints.push_back(arr[i].m_point);

  RemoveCoincidentOrderedPoints( polygonPoints, 1.0e-8 );
  //polygonPoints is now ordered in counter-clockwise order!
  return;
}

/**
 * @brief Remove redundant points
 * @author Scott Johnson
 *
 * @param[in,out] polygonPoints List of ordered points coming in and going out
 *(less the redundant points)
 */
template< int T_size >
void RemoveCoincidentOrderedPoints(array<R1TensorT<T_size> >& polygonPoints, const realT tol)
{
  if(polygonPoints.size() < 2)
    return;

  const realT tt = tol * tol;
  typename array<R1TensorT<T_size> >::iterator iter = polygonPoints.end();
  --iter;

  R1TensorT<T_size> v(*iter);
  iter = polygonPoints.begin();
  while(iter != polygonPoints.end())
  {
    v -= *iter;
    const bool rm = Dot(v,v) < tt;
    v = *iter;
    if(rm)
      iter = polygonPoints.erase(iter);
    else
      ++iter;
  }
}

template void RemoveCoincidentOrderedPoints<2>(array<R1TensorT<2> >& polygonPoints, const realT tol);
template void RemoveCoincidentOrderedPoints<3>(array<R1TensorT<3> >& polygonPoints, const realT tol);

///Calculate the area and volume of intersection of 2D polygons quickly
/**
 * @author E. Herbold
 * Translated to C++ by S. Johnson <- blame errors here
 * This does not calculate the overlap region geometry, just the area
 * NOTE: POINTS DEFINED ANTI-CLOCKWISE AROUND THE BOUNDARY
 * @param[in] a1 Local planar coordinates of the projected face 1
 * @param[in] a2 Local planar coordinates of the projected face 2
 * @return Area of intersection
 */
realT Intersection_2DPolygons(const array<R1TensorT<2> >& a1,
                              const array<R1TensorT<2> >& a2)
{
  if (a2.size() < 3 || a1.size() < 3)
    return 0.0;

  realT qy = 0.0, area = 0.0;

  //************************************************************************
  //CALCULATE ORIGIN SHIFT TO AVOID ROUNDOFF ERRORS
  R1TensorT<2> x0 = static_cast< R1TensorT<2> >( std::numeric_limits<realT>::max() );

  //get extrema and leave if no overlap
  {
    R1TensorT<2> xmin1 = static_cast< R1TensorT<2> >(std::numeric_limits<realT>::max());
    R1TensorT<2> xmax1 = static_cast< R1TensorT<2> >(-std::numeric_limits<realT>::max());
    for (localIndex i = 0 ; i < a1.size() ; ++i)
    {
      xmin1.SetMin(a1[i]);
      xmax1.SetMax(a1[i]);
      x0.SetMin(a1[i]);
    }

    R1TensorT<2> xmin2 = static_cast< R1TensorT<2> >(std::numeric_limits<realT>::max());
    R1TensorT<2> xmax2 = static_cast< R1TensorT<2> >(-std::numeric_limits<realT>::max());
    for (localIndex i = 0 ; i < a2.size() ; ++i)
    {
      xmin2.SetMin(a2[i]);
      xmax2.SetMax(a2[i]);
      x0.SetMin(a2[i]);
    }
    if (xmin1(0) > xmax2(0) || xmin2(0) > xmax1(0) || xmin1(1) > xmax2(1) || xmin2(1) > xmax1(1))
      return 0.0;
  }

  //loop over a1
  R1TensorT<2> x1, xp1, dx1;
  R1TensorT<2> x2, xp2, dx2;
  for (localIndex i1 = 0 ; i1 < a1.size() ; ++i1)
  {
    x1 = a1[i1];
    x1 -= x0;
    xp1 = a1[(i1 == (a1.size() - 1)) ? 0 : (i1 + 1)];
    xp1 -= x0;
    dx1 = xp1;
    dx1 -= x1;
    if (isEqual(dx1(0), 0.0))
      continue;
    const realT slope1 = dx1(1) / dx1(0);

    //loop over a2
    for (localIndex i2 = 0 ; i2 < a2.size() ; ++i2)
    {
      x2 = a2[i2];
      x2 -= x0;
      xp2 = a2[(i2 == (a2.size() - 1)) ? 0 : (i2 + 1)];
      xp2 -= x0;
      dx2 = xp2;
      dx2 -= x2;
      if (isEqual(dx2(0), 0.0))
        continue;
      const realT slope2 = dx2(1) / dx2(0);

      const realT ss = dx1(0) * dx2(0);

      //CALCULATE LEFT AND RIGHT COORDINATES OF OVERLAP (X)
      realT xl, xr;
      {
        xl = x1(0) < xp1(0) ? x1(0) : xp1(0);
        realT tmp = x2(0) < xp2(0) ? x2(0) : xp2(0);
        xl = xl > tmp ? xl : tmp;

        xr = x1(0) > xp1(0) ? x1(0) : xp1(0);
        tmp = x2(0) > xp2(0) ? x2(0) : xp2(0);
        xr = xr < tmp ? xr : tmp;

        if (xl >= xr)
          continue;
      }

      //CALCULATE LEFT AND RIGHT COORDINATES OF OVERLAP ()
      realT yl, yr;
      {
        const realT yl1 = x1(1) + (xl - x1(0)) * slope1;
        const realT yl2 = x2(1) + (xl - x2(0)) * slope2;
        const realT yr1 = x1(1) + (xr - x1(0)) * slope1;
        const realT yr2 = x2(1) + (xr - x2(0)) * slope2;
        yl = yl1 < yl2 ? yl1 : yl2;
        yr = yr1 < yr2 ? yr1 : yr2;
      }

      //CHECK WHETHER LINES INTERSECT
      const realT dslope = slope1 - slope2;

      if (!isEqual(dslope, 0.0))
      {
        const realT xm = (x2(1) - x1(1) + slope1 * x1(0) - slope2 * x2(0)) / dslope;
        const realT ym = x1(1) + slope1 * (xm - x1(0));
        if (xm > xl && xm < xr)
        {
          //LINES INTERSECT, CASE II
          const realT area1 = 0.5 * fabs((yl + ym) * (xm - xl)) * (ss < 0 ? -1.0 : 1.0);
          const realT area2 = 0.5 * fabs((ym + yr) * (xr - xm)) * (ss < 0 ? -1.0 : 1.0);
          area += area1;
          area += area2;

          if ((yl + ym) > 0)
            qy += (ym + yl * yl / (yl + ym)) * area1 / 3.0;
          if ((ym + yr) > 0)
            qy += (yr + ym * ym / (ym + yr)) * area2 / 3.0;
          continue;
        }
      }

      //LINES DO NOT INTERSECT, CASE I
      {
        const realT area1 = 0.5 * fabs((xr - xl) * (yr + yl)) * (ss < 0 ? -1.0 : 1.0);
        area += area1;
        if ((yl + yr) > 0)
          qy += (yr + yl * yl / (yl + yr)) * area1 / 3.0;
      }
    }
  }

  if (area < 0)
    throw GPException("Area less than 0");
  return area;
}

/**
 * @brief Determine the area of intersection and polygon of intersection of two
 * coplanar, convex polygons
 * @author Scott Johnson
 *
 * Description: this is an adaptation of the method proposed by
 * --unlimited release--
 * Didonato (1993) "An algorithm to find the intersection of two convex
 * polygons"
 * AD-A274 722 NSWC (Dahlgren Div., Dahlgren, VA 22448)
 * --unlimited release--
 * The same general workflow is used, but the method of determining
 * intersections and inclusions is different
 * as is the method of ordering the resultant point loop to define the convex
 * polygon of overlap
 * This is probably publishable as its own paper, and I will update the citation
 * when I have a chance to write it up
 *
 * @param[in] a1 Local planar coordinates of the projected face 1
 * @param[in] a2 Local planar coordinates of the projected face 2
 * @param[out] aa Local planar coordinates of the polygon of intersection
 * @return Area of intersection
 *
 * NOTE: for this to work, both a1 and a2 have to be listed in the same loop
 * direction (either CW or CCW)
 * THERE IS NO CHECKING FOR THIS CONDITION, SO BE CAREFUL
 */
realT Intersection_2DPolygons(const array<R1TensorT<2> >& a1,
                              const array<R1TensorT<2> >& a2,
                              const realT& positionTolerance,
                              array<R1TensorT<2> >& polygonPoints)
{
  //make sure the polygonPoints have been cleared
  polygonPoints.clear();

  //(1) first, if either polygon is 1D or less, there is, by definition, no
  // contact area
  if(a1.size() < 3 || a2.size() < 3)
    return 0.0;

  //(2) second, if the bounding boxes in local frame don't overlap, there is no
  // contact area
  {
    R1TensorT<2> min1 = static_cast< R1TensorT<2> >(std::numeric_limits<realT>::max());
    R1TensorT<2> min2 = static_cast< R1TensorT<2> >(std::numeric_limits<realT>::max());
    R1TensorT<2> max1 = static_cast< R1TensorT<2> >(-std::numeric_limits<realT>::max());
    R1TensorT<2> max2 = static_cast< R1TensorT<2> >(-std::numeric_limits<realT>::max());

    for(array<R1TensorT<2> >::size_type i = 0 ; i < a1.size() ; ++i)
    {
      min1.SetMin(a1[i]);
      max1.SetMax(a1[i]);
    }
    for(array<R1TensorT<2> >::size_type i = 0 ; i < a2.size() ; ++i)
    {
      min2.SetMin(a2[i]);
      max2.SetMax(a2[i]);
    }

    if (min2(0) > max1(0) || min1(0) > max2(0) || min2(1) > max1(1) || min1(1) > max2(1))
      return 0.0;
  }

  bool ret = false;

  //(3) Are all of the points of 1 contained in 2? -> YES, then the overlap is 1
  // RETURN
  //here is the contact problem
  {
    ret = true;
    for(array<R1TensorT<2> >::size_type i = 0 ; i < a1.size() ; ++i)
    {
      if(!PointInFace(a2, a1[i]))
        ret = false;
      else
        polygonPoints.push_back(a1[i]);
    }
  }

  //(4) Are all of the points of 2 contained in 1? -> YES, then the overlap is 2
  // RETURN
  if(!ret)
  {
    ret = polygonPoints.size() == 0;
    for(array<R1TensorT<2> >::size_type i = 0 ; i < a2.size() ; ++i)
    {
      if(!PointInFace(a1, a2[i]))
        ret = false;
      else
        polygonPoints.push_back(a2[i]);
    }
  }

  //(5) Are there intersections between 1 and 2
  if(!ret)
    Intersections_2DPolygons(a1, a2, polygonPoints);

  //(6) order the points ... if there are any
  OrderPoints_2DPolygon(polygonPoints, positionTolerance);

  const realT area = polygonPoints.size() < 3 ? 0.0 : Area_2DPolygon(polygonPoints);
  return area;
}   //PolygonPairIntersection


/*
 * A version of the above with modifications and bug-fixes as were apparent to
 * me when called from XfemManager
 * There is a lot of machinery in the original routines that is specific to
 * contact detection
 * For my purposes, those were unnecessary and were also resulting in errors.
 */
realT Overlap2DPolygons( const array<R1TensorT<2> >& a1,
                         const array<R1TensorT<2> >& a2,
                         const realT& positionTolerance,
                         array<R1TensorT<2> >& polygonPoints)
{
  //make sure the polygonPoints have been cleared
  polygonPoints.clear();

  //(1) first, if either polygon is 1D or less, there is, by definition, no
  // contact area
  if(a1.size() < 3 || a2.size() < 3)
    return 0.0;

  //(2) second, if the bounding boxes in local frame don't overlap, there is no
  // contact area
  {
    R1TensorT<2> min1 = static_cast< R1TensorT<2> >(std::numeric_limits<realT>::max());
    R1TensorT<2> min2 = static_cast< R1TensorT<2> >(std::numeric_limits<realT>::max());
    R1TensorT<2> max1 = static_cast< R1TensorT<2> >(-std::numeric_limits<realT>::max());
    R1TensorT<2> max2 = static_cast< R1TensorT<2> >(-std::numeric_limits<realT>::max());

    for(array<R1TensorT<2> >::size_type i = 0 ; i < a1.size() ; ++i)
    {
      min1.SetMin(a1[i]);
      max1.SetMax(a1[i]);
    }
    for(array<R1TensorT<2> >::size_type i = 0 ; i < a2.size() ; ++i)
    {
      min2.SetMin(a2[i]);
      max2.SetMax(a2[i]);
    }

    if (min2(0) > max1(0) || min1(0) > max2(0) || min2(1) > max1(1) || min1(1) > max2(1))
      return 0.0;
  }

  bool allInOne = true;
  for(array<R1TensorT<2> >::size_type i = 0 ; i < a1.size() ; ++i)
  {
    if(!PointInFace(a2, a1[i]))
    {
      allInOne = false;
      polygonPoints.clear();
      break;
    }
    else
      polygonPoints.push_back(a1[i]);
  }

  bool allInTwo = true;
  if(!allInOne)
  {
    for(array<R1TensorT<2> >::size_type i = 0 ; i < a2.size() ; ++i)
    {
      if(!PointInFace(a1, a2[i]))
      {
        allInTwo = false;
        polygonPoints.clear();
        break;
      }
      else
        polygonPoints.push_back(a2[i]);
    }
  }

  if(!allInOne && !allInTwo)
  {
    bool ret = false;

    //(3) Are all of the points of 1 contained in 2? -> YES, then the overlap is
    // 1 RETURN
    //here is the contact problem
    {
      ret = true;
      for(array<R1TensorT<2> >::size_type i = 0 ; i < a1.size() ; ++i)
      {
        if(!PointInFace(a2, a1[i]))
          ret = false;
        else
          polygonPoints.push_back(a1[i]);
      }
    }

    //(4) Are all of the points of 2 contained in 1? -> YES, then the overlap is
    // 2 RETURN
    if(!ret)
    {
      ret = polygonPoints.size() == 0;
      for(array<R1TensorT<2> >::size_type i = 0 ; i < a2.size() ; ++i)
      {
        if(!PointInFace(a1, a2[i]))
          ret = false;
        else
          polygonPoints.push_back(a2[i]);
      }
    }

    //(5) Are there intersections between 1 and 2
    if(!ret)
      Intersections2DPolygon(a1, a2, polygonPoints);
  }

  //(6) order the points ... if there are any
  OrderPoints_2DPolygon(polygonPoints, positionTolerance);

  const realT area = polygonPoints.size() < 3 ? 0.0 : Area_2DPolygon(polygonPoints);
  return area;
}  //PolygonPairIntersection

void Intersections2DPolygon( const array<R1TensorT<2> >& a1,
                             const array<R1TensorT<2> >& a2,
                             array<R1TensorT<2> >& intersections)
{
  //note: this is an O(N^2) algorithm ... shouldn't matter since typical faces
  // have
  //limited numbers of segments
  //intersections.clear();

  //identify intersections
  for (array<R1TensorT<2> >::size_type i1 = 1 ; i1 <= a1.size() ; ++i1)
  {
    array<R1TensorT<2> >::size_type ii1 = i1 == a1.size() ? 0 : i1;
    for (array<R1TensorT<2> >::size_type i2 = 1 ; i2 <= a2.size() ; ++i2)
    {
      array<R1TensorT<2> >::size_type ii2 = i2 == a2.size() ? 0 : i2;

      //find intersection of segments 1 and 2, otherwise go to the next
      // iteration
      R1TensorT<2> intersection;
      localIndex contactType = LineIntersection(a1[i1-1],a1[ii1],a2[i2-1],a2[ii2],
                                                intersection, 0.0);
      if(contactType < 3)   // no contact
      {
        intersections.push_back(intersection);
      }
    }  //foreach i2
  }  //foreach i1

  //return intersections.size();
}

/*
 * The routine below will rule out segment to segment overlap cases
 * This is being used only in calls to Intersection_2DPolygons by XfemManager to
 * calculate overlapping area and intersection points between polygons
 */
localIndex LineIntersection( const R1TensorT<2>& uv1,
                             const R1TensorT<2>& uv1p,
                             const R1TensorT<2>& uv2,
                             const R1TensorT<2>& uv2p,
                             R1TensorT<2>& intersection,
                             const realT small)
{
  //segment 1 relative to its first node
  R1TensorT<2> uvl1(uv1p);
  uvl1 -= uv1;
  const realT dd_uvl1 = Dot(uvl1, uvl1);

  R1TensorT<2> uv21(uv2);
  uv21 -= uv1;
  R1TensorT<2> uv2p1(uv2p);
  uv2p1 -= uv1;

  //this is a point-to-segment contact where segment 1 is a point
  if( isZero(dd_uvl1) )
  {
    intersection = uv1;
    uv21.Normalize();
    uv2p1.Normalize();
    const realT tt = Dot(uv21, uv2p1);
    return ( isEqual(tt,0) || isEqual(tt,-1) ) ? 1 : 3;   // point or no contact
  }

  //normal to segment 1
  R1TensorT<2> normal1;
  {
    normal1(0) = -uvl1(1);
    normal1(1) = uvl1(0);
    normal1.Normalize();
  }

  //get the relative positions of segment 2's nodes to
  //segment 1's first node projected along the normal to 1
  realT dnuv21 = 0., dnuv2p1 = 0., dd = 0.0;
  {
    dnuv21 = Dot(normal1, uv21);
    dnuv2p1 = Dot(normal1, uv2p1);

    //if the nodes of 2 projected along 1's normal
    //do not lie on opposite sides of 1, then there
    //is no intersection
    if((dnuv21 * dnuv2p1) > 0)
      return 3;   // no contact

    //get the position of the intersection along line 2
    //as a percentage of the segment length
    dd = dnuv2p1 - dnuv21;

    //if the normal projection of segment 2 has no length, the segments are
    // parallel
    //now, we need to check cases and see whether the segments overlap
    if( std::fabs(dd)<=1e-12 )
    {
      return 3;
    }  //END: dd == 0
  }

  //calculate intersection
  dd = fabs(dnuv21/dd);
  intersection = uv2p;
  intersection -= uv2;
  intersection *= dd;
  intersection += uv2;

  //project intersection onto line segment 1
  {
    R1TensorT<2> tmp(intersection);
    tmp -= uv1;
    dd = Dot(uvl1, tmp);
    dd /= dd_uvl1;

    //dd now holds the squared distance between node 1's first node and the
    // intersection
    //divided by the squared length of segment 1, so now check to make sure the
    //(normalized) distance (squared) is within tolerance ... i.e., the
    // intersection point
    //is not allowed to be too close to either end of segment 1
    return (dd > small && dd < (1. -  small)) ? 0 : 3;   // skew or no contact
  }
}

/**
 * @brief Find two orthogonal unit vectors that lie in a given plane
 * @author Scott Johnson
 * Determine a set of unit vectors in the plane given the normal
 * @param[in] normal Normal to the plane
 * @param[out] e1 first unit vector
 * @param[out] e2 second unit vector
 * @return Flag whether the given normal is pathologically zero
 */
bool VectorsInPlane(const R1Tensor& normal, R1Tensor& e1, R1Tensor& e2)
{
  if( isZero(Dot(normal,normal)) )
    return false;
  R1Tensor tmp = static_cast<R1Tensor>(0.0);
  tmp(0) = 1.0;
  e1.Cross(normal,tmp);
  if( isZero(Dot(e1,e1)) )
  {
    tmp(1) = 1.0;
    tmp(0) = 0.0;
    e1.Cross(normal,tmp);
  }
  e1.Normalize();
  e2.Cross(e1,normal);
  return true;
}


/**
 * @brief Project a list of points onto a plane
 * @author Scott Johnson
 * @param[in] pointsEntries list of indices in the point references and
 * displacements list
 * @param[in] pointReferences list of reference positions
 * @param[in] pointDisplacements list of displacements of the point from its
 * reference position
 * @param[in] e1 unit vector in the plane
 * @param[in] e2 unit vector in the plane orthogonal to e1
 * @param[out] localPoints coordinates of the given points in the e1 and e2
 * coordinates
 * @param[out] min minima of the point list in the e1 and e2 coordinates
 * @param[out] max maxima of the point list in the e1 and e2 coordinates
 * @return return
 */
void CartesianPointsProjectedToPlanarPoints(const lArray1d& pointsEntries,
                                            const array<R1Tensor>& pointReferences,
                                            const array<R1Tensor>& pointDisplacements,
                                            const R1Tensor& e1,
                                            const R1Tensor& e2,
                                            array<R1TensorT<2> >& localPoints,
                                            R1TensorT<2>& min,
                                            R1TensorT<2>& max)
{
  localPoints.clear();
  localPoints.reserve(pointsEntries.size());
  min = std::numeric_limits<realT>::max();
  max = -std::numeric_limits<realT>::max();
  for(lArray1d::size_type a = 0 ; a < pointsEntries.size() ; ++a)
  {
    const localIndex b = pointsEntries[a];
    R1TensorT<2> tmp;
    R1Tensor gtmp(pointReferences[b]);
    gtmp += pointDisplacements[b];
    for(localIndex i = 0 ; i < 2 ; ++i)
    {
      tmp(i) = Dot(gtmp, i == 0 ? e1 : e2);
      if(tmp(i) < min(i))
        min(i) = tmp(i);
      if(tmp(i) > max(i))
        max(i) = tmp(i);
    }
    localPoints.push_back(tmp);
  }
}

/**
 * @brief Project a list of points onto a plane
 * @author Scott Johnson
 * @param[in] points list of point positions in Cartesian space
 * @param[in] e1 unit vector in the plane
 * @param[in] e2 unit vector in the plane orthogonal to e1
 * @param[out] localPoints coordinates of the given points in the e1 and e2
 * coordinates
 * @param[out] min minima of the point list in the e1 and e2 coordinates
 * @param[out] max maxima of the point list in the e1 and e2 coordinates
 * @return return
 */
void CartesianPointsProjectedToPlanarPoints(  const array<R1Tensor>& points,
                                              const R1Tensor& e1,
                                              const R1Tensor& e2,
                                              array<R1TensorT<2> >& localPoints,
                                              R1TensorT<2>& min,
                                              R1TensorT<2>& max)
{
  localPoints.clear();
  localPoints.reserve(points.size());
  min = std::numeric_limits<realT>::max();
  max = -std::numeric_limits<realT>::max();
  for(lArray1d::size_type a = 0 ; a < points.size() ; ++a)
  {
    R1TensorT<2> tmp;
    tmp(0) = Dot(points[a], e1);
    tmp(1) = Dot(points[a], e2);
    max.SetMax(tmp);
    min.SetMin(tmp);
    localPoints.push_back(tmp);
  }
}

/**
 * @brief Returns the signed distance to the plane
 * @author Scott Johnson
 * @param[in] normal plane normal
 * @param[in] pointOnPlane any point on the plane
 * @param[in] pointToQuery point in space to query for distance to the plane
 * @return Signed distance to the plane (along normal)
 */
realT DistanceToPlane(const R1Tensor& normal,
                      const R1Tensor& pointOnPlane,
                      const R1Tensor& pointToQuery)
{
  R1Tensor tmp(pointToQuery);
  tmp -= pointOnPlane;
  return Dot(tmp, normal);
}

/**
 * @brief Clip a 3D convex polygon by a plane
 * @author Scott Johnson
 * Provides both the inside and outside surfaces as a tesselation
 * @param[in] normal plane normal
 * @param[in] pointOnPlane a point that lies on the plane and is inside the
 * projection of the polygon on the plane
 * @param[in] points ordered list of points that describes the boundary of the
 * polygon
 * @param[in] e1 a unit vector in the plane
 * @param[in] e2 a unit vector in the plane that is orthogonal to e1
 * @param[in] penetration list of signed distances to the common plane
 * @param[out] inContact list of triangles for each face overlapping the plane
 *(0 points ... or 2 points with the pointOnPlane implied as the terminus)
 * @param[out] outContact list of triangles for each face not overlapping the
 * plane (0 points ... or 2 points with the pointOnPlane implied as the
 * terminus)
 * @param[in] reverse Option to specify that both the polygon loop order should
 * be reversed and the normal inverted
 */
void ClipByPlane_3DPolygon(const R1Tensor& normal,
                           const R1Tensor& pointOnPlane,
                           const array<R1Tensor>& points,
                           const R1Tensor& e1,
                           const R1Tensor& e2,
                           const array<real64>& penetration,
                           array<array<R1TensorT<2> > >& inContact,
                           array<array<R1TensorT<2> > >& outContact,
                           const bool reverse)
{
  if(inContact.size() != points.size())
  {
    inContact.clear();
    inContact.resize(points.size());
  }
  if(outContact.size() != points.size())
  {
    outContact.clear();
    outContact.resize(points.size());
  }

  const int num_nod = static_cast<int>(points.size());
  const realT dfct = reverse ? 1.0 : -1.0;
  for ( int iii = 0 ; iii < num_nod ; ++iii)
  {
    const int n = reverse ? (num_nod - iii - 1) : iii;
    const int nplus1 = reverse ? (iii == (num_nod-1) ? (num_nod-1) : n-1) : (n == (num_nod-1) ? 0 : n + 1);

    //find intersection of straight segment [n,n+1] with the common plane
    //only search for an intersection point if the nodal positions
    //lie on opposite sides of the common plane
    if (penetration[n] * penetration[nplus1] < 0.0)
    {
      //get the fraction of the distance along the line segment at which the
      // line segment intersects the common plane
      const realT prod = fabs(penetration[n]) / (fabs(penetration[n]) + fabs(penetration[nplus1]));

      //set dxr to the point of intersection
      R1Tensor dxr(points[nplus1]);
      dxr -= points[n];
      dxr *= prod;
      dxr += points[n];
      dxr -= pointOnPlane;

      R1TensorT<2> a;
      a(0) = Dot(e1, dxr);
      a(1) = Dot(e2, dxr);
      inContact[n].push_back(a);
      outContact[n].push_back(a);
    }

    //node nplus1 penetrated
    {
      R1Tensor dxr(normal);
      dxr *= dfct * penetration[nplus1];
      dxr += points[nplus1];
      dxr -= pointOnPlane;

      // set the penetrating nodes' e1-e2 coordinate vector position
      // set the local nodal index
      // increment the array position
      R1TensorT<2> a;
      a(0) = Dot(e1, dxr);
      a(1) = Dot(e2, dxr);
      if (penetration[nplus1] > 0)
      {
        inContact[n].push_back(a);
        inContact[nplus1].insert(inContact[nplus1].begin(),a);
      }
      else
      {
        outContact[n].push_back(a);
        outContact[nplus1].insert(outContact[nplus1].begin(),a);
      }
    }
  }
}

/**
 * @brief Clip a 3D convex polygon by a plane: only provides an approximation of
 * the in and out of contact areas
 * @author Scott Johnson
 * Provides both the inside and outside surfaces as a tesselation
 * @param[in] normal plane normal
 * @param[in] pointOnPlane a point that lies on the plane and is inside the
 * projection of the polygon on the plane
 * @param[in] points ordered list of points that describes the boundary of the
 * polygon
 * @param[in] e1 a unit vector in the plane
 * @param[in] e2 a unit vector in the plane that is orthogonal to e1
 * @param[in] penetration list of signed normal distances
 * @param[out] inContact list of points along the edge that are in contact or
 * intersection with the common plane
 * @param[out] outContact list of points along the edge that are out of contact
 * or are intersecting with the common plane
 * @param[out] inDistance list of signed distances to the common plane for each
 * of the in contact points
 * @param[out] outDistance list of signed distances to the common plane for each
 * of the out of contact points
 * @param[in] reverse Option to specify that both the polygon loop order should
 * be reversed and the normal inverted
 */
void ClipByPlane_3DPolygon(const R1Tensor& normal,
                           const R1Tensor& pointOnPlane,
                           const array<R1Tensor>& points,
                           const R1Tensor& e1,
                           const R1Tensor& e2,
                           const array<real64>& penetration,
                           array<R1TensorT<2> >& inContact,
                           array<R1TensorT<2> >& outContact,
                           array<real64>& inDistance,
                           array<real64>& outDistance,
                           const bool reverse)
{
  inContact.clear();
  outContact.clear();
  inDistance.clear();
  outDistance.clear();

  const int num_nod = static_cast<int>(points.size());
  const realT dfct = reverse ? 1.0 : -1.0;
  for ( int iii = 0 ; iii < num_nod ; ++iii)
  {
    const int n = reverse ? (num_nod - iii - 1) : iii;
    const int nplus1 = reverse ? (iii == (num_nod-1) ? (num_nod-1) : n-1) : (n == (num_nod-1) ? 0 : n + 1);

    //node n penetrated
    {
      R1Tensor dxr(normal);
      dxr *= dfct * penetration[n];
      dxr += points[n];
      dxr -= pointOnPlane;

      // set the penetrating nodes' e1-e2 coordinate vector position
      // set the local nodal index
      // increment the array position
      R1TensorT<2> a;
      a(0) = Dot(e1, dxr);
      a(1) = Dot(e2, dxr);
      if (penetration[n] > 0)
      {
        inContact.push_back(a);
        inDistance.push_back(penetration[n]);
      }
      else
      {
        outContact.push_back(a);
        outDistance.push_back(penetration[n]);
      }
    }

    //find intersection of straight segment [n,n+1] with the common plane
    //only search for an intersection point if the nodal positions
    //lie on opposite sides of the common plane
    if (penetration[n] * penetration[nplus1] < 0.0)
    {
      //get the fraction of the distance along the line segment at which the
      // line segment intersects the common plane
      const realT prod = fabs(penetration[n]) / (fabs(penetration[n]) + fabs(penetration[nplus1]));

      //set dxr to the point of intersection
      R1Tensor dxr(points[nplus1]);
      dxr -= points[n];
      dxr *= prod;
      dxr += points[n];
      dxr -= pointOnPlane;

      R1TensorT<2> a;
      a(0) = Dot(e1, dxr);
      a(1) = Dot(e2, dxr);
      inContact.push_back(a);
      inDistance.push_back(0.0);
      outContact.push_back(a);
      outDistance.push_back(0.0);
    }
  }
}

/**
 * @brief Clip a 3D convex polygon by a plane: only provides an approximation of
 * the in and out of contact areas
 * @author Scott Johnson
 * Provides both the inside and outside surfaces as a tesselation
 * @param[in] normal plane normal
 * @param[in] pointOnPlane a point that lies on the plane and is inside the
 * projection of the polygon on the plane
 * @param[in] points ordered list of points that describes the boundary of the
 * polygon
 * @param[out] inContact list of points along the edge that are in contact or
 * intersection with the common plane
 * @param[out] outContact list of points along the edge that are out of contact
 * or are intersecting with the common plane
 * @param[in] reverse Option to specify that both the polygon loop order should
 * be reversed and the normal inverted
 */
void ClipByPlane_3DPolygon(const R1Tensor& normal, const R1Tensor& pointOnPlane,
                           const array<R1Tensor>& points,
                           array<array<R1TensorT<2> > >& inContact,
                           array<array<R1TensorT<2> > >& outContact, const bool reverse)
{
  //get e1, e2
  R1Tensor e1, e2;
  VectorsInPlane(normal, e1, e2);

  //get normal penetration distances
  array<real64> penetration(points.size(), 0.0);
  array<real64>::size_type a = 0;
  for (array<R1Tensor>::const_iterator it = points.begin() ; it != points.end() ; ++it, ++a)
    penetration[a] = (reverse ? -1.0 : 1.0) * DistanceToPlane(normal, pointOnPlane, *it);

  //call the base function
  ClipByPlane_3DPolygon(normal, pointOnPlane, points, e1, e2, penetration, inContact, outContact,
                        reverse);
}



///Get the strike vector for the given up and normal vectors
/**
   @param[in] nt Unit normal to the fault surface (not checked to be unit)
   @param[in] upt Unit vector in the direction of vertical
   @param[out] strikeVec Unit vector along the strike direction
 */
void Strike( const R1Tensor& nt,
             const R1Tensor& upt,
             R1Tensor& strikeVec)
{
  const realT dot = Dot(upt, nt);
  const realT dot2 = dot * dot;
  if(isEqual(dot2, 1.0))
    throw GPException("The plane normal is aligned with up, so plane frame transformation is ill-defined.");

  strikeVec = nt;
  R1Tensor tmp(upt);
  tmp *= -dot;
  strikeVec += tmp;
  tmp = strikeVec;
  strikeVec.Cross(tmp, upt);
  strikeVec.Normalize();
}

///Get the strike scalar for the given unit strike vector (both in the plane of
// the ground)
/**
   @param[in] north Unit vector in the direction of north
   @param[in] strike Unit vector in the strike direction
   @return Strike scalar
 */
realT StrikeAngle( const R1Tensor& nt,
                   const R1Tensor& upt,
                   const R1Tensor& strikeVec)
{
  const realT dot = Dot(nt, strikeVec);
  if(isEqual(dot,-1.0))
    return 0.5 * M_PI;
  else if(isEqual(dot, 1.0))
    return 0;

  R1Tensor east;
  east.Cross(nt, upt);
  const realT edot = Dot(east, strikeVec);
  return edot < 0 ? 2*M_PI - acos(dot) : acos(dot);
}

void Strike( const R1Tensor& east,
             const R1Tensor& north,
             const realT strikeAngle,
             R1Tensor& strike)
{
  const realT efct = sin(strikeAngle);
  const realT nfct = cos(strikeAngle);
  R1Tensor tmp(east);
  tmp *= efct;
  strike = tmp;
  tmp = north;
  tmp *= nfct;
  strike += tmp;
}

void Dip( const R1Tensor& up,
          const R1Tensor& strike,
          const realT dipAngle,
          R1Tensor& dip)
{
  const realT ufct = sin(dipAngle);
  const realT nfct = cos(dipAngle);
  R1Tensor tmp;
  tmp.Cross(strike, up);
  tmp *= nfct;
  dip = tmp;
  tmp = up;
  tmp *= ufct;
  dip += tmp;
}


///Get the transform matrix from plane frame to space frame
/**
   @param[in] normalToPlane Vector normal to the fault surface
   @param[in] up Vector in the direction of vertical
   @param[out] T Transform tensor
 */
void TransformPlaneFrame(const R1Tensor& normalToPlane,
                         const R1Tensor& up,
                         R2Tensor& T)
{
  T = 0.0;

  R1Tensor upt(up);
  upt.Normalize();

  R1Tensor nt(normalToPlane);
  nt.Normalize();

  //GET STRIKE
  R1Tensor strikeVec;
  Strike(nt, upt, strikeVec);

  //GET DIP
  R1Tensor dipVec;
  dipVec.Cross(strikeVec, nt);

  T.FillColumn(0, strikeVec);
  T.FillColumn(1, nt);
  T.FillColumn(2, dipVec);
}

///Get the transform matrix from strike frame to space frame
/**
   @param[in] normalToPlane Vector normal to the fault surface
   @param[in] up Vector in the direction of vertical
   @param[out] T Transform tensor
 */
void TransformStrikeFrame(const R1Tensor& normalToPlane,
                          const R1Tensor& up,
                          R2Tensor& T)
{

  R1Tensor upt(up);
  upt.Normalize();

  R1Tensor nt(normalToPlane);
  nt.Normalize();

  //GET STRIKE
  R1Tensor strikeVec;
  Strike(nt, upt, strikeVec);

  //GET NORMAL PROJECTED TO THE GROUND
  R1Tensor normVec(nt);
  {
    R1Tensor tmp(upt);
    tmp *= Dot(nt, upt);
    normVec -= tmp;
    normVec.Normalize();
  }

  T.FillColumn(0, strikeVec);
  T.FillColumn(1, normVec);
  T.FillColumn(2, upt);
}

///Get the transform matrix from strike frame to space frame
/**
   @param[in] strike Vector along the strike
   @param[in] up Vector in the direction of vertical
   @param[out] T Transform tensor
 */
void TransformStrikeFrame2(const R1Tensor& strike,
                           const R1Tensor& up,
                           R2Tensor& T)
{
  R1Tensor upt(up);
  upt.Normalize();

  R1Tensor st(strike);
  st.Normalize();

  //GET normal
  R1Tensor nt;
  nt.Cross(st, upt);

  T.FillColumn(0, st);
  T.FillColumn(1, nt);
  T.FillColumn(2, upt);
}

///Get the transform matrix from WPP frame to space frame
/**
   @param[in] north Vector in the direction of north
   @param[in] up Vector in the direction of vertical
   @param[out] T Transform tensor
 */
void TransformWPPFrame(const R1Tensor& north,
                       const R1Tensor& up,
                       R2Tensor& T)
{

  R1Tensor upt(up);
  upt.Normalize();

  R1Tensor northt(north);
  northt.Normalize();

  R1Tensor east;
  east.Cross(northt, up);
  upt *= -1.0;

  T.FillColumn(0, northt);
  T.FillColumn(1, east);
  T.FillColumn(2, upt);
}

///Transform the R2 tensor from a frame to another
/**
   @param[in] transform Transform matrix TO prime frame
   @param[in] tensor Current matrix
   @param[out] tensorPrime Matrix transformed to the prime frame
 */
void TransformTensorFrame(const R2Tensor& transform,
                          const R2Tensor& tensor,
                          R2Tensor& tensorPrime)
{
  R2Tensor tmp;
  tmp.AijBjk(transform, tensor);
  tensorPrime.AijBkj(tmp, transform);
}

///Transform the R2 tensor from a frame to another
/**
   @param[in] transform Transform matrix FROM prime frame
   @param[in] tensor Current matrix
   @param[out] tensorPrime Matrix transformed to the prime frame
 */
void TransformTensorFrameTranspose(const R2Tensor& transform,
                                   const R2Tensor& tensor,
                                   R2Tensor& tensorPrime)
{
  //(transform * tensor)
  R2Tensor tmp;
  tmp.AjiBjk(transform, tensor);
  //((transform * tensor) * transform^T)
  tensorPrime.AijBjk(tmp, transform);
}

///Calculate rotation matrix for given strike/dip
/**
   @param[in] strike angle (degrees)
   @param[in] dip angle (degrees)
   @param[out] rotation matrix
 */
void StrikeDipRotationMatrix(const realT strike,
                             const realT dip,
                             bool invert,
                             R2Tensor& rotationMatrix)
{
  realT stheta = sin(strike);
  realT ctheta = cos(strike);
  realT sdelta = sin(dip);
  realT cdelta = cos(dip);
  if (invert)
  {
    // x-north, y-east, z-down
    rotationMatrix(0,0) = ctheta;
    rotationMatrix(0,1) = stheta;
    rotationMatrix(0,2) = 0.0;
    rotationMatrix(1,0) = -stheta*cdelta;
    rotationMatrix(1,1) = ctheta*cdelta;
    rotationMatrix(1,2) = -sdelta;
    rotationMatrix(2,0) = -stheta*sdelta;
    rotationMatrix(2,1) = ctheta*sdelta;
    rotationMatrix(2,2) = cdelta;

    // x-north, y-down, z-east
    // rotationMatrix(0,0) = cdelta*ctheta;
    // rotationMatrix(0,1) = -sdelta;
    // rotationMatrix(0,2) = -cdelta*stheta;
    // rotationMatrix(1,0) = sdelta*ctheta;
    // rotationMatrix(1,1) = cdelta;
    // rotationMatrix(1,2) = -sdelta*stheta;
    // rotationMatrix(2,0) = stheta;
    // rotationMatrix(2,1) = 0.0;
    // rotationMatrix(2,2) = ctheta;
  }
  else
  {
    rotationMatrix(0,0) = ctheta;
    rotationMatrix(1,0) = stheta;
    rotationMatrix(2,0) = 0.0;
    rotationMatrix(0,1) = -stheta*cdelta;
    rotationMatrix(1,1) = ctheta*cdelta;
    rotationMatrix(2,1) = -sdelta;
    rotationMatrix(0,2) = -stheta*sdelta;
    rotationMatrix(1,2) = ctheta*sdelta;
    rotationMatrix(2,2) = cdelta;

    // rotationMatrix(0,0) = cdelta*ctheta;
    // rotationMatrix(1,0) = -sdelta;
    // rotationMatrix(2,0) = -cdelta*stheta;
    // rotationMatrix(0,1) = sdelta*ctheta;
    // rotationMatrix(1,1) = cdelta;
    // rotationMatrix(2,1) = -sdelta*stheta;
    // rotationMatrix(0,2) = stheta;
    // rotationMatrix(1,2) = 0.0;
    // rotationMatrix(2,2) = ctheta;
  }
}


} /* namespace GeometryUtilities */
