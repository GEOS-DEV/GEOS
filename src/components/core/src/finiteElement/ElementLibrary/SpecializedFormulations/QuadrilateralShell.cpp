// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * @file QuadrilateralShell.cpp
 * @author Fu, Pengcheng
 * @date July 12, 2012
 */

#include "QuadrilateralShell.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


QuadrilateralShell::QuadrilateralShell():
  FiniteElement<3>(1,4,0)
{
  m_nodeOrdering.resize(4);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;
  m_nodeOrdering[3] = 3;

}

QuadrilateralShell::~QuadrilateralShell()
{
  // TODO Auto-generated destructor stub
}

//Based on the formulation of UniformStrainQuadrilateral
//First has to map the quadrilateral from the 3D space to a 2D plane
void QuadrilateralShell::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{
  assert(mapped_support_points.size() == n_dofs);

  const unsigned int q = 0;



  R1Tensor x[4] = { mapped_support_points[m_nodeOrdering[0]],
                    mapped_support_points[m_nodeOrdering[1]],
                    mapped_support_points[m_nodeOrdering[2]],
                    mapped_support_points[m_nodeOrdering[3]] };

  //Transform

  R1Tensor x12 = x[1];
  x12 -= x[0];
  realT l12 = x12.Normalize();

  R1Tensor x13 = x[2];
  x13 -= x[0];
  realT l13 = x13.Normalize();
  realT cos213 = Dot(x13, x12);

  R1Tensor x14 = x[3];
  x14 -= x[0];
  realT l14 = x14.Normalize();
  realT cos214 = Dot(x14, x12);

  x[0] *= 0.0;

  x[1][0] = l12;
  x[1][1] = 0.0;
  x[1][2] = 0.0;

  x[2][0] = l13 * cos213;
  x[2][1] = l13 * sqrt(1.0 - pow(cos213,2));
  x[2][2] = 0.0;

  x[3][0] = l14 * cos214;
  x[3][1] = l14 * sqrt(1.0 - pow(cos214,2));
  x[3][2] = 0.0;

  // Now x is the coordinates in the 2D plane.

  R1Tensor x4_x2 = x[3];
  x4_x2 -= x[1];
  x4_x2 *= 0.5;

  R1Tensor x3_x1 = x[2];
  x3_x1 -= x[0];
  x3_x1 *= 0.5;

  R1Tensor b[4];

  b[0][0] = -x4_x2[1];
  b[0][1] =  x4_x2[0];
  b[0][2] =  1.0;

  b[1][0] =  x3_x1[1];
  b[1][1] = -x3_x1[0];
  b[1][2] =  1.0;

  b[2][0] =  x4_x2[1];
  b[2][1] = -x4_x2[0];
  b[2][2] =  1.0;

  b[3][0] = -x3_x1[1];
  b[3][1] =  x3_x1[0];
  b[3][2] =  1.0;



  data[q].jacobian_determinant = 2.0 * (x3_x1[0]*x4_x2[1] - x3_x1[1]*x4_x2[0]);

//  std::cout<<"data[q].jacobian_determinant =
// "<<data[q].jacobian_determinant<<std::endl;


  for( int a=0 ; a<4 ; ++a )
  {
    data[q].mapped_gradients[a](0) = b[m_nodeOrdering[a]][0] / data[q].jacobian_determinant;
    data[q].mapped_gradients[a](1) = b[m_nodeOrdering[a]][1] / data[q].jacobian_determinant;
    data[q].mapped_gradients[a](2) = 0.0;

//    std::cout<<"data["<<q<<"].mapped_gradients["<<a<<"] =
// "<<data[q].mapped_gradients[a]<<std::endl;
  }

}
