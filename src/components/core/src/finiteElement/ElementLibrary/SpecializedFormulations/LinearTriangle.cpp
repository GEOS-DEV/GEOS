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
 * @file LinearTriangle.cpp
 * @author Fu, Pengcheng
 * @date July 3, 2012
 */

#include "LinearTriangle.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


LinearTriangle::LinearTriangle():
  FiniteElement<2>(1,3,0)
{
  m_nodeOrdering.resize(3);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;

}

LinearTriangle::~LinearTriangle()
{
  // TODO Auto-generated destructor stub
}


/**
 * Reinitialize the finite element basis on a particular element.
 * We use the coordinates of the support points in real space to
 * construct the forward mapping from the parent coordinate system.  The
 * support points are assumed to follow a lexicographic ordering:
 * On the parent element, we loop over the x-coordinate fastest,
 * the y, then z (depending on the desired spatial dimension of the
 * element).
 */

void LinearTriangle::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{

  assert(mapped_support_points.size() == n_dofs);

  //See Chapter 15 of Int. FEM of U Colorado by Carlos Felippa for detailed
  // formulation.
  //http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/
  //Accessed in July 2012

  const std::vector<R1TensorT<3> >& X = mapped_support_points;

  realT V;
  const realT half = 1.0 / 2.0;



  V = (X[1][0] * X[2][1] - X[2][0] * X[1][1]) + (X[2][0] * X[0][1] - X[0][0] * X[2][1]) + (X[0][0] * X[1][1] - X[1][0] * X[0][1]);
  V *= half;
  data[0].jacobian_determinant = V;

  data[0].mapped_gradients[0](0) = (X[1][1] - X[2][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[0](1) = (X[2][0] - X[1][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](0) = (X[2][1] - X[0][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](1) = (X[0][0] - X[2][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](0) = (X[0][1] - X[1][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](1) = (X[1][0] - X[0][0]) * half / data[0].jacobian_determinant;



}
