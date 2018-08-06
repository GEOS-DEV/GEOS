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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TriangleShell.cpp
 * @author Fu, Pengcheng
 * @date August 3, 2012
 */

#include "TriangleShell.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


TriangleShell::TriangleShell():
  FiniteElement<3>(1,3,0)
{
  m_nodeOrdering.resize(3);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;

}

TriangleShell::~TriangleShell()
{
  // TODO Auto-generated destructor stub
}

void TriangleShell::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{

  assert(mapped_support_points.size() == n_dofs);
  //  printf("n=%d\n", n_dofs);

  //See Chapter 15 of Int. FEM of U Colorado by Carlos Felippa for detailed
  // formulation.
  //http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/
  //Accessed in July 2012

  const std::vector<R1TensorT<3> >& Y = mapped_support_points;

  // transform

  realT l1 = sqrt(pow(Y[1][0] - Y[0][0], 2.0) + pow(Y[1][1] - Y[0][1], 2.0) + pow(Y[1][2] - Y[0][2], 2.0));
  realT l2 = sqrt(pow(Y[2][0] - Y[0][0], 2.0) + pow(Y[2][1] - Y[0][1], 2.0) + pow(Y[2][2] - Y[0][2], 2.0));

  realT theta = ((Y[1][0] - Y[0][0]) * (Y[2][0] - Y[0][0]) + (Y[1][1] - Y[0][1]) * (Y[2][1] - Y[0][1]) + (Y[1][2] - Y[0][2]) * (Y[2][2] - Y[0][2])) / l1 / l2;


  std::vector<R1TensorT<3> > X(n_dofs);

  X[0][0] = 0.0;
  X[0][1] = 0.0;

  X[1][0] = l1;
  X[1][1] = 0.0;

  X[2][0] = l2 * theta;
  X[2][1] = l2 * sqrt(1.0 - theta * theta);

  realT V;
  const realT half = 1.0 / 2.0;

  V = (X[1][0] * X[2][1] - X[2][0] * X[1][1]) + (X[2][0] * X[0][1] - X[0][0] * X[2][1]) + (X[0][0] * X[1][1] - X[1][0] * X[0][1]);
  V *= half;
  data[0].jacobian_determinant = V;

  data[0].mapped_gradients[0](0) = (X[1][1] - X[2][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[0](1) = (X[2][0] - X[1][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[0](2) = 0.0;

  data[0].mapped_gradients[1](0) = (X[2][1] - X[0][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](1) = (X[0][0] - X[2][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](2) = 0.0;


  data[0].mapped_gradients[2](0) = (X[0][1] - X[1][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](1) = (X[1][0] - X[0][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](2) = 0.0;



}
