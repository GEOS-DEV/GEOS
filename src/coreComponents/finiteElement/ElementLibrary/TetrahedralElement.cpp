/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file Tetrahedron.cpp
 */

#include "TetrahedralElement.hpp"

#include "finiteElement/basis/LagrangeBasis.hpp"
#include "finiteElement/quadrature/GaussQuadrature.hpp"

namespace geosx
{

TetrahedralElement::TetrahedralElement( BasisBase const & basis,
                                        QuadratureBase const & quadrature,
                                        const int num_zero_energy_modes ):
  FiniteElement<3>( basis, quadrature, num_zero_energy_modes )
{
  m_nodeOrdering.resize(4);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;
  m_nodeOrdering[3] = 3;

}

TetrahedralElement::~TetrahedralElement()
{}



void TetrahedralElement::reinit( array1d<R1TensorT<3> > const & mapped_support_points )
{

  assert(mapped_support_points.size() == n_dofs);

  const unsigned int q = 0;

  array1d<R1TensorT<3> > const & X = mapped_support_points;

  realT V, a[4], b[4], c[4];
  const realT sixth = 1.0 / 6.0;


  //a[1]=(y4-y2)*(z3-z2)-(y3-y2)*(z4-z2)
  a[0] = (X[3][1] - X[1][1]) * (X[2][2] - X[1][2]) - (X[2][1] - X[1][1]) * (X[3][2] - X[1][2]);
  //      b[1]=x32 z42 - x42 z32
  b[0] = (X[2][0] - X[1][0]) * (X[3][2] - X[1][2]) - (X[3][0] - X[1][0]) * (X[2][2] - X[1][2]);
  //      c[1]=x42 y32 - x32 y42
  c[0] = (X[3][0] - X[1][0]) * (X[2][1] - X[1][1]) - (X[2][0] - X[1][0]) * (X[3][1] - X[1][1]);
  //      a[2]=y31 z43 - y34 z13
  a[1] = (X[2][1] - X[0][1]) * (X[3][2] - X[2][2]) - (X[2][1] - X[3][1]) * (X[0][2] - X[2][2]);
  //      b[2]=x43 z31 - x13 z34
  b[1] = (X[3][0] - X[2][0]) * (X[2][2] - X[0][2]) - (X[0][0] - X[2][0]) * (X[2][2] - X[3][2]);
  //      c[2]=x31 y43 - x34 y13
  c[1] = (X[2][0] - X[0][0]) * (X[3][1] - X[2][1]) - (X[2][0] - X[3][0]) * (X[0][1] - X[2][1]);
  //      a[3]=y24 z14 - y14 z24
  a[2] = (X[1][1] - X[3][1]) * (X[0][2] - X[3][2]) - (X[0][1] - X[3][1]) * (X[1][2] - X[3][2]);
  //      b[3]=x14 z24 - x24 z14
  b[2] = (X[0][0] - X[3][0]) * (X[1][2] - X[3][2]) - (X[1][0] - X[3][0]) * (X[0][2] - X[3][2]);
  //      c[3]=x24 y14 - x14 y24
  c[2] = (X[1][0] - X[3][0]) * (X[0][1] - X[3][1]) - (X[0][0] - X[3][0]) * (X[1][1] - X[3][1]);
  //      a[4]=y13 z21 - y12 z31
  a[3] = (X[0][1] - X[2][1]) * (X[1][2] - X[0][2]) - (X[0][1] - X[1][1]) * (X[2][2] - X[0][2]);
  //      b[4]=x21 z13 - x31 z12
  b[3] = (X[1][0] - X[0][0]) * (X[0][2] - X[2][2]) - (X[2][0] - X[0][0]) * (X[0][2] - X[1][2]);
  //      c[4]=x13 y21 - x12 y31
  c[3] = (X[0][0] - X[2][0]) * (X[1][1] - X[0][1]) - (X[0][0] - X[1][0]) * (X[2][1] - X[0][1]);
  //
  //6V=x21 (y23 z34 - y34 z23) + x32 (y34 z12 - y12 z34) + x43 (y12 z23 - y23
  // z12),
  V = (X[1][0] - X[0][0]) * ((X[1][1] - X[2][1]) * (X[2][2] - X[3][2]) - (X[2][1] - X[3][1]) * (X[1][2] - X[2][2])) + (X[2][0] - X[1][0]) *
      ((X[2][1] - X[3][1]) * (X[0][2] - X[1][2]) - (X[0][1] - X[1][1]) * (X[2][2] - X[3][2])) + (X[3][0] - X[2][0]) *
      ((X[0][1] - X[1][1]) * (X[1][2] - X[2][2]) - (X[1][1] - X[2][1]) * (X[0][2] - X[1][2]));
  V *= sixth;

  data[q].jacobian_determinant = V;
  for( int iNd=0 ; iNd<4 ; ++iNd )
  {
    data[q].mapped_gradients[iNd](0) = a[iNd] * sixth / data[q].jacobian_determinant;
    data[q].mapped_gradients[iNd](1) = b[iNd] * sixth / data[q].jacobian_determinant;
    data[q].mapped_gradients[iNd](2) = c[iNd] * sixth / data[q].jacobian_determinant;
  }


}

}
