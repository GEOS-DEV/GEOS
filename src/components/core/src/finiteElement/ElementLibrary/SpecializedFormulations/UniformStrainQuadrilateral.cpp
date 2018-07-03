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

/**
 * @file UniformStrainQuadrilateral.cpp
 * @author settgast1
 * @date Jun 6, 2011
 */

#include "UniformStrainQuadrilateral.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


static void
CalcFBHourForce( const R1Tensor vel[4],
                 const realT gamma[4],
                 const realT& dampcoef,
                 const realT& stiffcoef,
                 const realT& rho,
                 const realT& modulus,
                 const realT& area,
                 const array<R1Tensor>& dNdx,
                 const realT& dt,
                 array<R1Tensor>& Qstiffness,
                 R1Tensor hgforce[4] );


UniformStrainQuadrilateral::UniformStrainQuadrilateral():
  FiniteElement<3>(1,4,1)
{
  m_nodeOrdering.resize(4);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 3;
  m_nodeOrdering[3] = 2;


}

UniformStrainQuadrilateral::~UniformStrainQuadrilateral()
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

void UniformStrainQuadrilateral::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{
  assert(mapped_support_points.size() == n_dofs);

  const unsigned int q = 0;



  const R1Tensor x[4] = { mapped_support_points[m_nodeOrdering[0]],
                          mapped_support_points[m_nodeOrdering[1]],
                          mapped_support_points[m_nodeOrdering[2]],
                          mapped_support_points[m_nodeOrdering[3]] };


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

void UniformStrainQuadrilateral::zero_energy_mode_control( const array<R1Tensor>& dNdx,
                                                           const realT& volume,
                                                           const array<R1Tensor>& x0,
                                                           const array<R1Tensor>& v0,
                                                           const realT& dampcoef,
                                                           const realT& stiffcoef,
                                                           const realT& rho,
                                                           const realT& modulus,
                                                           const realT& dt,
                                                           array<R1Tensor>& Qstiffness,
                                                           array<R1Tensor>& force )
{

  const R1Tensor x[4] = { x0[m_nodeOrdering[0]],
                          x0[m_nodeOrdering[1]],
                          x0[m_nodeOrdering[2]],
                          x0[m_nodeOrdering[3]] };

  const R1Tensor vel[4] = { v0[m_nodeOrdering[0]],
                            v0[m_nodeOrdering[1]],
                            v0[m_nodeOrdering[2]],
                            v0[m_nodeOrdering[3]] };

  R1Tensor hgforce[4];

  const realT inv4A = 1.0 / ( 4.0 * volume );

/*
   const realT gamma[4] = { inv4A * ( x2*(y3-y4) + x3*(y4-y2) + x4*(y2-y3) ),
                           inv4A * ( x3*(y1-y4) + x4*(y3-y1) + x1*(y4-y3) ),
                           inv4A * ( x4*(y1-y2) + x1*(y2-y4) + x2*(y4-y1) ),
                           inv4A * ( x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1) ) };
 */

  const realT gamma[4] = { inv4A * ( x[1][0]*(x[2][1]-x[3][1]) + x[2][0]*(x[3][1]-x[1][1]) + x[3][0]*(x[1][1]-x[2][1]) ),
                           inv4A * ( x[2][0]*(x[0][1]-x[3][1]) + x[3][0]*(x[2][1]-x[0][1]) + x[0][0]*(x[3][1]-x[2][1]) ),
                           inv4A * ( x[3][0]*(x[0][1]-x[1][1]) + x[0][0]*(x[1][1]-x[3][1]) + x[1][0]*(x[3][1]-x[0][1]) ),
                           inv4A * ( x[0][0]*(x[2][1]-x[1][1]) + x[1][0]*(x[0][1]-x[2][1]) + x[2][0]*(x[1][1]-x[0][1]) ) };

  CalcFBHourForce( vel, gamma, dampcoef, stiffcoef, rho, modulus, volume, dNdx, dt, Qstiffness, hgforce );

  for( int a=0 ; a<4 ; ++a )
  {
    force[a] = hgforce[m_nodeOrdering[a]];
  }

}


void
CalcFBHourForce( const R1Tensor vel[4],
                 const realT gamma[4],
                 const realT& dampcoef,
                 const realT& stiffcoef,
                 const realT& rho,
                 const realT& modulus,
                 const realT& area,
                 const array<R1Tensor>& dNdx,
                 const realT& dt,
                 array<R1Tensor>& Qstiffness,
                 R1Tensor hgforce[4] )
{


  R1Tensor q, Q;
  R1Tensor temp;

  const realT BB = Dot(dNdx[0],dNdx[0]) + Dot(dNdx[1],dNdx[1]) + Dot(dNdx[2],dNdx[2]) + Dot(dNdx[3],dNdx[3]);

  const realT Cdamp  = dampcoef * sqrt( rho*modulus*BB / 6.0 ) * area;

  const realT Cstiff = stiffcoef * modulus * BB * area / 3.0;


  q = 0.0;
  for( int a=0 ; a<4 ; ++a )
  {
    temp = vel[a];
    temp *= gamma[a];
    q += temp;
  }

  q *= 0.5;

  temp = q;
  temp *= Cstiff * dt;
  Qstiffness[0] += temp;

  q *= Cdamp;
  q += Qstiffness[0];
  q *= 0.5;

  for( int a=0 ; a<4 ; ++a )
  {
    temp  = q;
    temp *= -gamma[a];
    hgforce[a] = temp;
  }
}
