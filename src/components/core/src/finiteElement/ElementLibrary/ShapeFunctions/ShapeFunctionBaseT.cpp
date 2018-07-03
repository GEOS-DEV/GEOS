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
 * File: ShapeFunctionBaseT.h
 * Shape function Class
 *
 * created : RRS (09/14/2010)
 */

#include "ShapeFunctionBaseT.h"

ShapeFunctionBaseT::ShapeT() :
  rElement(0),
  pNodes(0)
{}

ShapeFunctionBaseT::~ShapeFunctionBaseT(void)
{}



void ShapeFunctionBaseT::CalculateJacobian(const int elem)
{
  J = 0.0;
  for( int i=1 ; i<=nsdof ; ++i )
    for( int j=1 ; j<=nsdof ; ++j )
      for( int a=1 ; a<=rElement->NumNodeElem() ; ++a )
        J(i,j) += pNodes->Xref( rElement->Connectivity(elem,a),i)
                  * dNdXi(a)(j);

}

void ContinuumShapeT::Calc_Shape_Deriv(const realT fac)
{

  realT factor = 1.0 / pow(static_cast<realT>(2),nsdof);

  // loop over "nodal indexed" shape functions and calculate the
  // shape function derivitives wrt the local coordinates
  for( int a=1 ; a<=rElement->NumNodeElem() ; ++a )
  {
    for( int i=1 ; i<=nsdof ; ++i )
    {
      dNdXi(a)(i) = factor * rElement->Xi_node(a,i);
      for( int j=1 ; j<=nsdof ; j++ )
        if( i!=j )
          dNdXi(a)(i) *= ( 1 + ip_coord_fac * rElement->Xi_node(ip,j) * rElement->Xi_node(a,j));

    }
    //std::cout<<a<<' '<<dNdXi(a)<<std::endl;
  }

  CalculateJacobian(elem);

  m_detJ(elem,ip) = J.Det();

  dNdX.AijBi(Jinv,dNdXi);


}


inline realT ShapeFunctionBaseT::Shape( const R1Tensor& Xi,
                                        const R1Tensor& Xi_node )
{
  realT N = 1.0 / pow(2.0,nsdof);

  for( int i=1 ; i<=nsdof ; ++i )
    N *= ( 1 +  Xi(i) * Xi_node(i) );

  return N;
}
