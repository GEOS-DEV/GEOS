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
 * @file IntegrationPointsT.cpp
 * @author settgast1
 * @date Dec 6, 2010
 */

#include "IntegrationRuleT.h"

IntegrationRuleT::IntegrationRuleT()
{
  // TODO Auto-generated constructor stub

}

IntegrationRuleT::~IntegrationRuleT()
{
  // TODO Auto-generated destructor stub
}


void IntegrationRuleT::CalculateShapeFunctionDerivatives( const array<R1Tensor>& X,
                                                          array<R1Tensor>& dNdX,
                                                          realT& detJ )
{
  const realT x[8] = { X(0)(0),
                       X(1)(0),
                       X(2)(0),
                       X(3)(0),
                       X(4)(0),
                       X(5)(0),
                       X(6)(0),
                       X(7)(0) };

  const realT y[8] = { X(0)(1),
                       X(1)(1),
                       X(2)(1),
                       X(3)(1),
                       X(4)(1),
                       X(5)(1),
                       X(6)(1),
                       X(7)(1) };

  const realT z[8] = { X(0)(2),
                       X(1)(2),
                       X(2)(2),
                       X(3)(2),
                       X(4)(2),
                       X(5)(2),
                       X(6)(2),
                       X(7)(2) };

  realT b[3][8];

  CalculateShapeFunctionDerivative(  y, z, b[0] );
  CalculateShapeFunctionDerivative(  z, x, b[1] );
  CalculateShapeFunctionDerivative(  x, y, b[2] );

  detJ = 0.0;
  for( int a=0 ; a<8 ; ++a )
  {
    detJ += x[a]*b[0][a] + y[a]*b[1][a] + z[a]*b[2][a];
  }
  detJ /= 3.0;

  for( int a=0 ; a<8 ; ++a )
  {
    dNdX(a)(0) = b[0][a] / detJ;
    dNdX(a)(1) = b[1][a] / detJ;
    dNdX(a)(2) = b[2][a] / detJ;
  }

}
