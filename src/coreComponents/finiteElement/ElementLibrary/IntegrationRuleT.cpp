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
 * @file IntegrationPointsT.cpp
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


void IntegrationRuleT::CalculateShapeFunctionDerivatives( const array1d< R1Tensor > & X,
                                                          array1d< R1Tensor > & dNdX,
                                                          realT & detJ )
{
  const realT x[8] = { X( 0 )( 0 ),
                       X( 1 )( 0 ),
                       X( 2 )( 0 ),
                       X( 3 )( 0 ),
                       X( 4 )( 0 ),
                       X( 5 )( 0 ),
                       X( 6 )( 0 ),
                       X( 7 )( 0 ) };

  const realT y[8] = { X( 0 )( 1 ),
                       X( 1 )( 1 ),
                       X( 2 )( 1 ),
                       X( 3 )( 1 ),
                       X( 4 )( 1 ),
                       X( 5 )( 1 ),
                       X( 6 )( 1 ),
                       X( 7 )( 1 ) };

  const realT z[8] = { X( 0 )( 2 ),
                       X( 1 )( 2 ),
                       X( 2 )( 2 ),
                       X( 3 )( 2 ),
                       X( 4 )( 2 ),
                       X( 5 )( 2 ),
                       X( 6 )( 2 ),
                       X( 7 )( 2 ) };

  realT b[3][8];

  CalculateShapeFunctionDerivative( y, z, b[0] );
  CalculateShapeFunctionDerivative( z, x, b[1] );
  CalculateShapeFunctionDerivative( x, y, b[2] );

  detJ = 0.0;
  for( int a=0; a<8; ++a )
  {
    detJ += x[a]*b[0][a] + y[a]*b[1][a] + z[a]*b[2][a];
  }
  detJ /= 3.0;

  for( int a=0; a<8; ++a )
  {
    dNdX( a )( 0 ) = b[0][a] / detJ;
    dNdX( a )( 1 ) = b[1][a] / detJ;
    dNdX( a )( 2 ) = b[2][a] / detJ;
  }

}
