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
 * File: ShapeFunctionBaseT.h
 * Shape function Class
 *
 * created : RRS (09/14/2010)
 */

#include "ShapeFunctionBaseT.h"

ShapeFunctionBaseT::ShapeT() :
  rElement( 0 ),
  pNodes( 0 )
{}

ShapeFunctionBaseT::~ShapeFunctionBaseT( void )
{}



void ShapeFunctionBaseT::CalculateJacobian( const int elem )
{
  J = 0.0;
  for( int i=1; i<=nsdof; ++i )
    for( int j=1; j<=nsdof; ++j )
      for( int a=1; a<=rElement->NumNodeElem(); ++a )
        J( i, j ) += pNodes->Xref( rElement->Connectivity( elem, a ), i )
                     * dNdXi( a )( j );

}

void ContinuumShapeT::Calc_Shape_Deriv( const realT fac )
{

  realT factor = 1.0 / pow( static_cast< realT >(2), nsdof );

  // loop over "nodal indexed" shape functions and calculate the
  // shape function derivitives wrt the local coordinates
  for( int a=1; a<=rElement->NumNodeElem(); ++a )
  {
    for( int i=1; i<=nsdof; ++i )
    {
      dNdXi( a )( i ) = factor * rElement->Xi_node( a, i );
      for( int j=1; j<=nsdof; j++ )
        if( i!=j )
          dNdXi( a )( i ) *= ( 1 + ip_coord_fac * rElement->Xi_node( ip, j ) * rElement->Xi_node( a, j ));

    }
  }

  CalculateJacobian( elem );

  m_detJ( elem, ip ) = J.Det();

  dNdX.AijBi( Jinv, dNdXi );


}


inline realT ShapeFunctionBaseT::Shape( const R1Tensor & Xi,
                                        const R1Tensor & Xi_node )
{
  realT N = 1.0 / pow( 2.0, nsdof );

  for( int i=1; i<=nsdof; ++i )
    N *= ( 1 +  Xi( i ) * Xi_node( i ) );

  return N;
}
