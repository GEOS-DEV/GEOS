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
 * @file LinearTetrahedronShapeFunctionKernel.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_LINEARTETRAHEDRONSHAPEFUNCTIONKERNEL
#define GEOSX_CORE_FINITEELEMENT_LINEARTETRAHEDRONSHAPEFUNCTIONKERNEL

#include "FiniteElementShapeFunctionKernelBase.hpp"


namespace geosx
{

using LinearTetrahedronShapeFunctionKernelBaseClass = FiniteElementShapeFunctionKernelBase< 4, 1 >;

class LinearTetrahedronShapeFunctionKernel : public LinearTetrahedronShapeFunctionKernelBaseClass
{
public:
  using BaseClass = LinearTetrahedronShapeFunctionKernelBaseClass;

  using BaseClass::numNodes;
  using BaseClass::numQuadraturePoints;
  constexpr static real64 parentVolume = 1.0 / 6.0 ;
  constexpr static real64 weight = parentVolume / numQuadraturePoints; // point weight for 1-point (r = 1/4, s = 1/4, t = 1/4) formula

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void shapeFunctionValues( localIndex const q,
                                   real64 N[numNodes] )
  {
    GEOSX_UNUSED_VAR( q );
    for( localIndex a=0; a<numNodes; ++a )
    {
      N[a] = 0.25;
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivatives( localIndex const q,
                                          real64 const (&X)[numNodes][3],
                                          real64 (& dNdX)[numNodes][3] )
  {
    GEOSX_UNUSED_VAR( q );

    // Jacobian determinant of the isoparametric mapping
    real64 detJ =   ( X[1][0] - X[0][0] )*( ( X[2][1] - X[0][1] )*( X[3][2] - X[0][2] ) - ( X[3][1] - X[0][1] )*( X[2][2] - X[0][2] ) )
                  + ( X[1][1] - X[0][1] )*( ( X[3][0] - X[0][0] )*( X[2][2] - X[0][2] ) - ( X[2][0] - X[0][0] )*( X[3][2] - X[0][2] ) )
                  + ( X[1][2] - X[0][2] )*( ( X[2][0] - X[0][0] )*( X[3][1] - X[0][1] ) - ( X[3][0] - X[0][0] )*( X[2][1] - X[0][1] ) );

    // Shape function derivatives
    dNdX[0][0] =  X[1][1]*( X[3][2] - X[2][2] ) - X[2][1]*( X[3][2] - X[1][2] ) + X[3][1]*( X[2][2] - X[1][2] );
    dNdX[0][1] = -X[1][0]*( X[3][2] - X[2][2] ) + X[2][0]*( X[3][2] - X[1][2] ) - X[3][0]*( X[2][2] - X[1][2] );
    dNdX[0][2] =  X[1][0]*( X[3][1] - X[2][1] ) - X[2][0]*( X[3][1] - X[1][1] ) + X[3][0]*( X[2][1] - X[1][1] );

    dNdX[1][0] = -X[0][1]*( X[3][2] - X[2][2] ) + X[2][1]*( X[3][2] - X[0][2] ) - X[3][1]*( X[2][2] - X[0][2] );
    dNdX[1][1] =  X[0][0]*( X[3][2] - X[2][2] ) - X[2][0]*( X[3][2] - X[0][2] ) + X[3][0]*( X[2][2] - X[0][2] );
    dNdX[1][2] = -X[0][0]*( X[3][1] - X[2][1] ) + X[2][0]*( X[3][1] - X[0][1] ) - X[3][0]*( X[2][1] - X[0][1] );

    dNdX[2][0] =  X[0][1]*( X[3][2] - X[1][2] ) - X[1][1]*( X[3][2] - X[0][2] ) + X[3][1]*( X[1][2] - X[0][2] );
    dNdX[2][1] = -X[0][0]*( X[3][2] - X[1][2] ) + X[1][0]*( X[3][2] - X[0][2] ) - X[3][0]*( X[1][2] - X[0][2] );
    dNdX[2][2] =  X[0][0]*( X[3][1] - X[1][1] ) - X[1][0]*( X[3][1] - X[0][1] ) + X[3][0]*( X[1][1] - X[0][1] );

    dNdX[3][0] = -X[0][1]*( X[2][2] - X[1][2] ) + X[1][1]*( X[2][2] - X[0][2] ) - X[2][1]*( X[1][2] - X[0][2] );
    dNdX[3][1] =  X[0][0]*( X[2][2] - X[1][2] ) - X[1][0]*( X[2][2] - X[0][2] ) + X[2][0]*( X[1][2] - X[0][2] );
    dNdX[3][2] = -X[0][0]*( X[2][1] - X[1][1] ) + X[1][0]*( X[2][1] - X[0][1] ) - X[2][0]*( X[1][1] - X[0][1] );

    real64 factor = 1.0 / ( detJ );
    for( int i = 0; i < numNodes; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        dNdX[i][j] *= factor;
      }
    }

    // Return determinant times the weight (i.e. for 1-point formula the volume of the tetrahedron)
    return detJ * weight;
  }

private:

};

}

#endif //GEOSX_CORE_FINITEELEMENT_LINEARTETRAHEDRONSHAPEFUNCTIONKERNEL
