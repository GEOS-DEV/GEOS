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
 * @file BiLinearQuadrilateralShapeShapeFunctionKernel.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_BILINEARQUADRILATERALFACE
#define GEOSX_CORE_FINITEELEMENT_BILINEARQUADRILATERALFACE

#include "FiniteElementShapeFunctionKernelBase.hpp"


namespace geosx
{

using BiLinearQuadrilateralFaceBaseClass = FiniteElementShapeFunctionKernelBase< 4, 4 >;

class BiLinearQuadrilateralFaceShapeFunctionKernel : public BiLinearQuadrilateralFaceBaseClass
{
public:
  using BaseClass = BiLinearQuadrilateralFaceBaseClass;

  using BaseClass::numNodes;
  using BaseClass::numQuadraturePoints;
  constexpr static real64 parentArea = 4.0;
  constexpr static real64 weight = parentArea / numQuadraturePoints;
  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void shapeFunctionValues( localIndex const q,
                                   real64 (& N)[numNodes] )
  {
    for( localIndex a=0; a<numNodes; ++a )
    {
      N[a] = 0.25 *
             ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
             ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) );
    }
  }

  GEOSX_HOST_DEVICE
  static real64 JxW( localIndex const q,
                     real64 const (&X)[numNodes][3] );

private:
  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static T linearMap( T const i, T const j )
  {
    return i + 2 * j;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static int basisIndex0( localIndex const a )
  {
    return (a & 1);
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static int basisIndex1( localIndex const a )
  {
    return ( a & 2 ) >> 1;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return -1.0 + 2.0 * (a & 1);
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return -1.0 + ( a & 2 );
  }

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 BiLinearQuadrilateralFaceShapeFunctionKernel::JxW( localIndex const q,
                                                          real64 const (&X)[numNodes][3] )
{
  real64 dXdXi[3][2] = {{0}};

  real64 const quadratureCoords[2] = { quadratureFactor *parentCoords0( q ),
                                       quadratureFactor *parentCoords1( q ) };

  real64 const psi0[2] = { 0.5*( 1.0 - quadratureCoords[0] ),
                           0.5*( 1.0 + quadratureCoords[0] ) };
  real64 const psi1[2] = { 0.5*( 1.0 - quadratureCoords[1] ),
                           0.5*( 1.0 + quadratureCoords[1] ) };
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

  for( int a=0; a<2; ++a )
  {
    for( int b=0; b<2; ++b )
    {
      real64 const dNdXi[2] = { dpsi[a] * psi1[b],
                                psi0[a] * dpsi[b] };

      localIndex const nodeIndex = linearMap( a, b );

      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 2; ++j )
        {
          dXdXi[i][j] = dXdXi[i][j] + dNdXi[ j ] * X[nodeIndex][i];
        }
      }
    }
  }

  real64 const detJ = pow( dXdXi[1][0] * dXdXi[2][1] - dXdXi[2][0] * dXdXi[1][1], 2.0 )
                      + pow( dXdXi[2][0] * dXdXi[0][1] - dXdXi[0][0] * dXdXi[2][1], 2.0 )
                      + pow( dXdXi[0][0] * dXdXi[1][1] - dXdXi[1][0] * dXdXi[0][1], 2.0 );

  return sqrt( detJ ) * weight;
}

}

#endif //GEOSX_CORE_FINITEELEMENT_BILINEARQUADRILATERALFACE
