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
  constexpr static real64 parentVolume = 4.0;
  constexpr static real64 weight = parentVolume / numQuadraturePoints;
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

//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  constexpr static int basisIndex0( localIndex const a )
//  {
//    return (a & 1);
//  }
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  constexpr static int basisIndex1( localIndex const a )
//  {
//    return ( a & 2 ) >> 1;
//  }
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  constexpr static int basisIndex2( localIndex const a )
//  {
//    return ( a & 4 ) >> 2;
//  }


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
