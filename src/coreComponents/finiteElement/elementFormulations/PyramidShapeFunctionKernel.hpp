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
 * @file PyramidShapeFunctionKernel.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_PYRAMIDSHAPEFUNCTIONKERNEL
#define GEOSX_CORE_FINITEELEMENT_PYRAMIDSHAPEFUNCTIONKERNEL

#include "FiniteElementBase.hpp"


namespace geosx
{
namespace finiteElement
{


class PyramidShapeFunctionKernel : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 5;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 5;

  //constexpr static real64 parentVolume = 8.0 / 3.0;
  constexpr static real64 weight = 81.0 / 100.0;
  constexpr static real64 weightDelta  = 125.0 / 27.0 - weight;
  constexpr static real64 quadratureCrossSectionCoord = 0.584237394672177;
  constexpr static real64 quadratureLongitudinalCoordNeg = -2.0 / 3.0;
  constexpr static real64 quadratureLongitudinalCoordDelta = 16.0 / 15.0;

  virtual ~PyramidShapeFunctionKernel() override final
  {}

  virtual localIndex getNumQuadraturePoints() const override final
  {
    return numQuadraturePoints;
  }

  virtual localIndex getNumSupportPoints() const override final
  {
    return numNodes;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void shapeFunctionValues( localIndex const q,
                                   real64 N[numNodes] )
  {
    real64 const xi[3] = { quadratureParentCoords0( q ),
                           quadratureParentCoords1( q ),
                           quadratureParentCoords2( q ) };

    N[0] = 0.125*( 1.0 - xi[0] ) * ( 1.0 - xi[1] ) * ( 1.0 - xi[2] );
    N[1] = 0.125*( 1.0 + xi[0] ) * ( 1.0 - xi[1] ) * ( 1.0 - xi[2] );
    N[2] = 0.125*( 1.0 - xi[0] ) * ( 1.0 + xi[1] ) * ( 1.0 - xi[2] );
    N[3] = 0.125*( 1.0 + xi[0] ) * ( 1.0 + xi[1] ) * ( 1.0 - xi[2] );
    N[4] = 0.5*( 1.0 + xi[2] );
  }

  GEOSX_HOST_DEVICE
  static real64 shapeFunctionDerivatives( localIndex const q,
                                          real64 const (&X)[numNodes][3],
                                          real64 ( &dNdX )[numNodes][3] );

private:
  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static T linearMap( T const i, T const j, T const k )
  {
    return i + 2 * j + 4 * k;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return -1.0 + 2.0 * ( a & 1 ) + 0.25 * ( a & 4 );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return -1.0 + ( a & 2 ) + 0.25 * ( a & 4 );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return -1.0 + 0.5 * ( a & 4 );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords0( localIndex const q )
  {
    return ( -1.0 + 2.0 * ( q & 1 ) + 0.25 * ( q & 4 ) ) * quadratureCrossSectionCoord;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords1( localIndex const q )
  {
    return ( -1.0 + ( q & 2 ) + 0.25 * ( q & 4 ) ) * quadratureCrossSectionCoord;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords2( localIndex const q )
  {
    return quadratureLongitudinalCoordNeg + 0.25 * ( q & 4 ) * quadratureLongitudinalCoordDelta;
  }

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 PyramidShapeFunctionKernel::shapeFunctionDerivatives( localIndex const q,
                                                             real64 const (&X)[numNodes][3],
                                                             real64 (& dNdX)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  real64 const quadratureCoords[3] = { quadratureParentCoords0( q ),
                                       quadratureParentCoords1( q ),
                                       quadratureParentCoords2( q ) };

  real64 const psi0[2] = { 0.5*( 1.0 - quadratureCoords[0] ),
                           0.5*( 1.0 + quadratureCoords[0] ) };
  real64 const psi1[2] = { 0.5*( 1.0 - quadratureCoords[1] ),
                           0.5*( 1.0 + quadratureCoords[1] ) };
  real64 const psi2 = 0.5*( 1.0 - quadratureCoords[2]);
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

  for( localIndex a=0; a<2; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2,
                                psi0[a] * dpsi[b] * psi2,
                                psi0[a] * psi1[b] * dpsi[0] };
      localIndex const nodeIndex = linearMap( a, b, LvArray::integerConversion< localIndex >( 0 ) );
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          J[i][j] = J[i][j] + dNdXi[ j ] * X[nodeIndex][i];
        }
      }
    }
  }

  {
    localIndex const nodeIndex = 4;

    for( int i = 0; i < 3; ++i )
    {
      J[i][2] = J[i][2] + dpsi[1] * X[nodeIndex][i];
    }
  }

  real64 const detJ = inverse( J );

  for( localIndex a=0; a<2; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2,
                                psi0[a] * dpsi[b] * psi2,
                                psi0[a] * psi1[b] * dpsi[0] };
      localIndex const nodeIndex = linearMap( a, b, LvArray::integerConversion< localIndex >( 0 ) );
      for( int i = 0; i < 3; ++i )
      {
        dNdX[nodeIndex][i] = 0.0;
        for( int j = 0; j < 3; ++j )
        {
          dNdX[nodeIndex][i] = dNdX[nodeIndex][i] + dNdXi[ j ] * J[j][i];
        }
      }
    }
  }

  {
    real64 const dNdXi[3] = { 0.0,
                              0.0,
                              dpsi[1] };
    localIndex const nodeIndex = 4;
    for( int i = 0; i < 3; ++i )
    {
      dNdX[nodeIndex][i] = dNdX[nodeIndex][i] + dNdXi[2] * J[2][i];
    }
  }

  // Return determinant times the weight (i.e. for 1-point formula the volume of the tetrahedron)
  return detJ * ( weight + 0.25 * ( q & 4 ) * weightDelta );
}


}
}
#endif //GEOSX_CORE_FINITEELEMENT_PYRAMIDSHAPEFUNCTIONKERNEL
