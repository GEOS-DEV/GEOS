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
 * @file BiLinearWedgeShapeFunctionKernel.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_BILINEARWEDGESHAPEFUNCTIONKERNEL
#define GEOSX_CORE_FINITEELEMENT_BILINEARWEDGESHAPEFUNCTIONKERNEL

#include "FiniteElementBase.hpp"


namespace geosx
{
namespace finiteElement
{


class BiLinearWedgeShapeFunctionKernel : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 6;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 6;

  constexpr static real64 parentVolume = 1.0;
  constexpr static real64 weight = parentVolume / numQuadraturePoints;
  constexpr static real64 quadratureCrossSectionCoord = 1.0 / 6.0;
  constexpr static real64 quadratureLongitudinalCoord = 1.0 / 1.732050807568877293528;


  virtual ~BiLinearWedgeShapeFunctionKernel() override final
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
                                   real64 (&N)[numNodes] )
  {
    real64 const xi[3] = { quadratureParentCoords0( q ),
                           quadratureParentCoords1( q ),
                           quadratureParentCoords2( q ) };

    N[0] = 0.5*( 1.0 - xi[0] - xi[1] ) * ( 1.0 - xi[2] );
    N[1] = 0.5*( 1.0 - xi[0] - xi[1] ) * ( 1.0 + xi[2] );
    N[2] = 0.5*( xi[0] ) * ( 1.0 - xi[2] );
    N[3] = 0.5*( xi[0] ) * ( 1.0 + xi[2] );
    N[4] = 0.5*( xi[1] ) * ( 1.0 - xi[2] );
    N[5] = 0.5*( xi[1] ) * ( 1.0 + xi[2] );
  }

  GEOSX_HOST_DEVICE
  static real64 shapeFunctionDerivatives( localIndex const q,
                                          real64 const (&X)[numNodes][3],
                                          real64 ( &dNdX )[numNodes][3] );

private:
  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static T linearMap( T const indexT, T const indexL )
  {
    return 2 * indexT + indexL;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return 0.5 * ( a & 2 );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return 0.25 * ( a & 4 );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return -1.0 + 2 * ( a & 1 );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords0( localIndex const q )
  {
    return ( 0.25 * ( q & 2 ) + quadratureCrossSectionCoord );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords1( localIndex const q )
  {
    return ( 0.125 * ( q & 4 ) + quadratureCrossSectionCoord );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords2( localIndex const q )
  {
    return parentCoords2( q ) * quadratureLongitudinalCoord;
  }


};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 BiLinearWedgeShapeFunctionKernel::shapeFunctionDerivatives( localIndex const q,
                                                                   real64 const (&X)[numNodes][3],
                                                                   real64 (& dNdX)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  real64 const quadratureCoords[3] = { quadratureParentCoords0( q ),
                                       quadratureParentCoords1( q ),
                                       quadratureParentCoords2( q ) };

  real64 const psiTRI[3] = { 1.0 - quadratureCoords[0]  - quadratureCoords[1],
                             quadratureCoords[0],
                             quadratureCoords[1] };
  real64 const psiLIN[2] = { 0.5 - 0.5*quadratureCoords[2],
                             0.5 + 0.5*quadratureCoords[2] };
  constexpr real64 dpsiTRI[2][3] = { { -1.0, 1.0, 0.0 },
    { -1.0, 0.0, 1.0 }
  };
  constexpr real64 dpsiLIN[2] = { -0.5, 0.5 };

  for( localIndex a=0; a<3; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsiTRI[0][a] * psiLIN[b],
                                dpsiTRI[1][a] * psiLIN[b],
                                psiTRI[a] * dpsiLIN[b] };
      localIndex const nodeIndex = linearMap( a, b );
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          J[i][j] = J[i][j] + dNdXi[ j ] * X[nodeIndex][i];
        }
      }
    }
  }

  real64 const detJ = inverse( J );


  for( localIndex a=0; a<3; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsiTRI[0][a] * psiLIN[b],
                                dpsiTRI[1][a] * psiLIN[b],
                                psiTRI[a] * dpsiLIN[b] };
      localIndex const nodeIndex = linearMap( a, b );
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

  // Return determinant times the weight (i.e. for 1-point formula the volume of the tetrahedron)
  return detJ * weight;
}

}
}

#endif //GEOSX_CORE_FINITEELEMENT_BILINEARWEDGESHAPEFUNCTIONKERNEL
