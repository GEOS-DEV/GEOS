

#ifndef GEOSX_CORE_FINITEELEMENT_TRILINEARHEXAHEDRON
#define GEOSX_CORE_FINITEELEMENT_TRILINEARHEXAHEDRON

#include "FiniteElementShapeFunctionKernelBase.hpp"


namespace geosx
{

using TrilinearHexahedronBaseClass = FiniteElementShapeFunctionKernelBase< 8, 8 >;

class TrilinearHexahedronShapeFunctionKernel : public TrilinearHexahedronBaseClass
{
public:
  using BaseClass = TrilinearHexahedronBaseClass;

  using BaseClass::numNodes;
  using BaseClass::numQuadraturePoints;
  constexpr static real64 parentVolume = 8.0;
  constexpr static real64 weight = parentVolume / numQuadraturePoints;
  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void shapeFunctionValues( localIndex const q,
                                   real64 (&N)[numNodes] )
  {
    for( localIndex a=0; a<numNodes; ++a )
    {
      N[a] = 0.125 *
             ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
             ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) ) *
             ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static void
  parentShapeFunctionDerivatives( localIndex const q,
                                  localIndex const a,
                                  real64 (& dNdXi)[3] )
  {
    dNdXi[0] = 0.125 *
               parentCoords0( a ) *
               ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) ) *
               ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
    dNdXi[1] = 0.125 *
               ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
               parentCoords1( a ) *
               ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
    dNdXi[2] = 0.125 *
               ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
               ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) ) *
               parentCoords2( a );
  }



  GEOSX_HOST_DEVICE
  static real64 shapeFunctionDerivatives( localIndex const q,
                                          real64 const (&X)[numNodes][3],
                                          real64 (& dNdX)[numNodes][3] );

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
  constexpr static int basisIndex2( localIndex const a )
  {
    return ( a & 4 ) >> 2;
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

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return -1.0 + 0.5 * ( a & 4 );
  }


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 oneDimensionalShape( real64 const parentCoord,
                                               real64 const coord )
  {
    return 0.5 * ( 1.0 + parentCoord * coord );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 oneDimensionalShapeHalf( real64 const halfParentCoord,
                                                   real64 const coord )
  {
    return 0.5 + halfParentCoord * coord;
  }

};


#define SUM_FACTORIZATION
#if defined(SUM_FACTORIZATION)
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 TrilinearHexahedronShapeFunctionKernel::shapeFunctionDerivatives( localIndex const q,
                                                           real64 const (&X)[numNodes][3],
                                                           real64 (& dNdX)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

#define PD 3

#if PD==1
  real64 const quadratureCoords[3] = { quadratureFactor *parentCoords0( q ),
                                       quadratureFactor *parentCoords1( q ),
                                       quadratureFactor *parentCoords2( q ) };

  real64 const psi0[2] = { oneDimensionalShapeHalf( -0.5, quadratureCoords[0] ),
                           oneDimensionalShapeHalf( 0.5, quadratureCoords[0] ) };
  real64 const psi1[2] = { oneDimensionalShapeHalf( -0.5, quadratureCoords[1] ),
                           oneDimensionalShapeHalf( 0.5, quadratureCoords[1] ) };
  real64 const psi2[2] = { oneDimensionalShapeHalf( -0.5, quadratureCoords[2] ),
                           oneDimensionalShapeHalf( 0.5, quadratureCoords[2] ) };
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

#elif PD==2
  constexpr static real64 linearBasisAtQuadrature[2] = { 0.5 * ( 1 + quadratureFactor ),
                                                         0.5 * ( 1 - quadratureFactor ) };
  real64 const psi0[2] = { linearBasisAtQuadrature[basisIndex0(q)], linearBasisAtQuadrature[!basisIndex0(q)] };
  real64 const psi1[2] = { linearBasisAtQuadrature[basisIndex1(q)], linearBasisAtQuadrature[!basisIndex1(q)] };
  real64 const psi2[2] = { linearBasisAtQuadrature[basisIndex2(q)], linearBasisAtQuadrature[!basisIndex2(q)] };
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

#elif PD==3
  int const qa = basisIndex0(q);
  int const qb = basisIndex1(q);
  int const qc = basisIndex2(q);

  constexpr static real64 linearBasisAtQuadrature[2] = { 0.5 + 0.5 * quadratureFactor,
                                                         0.5 - 0.5 * quadratureFactor };
  constexpr static real64 psiProduct[3] = { 0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[0],
                                            0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[1],
                                            0.5 * linearBasisAtQuadrature[1]*linearBasisAtQuadrature[1] };

  constexpr static int dpsi[2] = { -1, 1 };

#endif


  for( int a=0; a<2; ++a )
  {
#if PD==3
    int const qaa = abs(a-qa);
#endif
    for( int b=0; b<2; ++b )
    {
#if PD==3
      int const qbb = abs(b-qb);
#endif
      for( int c=0; c<2; ++c )
      {
#if PD<3
        real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2[c],
                                  psi0[a] * dpsi[b] * psi2[c],
                                  psi0[a] * psi1[b] * dpsi[c] };
#else
        int const qcc = abs(c-qc);
        real64 const dNdXi[3] = { dpsi[a] * psiProduct[ qbb + qcc ],
                                   dpsi[b] * psiProduct[ qaa + qcc ],
                                   dpsi[c] * psiProduct[ qaa + qbb ] };

#endif

        localIndex const nodeIndex = linearMap( a, b, c );

        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            J[i][j] = J[i][j] + dNdXi[ j ] * X[nodeIndex][i];
          }
        }
      }
    }
  }

  real64 const detJ = inverse( J );


  for( int a=0; a<2; ++a )
  {
#if PD==3
    int const qaa = abs(a-qa);
#endif
    for( int b=0; b<2; ++b )
    {
#if PD==3
      int const qbb = abs(b-qb);
#endif
      for( int c=0; c<2; ++c )
      {
#if PD<3
        real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2[c],
                                  psi0[a] * dpsi[b] * psi2[c],
                                  psi0[a] * psi1[b] * dpsi[c] };
#else
        int const qcc = abs(c-qc);
        real64 const dNdXi[3] = { dpsi[a] * psiProduct[ qbb + qcc ],
                                  dpsi[b] * psiProduct[ qaa + qcc ],
                                  dpsi[c] * psiProduct[ qaa + qbb ] };
#endif
        localIndex const nodeIndex = linearMap( a, b, c );
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
  }

  return detJ * weight;
}

#else

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
static real64 shapeFunctionDerivatives( localIndex const q,
                                        real64 const (&X)[numNodes][3],
                                        real64 (& dNdX)[numNodes][3] )
{
  real64 J[3][3] = {{0}};
  real64 dNdXi[3];

  for( localIndex a=0; a<numNodes; ++a )
  {
    parentShapeFunctionDerivatives( q, a, dNdXi );
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        J[i][j] = J[i][j] + X[a][i] * dNdXi[ j ];
      }
    }
  }

  real64 const invDetJ = inverse( J, &(dNdX[0][0]) );

  for( localIndex a=0; a<numNodes; ++a )
  {
    parentShapeFunctionDerivatives( q, a, dNdXi );
    for( int i = 0; i < 3; ++i )
    {
      dNdX[a][i] = 0.0;
      for( int j = 0; j < 3; ++j )
      {
        dNdX[a][i] = dNdX[a][i] + dNdXi[ j ] * J[j][i];
      }
    }
  }
  return 1.0 / invDetJ;
}
#endif


}

#endif //FINITE_ELEMENT_SHAPE_KERNEL
