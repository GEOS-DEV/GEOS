

#ifndef FINITE_ELEMENT_SHAPE_KERNEL
#define FINITE_ELEMENT_SHAPE_KERNEL

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"


namespace geosx
{

class FiniteElementShapeKernel
{
public:
  constexpr static localIndex numNodes = 8;
  constexpr static real64 weight = 8.0 / numNodes;
  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

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


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  parentShapeFunctionValue( localIndex const a,
                            real64 const xi0,
                            real64 const xi1,
                            real64 const xi2 )
  {
    return 0.125 *
           ( 1 + xi0*parentCoords0( a ) ) *
           ( 1 + xi1*parentCoords1( a ) ) *
           ( 1 + xi2*parentCoords2( a ) );
  }


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void shapeFunctionValues( localIndex const q,
                                   real64 N[numNodes] )
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



  #define SUM_FACTORIZATION
#if defined(SUM_FACTORIZATION)
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivatives( localIndex const q,
                                          real64 const (&X)[numNodes][3],
                                          real64 (& dNdX)[numNodes][3] )
  {
    real64 J[3][3] = {{0}};

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



    for( localIndex a=0; a<2; ++a )
    {
      for( localIndex b=0; b<2; ++b )
      {
        for( localIndex c=0; c<2; ++c )
        {
          real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2[c],
                                    psi0[a] * dpsi[b] * psi2[c],
                                    psi0[a] * psi1[b] * dpsi[c] };
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

    real64 const detJ = inverse( J, &(dNdX[0][0]) );


    for( localIndex a=0; a<2; ++a )
    {
      for( localIndex b=0; b<2; ++b )
      {
        for( localIndex c=0; c<2; ++c )
        {
          real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2[c],
                                    psi0[a] * dpsi[b] * psi2[c],
                                    psi0[a] * psi1[b] * dpsi[c] };
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

    return detJ;
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

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 inverse( real64 (& J)[3][3], real64 * GEOSX_RESTRICT const scratch )
  {
    scratch[0] = J[1][1]*J[2][2] - J[1][2]*J[2][1];
    scratch[1] = J[0][2]*J[2][1] - J[0][1]*J[2][2];
    scratch[2] = J[0][1]*J[1][2] - J[0][2]*J[1][1];
    scratch[3] = J[1][2]*J[2][0] - J[1][0]*J[2][2];
    scratch[4] = J[0][0]*J[2][2] - J[0][2]*J[2][0];
    scratch[5] = J[0][2]*J[1][0] - J[0][0]*J[1][2];
    scratch[6] = J[1][0]*J[2][1] - J[1][1]*J[2][0];
    scratch[7] = J[0][1]*J[2][0] - J[0][0]*J[2][1];
    scratch[8] = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    scratch[9] =  J[0][0] * scratch[0] + J[1][0] * scratch[1] + J[2][0] * scratch[2];
    scratch[10] = 1.0 / scratch[9];

    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        J[i][j] = scratch[3*i+j] * scratch[10];
      }
    }

    return scratch[9];
  }



//  typedef real64 const TypeOfdNdX[3][numNodes];
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  real64 shapeFunctionDerivatives( localIndex , localIndex , localIndex a, localIndex i )
//  {
//    return m_dNdX[i][a];
//  }
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  TypeOfdNdX const & shapeFunctionDerivatives(  ) const
//  {
//    return m_dNdX;
//  }


private:
//  real64 m_dNdX[3][numNodes];
//  real64 m_detJ;

};

}

#endif //FINITE_ELEMENT_SHAPE_KERNEL
