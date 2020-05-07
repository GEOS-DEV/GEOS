

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



  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords( localIndex const i, localIndex const a )
  {
    constexpr real64 pCoords[3][8] = {
      { -1, 1, -1, 1, -1, 1, -1, 1 },
      { -1, -1, 1, 1, -1, -1, 1, 1 },
      { -1, -1, -1, -1, 1, 1, 1, 1 } };
    return pCoords[i][a];
  }

  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static bool isOdd( T const a )
  {
    return a % 2;
  }

  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static bool isEven( T const a )
  {
    return !(a % 2);
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return isEven( a ) ? -1 : 1;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return isEven( a/2 ) ? -1 : 1;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return isEven( a/4 ) ? -1 : 1;
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
           ( 1 + xi0*parentCoords( 0, a ) ) *
           ( 1 + xi1*parentCoords( 1, a ) ) *
           ( 1 + xi2*parentCoords( 2, a ) );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static void
  parentShapeFunctionDerivatives( localIndex const q,
                                  localIndex const a,
                                  real64 (& dNdXi)[3] )
  {
//    const static real64 pCoords[3][8] = { { -1,  1, -1,  1, -1,  1, -1,  1 },
//                                          { -1, -1,  1,  1, -1, -1,  1,  1 },
//                                          { -1, -1, -1, -1,  1,  1,  1,  1 } };
//
//    dNdXi[0] = 0.125 * pCoords[0][a] *
//                       ( 1 + quadratureFactor*pCoords[1][q]*pCoords[1][a] ) *
//                       ( 1 + quadratureFactor*pCoords[2][q]*pCoords[2][a] ) ;
//    dNdXi[1] = 0.125 * ( 1 + quadratureFactor*pCoords[0][q]*pCoords[0][a] ) *
//                       pCoords[1][a] *
//                       ( 1 + quadratureFactor*pCoords[2][q]*pCoords[2][a] ) ;
//    dNdXi[2] = 0.125 * ( 1 + quadratureFactor*pCoords[0][q]*pCoords[0][a] ) *
//                       ( 1 + quadratureFactor*pCoords[1][q]*pCoords[1][a] ) *
//                       pCoords[2][a] ;

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
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi0( localIndex const q, localIndex const a )
  {
    return 0.125 *
           parentCoords0( a ) *
           ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a )) *
           ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi1( localIndex const q, localIndex const a )
  {
    return 0.125 *
           ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
           parentCoords1( a ) *
           ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
  }


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi2( localIndex const q, localIndex const a )
  {
    return 0.125 *
           ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
           ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) ) *
           parentCoords2( a );
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
  static real64 shapeFunctionDerivatives( localIndex const q,
                                          real64 const (&X)[numNodes][3],
                                          real64 (& dNdX)[numNodes][3] )
  {
    real64 J[3][3] = {{0}};

//#define STOREPARENTSHAPEDER
#ifdef STOREPARENTSHAPEDER
    constexpr static real64 dNdXi[8][8][3] = {
      {
        { dNdXi0( 0, 0 ), dNdXi1( 0, 0 ), dNdXi2( 0, 0 ) },
        { dNdXi0( 0, 1 ), dNdXi1( 0, 1 ), dNdXi2( 0, 1 ) },
        { dNdXi0( 0, 2 ), dNdXi1( 0, 2 ), dNdXi2( 0, 2 ) },
        { dNdXi0( 0, 3 ), dNdXi1( 0, 3 ), dNdXi2( 0, 3 ) },
        { dNdXi0( 0, 4 ), dNdXi1( 0, 4 ), dNdXi2( 0, 4 ) },
        { dNdXi0( 0, 5 ), dNdXi1( 0, 5 ), dNdXi2( 0, 5 ) },
        { dNdXi0( 0, 6 ), dNdXi1( 0, 6 ), dNdXi2( 0, 6 ) },
        { dNdXi0( 0, 7 ), dNdXi1( 0, 7 ), dNdXi2( 0, 7 ) }
      },
      {
        { dNdXi0( 1, 0 ), dNdXi1( 1, 0 ), dNdXi2( 1, 0 ) },
        { dNdXi0( 1, 1 ), dNdXi1( 1, 1 ), dNdXi2( 1, 1 ) },
        { dNdXi0( 1, 2 ), dNdXi1( 1, 2 ), dNdXi2( 1, 2 ) },
        { dNdXi0( 1, 3 ), dNdXi1( 1, 3 ), dNdXi2( 1, 3 ) },
        { dNdXi0( 1, 4 ), dNdXi1( 1, 4 ), dNdXi2( 1, 4 ) },
        { dNdXi0( 1, 5 ), dNdXi1( 1, 5 ), dNdXi2( 1, 5 ) },
        { dNdXi0( 1, 6 ), dNdXi1( 1, 6 ), dNdXi2( 1, 6 ) },
        { dNdXi0( 1, 7 ), dNdXi1( 1, 7 ), dNdXi2( 1, 7 ) }
      },
      {
        { dNdXi0( 2, 0 ), dNdXi1( 2, 0 ), dNdXi2( 2, 0 ) },
        { dNdXi0( 2, 1 ), dNdXi1( 2, 1 ), dNdXi2( 2, 1 ) },
        { dNdXi0( 2, 2 ), dNdXi1( 2, 2 ), dNdXi2( 2, 2 ) },
        { dNdXi0( 2, 3 ), dNdXi1( 2, 3 ), dNdXi2( 2, 3 ) },
        { dNdXi0( 2, 4 ), dNdXi1( 2, 4 ), dNdXi2( 2, 4 ) },
        { dNdXi0( 2, 5 ), dNdXi1( 2, 5 ), dNdXi2( 2, 5 ) },
        { dNdXi0( 2, 6 ), dNdXi1( 2, 6 ), dNdXi2( 2, 6 ) },
        { dNdXi0( 2, 7 ), dNdXi1( 2, 7 ), dNdXi2( 2, 7 ) }
      },
      {
        { dNdXi0( 3, 0 ), dNdXi1( 3, 0 ), dNdXi2( 3, 0 ) },
        { dNdXi0( 3, 1 ), dNdXi1( 3, 1 ), dNdXi2( 3, 1 ) },
        { dNdXi0( 3, 2 ), dNdXi1( 3, 2 ), dNdXi2( 3, 2 ) },
        { dNdXi0( 3, 3 ), dNdXi1( 3, 3 ), dNdXi2( 3, 3 ) },
        { dNdXi0( 3, 4 ), dNdXi1( 3, 4 ), dNdXi2( 3, 4 ) },
        { dNdXi0( 3, 5 ), dNdXi1( 3, 5 ), dNdXi2( 3, 5 ) },
        { dNdXi0( 3, 6 ), dNdXi1( 3, 6 ), dNdXi2( 3, 6 ) },
        { dNdXi0( 3, 7 ), dNdXi1( 3, 7 ), dNdXi2( 3, 7 ) }
      },
      {
        { dNdXi0( 4, 0 ), dNdXi1( 4, 0 ), dNdXi2( 4, 0 ) },
        { dNdXi0( 4, 1 ), dNdXi1( 4, 1 ), dNdXi2( 4, 1 ) },
        { dNdXi0( 4, 2 ), dNdXi1( 4, 2 ), dNdXi2( 4, 2 ) },
        { dNdXi0( 4, 3 ), dNdXi1( 4, 3 ), dNdXi2( 4, 3 ) },
        { dNdXi0( 4, 4 ), dNdXi1( 4, 4 ), dNdXi2( 4, 4 ) },
        { dNdXi0( 4, 5 ), dNdXi1( 4, 5 ), dNdXi2( 4, 5 ) },
        { dNdXi0( 4, 6 ), dNdXi1( 4, 6 ), dNdXi2( 4, 6 ) },
        { dNdXi0( 4, 7 ), dNdXi1( 4, 7 ), dNdXi2( 4, 7 ) }
      },
      {
        { dNdXi0( 5, 0 ), dNdXi1( 5, 0 ), dNdXi2( 5, 0 ) },
        { dNdXi0( 5, 1 ), dNdXi1( 5, 1 ), dNdXi2( 5, 1 ) },
        { dNdXi0( 5, 2 ), dNdXi1( 5, 2 ), dNdXi2( 5, 2 ) },
        { dNdXi0( 5, 3 ), dNdXi1( 5, 3 ), dNdXi2( 5, 3 ) },
        { dNdXi0( 5, 4 ), dNdXi1( 5, 4 ), dNdXi2( 5, 4 ) },
        { dNdXi0( 5, 5 ), dNdXi1( 5, 5 ), dNdXi2( 5, 5 ) },
        { dNdXi0( 5, 6 ), dNdXi1( 5, 6 ), dNdXi2( 5, 6 ) },
        { dNdXi0( 5, 7 ), dNdXi1( 5, 7 ), dNdXi2( 5, 7 ) }
      },
      {
        { dNdXi0( 6, 0 ), dNdXi1( 6, 0 ), dNdXi2( 6, 0 ) },
        { dNdXi0( 6, 1 ), dNdXi1( 6, 1 ), dNdXi2( 6, 1 ) },
        { dNdXi0( 6, 2 ), dNdXi1( 6, 2 ), dNdXi2( 6, 2 ) },
        { dNdXi0( 6, 3 ), dNdXi1( 6, 3 ), dNdXi2( 6, 3 ) },
        { dNdXi0( 6, 4 ), dNdXi1( 6, 4 ), dNdXi2( 6, 4 ) },
        { dNdXi0( 6, 5 ), dNdXi1( 6, 5 ), dNdXi2( 6, 5 ) },
        { dNdXi0( 6, 6 ), dNdXi1( 6, 6 ), dNdXi2( 6, 6 ) },
        { dNdXi0( 6, 7 ), dNdXi1( 6, 7 ), dNdXi2( 6, 7 ) }
      },
      {
        { dNdXi0( 7, 0 ), dNdXi1( 7, 0 ), dNdXi2( 7, 0 ) },
        { dNdXi0( 7, 1 ), dNdXi1( 7, 1 ), dNdXi2( 7, 1 ) },
        { dNdXi0( 7, 2 ), dNdXi1( 7, 2 ), dNdXi2( 7, 2 ) },
        { dNdXi0( 7, 3 ), dNdXi1( 7, 3 ), dNdXi2( 7, 3 ) },
        { dNdXi0( 7, 4 ), dNdXi1( 7, 4 ), dNdXi2( 7, 4 ) },
        { dNdXi0( 7, 5 ), dNdXi1( 7, 5 ), dNdXi2( 7, 5 ) },
        { dNdXi0( 7, 6 ), dNdXi1( 7, 6 ), dNdXi2( 7, 6 ) },
        { dNdXi0( 7, 7 ), dNdXi1( 7, 7 ), dNdXi2( 7, 7 ) }
      }
    };

#define DNDXI( q, a, j ) dNdXi[q][a][j]
#define CALC_DNDXI( q, a )
#else
    real64 dNdXi[3];
#define CALC_DNDXI( q, a ) parentShapeFunctionDerivatives( q, a, dNdXi )
#define DNDXI( q, a, i ) dNdXi[i]
#endif


    for( localIndex a=0; a<numNodes; ++a )
    {
      CALC_DNDXI( q, a );
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          J[i][j] = J[i][j] + X[a][i] * DNDXI( q, a, j );
        }
      }
    }

    real64 const invDetJ = inverse( J, &(dNdX[0][0]) );

    for( localIndex a=0; a<numNodes; ++a )
    {
      CALC_DNDXI( q, a );
      for( int i = 0; i < 3; ++i )
      {
        dNdX[a][i] = 0.0;
        for( int j = 0; j < 3; ++j )
        {
          dNdX[a][i] = dNdX[a][i] + DNDXI( q, a, j ) * J[j][i];
        }
      }
    }

    return 1.0 / invDetJ;
#undef CALC_DNDXI
#undef DNDXI
  }

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

    scratch[9] = 1 / ( J[0][0] * scratch[0] + J[1][0] * scratch[1] + J[2][0] * scratch[2] );

    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        J[i][j] = scratch[3*i+j] * scratch[9];
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
