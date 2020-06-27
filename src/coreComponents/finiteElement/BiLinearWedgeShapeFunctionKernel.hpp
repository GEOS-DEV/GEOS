

#ifndef SIX_NODE_WEDGE_SHAPE_FUNCTION_KERNEL
#define SIX_NODE_WEDGE_SHAPE_FUNCTION_KERNEL

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"


namespace geosx
{

class BiLinearWedgeShapeFunctionKernel
{
public:
  constexpr static localIndex numNodes = 6;
  constexpr static real64 weight = 1.0 / 6.0;
  constexpr static real64 quadratureCrossSectionCoord = 1.0 / 6.0;
  constexpr static real64 quadratureLongitudinalCoord = 1.0 / 1.732050807568877293528;

  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static T linearMap( T const i, T const j )
  {
    return i + 3 * j;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords0( localIndex const q )
  {
    return ( 1.0 + 3.0*static_cast< double >( ( ( q % 3 ) == 1 ) ) ) * quadratureCrossSectionCoord;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords1( localIndex const q )
  {
    return ( 1.0 + 3.0*static_cast< double >( ( ( q % 3 ) == 2 ) ) ) * quadratureCrossSectionCoord;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords2( localIndex const q )
  {
    return ( -1.0 + 2 * static_cast< double >( q / 3 ) ) * quadratureLongitudinalCoord;
  }

//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  constexpr static real64 oneDimensionalShape( real64 const parentCoord,
//                                               real64 const coord )
//  {
//    return 0.5 * ( 1.0 + parentCoord * coord );
//  }
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  constexpr static real64 oneDimensionalShapeHalf( real64 const halfParentCoord,
//                                                   real64 const coord )
//  {
//    return 0.5 + halfParentCoord * coord;
//  }
//
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  constexpr static real64
//  parentShapeFunctionValue( localIndex const a,
//                            real64 const xi0,
//                            real64 const xi1,
//                            real64 const xi2 )
//  {
//    return 0.125 *
//           ( 1 + xi0*parentCoords0( a ) ) *
//           ( 1 + xi1*parentCoords1( a ) ) *
//           ( 1 + xi2*parentCoords2( a ) );
//  }


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void shapeFunctionValues( localIndex const q,
                                   real64 N[numNodes] )
  {
    real64 const xi[3] = { quadratureParentCoords0( q ),
                           quadratureParentCoords1( q ),
                           quadratureParentCoords2( q ) };

    N[0] = 0.5*( 1.0 - xi[0] - xi[1] ) * ( 1.0 - xi[2] );
    N[1] = 0.5*( xi[0] ) * ( 1.0 - xi[2] );
    N[2] = 0.5*( xi[1] ) * ( 1.0 - xi[2] );
    N[3] = 0.5*( 1.0 - xi[0] - xi[1] ) * ( 1.0 + xi[2] );
    N[4] = 0.5*( xi[0] ) * ( 1.0 + xi[2] );
    N[5] = 0.5*( xi[1] ) * ( 1.0 + xi[2] );
  }

//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  constexpr static void
//  parentShapeFunctionDerivatives( localIndex const q,
//                                  localIndex const a,
//                                  real64 (& dNdXi)[3] )
//  {
//    GEOSX_UNUSED_VAR( a );
//    real64 const xi[3] = { quadratureParentCoords2( q ),
//                           quadratureParentCoords2( q ),
//                           quadratureParentCoords2( q ) };
//    real64 c0 = static_cast< double >( ( q % 3 ) == 0 );
//    real64 c1 = -c0 + static_cast< double >( ( q % 3 ) == 1 );
//    real64 c2 = -c0 + static_cast< double >( ( q % 3 ) == 2 );
//    real64 c3 = ( - 1.0 + 2 * static_cast< double >( q / 3 ) );
//    dNdXi[0] = 0.5*c1*( 1 + c3*xi[2] );
//    dNdXi[1] = 0.5*c2*( 1 + c3*xi[2] );
//    dNdXi[2] = 0.5*( c0 + c1*xi[0] + c2*xi[3] )*c3;
//  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivatives( localIndex const q,
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

    real64 const detJ = inverse( J, &(dNdX[0][0]) );


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

#endif //SIX_NODE_WEDGE_SHAPE_FUNCTION_KERNEL
