

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
    constexpr  real64 pCoords[3][8] = { { -1,  1, -1,  1, -1,  1, -1,  1 },
                                        { -1, -1,  1,  1, -1, -1,  1,  1 },
                                        { -1, -1, -1, -1,  1,  1,  1,  1 } };
    return pCoords[i][a];
  }

  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static bool isOdd( T const a )
  {
    return a % 2;
  };

  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static bool isEven( T const a )
  {
    return !(a % 2);
  };

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return isEven(a) ? -1 : 1;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return isEven(a/2) ? -1 : 1;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return isEven(a/4) ? -1 : 1;
  }


  // GEOSX_HOST_DEVICE
  // GEOSX_FORCE_INLINE
  // static real64
  // parentShapeFunctionValue( localIndex const a,
  //                           real64 const xi0,
  //                           real64 const xi1,
  //                           real64 const xi2 )
  // {
  //   return 0.125 * ( 1 + xi0*parentCoords(0,a) ) * ( 1 + xi1*parentCoords(1,a) ) * ( 1 + xi2*parentCoords(2,a) ) ;
  // }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void parentShapeFunctionDerivatives( localIndex const a,
                                              real64 const xi0,
                                              real64 const xi1,
                                              real64 const xi2,
                                              real64 (&dNdXi)[3] )
  {
    real64 const pC_0_a = parentCoords(0,a);
    real64 const pC_1_a = parentCoords(1,a);
    real64 const pC_2_a = parentCoords(2,a);

    dNdXi[0] = 0.125 * pC_0_a * ( 1 + xi1*pC_1_a ) * ( 1 + xi2*pC_2_a ) ;
    dNdXi[1] = 0.125 * ( 1 + xi0*pC_0_a ) * pC_1_a * ( 1 + xi2*pC_2_a ) ;
    dNdXi[2] = 0.125 * ( 1 + xi0*pC_0_a ) * ( 1 + xi1*pC_1_a ) * pC_2_a ;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi0( localIndex const q, localIndex const a )
  {
    return 0.125 * parentCoords0(a) *
                   ( 1 + quadratureFactor*parentCoords1(q)*parentCoords1(a)) *
                   ( 1 + quadratureFactor*parentCoords2(q)*parentCoords2(a) ) ;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi1( localIndex const q, localIndex const a )
  {
    return 0.125 * ( 1 + quadratureFactor*parentCoords0(q)*parentCoords0(a) ) *
                   parentCoords1(a) *
                   ( 1 + quadratureFactor*parentCoords2(q)*parentCoords2(a) ) ;
  }


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi2( localIndex const q, localIndex const a )
  {
    return 0.125 * ( 1 + quadratureFactor*parentCoords0(q)*parentCoords0(a) ) *
                   ( 1 + quadratureFactor*parentCoords1(q)*parentCoords1(a) ) *
                   parentCoords2(a) ;
  }




  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivatives( localIndex const k,
                                          localIndex const q,
                                          arrayView2d<localIndex const> const & elemsToNodes,
                                          arrayView1d<R1Tensor const> const & X,
                                          real64 (&dNdX)[3][numNodes] )
  {
    real64 J[3][3] = {{0}};

    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      for ( int i = 0; i < 3; ++i )
      {
        J[i][0] += X[elemsToNodes(k, a)][i] * dNdXi0(q,a);
        J[i][1] += X[elemsToNodes(k, a)][i] * dNdXi1(q,a);
        J[i][2] += X[elemsToNodes(k, a)][i] * dNdXi2(q,a);
//        #pragma unroll
//        for ( int j = 0; j < 3; ++j )
//        {
//          J[i][j] += X[elemsToNodes(k, a)][i] * dNdXi[q][a][j];
//        }
      }
    }

    real64 const detJ = inverse( J );

    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      for ( int i = 0; i < 3; ++i )
      {
//        dNdX[i][a] = J[0][i] * dNdXi[q][a][0] +
//                     J[1][i] * dNdXi[q][a][1] +
//                     J[2][i] * dNdXi[q][a][2];
        dNdX[i][a] = J[0][i] * dNdXi0(q,a) +
                     J[1][i] * dNdXi1(q,a) +
                     J[2][i] * dNdXi2(q,a);
      }
    }

    return detJ;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 inverse( real64 (&J)[3][3] )
  {
    real64 const o1 = J[1][1]*J[2][2] - J[1][2]*J[2][1];
    real64 const o2 = J[0][2]*J[2][1] - J[0][1]*J[2][2];
    real64 const o3 = J[0][1]*J[1][2] - J[0][2]*J[1][1];
    real64 const o4 = J[1][2]*J[2][0] - J[1][0]*J[2][2];
    real64 const o5 = J[0][0]*J[2][2] - J[0][2]*J[2][0];
    real64 const o6 = J[0][2]*J[1][0] - J[0][0]*J[1][2];
    real64 const o7 = J[1][0]*J[2][1] - J[1][1]*J[2][0];
    real64 const o8 = J[0][1]*J[2][0] - J[0][0]*J[2][1];
    real64 const o9 = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    real64 const detJ = J[0][0] * o1 + J[1][0] * o2 + J[2][0] * o3 ;

    J[0][0] = o1 / detJ;
    J[0][1] = o2 / detJ;
    J[0][2] = o3 / detJ;
    J[1][0] = o4 / detJ;
    J[1][1] = o5 / detJ;
    J[1][2] = o6 / detJ;
    J[2][0] = o7 / detJ;
    J[2][1] = o8 / detJ;
    J[2][2] = o9 / detJ;

    return detJ;
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
