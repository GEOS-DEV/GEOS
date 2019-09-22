

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

  template< int q, int a >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi0()
  {
    constexpr real64 xi1 = quadratureFactor*parentCoords(1,q);
    constexpr real64 xi2 = quadratureFactor*parentCoords(2,q);

    real64 const pC_0_a = parentCoords(0,a);
    real64 const pC_1_a = parentCoords(1,a);
    real64 const pC_2_a = parentCoords(2,a);

    return 0.125 * pC_0_a * ( 1 + xi1*pC_1_a ) * ( 1 + xi2*pC_2_a ) ;
  }

  template< int q, int a >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi1()
  {
    constexpr real64 xi0 = quadratureFactor*parentCoords(0,q);
    constexpr real64 xi2 = quadratureFactor*parentCoords(2,q);

    real64 const pC_0_a = parentCoords(0,a);
    real64 const pC_1_a = parentCoords(1,a);
    real64 const pC_2_a = parentCoords(2,a);

    return 0.125 * ( 1 + xi0*pC_0_a ) * pC_1_a * ( 1 + xi2*pC_2_a ) ;
  }


  template< int q, int a >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 dNdXi2()
  {
    constexpr real64 xi0 = quadratureFactor*parentCoords(0,q);
    constexpr real64 xi1 = quadratureFactor*parentCoords(1,q);

    real64 const pC_0_a = parentCoords(0,a);
    real64 const pC_1_a = parentCoords(1,a);
    real64 const pC_2_a = parentCoords(2,a);

    return 0.125 * ( 1 + xi0*pC_0_a ) * ( 1 + xi1*pC_1_a ) * pC_2_a ;
  }




  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivativesImpl( localIndex const k,
                                          localIndex const q,
                                          arrayView2d<localIndex const> const & elemsToNodes,
//                                          real64 const xi0,
//                                          real64 const xi1,
//                                          real64 const xi2,
                                          arrayView1d<R1Tensor const> const & X,
                                          real64 (&dNdX)[3][numNodes] )
  {
    real64 J[3][3] = {{0}};
    constexpr static real64 dNdXi[8][8][3] = { { { dNdXi0<0,0>(), dNdXi1<0,0>(), dNdXi2<0,0>() },
                                          { dNdXi0<0,1>(), dNdXi1<0,1>(), dNdXi2<0,1>() },
                                          { dNdXi0<0,2>(), dNdXi1<0,2>(), dNdXi2<0,2>() },
                                          { dNdXi0<0,3>(), dNdXi1<0,3>(), dNdXi2<0,3>() },
                                          { dNdXi0<0,4>(), dNdXi1<0,4>(), dNdXi2<0,4>() },
                                          { dNdXi0<0,5>(), dNdXi1<0,5>(), dNdXi2<0,5>() },
                                          { dNdXi0<0,6>(), dNdXi1<0,6>(), dNdXi2<0,6>() },
                                          { dNdXi0<0,7>(), dNdXi1<0,7>(), dNdXi2<0,7>() }
                                        },
                                        { { dNdXi0<1,0>(), dNdXi1<1,0>(), dNdXi2<1,0>() },
                                          { dNdXi0<1,1>(), dNdXi1<1,1>(), dNdXi2<1,1>() },
                                          { dNdXi0<1,2>(), dNdXi1<1,2>(), dNdXi2<1,2>() },
                                          { dNdXi0<1,3>(), dNdXi1<1,3>(), dNdXi2<1,3>() },
                                          { dNdXi0<1,4>(), dNdXi1<1,4>(), dNdXi2<1,4>() },
                                          { dNdXi0<1,5>(), dNdXi1<1,5>(), dNdXi2<1,5>() },
                                          { dNdXi0<1,6>(), dNdXi1<1,6>(), dNdXi2<1,6>() },
                                          { dNdXi0<1,7>(), dNdXi1<1,7>(), dNdXi2<1,7>() }
                                        },
                                        { { dNdXi0<2,0>(), dNdXi1<2,0>(), dNdXi2<2,0>() },
                                          { dNdXi0<2,1>(), dNdXi1<2,1>(), dNdXi2<2,1>() },
                                          { dNdXi0<2,2>(), dNdXi1<2,2>(), dNdXi2<2,2>() },
                                          { dNdXi0<2,3>(), dNdXi1<2,3>(), dNdXi2<2,3>() },
                                          { dNdXi0<2,4>(), dNdXi1<2,4>(), dNdXi2<2,4>() },
                                          { dNdXi0<2,5>(), dNdXi1<2,5>(), dNdXi2<2,5>() },
                                          { dNdXi0<2,6>(), dNdXi1<2,6>(), dNdXi2<2,6>() },
                                          { dNdXi0<2,7>(), dNdXi1<2,7>(), dNdXi2<2,7>() }
                                        },
                                        { { dNdXi0<3,0>(), dNdXi1<3,0>(), dNdXi2<3,0>() },
                                          { dNdXi0<3,1>(), dNdXi1<3,1>(), dNdXi2<3,1>() },
                                          { dNdXi0<3,2>(), dNdXi1<3,2>(), dNdXi2<3,2>() },
                                          { dNdXi0<3,3>(), dNdXi1<3,3>(), dNdXi2<3,3>() },
                                          { dNdXi0<3,4>(), dNdXi1<3,4>(), dNdXi2<3,4>() },
                                          { dNdXi0<3,5>(), dNdXi1<3,5>(), dNdXi2<3,5>() },
                                          { dNdXi0<3,6>(), dNdXi1<3,6>(), dNdXi2<3,6>() },
                                          { dNdXi0<3,7>(), dNdXi1<3,7>(), dNdXi2<3,7>() }
                                        },
                                        { { dNdXi0<4,0>(), dNdXi1<4,0>(), dNdXi2<4,0>() },
                                          { dNdXi0<4,1>(), dNdXi1<4,1>(), dNdXi2<4,1>() },
                                          { dNdXi0<4,2>(), dNdXi1<4,2>(), dNdXi2<4,2>() },
                                          { dNdXi0<4,3>(), dNdXi1<4,3>(), dNdXi2<4,3>() },
                                          { dNdXi0<4,4>(), dNdXi1<4,4>(), dNdXi2<4,4>() },
                                          { dNdXi0<4,5>(), dNdXi1<4,5>(), dNdXi2<4,5>() },
                                          { dNdXi0<4,6>(), dNdXi1<4,6>(), dNdXi2<4,6>() },
                                          { dNdXi0<4,7>(), dNdXi1<4,7>(), dNdXi2<4,7>() }
                                        },
                                        { { dNdXi0<5,0>(), dNdXi1<5,0>(), dNdXi2<5,0>() },
                                          { dNdXi0<5,1>(), dNdXi1<5,1>(), dNdXi2<5,1>() },
                                          { dNdXi0<5,2>(), dNdXi1<5,2>(), dNdXi2<5,2>() },
                                          { dNdXi0<5,3>(), dNdXi1<5,3>(), dNdXi2<5,3>() },
                                          { dNdXi0<5,4>(), dNdXi1<5,4>(), dNdXi2<5,4>() },
                                          { dNdXi0<5,5>(), dNdXi1<5,5>(), dNdXi2<5,5>() },
                                          { dNdXi0<5,6>(), dNdXi1<5,6>(), dNdXi2<5,6>() },
                                          { dNdXi0<5,7>(), dNdXi1<5,7>(), dNdXi2<5,7>() }
                                        },
                                        { { dNdXi0<6,0>(), dNdXi1<6,0>(), dNdXi2<6,0>() },
                                          { dNdXi0<6,1>(), dNdXi1<6,1>(), dNdXi2<6,1>() },
                                          { dNdXi0<6,2>(), dNdXi1<6,2>(), dNdXi2<6,2>() },
                                          { dNdXi0<6,3>(), dNdXi1<6,3>(), dNdXi2<6,3>() },
                                          { dNdXi0<6,4>(), dNdXi1<6,4>(), dNdXi2<6,4>() },
                                          { dNdXi0<6,5>(), dNdXi1<6,5>(), dNdXi2<6,5>() },
                                          { dNdXi0<6,6>(), dNdXi1<6,6>(), dNdXi2<6,6>() },
                                          { dNdXi0<6,7>(), dNdXi1<6,7>(), dNdXi2<6,7>() }
                                        },
                                        { { dNdXi0<7,0>(), dNdXi1<7,0>(), dNdXi2<7,0>() },
                                          { dNdXi0<7,1>(), dNdXi1<7,1>(), dNdXi2<7,1>() },
                                          { dNdXi0<7,2>(), dNdXi1<7,2>(), dNdXi2<7,2>() },
                                          { dNdXi0<7,3>(), dNdXi1<7,3>(), dNdXi2<7,3>() },
                                          { dNdXi0<7,4>(), dNdXi1<7,4>(), dNdXi2<7,4>() },
                                          { dNdXi0<7,5>(), dNdXi1<7,5>(), dNdXi2<7,5>() },
                                          { dNdXi0<7,6>(), dNdXi1<7,6>(), dNdXi2<7,6>() },
                                          { dNdXi0<7,7>(), dNdXi1<7,7>(), dNdXi2<7,7>() }
                                        } };


    for( localIndex a=0 ; a<numNodes ; ++a )
    {
//      parentShapeFunctionDerivatives( a, xi0, xi1, xi2, dNdXi );
      
      #pragma unroll
      for ( int i = 0; i < 3; ++i )
      {
        #pragma unroll
        for ( int j = 0; j < 3; ++j )
        {
          J[i][j] += X[elemsToNodes(k, a)][i] * dNdXi[q][a][j];
        }
      }
    }

    real64 const detJ = inverse( J );

    #pragma unroll
    for( localIndex a=0 ; a<numNodes ; ++a )
    {
//      parentShapeFunctionDerivatives( a, xi0, xi1, xi2, dNdXi );

      #pragma unroll
      for ( int i = 0; i < 3; ++i )
      {
//        dNdX[i][a] = J[0][i] * dNdXi[0] + J[1][i] * dNdXi[1] + J[2][i] * dNdXi[2];
        dNdX[i][a] = J[0][i] * dNdXi[q][a][0] + J[1][i] * dNdXi[q][a][1] + J[2][i] * dNdXi[q][a][2];
      }
    }

    return detJ;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivatives( localIndex const k,
                                          localIndex const q,
                                          arrayView2d<localIndex const> const & elemsToNodes,
                                          arrayView1d<R1Tensor const> const & X,
                                          real64 (&dNdXi)[3][numNodes] )
  {
    return shapeFunctionDerivativesImpl( k,
                                     q,
                                     elemsToNodes,
//                                     quadratureFactor*parentCoords(0,q),
//                                     quadratureFactor*parentCoords(1,q),
//                                     quadratureFactor*parentCoords(2,q),
                                     X,
                                     dNdXi );
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
