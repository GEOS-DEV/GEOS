

#ifndef FINITE_ELEMENT_SHAPE_KERNEL
#define FINITE_ELEMENT_SHAPE_KERNEL

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "physicsSolvers/solidMechanics/KernelMacros.hpp"

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
  static real64 parentCoords( localIndex const i, localIndex const a )
  {
    constexpr static  real64 pCoords[3][8] = { { -1,  1, -1,  1, -1,  1, -1,  1 },
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



  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivatives( localIndex const k,
                                          localIndex const q,
                                          arrayView2d<localIndex const> const & elemsToNodes,
                                          real64 const xi0,
                                          real64 const xi1,
                                          real64 const xi2,
                                        #if STORE_NODE_DATA_LOCALLY
                                          real64 const (&x_local)[3][numNodes],
                                        #else
                                          arrayView1d<R1Tensor const> const & X,
                                        #endif
                                          real64 (&dNdX)[3][numNodes] )
  {
    real64 J[3][3] = {{0}};
    real64 dNdXi[3];

    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      parentShapeFunctionDerivatives( a, xi0, xi1, xi2, dNdXi );
      
      #pragma unroll
      for ( int i = 0; i < 3; ++i )
      {
        real64 const pos_k_a_i = POSITION_ACCESSOR(k, a, i);

        #pragma unroll
        for ( int j = 0; j < 3; ++j )
        {
          J[i][j] += pos_k_a_i * dNdXi[j];
        }
      }
    }

    real64 const detJ = inverse( J );

    #pragma unroll
    for( localIndex a=0 ; a<numNodes ; ++a )
    {
      parentShapeFunctionDerivatives( a, xi0, xi1, xi2, dNdXi );

      #pragma unroll
      for ( int i = 0; i < 3; ++i )
      {
        dNdX[i][a] = J[0][i] * dNdXi[0] + J[1][i] * dNdXi[1] + J[2][i] * dNdXi[2];
      }
    }

    return detJ;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 shapeFunctionDerivatives( localIndex const k,
                                          localIndex const q,
                                          arrayView2d<localIndex const> const & elemsToNodes,
                                        #if STORE_NODE_DATA_LOCALLY
                                          real64 const (&X)[3][numNodes],
                                        #else
                                          arrayView1d<R1Tensor const> const & X,
                                        #endif
                                          real64 (&dNdXi)[3][numNodes] )
  {
    return shapeFunctionDerivatives( k,
                                     q,
                                     elemsToNodes,
                                     quadratureFactor*parentCoords(0,q),
                                     quadratureFactor*parentCoords(1,q),
                                     quadratureFactor*parentCoords(2,q),
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
    real64 const o10 = 1.0 / detJ;

    J[0][0] = o1 * o10;
    J[0][1] = o2 * o10;
    J[0][2] = o3 * o10;
    J[1][0] = o4 * o10;
    J[1][1] = o5 * o10;
    J[1][2] = o6 * o10;
    J[2][0] = o7 * o10;
    J[2][1] = o8 * o10;
    J[2][2] = o9 * o10;

    return detJ;
  }

};

}

#endif //FINITE_ELEMENT_SHAPE_KERNEL
