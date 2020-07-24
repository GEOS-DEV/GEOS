

#ifndef GEOSX_CORE_FINITEELEMENT_FINITEELEMENTSHAPEFUNCTIONKERNELBASE
#define GEOSX_CORE_FINITEELEMENT_FINITEELEMENTSHAPEFUNCTIONKERNELBASE

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"


namespace geosx
{

/**
 * @class FiniteElementShapeFunctionKernelBase
 * @brief Base class for the finite element kernels.
 */
class FiniteElementShapeFunctionKernelBase
{
public:

  virtual localIndex getNumQuadraturePoints() const = 0;
  virtual localIndex getNumSupportPoints() const = 0;

  virtual ~FiniteElementShapeFunctionKernelBase()
  {}

  /**
   * @brief Computes the inverse of a 3x3 c-array.
   * @param J The array to invert...which is also used to store the inverse.
   * @return The determinant of @p J.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 inverse( real64 (& J)[3][3] )
  {
    real64 const temp[3][3] =
    { { J[1][1]*J[2][2] - J[1][2]*J[2][1], J[0][2]*J[2][1] - J[0][1]*J[2][2], J[0][1]*J[1][2] - J[0][2]*J[1][1] },
      { J[1][2]*J[2][0] - J[1][0]*J[2][2], J[0][0]*J[2][2] - J[0][2]*J[2][0], J[0][2]*J[1][0] - J[0][0]*J[1][2] },
      { J[1][0]*J[2][1] - J[1][1]*J[2][0], J[0][1]*J[2][0] - J[0][0]*J[2][1], J[0][0]*J[1][1] - J[0][1]*J[1][0] } };

    real64 const det =  J[0][0] * temp[0][0] + J[1][0] * temp[0][1] + J[2][0] * temp[0][2];
    real64 const invDet = 1.0 / det;

    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        J[i][j] = temp[i][j] * invDet;
      }
    }
    return det;
  }

  /**
   * @brief Calculate the determinant of a 3x3 c-array.
   * @param J The input array.
   * @return The determinant of @p J
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 detJ( real64 const (&J)[3][3] )
  {
    return J[0][0] * ( J[1][1]*J[2][2] - J[1][2]*J[2][1] ) +
           J[1][0] * ( J[0][2]*J[2][1] - J[0][1]*J[2][2] ) +
           J[2][0] * ( J[0][0]*J[1][1] - J[0][1]*J[1][0] );
  }


//TODO we want to keep views and provide interfaces to this data here for cases
//     where we pre-compute the shape function derivatives...maybe...tbd.
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  real64 shapeFunctionDerivatives( localIndex k, localIndex q, real const (&X)[numNodes][3] , real64 (& dNdX)[numNodes][3] )
//  {
//    for( int a=0 ; a<numNodes; ++a )
//    {
//      for( int i=0; i<3; ++a )
//      {
//        dNdX[a][i] = m_dNdX( k, q, a, i);
//      }
//    }
//    return m_detJ( k, q );
//  }
//private:
//  arrayView4d< real64 const > const m_dNdX;
//  arrayView2d< real64 const > const m_detJ;


};

}

#endif //GEOSX_CORE_FINITEELEMENT_FINITEELEMENTSHAPEFUNCTIONKERNELBASE
