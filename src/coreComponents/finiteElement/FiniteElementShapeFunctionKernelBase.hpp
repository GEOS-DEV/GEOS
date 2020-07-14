

#ifndef GEOSX_CORE_FINITEELEMENT_FINITEELEMENTSHAPEFUNCTIONKERNELBASE
#define GEOSX_CORE_FINITEELEMENT_FINITEELEMENTSHAPEFUNCTIONKERNELBASE

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"


namespace geosx
{

template< int NUM_SUPPORT_POINTS, int NUM_QUADRATURE_POINTS >
class FiniteElementShapeFunctionKernelBase
{
public:
  constexpr static localIndex numNodes = NUM_SUPPORT_POINTS;
  constexpr static localIndex numQuadraturePoints = NUM_QUADRATURE_POINTS;


//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  static void shapeFunctionValues( localIndex const q,
//                                   real64 N[numNodes] ) {}
//
//
//
//  GEOSX_HOST_DEVICE
//  GEOSX_FORCE_INLINE
//  static real64 shapeFunctionDerivatives( localIndex const q,
//                                          real64 const (&X)[numNodes][3],
//                                          real64 (& dNdX)[numNodes][3] ) {}


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 inverse( real64 (& J)[3][3] )
  {
    real64 const temp[3][3] = { { J[1][1]*J[2][2] - J[1][2]*J[2][1], J[0][2]*J[2][1] - J[0][1]*J[2][2], J[0][1]*J[1][2] - J[0][2]*J[1][1] },
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

};

}

#endif //GEOSX_CORE_FINITEELEMENT_FINITEELEMENTSHAPEFUNCTIONKERNELBASE
