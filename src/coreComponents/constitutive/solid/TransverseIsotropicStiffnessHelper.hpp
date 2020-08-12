/*
 * BTransposeDB_Helper.hpp
 *
 *  Created on: Aug 11, 2020
 *      Author: settgast
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_TRANSVERSEISOTROPICSTIFFNESSHELPER_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_TRANSVERSEISOTROPICSTIFFNESSHELPER_HPP_

#include "StiffnessHelperBase.hpp"

namespace geosx
{
namespace constitutive
{


struct TransverseIsotropicStiffnessHelper : public StiffnessHelperBase
{
  template< int NUM_SUPPORT_POINTS,
            int NUM_DOF_PER_SUPPORT_POINT,
            typename BASIS_GRADIENT >
  void BTDB( BASIS_GRADIENT const & dNdX,
             real64 const & detJ,
             real64 ( &localJacobian )[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT] );

  template< int NUM_SUPPORT_POINTS,
            int NUM_DOF_PER_SUPPORT_POINT,
            typename BASIS_GRADIENT >
  void diagBTDB( BASIS_GRADIENT const & dNdX,
                 real64 const & detJ,
                 real64 ( &diagLocalJacobian )[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT] );

  void setParams( real64 (& C)[6][6] )
  {
    m_c11 = C[0][0];
    m_c13 = C[0][2];
    m_c33 = C[2][2];
    m_c44 = C[3][3];
    m_c66 = C[5][5];
  }

  real64 m_c11;
  real64 m_c13;
  real64 m_c33;
  real64 m_c44;
  real64 m_c66;
};


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT >
void TransverseIsotropicStiffnessHelper::BTDB( BASIS_GRADIENT const & dNdX,
                                               real64 const & detJ,
                                               real64 (& localJacobian)[NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT] )
{
  real64 const c12 = (m_c11 - 2 * m_c66) * detJ;
  StiffnessHelperBase::BTDB< NUM_SUPPORT_POINTS,
                             NUM_DOF_PER_SUPPORT_POINT >( dNdX,
                                                          localJacobian,
                                                          [ c11 = this->m_c11 * detJ,
                                                            c13 = this->m_c13 * detJ,
                                                            c33 = this->m_c33 * detJ,
                                                            c44 = this->m_c44 * detJ,
                                                            c66 = this->m_c66 * detJ,
                                                            c12 ]
                                                            ( int const a,
                                                            int const b,
                                                            real64 const (&dNdXadNdXb)[3][3],
                                                            real64 (& localJacobian)[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT] )
  {
    localJacobian[a*3+0][b*3+0] = localJacobian[a*3+0][b*3+0] - c11 * dNdXadNdXb[0][0] - c66 * dNdXadNdXb[1][1] - c44 * dNdXadNdXb[2][2];
    localJacobian[a*3+0][b*3+1] = localJacobian[a*3+0][b*3+1] - c12 * dNdXadNdXb[0][1] - c66 * dNdXadNdXb[1][0];
    localJacobian[a*3+0][b*3+2] = localJacobian[a*3+0][b*3+2] - c13 * dNdXadNdXb[0][2] - c44 * dNdXadNdXb[2][0];
    localJacobian[a*3+1][b*3+0] = localJacobian[a*3+1][b*3+0] - c66 * dNdXadNdXb[0][1] - c12 * dNdXadNdXb[1][0];
    localJacobian[a*3+1][b*3+1] = localJacobian[a*3+1][b*3+1] - c66 * dNdXadNdXb[0][0] - c11 * dNdXadNdXb[1][1] - c44 * dNdXadNdXb[2][2];
    localJacobian[a*3+1][b*3+2] = localJacobian[a*3+1][b*3+2] - c13 * dNdXadNdXb[1][2] - c44 * dNdXadNdXb[2][1];
    localJacobian[a*3+2][b*3+0] = localJacobian[a*3+2][b*3+0] - c44 * dNdXadNdXb[0][2] - c13 * dNdXadNdXb[2][0];
    localJacobian[a*3+2][b*3+1] = localJacobian[a*3+2][b*3+1] - c44 * dNdXadNdXb[1][2] - c13 * dNdXadNdXb[2][1];
    localJacobian[a*3+2][b*3+2] = localJacobian[a*3+2][b*3+2] - c44 * dNdXadNdXb[0][0] - c44 * dNdXadNdXb[1][1] - c33 * dNdXadNdXb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT >
void TransverseIsotropicStiffnessHelper::diagBTDB( BASIS_GRADIENT const & dNdX,
                                                   real64 const & detJ,
                                                   real64 (& diagLocalJacobian)[NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT] )
{
  StiffnessHelperBase::diagBTDB< NUM_SUPPORT_POINTS,
                                 NUM_DOF_PER_SUPPORT_POINT >( dNdX,
                                                              diagLocalJacobian,
                                                              [ c11 = this->m_c11 * detJ,
                                                                c33 = this->m_c33 * detJ,
                                                                c44 = this->m_c44 * detJ,
                                                                c66 = this->m_c66 * detJ ]
                                                                ( const int a,
                                                                const real64 (& dNdXdNdX)[3],
                                                                real64 (& diagLocalJacobian)[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT] )
  {
    diagLocalJacobian[a*3+0] = diagLocalJacobian[a*3+0] - c11 * dNdXdNdX[0] - c66 * dNdXdNdX[1] - c44 * dNdXdNdX[2];
    diagLocalJacobian[a*3+1] = diagLocalJacobian[a*3+1] - c66 * dNdXdNdX[0] - c11 * dNdXdNdX[1] - c44 * dNdXdNdX[2];
    diagLocalJacobian[a*3+2] = diagLocalJacobian[a*3+2] - c44 * dNdXdNdX[0] - c44 * dNdXdNdX[1] - c33 * dNdXdNdX[2];
  } );
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_TRANSVERSEISOTROPICSTIFFNESSHELPER_HPP_ */
