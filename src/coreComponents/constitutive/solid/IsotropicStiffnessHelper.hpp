/*
 * BTransposeDB_Helper.hpp
 *
 *  Created on: Aug 11, 2020
 *      Author: settgast
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_ISOTROPICSTIFFNESSHELPER_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_ISOTROPICSTIFFNESSHELPER_HPP_

#include "StiffnessHelperBase.hpp"

namespace geosx
{
namespace constitutive
{


struct IsotropicStiffnessHelper : public StiffnessHelperBase
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
    m_lambda = C[1][0];
    m_shearModulus = C[5][5];
  }

  real64 m_lambda;
  real64 m_shearModulus;
};


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT >
void IsotropicStiffnessHelper::BTDB( BASIS_GRADIENT const & dNdX,
                                     real64 const & detJ,
                                     real64 (& localJacobian)[NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJ;
  StiffnessHelperBase::BTDB< NUM_SUPPORT_POINTS,
                             NUM_DOF_PER_SUPPORT_POINT >( dNdX,
                                                          localJacobian,
                                                          [ lambda = this->m_lambda * detJ,
                                                            G = this->m_shearModulus * detJ,
                                                            lambda2G ]
                                                            ( int const a,
                                                            int const b,
                                                            real64 const (&dNdXadNdXb)[3][3],
                                                            real64 (& localJacobian)[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT] )
  {
    localJacobian[a*3+0][b*3+0] = localJacobian[a*3+0][b*3+0] - dNdXadNdXb[1][1] * G - dNdXadNdXb[2][2] * G - dNdXadNdXb[0][0] * lambda2G;
    localJacobian[a*3+0][b*3+1] = localJacobian[a*3+0][b*3+1] - dNdXadNdXb[1][0] * G - dNdXadNdXb[0][1] * lambda;
    localJacobian[a*3+0][b*3+2] = localJacobian[a*3+0][b*3+2] - dNdXadNdXb[2][0] * G - dNdXadNdXb[0][2] * lambda;
    localJacobian[a*3+1][b*3+0] = localJacobian[a*3+1][b*3+0] - dNdXadNdXb[0][1] * G - dNdXadNdXb[1][0] * lambda;
    localJacobian[a*3+1][b*3+1] = localJacobian[a*3+1][b*3+1] - dNdXadNdXb[0][0] * G - dNdXadNdXb[2][2] * G - dNdXadNdXb[1][1] * lambda2G;
    localJacobian[a*3+1][b*3+2] = localJacobian[a*3+1][b*3+2] - dNdXadNdXb[2][1] * G - dNdXadNdXb[1][2] * lambda;
    localJacobian[a*3+2][b*3+0] = localJacobian[a*3+2][b*3+0] - dNdXadNdXb[0][2] * G - dNdXadNdXb[2][0] * lambda;
    localJacobian[a*3+2][b*3+1] = localJacobian[a*3+2][b*3+1] - dNdXadNdXb[1][2] * G - dNdXadNdXb[2][1] * lambda;
    localJacobian[a*3+2][b*3+2] = localJacobian[a*3+2][b*3+2] - dNdXadNdXb[0][0] * G - dNdXadNdXb[1][1] * G - dNdXadNdXb[2][2] * lambda2G;
  } );
}


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT >
void IsotropicStiffnessHelper::diagBTDB( BASIS_GRADIENT const & dNdX,
                                         real64 const & detJ,
                                         real64 (& diagLocalJacobian)[NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJ;

  StiffnessHelperBase::diagBTDB< NUM_SUPPORT_POINTS,
                                 NUM_DOF_PER_SUPPORT_POINT >( dNdX,
                                                              diagLocalJacobian,
                                                              [ lambda = m_lambda * detJ,
                                                                G = m_shearModulus * detJ,
                                                                lambda2G ]
                                                                ( const int a,
                                                                const real64 (& dNdXdNdX)[3],
                                                                real64 (& diagLocalJacobian)[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT] )
  {
    diagLocalJacobian[ a*3+0 ] = diagLocalJacobian[ a*3+0 ] - dNdXdNdX[0] * lambda2G - dNdXdNdX[1] * G - dNdXdNdX[2] * G;
    diagLocalJacobian[ a*3+1 ] = diagLocalJacobian[ a*3+1 ] - dNdXdNdX[1] * lambda2G - dNdXdNdX[0] * G - dNdXdNdX[2] * G;
    diagLocalJacobian[ a*3+2 ] = diagLocalJacobian[ a*3+2 ] - dNdXdNdX[2] * lambda2G - dNdXdNdX[0] * G - dNdXdNdX[1] * G;
  } );
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_ISOTROPICSTIFFNESSHELPER_HPP_ */
