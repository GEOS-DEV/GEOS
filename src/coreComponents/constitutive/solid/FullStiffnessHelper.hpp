/*
 * BTransposeDB_Helper.hpp
 *
 *  Created on: Aug 11, 2020
 *      Author: settgast
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_FULLSTIFFNESSHELPER_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_FULLSTIFFNESSHELPER_HPP_

#include "StiffnessHelperBase.hpp"

namespace geosx
{
namespace constitutive
{


struct FullStiffnessHelper : public StiffnessHelperBase
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
    for( int i=0; i<6; ++i )
    {
      for( int j=0; j<6; ++j )
      {
        m_c[i][j] = C[i][j];
      }
    }
  }

  real64 m_c[6][6];
};


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT >
void FullStiffnessHelper::BTDB( BASIS_GRADIENT const & dNdX,
                                real64 const & detJ,
                                real64 (& localJacobian)[NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT] )
{
  GEOSX_UNUSED_VAR( detJ );
  StiffnessHelperBase::BTDB< NUM_SUPPORT_POINTS,
                             NUM_DOF_PER_SUPPORT_POINT >( dNdX,
                                                          localJacobian,
                                                          [ c = this->m_c ]
                                                            ( int const a,
                                                            int const b,
                                                            real64 const (&dNdXadNdXb)[3][3],
                                                            real64 (& localJacobian)[NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT][NUM_SUPPORT_POINTS*NUM_DOF_PER_SUPPORT_POINT] )
  {
    localJacobian[a*3+0][b*3+0] = localJacobian[a*3+0][b*3+0] -
                                  c[0][0] * dNdXadNdXb[0][0] - c[0][5] * dNdXadNdXb[0][1] - c[0][4] * dNdXadNdXb[0][2] -
                                  c[5][0] * dNdXadNdXb[1][0] - c[5][5] * dNdXadNdXb[1][1] - c[5][4] * dNdXadNdXb[1][2] -
                                  c[4][0] * dNdXadNdXb[2][0] - c[4][5] * dNdXadNdXb[2][1] - c[4][4] * dNdXadNdXb[2][2];
    localJacobian[a*3+0][b*3+1] = localJacobian[a*3+0][b*3+1] -
                                  c[0][5] * dNdXadNdXb[0][0] - c[0][1] * dNdXadNdXb[0][1] - c[0][3] * dNdXadNdXb[0][2] -
                                  c[5][5] * dNdXadNdXb[1][0] - c[5][1] * dNdXadNdXb[1][1] - c[5][3] * dNdXadNdXb[1][2] -
                                  c[4][5] * dNdXadNdXb[2][0] - c[4][1] * dNdXadNdXb[2][1] - c[4][3] * dNdXadNdXb[2][2];
    localJacobian[a*3+0][b*3+2] = localJacobian[a*3+0][b*3+2] -
                                  c[0][4] * dNdXadNdXb[0][0] - c[0][3] * dNdXadNdXb[0][1] - c[0][2] * dNdXadNdXb[0][2] -
                                  c[5][4] * dNdXadNdXb[1][0] - c[5][3] * dNdXadNdXb[1][1] - c[5][2] * dNdXadNdXb[1][2] -
                                  c[4][4] * dNdXadNdXb[2][0] - c[4][3] * dNdXadNdXb[2][1] - c[4][2] * dNdXadNdXb[2][2];
    localJacobian[a*3+1][b*3+0] = localJacobian[a*3+1][b*3+0] -
                                  c[5][0] * dNdXadNdXb[0][0] - c[5][5] * dNdXadNdXb[0][1] - c[5][4] * dNdXadNdXb[0][2] -
                                  c[1][0] * dNdXadNdXb[1][0] - c[1][5] * dNdXadNdXb[1][1] - c[1][4] * dNdXadNdXb[1][2] -
                                  c[3][0] * dNdXadNdXb[2][0] - c[3][5] * dNdXadNdXb[2][1] - c[3][4] * dNdXadNdXb[2][2];
    localJacobian[a*3+1][b*3+1] = localJacobian[a*3+1][b*3+1] -
                                  c[5][5] * dNdXadNdXb[0][0] - c[5][1] * dNdXadNdXb[0][1] - c[5][3] * dNdXadNdXb[0][2] -
                                  c[1][5] * dNdXadNdXb[1][0] - c[1][1] * dNdXadNdXb[1][1] - c[1][3] * dNdXadNdXb[1][2] -
                                  c[3][5] * dNdXadNdXb[2][0] - c[3][1] * dNdXadNdXb[2][1] - c[3][3] * dNdXadNdXb[2][2];
    localJacobian[a*3+1][b*3+2] = localJacobian[a*3+1][b*3+2] -
                                  c[5][4] * dNdXadNdXb[0][0] - c[5][3] * dNdXadNdXb[0][1] - c[5][2] * dNdXadNdXb[0][2] -
                                  c[1][4] * dNdXadNdXb[1][0] - c[1][3] * dNdXadNdXb[1][1] - c[1][2] * dNdXadNdXb[1][2] -
                                  c[3][4] * dNdXadNdXb[2][0] - c[3][3] * dNdXadNdXb[2][1] - c[3][2] * dNdXadNdXb[2][2];
    localJacobian[a*3+2][b*3+0] = localJacobian[a*3+2][b*3+0] -
                                  c[4][0] * dNdXadNdXb[0][0] - c[4][5] * dNdXadNdXb[0][1] - c[4][4] * dNdXadNdXb[0][2] -
                                  c[3][0] * dNdXadNdXb[1][0] - c[3][5] * dNdXadNdXb[1][1] - c[3][4] * dNdXadNdXb[1][2] -
                                  c[2][0] * dNdXadNdXb[2][0] - c[2][5] * dNdXadNdXb[2][1] - c[2][4] * dNdXadNdXb[2][2];
    localJacobian[a*3+2][b*3+1] = localJacobian[a*3+2][b*3+1] -
                                  c[4][5] * dNdXadNdXb[0][0] - c[4][1] * dNdXadNdXb[0][1] - c[4][3] * dNdXadNdXb[0][2] -
                                  c[3][5] * dNdXadNdXb[1][0] - c[3][1] * dNdXadNdXb[1][1] - c[3][3] * dNdXadNdXb[1][2] -
                                  c[2][5] * dNdXadNdXb[2][0] - c[2][1] * dNdXadNdXb[2][1] - c[2][3] * dNdXadNdXb[2][2];
    localJacobian[a*3+2][b*3+2] = localJacobian[a*3+2][b*3+2] -
                                  c[4][4] * dNdXadNdXb[0][0] - c[4][3] * dNdXadNdXb[0][1] - c[4][2] * dNdXadNdXb[0][2] -
                                  c[3][4] * dNdXadNdXb[1][0] - c[3][3] * dNdXadNdXb[1][1] - c[3][2] * dNdXadNdXb[1][2] -
                                  c[2][4] * dNdXadNdXb[2][0] - c[2][3] * dNdXadNdXb[2][1] - c[2][2] * dNdXadNdXb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          int NUM_DOF_PER_SUPPORT_POINT,
          typename BASIS_GRADIENT >
void FullStiffnessHelper::diagBTDB( BASIS_GRADIENT const & dNdX,
                                    real64 const & detJ,
                                    real64 (& diagLocalJacobian)[NUM_SUPPORT_POINTS *NUM_DOF_PER_SUPPORT_POINT] )
{
  GEOSX_UNUSED_VAR( detJ );

  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    real64 const dNdXadNdXa[3][3] = { { dNdX[a][0] * dNdX[a][0], dNdX[a][0] * dNdX[a][1], dNdX[a][0] * dNdX[a][2] },
      { dNdX[a][1] * dNdX[a][0], dNdX[a][1] * dNdX[a][1], dNdX[a][1] * dNdX[a][2] },
      { dNdX[a][2] * dNdX[a][0], dNdX[a][2] * dNdX[a][1], dNdX[a][2] * dNdX[a][2] } };
    diagLocalJacobian[a*3+0] = diagLocalJacobian[a*3+0] -
                               m_c[0][0] * dNdXadNdXa[0][0] - m_c[0][5] * dNdXadNdXa[0][1] - m_c[0][4] * dNdXadNdXa[0][2] -
                               m_c[5][0] * dNdXadNdXa[1][0] - m_c[5][5] * dNdXadNdXa[1][1] - m_c[5][4] * dNdXadNdXa[1][2] -
                               m_c[4][0] * dNdXadNdXa[2][0] - m_c[4][5] * dNdXadNdXa[2][1] - m_c[4][4] * dNdXadNdXa[2][2];
    diagLocalJacobian[a*3+1] = diagLocalJacobian[a*3+1] -
                               m_c[5][5] * dNdXadNdXa[0][0] - m_c[5][1] * dNdXadNdXa[0][1] - m_c[5][3] * dNdXadNdXa[0][2] -
                               m_c[1][5] * dNdXadNdXa[1][0] - m_c[1][1] * dNdXadNdXa[1][1] - m_c[1][3] * dNdXadNdXa[1][2] -
                               m_c[3][5] * dNdXadNdXa[2][0] - m_c[3][1] * dNdXadNdXa[2][1] - m_c[3][3] * dNdXadNdXa[2][2];
    diagLocalJacobian[a*3+2] = diagLocalJacobian[a*3+2] -
                               m_c[4][4] * dNdXadNdXa[0][0] - m_c[4][3] * dNdXadNdXa[0][1] - m_c[4][2] * dNdXadNdXa[0][2] -
                               m_c[3][4] * dNdXadNdXa[1][0] - m_c[3][3] * dNdXadNdXa[1][1] - m_c[3][2] * dNdXadNdXa[1][2] -
                               m_c[2][4] * dNdXadNdXa[2][0] - m_c[2][3] * dNdXadNdXa[2][1] - m_c[2][2] * dNdXadNdXa[2][2];
  }
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_FULLSTIFFNESSHELPER_HPP_ */
