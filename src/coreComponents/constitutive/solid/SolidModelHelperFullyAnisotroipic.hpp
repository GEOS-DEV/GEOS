/*
 * BTransposeDB_Helper.hpp
 *
 *  Created on: Aug 11, 2020
 *      Author: settgast
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERANISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERANISOTROPIC_HPP_

#include "SolidModelHelperBase.hpp"

namespace geosx
{
namespace constitutive
{


struct SolidModelHelperFullyAnisotroipic : public SolidModelHelperBase
{
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void BTDB( BASIS_GRADIENT const & gradN,
             real64 const & detJxW,
             real64 ( &localStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 const & detJxW,
                 real64 ( &diagLocalStiffness )[NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 const & detJxW,
                       real64 ( &diagSumLocalStiffness )[NUM_SUPPORT_POINTS*3] );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
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
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperFullyAnisotroipic::BTDB( BASIS_GRADIENT const & gradN,
                                              real64 const & detJxW,
                                              real64 (& localStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };


      localStiffness[a*3+0][b*3+0] = localStiffness[a*3+0][b*3+0] -
                                     ( m_c[0][0] * gradNa_gradNb[0][0] - m_c[0][5] * gradNa_gradNb[0][1] - m_c[0][4] * gradNa_gradNb[0][2] -
                                       m_c[5][0] * gradNa_gradNb[1][0] - m_c[5][5] * gradNa_gradNb[1][1] - m_c[5][4] * gradNa_gradNb[1][2] -
                                       m_c[4][0] * gradNa_gradNb[2][0] - m_c[4][5] * gradNa_gradNb[2][1] - m_c[4][4] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+0][b*3+1] = localStiffness[a*3+0][b*3+1] -
                                     (m_c[0][5] * gradNa_gradNb[0][0] - m_c[0][1] * gradNa_gradNb[0][1] - m_c[0][3] * gradNa_gradNb[0][2] -
                                      m_c[5][5] * gradNa_gradNb[1][0] - m_c[5][1] * gradNa_gradNb[1][1] - m_c[5][3] * gradNa_gradNb[1][2] -
                                      m_c[4][5] * gradNa_gradNb[2][0] - m_c[4][1] * gradNa_gradNb[2][1] - m_c[4][3] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+0][b*3+2] = localStiffness[a*3+0][b*3+2] -
                                     (m_c[0][4] * gradNa_gradNb[0][0] - m_c[0][3] * gradNa_gradNb[0][1] - m_c[0][2] * gradNa_gradNb[0][2] -
                                      m_c[5][4] * gradNa_gradNb[1][0] - m_c[5][3] * gradNa_gradNb[1][1] - m_c[5][2] * gradNa_gradNb[1][2] -
                                      m_c[4][4] * gradNa_gradNb[2][0] - m_c[4][3] * gradNa_gradNb[2][1] - m_c[4][2] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+1][b*3+0] = localStiffness[a*3+1][b*3+0] -
                                     (m_c[5][0] * gradNa_gradNb[0][0] - m_c[5][5] * gradNa_gradNb[0][1] - m_c[5][4] * gradNa_gradNb[0][2] -
                                      m_c[1][0] * gradNa_gradNb[1][0] - m_c[1][5] * gradNa_gradNb[1][1] - m_c[1][4] * gradNa_gradNb[1][2] -
                                      m_c[3][0] * gradNa_gradNb[2][0] - m_c[3][5] * gradNa_gradNb[2][1] - m_c[3][4] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+1][b*3+1] = localStiffness[a*3+1][b*3+1] -
                                     (m_c[5][5] * gradNa_gradNb[0][0] - m_c[5][1] * gradNa_gradNb[0][1] - m_c[5][3] * gradNa_gradNb[0][2] -
                                      m_c[1][5] * gradNa_gradNb[1][0] - m_c[1][1] * gradNa_gradNb[1][1] - m_c[1][3] * gradNa_gradNb[1][2] -
                                      m_c[3][5] * gradNa_gradNb[2][0] - m_c[3][1] * gradNa_gradNb[2][1] - m_c[3][3] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+1][b*3+2] = localStiffness[a*3+1][b*3+2] -
                                     (m_c[5][4] * gradNa_gradNb[0][0] - m_c[5][3] * gradNa_gradNb[0][1] - m_c[5][2] * gradNa_gradNb[0][2] -
                                      m_c[1][4] * gradNa_gradNb[1][0] - m_c[1][3] * gradNa_gradNb[1][1] - m_c[1][2] * gradNa_gradNb[1][2] -
                                      m_c[3][4] * gradNa_gradNb[2][0] - m_c[3][3] * gradNa_gradNb[2][1] - m_c[3][2] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+2][b*3+0] = localStiffness[a*3+2][b*3+0] -
                                     (m_c[4][0] * gradNa_gradNb[0][0] - m_c[4][5] * gradNa_gradNb[0][1] - m_c[4][4] * gradNa_gradNb[0][2] -
                                      m_c[3][0] * gradNa_gradNb[1][0] - m_c[3][5] * gradNa_gradNb[1][1] - m_c[3][4] * gradNa_gradNb[1][2] -
                                      m_c[2][0] * gradNa_gradNb[2][0] - m_c[2][5] * gradNa_gradNb[2][1] - m_c[2][4] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+2][b*3+1] = localStiffness[a*3+2][b*3+1] -
                                     (m_c[4][5] * gradNa_gradNb[0][0] - m_c[4][1] * gradNa_gradNb[0][1] - m_c[4][3] * gradNa_gradNb[0][2] -
                                      m_c[3][5] * gradNa_gradNb[1][0] - m_c[3][1] * gradNa_gradNb[1][1] - m_c[3][3] * gradNa_gradNb[1][2] -
                                      m_c[2][5] * gradNa_gradNb[2][0] - m_c[2][1] * gradNa_gradNb[2][1] - m_c[2][3] * gradNa_gradNb[2][2] ) * detJxW;
      localStiffness[a*3+2][b*3+2] = localStiffness[a*3+2][b*3+2] -
                                     (m_c[4][4] * gradNa_gradNb[0][0] - m_c[4][3] * gradNa_gradNb[0][1] - m_c[4][2] * gradNa_gradNb[0][2] -
                                      m_c[3][4] * gradNa_gradNb[1][0] - m_c[3][3] * gradNa_gradNb[1][1] - m_c[3][2] * gradNa_gradNb[1][2] -
                                      m_c[2][4] * gradNa_gradNb[2][0] - m_c[2][3] * gradNa_gradNb[2][1] - m_c[2][2] * gradNa_gradNb[2][2] ) * detJxW;
    }
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperFullyAnisotroipic::diagBTDB( BASIS_GRADIENT const & gradN,
                                                  real64 const & detJxW,
                                                  real64 (& diagLocalStiffness)[NUM_SUPPORT_POINTS *3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    real64 const gradN_gradN[3][3] =
    { { gradN[a][0] * gradN[a][0], gradN[a][0] * gradN[a][1], gradN[a][0] * gradN[a][2] },
      { gradN[a][1] * gradN[a][0], gradN[a][1] * gradN[a][1], gradN[a][1] * gradN[a][2] },
      { gradN[a][2] * gradN[a][0], gradN[a][2] * gradN[a][1], gradN[a][2] * gradN[a][2] } };
    diagLocalStiffness[a*3+0] = diagLocalStiffness[a*3+0] -
                                (m_c[0][0] * gradN_gradN[0][0] - m_c[0][5] * gradN_gradN[0][1] - m_c[0][4] * gradN_gradN[0][2] -
                                 m_c[5][0] * gradN_gradN[1][0] - m_c[5][5] * gradN_gradN[1][1] - m_c[5][4] * gradN_gradN[1][2] -
                                 m_c[4][0] * gradN_gradN[2][0] - m_c[4][5] * gradN_gradN[2][1] - m_c[4][4] * gradN_gradN[2][2] ) * detJxW;
    diagLocalStiffness[a*3+1] = diagLocalStiffness[a*3+1] -
                                (m_c[5][5] * gradN_gradN[0][0] - m_c[5][1] * gradN_gradN[0][1] - m_c[5][3] * gradN_gradN[0][2] -
                                 m_c[1][5] * gradN_gradN[1][0] - m_c[1][1] * gradN_gradN[1][1] - m_c[1][3] * gradN_gradN[1][2] -
                                 m_c[3][5] * gradN_gradN[2][0] - m_c[3][1] * gradN_gradN[2][1] - m_c[3][3] * gradN_gradN[2][2] ) * detJxW;
    diagLocalStiffness[a*3+2] = diagLocalStiffness[a*3+2] -
                                (m_c[4][4] * gradN_gradN[0][0] - m_c[4][3] * gradN_gradN[0][1] - m_c[4][2] * gradN_gradN[0][2] -
                                 m_c[3][4] * gradN_gradN[1][0] - m_c[3][3] * gradN_gradN[1][1] - m_c[3][2] * gradN_gradN[1][2] -
                                 m_c[2][4] * gradN_gradN[2][0] - m_c[2][3] * gradN_gradN[2][1] - m_c[2][2] * gradN_gradN[2][2] ) * detJxW;
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperFullyAnisotroipic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                        real64 const & detJxW,
                                                        real64 ( & diagSumLocalStiffness )[NUM_SUPPORT_POINTS*3] )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };


      diagSumLocalStiffness[a*3+0] = diagSumLocalStiffness[a*3+0] -
                                     ( m_c[0][0] * gradNa_gradNb[0][0] - m_c[0][5] * gradNa_gradNb[0][1] - m_c[0][4] * gradNa_gradNb[0][2] -
                                       m_c[5][0] * gradNa_gradNb[1][0] - m_c[5][5] * gradNa_gradNb[1][1] - m_c[5][4] * gradNa_gradNb[1][2] -
                                       m_c[4][0] * gradNa_gradNb[2][0] - m_c[4][5] * gradNa_gradNb[2][1] - m_c[4][4] * gradNa_gradNb[2][2] -
                                       m_c[0][5] * gradNa_gradNb[0][0] - m_c[0][1] * gradNa_gradNb[0][1] - m_c[0][3] * gradNa_gradNb[0][2] -
                                       m_c[5][5] * gradNa_gradNb[1][0] - m_c[5][1] * gradNa_gradNb[1][1] - m_c[5][3] * gradNa_gradNb[1][2] -
                                       m_c[4][5] * gradNa_gradNb[2][0] - m_c[4][1] * gradNa_gradNb[2][1] - m_c[4][3] * gradNa_gradNb[2][2] -
                                       m_c[0][4] * gradNa_gradNb[0][0] - m_c[0][3] * gradNa_gradNb[0][1] - m_c[0][2] * gradNa_gradNb[0][2] -
                                       m_c[5][4] * gradNa_gradNb[1][0] - m_c[5][3] * gradNa_gradNb[1][1] - m_c[5][2] * gradNa_gradNb[1][2] -
                                       m_c[4][4] * gradNa_gradNb[2][0] - m_c[4][3] * gradNa_gradNb[2][1] - m_c[4][2] * gradNa_gradNb[2][2] ) * detJxW;
      diagSumLocalStiffness[a*3+1] = diagSumLocalStiffness[a*3+1] -
                                     (m_c[5][0] * gradNa_gradNb[0][0] - m_c[5][5] * gradNa_gradNb[0][1] - m_c[5][4] * gradNa_gradNb[0][2] -
                                      m_c[1][0] * gradNa_gradNb[1][0] - m_c[1][5] * gradNa_gradNb[1][1] - m_c[1][4] * gradNa_gradNb[1][2] -
                                      m_c[3][0] * gradNa_gradNb[2][0] - m_c[3][5] * gradNa_gradNb[2][1] - m_c[3][4] * gradNa_gradNb[2][2] -
                                      m_c[5][5] * gradNa_gradNb[0][0] - m_c[5][1] * gradNa_gradNb[0][1] - m_c[5][3] * gradNa_gradNb[0][2] -
                                      m_c[1][5] * gradNa_gradNb[1][0] - m_c[1][1] * gradNa_gradNb[1][1] - m_c[1][3] * gradNa_gradNb[1][2] -
                                      m_c[3][5] * gradNa_gradNb[2][0] - m_c[3][1] * gradNa_gradNb[2][1] - m_c[3][3] * gradNa_gradNb[2][2] -
                                      m_c[5][4] * gradNa_gradNb[0][0] - m_c[5][3] * gradNa_gradNb[0][1] - m_c[5][2] * gradNa_gradNb[0][2] -
                                      m_c[1][4] * gradNa_gradNb[1][0] - m_c[1][3] * gradNa_gradNb[1][1] - m_c[1][2] * gradNa_gradNb[1][2] -
                                      m_c[3][4] * gradNa_gradNb[2][0] - m_c[3][3] * gradNa_gradNb[2][1] - m_c[3][2] * gradNa_gradNb[2][2] ) * detJxW;
      diagSumLocalStiffness[a*3+2] = diagSumLocalStiffness[a*3+2] -
                                     (m_c[4][0] * gradNa_gradNb[0][0] - m_c[4][5] * gradNa_gradNb[0][1] - m_c[4][4] * gradNa_gradNb[0][2] -
                                      m_c[3][0] * gradNa_gradNb[1][0] - m_c[3][5] * gradNa_gradNb[1][1] - m_c[3][4] * gradNa_gradNb[1][2] -
                                      m_c[2][0] * gradNa_gradNb[2][0] - m_c[2][5] * gradNa_gradNb[2][1] - m_c[2][4] * gradNa_gradNb[2][2] -
                                      m_c[4][5] * gradNa_gradNb[0][0] - m_c[4][1] * gradNa_gradNb[0][1] - m_c[4][3] * gradNa_gradNb[0][2] -
                                      m_c[3][5] * gradNa_gradNb[1][0] - m_c[3][1] * gradNa_gradNb[1][1] - m_c[3][3] * gradNa_gradNb[1][2] -
                                      m_c[2][5] * gradNa_gradNb[2][0] - m_c[2][1] * gradNa_gradNb[2][1] - m_c[2][3] * gradNa_gradNb[2][2] -
                                      m_c[4][4] * gradNa_gradNb[0][0] - m_c[4][3] * gradNa_gradNb[0][1] - m_c[4][2] * gradNa_gradNb[0][2] -
                                      m_c[3][4] * gradNa_gradNb[1][0] - m_c[3][3] * gradNa_gradNb[1][1] - m_c[3][2] * gradNa_gradNb[1][2] -
                                      m_c[2][4] * gradNa_gradNb[2][0] - m_c[2][3] * gradNa_gradNb[2][1] - m_c[2][2] * gradNa_gradNb[2][2] ) * detJxW;
    }
  }
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERANISOTROPIC_HPP_ */
