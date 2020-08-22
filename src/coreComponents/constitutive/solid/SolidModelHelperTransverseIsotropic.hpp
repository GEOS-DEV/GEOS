/*
 * BTransposeDB_Helper.hpp
 *
 *  Created on: Aug 11, 2020
 *      Author: settgast
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERTRANSVERSEISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERTRANSVERSEISOTROPIC_HPP_

#include "SolidModelHelperBase.hpp"

namespace geosx
{
namespace constitutive
{


struct SolidModelHelperTransverseIsotropic : public SolidModelHelperBase
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
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperTransverseIsotropic::BTDB( BASIS_GRADIENT const & gradN,
                                                real64 const & detJxW,
                                                real64 (& localStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const c12 = (m_c11 - 2 * m_c66) * detJxW;
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c13 = this->m_c13 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;

  SolidModelHelperBase::BTDB< NUM_SUPPORT_POINTS >( gradN,
                                                    localStiffness,
                                                    [ c11,
                                                      c13,
                                                      c33,
                                                      c44,
                                                      c66,
                                                      c12 ] GEOSX_HOST_DEVICE
                                                      ( int const a,
                                                      int const b,
                                                      real64 const (&gradNa_gradNb)[3][3],
                                                      real64 (& localStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    localStiffness[a*3+0][b*3+0] = localStiffness[a*3+0][b*3+0] - c11 * gradNa_gradNb[0][0] - c66 * gradNa_gradNb[1][1] - c44 * gradNa_gradNb[2][2];
    localStiffness[a*3+0][b*3+1] = localStiffness[a*3+0][b*3+1] - c12 * gradNa_gradNb[0][1] - c66 * gradNa_gradNb[1][0];
    localStiffness[a*3+0][b*3+2] = localStiffness[a*3+0][b*3+2] - c13 * gradNa_gradNb[0][2] - c44 * gradNa_gradNb[2][0];
    localStiffness[a*3+1][b*3+0] = localStiffness[a*3+1][b*3+0] - c66 * gradNa_gradNb[0][1] - c12 * gradNa_gradNb[1][0];
    localStiffness[a*3+1][b*3+1] = localStiffness[a*3+1][b*3+1] - c66 * gradNa_gradNb[0][0] - c11 * gradNa_gradNb[1][1] - c44 * gradNa_gradNb[2][2];
    localStiffness[a*3+1][b*3+2] = localStiffness[a*3+1][b*3+2] - c13 * gradNa_gradNb[1][2] - c44 * gradNa_gradNb[2][1];
    localStiffness[a*3+2][b*3+0] = localStiffness[a*3+2][b*3+0] - c44 * gradNa_gradNb[0][2] - c13 * gradNa_gradNb[2][0];
    localStiffness[a*3+2][b*3+1] = localStiffness[a*3+2][b*3+1] - c44 * gradNa_gradNb[1][2] - c13 * gradNa_gradNb[2][1];
    localStiffness[a*3+2][b*3+2] = localStiffness[a*3+2][b*3+2] - c44 * gradNa_gradNb[0][0] - c44 * gradNa_gradNb[1][1] - c33 * gradNa_gradNb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperTransverseIsotropic::diagBTDB( BASIS_GRADIENT const & gradN,
                                                    real64 const & detJxW,
                                                    real64 (& diagLocalStiffness)[NUM_SUPPORT_POINTS *3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;
  SolidModelHelperBase::diagBTDB< NUM_SUPPORT_POINTS,
                                  3 >( gradN,
                                       diagLocalStiffness,
                                       [ c11,
                                         c33,
                                         c44,
                                         c66 ] GEOSX_HOST_DEVICE
                                         ( const int a,
                                         const real64 (& gradN_gradN)[3],
                                         real64 (& diagLocalStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagLocalStiffness[a*3+0] = diagLocalStiffness[a*3+0] - c11 * gradN_gradN[0] - c66 * gradN_gradN[1] - c44 * gradN_gradN[2];
    diagLocalStiffness[a*3+1] = diagLocalStiffness[a*3+1] - c66 * gradN_gradN[0] - c11 * gradN_gradN[1] - c44 * gradN_gradN[2];
    diagLocalStiffness[a*3+2] = diagLocalStiffness[a*3+2] - c44 * gradN_gradN[0] - c44 * gradN_gradN[1] - c33 * gradN_gradN[2];
  } );
}

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperTransverseIsotropic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                          real64 const & detJxW,
                                                          real64 ( & diagSumLocalStiffness )[NUM_SUPPORT_POINTS*3] )
{
  real64 const c12 = (m_c11 - 2 * m_c66) * detJxW;
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c13 = this->m_c13 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;

  SolidModelHelperBase::diagRowSumBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                              diagSumLocalStiffness,
                                                              [ c11,
                                                                c13,
                                                                c33,
                                                                c44,
                                                                c66,
                                                                c12 ] GEOSX_HOST_DEVICE
                                                                ( int const a,
                                                                real64 const (&gradNa_gradNb)[3][3],
                                                                real64 (& diagSumLocalStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagSumLocalStiffness[a*3+0] = diagSumLocalStiffness[a*3+0] -
                                   c11 * gradNa_gradNb[0][0] -
                                   c66 * gradNa_gradNb[1][1] -
                                   c44 * gradNa_gradNb[2][2] -
                                   c12 * gradNa_gradNb[0][1] -
                                   c66 * gradNa_gradNb[1][0] -
                                   c13 * gradNa_gradNb[0][2] -
                                   c44 * gradNa_gradNb[2][0];
    diagSumLocalStiffness[a*3+1] = diagSumLocalStiffness[a*3+1] -
                                   c66 * gradNa_gradNb[0][1] -
                                   c12 * gradNa_gradNb[1][0] -
                                   c66 * gradNa_gradNb[0][0] -
                                   c11 * gradNa_gradNb[1][1] -
                                   c44 * gradNa_gradNb[2][2] -
                                   c13 * gradNa_gradNb[1][2] -
                                   c44 * gradNa_gradNb[2][1];
    diagSumLocalStiffness[a*3+2] = diagSumLocalStiffness[a*3+2] -
                                   c44 * gradNa_gradNb[0][2] -
                                   c13 * gradNa_gradNb[2][0] -
                                   c44 * gradNa_gradNb[1][2] -
                                   c13 * gradNa_gradNb[2][1] -
                                   c44 * gradNa_gradNb[0][0] -
                                   c44 * gradNa_gradNb[1][1] -
                                   c33 * gradNa_gradNb[2][2];
  } );
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERTRANSVERSEISOTROPIC_HPP_ */
