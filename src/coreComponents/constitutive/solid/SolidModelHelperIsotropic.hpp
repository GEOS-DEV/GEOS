/*
 * BTransposeDB_Helper.hpp
 *
 *  Created on: Aug 11, 2020
 *      Author: settgast
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERISOTROPIC_HPP_

#include "SolidModelHelperBase.hpp"

namespace geosx
{
namespace constitutive
{


struct SolidModelHelperIsotropic : public SolidModelHelperBase
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
    m_lambda = C[1][0];
    m_shearModulus = C[5][5];
  }

  real64 m_lambda;
  real64 m_shearModulus;
};


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperIsotropic::BTDB( BASIS_GRADIENT const & gradN,
                                      real64 const & detJxW,
                                      real64 (& localStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJxW;
  real64 const lambda = this->m_lambda * detJxW;
  real64 const G = this->m_shearModulus * detJxW;
  SolidModelHelperBase::BTDB< NUM_SUPPORT_POINTS >( gradN,
                                                    localStiffness,
                                                    [ lambda,
                                                      G,
                                                      lambda2G ] GEOSX_HOST_DEVICE
                                                      ( int const a,
                                                      int const b,
                                                      real64 const (&gradNa_gradNb)[3][3],
                                                      real64 (& localStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    localStiffness[a*3+0][b*3+0] = localStiffness[a*3+0][b*3+0] - gradNa_gradNb[1][1] * G - gradNa_gradNb[2][2] * G - gradNa_gradNb[0][0] * lambda2G;
    localStiffness[a*3+0][b*3+1] = localStiffness[a*3+0][b*3+1] - gradNa_gradNb[1][0] * G - gradNa_gradNb[0][1] * lambda;
    localStiffness[a*3+0][b*3+2] = localStiffness[a*3+0][b*3+2] - gradNa_gradNb[2][0] * G - gradNa_gradNb[0][2] * lambda;
    localStiffness[a*3+1][b*3+0] = localStiffness[a*3+1][b*3+0] - gradNa_gradNb[0][1] * G - gradNa_gradNb[1][0] * lambda;
    localStiffness[a*3+1][b*3+1] = localStiffness[a*3+1][b*3+1] - gradNa_gradNb[0][0] * G - gradNa_gradNb[2][2] * G - gradNa_gradNb[1][1] * lambda2G;
    localStiffness[a*3+1][b*3+2] = localStiffness[a*3+1][b*3+2] - gradNa_gradNb[2][1] * G - gradNa_gradNb[1][2] * lambda;
    localStiffness[a*3+2][b*3+0] = localStiffness[a*3+2][b*3+0] - gradNa_gradNb[0][2] * G - gradNa_gradNb[2][0] * lambda;
    localStiffness[a*3+2][b*3+1] = localStiffness[a*3+2][b*3+1] - gradNa_gradNb[1][2] * G - gradNa_gradNb[2][1] * lambda;
    localStiffness[a*3+2][b*3+2] = localStiffness[a*3+2][b*3+2] - gradNa_gradNb[0][0] * G - gradNa_gradNb[1][1] * G - gradNa_gradNb[2][2] * lambda2G;
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperIsotropic::diagBTDB( BASIS_GRADIENT const & gradN,
                                          real64 const & detJxW,
                                          real64 (& diagLocalStiffness)[NUM_SUPPORT_POINTS *3] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJxW;
  real64 const lambda = this->m_lambda * detJxW;
  real64 const G = this->m_shearModulus * detJxW;

  SolidModelHelperBase::diagBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                        diagLocalStiffness,
                                                        [ lambda,
                                                          G,
                                                          lambda2G ] GEOSX_HOST_DEVICE
                                                          ( const int a,
                                                          const real64 (& gradN_gradN)[3],
                                                          real64 (& diagLocalStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagLocalStiffness[ a*3+0 ] = diagLocalStiffness[ a*3+0 ] - gradN_gradN[0] * lambda2G - gradN_gradN[1] * G - gradN_gradN[2] * G;
    diagLocalStiffness[ a*3+1 ] = diagLocalStiffness[ a*3+1 ] - gradN_gradN[1] * lambda2G - gradN_gradN[0] * G - gradN_gradN[2] * G;
    diagLocalStiffness[ a*3+2 ] = diagLocalStiffness[ a*3+2 ] - gradN_gradN[2] * lambda2G - gradN_gradN[0] * G - gradN_gradN[1] * G;
  } );
}



template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelHelperIsotropic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                real64 const & detJxW,
                                                real64 ( & diagSumLocalStiffness )[NUM_SUPPORT_POINTS*3] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJxW;
  real64 const lambda = this->m_lambda * detJxW;
  real64 const G = this->m_shearModulus * detJxW;
  SolidModelHelperBase::diagRowSumBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                              diagSumLocalStiffness,
                                                              [ lambda,
                                                                G,
                                                                lambda2G ] GEOSX_HOST_DEVICE
                                                                ( int const a,
                                                                real64 const (&gradNa_gradNb)[3][3],
                                                                real64 (& diagSumLocalStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagSumLocalStiffness[a*3+0] = diagSumLocalStiffness[a*3+0] -
                                   gradNa_gradNb[1][1] * G -
                                   gradNa_gradNb[2][2] * G -
                                   gradNa_gradNb[0][0] * lambda2G -
                                   gradNa_gradNb[1][0] * G -
                                   gradNa_gradNb[0][1] * lambda -
                                   gradNa_gradNb[2][0] * G -
                                   gradNa_gradNb[0][2] * lambda;
    diagSumLocalStiffness[a*3+1] = diagSumLocalStiffness[a*3+1] -
                                   gradNa_gradNb[0][1] * G -
                                   gradNa_gradNb[1][0] * lambda -
                                   gradNa_gradNb[0][0] * G -
                                   gradNa_gradNb[2][2] * G -
                                   gradNa_gradNb[1][1] * lambda2G -
                                   gradNa_gradNb[2][1] * G -
                                   gradNa_gradNb[1][2] * lambda;
    diagSumLocalStiffness[a*3+2] = diagSumLocalStiffness[a*3+2] -
                                   gradNa_gradNb[0][2] * G -
                                   gradNa_gradNb[2][0] * lambda -
                                   gradNa_gradNb[1][2] * G -
                                   gradNa_gradNb[2][1] * lambda -
                                   gradNa_gradNb[0][0] * G -
                                   gradNa_gradNb[1][1] * G -
                                   gradNa_gradNb[2][2] * lambda2G;
  } );
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELHELPERISOTROPIC_HPP_ */
