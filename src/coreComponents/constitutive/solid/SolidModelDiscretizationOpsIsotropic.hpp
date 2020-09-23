/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidModelHelperIsotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSISOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSISOTROPIC_HPP_

#include "SolidModelDiscretizationOps.hpp"

namespace geosx
{
namespace constitutive
{


struct SolidModelDiscretizationOpsIsotropic : public SolidModelDiscretizationOps
{
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void upperBTDB( BASIS_GRADIENT const & gradN,
                  real64 const & detJxW,
                  real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 const & detJxW,
                 real64 ( &diagElementStiffness )[NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 const & detJxW,
                       real64 ( &diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void scaleParams( real64 const scale )
  {
    m_lambda *= scale;
    m_shearModulus *= scale;
  }

  real64 m_lambda;
  real64 m_shearModulus;
};

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOpsIsotropic::upperBTDB( BASIS_GRADIENT const & gradN,
                                                      real64 const & detJxW,
                                                      real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJxW;
  real64 const lambda = this->m_lambda * detJxW;
  real64 const G = this->m_shearModulus * detJxW;
  SolidModelDiscretizationOps::upperBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                elementStiffness,
                                                                [ lambda,
                                                                  G,
                                                                  lambda2G ] GEOSX_HOST_DEVICE
                                                                  ( int const a,
                                                                  int const b,
                                                                  real64 const (&gradNa_gradNb)[3][3],
                                                                  real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    elementStiffness[a*3+0][b*3+0] = elementStiffness[a*3+0][b*3+0] + gradNa_gradNb[1][1] * G + gradNa_gradNb[2][2] * G + gradNa_gradNb[0][0] * lambda2G;
    elementStiffness[a*3+0][b*3+1] = elementStiffness[a*3+0][b*3+1] + gradNa_gradNb[1][0] * G + gradNa_gradNb[0][1] * lambda;
    elementStiffness[a*3+0][b*3+2] = elementStiffness[a*3+0][b*3+2] + gradNa_gradNb[2][0] * G + gradNa_gradNb[0][2] * lambda;
    elementStiffness[a*3+1][b*3+0] = elementStiffness[a*3+1][b*3+0] + gradNa_gradNb[0][1] * G + gradNa_gradNb[1][0] * lambda;
    elementStiffness[a*3+1][b*3+1] = elementStiffness[a*3+1][b*3+1] + gradNa_gradNb[0][0] * G + gradNa_gradNb[2][2] * G + gradNa_gradNb[1][1] * lambda2G;
    elementStiffness[a*3+1][b*3+2] = elementStiffness[a*3+1][b*3+2] + gradNa_gradNb[2][1] * G + gradNa_gradNb[1][2] * lambda;
    elementStiffness[a*3+2][b*3+0] = elementStiffness[a*3+2][b*3+0] + gradNa_gradNb[0][2] * G + gradNa_gradNb[2][0] * lambda;
    elementStiffness[a*3+2][b*3+1] = elementStiffness[a*3+2][b*3+1] + gradNa_gradNb[1][2] * G + gradNa_gradNb[2][1] * lambda;
    elementStiffness[a*3+2][b*3+2] = elementStiffness[a*3+2][b*3+2] + gradNa_gradNb[0][0] * G + gradNa_gradNb[1][1] * G + gradNa_gradNb[2][2] * lambda2G;
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOpsIsotropic::diagBTDB( BASIS_GRADIENT const & gradN,
                                                     real64 const & detJxW,
                                                     real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS *3] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJxW;
  real64 const lambda = this->m_lambda * detJxW;
  real64 const G = this->m_shearModulus * detJxW;

  SolidModelDiscretizationOps::diagBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                               diagElementStiffness,
                                                               [ lambda,
                                                                 G,
                                                                 lambda2G ] GEOSX_HOST_DEVICE
                                                                 ( const int a,
                                                                 const real64 (& gradN_gradN)[3],
                                                                 real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagElementStiffness[ a*3+0 ] = diagElementStiffness[ a*3+0 ] + gradN_gradN[0] * lambda2G + gradN_gradN[1] * G + gradN_gradN[2] * G;
    diagElementStiffness[ a*3+1 ] = diagElementStiffness[ a*3+1 ] + gradN_gradN[1] * lambda2G + gradN_gradN[0] * G + gradN_gradN[2] * G;
    diagElementStiffness[ a*3+2 ] = diagElementStiffness[ a*3+2 ] + gradN_gradN[2] * lambda2G + gradN_gradN[0] * G + gradN_gradN[1] * G;
  } );
}



template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOpsIsotropic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                           real64 const & detJxW,
                                                           real64 ( & diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] )
{
  real64 const lambda2G = ( 2 * m_shearModulus + m_lambda ) * detJxW;
  real64 const lambda = this->m_lambda * detJxW;
  real64 const G = this->m_shearModulus * detJxW;
  SolidModelDiscretizationOps::diagRowSumBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                     diagSumElementStiffness,
                                                                     [ lambda,
                                                                       G,
                                                                       lambda2G ] GEOSX_HOST_DEVICE
                                                                       ( int const a,
                                                                       real64 const (&gradNa_gradNb)[3][3],
                                                                       real64 (& diagSumElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagSumElementStiffness[a*3+0] = diagSumElementStiffness[a*3+0] +
                                     gradNa_gradNb[1][1] * G +
                                     gradNa_gradNb[2][2] * G +
                                     gradNa_gradNb[0][0] * lambda2G +
                                     gradNa_gradNb[1][0] * G +
                                     gradNa_gradNb[0][1] * lambda +
                                     gradNa_gradNb[2][0] * G +
                                     gradNa_gradNb[0][2] * lambda;
    diagSumElementStiffness[a*3+1] = diagSumElementStiffness[a*3+1] +
                                     gradNa_gradNb[0][1] * G +
                                     gradNa_gradNb[1][0] * lambda +
                                     gradNa_gradNb[0][0] * G +
                                     gradNa_gradNb[2][2] * G +
                                     gradNa_gradNb[1][1] * lambda2G +
                                     gradNa_gradNb[2][1] * G +
                                     gradNa_gradNb[1][2] * lambda;
    diagSumElementStiffness[a*3+2] = diagSumElementStiffness[a*3+2] +
                                     gradNa_gradNb[0][2] * G +
                                     gradNa_gradNb[2][0] * lambda +
                                     gradNa_gradNb[1][2] * G +
                                     gradNa_gradNb[2][1] * lambda +
                                     gradNa_gradNb[0][0] * G +
                                     gradNa_gradNb[1][1] * G +
                                     gradNa_gradNb[2][2] * lambda2G;
  } );
}
#if __GNUC__
#pragma GCC diagnostic pop
#endif

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSISOTROPIC_HPP_ */
