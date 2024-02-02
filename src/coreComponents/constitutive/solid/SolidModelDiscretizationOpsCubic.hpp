/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidModelHelperCubic.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSCUBIC_HPP_
#define GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSCUBIC_HPP_

#include "SolidModelDiscretizationOps.hpp"

namespace geos
{
namespace constitutive
{

struct SolidModelDiscretizationOpsCubic : public SolidModelDiscretizationOps
{
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void BTDB( BASIS_GRADIENT const & gradN,
             real64 const & detJxW,
             real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void upperBTDB( BASIS_GRADIENT const & gradN,
                  real64 const & detJxW,
                  real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 const & detJxW,
                 real64 ( &diagElementStiffness )[NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOS_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 const & detJxW,
                       real64 ( &diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] );

  GEOS_HOST_DEVICE
  inline
  void scaleParams( real64 const scale )
  {
    m_c11 *= scale;
    m_c12 *= scale;
    m_c44 *= scale;
  }

  real64 m_c11;
  real64 m_c12;
  real64 m_c44;

};

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsCubic::
  BTDB( BASIS_GRADIENT const & gradN,
        real64 const & detJxW,
        real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c12 = this->m_c12 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;

  SolidModelDiscretizationOps::BTDB< NUM_SUPPORT_POINTS >( gradN,
                                                           elementStiffness,
                                                           [ c11,
                                                             c12,
                                                             c44 ] GEOS_HOST_DEVICE
                                                             ( int const a,
                                                             int const b,
                                                             real64 const (&gradNa_gradNb)[3][3],
                                                             real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    elementStiffness[a*3+0][b*3+0] += c11 * gradNa_gradNb[0][0] + c44 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+0][b*3+1] += c12 * gradNa_gradNb[0][1] + c44 * gradNa_gradNb[1][0];
    elementStiffness[a*3+0][b*3+2] += c12 * gradNa_gradNb[0][2] + c44 * gradNa_gradNb[2][0];
    elementStiffness[a*3+1][b*3+0] += c44 * gradNa_gradNb[0][1] + c12 * gradNa_gradNb[1][0];
    elementStiffness[a*3+1][b*3+1] += c44 * gradNa_gradNb[0][0] + c11 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+1][b*3+2] += c12 * gradNa_gradNb[1][2] + c44 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+0] += c44 * gradNa_gradNb[0][2] + c12 * gradNa_gradNb[2][0];
    elementStiffness[a*3+2][b*3+1] += c44 * gradNa_gradNb[1][2] + c12 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+2] += c44 * gradNa_gradNb[0][0] + c44 * gradNa_gradNb[1][1] + c11 * gradNa_gradNb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsCubic::
  upperBTDB( BASIS_GRADIENT const & gradN,
             real64 const & detJxW,
             real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c12 = this->m_c12 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;

  SolidModelDiscretizationOps::upperBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                elementStiffness,
                                                                [ c11,
                                                                  c12,
                                                                  c44 ] GEOS_HOST_DEVICE
                                                                  ( int const a,
                                                                  int const b,
                                                                  real64 const (&gradNa_gradNb)[3][3],
                                                                  real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    elementStiffness[a*3+0][b*3+0] += c11 * gradNa_gradNb[0][0] + c44 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+0][b*3+1] += c12 * gradNa_gradNb[0][1] + c44 * gradNa_gradNb[1][0];
    elementStiffness[a*3+0][b*3+2] += c12 * gradNa_gradNb[0][2] + c44 * gradNa_gradNb[2][0];
    elementStiffness[a*3+1][b*3+0] += c44 * gradNa_gradNb[0][1] + c12 * gradNa_gradNb[1][0];
    elementStiffness[a*3+1][b*3+1] += c44 * gradNa_gradNb[0][0] + c11 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+1][b*3+2] += c12 * gradNa_gradNb[1][2] + c44 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+0] += c44 * gradNa_gradNb[0][2] + c12 * gradNa_gradNb[2][0];
    elementStiffness[a*3+2][b*3+1] += c44 * gradNa_gradNb[1][2] + c12 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+2] += c44 * gradNa_gradNb[0][0] + c44 * gradNa_gradNb[1][1] + c11 * gradNa_gradNb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsCubic::diagBTDB( BASIS_GRADIENT const & gradN,
                                                       real64 const & detJxW,
                                                       real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS *3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  SolidModelDiscretizationOps::diagBTDB< NUM_SUPPORT_POINTS,
                                         3 >( gradN,
                                              diagElementStiffness,
                                              [ c11,
                                                c44 ] GEOS_HOST_DEVICE
                                                ( const int a,
                                                const real64 (& gradN_gradN)[3],
                                                real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagElementStiffness[a*3+0] += c11 * gradN_gradN[0] + c44 * gradN_gradN[1] + c44 * gradN_gradN[2];
    diagElementStiffness[a*3+1] += c44 * gradN_gradN[0] + c11 * gradN_gradN[1] + c44 * gradN_gradN[2];
    diagElementStiffness[a*3+2] += c44 * gradN_gradN[0] + c44 * gradN_gradN[1] + c11 * gradN_gradN[2];
  } );
}

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOS_HOST_DEVICE
inline
void SolidModelDiscretizationOpsCubic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                             real64 const & detJxW,
                                                             real64 ( & diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c12 = this->m_c12 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;

  SolidModelDiscretizationOps::diagRowSumBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                     diagSumElementStiffness,
                                                                     [ c11,
                                                                       c12,
                                                                       c44 ] GEOS_HOST_DEVICE
                                                                       ( int const a,
                                                                       real64 const (&gradNa_gradNb)[3][3],
                                                                       real64 (& diagSumElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagSumElementStiffness[a*3+0] = diagSumElementStiffness[a*3+0] +
                                     c11 * gradNa_gradNb[0][0] +
                                     c44 * gradNa_gradNb[1][1] +
                                     c44 * gradNa_gradNb[2][2] +
                                     c12 * gradNa_gradNb[0][1] +
                                     c44 * gradNa_gradNb[1][0] +
                                     c12 * gradNa_gradNb[0][2] +
                                     c44 * gradNa_gradNb[2][0];

    diagSumElementStiffness[a*3+1] = diagSumElementStiffness[a*3+1] +
                                     c44 * gradNa_gradNb[0][1] +
                                     c12 * gradNa_gradNb[1][0] +
                                     c44 * gradNa_gradNb[0][0] +
                                     c11 * gradNa_gradNb[1][1] +
                                     c44 * gradNa_gradNb[2][2] +
                                     c12 * gradNa_gradNb[1][2] +
                                     c44 * gradNa_gradNb[2][1];

    diagSumElementStiffness[a*3+2] = diagSumElementStiffness[a*3+2] +
                                     c44 * gradNa_gradNb[0][2] +
                                     c12 * gradNa_gradNb[2][0] +
                                     c44 * gradNa_gradNb[1][2] +
                                     c12 * gradNa_gradNb[2][1] +
                                     c44 * gradNa_gradNb[0][0] +
                                     c44 * gradNa_gradNb[1][1] +
                                     c11 * gradNa_gradNb[2][2];
  } );
}

#if __GNUC__
#pragma GCC diagnostic pop
#endif

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSCUBIC_HPP_ */
