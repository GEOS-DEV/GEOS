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
 * @file SolidModelHelperBase.hpp
 */


#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPS_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPS_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{
namespace constitutive
{


struct SolidModelDiscretizationOps
{
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT,
            typename CBF >
  GEOSX_HOST_DEVICE
  void BTDB( BASIS_GRADIENT const & gradN,
             real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3],
             CBF && callbackFunction );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT,
            typename CBF >
  GEOSX_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 ( &diagElementStiffness )[NUM_SUPPORT_POINTS*3],
                 CBF && callbackFunction );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT,
            typename CBF >
  GEOSX_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 ( &diagSumElementStiffness )[NUM_SUPPORT_POINTS*3],
                       CBF && callbackFunction );

};


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT,
          typename CBF >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOps::BTDB( BASIS_GRADIENT const & gradN,
                                 real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3],
                                 CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };
      callbackFunction( a, b, gradNa_gradNb, elementStiffness );
    }
  }

  for( int row=0; row<NUM_SUPPORT_POINTS*3; ++row )
  {
    for( int col=0; col<row; ++col )
    {
      elementStiffness[row][col] = elementStiffness[col][row];
    }
  }
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT,
          typename CBF >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOps::diagBTDB( BASIS_GRADIENT const & gradN,
                                     real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS*3],
                                     CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    real64 const gradN_gradN[3] = { gradN[a][0] * gradN[a][0], gradN[a][1] * gradN[a][1], gradN[a][2] * gradN[a][2] };

    callbackFunction( a, gradN_gradN, diagElementStiffness );
  }
}

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT,
          typename CBF >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOps::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                           real64 ( & diagSumElementStiffness )[NUM_SUPPORT_POINTS*3],
                                           CBF && callbackFunction )
{
  for( int a=0; a<NUM_SUPPORT_POINTS; ++a )
  {
    for( int b=a; b<NUM_SUPPORT_POINTS; ++b )
    {
      real64 const gradNa_gradNb[3][3] =
      { { gradN[a][0] * gradN[b][0], gradN[a][0] * gradN[b][1], gradN[a][0] * gradN[b][2] },
        { gradN[a][1] * gradN[b][0], gradN[a][1] * gradN[b][1], gradN[a][1] * gradN[b][2] },
        { gradN[a][2] * gradN[b][0], gradN[a][2] * gradN[b][1], gradN[a][2] * gradN[b][2] } };
      callbackFunction( a, b, gradNa_gradNb, diagSumElementStiffness );
    }
  }
}

}
}


#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPS_HPP_ */
