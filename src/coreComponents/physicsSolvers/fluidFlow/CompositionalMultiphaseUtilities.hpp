/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseUtilities.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEUTILITIES_H_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEUTILITIES_H_

//#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
//#include "LvArray/src/limits.hpp"

namespace geosx
{

namespace CompositionalMultiphaseUtilities
{

template< typename VEC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyBlockLinearCombination( localIndex const M,
                                  localIndex const a,
                                  VEC && v )
{
  for( localIndex i = 0; i < a; ++i )
  {
    localIndex const ind = i * M + M - 1;
    real64 tmp = v[ind];
    for( int j = ind - 1; j >= i * M; --j )
    {
      v[j+1] = v[j];
      tmp += v[j];
    }
    v[i*M] = tmp;
  }
}

template< typename MATRIX, typename VEC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyBlockLinearCombination( localIndex const M,
                                  localIndex const N,
                                  localIndex const a,
                                  localIndex const b,
                                  MATRIX && mat,
                                  VEC && work )
{
  for( localIndex k = 0; k < a; ++k )
  {
    localIndex const ind = k * M + M - 1;
//    copy( N * b, mat[ind], work );
    for( localIndex j = 0; j < b * N; ++j )
    {
      work[j] = mat[ind][j];
    }
    for( localIndex i = ind - 1; i >= k * M; --i )
    {
      for( localIndex j = 0; j < b * N; ++j )
      {
        mat[i+1][j] = mat[i][j];
        work[j] += mat[i][j];
      }
    }
    for( localIndex j = 0; j < b * N; ++j )
    {
      mat[k*M][j] = work[j];
    }
  }
}

} // namespace CompositionalMultiphaseUtilities

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASE_HPP_
