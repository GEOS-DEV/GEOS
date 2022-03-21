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

#include "common/DataTypes.hpp"

namespace geosx
{

namespace compositionalMultiphaseUtilities
{

/**
 * @brief Shifts all elements one position ahead and replaces the first element
 * with the sum of all elements in each block of the input block one-dimensional
 * array of values.
 *
 * @tparam VEC type of one-dimensional array of values
 * @param N block size
 * @param NB number of blocks
 * @param v block one-dimensional array of values
 */
template< typename VEC >
GEOSX_HOST_DEVICE
void shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( integer const N,
                                                                integer const NB,
                                                                VEC && v )
{
  for( integer i = 0; i < NB; ++i )
  {
    integer const ind = i * N + N - 1;
    real64 tmp = v[ind];
    for( int j = ind - 1; j >= i * N; --j )
    {
      v[j+1] = v[j];
      tmp += v[j];
    }
    v[i*N] = tmp;
  }
}

/**
 * @brief Shifts all elements one position ahead and replaces the first element
 * with the sum of all elements in the input one-dimensional array of values.
 *
 * @tparam VEC type of one-dimensional array of values
 * @param N vector size
 * @param v one-dimensional array of values
 */
template< typename VEC >
GEOSX_HOST_DEVICE
void shiftElementsAheadByOneAndReplaceFirstElementWithSum( integer const N,
                                                           VEC && v )
{
  shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( N, 1, v );
}

/**
 * @brief Shifts all rows one position ahead and replaces the first row
 * with the sum of all rows in each block of the input block row
 * two-dimensional array of values.
 *
 * @tparam MATRIX type of two-dimensional array of values
 * @tparam VEC type of one-dimensional array of values
 * @param M block number of rows
 * @param N block number of columns
 * @param NB number of row blocks
 * @param mat block row two-dimensional array of values
 * @param work one-dimensional working array of values of size N
 */
template< typename MATRIX, typename VEC >
GEOSX_HOST_DEVICE
void shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( integer const M,
                                                              integer const N,
                                                              integer const NB,
                                                              MATRIX && mat,
                                                              VEC && work )
{
  for( integer k = 0; k < NB; ++k )
  {
    integer const ind = k * M + M - 1;
    for( integer j = 0; j < N; ++j )
    {
      work[j] = mat[ind][j];
    }
    for( integer i = ind - 1; i >= k * M; --i )
    {
      for( integer j = 0; j < N; ++j )
      {
        mat[i+1][j] = mat[i][j];
        work[j] += mat[i][j];
      }
    }
    for( integer j = 0; j < N; ++j )
    {
      mat[k*M][j] = work[j];
    }
  }
}

/**
 * @brief Shifts all rows one position ahead and replaces the first row
 * with the sum of all rows in the input block row two-dimensional array
 * of values.
 *
 * @tparam MATRIX type of two-dimensional array of values
 * @tparam VEC type of one-dimensional array of values
 * @param M number of rows
 * @param N number of columns
 * @param mat block row two-dimensional array of values
 * @param work one-dimensional working array of values of size N
 */
template< typename MATRIX, typename VEC >
GEOSX_HOST_DEVICE
void shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( integer const M,
                                                         integer const N,
                                                         MATRIX && mat,
                                                         VEC && work )
{
  shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( M, N, 1, mat, work );
}

} // namespace compositionalMultiphaseUtilities

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEUTILITIES_HPP_
