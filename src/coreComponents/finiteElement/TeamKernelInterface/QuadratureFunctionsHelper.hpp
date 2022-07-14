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
 * @file QuadratureFunctionsHelper.hpp
 */

#ifndef GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_
#define GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
// namespace finiteElement
// {

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 determinant( real64 const (& J)[3][3] )
{
  real64 const detJ = J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2])
                    - J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2])
                    + J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
  return detJ;
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void adjugate( real64 const (& J)[3][3], real64 (& AdjJ)[3][3] )
{
  AdjJ[0][0] = (J[1][1] * J[2][2]) - (J[1][2] * J[2][1]);
  AdjJ[0][1] = (J[2][1] * J[0][2]) - (J[0][1] * J[2][2]);
  AdjJ[0][2] = (J[0][1] * J[1][2]) - (J[1][1] * J[0][2]);
  AdjJ[1][0] = (J[2][0] * J[1][2]) - (J[1][0] * J[2][2]);
  AdjJ[1][1] = (J[0][0] * J[2][2]) - (J[0][2] * J[2][0]);
  AdjJ[1][2] = (J[1][0] * J[0][2]) - (J[0][0] * J[1][2]);
  AdjJ[2][0] = (J[1][0] * J[2][1]) - (J[2][0] * J[1][1]);
  AdjJ[2][1] = (J[2][0] * J[0][1]) - (J[0][0] * J[2][1]);
  AdjJ[2][2] = (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);
}

template < typename QField >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalLoad( localIndex const quad_x,
                 localIndex const quad_y,
                 localIndex const quad_z,
                 QField const & q_field,
                 real64 & q_value )
{
  q_value = q_field[quad_x][quad_y][quad_z];
}

template < typename QField, size_t num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalLoad( localIndex const quad_x,
                 localIndex const quad_y,
                 localIndex const quad_z,
                 QField const & q_field,
                 real64 (& q_value)[num_comp] )
{
  for (size_t c = 0; c < num_comp; c++)
  {
    q_value[c] = q_field[quad_x][quad_y][quad_z][c];
  }
}

template < typename QField, size_t num_comp_x, size_t num_comp_y >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalLoad( localIndex const quad_x,
                 localIndex const quad_y,
                 localIndex const quad_z,
                 QField const & q_field,
                 real64 (& q_value)[num_comp_x][num_comp_y] )
{
  for (size_t c_x = 0; c_x < num_comp_x; c_x++)
  {
    for (size_t c_y = 0; c_y < num_comp_y; c_y++)
    {
      q_value[c_x][c_y] = q_field[quad_x][quad_y][quad_z][c_x][c_y];
    }
  }
}

template < typename QField >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( localIndex const quad_x,
                  localIndex const quad_y,
                  localIndex const quad_z,
                  real64 const & q_value,
                  QField & q_field )
{
  q_field[quad_x][quad_y][quad_z] = q_value;
}

template < typename QField, size_t num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( localIndex const quad_x,
                  localIndex const quad_y,
                  localIndex const quad_z,
                  real64 const (& q_value)[num_comp],
                  QField & q_field )
{
  for (size_t c = 0; c < num_comp; c++)
  {
    q_field[quad_x][quad_y][quad_z][c] = q_value[c];
  }
}

template < typename QField, size_t num_comp_x, size_t num_comp_y >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( localIndex const quad_x,
                  localIndex const quad_y,
                  localIndex const quad_z,
                  real64 const (& q_value)[num_comp_x][num_comp_y],
                  QField & q_field )
{
  for (size_t c_x = 0; c_x < num_comp_x; c_x++)
  {
    for (size_t c_y = 0; c_y < num_comp_y; c_y++)
    {
      q_field[quad_x][quad_y][quad_z][c_x][c_y] = q_value[c_x][c_y];
    }
  }
}

// } // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONSHELPER_HPP_ */
