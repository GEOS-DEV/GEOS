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
 * @file qLocalWrite.hpp
 */

#ifndef GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONS_WRITE_HPP_
#define GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONS_WRITE_HPP_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "finiteElement/TeamKernelInterface/common.hpp"
#include "tensor/tensor_types.hpp"

namespace geosx
{

// Non-distributed/Shared version
template < localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const & q_value,
                  real64 (& q_field)[num_quads_1d][num_quads_1d][num_quads_1d] )
{
  q_field[quad_index.x][quad_index.y][quad_index.z] = q_value;
}

template < localIndex num_quads_1d, localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const (& q_value)[num_comp],
                  real64 (& q_field)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp] )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    q_field[quad_index.x][quad_index.y][quad_index.z][c] = q_value[c];
  }
}

template < localIndex num_quads_1d, localIndex num_comp_x, localIndex num_comp_y >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const (& q_value)[num_comp_x][num_comp_y],
                  real64 (& q_field)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp_x][num_comp_y] )
{
  for (localIndex c_x = 0; c_x < num_comp_x; c_x++)
  {
    for (localIndex c_y = 0; c_y < num_comp_y; c_y++)
    {
      q_field[quad_index.x][quad_index.y][quad_index.z][c_x][c_y] = q_value[c_x][c_y];
    }
  }
}

// 2D distributed (x ,y)
template < localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const & q_value,
                  tensor::Static2dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d > & q_field )
{
  q_field( quad_index.x, quad_index.y, quad_index.z ) = q_value;
}

template < localIndex num_quads_1d, localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const (& q_value)[num_comp],
                  tensor::Static2dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp > & q_field )
{
  for (localIndex c = 0; c < num_comp; c++)
  {
    q_field( quad_index.x, quad_index.y, quad_index.z, c ) = q_value[c];
  }
}

template < localIndex num_quads_1d, localIndex num_comp_x, localIndex num_comp_y >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const (& q_value)[num_comp_x][num_comp_y],
                  tensor::Static2dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp_x, num_comp_y > & q_field )
{
  for (localIndex c_x = 0; c_x < num_comp_x; c_x++)
  {
    for (localIndex c_y = 0; c_y < num_comp_y; c_y++)
    {
      q_field( quad_index.x, quad_index.y, quad_index.z, c_x, c_y ) = q_value[c_x][c_y];
    }
  }
}

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONS_WRITE_HPP_ */
