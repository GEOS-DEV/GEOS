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
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/common.hpp"
#include "tensor/tensor_types.hpp"

namespace geosx
{

/**
 * @brief Write a scalar field value local to the quadrature point.
 * 
 * @tparam Tensor The type of tensor.
 * @param quad_index The index of the quadrature point.
 * @param q_value The field local to the quadrature point.
 * @param q_field The field at all quadrature points.
 */
template < typename Tensor,
           std::enable_if_t<
             tensor::get_tensor_rank< Tensor > == 3,
             bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const & q_value,
                  Tensor & q_field )
{
  q_field( quad_index.x, quad_index.y, quad_index.z ) = q_value;
}

/**
 * @brief Write vector field values local to the quadrature point.
 * 
 * @tparam Tensor The type of tensor.
 * @param quad_index The index of the quadrature point.
 * @param q_value The field local to the quadrature point.
 * @param q_field The field at all quadrature points.
 */
template < typename Tensor,
           localIndex num_comp,
           std::enable_if_t<
             tensor::get_tensor_rank< Tensor > == 4 &&
             tensor::get_tensor_size< 3, Tensor > == num_comp,
             bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const (& q_value)[num_comp],
                  Tensor & q_field )
{
  #pragma unroll
  for (localIndex c = 0; c < num_comp; c++)
  {
    q_field( quad_index.x, quad_index.y, quad_index.z, c ) = q_value[c];
  }
}

/**
 * @brief Write rank 2 field values local to the quadrature point.
 * 
 * @tparam Tensor The type of tensor.
 * @param quad_index The index of the quadrature point.
 * @param q_value The field local to the quadrature point.
 * @param q_field The field at all quadrature points.
 */
template < typename Tensor,
           localIndex num_comp_x,
           localIndex num_comp_y,
           std::enable_if_t<
             tensor::get_tensor_rank< Tensor > == 5 &&
             tensor::get_tensor_size<3, Tensor> == num_comp_x &&
             tensor::get_tensor_size<4, Tensor> == num_comp_y,
             bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void qLocalWrite( TensorIndex const & quad_index,
                  real64 const (& q_value)[num_comp_x][num_comp_y],
                  Tensor & q_field )
{
  #pragma unroll
  for (localIndex c_x = 0; c_x < num_comp_x; c_x++)
  {
    #pragma unroll
    for (localIndex c_y = 0; c_y < num_comp_y; c_y++)
    {
      q_field( quad_index.x, quad_index.y, quad_index.z, c_x, c_y ) = q_value[c_x][c_y];
    }
  }
}

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_QUADRATUREFUNCTIONS_WRITE_HPP_ */
