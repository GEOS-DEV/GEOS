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
 * @file writeAddField.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_WRITE_ADD_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_WRITE_ADD_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_types.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/common.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

// Non-distributed/Shared version
template < typename StackVariables,
           typename Field,
           localIndex stride_x, localIndex stride_y, localIndex stride_z >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void writeAddField( StackVariables & stack,
                    real64 const (& local_field)[stride_x][stride_y][stride_z],
                    Field & field )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  loop<thread_x> (ctx, RangeSegment(0, stride_x), [&] (localIndex ind_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, stride_y), [&] (localIndex ind_y)
    {
      for (localIndex ind_z = 0; ind_z < stride_z; ind_z++)
      {
        localIndex const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        localIndex const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        RAJA::atomicAdd( RAJA::auto_atomic{},
                         &field[ global_node_index ],
                         local_field[ ind_x ][ ind_y ][ ind_z ]);
      }
    });
  });
}

template < typename StackVariables,
           typename Field,
           localIndex stride_x, localIndex stride_y, localIndex stride_z, localIndex dim >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void writeAddField( StackVariables & stack,
                    real64 const (& local_field)[stride_x][stride_y][stride_z][dim], 
                    Field & field )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  loop<thread_x> (ctx, RangeSegment(0, stride_x), [&] (localIndex ind_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, stride_y), [&] (localIndex ind_y)
    {
      for (localIndex ind_z = 0; ind_z < stride_z; ind_z++)
      {
        localIndex const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        localIndex const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        for (localIndex d = 0; d < dim; d++)
        {
          RAJA::atomicAdd( RAJA::auto_atomic{},
                           &field( global_node_index, d ),
                           local_field[ ind_x ][ ind_y ][ ind_z ][ d ] );
        }
      }
    });
  });
}

// Generic version
template < typename StackVariables,
           typename EtoNMap,
           typename Field,
           typename Tensor,
           std::enable_if_t< tensor::get_tensor_rank< Tensor > == 3, bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void writeAddField( StackVariables & stack,
                    EtoNMap const & m_elemsToNodes,
                    Tensor const & local_field,
                    Field & field )
{
  constexpr localIndex stride_x = tensor::get_tensor_size<0, Tensor>;
  constexpr localIndex stride_y = tensor::get_tensor_size<1, Tensor>;
  constexpr localIndex stride_z = tensor::get_tensor_size<2, Tensor>;

  loop3D( stack, stride_x, stride_y, stride_z,
          [&]( localIndex ind_x, localIndex ind_y, localIndex ind_z){
    localIndex const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
    localIndex const global_node_index = m_elemsToNodes( stack.element_index, local_node_index );
    RAJA::atomicAdd( RAJA::auto_atomic{},
                     &field[ global_node_index ],
                     local_field( ind_x, ind_y, ind_z ) );
  } );
}

template < typename StackVariables,
           typename EtoNMap,
           typename Field,
           typename Tensor,
           std::enable_if_t< tensor::get_tensor_rank< Tensor > == 4, bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void writeAddField( StackVariables & stack,
                    EtoNMap const & m_elemsToNodes,
                    Tensor const & local_field,
                    Field & field )
{
  constexpr localIndex stride_x = tensor::get_tensor_size<0, Tensor>;
  constexpr localIndex stride_y = tensor::get_tensor_size<1, Tensor>;
  constexpr localIndex stride_z = tensor::get_tensor_size<2, Tensor>;
  constexpr localIndex dim = tensor::get_tensor_size<3, Tensor>;

  loop3D( stack, stride_x, stride_y, stride_z,
          [&]( localIndex ind_x, localIndex ind_y, localIndex ind_z){
    localIndex const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
    localIndex const global_node_index = m_elemsToNodes( stack.element_index, local_node_index );
    #pragma unroll
    for (localIndex d = 0; d < dim; d++)
    {
      RAJA::atomicAdd( RAJA::auto_atomic{},
                       &field( global_node_index, d ),
                       local_field( ind_x, ind_y, ind_z, d ) );
    }
  } );
}

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_WRITE_ADD_HPP_ */
