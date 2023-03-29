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
 * @file interpolateAtQuadraturePoints/distributed_3d.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_3D_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_3D_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_types.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

namespace impl
{

template < localIndex... Sizes >
using SharedTensor = tensor::StaticPointerDTensor< Sizes... >;

// 3D Threaded version using RAJA teams
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints(
  StackVariables & stack,
  real64 const (& basis)[num_dofs_1d][num_quads_1d],
  tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > const & dofs,
  tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d > & q_values )
{
  RAJA::LaunchContext & ctx = stack.ctx;

  // Load in registers values for (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_dofs_1d], By[num_dofs_1d], Bz[num_dofs_1d];
  localIndex quad_x = stack.tidx;
  if ( quad_x < num_quads_1d )
  {
    #pragma unroll
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      Bx[dof_x] = basis[dof_x][quad_x];
    }
  }
  localIndex quad_y = stack.tidy;
  if ( quad_y < num_quads_1d )
  {
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      By[dof_y] = basis[dof_y][quad_y];
    }
  }
  localIndex quad_z = stack.tidz;
  if ( quad_z < num_quads_1d )
  {
    #pragma unroll
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      Bz[dof_z] = basis[dof_z][quad_z];
    }
  }

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > u( stack.shared_mem[3] );
  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z){
    u( dof_x, dof_y, dof_z ) = dofs( dof_x, dof_y, dof_z );
  } );

  ctx.teamSync();

  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z){
    real64 bbbu = 0.0;
    #pragma unroll
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      real64 const bz = Bz[dof_z];
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const by = By[dof_y];
        #pragma unroll
        for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
        {
          real64 const bx = Bx[dof_x];
          real64 const val = u( dof_x, dof_y, dof_z );
          bbbu += bx * by * bz * val;
        }
      }
    }
    q_values( quad_x, quad_y, quad_z ) = bbbu;
  } );
}

// 3D threaded vector case
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints(
  StackVariables & stack,
  real64 const (& basis)[num_dofs_1d][num_quads_1d],
  tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & dofs,
  tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp > & q_values )
{
  RAJA::LaunchContext & ctx = stack.ctx;

  // Load in registers values for (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_dofs_1d], By[num_dofs_1d], Bz[num_dofs_1d];
  localIndex quad_x = stack.tidx;
  if ( quad_x < num_quads_1d )
  {
    #pragma unroll
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      Bx[dof_x] = basis[dof_x][quad_x];
    }
  }
  localIndex quad_y = stack.tidy;
  if ( quad_y < num_quads_1d )
  {
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      By[dof_y] = basis[dof_y][quad_y];
    }
  }
  localIndex quad_z = stack.tidz;
  if ( quad_z < num_quads_1d )
  {
    #pragma unroll
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      Bz[dof_z] = basis[dof_z][quad_z];
    }
  }

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > u( stack.shared_mem[3] );
  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z){
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      u( dof_x, dof_y, dof_z, comp ) = dofs( dof_x, dof_y, dof_z, comp );
    }
  } );

  ctx.teamSync();

  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z){
    real64 bbbu[ num_comp ];
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      bbbu[ comp ] = 0.0;
    }
    #pragma unroll
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      real64 const bz = Bz[dof_z];
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const by = By[dof_y];
        #pragma unroll
        for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
        {
          real64 const bx = Bx[dof_x];
          real64 const b = bx * by * bz;
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            real64 const val = u( dof_x, dof_y, dof_z, comp );
            bbbu[ comp ] += b * val;
          }
        }
      }
    }
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      q_values( quad_x, quad_y, quad_z, comp ) = bbbu[ comp ];
    }
  } );
}

} // namespace impl

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_3D_HPP_ */
