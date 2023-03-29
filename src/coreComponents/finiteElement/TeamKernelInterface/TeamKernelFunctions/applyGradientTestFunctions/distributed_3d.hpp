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
 * @file applyGradientTestFunctions/distributed_3d.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_3D_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_3D_HPP_

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

// Non-distributed/Shared
template < localIndex... Sizes >
using SharedTensor = tensor::StaticPointerDTensor< Sizes... >;

///////////////////
// 3D distributed

template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyGradientTestFunctions(
  StackVariables & stack,
  real64 const (& basis)[num_dofs_1d][num_quads_1d],
  real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
  tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, 3 > const & q_values,
  tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > & dofs )
{
  RAJA::LaunchContext & ctx = stack.ctx;

  // Load in registers valuesfor (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_quads_1d], By[num_quads_1d], Bz[num_quads_1d];
  real64 Gx[num_quads_1d], Gy[num_quads_1d], Gz[num_quads_1d];
  localIndex dof_x = stack.tidx;
  if ( dof_x < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
    {
      Bx[quad_x] = basis[dof_x][quad_x];
      Gx[quad_x] = basis_gradient[dof_x][quad_x];
    }
  }
  localIndex dof_y = stack.tidy;
  if ( dof_y < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      By[quad_y] = basis[dof_y][quad_y];
      Gy[quad_y] = basis_gradient[dof_y][quad_y];
    }
  }
  localIndex dof_z = stack.tidz;
  if ( dof_z < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      Bz[quad_z] = basis[dof_z][quad_z];
      Gz[quad_z] = basis_gradient[dof_z][quad_z];
    }
  }

  SharedTensor< num_quads_1d, num_quads_1d, num_quads_1d, 3 > Du( stack.shared_mem[0] );
  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z ){
    #pragma unroll
    for (localIndex c = 0; c < 3; c++)
    {
      Du( quad_x, quad_y, quad_z, c ) = q_values( quad_x, quad_y, quad_z, c );
    }
  } );

  ctx.teamSync();

  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z ){
    real64 v = 0.0;
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      real64 const bz = Bz[ quad_z ];
      real64 const gz = Gz[ quad_z ];
      #pragma unroll
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const by = By[ quad_y ];
        real64 const gy = Gy[ quad_y ];
        #pragma unroll
        for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
        {
          real64 const bx = Bx[ quad_x ];
          real64 const gx = Gx[ quad_x ];
          real64 const valx = Du( quad_x, quad_y, quad_z, 0 ); // Avoid bank conflicts?
          real64 const valy = Du( quad_x, quad_y, quad_z, 1 );
          real64 const valz = Du( quad_x, quad_y, quad_z, 2 );
          v += gx * by * bz * valx;
          v += bx * gy * bz * valy;
          v += bx * by * gz * valz;
        }
      }
    }
    dofs( dof_x, dof_y, dof_z ) = v;
  } );
}

// Vector version
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyGradientTestFunctions(
  StackVariables & stack,
  real64 const (& basis)[num_dofs_1d][num_quads_1d],
  real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
  tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, 3 > const & q_values,
  tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & dofs )
{
  RAJA::LaunchContext & ctx = stack.ctx;

  // Load in registers values for (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_quads_1d], By[num_quads_1d], Bz[num_quads_1d];
  real64 Gx[num_quads_1d], Gy[num_quads_1d], Gz[num_quads_1d];
  localIndex dof_x = stack.tidx;
  if ( dof_x < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
    {
      Bx[quad_x] = basis[dof_x][quad_x];
      Gx[quad_x] = basis_gradient[dof_x][quad_x];
    }
  }
  localIndex dof_y = stack.tidy;
  if ( dof_y < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      By[quad_y] = basis[dof_y][quad_y];
      Gy[quad_y] = basis_gradient[dof_y][quad_y];
    }
  }
  localIndex dof_z = stack.tidz;
  if ( dof_z < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      Bz[quad_z] = basis[dof_z][quad_z];
      Gz[quad_z] = basis_gradient[dof_z][quad_z];
    }
  }

  SharedTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, 3 > Du( stack.shared_mem[0] );
  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z){
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      #pragma unroll
      for (localIndex dim = 0; dim < 3; dim++)
      {
        Du( quad_x, quad_y, quad_z, comp, dim ) = q_values( quad_x, quad_y, quad_z, comp, dim );
      }
    }
  } );

  ctx.teamSync();

  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z){
    real64 v[ num_comp ];
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      v[ comp ] = 0.0;
    }
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      real64 const bz = Bz[quad_z];
      real64 const gz = Gz[quad_z];
      #pragma unroll
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const by = By[quad_y];
        real64 const gy = Gy[quad_y];
        #pragma unroll
        for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
        {
          real64 const bx = Bx[quad_x];
          real64 const gx = Gx[quad_x];
          real64 const dx = gx * by * bz;
          real64 const dy = bx * gy * bz;
          real64 const dz = bx * by * gz;
          // localIndex const srcLane = quad_x + num_quads_1d * ( quad_y + num_quads_1d * ( quad_z + num_quads_1d * stack.batch_index ) );
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            real64 const valx = Du( quad_x, quad_y, quad_z, comp, 0 );
            // real64 const valx = __shfl_sync(0xffffffff, q_values( quad_x, quad_y, quad_z, comp, 0 ), srcLane);
            real64 const valy = Du( quad_x, quad_y, quad_z, comp, 1 );
            // real64 const valy = __shfl_sync(0xffffffff, q_values( quad_x, quad_y, quad_z, comp, 1 ), srcLane);
            real64 const valz = Du( quad_x, quad_y, quad_z, comp, 2 );
            // real64 const valz = __shfl_sync(0xffffffff, q_values( quad_x, quad_y, quad_z, comp, 2 ), srcLane);
            v[ comp ] += dx * valx;
            v[ comp ] += dy * valy;
            v[ comp ] += dz * valz;
          }
        }
      }
    }
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      dofs( dof_x, dof_y, dof_z, comp ) = v[ comp ];
    }
  } );
}

} // namespace impl

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_3D_HPP_ */
