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
 * @file interpolateAtQuadraturePoints.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_HPP_

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

template < localIndex... Sizes >
using SharedTensor = tensor::StaticPointerDTensor< Sizes... >;

// Non-distributed/2Dthreaded-Shared

/** @brief 3D interpolation operator using 2D block of threads with RAJA teams
 * 
 * @param[in] ctx RAJA team context.
 * @param[in] basis Two dimension array containing 1D shape functions evaluated
 *                  at quadrature points. Assumed to be in shared memory.
 * @param[in] dofs Three dimension array containing the dofs associated to the
 *                 shape functions. Assumed to be in shared memory.
 * @param[out] quads Three dimension array containing the values of the finite
 *                   element field evaluated at the quadrature points described
 *                   by the @a basis. Assumed to be in shared memory.
 * */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d>
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d],
                                    real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d > Bu( stack.shared_mem[0] );
  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      real64 res[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * dofs[dof_x][dof_y][dof_z];
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        Bu( quad_x, dof_y, dof_z ) = res[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d > BBu( stack.shared_mem[1] );

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
    {
      real64 res[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * Bu( quad_x, dof_y, dof_z );
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        BBu( quad_x, quad_y, dof_z ) = res[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        val[dof_z] = BBu( quad_x, quad_y, dof_z );
      }
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 res = 0.0;
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res += basis[dof_z][quad_z] * val[dof_z];
        }
        q_values[quad_x][quad_y][quad_z] = res;
      }
    });
  });

  ctx.teamSync();
}

/** @brief 3D interpolation operator using 2D block of threads with RAJA teams
 * 
 * @param[in] ctx RAJA team context.
 * @param[in] basis Two dimension array containing 1D shape functions evaluated
 *                  at quadrature points. Assumed to be in shared memory.
 * @param[in] dofs Three dimension array containing the dofs associated to the
 *                 shape functions. Assumed to be in shared memory.
 * @param[out] quads Three dimension array containing the values of the finite
 *                   element field evaluated at the quadrature points described
 *                   by the @a basis. Assumed to be in shared memory.
 * */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp],
                                    real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d, num_comp > Bu( stack.shared_mem[0] );

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      real64 res[num_dofs_1d][num_comp];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[dof_z][comp] = 0.0;
        }
      }
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] += b * dofs[dof_x][dof_y][dof_z][comp];
          }
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          Bu( quad_x, dof_y, dof_z, comp ) = res[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d, num_comp > BBu( stack.shared_mem[1] );

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
    {
      real64 res[num_dofs_1d][num_comp];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[dof_z][comp] = 0.0;
        }
      }
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] += b * Bu( quad_x, dof_y, dof_z, comp );
          }
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          BBu( quad_x, quad_y, dof_z, comp ) = res[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_dofs_1d][num_comp];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          val[dof_z][comp] = BBu( quad_x, quad_y, dof_z, comp );
        }
      }
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 res[num_comp];
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[comp] = 0.0;
        }
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[comp] += b * val[dof_z][comp];
          }
        }
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          q_values[quad_x][quad_y][quad_z][comp] = res[comp];
        }
      }
    });
  });

  ctx.teamSync();
}

// Stack tensor

/** @brief 3D interpolation operator using 2D block of threads with RAJA teams
 * 
 * @param[in] ctx RAJA team context.
 * @param[in] basis Two dimension array containing 1D shape functions evaluated
 *                  at quadrature points. Assumed to be in shared memory.
 * @param[in] dofs Three dimension array containing the dofs associated to the
 *                 shape functions. Assumed to be in shared memory.
 * @param[out] quads Three dimension array containing the values of the finite
 *                   element field evaluated at the quadrature points described
 *                   by the @a basis. Assumed to be in shared memory.
 * */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d>
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > const & dofs,
                                    tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_quads_1d > & q_values )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  tensor::StaticDTensor< num_quads_1d, num_dofs_1d, num_dofs_1d > Bu;
  #pragma unroll
  for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
  {
    #pragma unroll
    for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
    {
      real64 res[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * dofs( dof_x, dof_y, dof_z );
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        Bu( quad_x, dof_y, dof_z ) = res[dof_z];
      }
    }
  }

  // Contraction on the second dimension
  tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_dofs_1d > BBu;
  #pragma unroll
  for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
  {
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      real64 res[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * Bu( quad_x, dof_y, dof_z );
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        BBu( quad_x, quad_y, dof_z ) = res[dof_z];
      }
    }
  }

  // Contraction on the third dimension
  #pragma unroll
  for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
  {
    #pragma unroll
    for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        val[dof_z] = BBu( quad_x, quad_y, dof_z );
      }
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 res = 0.0;
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res += basis[dof_z][quad_z] * val[dof_z];
        }
        q_values( quad_x, quad_y, quad_z ) = res;
      }
    }
  }
}

/** @brief 3D interpolation operator using 2D block of threads with RAJA teams
 * 
 * @param[in] ctx RAJA team context.
 * @param[in] basis Two dimension array containing 1D shape functions evaluated
 *                  at quadrature points. Assumed to be in shared memory.
 * @param[in] dofs Three dimension array containing the dofs associated to the
 *                 shape functions. Assumed to be in shared memory.
 * @param[out] quads Three dimension array containing the values of the finite
 *                   element field evaluated at the quadrature points described
 *                   by the @a basis. Assumed to be in shared memory.
 * */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & dofs,
                                    tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp > & q_values )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  #pragma unroll
  for (localIndex c = 0; c < num_comp; c++)
  {
    // Contraction on the first dimension
    tensor::StaticDTensor< num_quads_1d, num_dofs_1d, num_dofs_1d > Bu;
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      #pragma unroll
      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 res[num_dofs_1d][num_comp];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] = 0.0;
          }
        }
        #pragma unroll
        for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
        {
          real64 const b = basis[dof_x][quad_x];
          #pragma unroll
          for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
          {
            #pragma unroll
            for (localIndex comp = 0; comp < num_comp; comp++)
            {
              res[dof_z][comp] += b * dofs( dof_x, dof_y, dof_z, comp );
            }
          }
        }
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            Bu( quad_x, dof_y, dof_z, comp ) = res[dof_z][comp];
          }
        }
      }
    }

    // Contraction on the second dimension
    tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_dofs_1d > BBu;
    #pragma unroll
    for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
    {
      #pragma unroll
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 res[num_dofs_1d][num_comp];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] = 0.0;
          }
        }
        #pragma unroll
        for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
        {
          real64 const b = basis[dof_y][quad_y];
          #pragma unroll
          for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
          {
            #pragma unroll
            for (localIndex comp = 0; comp < num_comp; comp++)
            {
              res[dof_z][comp] += b * Bu( quad_x, dof_y, dof_z, comp );
            }
          }
        }
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            BBu( quad_x, quad_y, dof_z, comp ) = res[dof_z][comp];
          }
        }
      }
    }

    // Contraction on the third dimension
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      #pragma unroll
      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        // Cache values in registers to read them only once from shared
        real64 val[num_dofs_1d][num_comp];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            val[dof_z][comp] = BBu( quad_x, quad_y, dof_z, comp );
          }
        }
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 res[num_comp];
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[comp] = 0.0;
          }
          #pragma unroll
          for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
          {
            real64 const b = basis[dof_z][quad_z];
            #pragma unroll
            for (localIndex comp = 0; comp < num_comp; comp++)
            {
              res[comp] += b * val[dof_z][comp];
            }
          }
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            q_values( quad_x, quad_y, quad_z, comp ) = res[comp];
          }
        }
      }
    }
  }
}

// 2D distributed

/** @brief 3D interpolation operator using 2D block of threads with RAJA teams
 * 
 * @param[in] ctx RAJA team context.
 * @param[in] basis Two dimension array containing 1D shape functions evaluated
 *                  at quadrature points. Assumed to be in shared memory.
 * @param[in] dofs Three dimension array containing the dofs associated to the
 *                 shape functions. Assumed to be in shared memory.
 * @param[out] quads Three dimension array containing the values of the finite
 *                   element field evaluated at the quadrature points described
 *                   by the @a basis. Assumed to be in shared memory.
 * */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d>
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    tensor::Static2dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > const & dofs,
                                    tensor::Static2dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d > & q_values )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > u( stack.shared_mem[3] );
  #pragma unroll
  for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
    {
      loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
      {
        u( dof_x, dof_y, dof_z ) = dofs( dof_x, dof_y, dof_z );
      } );
    } );
  }

  ctx.teamSync();

  // Contraction on the first dimension
  SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d > Bu( stack.shared_mem[0] );
  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      real64 res[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * u( dof_x, dof_y, dof_z );
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        Bu( quad_x, dof_y, dof_z ) = res[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d > BBu( stack.shared_mem[1] );

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
    {
      real64 res[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * Bu( quad_x, dof_y, dof_z );
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        BBu( quad_x, quad_y, dof_z ) = res[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_dofs_1d];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        val[dof_z] = BBu( quad_x, quad_y, dof_z );
      }
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 res = 0.0;
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res += basis[dof_z][quad_z] * val[dof_z];
        }
        q_values( quad_x, quad_y, quad_z ) = res;
      }
    });
  });

  // ctx.teamSync();
}

/** @brief 3D interpolation operator using 2D block of threads with RAJA teams
 * 
 * @param[in] ctx RAJA team context.
 * @param[in] basis Two dimension array containing 1D shape functions evaluated
 *                  at quadrature points. Assumed to be in shared memory.
 * @param[in] dofs Three dimension array containing the dofs associated to the
 *                 shape functions. Assumed to be in shared memory.
 * @param[out] quads Three dimension array containing the values of the finite
 *                   element field evaluated at the quadrature points described
 *                   by the @a basis. Assumed to be in shared memory.
 * */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    tensor::Static2dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & dofs,
                                    tensor::Static2dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp > & q_values )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > u( stack.shared_mem[3] );
  #pragma unroll
  for (localIndex c = 0; c < num_comp; c++)
  {
    #pragma unroll
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
      {
        loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
        {
          u( dof_x, dof_y, dof_z, c ) = dofs( dof_x, dof_y, dof_z, c );
        } );
      } );
    }
  }

  ctx.teamSync();

  // Contraction on the first dimension
  SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d, num_comp > Bu( stack.shared_mem[0] );

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      real64 res[num_dofs_1d][num_comp];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[dof_z][comp] = 0.0;
        }
      }
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] += b * u( dof_x, dof_y, dof_z, comp );
          }
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          Bu( quad_x, dof_y, dof_z, comp ) = res[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d, num_comp > BBu( stack.shared_mem[1] );

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
    {
      real64 res[num_dofs_1d][num_comp];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[dof_z][comp] = 0.0;
        }
      }
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] += b * Bu( quad_x, dof_y, dof_z, comp );
          }
        }
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          BBu( quad_x, quad_y, dof_z, comp ) = res[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_dofs_1d][num_comp];
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          val[dof_z][comp] = BBu( quad_x, quad_y, dof_z, comp );
        }
      }
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 res[num_comp];
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[comp] = 0.0;
        }
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[comp] += b * val[dof_z][comp];
          }
        }
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          q_values( quad_x, quad_y, quad_z, comp ) = res[comp];
        }
      }
    });
  });

  // ctx.teamSync();
}

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
  LaunchContext & ctx = stack.ctx;

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
  LaunchContext & ctx = stack.ctx;

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

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_HPP_ */
