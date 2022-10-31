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
 * @file interpolateGradientAtQuadraturePoints.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_GRAD_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_GRAD_HPP_

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
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d],
                                            real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][3] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d >
    Bu( stack.shared_mem[0] ),
    Gu( stack.shared_mem[1] );

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      real64 bu[num_dofs_1d];
      real64 gu[num_dofs_1d];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bu[dof_z] = 0.0;
        gu[dof_z] = 0.0;
      }
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const val = dofs[dof_x][dof_y][dof_z];
          bu[dof_z] += b * val; // assumes dofs in shared
          gu[dof_z] += g * val; // assumes dofs in shared
        }
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        Bu( quad_x, dof_y, dof_z ) = bu[dof_z];
        Gu( quad_x, dof_y, dof_z ) = gu[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d >
    BBu( stack.shared_mem[2] ),
    BGu( stack.shared_mem[3] ),
    GBu( stack.shared_mem[4] );

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
    {
      real64 bbu[num_dofs_1d];
      real64 bgu[num_dofs_1d];
      real64 gbu[num_dofs_1d];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bbu[dof_z] = 0.0;
        bgu[dof_z] = 0.0;
        gbu[dof_z] = 0.0;
      }
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const bu = Bu( quad_x, dof_y, dof_z );
          real64 const gu = Gu( quad_x, dof_y, dof_z );
          bbu[dof_z] += b * bu;
          bgu[dof_z] += b * gu;
          gbu[dof_z] += g * bu;
        }
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        BBu( quad_x, quad_y, dof_z ) = bbu[dof_z];
        BGu( quad_x, quad_y, dof_z ) = bgu[dof_z];
        GBu( quad_x, quad_y, dof_z ) = gbu[dof_z];
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
      real64 bbu[num_dofs_1d];
      real64 bgu[num_dofs_1d];
      real64 gbu[num_dofs_1d];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bbu[dof_z] = BBu( quad_x, quad_y, dof_z );
        bgu[dof_z] = BGu( quad_x, quad_y, dof_z );
        gbu[dof_z] = GBu( quad_x, quad_y, dof_z );
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 bbgu = 0.0;
        real64 bgbu = 0.0;
        real64 gbbu = 0.0;
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          bbgu += b * bgu[dof_z];
          bgbu += b * gbu[dof_z];
          gbbu += g * bbu[dof_z];
        }
        q_values[quad_x][quad_y][quad_z][0] = bbgu;
        q_values[quad_x][quad_y][quad_z][1] = bgbu;
        q_values[quad_x][quad_y][quad_z][2] = gbbu;
      }
    });
  });

  ctx.teamSync();
}

// 3D Threaded version using RAJA teams
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp],
                                            real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][3] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d, num_comp >
    Bu( stack.shared_mem[0] ),
    Gu( stack.shared_mem[1] );

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      real64 bu[num_dofs_1d][num_comp];
      real64 gu[num_dofs_1d][num_comp];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          bu[dof_z][comp] = 0.0;
          gu[dof_z][comp] = 0.0;
        }
      }
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            real64 const val = dofs[dof_x][dof_y][dof_z][comp];
            bu[dof_z][comp] += b * val;
            gu[dof_z][comp] += g * val;
          }
        }
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          Bu( quad_x, dof_y, dof_z, comp ) = bu[dof_z][comp];
          Gu( quad_x, dof_y, dof_z, comp ) = gu[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d, num_comp >
    BBu( stack.shared_mem[2] ),
    BGu( stack.shared_mem[3] ),
    GBu( stack.shared_mem[4] );

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
    {
      real64 bbu[num_dofs_1d][num_comp];
      real64 bgu[num_dofs_1d][num_comp];
      real64 gbu[num_dofs_1d][num_comp];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          bbu[dof_z][comp] = 0.0;
          bgu[dof_z][comp] = 0.0;
          gbu[dof_z][comp] = 0.0;
        }
      }
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            real64 const bu = Bu( quad_x, dof_y, dof_z, comp );
            real64 const gu = Gu( quad_x, dof_y, dof_z, comp );
            bbu[dof_z][comp] += b * bu;
            bgu[dof_z][comp] += b * gu;
            gbu[dof_z][comp] += g * bu;
          }
        }
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          BBu( quad_x, quad_y, dof_z, comp ) = bbu[dof_z][comp];
          BGu( quad_x, quad_y, dof_z, comp ) = bgu[dof_z][comp];
          GBu( quad_x, quad_y, dof_z, comp ) = gbu[dof_z][comp];
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
      real64 bbu[num_dofs_1d][num_comp];
      real64 bgu[num_dofs_1d][num_comp];
      real64 gbu[num_dofs_1d][num_comp];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          bbu[dof_z][comp] = BBu( quad_x, quad_y, dof_z, comp );
          bgu[dof_z][comp] = BGu( quad_x, quad_y, dof_z, comp );
          gbu[dof_z][comp] = GBu( quad_x, quad_y, dof_z, comp );
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 bbgu[num_comp];
        real64 bgbu[num_comp];
        real64 gbbu[num_comp];
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          bbgu[comp] = 0.0;
          bgbu[comp] = 0.0;
          gbbu[comp] = 0.0;
        }
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            bbgu[comp] += b * bgu[dof_z][comp];
            bgbu[comp] += b * gbu[dof_z][comp];
            gbbu[comp] += g * bbu[dof_z][comp];
          }
        }
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          q_values[quad_x][quad_y][quad_z][comp][0] = bbgu[comp];
          q_values[quad_x][quad_y][quad_z][comp][1] = bgbu[comp];
          q_values[quad_x][quad_y][quad_z][comp][2] = gbbu[comp];
        }
      }
    });
  });

  ctx.teamSync();
}

// 2Dthreaded
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            tensor::Static2dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > const & dofs,
                                            tensor::Static2dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, 3 > & q_values )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > u( stack.shared_mem[3] );
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
  SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d >
    Bu( stack.shared_mem[0] ),
    Gu( stack.shared_mem[1] );

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      real64 bu[num_dofs_1d];
      real64 gu[num_dofs_1d];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bu[dof_z] = 0.0;
        gu[dof_z] = 0.0;
      }
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const val = u( dof_x, dof_y, dof_z );
          bu[dof_z] += b * val; // assumes dofs in shared
          gu[dof_z] += g * val; // assumes dofs in shared
        }
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        Bu( quad_x, dof_y, dof_z ) = bu[dof_z];
        Gu( quad_x, dof_y, dof_z ) = gu[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d >
    BBu( stack.shared_mem[2] ),
    BGu( stack.shared_mem[3] ),
    GBu( stack.shared_mem[4] );

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
    {
      real64 bbu[num_dofs_1d];
      real64 bgu[num_dofs_1d];
      real64 gbu[num_dofs_1d];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bbu[dof_z] = 0.0;
        bgu[dof_z] = 0.0;
        gbu[dof_z] = 0.0;
      }
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const bu = Bu( quad_x, dof_y, dof_z );
          real64 const gu = Gu( quad_x, dof_y, dof_z );
          bbu[dof_z] += b * bu;
          bgu[dof_z] += b * gu;
          gbu[dof_z] += g * bu;
        }
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        BBu( quad_x, quad_y, dof_z ) = bbu[dof_z];
        BGu( quad_x, quad_y, dof_z ) = bgu[dof_z];
        GBu( quad_x, quad_y, dof_z ) = gbu[dof_z];
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
      real64 bbu[num_dofs_1d];
      real64 bgu[num_dofs_1d];
      real64 gbu[num_dofs_1d];
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bbu[dof_z] = BBu( quad_x, quad_y, dof_z );
        bgu[dof_z] = BGu( quad_x, quad_y, dof_z );
        gbu[dof_z] = GBu( quad_x, quad_y, dof_z );
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 bbgu = 0.0;
        real64 bgbu = 0.0;
        real64 gbbu = 0.0;
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          bbgu += b * bgu[dof_z];
          bgbu += b * gbu[dof_z];
          gbbu += g * bbu[dof_z];
        }
        q_values( quad_x, quad_y, quad_z, 0 ) = bbgu;
        q_values( quad_x, quad_y, quad_z, 1 ) = bgbu;
        q_values( quad_x, quad_y, quad_z, 2 ) = gbbu;
      }
    });
  });

  // ctx.teamSync();
}

// 2D Threaded version using RAJA teams
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            tensor::Static2dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & dofs,
                                            tensor::Static2dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, 3 > & q_values )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > u( stack.shared_mem[3] );
  for (localIndex comp = 0; comp < num_comp; comp++)
  {
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
      {
        loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
        {
          u( dof_x, dof_y, dof_z ) = dofs( dof_x, dof_y, dof_z, comp );
        } );
      } );
    }

    ctx.teamSync();

    // Contraction on the first dimension
    SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d >
      Bu( stack.shared_mem[0] ),
      Gu( stack.shared_mem[1] );

    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
    {
      loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
      {
        real64 bu[num_dofs_1d];
        real64 gu[num_dofs_1d];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          bu[dof_z] = 0.0;
          gu[dof_z] = 0.0;
        }
        for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
        {
          real64 const b = basis[dof_x][quad_x];
          real64 const g = basis_gradient[dof_x][quad_x];
          for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
          {
            real64 const val = u( dof_x, dof_y, dof_z );
            bu[dof_z] += b * val;
            gu[dof_z] += g * val;
          }
        }
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
            Bu( quad_x, dof_y, dof_z ) = bu[dof_z];
            Gu( quad_x, dof_y, dof_z ) = gu[dof_z];
        }
      });
    });

    ctx.teamSync();

    // Contraction on the second dimension
    SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d >
      BBu( stack.shared_mem[2] ),
      BGu( stack.shared_mem[3] ),
      GBu( stack.shared_mem[4] );

    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
    {
      loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
      {
        real64 bbu[num_dofs_1d];
        real64 bgu[num_dofs_1d];
        real64 gbu[num_dofs_1d];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          bbu[dof_z] = 0.0;
          bgu[dof_z] = 0.0;
          gbu[dof_z] = 0.0;
        }
        for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
        {
          real64 const b = basis[dof_y][quad_y];
          real64 const g = basis_gradient[dof_y][quad_y];
          for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
          {
            real64 const bu = Bu( quad_x, dof_y, dof_z );
            real64 const gu = Gu( quad_x, dof_y, dof_z );
            bbu[dof_z] += b * bu;
            bgu[dof_z] += b * gu;
            gbu[dof_z] += g * bu;
          }
        }
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
            BBu( quad_x, quad_y, dof_z ) = bbu[dof_z];
            BGu( quad_x, quad_y, dof_z ) = bgu[dof_z];
            GBu( quad_x, quad_y, dof_z ) = gbu[dof_z];
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
        real64 bbu[num_dofs_1d];
        real64 bgu[num_dofs_1d];
        real64 gbu[num_dofs_1d];
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
            bbu[dof_z] = BBu( quad_x, quad_y, dof_z );
            bgu[dof_z] = BGu( quad_x, quad_y, dof_z );
            gbu[dof_z] = GBu( quad_x, quad_y, dof_z );
        }
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 bbgu = 0.0;
          real64 bgbu = 0.0;
          real64 gbbu = 0.0; 
          for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
          {
            real64 const b = basis[dof_z][quad_z];
            real64 const g = basis_gradient[dof_z][quad_z];
            bbgu += b * bgu[dof_z];
            bgbu += b * gbu[dof_z];
            gbbu += g * bbu[dof_z];
          }
          q_values( quad_x, quad_y, quad_z, comp, 0 ) = bbgu;
          q_values( quad_x, quad_y, quad_z, comp, 1 ) = bgbu;
          q_values( quad_x, quad_y, quad_z, comp, 2 ) = gbbu;
        }
      });
    });

    ctx.teamSync();
  }
}

// 3D Threaded version using RAJA teams
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > const & dofs,
                                            tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, 3 > & q_values )
{
  LaunchContext & ctx = stack.ctx;

  // Load in registers values for (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_dofs_1d], By[num_dofs_1d], Bz[num_dofs_1d];
  real64 Gx[num_dofs_1d], Gy[num_dofs_1d], Gz[num_dofs_1d];
  localIndex quad_x = stack.tidx;
  if ( quad_x < num_quads_1d )
  {
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      Bx[dof_x] = basis[dof_x][quad_x];
      Gx[dof_x] = basis_gradient[dof_x][quad_x];
    }
  }
  localIndex quad_y = stack.tidy;
  if ( quad_y < num_quads_1d )
  {
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      By[dof_y] = basis[dof_y][quad_y];
      Gy[dof_y] = basis_gradient[dof_y][quad_y];
    }
  }
  localIndex quad_z = stack.tidz;
  if ( quad_z < num_quads_1d )
  {
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      Bz[dof_z] = basis[dof_z][quad_z];
      Gz[dof_z] = basis_gradient[dof_z][quad_z];
    }
  }

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > u( stack.shared_mem[0] );
  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z){
    u( dof_x, dof_y, dof_z ) = dofs( dof_x, dof_y, dof_z );
  } );

  ctx.teamSync();

  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z){
    real64 bbgu = 0.0;
    real64 bgbu = 0.0;
    real64 gbbu = 0.0;
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      real64 const bz = Bz[dof_z];
      real64 const gz = Gz[dof_z];
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const by = By[dof_y];
        real64 const gy = Gy[dof_y];
        for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
        {
          real64 const bx = Bx[dof_x];
          real64 const gx = Gx[dof_x];
          real64 const val = u( dof_x, dof_y, dof_z );
          bbgu += gx * by * bz * val;
          bgbu += bx * gy * bz * val;
          gbbu += bx * by * gz * val;
        }
      }
    }
    q_values( quad_x, quad_y, quad_z, 0 ) = bbgu;
    q_values( quad_x, quad_y, quad_z, 1 ) = bgbu;
    q_values( quad_x, quad_y, quad_z, 2 ) = gbbu;
  } );
}

// 3D threaded vector case
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & dofs,
                                            tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, 3 > & q_values )
{
  LaunchContext & ctx = stack.ctx;

  // Load in registers values for (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_dofs_1d], By[num_dofs_1d], Bz[num_dofs_1d];
  real64 Gx[num_dofs_1d], Gy[num_dofs_1d], Gz[num_dofs_1d];
  localIndex quad_x = stack.tidx;
  if ( quad_x < num_quads_1d )
  {
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      Bx[dof_x] = basis[dof_x][quad_x];
      Gx[dof_x] = basis_gradient[dof_x][quad_x];
    }
  }
  localIndex quad_y = stack.tidy;
  if ( quad_y < num_quads_1d )
  {
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      By[dof_y] = basis[dof_y][quad_y];
      Gy[dof_y] = basis_gradient[dof_y][quad_y];
    }
  }
  localIndex quad_z = stack.tidz;
  if ( quad_z < num_quads_1d )
  {
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      Bz[dof_z] = basis[dof_z][quad_z];
      Gz[dof_z] = basis_gradient[dof_z][quad_z];
    }
  }

  SharedTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > u( stack.shared_mem[0] );
  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z){
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      u( dof_x, dof_y, dof_z, comp ) = dofs( dof_x, dof_y, dof_z, comp );
    }
  } );

  ctx.teamSync();

  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z){
    real64 bbgu[ num_comp ];
    real64 bgbu[ num_comp ];
    real64 gbbu[ num_comp ];
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      bbgu[ comp ] = 0.0;
      bgbu[ comp ] = 0.0;
      gbbu[ comp ] = 0.0;
    }
    for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    {
      real64 const bz = Bz[dof_z];
      real64 const gz = Gz[dof_z];
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const by = By[dof_y];
        real64 const gy = Gy[dof_y];
        for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
        {
          real64 const bx = Bx[dof_x];
          real64 const gx = Gx[dof_x];
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            real64 const val = u( dof_x, dof_y, dof_z, comp );
            bbgu[ comp ] += gx * by * bz * val;
            bgbu[ comp ] += bx * gy * bz * val;
            gbbu[ comp ] += bx * by * gz * val;
          }
        }
      }
    }
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      q_values( quad_x, quad_y, quad_z, comp, 0 ) = bbgu[ comp ];
      q_values( quad_x, quad_y, quad_z, comp, 1 ) = bgbu[ comp ];
      q_values( quad_x, quad_y, quad_z, comp, 2 ) = gbbu[ comp ];
    }
  } );
}

// 3D Threaded version using RAJA teams
// template < typename StackVariables,
//            localIndex num_comp >
// GEOSX_HOST_DEVICE
// GEOSX_FORCE_INLINE
// void interpolateGradientAtQuadraturePoints( StackVariables & stack,
//                                             real64 const (& basis)[2][2],
//                                             real64 const (& basis_gradient)[2][2],
//                                             real64 const (& dofs)[2][2][2][num_comp],
//                                             real64 (& q_values)[2][2][2][num_comp][3] )
// {
//   using RAJA::RangeSegment;
//   LaunchContext & ctx = stack.ctx;
//   constexpr localIndex num_dofs_1d = 2;
//   constexpr localIndex num_quads_1d = 2;

//   for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
//   {
//     loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
//     {
//       loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
//       {
//         real64 bbgu = 0.0;
//         real64 bgbu = 0.0;
//         real64 gbbu = 0.0;
//         for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//         {
//           real64 const bz = basis[dof_z][quad_z];
//           real64 const gz = basis_gradient[dof_z][quad_z];
//           for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
//           {
//             real64 const by = basis[dof_y][quad_y];
//             real64 const gy = basis_gradient[dof_y][quad_y];
//             for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
//             {
//               real64 const bx = basis[dof_x][quad_x];
//               real64 const gx = basis_gradient[dof_x][quad_x];
//               for (localIndex comp = 0; comp < num_comp; comp++)
//               {
//                 real64 const val = dofs[dof_x][dof_y][dof_z][comp];
//                 bbgu = gx * by * bz * val;
//                 bgbu = bx * gy * bz * val;
//                 gbbu = bx * by * gz * val;
//               }
//             }
//           }
//         }
//         q_values[quad_x][quad_y][quad_z][0] = bbgu;
//         q_values[quad_x][quad_y][quad_z][1] = bgbu;
//         q_values[quad_x][quad_y][quad_z][2] = gbbu;
//       });
//     });
//   }

//   // Contraction on the first dimension
//   SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d, num_comp >
//     Bu( stack.shared_mem[0] ),
//     Gu( stack.shared_mem[1] );

//   loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
//   {
//     loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
//     {
//       real64 bu[num_dofs_1d][num_comp];
//       real64 gu[num_dofs_1d][num_comp];
//       for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//       {
//         for (localIndex comp = 0; comp < num_comp; comp++)
//         {
//           bu[dof_z][comp] = 0.0;
//           gu[dof_z][comp] = 0.0;
//         }
//       }
//       for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
//       {
//         real64 const b = basis[dof_x][quad_x];
//         real64 const g = basis_gradient[dof_x][quad_x];
//         for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//         {
//           for (localIndex comp = 0; comp < num_comp; comp++)
//           {
//             real64 const val = dofs[dof_x][dof_y][dof_z][comp];
//             bu[dof_z][comp] += b * val;
//             gu[dof_z][comp] += g * val;
//           }
//         }
//       }
//       for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//       {
//         for (localIndex comp = 0; comp < num_comp; comp++)
//         {
//           Bu( quad_x, dof_y, dof_z, comp ) = bu[dof_z][comp];
//           Gu( quad_x, dof_y, dof_z, comp ) = gu[dof_z][comp];
//         }
//       }
//     });
//   });

//   ctx.teamSync();

//   // Contraction on the second dimension
//   SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d, num_comp >
//     BBu( stack.shared_mem[2] ),
//     BGu( stack.shared_mem[3] ),
//     GBu( stack.shared_mem[4] );

//   loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
//   {
//     loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
//     {
//       real64 bbu[num_dofs_1d][num_comp];
//       real64 bgu[num_dofs_1d][num_comp];
//       real64 gbu[num_dofs_1d][num_comp];
//       for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//       {
//         for (localIndex comp = 0; comp < num_comp; comp++)
//         {
//           bbu[dof_z][comp] = 0.0;
//           bgu[dof_z][comp] = 0.0;
//           gbu[dof_z][comp] = 0.0;
//         }
//       }
//       for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
//       {
//         real64 const b = basis[dof_y][quad_y];
//         real64 const g = basis_gradient[dof_y][quad_y];
//         for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//         {
//           for (localIndex comp = 0; comp < num_comp; comp++)
//           {
//             real64 const bu = Bu( quad_x, dof_y, dof_z, comp );
//             real64 const gu = Gu( quad_x, dof_y, dof_z, comp );
//             bbu[dof_z][comp] += b * bu;
//             bgu[dof_z][comp] += b * gu;
//             gbu[dof_z][comp] += g * bu;
//           }
//         }
//       }
//       for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//       {
//         for (localIndex comp = 0; comp < num_comp; comp++)
//         {
//           BBu( quad_x, quad_y, dof_z, comp ) = bbu[dof_z][comp];
//           BGu( quad_x, quad_y, dof_z, comp ) = bgu[dof_z][comp];
//           GBu( quad_x, quad_y, dof_z, comp ) = gbu[dof_z][comp];
//         }
//       }
//     });
//   });

//   ctx.teamSync();

//   // Contraction on the third dimension
//   loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
//   {
//     loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
//     {
//       // Cache values in registers to read them only once from shared
//       real64 bbu[num_dofs_1d][num_comp];
//       real64 bgu[num_dofs_1d][num_comp];
//       real64 gbu[num_dofs_1d][num_comp];
//       for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//       {
//         for (localIndex comp = 0; comp < num_comp; comp++)
//         {
//           bbu[dof_z][comp] = BBu( quad_x, quad_y, dof_z, comp );
//           bgu[dof_z][comp] = BGu( quad_x, quad_y, dof_z, comp );
//           gbu[dof_z][comp] = GBu( quad_x, quad_y, dof_z, comp );
//         }
//       }
//       for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
//       {
//         real64 bbgu[num_comp];
//         real64 bgbu[num_comp];
//         real64 gbbu[num_comp];
//         for (localIndex comp = 0; comp < num_comp; comp++)
//         {
//           bbgu[comp] = 0.0;
//           bgbu[comp] = 0.0;
//           gbbu[comp] = 0.0;
//         }
//         for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
//         {
//           real64 const b = basis[dof_z][quad_z];
//           real64 const g = basis_gradient[dof_z][quad_z];
//           for (localIndex comp = 0; comp < num_comp; comp++)
//           {
//             bbgu[comp] += b * bgu[dof_z][comp];
//             bgbu[comp] += b * gbu[dof_z][comp];
//             gbbu[comp] += g * bbu[dof_z][comp];
//           }
//         }
//         for (localIndex comp = 0; comp < num_comp; comp++)
//         {
//           q_values[quad_x][quad_y][quad_z][comp][0] = bbgu[comp];
//           q_values[quad_x][quad_y][quad_z][comp][1] = bgbu[comp];
//           q_values[quad_x][quad_y][quad_z][comp][2] = gbbu[comp];
//         }
//       }
//     });
//   });

//   ctx.teamSync();
// }

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_GRAD_HPP_ */
