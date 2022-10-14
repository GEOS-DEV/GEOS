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

// 3D Threaded version using RAJA teams
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

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_GRAD_HPP_ */
