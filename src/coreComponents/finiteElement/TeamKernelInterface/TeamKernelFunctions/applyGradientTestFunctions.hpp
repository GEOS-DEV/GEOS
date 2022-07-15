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
 * @file applyGradientTestFunctions.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_HPP_

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
void applyGradientTestFunctions( StackVariables & stack,
                                 real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                 real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                 real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][3],
                                 real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  SharedTensor< num_dofs_1d, num_quads_1d, num_quads_1d >
    Gqx( stack.shared_mem[0] ),
    Bqy( stack.shared_mem[1] ),
    Bqz( stack.shared_mem[2] );

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    {
      real64 gqx[num_quads_1d];
      real64 bqy[num_quads_1d];
      real64 bqz[num_quads_1d];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        gqx[quad_z] = 0.0;
        bqy[quad_z] = 0.0;
        bqz[quad_z] = 0.0;
      }
      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          // assumes quads in shared
          real64 const qx = q_values[quad_x][quad_y][quad_z][0];
          real64 const qy = q_values[quad_x][quad_y][quad_z][1];
          real64 const qz = q_values[quad_x][quad_y][quad_z][2];
          gqx[quad_z] += g * qx;
          bqy[quad_z] += b * qy;
          bqz[quad_z] += b * qz;
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        Gqx( dof_x, quad_y, quad_z ) = gqx[quad_z];
        Bqy( dof_x, quad_y, quad_z ) = bqy[quad_z];
        Bqz( dof_x, quad_y, quad_z ) = bqz[quad_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_dofs_1d, num_dofs_1d, num_quads_1d >
    BGqx( stack.shared_mem[3] ),
    GBqy( stack.shared_mem[4] ),
    BBqz( stack.shared_mem[5] );

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
    {
      real64 bgqx[num_quads_1d];
      real64 gbqy[num_quads_1d];
      real64 bbqz[num_quads_1d];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        bgqx[quad_z] = 0.0;
        gbqy[quad_z] = 0.0;
        bbqz[quad_z] = 0.0;
      }
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 const gqx = Gqx( dof_x, quad_y, quad_z );
          real64 const bqy = Bqy( dof_x, quad_y, quad_z );
          real64 const bqz = Bqz( dof_x, quad_y, quad_z );
          bgqx[quad_z] += b * gqx;
          gbqy[quad_z] += g * bqy;
          bbqz[quad_z] += b * bqz;
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        BGqx( dof_x, dof_y, quad_z ) = bgqx[quad_z];
        GBqy( dof_x, dof_y, quad_z ) = gbqy[quad_z];
        BBqz( dof_x, dof_y, quad_z ) = bbqz[quad_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    {
      // Cache values in registers to read them only once from shared
      real64 bgqx[num_quads_1d];
      real64 gbqy[num_quads_1d];
      real64 bbqz[num_quads_1d];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        bgqx[quad_z] = BGqx( dof_x, dof_y, quad_z );
        gbqy[quad_z] = GBqy( dof_x, dof_y, quad_z );
        bbqz[quad_z] = BBqz( dof_x, dof_y, quad_z );
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res = 0.0;
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          res += b * bgqx[quad_z] + b * gbqy[quad_z] + g * bbqz[quad_z];
        }
        dofs[dof_x][dof_y][dof_z] = res;
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
void applyGradientTestFunctions( StackVariables & stack,
                                 real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                 real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                 real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][3],
                                 real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  SharedTensor< num_dofs_1d, num_quads_1d, num_quads_1d, num_comp >
    Gqx( stack.shared_mem[0] ),
    Bqy( stack.shared_mem[1] ),
    Bqz( stack.shared_mem[2] );

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    {
      real64 gqx[num_quads_1d][num_comp];
      real64 bqy[num_quads_1d][num_comp];
      real64 bqz[num_quads_1d][num_comp];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          gqx[quad_z][comp] = 0.0;
          bqy[quad_z][comp] = 0.0;
          bqz[quad_z][comp] = 0.0;
        }
      }
      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            // assumes quads in shared
            real64 const qx = q_values[quad_x][quad_y][quad_z][comp][0];
            real64 const qy = q_values[quad_x][quad_y][quad_z][comp][1];
            real64 const qz = q_values[quad_x][quad_y][quad_z][comp][2];
            gqx[quad_z][comp] += g * qx;
            bqy[quad_z][comp] += b * qy;
            bqz[quad_z][comp] += b * qz;
          }
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          Gqx( dof_x, quad_y, quad_z, comp ) = gqx[quad_z][comp];
          Bqy( dof_x, quad_y, quad_z, comp ) = bqy[quad_z][comp];
          Bqz( dof_x, quad_y, quad_z, comp ) = bqz[quad_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_dofs_1d, num_dofs_1d, num_quads_1d, num_comp >
    BGqx( stack.shared_mem[3] ),
    GBqy( stack.shared_mem[4] ),
    BBqz( stack.shared_mem[5] );

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
    {
      real64 bgqx[num_quads_1d][num_comp];
      real64 gbqy[num_quads_1d][num_comp];
      real64 bbqz[num_quads_1d][num_comp];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          bgqx[quad_z][comp] = 0.0;
          gbqy[quad_z][comp] = 0.0;
          bbqz[quad_z][comp] = 0.0;
        }
      }
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            real64 const gqx = Gqx( dof_x, quad_y, quad_z, comp );
            real64 const bqy = Bqy( dof_x, quad_y, quad_z, comp );
            real64 const bqz = Bqz( dof_x, quad_y, quad_z, comp );
            bgqx[quad_z][comp] += b * gqx;
            gbqy[quad_z][comp] += g * bqy;
            bbqz[quad_z][comp] += b * bqz;
          }
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          BGqx( dof_x, dof_y, quad_z, comp ) = bgqx[quad_z][comp];
          GBqy( dof_x, dof_y, quad_z, comp ) = gbqy[quad_z][comp];
          BBqz( dof_x, dof_y, quad_z, comp ) = bbqz[quad_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    {
      // Cache values in registers to read them only once from shared
      real64 bgqx[num_quads_1d][num_comp];
      real64 gbqy[num_quads_1d][num_comp];
      real64 bbqz[num_quads_1d][num_comp];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          bgqx[quad_z][comp] = BGqx( dof_x, dof_y, quad_z, comp );
          gbqy[quad_z][comp] = GBqy( dof_x, dof_y, quad_z, comp );
          bbqz[quad_z][comp] = BBqz( dof_x, dof_y, quad_z, comp );
        }
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res[num_comp];
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[comp] = 0.0;
        }
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[comp] += b * bgqx[quad_z][comp]
                       + b * gbqy[quad_z][comp]
                       + g * bbqz[quad_z][comp];
          }
        }
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          dofs[dof_x][dof_y][dof_z][comp] = res[comp];
        }
      }
    });
  });

  ctx.teamSync();
}

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_HPP_ */
