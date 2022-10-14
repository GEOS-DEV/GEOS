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
 * @file applyTestFunctions.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_HPP_

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

/**
 * @brief "Apply" test functions.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] q_values Values of the "qFunction" at quadrature points.
 * @return Contribution of the q_values to the degrees of freedom.
 */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyTestFunctions( StackVariables & stack,
                         real64 const (& basis)[num_dofs_1d][num_quads_1d],
                         real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d],
                         real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;
  
  // Contraction on the first dimension
  SharedTensor< num_dofs_1d, num_quads_1d, num_quads_1d > Bu( stack.shared_mem[0] );

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    {
      real64 res[num_quads_1d];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        res[quad_z] = 0.0;
      }
      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] += b * q_values[quad_x][quad_y][quad_z]; // assumes quads in shared
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        Bu( dof_x, quad_y, quad_z ) = res[quad_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_dofs_1d, num_dofs_1d, num_quads_1d > BBu( stack.shared_mem[1] );

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
    {
      real64 res[num_quads_1d];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        res[quad_z] = 0.0;
      }
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] += b * Bu( dof_x, quad_y, quad_z );
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        BBu( dof_x, dof_y, quad_z ) = res[quad_z];
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
      real64 val[num_quads_1d];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        val[quad_z] = BBu( dof_x, dof_y, quad_z );
      }
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res = 0.0;
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res += basis[dof_z][quad_z] * val[quad_z];
        }
        dofs[dof_x][dof_y][dof_z] = res;
      }
    });
  });

  ctx.teamSync();
}

/**
 * @brief "Apply" test functions.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] q_values Values of the "qFunction" at quadrature points.
 * @return Contribution of the q_values to the degrees of freedom.
 */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyTestFunctions( StackVariables & stack,
                         real64 const (& basis)[num_dofs_1d][num_quads_1d],
                         real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp],
                         real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;
  
  // Contraction on the first dimension
  SharedTensor< num_dofs_1d, num_quads_1d, num_quads_1d, num_comp > Bu( stack.shared_mem[0] );

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
    {
      real64 res[num_quads_1d][num_comp];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[quad_z][comp] = 0.0;
        }
      }
      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[quad_z][comp] += b * q_values[quad_x][quad_y][quad_z][comp];
          }
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          Bu( dof_x, quad_y, quad_z, comp ) = res[quad_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  SharedTensor< num_dofs_1d, num_dofs_1d, num_quads_1d, num_comp > BBu( stack.shared_mem[1] );

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
    {
      real64 res[num_quads_1d][num_comp];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          res[quad_z][comp] = 0.0;
        }
      }
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[quad_z][comp] += b * Bu( dof_x, quad_y, quad_z, comp );
          }
        }
      }
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          BBu( dof_x, dof_y, quad_z, comp ) = res[quad_z][comp];
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
      real64 val[num_quads_1d][num_comp];
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          val[quad_z][comp] = BBu( dof_x, dof_y, quad_z, comp );
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
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            res[comp] += b * val[quad_z][comp];
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

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_HPP_ */
