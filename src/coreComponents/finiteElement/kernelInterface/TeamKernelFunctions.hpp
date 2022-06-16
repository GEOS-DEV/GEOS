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
 * @file TeamKernelFunction.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

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
           size_t num_dofs_1d,
           size_t num_quads_1d>
GEOSX_HOST_DEVICE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d],
                                    real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Bu[num_quads_1d][num_dofs_1d][num_dofs_1d];

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      real64 res[num_dofs_1d];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      for (size_t dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * dofs[dof_x][dof_y][dof_z];
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        Bu[quad_x][dof_y][dof_z] = res[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BBu[num_quads_1d][num_quads_1d][num_dofs_1d];

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
    {
      real64 res[num_dofs_1d];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        res[dof_z] = 0.0;
      }
      for (size_t dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          res[dof_z] += b * Bu[quad_x][dof_y][dof_z];
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        BBu[quad_x][quad_y][dof_z] = res[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_dofs_1d];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        val[dof_z] = BBu[quad_x][quad_y][dof_z];
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 res = 0.0;
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
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
           size_t num_dofs_1d,
           size_t num_quads_1d,
           size_t num_comp >
GEOSX_HOST_DEVICE
void interpolateAtQuadraturePoints( StackVariables & stack,
                                    real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                    real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp],
                                    real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Bu[num_quads_1d][num_dofs_1d][num_dofs_1d][num_comp];

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      real64 res[num_dofs_1d][num_comp];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          res[dof_z][comp] = 0.0;
        }
      }
      for (size_t dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] += b * dofs[dof_x][dof_y][dof_z][comp];
          }
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          Bu[quad_x][dof_y][dof_z][comp] = res[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BBu[num_quads_1d][num_quads_1d][num_dofs_1d][num_comp];

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
    {
      real64 res[num_dofs_1d][num_comp];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          res[dof_z][comp] = 0.0;
        }
      }
      for (size_t dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            res[dof_z][comp] += b * Bu[quad_x][dof_y][dof_z][comp];
          }
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          BBu[quad_x][quad_y][dof_z][comp] = res[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_dofs_1d][num_comp];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          val[dof_z][comp] = BBu[quad_x][quad_y][dof_z][comp];
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 res[num_comp];
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          res[comp] = 0.0;
        }
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            res[comp] += b * val[dof_z][comp];
          }
        }
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          q_values[quad_x][quad_y][quad_z][comp] = res[comp];
        }
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
           size_t num_dofs_1d,
           size_t num_quads_1d >
GEOSX_HOST_DEVICE
void applyTestFunctions( StackVariables & stack,
                         real64 const (& basis)[num_dofs_1d][num_quads_1d],
                         real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d],
                         real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;
  
  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Bu[num_dofs_1d][num_quads_1d][num_quads_1d];

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      real64 res[num_quads_1d];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        res[quad_z] = 0.0;
      }
      for (size_t quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] += b * q_values[quad_x][quad_y][quad_z]; // assumes quads in shared
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        Bu[dof_x][quad_y][quad_z] = res[quad_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BBu[num_dofs_1d][num_dofs_1d][num_quads_1d];

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
    {
      real64 res[num_quads_1d];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        res[quad_z] = 0.0;
      }
      for (size_t quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] += b * Bu[dof_x][quad_y][quad_z];
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        BBu[dof_x][dof_y][quad_z] = res[quad_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_quads_1d];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        val[quad_z] = BBu[dof_x][dof_y][quad_z];
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res = 0.0;
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
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
           size_t num_dofs_1d,
           size_t num_quads_1d,
           size_t num_comp >
GEOSX_HOST_DEVICE
void applyTestFunctions( StackVariables & stack,
                         real64 const (& basis)[num_dofs_1d][num_quads_1d],
                         real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp],
                         real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;
  
  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Bu[num_dofs_1d][num_quads_1d][num_quads_1d][num_comp];

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      real64 res[num_quads_1d][num_comp];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          res[quad_z][comp] = 0.0;
        }
      }
      for (size_t quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            res[quad_z][comp] += b * q_values[quad_x][quad_y][quad_z][comp];
          }
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          Bu[dof_x][quad_y][quad_z][comp] = res[quad_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BBu[num_dofs_1d][num_dofs_1d][num_quads_1d][num_comp];

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
    {
      real64 res[num_quads_1d][num_comp];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          res[quad_z][comp] = 0.0;
        }
      }
      for (size_t quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            res[quad_z][comp] += b * Bu[dof_x][quad_y][quad_z][comp];
          }
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          BBu[dof_x][dof_y][quad_z][comp] = res[quad_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_quads_1d][num_comp];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          val[quad_z][comp] = BBu[dof_x][dof_y][quad_z][comp];
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res[num_comp];
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          res[comp] = 0.0;
        }
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 const b = basis[dof_z][quad_z]; 
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            res[comp] += b * val[quad_z][comp];
          }
        }
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          dofs[dof_x][dof_y][dof_z][comp] = res[comp];
        }
      }
    });
  });

  ctx.teamSync();
}

// 3D Threaded version using RAJA teams
template < typename StackVariables,
           size_t num_dofs_1d,
           size_t num_quads_1d>
GEOSX_HOST_DEVICE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d],
                                            real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][3] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Bu[num_quads_1d][num_dofs_1d][num_dofs_1d],
                          Gu[num_quads_1d][num_dofs_1d][num_dofs_1d];

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      real64 bu[num_dofs_1d];
      real64 gu[num_dofs_1d];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bu[dof_z] = 0.0;
        gu[dof_z] = 0.0;
      }
      for (size_t dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const val = dofs[dof_x][dof_y][dof_z];
          bu[dof_z] += b * val; // assumes dofs in shared
          gu[dof_z] += g * val; // assumes dofs in shared
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        Bu[quad_x][dof_y][dof_z] = bu[dof_z];
        Gu[quad_x][dof_y][dof_z] = gu[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BBu[num_quads_1d][num_quads_1d][num_dofs_1d],
                          BGu[num_quads_1d][num_quads_1d][num_dofs_1d],
                          GBu[num_quads_1d][num_quads_1d][num_dofs_1d];

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
    {
      real64 bbu[num_dofs_1d];
      real64 bgu[num_dofs_1d];
      real64 gbu[num_dofs_1d];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bbu[dof_z] = 0.0;
        bgu[dof_z] = 0.0;
        gbu[dof_z] = 0.0;
      }
      for (size_t dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const bu = Bu[quad_x][dof_y][dof_z];
          real64 const gu = Gu[quad_x][dof_y][dof_z];
          bbu[dof_z] += b * bu;
          bgu[dof_z] += b * gu;
          gbu[dof_z] += g * bu;
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        BBu[quad_x][quad_y][dof_z] = bbu[dof_z];
        BGu[quad_x][quad_y][dof_z] = bgu[dof_z];
        GBu[quad_x][quad_y][dof_z] = gbu[dof_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 bbu[num_dofs_1d];
      real64 bgu[num_dofs_1d];
      real64 gbu[num_dofs_1d];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        bbu[dof_z] = BBu[quad_x][quad_y][dof_z];
        bgu[dof_z] = BGu[quad_x][quad_y][dof_z];
        gbu[dof_z] = GBu[quad_x][quad_y][dof_z];
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 bbgu = 0.0;
        real64 bgbu = 0.0;
        real64 gbbu = 0.0;
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
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
           size_t num_dofs_1d,
           size_t num_quads_1d,
           size_t num_comp >
GEOSX_HOST_DEVICE
void interpolateGradientAtQuadraturePoints( StackVariables & stack,
                                            real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                            real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                            real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp],
                                            real64 (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][3] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Bu[num_quads_1d][num_dofs_1d][num_dofs_1d][num_comp],
                          Gu[num_quads_1d][num_dofs_1d][num_dofs_1d][num_comp];

  loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      real64 bu[num_dofs_1d][num_comp];
      real64 gu[num_dofs_1d][num_comp];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          bu[dof_z][comp] = 0.0;
          gu[dof_z][comp] = 0.0;
        }
      }
      for (size_t dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            real64 const val = dofs[dof_x][dof_y][dof_z][comp];
            bu[dof_z][comp] += b * val;
            gu[dof_z][comp] += g * val;
          }
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          Bu[quad_x][dof_y][dof_z][comp] = bu[dof_z][comp];
          Gu[quad_x][dof_y][dof_z][comp] = gu[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BBu[num_quads_1d][num_quads_1d][num_dofs_1d][num_comp],
                          BGu[num_quads_1d][num_quads_1d][num_dofs_1d][num_comp],
                          GBu[num_quads_1d][num_quads_1d][num_dofs_1d][num_comp];

  loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
    {
      real64 bbu[num_dofs_1d][num_comp];
      real64 bgu[num_dofs_1d][num_comp];
      real64 gbu[num_dofs_1d][num_comp];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          bbu[dof_z][comp] = 0.0;
          bgu[dof_z][comp] = 0.0;
          gbu[dof_z][comp] = 0.0;
        }
      }
      for (size_t dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            real64 const bu = Bu[quad_x][dof_y][dof_z][comp];
            real64 const gu = Gu[quad_x][dof_y][dof_z][comp];
            bbu[dof_z][comp] += b * bu;
            bgu[dof_z][comp] += b * gu;
            gbu[dof_z][comp] += g * bu;
          }
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          BBu[quad_x][quad_y][dof_z][comp] = bbu[dof_z][comp];
          BGu[quad_x][quad_y][dof_z][comp] = bgu[dof_z][comp];
          GBu[quad_x][quad_y][dof_z][comp] = gbu[dof_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_quads_1d), [&] (size_t quad_x)
    {
      // Cache values in registers to read them only once from shared
      real64 bbu[num_dofs_1d][num_comp];
      real64 bgu[num_dofs_1d][num_comp];
      real64 gbu[num_dofs_1d][num_comp];
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          bbu[dof_z][comp] = BBu[quad_x][quad_y][dof_z][comp];
          bgu[dof_z][comp] = BGu[quad_x][quad_y][dof_z][comp];
          gbu[dof_z][comp] = GBu[quad_x][quad_y][dof_z][comp];
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        real64 bbgu[num_comp];
        real64 bgbu[num_comp];
        real64 gbbu[num_comp];
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          bbgu[comp] = 0.0;
          bgbu[comp] = 0.0;
          gbbu[comp] = 0.0;
        }
        for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            bbgu[comp] += b * bgu[dof_z][comp];
            bgbu[comp] += b * gbu[dof_z][comp];
            gbbu[comp] += g * bbu[dof_z][comp];
          }
        }
        for (size_t comp = 0; comp < num_comp; comp++)
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

// 3D Threaded version using RAJA teams
template < typename StackVariables,
           size_t num_dofs_1d,
           size_t num_quads_1d>
void applyGradientTestFunctions( StackVariables & stack,
                                 real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                 real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                 real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][3],
                                 real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Gqx[num_dofs_1d][num_quads_1d][num_quads_1d],
                          Bqy[num_dofs_1d][num_quads_1d][num_quads_1d],
                          Bqz[num_dofs_1d][num_quads_1d][num_quads_1d];

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      real64 gqx[num_quads_1d];
      real64 bqy[num_quads_1d];
      real64 bqz[num_quads_1d];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        gqx[quad_z] = 0.0;
        bqy[quad_z] = 0.0;
        bqz[quad_z] = 0.0;
      }
      for (size_t quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
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
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        Gqx[dof_x][quad_y][quad_z] = gqx[quad_z];
        Bqy[dof_x][quad_y][quad_z] = bqy[quad_z];
        Bqz[dof_x][quad_y][quad_z] = bqz[quad_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BGqx[num_dofs_1d][num_dofs_1d][num_quads_1d],
                          GBqy[num_dofs_1d][num_dofs_1d][num_quads_1d],
                          BBqz[num_dofs_1d][num_dofs_1d][num_quads_1d];

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
    {
      real64 bgqx[num_quads_1d];
      real64 gbqy[num_quads_1d];
      real64 bbqz[num_quads_1d];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        bgqx[quad_z] = 0.0;
        gbqy[quad_z] = 0.0;
        bbqz[quad_z] = 0.0;
      }
      for (size_t quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 const gqx = Gqx[dof_x][quad_y][quad_z];
          real64 const bqy = Bqy[dof_x][quad_y][quad_z];
          real64 const bqz = Bqz[dof_x][quad_y][quad_z];
          bgqx[quad_z] += b * gqx;
          gbqy[quad_z] += g * bqy;
          bbqz[quad_z] += b * bqz;
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        BGqx[dof_x][dof_y][quad_z] = bgqx[quad_z];
        GBqy[dof_x][dof_y][quad_z] = gbqy[quad_z];
        BBqz[dof_x][dof_y][quad_z] = bbqz[quad_z];
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      // Cache values in registers to read them only once from shared
      real64 bgqx[num_quads_1d];
      real64 gbqy[num_quads_1d];
      real64 bbqz[num_quads_1d];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        bgqx[quad_z] = BGqx[dof_x][dof_y][quad_z];
        gbqy[quad_z] = GBqy[dof_x][dof_y][quad_z];
        bbqz[quad_z] = BBqz[dof_x][dof_y][quad_z];
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res = 0.0;
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
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
           size_t num_dofs_1d,
           size_t num_quads_1d,
           size_t num_comp >
void applyGradientTestFunctions( StackVariables & stack,
                                 real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                 real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                 real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][num_comp][3],
                                 real64 (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  RAJA_TEAM_SHARED real64 Gqx[num_dofs_1d][num_quads_1d][num_quads_1d][num_comp],
                          Bqy[num_dofs_1d][num_quads_1d][num_quads_1d][num_comp],
                          Bqz[num_dofs_1d][num_quads_1d][num_quads_1d][num_comp];

  loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (size_t quad_y)
  {
    loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      real64 gqx[num_quads_1d][num_comp];
      real64 bqy[num_quads_1d][num_comp];
      real64 bqz[num_quads_1d][num_comp];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          gqx[quad_z][comp] = 0.0;
          bqy[quad_z][comp] = 0.0;
          bqz[quad_z][comp] = 0.0;
        }
      }
      for (size_t quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        real64 const g = basis_gradient[dof_x][quad_x];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
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
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          Gqx[dof_x][quad_y][quad_z][comp] = gqx[quad_z][comp];
          Bqy[dof_x][quad_y][quad_z][comp] = bqy[quad_z][comp];
          Bqz[dof_x][quad_y][quad_z][comp] = bqz[quad_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the second dimension
  RAJA_TEAM_SHARED real64 BGqx[num_dofs_1d][num_dofs_1d][num_quads_1d][num_comp],
                          GBqy[num_dofs_1d][num_dofs_1d][num_quads_1d][num_comp],
                          BBqz[num_dofs_1d][num_dofs_1d][num_quads_1d][num_comp];

  loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
  {
    loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
    {
      real64 bgqx[num_quads_1d][num_comp];
      real64 gbqy[num_quads_1d][num_comp];
      real64 bbqz[num_quads_1d][num_comp];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          bgqx[quad_z][comp] = 0.0;
          gbqy[quad_z][comp] = 0.0;
          bbqz[quad_z][comp] = 0.0;
        }
      }
      for (size_t quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        real64 const g = basis_gradient[dof_y][quad_y];
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            real64 const gqx = Gqx[dof_x][quad_y][quad_z][comp];
            real64 const bqy = Bqy[dof_x][quad_y][quad_z][comp];
            real64 const bqz = Bqz[dof_x][quad_y][quad_z][comp];
            bgqx[quad_z][comp] += b * gqx;
            gbqy[quad_z][comp] += g * bqy;
            bbqz[quad_z][comp] += b * bqz;
          }
        }
      }
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          BGqx[dof_x][dof_y][quad_z][comp] = bgqx[quad_z][comp];
          GBqy[dof_x][dof_y][quad_z][comp] = gbqy[quad_z][comp];
          BBqz[dof_x][dof_y][quad_z][comp] = bbqz[quad_z][comp];
        }
      }
    });
  });

  ctx.teamSync();

  // Contraction on the third dimension
  RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_y)
  {
    RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, num_dofs_1d), [&] (size_t dof_x)
    {
      // Cache values in registers to read them only once from shared
      real64 bgqx[num_quads_1d][num_comp];
      real64 gbqy[num_quads_1d][num_comp];
      real64 bbqz[num_quads_1d][num_comp];
      for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          bgqx[quad_z][comp] = BGqx[dof_x][dof_y][quad_z][comp];
          gbqy[quad_z][comp] = GBqy[dof_x][dof_y][quad_z][comp];
          bbqz[quad_z][comp] = BBqz[dof_x][dof_y][quad_z][comp];
        }
      }
      for (size_t dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res[num_comp];
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          res[comp] = 0.0;
        }
        for (size_t quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          for (size_t comp = 0; comp < num_comp; comp++)
          {
            res[comp] += b * bgqx[quad_z][comp]
                       + b * gbqy[quad_z][comp]
                       + g * bbqz[quad_z][comp];
          }
        }
        for (size_t comp = 0; comp < num_comp; comp++)
        {
          dofs[dof_x][dof_y][dof_z][comp] = res[comp];
        }
      }
    });
  });

  ctx.teamSync();
}

template < typename StackVariables,
           typename Field,
           size_t stride_x, size_t stride_y, size_t stride_z >
void readField( StackVariables & stack,
                Field & field,
                real64 (& local_field)[stride_x][stride_y][stride_z] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, stride_x), [&] (size_t ind_x)
  {
    RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, stride_y), [&] (size_t ind_y)
    {
      for (size_t ind_z = 0; ind_z < stride_z; ind_z++)
      {
        size_t const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        size_t const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        local_field[ ind_x ][ ind_y ][ ind_z ] = field[ global_node_index ];
      }
    });
  });
}

template < typename StackVariables,
           typename Field,
           size_t stride_x, size_t stride_y, size_t stride_z, size_t dim >
void readField( StackVariables & stack,
                Field & field,
                real64 (& local_field)[stride_x][stride_y][stride_z][dim] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, stride_x), [&] (size_t ind_x)
  {
    RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, stride_y), [&] (size_t ind_y)
    {
      for (size_t ind_z = 0; ind_z < stride_z; ind_z++)
      {
        size_t const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        size_t const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        for (size_t d = 0; d < dim; d++)
        {
          local_field[ ind_x ][ ind_y ][ ind_z ][ d ] = field[ global_node_index ][ d ];
        }
      }
    });
  });
}

template < typename StackVariables,
           typename Field,
           size_t stride_x, size_t stride_y, size_t stride_z >
void writeField( StackVariables & stack,
                 real64 const (& local_field)[stride_x][stride_y][stride_z],
                 Field & field )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, stride_x), [&] (size_t ind_x)
  {
    RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, stride_y), [&] (size_t ind_y)
    {
      for (size_t ind_z = 0; ind_z < stride_z; ind_z++)
      {
        size_t const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        size_t const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        field[ global_node_index ] = local_field[ ind_x ][ ind_y ][ ind_z ];
      }
    });
  });
}

template < typename StackVariables,
           typename Field,
           size_t stride_x, size_t stride_y, size_t stride_z, size_t dim >
void writeField( StackVariables & stack,
                 real64 const (& local_field)[stride_x][stride_y][stride_z][dim], 
                 Field & field )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, stride_x), [&] (size_t ind_x)
  {
    RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, stride_y), [&] (size_t ind_y)
    {
      for (size_t ind_z = 0; ind_z < stride_z; ind_z++)
      {
        size_t const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        size_t const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        for (size_t d = 0; d < dim; d++)
        {
          field[ global_node_index ][ d ] = local_field[ ind_x ][ ind_y ][ ind_z ][ d ];
        }
      }
    });
  });
}

template < typename StackVariables,
           typename Field,
           size_t stride_x, size_t stride_y, size_t stride_z >
void writeAddField( StackVariables & stack,
                    real64 const (& local_field)[stride_x][stride_y][stride_z],
                    Field & field )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, stride_x), [&] (size_t ind_x)
  {
    RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, stride_y), [&] (size_t ind_y)
    {
      for (size_t ind_z = 0; ind_z < stride_z; ind_z++)
      {
        size_t const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        size_t const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        atomicAdd( field[ global_node_index ], local_field[ ind_x ][ ind_y ][ ind_z ]);
      }
    });
  });
}

template < typename StackVariables,
           typename Field,
           size_t stride_x, size_t stride_y, size_t stride_z, size_t dim >
void writeAddField( StackVariables & stack,
                    real64 const (& local_field)[stride_x][stride_y][stride_z][dim], 
                    Field & field )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  RAJA::expt::loop<thread_x> (ctx, RAJA::RangeSegment(0, stride_x), [&] (size_t ind_x)
  {
    RAJA::expt::loop<thread_y> (ctx, RAJA::RangeSegment(0, stride_y), [&] (size_t ind_y)
    {
      for (size_t ind_z = 0; ind_z < stride_z; ind_z++)
      {
        size_t const local_node_index = ind_x + stride_x * ( ind_y + stride_y * ind_z );
        size_t const global_node_index = stack.kernelComponent.m_elemsToNodes( stack.element_index, local_node_index );
        for (size_t d = 0; d < dim; d++)
        {
          atomicAdd( field[ global_node_index ][ d ], local_field[ ind_x ][ ind_y ][ ind_z ][ d ] );
        }
      }
    });
  });
}

} // namespace finiteElement
} // namespace geosx



#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_HPP_ */
