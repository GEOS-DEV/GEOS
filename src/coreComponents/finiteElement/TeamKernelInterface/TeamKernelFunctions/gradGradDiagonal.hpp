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
 * @file gradGradDiagonal.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_GRADGRADDIAG_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_GRADGRADDIAG_HPP_

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

template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void computeGradGradLocalDiagonal( StackVariables & stack,
                                   real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                   real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                   real64 const (& q_values)[num_quads_1d][num_quads_1d][num_quads_1d][3][3],
                                   real64 (& diag)[num_dofs_1d][num_dofs_1d][num_dofs_1d] )
{
  using RAJA::RangeSegment;
  LaunchContext & ctx = stack.ctx;

  constexpr localIndex DIM = 3;

  for (int i = 0; i < DIM; ++i)
  {
    for (int j = 0; j < DIM; ++j)
    {
      // first tensor contraction, along z direction
      SharedTensor< num_quads_1d, num_quads_1d, num_dofs_1d > QQD( stack.shared_mem[0] );

      loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
      {
        loop<thread_y> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_y)
        {
          for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
          {
            QQD( quad_x, quad_y, dof_z ) = 0.0;
            for (int quad_z = 0; quad_z < num_quads_1d; ++quad_z)
            {
              const double Bz = basis[dof_z][quad_z];
              const double Gz = basis_gradient[dof_z][quad_z];
              const double L = i==2 ? Gz : Bz;
              const double R = j==2 ? Gz : Bz;
              const double O = q_values[quad_x][quad_y][quad_z][i][j];
              QQD( quad_x, quad_y, dof_z ) += L * O * R;
            }
          }
        } );
      } );
      
      ctx.teamSync();

      // second tensor contraction, along y direction
      SharedTensor< num_quads_1d, num_dofs_1d, num_dofs_1d > QDD( stack.shared_mem[1] );

      loop<thread_x> (ctx, RangeSegment(0, num_quads_1d), [&] (localIndex quad_x)
      {
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
          {
            QDD( quad_x, dof_y, dof_z ) = 0.0;
            for (int quad_y = 0; quad_y < num_quads_1d; ++quad_y)
            {
              const double By = basis[dof_y][quad_y];
              const double Gy = basis_gradient[dof_y][quad_y];
              const double L = i==1 ? Gy : By;
              const double R = j==1 ? Gy : By;
              QDD( quad_x, dof_y, dof_z ) += L * QQD( quad_x, quad_y, dof_z ) * R;
            }
          } );
        }
      } );
      
      ctx.teamSync();

      // third tensor contraction, along x direction
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        loop<thread_y> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_y)
        {
          loop<thread_x> (ctx, RangeSegment(0, num_dofs_1d), [&] (localIndex dof_x)
          {
            for (int quad_x = 0; quad_x < num_quads_1d; ++quad_x)
            {
              const double Bx = basis[dof_x][quad_x];
              const double Gx = basis_gradient[dof_x][quad_x];
              const double L = i==0 ? Gx : Bx;
              const double R = j==0 ? Gx : Bx;
              diag[dof_x][dof_y][dof_z] += L * QDD( quad_x, dof_y, dof_z ) * R;
            }
          } );
        } );
      }
      
      ctx.teamSync();
    }
  }
}

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_GRADGRADDIAG_HPP_ */
