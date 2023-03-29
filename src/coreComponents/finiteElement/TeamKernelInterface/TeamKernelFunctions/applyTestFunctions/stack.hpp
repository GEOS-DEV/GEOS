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
 * @file applyTestFunctions/stack.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_STACK_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_STACK_HPP_

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

// Stack tensor

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
                         tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_quads_1d > const & q_values,
                         tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > & dofs )
{
  using RAJA::RangeSegment;
  RAJA::LaunchContext & ctx = stack.ctx;

  // Contraction on the first dimension
  tensor::StaticDTensor< num_dofs_1d, num_quads_1d, num_quads_1d > Bu;
  #pragma unroll
  for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
  {
    #pragma unroll
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      real64 res[num_quads_1d];
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        res[quad_z] = 0.0;
      }
      #pragma unroll
      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
      {
        real64 const b = basis[dof_x][quad_x];
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] += b * q_values( quad_x, quad_y, quad_z ); // assumes quads in shared
        }
      }
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        Bu( dof_x, quad_y, quad_z ) = res[quad_z];
      }
    }
  }

  // Contraction on the second dimension
  tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_quads_1d > BBu;
  #pragma unroll
  for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
  {
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      real64 res[num_quads_1d];
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        res[quad_z] = 0.0;
      }
      #pragma unroll
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const b = basis[dof_y][quad_y];
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] += b * Bu( dof_x, quad_y, quad_z );
        }
      }
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        BBu( dof_x, dof_y, quad_z ) = res[quad_z];
      }
    }
  }

  // Contraction on the third dimension
  #pragma unroll
  for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
  {
    #pragma unroll
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      // Cache values in registers to read them only once from shared
      real64 val[num_quads_1d];
      #pragma unroll
      for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
      {
        val[quad_z] = BBu( dof_x, dof_y, quad_z );
      }
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res = 0.0;
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res += basis[dof_z][quad_z] * val[quad_z];
        }
        dofs( dof_x, dof_y, dof_z ) = res;
      }
    }
  }
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
                         tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp > const & q_values,
                         tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & dofs )
{
  using RAJA::RangeSegment;
  RAJA::LaunchContext & ctx = stack.ctx;

  #pragma unroll
  for (localIndex comp = 0; comp < num_comp; comp++)
  {
    // Contraction on the first dimension
    tensor::StaticDTensor< num_dofs_1d, num_quads_1d, num_quads_1d > Bu;
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 res[num_quads_1d];
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] = 0.0;
        }
        #pragma unroll
        for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
        {
          real64 const b = basis[dof_x][quad_x];
          #pragma unroll
          for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
          {
            res[quad_z] += b * q_values( quad_x, quad_y, quad_z, comp );
          }
        }
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          Bu( dof_x, quad_y, quad_z ) = res[quad_z];
        }
      }
    }

    // Contraction on the second dimension
    tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_quads_1d > BBu;
    #pragma unroll
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 res[num_quads_1d];
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          res[quad_z] = 0.0;
        }
        #pragma unroll
        for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
        {
          real64 const b = basis[dof_y][quad_y];
          #pragma unroll
          for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
          {
            res[quad_z] += b * Bu( dof_x, quad_y, quad_z );
          }
        }
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          BBu( dof_x, dof_y, quad_z ) = res[quad_z];
        }
      }
    }

    // Contraction on the third dimension
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        // Cache values in registers to read them only once from shared
        real64 val[num_quads_1d];
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
            val[quad_z] = BBu( dof_x, dof_y, quad_z );
        }
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 res = 0.0;
          #pragma unroll
          for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
          {
            real64 const b = basis[dof_z][quad_z];
            res += b * val[quad_z];
          }
          dofs( dof_x, dof_y, dof_z, comp ) = res;
        }
      }
    }
  }
}

} // namespace impl

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_STACK_HPP_ */
