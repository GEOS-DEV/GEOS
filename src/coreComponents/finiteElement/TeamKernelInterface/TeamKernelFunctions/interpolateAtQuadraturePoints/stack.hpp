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
 * @file interpolateAtQuadraturePoints/stack.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_STACK_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_STACK_HPP_

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
  RAJA::LaunchContext & ctx = stack.ctx;

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
  RAJA::LaunchContext & ctx = stack.ctx;

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

} // namespace impl

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_STACK_HPP_ */
