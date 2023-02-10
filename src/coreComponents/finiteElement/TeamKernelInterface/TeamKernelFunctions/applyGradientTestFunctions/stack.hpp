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
 * @file applyGradientTestFunctions/stack.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_STACK_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_STACK_HPP_

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

///////////////////
// Stack

template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyGradientTestFunctions(
  StackVariables & stack,
  real64 const (& basis)[num_dofs_1d][num_quads_1d],
  real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
  tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_quads_1d, 3 > const & q_values,
  tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > & dofs )
{
  // Contraction on the first dimension
  tensor::StaticDTensor< num_dofs_1d, num_quads_1d, num_quads_1d > Gqx, Bqy, Bqz;
  #pragma unroll
  for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
  {
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 gqx = 0.0;
        real64 bqy = 0.0;
        real64 bqz = 0.0;
        #pragma unroll
        for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
        {
          real64 const b = basis[ dof_x ][ quad_x ];
          real64 const g = basis_gradient[ dof_x ][ quad_x ];

          real64 const qx = q_values( quad_x, quad_y, quad_z, 0 );
          real64 const qy = q_values( quad_x, quad_y, quad_z, 1 );
          real64 const qz = q_values( quad_x, quad_y, quad_z, 2 );

          gqx = gqx + g * qx;
          bqy = bqy + b * qy;
          bqz = bqz + b * qz;
        }
        Gqx( dof_x, quad_y, quad_z ) = gqx;
        Bqy( dof_x, quad_y, quad_z ) = bqy;
        Bqz( dof_x, quad_y, quad_z ) = bqz;
      }
    }
  }

  // Contraction on the second dimension
  tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_quads_1d > BGqx, GBqy, BBqz;
  #pragma unroll
  for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
  {
    #pragma unroll
    for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    {
      #pragma unroll
      for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
      {
        real64 bgqx = 0.0;
        real64 gbqy = 0.0;
        real64 bbqz = 0.0;
        #pragma unroll
        for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
        {
          real64 const b = basis[dof_y][quad_y];
          real64 const g = basis_gradient[dof_y][quad_y];
          real64 const gqx = Gqx( dof_x, quad_y, quad_z );
          real64 const bqy = Bqy( dof_x, quad_y, quad_z );
          real64 const bqz = Bqz( dof_x, quad_y, quad_z );
          bgqx += b * gqx;
          gbqy += g * bqy;
          bbqz += b * bqz;
        }
        BGqx( dof_x, dof_y, quad_z ) = bgqx;
        GBqy( dof_x, dof_y, quad_z ) = gbqy;
        BBqz( dof_x, dof_y, quad_z ) = bbqz;
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
      #pragma unroll
      for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
      {
        real64 res = 0.0;
        #pragma unroll
        for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
        {
          real64 const bgqx = BGqx( dof_x, dof_y, quad_z );
          real64 const gbqy = GBqy( dof_x, dof_y, quad_z );
          real64 const bbqz = BBqz( dof_x, dof_y, quad_z );
          real64 const b = basis[dof_z][quad_z];
          real64 const g = basis_gradient[dof_z][quad_z];
          res = res + b * bgqx;
          res = res + b * gbqy;
          res = res + g * bbqz;
        }
        dofs( dof_x, dof_y, dof_z ) = res;
      }
    }
  }
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
  tensor::StaticDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp, 3 > const & q_values,
  tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & dofs )
{
  #pragma unroll
  for (localIndex comp = 0; comp < num_comp; comp++)
  {
    // Contraction on the first dimension
    tensor::StaticDTensor< num_dofs_1d, num_quads_1d, num_quads_1d > Gqx, Bqy, Bqz;
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      #pragma unroll
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        #pragma unroll
        for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
        {
          real64 gqx = 0.0;
          real64 bqy = 0.0;
          real64 bqz = 0.0;
          #pragma unroll
          for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
          {
            real64 const b = 1.0;//basis[dof_x][quad_x];
            real64 const g = 1.0;//basis_gradient[dof_x][quad_x];
            real64 const qx = q_values( quad_x, quad_y, quad_z, comp, 0 );
            real64 const qy = q_values( quad_x, quad_y, quad_z, comp, 1 );
            real64 const qz = q_values( quad_x, quad_y, quad_z, comp, 2 );
            gqx = gqx + g * qx;
            bqy = bqy + b * qy;
            bqz = bqz + b * qz;
          }
          Gqx( dof_x, quad_y, quad_z ) = gqx;
          Bqy( dof_x, quad_y, quad_z ) = bqy;
          Bqz( dof_x, quad_y, quad_z ) = bqz;
        }
      }
    }

    // Contraction on the second dimension
    tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_quads_1d > BGqx, GBqy, BBqz;
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        #pragma unroll
        for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
        {
          real64 bgqx = 0.0;
          real64 gbqy = 0.0;
          real64 bbqz = 0.0;
          #pragma unroll
          for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
          {
            real64 const b = 1.0;//basis[dof_y][quad_y];
            real64 const g = 1.0;//basis_gradient[dof_y][quad_y];
            real64 const gqx = Gqx( dof_x, quad_y, quad_z );
            real64 const bqy = Bqy( dof_x, quad_y, quad_z );
            real64 const bqz = Bqz( dof_x, quad_y, quad_z );
            bgqx = bgqx + b * gqx;
            gbqy = gbqy + g * bqy;
            bbqz = bbqz + b * bqz;
          }
          BGqx( dof_x, dof_y, quad_z ) = bgqx;
          GBqy( dof_x, dof_y, quad_z ) = gbqy;
          BBqz( dof_x, dof_y, quad_z ) = bbqz;
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
        #pragma unroll
        for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
        {
          real64 res = 0.0;
          #pragma unroll
          for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
          {
            real64 const b = 1.0;//basis[dof_z][quad_z];
            real64 const g = 1.0;//basis_gradient[dof_z][quad_z];
            real64 bgqx = BGqx( dof_x, dof_y, quad_z );
            real64 gbqy = GBqy( dof_x, dof_y, quad_z );
            real64 bbqz = BBqz( dof_x, dof_y, quad_z );
            res = res + b * bgqx;
            res = res + b * gbqy;
            res = res + g * bbqz;
          }
          dofs( dof_x, dof_y, dof_z, comp ) = res;
        }
      }
    }
    // #pragma unroll
    // for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    // {
    //   #pragma unroll
    //   for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
    //   {
    //     #pragma unroll
    //     for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
    //     {
    //       // real64 res = 0.0;
    //       // #pragma unroll
    //       // for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    //       // {
    //       //   real64 const b = basis[dof_z][quad_z];
    //       //   real64 const g = basis_gradient[dof_z][quad_z];
    //       //   real64 bgqx = BGqx( dof_x, dof_y, quad_z );
    //       //   real64 gbqy = GBqy( dof_x, dof_y, quad_z );
    //       //   real64 bbqz = BBqz( dof_x, dof_y, quad_z );
    //       //   res = res + b * bgqx;
    //       //   res = res + b * gbqy;
    //       //   res = res + g * bbqz;
    //       // }
    //       dofs( dof_x, dof_y, dof_z, comp ) = q_values( dof_x, dof_y, dof_z, comp, 0 )
    //                                         + q_values( dof_x, dof_y, dof_z, comp, 1 )
    //                                         + q_values( dof_x, dof_y, dof_z, comp, 2 );
    //     }
    //   }
    // }
  }
  // loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
  //         [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z){
  //   real64 v[ num_comp ];
  //   #pragma unroll
  //  for (localIndex comp = 0; comp < num_comp; comp++)
  //   {
  //     v[ comp ] = 0.0;
  //   }
  //   #pragma unroll
  //  for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
  //   {
  //     real64 const bz = 1.0;//Bz[quad_z];
  //     real64 const gz = 1.0;//Gz[quad_z];
  //     #pragma unroll
  //    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
  //     {
  //       real64 const by = 1.0;//By[quad_y];
  //       real64 const gy = 1.0;//Gy[quad_y];
  //       #pragma unroll
  //      for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
  //       {
  //         real64 const bx = 1.0; //Bx[quad_x];
  //         real64 const gx = 1.0; //Gx[quad_x];
  //         real64 const dx = gx * by * bz;
  //         real64 const dy = bx * gy * bz;
  //         real64 const dz = bx * by * gz;
  //         // localIndex const srcLane = quad_x + num_quads_1d * ( quad_y + num_quads_1d * ( quad_z + num_quads_1d * stack.batch_index ) );
  //         #pragma unroll
  //        for (localIndex comp = 0; comp < num_comp; comp++)
  //         {
  //           real64 const valx = q_values( quad_x, quad_y, quad_z, comp, 0 );
  //           // real64 const valx = __shfl_sync(0xffffffff, q_values( quad_x, quad_y, quad_z, comp, 0 ), srcLane);
  //           real64 const valy = q_values( quad_x, quad_y, quad_z, comp, 1 );
  //           // real64 const valy = __shfl_sync(0xffffffff, q_values( quad_x, quad_y, quad_z, comp, 1 ), srcLane);
  //           real64 const valz = q_values( quad_x, quad_y, quad_z, comp, 2 );
  //           // real64 const valz = __shfl_sync(0xffffffff, q_values( quad_x, quad_y, quad_z, comp, 2 ), srcLane);
  //           v[ comp ] += dx * valx;
  //           v[ comp ] += dy * valy;
  //           v[ comp ] += dz * valz;
  //         }
  //       }
  //     }
  //   }
  //   #pragma unroll
  //  for (localIndex comp = 0; comp < num_comp; comp++)
  //   {
  //     dofs( dof_x, dof_y, dof_z, comp ) = v[ comp ];
  //   }
  // } );
}

} // impl

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_GRAD_TEST_STACK_HPP_ */
