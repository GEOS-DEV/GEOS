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
 * @file interpolateGradientAtQuadraturePoint.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_GRAD_SINGLE_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_GRAD_SINGLE_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_types.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/common.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

// Stack tensor
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoint( StackVariables & stack,
                                           geosx::TensorIndex const & quad_index,
                                           real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                           real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                           real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d],
                                           real64 (& q_value)[3] )
{
  q_value[ 0 ] = 0.0;
  q_value[ 1 ] = 0.0;
  q_value[ 2 ] = 0.0;
  #pragma unroll
  for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
  {
    real64 const bz = basis[dof_z][quad_index.z];
    real64 const gz = basis_gradient[dof_z][quad_index.z];
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      real64 const by = basis[dof_y][quad_index.y];
      real64 const gy = basis_gradient[dof_y][quad_index.y];
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const bx = basis[dof_x][quad_index.x];
        real64 const gx = basis_gradient[dof_x][quad_index.x];
        real64 const val = dofs[ dof_x ][ dof_y ][ dof_z ];
        q_value[ 0 ] += gx * by * bz * val;
        q_value[ 1 ] += bx * gy * bz * val;
        q_value[ 2 ] += bx * by * gz * val;
      }
    }
  }
}

// Stack version using RAJA teams
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoint( StackVariables & stack,
                                           geosx::TensorIndex const & quad_index,
                                           real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                           real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                           real64 const (& dofs)[num_dofs_1d][num_dofs_1d][num_dofs_1d][num_comp],
                                           real64 (& q_value)[num_comp][3] )
{
  #pragma unroll
  for (localIndex comp = 0; comp < num_comp; comp++)
  {
    q_value[ comp ][ 0 ] = 0.0;
    q_value[ comp ][ 1 ] = 0.0;
    q_value[ comp ][ 2 ] = 0.0;
  }
  #pragma unroll
  for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
  {
    real64 const bz = basis[dof_z][quad_index.z];
    real64 const gz = basis_gradient[dof_z][quad_index.z];
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      real64 const by = basis[dof_y][quad_index.y];
      real64 const gy = basis_gradient[dof_y][quad_index.y];
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const bx = basis[dof_x][quad_index.x];
        real64 const gx = basis_gradient[dof_x][quad_index.x];
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          real64 const val = dofs[ dof_x ][ dof_y ][ dof_z ][ comp ];
          q_value[ comp ][ 0 ] += gx * by * bz * val;
          q_value[ comp ][ 1 ] += bx * gy * bz * val;
          q_value[ comp ][ 2 ] += bx * by * gz * val;
        }
      }
    }
  }
}

// Stack tensor
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoint( StackVariables & stack,
                                           geosx::TensorIndex const & quad_index,
                                           real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                           real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                           tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > const & dofs,
                                           real64 (& q_value)[3] )
{
  q_value[ 0 ] = 0.0;
  q_value[ 1 ] = 0.0;
  q_value[ 2 ] = 0.0;
  #pragma unroll
  for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
  {
    real64 const bz = basis[dof_z][quad_index.z];
    real64 const gz = basis_gradient[dof_z][quad_index.z];
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      real64 const by = basis[dof_y][quad_index.y];
      real64 const gy = basis_gradient[dof_y][quad_index.y];
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const bx = basis[dof_x][quad_index.x];
        real64 const gx = basis_gradient[dof_x][quad_index.x];
        real64 const val = dofs( dof_x, dof_y, dof_z );
        q_value[ 0 ] += gx * by * bz * val;
        q_value[ 1 ] += bx * gy * bz * val;
        q_value[ 2 ] += bx * by * gz * val;
      }
    }
  }
}

// Stack version using RAJA teams
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoint( StackVariables & stack,
                                          //  geosx::TensorIndex const & quad_index,
                                           localIndex const & qx,
                                           localIndex const & qy,
                                           localIndex const & qz,
                                           real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                           real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                           tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & dofs,
                                           real64 (& q_value)[num_comp][3] )
{
  // real64 const w = 0.57735026919;
  // real64 basis_qx[num_dofs_1d];
  // real64 basis_gradient_qx[num_dofs_1d];
  // real64 basis_qy[num_dofs_1d];
  // real64 basis_gradient_qy[num_dofs_1d];
  // real64 basis_qz[num_dofs_1d];
  // real64 basis_gradient_qz[num_dofs_1d];
  // #pragma unroll
  // for (localIndex i = 0; i < num_dofs_1d; i++)
  // {
  //   // basis_q[i] = finiteElement::LagrangeBasis1::value(
  //   //       i, w * finiteElement::LagrangeBasis1::parentSupportCoord( qx ) );
  //   // basis_gradient_q[i] = finiteElement::LagrangeBasis1::gradient(
  //   //       i, w * finiteElement::LagrangeBasis1::parentSupportCoord( qx ) );
  //   basis_qx[i] = 0.1 * (i );//+ qx);
  //   basis_gradient_qx[i] = 0.2 * (i );//+ qx);
  //   basis_qy[i] = 0.1 * (i );//+ qy);
  //   basis_gradient_qy[i] = 0.2 * (i );//+ qy);
  //   basis_qz[i] = 0.1 * (i );//+ qz);
  //   basis_gradient_qz[i] = 0.2 * (i );//+ qz);
  // }
  // constexpr real64 basis[num_dofs_1d][num_quads_1d] = { { 0.123123123, 123.123 }, {0.64563452, 45.23472} };
  // constexpr real64 basis_gradient[num_dofs_1d][num_quads_1d] = { { 0.643542346578, 0.4567812 }, {0.98765234, 0.3457823} };
  // constexpr real64 val0 = 0.3455756454;
  // constexpr real64 val1 = 0.8765643458;
  // constexpr real64 grad0 = 0.12343765432;
  // constexpr real64 grad1 = 0.9876524569;

  #pragma unroll
  for (localIndex comp = 0; comp < num_comp; comp++)
  {
    q_value[ comp ][ 0 ] = 0.0;
    q_value[ comp ][ 1 ] = 0.0;
    q_value[ comp ][ 2 ] = 0.0;
  }
  #pragma unroll
  for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
  {
    // real64 const bz = basis_qz[dof_z];
    // real64 const gz = basis_gradient_qz[dof_z];
    real64 const bz = basis[dof_z][qz];
    real64 const gz = basis_gradient[dof_z][qz];
    // real64 const bz = dof_z == qz ? val0 : val1;
    // real64 const gz = dof_z == qz ? grad0 : grad1;
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      // real64 const by = basis_qy[dof_y];
      // real64 const gy = basis_gradient_qy[dof_y];
      real64 const by = basis[dof_y][qy];
      real64 const gy = basis_gradient[dof_y][qy];
      // real64 const by = dof_y == qy ? val0 : val1;
      // real64 const gy = dof_y == qy ? grad0 : grad1;
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        // real64 const bx = basis_qx[dof_x];
        // real64 const gx = basis_gradient_qx[dof_x];
        real64 const bx = basis[dof_x][qx];
        real64 const gx = basis_gradient[dof_x][qx];
        // real64 const bx = dof_x == qx ? val0 : val1;
        // real64 const gx = dof_x == qx ? grad0 : grad1;
        real64 const Gx = gx * by * bz;
        real64 const Gy = bx * gy * bz;
        real64 const Gz = bx * by * gz;
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          real64 const val = dofs( dof_x, dof_y, dof_z, comp );
          q_value[ comp ][ 0 ] += Gx * val;
          q_value[ comp ][ 1 ] += Gy * val;
          q_value[ comp ][ 2 ] += Gz * val;
        }
      }
    }
  }
}

// Stack tensor
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoint( StackVariables & stack,
                                           geosx::TensorIndex const & quad_index,
                                           real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                           real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                           tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > const & dofs,
                                           tensor::StaticDTensor< 3 > & q_value )
{
  q_value( 0 ) = 0.0;
  q_value( 1 ) = 0.0;
  q_value( 2 ) = 0.0;
  #pragma unroll
  for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
  {
    real64 const bz = basis[dof_z][quad_index.z];
    real64 const gz = basis_gradient[dof_z][quad_index.z];
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      real64 const by = basis[dof_y][quad_index.y];
      real64 const gy = basis_gradient[dof_y][quad_index.y];
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const bx = basis[dof_x][quad_index.x];
        real64 const gx = basis_gradient[dof_x][quad_index.x];
        real64 const val = dofs( dof_x, dof_y, dof_z );
        q_value( 0 ) += gx * by * bz * val;
        q_value( 1 ) += bx * gy * bz * val;
        q_value( 2 ) += bx * by * gz * val;
      }
    }
  }
}

// Stack version using RAJA teams
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoint( StackVariables & stack,
                                           geosx::TensorIndex const & quad_index,
                                           real64 const (& basis)[num_dofs_1d][num_quads_1d],
                                           real64 const (& basis_gradient)[num_dofs_1d][num_quads_1d],
                                           tensor::StaticDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > const & dofs,
                                           tensor::StaticDTensor< num_comp, 3 > & q_value )
{
  #pragma unroll
  for (localIndex comp = 0; comp < num_comp; comp++)
  {
    q_value( comp, 0 ) = 0.0;
    q_value( comp, 1 ) = 0.0;
    q_value( comp, 2 ) = 0.0;
  }
  #pragma unroll
  for (localIndex dof_z = 0; dof_z < num_dofs_1d; dof_z++)
  {
    real64 const bz = basis[dof_z][quad_index.z];
    real64 const gz = basis_gradient[dof_z][quad_index.z];
    #pragma unroll
    for (localIndex dof_y = 0; dof_y < num_dofs_1d; dof_y++)
    {
      real64 const by = basis[dof_y][quad_index.y];
      real64 const gy = basis_gradient[dof_y][quad_index.y];
      #pragma unroll
      for (localIndex dof_x = 0; dof_x < num_dofs_1d; dof_x++)
      {
        real64 const bx = basis[dof_x][quad_index.x];
        real64 const gx = basis_gradient[dof_x][quad_index.x];
        #pragma unroll
        for (localIndex comp = 0; comp < num_comp; comp++)
        {
          real64 const val = dofs( dof_x, dof_y, dof_z, comp );
          q_value( comp, 0 ) += gx * by * bz * val;
          q_value( comp, 1 ) += bx * gy * bz * val;
          q_value( comp, 2 ) += bx * by * gz * val;
        }
      }
    }
  }
}

template < typename StackVariables,
           typename Basis,
           typename Dofs,
           typename QValues >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void interpolateGradientAtQuadraturePoint( StackVariables & stack,
                                          //  geosx::TensorIndex const & quad_index,
  localIndex const & qx,
  localIndex const & qy,
  localIndex const & qz,
                                           Basis const & basis,
                                           Dofs const & dofs,
                                           QValues & q_values )
{
  interpolateGradientAtQuadraturePoint( stack,
                                        // quad_index,
                                        qx, qy, qz,
                                        basis.getValuesAtQuadPts(),
                                        basis.getGradientValuesAtQuadPts(),
                                        dofs,
                                        q_values );
}

} // namespace finiteElement
} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_INTERP_GRAD_SINGLE_HPP_ */
