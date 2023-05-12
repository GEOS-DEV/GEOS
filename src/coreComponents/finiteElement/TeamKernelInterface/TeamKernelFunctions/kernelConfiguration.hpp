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
 * @file kernelConfiguration.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_KERNEL_CONFIG_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_KERNEL_CONFIG_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_types.hpp"
#include "tensor/tensor_traits.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/Basis.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/SharedMemBuffers.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

/**
 * @brief A type describing the threading strategy and other aspects about the kernel configuration in a forallElements.
 * 
 * @tparam threading_model The threading model to compute on each element.
 * @tparam num_threads_1d The number of threads per threading dimension, this is usually equal to the number of quadrature points in 1D.
 * @tparam batch_size The number of elements treated in a same block of threads, used to optimize the number of threads per block of threads.
 */
template < ThreadingModel threading_model, localIndex num_threads_1d, localIndex batch_size >
struct KernelConfiguration;

// @brief Serial kernel configuration for kernels using `forallElements`.
template < localIndex num_threads_1d, localIndex batchSize >
struct KernelConfiguration< ThreadingModel::Serial, num_threads_1d, batchSize >
{
  static constexpr ThreadingModel threading_model = ThreadingModel::Serial;
  static constexpr localIndex num_threads_x = 1;
  static constexpr localIndex num_threads_y = 1;
  static constexpr localIndex num_threads_z = batchSize;
  static constexpr localIndex batch_size = batchSize;

  template < typename T, localIndex... Dims >
  using Tensor = tensor::StaticTensor< T, Dims... >;

  template < localIndex num_dofs, localIndex num_quads >
  using Basis = geosx::stackVariables::template StackBasis< num_dofs, num_quads >;

  // Index inside a batch of elements
  localIndex batch_index;
  // RAJA launch context
  RAJA::LaunchContext & ctx;

  GEOSX_HOST_DEVICE
  KernelConfiguration( RAJA::LaunchContext & ctx )
  : ctx( ctx )
  { }
};

// Distributed 1D
template < localIndex num_threads_1d, localIndex batchSize >
struct KernelConfiguration< ThreadingModel::Distributed1D, num_threads_1d, batchSize >
{
  static constexpr ThreadingModel threading_model = ThreadingModel::Distributed1D;
  static constexpr localIndex num_threads_x = num_threads_1d;
  static constexpr localIndex num_threads_y = 1;
  static constexpr localIndex num_threads_z = batchSize;
  static constexpr localIndex batch_size = batchSize;

  template < typename T, localIndex... Dims >
  using Tensor = tensor::Static2dThreadTensor< T, Dims... >; // FIXME

  template < localIndex num_dofs, localIndex num_quads >
  using Basis = geosx::stackVariables::template SharedBasis< num_dofs, num_quads >;

  // Thread index x
  localIndex tidx;
  // Index inside a batch of elements
  localIndex batch_index;
  // RAJA launch context
  RAJA::LaunchContext & ctx;

  GEOSX_HOST_DEVICE
  KernelConfiguration( RAJA::LaunchContext & ctx )
  : ctx( ctx )
  {
    using RAJA::RangeSegment;
    loop<thread_x>( ctx, RangeSegment( 0, num_threads_1d ), [&]( localIndex const tid_x )
    {
      tidx = tid_x;
    } );
  }
};

// Distributed 2D
template < localIndex num_threads_1d, localIndex batchSize >
struct KernelConfiguration< ThreadingModel::Distributed2D, num_threads_1d, batchSize >
{
  static constexpr ThreadingModel threading_model = ThreadingModel::Distributed2D;
  static constexpr localIndex num_threads_x = num_threads_1d;
  static constexpr localIndex num_threads_y = num_threads_1d;
  static constexpr localIndex num_threads_z = batchSize;
  static constexpr localIndex batch_size = batchSize;

  template < typename T, localIndex... Dims >
  using Tensor = tensor::Static2dThreadTensor< T, Dims... >;

  template < localIndex num_dofs, localIndex num_quads >
  using Basis = geosx::stackVariables::template SharedBasis< num_dofs, num_quads >;

  // Thread index x
  localIndex tidx;
  // Thread index y
  localIndex tidy;
  // Index inside a batch of elements
  localIndex batch_index;
  // RAJA launch context
  RAJA::LaunchContext & ctx;
  // Shared memory buffers, using buffers allows to avoid using too much shared memory.
  static constexpr localIndex dim = 3;
  static constexpr localIndex buffer_size = num_threads_1d * num_threads_1d * num_threads_1d;
  static constexpr localIndex num_buffers = 2 * dim;
  geosx::stackVariables::SharedMemBuffers< buffer_size, num_buffers, batch_size > shared_mem;


  GEOSX_HOST_DEVICE
  KernelConfiguration( RAJA::LaunchContext & ctx )
  : ctx( ctx ), shared_mem( ctx )
  {
    using RAJA::RangeSegment;
    loop<thread_x>( ctx, RangeSegment( 0, num_threads_1d ), [&]( localIndex const tid_x )
    {
      tidx = tid_x;
    } );
    loop<thread_y>( ctx, RangeSegment( 0, num_threads_1d ), [&]( localIndex const tid_y )
    {
      tidy = tid_y;
    } );
  }
};

// Distributed 3D
template < localIndex num_threads_1d, localIndex batchSize >
struct KernelConfiguration< ThreadingModel::Distributed3D, num_threads_1d, batchSize >
{
  static constexpr ThreadingModel threading_model = ThreadingModel::Distributed3D;
  static constexpr localIndex num_threads_x = num_threads_1d * num_threads_1d * num_threads_1d;
  static constexpr localIndex num_threads_y = 1;
  static constexpr localIndex num_threads_z = batchSize;
  static constexpr localIndex batch_size = batchSize;

  template < typename T, localIndex... Dims >
  using Tensor = tensor::Static3dThreadTensor< T, Dims... >;

  template < localIndex num_dofs, localIndex num_quads >
  using Basis = geosx::stackVariables::template SharedBasis< num_dofs, num_quads >;

  // Thread index x
  localIndex tidx;
  // Thread index y
  localIndex tidy;
  // Thread index z
  localIndex tidz;
  // Index inside a batch of elements
  localIndex batch_index;
  // RAJA launch context
  RAJA::LaunchContext & ctx;
  // Shared memory buffers, using buffers allows to avoid using too much shared memory.
  static constexpr localIndex dim = 3;
  static constexpr localIndex buffer_size = num_threads_1d * num_threads_1d * num_threads_1d * dim * dim;
  static constexpr localIndex num_buffers = 1;
  geosx::stackVariables::SharedMemBuffers< buffer_size, num_buffers, batch_size > shared_mem;

  GEOSX_HOST_DEVICE
  KernelConfiguration( RAJA::LaunchContext & ctx )
  : ctx( ctx ), shared_mem( ctx )
  {
    using RAJA::RangeSegment;
    loop<thread_x>( ctx, RangeSegment( 0, num_threads_x ), [&]( localIndex const tid )
    {
      tidx = tid % num_threads_1d;
      tidy = ( tid % ( num_threads_1d * num_threads_1d ) ) / num_threads_1d;
      tidz = tid / ( num_threads_1d * num_threads_1d );
    } );
  }
};

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_KERNEL_CONFIG_HPP_ */
