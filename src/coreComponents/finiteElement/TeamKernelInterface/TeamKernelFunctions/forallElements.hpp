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
 * @file forallElements.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_FORALL_ELEMENTS_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_FORALL_ELEMENTS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_types.hpp"
#include "tensor/tensor_traits.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

// @brief Iterator on the elements of TeamKernelBase object. TODO descrine how to use this
template < typename KernelConfig, typename KernelComponents, typename Lambda >
void forallElements( localIndex const numElems, KernelComponents const & fields, Lambda && element_kernel )
{
  constexpr localIndex num_threads_x = KernelConfig::num_threads_x;
  constexpr localIndex num_threads_y = KernelConfig::num_threads_y;
  constexpr localIndex num_threads_z = KernelConfig::num_threads_z;
  constexpr localIndex batch_size = KernelConfig::batch_size;
  localIndex const num_batches = ( numElems + batch_size - 1 ) / batch_size;
  // localIndex const num_SM = 80;
  // localIndex const num_blocks = 64 * num_SM; //( numElems + batch_size - 1 ) / batch_size;
  localIndex const num_blocks = num_batches;

  launch< team_launch_policy >
  ( GEOSX_RAJA_DEVICE, Grid( Teams( num_blocks ), Threads( num_threads_x, num_threads_y, num_threads_z ) ),
  [=] GEOSX_HOST_DEVICE ( LaunchContext ctx )
  {
    using RAJA::RangeSegment;
    typename KernelComponents::template StackVariables<KernelConfig> stack( ctx );

    // Each block of threads treats "batch_size" elements.
    loop<team_x>( ctx, RangeSegment( 0, num_batches ), [&]( localIndex const & block_index )
    {
      // We batch elements over the z-thread dimension
      loop<thread_z>( ctx, RangeSegment( 0, batch_size ), [&]( localIndex const & batch_index )
      {
        localIndex const element_index = block_index * batch_size + batch_index;
        if ( element_index < numElems )
        {
          stack.batch_index = batch_index;
          stack.element_index = element_index;

          element_kernel( stack );
        }
      } );
    } );
  } );
  parallelDeviceSync();
}

// TODO document
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

  /// Index inside a batch of elements
  localIndex batch_index;
  /// RAJA launch context
  LaunchContext & ctx;

  GEOSX_HOST_DEVICE
  KernelConfiguration( LaunchContext & ctx )
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

  /// Thread index x
  localIndex tidx;
  /// Index inside a batch of elements
  localIndex batch_index;
  /// RAJA launch context
  LaunchContext & ctx;

  GEOSX_HOST_DEVICE
  KernelConfiguration( LaunchContext & ctx )
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

  /// Thread index x
  localIndex tidx;
  /// Thread index y
  localIndex tidy;
  /// Index inside a batch of elements
  localIndex batch_index;
  /// RAJA launch context
  LaunchContext & ctx;

  GEOSX_HOST_DEVICE
  KernelConfiguration( LaunchContext & ctx )
  : ctx( ctx )
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

  /// Thread index x
  localIndex tidx;
  /// Thread index y
  localIndex tidy;
  /// Thread index z
  localIndex tidz;
  /// Index inside a batch of elements
  localIndex batch_index;
  /// RAJA launch context
  LaunchContext & ctx;

  GEOSX_HOST_DEVICE
  KernelConfiguration( LaunchContext & ctx )
  : ctx( ctx )
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

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_FORALL_ELEMENTS_HPP_ */
