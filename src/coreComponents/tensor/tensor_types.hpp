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
 * @file tensor_types.hpp
 */

#ifndef GEOSX_TENSOR_TYPES_HPP_
#define GEOSX_TENSOR_TYPES_HPP_

#include "tensor.hpp"
/// The memory storage for the tensors' values
#include "containers/containers.hpp"
/// The rank N index mapping to linear memory index
#include "layouts/layouts.hpp"
/// Implementation details
#include "utilities/tensor_types_impl.hpp"

namespace geosx
{

namespace tensor
{

/****************************************************************************
 *              Behavioral tensor types.
 * This file defines types of tensors with different behaviors based on their
 * type of container, and their type of layout.
 * Their list is:
 *    - DynamicTensor,
 *    - StaticTensor,
 *    - Dynamic1dThreadTensor,
 *    - Static1dThreadTensor,
 *    - Dynamic2dThreadTensor,
 *    - Static2dThreadTensor,
 *    - Dynamic3dThreadTensor,
 *    - Static3dThreadTensor,
 *    - ThreadTensor,
 *    - MyDeviceTensor,
 *    - StaticPointerTensor,
 *    - ReadTensor.
 * */

/// Dynamically sized tensor
/** DyncamicTensor represent stack allocated tensors with sizes known at runtime.
    ex: `DynamicTensor<4,real64> u(2,3,4,5);',
    represents a dynamically sized stack allocated rank 4 tensor with respective
    dimensions 2,3,4, and 5. The memory stack allocated being
    `pow(DynamicMaxSize,Rank)'.
   */
template <int Rank, typename T, int MaxSize = pow(DynamicMaxSize,Rank)>
using DynamicTensor = TensorBase<StackContainer<T, MaxSize>,
                             DynamicLayout<Rank> >;

/// Helper type for DynamicTensor with real64.
template <int Rank, int MaxSize = pow(DynamicMaxSize,Rank)>
using DynamicDTensor = DynamicTensor<Rank,real64,MaxSize>;

/// Statically sized Tensor
/** StaticTensor represent stack allocated tensors with dimensions known at the
    compilation.
    ex: `StaticTensor<real64,2,3,4,5> u;',
    represents a stack allocated rank 4 tensor with compilation time dimensions
    2, 3, 4, and 5.
    These tensors have the propriety to be thread private on GPU.
   */
template <typename T, int... Sizes>
using StaticTensor = TensorBase<StackContainer<T, prod(Sizes...)>,
                            StaticLayout<Sizes...> >;

/// Helper type for StaticTensor with real64
template <int... Sizes>
using StaticDTensor = StaticTensor<real64,Sizes...>;

/// A Tensor dynamically distributed over threads
/** Dynamic1dThreadTensor represent stack allocated tensors whith dimensions
    known at runtime, their data is distributed over threads (ThreadIdx.x).
    Shared memory MUST be used to read data located in a different thread.
   */
template <int Rank, typename T, int BatchSize, int MaxSize = get_Dynamic1dThreadLayout_size(DynamicMaxSize,Rank)>
using Dynamic1dThreadTensor = TensorBase<
                              StackContainer<
                                 T,
                                 MaxSize
                              >,
                              Dynamic1dThreadLayout<Rank,BatchSize>
                           >;

/// Helper type for Dynamic1dThreadTensor with real64
template <int Rank, int BatchSize, int MaxSize = DynamicMaxSize>
using Dynamic1dThreadDTensor = Dynamic1dThreadTensor<Rank,real64,BatchSize,MaxSize>;

/// A Tensor statically distributed over threads
/** Static1dThreadTensor represent stack allocated tensors whith dimensions
    known at compilation, their data is distributed over threads (ThreadIdx.x).
    Shared memory MUST be used to read data located in a different thread.
   */
template <typename T, int BatchSize, int... Sizes>
using Static1dThreadTensor = TensorBase<StackContainer<
                                       T,
                                       get_Static1dThreadTensor_size(Sizes...)
                                    >,
                                    Static1dThreadLayout<BatchSize, Sizes...>
                             >;

/// Helper type for Static1dThreadTensor with real64
template <int BatchSize, int... Sizes>
using Static1dThreadDTensor = Static1dThreadTensor<real64,BatchSize,Sizes...>;

/// A Tensor dynamically distributed over a plane of threads
/** Dynamic2dThreadTensor represent stack allocated tensors whith dimensions
    known at runtime, their data is distributed over a plane of threads
    (ThreadIdx.x and ThreadIdx.y).
    Shared memory MUST be used to read data located in a different thread.
   */
template <int Rank, typename T, int BatchSize, int MaxSize = get_Dynamic2dThreadLayout_size(DynamicMaxSize,Rank)>
using Dynamic2dThreadTensor = TensorBase<
                              StackContainer<
                                 T,
                                 MaxSize
                              >,
                              Dynamic2dThreadLayout<Rank,BatchSize>
                           >;

/// Helper type for Dynamic2dThreadTensor with real64
template <int Rank, int BatchSize, int MaxSize = DynamicMaxSize>
using Dynamic2dThreadDTensor = Dynamic2dThreadTensor<Rank,real64,BatchSize,MaxSize>;

/// A Tensor statically distributed over a plane of threads
/** Static2dThreadTensor represent stack allocated tensors whith dimensions
    known at compilation, their data is distributed over a plane of threads
    (ThreadIdx.x and ThreadIdx.y).
    Shared memory MUST be used to read data located in a different thread.
   */
template <typename T, int BatchSize, int... Sizes>
using Static2dThreadTensor = TensorBase<StackContainer<
                                       T,
                                       get_Static2dThreadTensor_size(Sizes...)
                                    >,
                                    Static2dThreadLayout<BatchSize, Sizes...>
                             >;

/// Helper type for Static2dThreadTensor with real64
template <int BatchSize, int... Sizes>
using Static2dThreadDTensor = Static2dThreadTensor<real64,BatchSize,Sizes...>;

/// A Tensor dynamically distributed over a cube of threads
/** Dynamic3dThreadTensor represent stack allocated tensors whith dimensions
    known at runtime, their data is distributed over a cube of threads
    (ThreadIdx.x, ThreadIdx.y, and ThreadIdx.z).
    Shared memory MUST be used to read data located in a different thread.
   */
template <int Rank, typename T, int BatchSize, int MaxSize = get_Dynamic3dThreadLayout_size(DynamicMaxSize,Rank)>
using Dynamic3dThreadTensor = TensorBase<
                                 StackContainer<
                                    T,
                                    MaxSize
                                 >,
                                 Dynamic3dThreadLayout<Rank,BatchSize>
                              >;

/// Helper type for Dynamic3dThreadTensor with real64
template <int Rank, int BatchSize, int MaxSize = DynamicMaxSize>
using Dynamic3dThreadDTensor = Dynamic3dThreadTensor<Rank,real64,BatchSize,MaxSize>;

/// A Tensor statically distributed over a cube of threads
/** Static3dThreadTensor represent stack allocated tensors whith dimensions
    known at compilation, their data is distributed over a cube of threads
    (ThreadIdx.x, ThreadIdx.y, and ThreadIdx.z).
    Shared memory MUST be used to read data located in a different thread.
   */
template <typename T, int BatchSize, int... Sizes>
using Static3dThreadTensor = TensorBase<StackContainer<
                                       T,
                                       get_Static3dThreadTensor_size(Sizes...)
                                    >,
                                    Static3dThreadLayout<BatchSize, Sizes...>
                             >;

template <int BatchSize, int... Sizes>
using Static3dThreadDTensor = Static3dThreadTensor<real64,BatchSize,Sizes...>;

/// A threaded tensor type with Dim threaded dimensions
/** ThreadTensor<Dim> contains static_type, and dynamic_type corresponding to
    Static(Dim)ThreadTensor, and Dynamic(Dim)ThreadTensor.
*/
template <int Dim>
struct ThreadTensor;

/// Dim == 1
template <>
struct ThreadTensor<1>
{
   template <typename T, int BatchSize, int... Sizes>
   using static_type = Static2dThreadTensor<T,BatchSize,Sizes...>;

   template <int Rank, typename T, int BatchSize, int MaxSize = DynamicMaxSize>
   using dynamic_type = Dynamic2dThreadTensor<Rank,T,BatchSize,MaxSize>;
};

/// Dim == 2
template <>
struct ThreadTensor<2>
{
   template <typename T, int BatchSize, int... Sizes>
   using static_type = Static2dThreadTensor<T,BatchSize,Sizes...>;

   template <int Rank, typename T, int BatchSize, int MaxSize = DynamicMaxSize>
   using dynamic_type = Dynamic2dThreadTensor<Rank,T,BatchSize,MaxSize>;
};

/// Dim == 3
template <>
struct ThreadTensor<3>
{
   template <typename T, int BatchSize, int... Sizes>
   using static_type = Static3dThreadTensor<T,BatchSize,Sizes...>;

   template <int Rank, typename T, int BatchSize, int MaxSize = DynamicMaxSize>
   using dynamic_type = Dynamic3dThreadTensor<Rank,T,BatchSize,MaxSize>;
};

/// Statically sized threaded tensors where Dim is the number of threaded
/// dimensions.
template <int Dim, typename T, int BatchSize, int... Sizes>
using StaticThreadTensor = typename ThreadTensor<Dim>
                              ::template static_type<T,BatchSize,Sizes...>;

/// Dynamically sized threaded tensors where Dim is the number of threaded
/// dimensions.
template <int Dim, int Rank, typename T, int BatchSize, int MaxSize = DynamicMaxSize>
using DynamicThreadTensor = typename ThreadTensor<Dim>
                               ::template dynamic_type<Rank,T,BatchSize,MaxSize>;


/// A tensor using a read write access pointer and a dynamic data layout.
/// Note: Backward compatible if renamed in DeviceTensor.
template <int Rank, typename T>
using MyDeviceTensor = TensorBase<DeviceContainer<T>,
                              DynamicLayout<Rank> >;

template <int Rank>
using DeviceDTensor = MyDeviceTensor<Rank,real64>;

/// Dynamically sized tensors using a pointer as container.
template <int Rank, typename T>
using DynamicPointerTensor = TensorBase<DeviceContainer<T>,
                                    DynamicLayout<Rank> >;

template <int Rank>
using DynamicPointerDTensor = DynamicPointerTensor<Rank,real64>;

/// Statically sized tensor using a pointer container.
template <typename T, int... Sizes>
using StaticPointerTensor = TensorBase<DeviceContainer<T>,
                                   StaticLayout<Sizes...> >;

template <int... Sizes>
using StaticPointerDTensor = StaticPointerTensor<real64,Sizes...>;

/// A tensor using a read only const pointer and a dynamic data layout.
template <int Rank, typename T>
using ReadTensor = TensorBase<ReadContainer<T>,
                          DynamicLayout<Rank> >;

template <int Rank>
using ReadDTensor = ReadTensor<Rank,real64>;

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_TYPES_HPP_
