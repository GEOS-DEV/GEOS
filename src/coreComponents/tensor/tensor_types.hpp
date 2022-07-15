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

namespace geosx
{

namespace tensor
{

/****************************************************************************
 *              Behavioral tensor types.
 * This file defines types of tensors with different behaviors based on their
 * type of container, and their type of layout.
 * */

/// Statically sized Tensor
/** StaticTensor represent stack allocated tensors with dimensions known at the
    compilation.
    ex: `StaticTensor<real64,2,3,4,5> u;',
    represents a stack allocated rank 4 tensor with compilation time dimensions
    2, 3, 4, and 5.
    These tensors have the propriety to be thread private on GPU.
   */
template <typename T, int... Sizes>
using StaticTensor = TensorBase<
                        StackContainer<T, prod(Sizes...)>,
                        StaticLayout<Sizes...> >;

/// Helper type for StaticTensor with real64
template <int... Sizes>
using StaticDTensor = StaticTensor<real64,Sizes...>;

/// Statically sized tensor using a pointer container.
template <typename T, int... Sizes>
using StaticPointerTensor = TensorBase<
                              PointerContainer<T>,
                              StaticLayout<Sizes...> >;

template <int... Sizes>
using StaticPointerDTensor = StaticPointerTensor<real64,Sizes...>;

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_TYPES_HPP_
