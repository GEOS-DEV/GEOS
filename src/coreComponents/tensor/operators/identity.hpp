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
 * @file identity.hpp
 */

#ifndef GEOSX_TENSOR_IDENTITY
#define GEOSX_TENSOR_IDENTITY

#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace tensor
{

struct Identity { };

template <typename Tensor,
          std::enable_if_t<
             is_tensor<Tensor> &&
             !std::is_same<Tensor, ResultTensor<Tensor>>::value,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto operator*(const Identity &I, const Tensor &u)
{
   using Result = ResultTensor<Tensor>;
   return Result(u);
}

template <typename Tensor,
          std::enable_if_t<
             is_tensor<Tensor> &&
             std::is_same<Tensor, ResultTensor<Tensor>>::value,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto operator*(const Identity &I, const Tensor &u)
{
   return u;
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_IDENTITY
