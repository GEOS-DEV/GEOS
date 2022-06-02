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
 * @file KernelBase.hpp
 */

#ifndef GEOSX_TENSOR_SCALAR_MULT
#define GEOSX_TENSOR_SCALAR_MULT

#include "tensor/tensor_traits.hpp"
#include "tensor/utilities/foreach.hpp"
#include "tensor/utilities/get_layout.hpp"

namespace geosx
{

namespace tensor
{

template <typename Scalar,
          typename Tensor,
          std::enable_if_t<
             is_tensor<Tensor> &&
             std::is_same<Scalar, get_tensor_type<Tensor>>::value,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto operator*(const Scalar &a, const Tensor &u)
{
   using Res = ResultTensor<Tensor>;
   Res v(GetLayout(u));
   ForallDims<Tensor>::ApplyBinOp(u, v, [&](auto... idx)
   {
      v(idx...) = a * u(idx...);
   });
   return v;
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_SCALAR_MULT
