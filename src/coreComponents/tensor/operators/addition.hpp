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
 * @file addition.hpp
 */


#ifndef GEOSX_TENSOR_SCALAR_ADD
#define GEOSX_TENSOR_SCALAR_ADD

#include "tensor/tensor_traits.hpp"
#include "tensor/utilities/foreach.hpp"

namespace geosx
{

namespace tensor
{

template <typename Tensor,
          std::enable_if_t<
             is_tensor<Tensor>,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto operator+(const Tensor &u, const Tensor &v)
{
   using Res = ResultTensor<Tensor>;
   Res w(GetLayout(u));
   ForallDims<Tensor>::ApplyBinOp(u, v, [&](auto... idx)
   {
      w(idx...) = u(idx...) + v(idx...);
   });
   return w;
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_SCALAR_ADD
