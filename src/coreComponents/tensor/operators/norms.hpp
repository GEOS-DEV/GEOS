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
 * @file norms.hpp
 */

#ifndef GEOSX_TENSOR_NORMS
#define GEOSX_TENSOR_NORMS

// #include "../../../general/backends.hpp"
#include "tensor/tensor_traits.hpp"
#include "tensor/utilities/foreach.hpp"

namespace geosx
{

namespace tensor
{

/// Serial version
template <typename Tensor,
          std::enable_if_t<
             is_serial_tensor<Tensor>,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto SquaredNorm(const Tensor& t)
{
   using Scalar = get_tensor_type<Tensor>;

   Scalar norm = 0;

   ForallDims<Tensor>::Apply(t, [&](auto... idx){
      const Scalar& val = t(idx...);
      norm += val*val;
   });

   return norm;
}

/// Threaded version
template <typename Tensor,
          std::enable_if_t<
             !is_serial_tensor<Tensor>,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto SquaredNorm(const Tensor& t)
{
   using Scalar = get_tensor_type<Tensor>;

   GEOSX_SHARED Scalar res;
   if (GEOSX_THREAD_ID(x)==0 && GEOSX_THREAD_ID(y)==0 && GEOSX_THREAD_ID(z)==0)
   {
      res = 0.0;
   }
   GEOSX_SYNC_THREAD;
   Scalar norm = 0;

   ForallDims<Tensor>::Apply(t, [&](auto... idx){
      const Scalar& val = t(idx...);
      norm += val*val;
   });
   AtomicAdd(res, norm);
   GEOSX_SYNC_THREAD;

   return res;
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_NORMS
