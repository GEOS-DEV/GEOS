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
 * @file dot_product.hpp
 */

#ifndef GEOSX_TENSOR_MAT_MULT
#define GEOSX_TENSOR_MAT_MULT

// #include "../../../general/backends.hpp"
#include "tensor/tensor_traits.hpp"
#include "tensor/utilities/foreach.hpp"

namespace geosx
{

namespace tensor
{

/// Serial version
template <typename RHS,
          typename LHS,
          std::enable_if_t<
             get_tensor_rank<RHS> == get_tensor_rank<LHS> &&
             is_serial_tensor<RHS> &&
             is_serial_tensor<LHS>,
             bool> = true >
GEOSX_HOST_DEVICE inline
auto Dot(const LHS &lhs, const RHS &rhs)
{
   using Scalar = get_tensor_type<RHS>;

   Scalar res = 0;

   ForallDims<RHS>::ApplyBinOp(lhs, rhs, [&](auto... idx){
      res += lhs(idx...)*rhs(idx...);
   });

   return res;
}

/// Threaded version
template <typename RHS,
          typename LHS,
          std::enable_if_t<
             get_tensor_rank<RHS> == get_tensor_rank<LHS> &&
             (!is_serial_tensor<RHS> ||
              !is_serial_tensor<LHS>),
             bool> = true >
GEOSX_HOST_DEVICE inline
auto Dot(const LHS &lhs, const RHS &rhs)
{
   using Scalar = get_tensor_type<RHS>;

   GEOSX_SHARED Scalar res;
   if (GEOSX_THREAD_ID(x)==0 && GEOSX_THREAD_ID(y)==0 && GEOSX_THREAD_ID(z)==0)
   {
      res = 0.0;
   }
   GEOSX_SYNC_THREAD;

   Scalar loc_res = 0.0;
   ForallDims<RHS>::ApplyBinOp(lhs, rhs, [&](auto... idx){
      loc_res += lhs(idx...)*rhs(idx...);
   });
   AtomicAdd(res, loc_res);
   GEOSX_SYNC_THREAD;

   return res;
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_MAT_MULT
