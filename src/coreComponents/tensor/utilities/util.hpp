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
 * @file util.hpp
 */

#ifndef GEOSX_TENSOR_UTIL
#define GEOSX_TENSOR_UTIL

#include "error.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
// #include "../../../general/forall.hpp"

namespace geosx
{

namespace tensor
{

/// Getter for the N-th dimension value
template <int N, int... Dims>
struct get_value_v;

template <int Dim0, int... Dims>
struct get_value_v<0, Dim0, Dims...>
{
   static constexpr int value = Dim0;
};

template <int N, int Dim0, int... Dims>
struct get_value_v<N, Dim0, Dims...>
{
   static constexpr int value = get_value_v<N-1,Dims...>::value;
};

template <int N, int... Dims>
constexpr int get_value = get_value_v<N,Dims...>::value;

/// Get the last value
template <typename T> GEOSX_HOST_DEVICE
T GetLast(T first)
{
   return first;
}

template <typename T, typename... Ts> GEOSX_HOST_DEVICE
auto GetLast(T first, Ts... rest)
{
   GEOSX_UNUSED_VAR( first );
   return GetLast(rest...);
}

/// Does the same as sizeof...
template <int... Sizes>
constexpr int rank = sizeof...(Sizes);

/// IfThenElse
template <bool Cond, typename TrueType, typename FalseType>
struct IfThenElse_t
{
   using type = TrueType;
};

template <typename TrueType, typename FalseType>
struct IfThenElse_t<false, TrueType, FalseType>
{
   using type = FalseType;
};

template <bool Cond, typename TrueType, typename FalseType>
using IfThenElse = typename IfThenElse_t<Cond,TrueType,FalseType>::type;

template <typename... Args> GEOSX_HOST_DEVICE
void one_print(const char* msg, Args... vals)
{
#ifdef GEOSX_DEVICE_COMPILE
   if (GEOSX_THREAD_ID(x)==0 && GEOSX_THREAD_ID(y)==0 && GEOSX_THREAD_ID(z)==0)
   {
      printf(msg,vals...);
   }
#else
   printf(msg,vals...);
#endif
}

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_UTIL
