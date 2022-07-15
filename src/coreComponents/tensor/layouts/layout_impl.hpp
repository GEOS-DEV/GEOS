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
 * @file layout_impl.hpp
 */

#ifndef GEOSX_TENSOR_LAYOUT_IMPL_HPP_
#define GEOSX_TENSOR_LAYOUT_IMPL_HPP_

#include "tensor/utilities/error.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace tensor
{
/// Layout utility classes

//Compute the index inside a StaticLayout
template<int Cpt, int rank, int... Dims>
struct StaticIndex
{
   template <typename... Idx> GEOSX_HOST_DEVICE inline
   static int eval(int first, Idx... args)
   {
      constexpr int size = get_value<Cpt-1,Dims...>;
      GEOSX_ASSERT_KERNEL(
         first<size,
         "Index (%d) greater than the static size (%d).\n", first, size);
      return first + size * StaticIndex<Cpt+1, rank, Dims...>::eval(args...);
   }
};

template<int rank, int... Dims>
struct StaticIndex<rank,rank,Dims...>
{
   GEOSX_HOST_DEVICE inline
   static int eval(int first)
   {
#ifdef GEOSX_DEBUG
      constexpr int size = get_value<rank-1,Dims...>;
      GEOSX_ASSERT_KERNEL(
         size==Dynamic || first<size,
         "Index (%d) greater than the static size (%d).\n", first, size);
#endif
      return first;
   }
};

template<int... Dims>
struct StaticLayoutIndex
{
   template <typename... Idx> GEOSX_HOST_DEVICE inline
   static int eval(Idx... args)
   {
      return StaticIndex<1,sizeof...(Dims),Dims...>::eval(args...);
   }
};

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_LAYOUT_IMPL_HPP_
