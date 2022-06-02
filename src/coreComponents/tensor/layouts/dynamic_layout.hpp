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
 * @file dynamic_layout.hpp
 */

#ifndef GEOSX_DYNAMIC_LAYOUT_HPP_
#define GEOSX_DYNAMIC_LAYOUT_HPP_

#include "tensor/utilities/helper_constants.hpp"
#include "layout_impl.hpp"
#include "tensor/utilities/error.hpp"
#include "layout_traits.hpp"

namespace geosx
{

namespace tensor
{

template <int... Sizes>
class SizedDynamicLayout;

template <int Rank>
using DynamicLayout = instantiate<SizedDynamicLayout,int_repeat<Dynamic,Rank>>;

/// A dynamic layout with first index fastest
template <int... Sizes>
class SizedDynamicLayout
{
private:
   static constexpr int Rank = sizeof...(Sizes);
   int sizes[Rank];

public:
   template <typename... Dims> GEOSX_HOST_DEVICE
   SizedDynamicLayout(int arg0, Dims... args)
   {
      // TODO Add Sizes checks
      InitDynamicLayout<Rank>::result(sizes, arg0, args...);
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamicLayout(const Layout &rhs)
   {
      // TODO Add Sizes checks
      InitDynamicLayout<Rank>::result(sizes,rhs);
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   int index(Idx... idx) const
   {
      static_assert(
         Rank==sizeof...(Idx),
         "Wrong number of arguments.");
      return DynamicLayoutIndex<Rank>::eval(sizes, idx...);
   }

   template <int N> GEOSX_HOST_DEVICE inline
   int Size() const
   {
      static_assert(
         N>=0 && N<Rank,
         "Accessed size is higher than the rank of the Tensor.");
      return sizes[N];
   }
};

// get_layout_rank
template <int... Sizes>
struct get_layout_rank_v<SizedDynamicLayout<Sizes...>>
{
   static constexpr int value = sizeof...(Sizes);
};

// is_dynamic_layout
template <int... Sizes>
struct is_dynamic_layout_v<SizedDynamicLayout<Sizes...>>
{
   static constexpr bool value = true;
};

// is_serial_layout
template <int... Sizes>
struct is_serial_layout_v<SizedDynamicLayout<Sizes...>>
{
   static constexpr bool value = true;
};

// get_layout_size
template <int N, int... Sizes>
struct get_layout_size_v<N, SizedDynamicLayout<Sizes...>>
{
   static constexpr int value = get_value<N,Sizes...>;
};

// get_layout_sizes
template <int... Sizes>
struct get_layout_sizes_t<SizedDynamicLayout<Sizes...>>
{
   using type = int_list<Sizes...>;//int_repeat<Dynamic,Rank>;
};

// get_layout_capacity
template <int First, int... Rest>
struct get_layout_capacity_v<SizedDynamicLayout<First,Rest...>>
{
   static constexpr int value = (First==Dynamic? DynamicMaxSize : First) *
                                get_layout_capacity_v<SizedDynamicLayout<Rest...>>::value;
};

template <int Size>
struct get_layout_capacity_v<SizedDynamicLayout<Size>>
{
   static constexpr int value = (Size==Dynamic? DynamicMaxSize : Size);
};

// get_layout_result_type
template <int... Sizes>
struct get_layout_result_type_t<SizedDynamicLayout<Sizes...>>
{
   template <int... mySizes>
   using type = SizedDynamicLayout<mySizes...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_DYNAMIC_LAYOUT_HPP_
