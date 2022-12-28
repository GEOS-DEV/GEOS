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
 * @file static_layout.hpp
 */

#ifndef GEOSX_STATIC_LAYOUT_HPP_
#define GEOSX_STATIC_LAYOUT_HPP_

#include "layout_traits.hpp"
#include "layout_impl.hpp"

namespace geosx
{

namespace tensor
{

/// A static layout
template <int... Sizes>
class StaticLayout
{
public:
   GEOSX_HOST_DEVICE
   constexpr StaticLayout() { }

   template <typename... Dims> GEOSX_HOST_DEVICE
   constexpr StaticLayout(Dims... args)
   {
      GEOSX_UNUSED_VAR(args...);
      static_assert(
         sizeof...(Dims)==sizeof...(Sizes),
         "Static and dynamic sizes don't match.");
      // TODO verify that Dims == sizes in Debug mode
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   constexpr StaticLayout(const Layout& rhs)
   {
      GEOSX_UNUSED_VAR(rhs);
      // for (int i = 0; i < Rank; i++)
      // {
      //    GEOSX_ASSERT(Sizes...[i] == lhs.Size<i>());
      // }
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(Idx... idx) const
   {
      static_assert(
         sizeof...(Sizes)==sizeof...(Idx),
         "Wrong number of arguments.");
      return StaticLayoutIndex<Sizes...>::eval(idx...);
   }

   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      static_assert(
         N>=0 && N<sizeof...(Sizes),
         "Accessed size is higher than the rank of the Tensor.");
      return get_value<N,Sizes...>;
   }
};

// get_layout_rank
template <int... Dims>
struct get_layout_rank_v<StaticLayout<Dims...>>
{
   static constexpr int value = sizeof...(Dims);
};

// is_static_layout
template<int... Dims>
struct is_static_layout_v<StaticLayout<Dims...>>
{
   static constexpr bool value = true;
};

// is_serial_layout
template<int... Dims>
struct is_serial_layout_v<StaticLayout<Dims...>>
{
   static constexpr bool value = true;
};

// get_layout_size
template <int N, int... Dims>
struct get_layout_size_v<N, StaticLayout<Dims...>>
{
   static constexpr int value = get_value<N, Dims...>;
};

// get_layout_sizes
template <int... Dims>
struct get_layout_sizes_t<StaticLayout<Dims...>>
{
   using type = int_list<Dims...>;
};

// get_layout_capacity
template <int... Sizes>
struct get_layout_capacity_v<StaticLayout<Sizes...>>
{
   static constexpr int value = prod(Sizes...);
};

// get_layout_result_type
template <int... Sizes>
struct get_layout_result_type_t<StaticLayout<Sizes...>>
{
   template <int... Dims>
   using type = StaticLayout<Dims...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_STATIC_LAYOUT_HPP_
