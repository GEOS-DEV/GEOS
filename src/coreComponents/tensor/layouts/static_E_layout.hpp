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
 * @file static_E__layout.hpp
 */

#ifndef GEOSX_STATIC_E_LAYOUT_HPP_
#define GEOSX_STATIC_E_LAYOUT_HPP_

#include "layout_traits.hpp"
#include "restricted_layout.hpp"

namespace geosx
{

namespace tensor
{

/// A static layout with the last dimension dynamic (typically for element index)
template<int... Sizes>
class StaticELayout
{
private:
   const int last_size;

public:
   template <typename... Dims> GEOSX_HOST_DEVICE
   StaticELayout(int arg0, Dims... args)
   : last_size(GetLast(arg0, args...))
   {
      static_assert(sizeof...(Dims)==sizeof...(Sizes),
         "Static and dynamic sizes don't match.");
      // TODO verify that Dims == sizes in Debug mode
   }

   // GEOSX_HOST_DEVICE
   // constexpr StaticELayout(const StaticELayout& rhs)
   // : last_size(rhs.last_size)
   // {
   //    // for (int i = 0; i < Rank; i++)
   //    // {
   //    //    GEOSX_ASSERT(Sizes...[i] == lhs.Size<i>());
   //    // }
   // }

   // template <typename Layout> GEOSX_HOST_DEVICE
   // constexpr StaticELayout(const Layout& rhs)
   // : last_size(rhs.template Size<sizeof...(Sizes)>())
   // {
   //    // for (int i = 0; i < Rank; i++)
   //    // {
   //    //    GEOSX_ASSERT(Sizes...[i] == lhs.Size<i>());
   //    // }
   // }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   int index(Idx... idx) const
   {
      static_assert(
         sizeof...(Sizes)+1==sizeof...(Idx),
         "Wrong number of arguments.");
      GEOSX_ASSERT_KERNEL(
         GetLast(idx...) < last_size,
         "Last index (%d) is out of bound (%d).\n",
         GetLast(idx...), last_size);
      return StaticELayoutIndex<Sizes...>::eval(idx...);
   }

   template <int N> GEOSX_HOST_DEVICE inline
   int Size() const
   {
      static_assert(
         N>=0 && N<sizeof...(Sizes)+1,
         "Accessed size is higher than the rank of the Tensor.");
      return StaticELayoutSize<sizeof...(Sizes),N,Sizes...>::eval(last_size);
   }
};

// get_layout_rank
template <int... Dims>
struct get_layout_rank_v<StaticELayout<Dims...>>
{
   static constexpr int value = sizeof...(Dims)+1;
};

// is_static_layout
template<int... Dims>
struct is_static_layout_v<StaticELayout<Dims...>>
{
   static constexpr bool value = true;
};

// is_serial_layout
template<int... Dims>
struct is_serial_layout_v<StaticELayout<Dims...>>
{
   static constexpr bool value = true;
};

// get_layout_size
template <int N, int... Dims>
struct get_layout_size_v<N, StaticELayout<Dims...>>
{
   static constexpr int value = get_value<N, Dims...>;
};

// get_layout_sizes
template <int... Dims>
struct get_layout_sizes_t<StaticELayout<Dims...>>
{
   using type = int_list<Dims...,Dynamic>;
};

// get_layout_capacity
template <int... Sizes>
struct get_layout_capacity_v<
   RestrictedLayout<sizeof...(Sizes),StaticELayout<Sizes...>>>
{
   static constexpr int value = prod(Sizes...);
};

// get_layout_result_type
template <int... Sizes>
struct get_layout_result_type_t<StaticELayout<Sizes...>>
{
   template <int... Dims>
   using type = StaticLayout<Dims...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_STATIC_E_LAYOUT_HPP_
