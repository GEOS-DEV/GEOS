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
 * @file restricted_layout.hpp
 */

#ifndef GEOSX_RESTRICTED_LAYOUT_HPP_
#define GEOSX_RESTRICTED_LAYOUT_HPP_

#include "tensor/utilities/error.hpp"
#include "layout_traits.hpp"

namespace geosx
{

namespace tensor
{

// TODO possible to write directly in a generic way?
// TODO throw an error if N>8?
template <int N, typename Layout>
class RestrictedLayout;

template <typename Layout>
class RestrictedLayout<0,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<0>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<0>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(Idx... idx) const
   {
      return layout.index(i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   int Size() const
   {
      return layout.template Size<M+1>();
   }
};

template <typename Layout>
class RestrictedLayout<1,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<1>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<1>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, Idx... idx) const
   {
      return layout.index(idx0,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<1?M:M+1;
      return layout.template Size<Dim>();
   }
};

template <typename Layout>
class RestrictedLayout<2,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<2>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<2>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, int idx1, Idx... idx) const
   {
      return layout.index(idx0,idx1,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<2?M:M+1;
      return layout.template Size<Dim>();
   }
};

template <typename Layout>
class RestrictedLayout<3,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<3>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<3>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, int idx1, int idx2, Idx... idx) const
   {
      return layout.index(idx0,idx1,idx2,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<3?M:M+1;
      return layout.template Size<Dim>();
   }
};

template <typename Layout>
class RestrictedLayout<4,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<4>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<4>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, int idx1, int idx2, int idx3, Idx... idx) const
   {
      return layout.index(idx0,idx1,idx2,idx3,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<4?M:M+1;
      return layout.template Size<Dim>();
   }
};

template <typename Layout>
class RestrictedLayout<5,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<5>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<5>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, int idx1, int idx2, int idx3, int idx4, Idx... idx) const
   {
      return layout.index(idx0,idx1,idx2,idx3,idx4,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<5?M:M+1;
      return layout.template Size<Dim>();
   }
};

template <typename Layout>
class RestrictedLayout<6,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<6>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<6>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, int idx1, int idx2, int idx3, int idx4, int idx5, Idx... idx) const
   {
      return layout.index(idx0,idx1,idx2,idx3,idx4,idx5,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<6?M:M+1;
      return layout.template Size<Dim>();
   }
};

template <typename Layout>
class RestrictedLayout<7,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<7>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<7>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, int idx1, int idx2, int idx3, int idx4, int idx5, int idx6, Idx... idx) const
   {
      return layout.index(idx0,idx1,idx2,idx3,idx4,idx5,idx6,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<7?M:M+1;
      return layout.template Size<Dim>();
   }
};

template <typename Layout>
class RestrictedLayout<8,Layout>
{
private:
   const int i;
   const Layout &layout;

public:
   GEOSX_HOST_DEVICE
   RestrictedLayout(int i, const Layout &layout): i(i), layout(layout)
   {
      GEOSX_ASSERT_KERNEL(
         i<layout.template Size<8>(),
         "The RestrictedLayout index (%d) is out of bounds (%d).\n",
         i, layout.template Size<8>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   constexpr int index(int idx0, int idx1, int idx2, int idx3,int idx4, int idx5, int idx6, int idx7, Idx... idx) const
   {
      return layout.index(idx0,idx1,idx2,idx3,idx4,idx5,idx6,idx7,i,idx...);
   }

   template <int M> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Dim = M<8?M:M+1;
      return layout.template Size<Dim>();
   }
};

// get_layout_rank
template <int I, typename Layout>
struct get_layout_rank_v<RestrictedLayout<I,Layout>>
{
   static constexpr int value = get_layout_rank<Layout> - 1;
};

// is_dynamic_layout
template <int N, typename Layout>
struct is_dynamic_layout_v<RestrictedLayout<N,Layout>>
{
   static constexpr bool value = is_dynamic_layout<Layout>;
};

// is_static_layout
template <int N, typename Layout>
struct is_static_layout_v<RestrictedLayout<N,Layout>>
{
   static constexpr bool value = is_static_layout<Layout>;
};

// is_serial_layout
template <int N, typename Layout>
struct is_serial_layout_v<RestrictedLayout<N,Layout>>
{
   static constexpr bool value = is_serial_layout<Layout>;
};

// is_2d_threaded_layout
template <int N, typename Layout>
struct is_2d_threaded_layout_v<RestrictedLayout<N,Layout>>
{
   static constexpr bool value = is_2d_threaded_layout<Layout>;
};

// is_3d_threaded_layout
template <int N, typename Layout>
struct is_3d_threaded_layout_v<RestrictedLayout<N,Layout>>
{
   static constexpr bool value = is_3d_threaded_layout<Layout>;
};

// is_serial_layout_dim
template <int N, int R, typename Layout>
struct is_serial_layout_dim_v<RestrictedLayout<R,Layout>, N>
{
   static constexpr bool value = N<R?
                                 is_serial_layout_dim<Layout,N>:
                                 is_serial_layout_dim<Layout,N+1>;
};

// is_threaded_layout_dim
template <int N, int R, typename Layout>
struct is_threaded_layout_dim_v<RestrictedLayout<R,Layout>, N>
{
   static constexpr bool value = N<R?
                                 is_threaded_layout_dim<Layout,N>:
                                 is_threaded_layout_dim<Layout,N+1>;
};

// get_layout_size
template <int N, int I, typename Layout>
struct get_layout_size_v<N,RestrictedLayout<I,Layout>>
{
   static constexpr int value = get_layout_size<N+(N>=I),Layout>;
};

// get_layout_sizes
template <int N, typename Layout>
struct get_layout_sizes_t<RestrictedLayout<N,Layout>>
{
   using type = remove< N, get_layout_sizes<Layout> >;
};

// get_layout_batch_size
template <int N, typename Layout>
struct get_layout_batch_size_v<RestrictedLayout<N, Layout>>
{
   static constexpr int value = get_layout_batch_size<Layout>;
};

// get_layout_capacity
template <int N, typename Layout>
struct get_layout_capacity_v<RestrictedLayout<N,Layout>>
{
   static constexpr int capacity = get_layout_capacity<Layout>;
   static constexpr int sizeN = get_layout_size<N,Layout>;
   static constexpr int value = sizeN != Dynamic ?
                                         ( capacity / sizeN) :
                                         Dynamic;
};

// get_layout_result_type
template <int N, typename Layout>
struct get_layout_result_type_t< RestrictedLayout<N,Layout> >
{
   template <int... Sizes>
   using type = get_layout_result_type<Layout,Sizes...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_RESTRICTED_LAYOUT_HPP_
