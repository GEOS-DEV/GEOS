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
 * @file dynamic_2dthread_layout.hpp
 */

#ifndef GEOSX_DYNAMIC_2DTHREAD_LAYOUT_HPP_
#define GEOSX_DYNAMIC_2DTHREAD_LAYOUT_HPP_

#include "tensor/utilities/error.hpp"
#include "dynamic_layout.hpp"
#include "layout_traits.hpp"
#include "tensor/operators/get.hpp"

namespace geosx
{

namespace tensor
{

template <int BatchSize, int... Sizes>
class SizedDynamic2dThreadLayout;

template <int Rank, int BatchSize>
using Dynamic2dThreadLayout = instantiate<SizedDynamic2dThreadLayout,
                                          append<
                                             int_list<BatchSize>,
                                             int_repeat<Dynamic,Rank> > >;

template <int BatchSize, int FirstSize, int SecondSize, int... Sizes>
class SizedDynamic2dThreadLayout<BatchSize,FirstSize,SecondSize,Sizes...>
{
private:
   const int size0;
   const int size1;
   SizedDynamicLayout<Sizes...> layout;
public:
   template <typename... Ss> GEOSX_HOST_DEVICE inline
   SizedDynamic2dThreadLayout(int size0, int size1,  Ss... sizes)
   : size0(size0), size1(size1), layout(sizes...)
   {
      GEOSX_ASSERT_KERNEL(
         FirstSize==Dynamic || FirstSize==size0,
         "Compilation time (%d) and runtime sizes (%d) must be the same.\n",
         FirstSize, size0);
      GEOSX_ASSERT_KERNEL(
         SecondSize==Dynamic || SecondSize==size1,
         "Compilation time (%d) and runtime sizes (%d) must be the same.\n",
         SecondSize, size1);
      GEOSX_ASSERT_KERNEL(
         size0<=GEOSX_THREAD_SIZE(x),
         "The first dimension (%d) exceeds the number of x threads (%d).\n",
         size0, GEOSX_THREAD_SIZE(x));
      GEOSX_ASSERT_KERNEL(
         size1<=GEOSX_THREAD_SIZE(y),
         "The second dimension (%d) exceeds the number of x threads (%d).\n",
         size1, GEOSX_THREAD_SIZE(y));
      GEOSX_ASSERT_KERNEL(
         BatchSize==GEOSX_THREAD_SIZE(z),
         "The batchsize (%d) is not equal to the number of z threads (%d).\n",
         BatchSize, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamic2dThreadLayout(const Layout &rhs)
   : size0(rhs.template Size<0>()),
     size1(rhs.template Size<1>()),
     layout( Get<0>( 0, Get<0>( 0, rhs ) ) )
   {
      constexpr int Rank = sizeof...(Sizes) + 2;
      static_assert(
         Rank == get_layout_rank<Layout>,
         "Can't copy-construct a layout of different rank.");
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   int index(int idx0, int idx1, Idx... idx) const
   {
      GEOSX_UNUSED_VAR(idx0, idx1);
      GEOSX_ASSERT_KERNEL(
         idx0==GEOSX_THREAD_ID(x),
         "The first index (%d) must be equal to the x thread index (%d)"
         " when using SizedDynamic2dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx0, GEOSX_THREAD_ID(x));
      GEOSX_ASSERT_KERNEL(
         idx1==GEOSX_THREAD_ID(y),
         "The second index (%d) must be equal to the y thread index (%d)"
         " when using SizedDynamic2dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx1, GEOSX_THREAD_ID(y));
      return layout.index(idx...);
   }

   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Rank = sizeof...(Sizes) + 2;
      static_assert(
         N>=0 && N<Rank,
         "Accessed size is higher than the rank of the Tensor.");
      return DynamicBlockLayoutSize<N,Rank>::eval(size0,size1,layout);
   }
};

template <int BatchSize, int FirstSize>
class SizedDynamic2dThreadLayout<BatchSize,FirstSize>
{
private:
   const int size0;
public:
   GEOSX_HOST_DEVICE inline
   SizedDynamic2dThreadLayout(int size0)
   : size0(size0)
   {
      GEOSX_ASSERT_KERNEL(
         FirstSize==Dynamic || FirstSize==size0,
         "Compilation time (%d) and runtime sizes (%d) must be the same.\n",
         FirstSize, size0);
      GEOSX_ASSERT_KERNEL(
         size0<=GEOSX_THREAD_SIZE(x),
         "The first dimension (%d) exceeds the number of x threads (%d).\n",
         size0, GEOSX_THREAD_SIZE(x));
      GEOSX_ASSERT_KERNEL(
         BatchSize==GEOSX_THREAD_SIZE(z),
         "The batchsize (%d) is not equal to the number of z threads (%d).\n",
         BatchSize, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamic2dThreadLayout(const Layout &rhs)
   : size0(rhs.template Size<0>())
   {
      static_assert(
         1 == get_layout_rank<Layout>,
         "Can't copy-construct a layout of different rank.");
   }

   GEOSX_HOST_DEVICE inline
   int index(int idx) const
   {
      GEOSX_UNUSED_VAR(idx);
      GEOSX_ASSERT_KERNEL(
         idx==GEOSX_THREAD_ID(x),
         "The first index (%d) must be equal to the x thread index (%d)"
         " when using SizedDynamic2dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx, GEOSX_THREAD_ID(x));
      return 0;
   }

   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      static_assert(
         N==0,
         "Accessed size is higher than the rank of the Tensor.");
      return size0;
   }
};

template <int BatchSize, int FirstSize, int SecondSize>
class SizedDynamic2dThreadLayout<BatchSize,FirstSize,SecondSize>
{
private:
   const int size0;
   const int size1;
public:
   GEOSX_HOST_DEVICE inline
   SizedDynamic2dThreadLayout(int size0, int size1)
   : size0(size0), size1(size1)
   {
      GEOSX_ASSERT_KERNEL(
         FirstSize==Dynamic || FirstSize==size0,
         "Compilation time (%d) and runtime sizes (%d) must be the same.\n",
         FirstSize, size0);
      GEOSX_ASSERT_KERNEL(
         SecondSize==Dynamic || SecondSize==size1,
         "Compilation time (%d) and runtime sizes (%d) must be the same.\n",
         SecondSize, size1);
      GEOSX_ASSERT_KERNEL(
         size0<=GEOSX_THREAD_SIZE(x),
         "The first dimension (%d) exceeds the number of x threads (%d).\n",
         size0, GEOSX_THREAD_SIZE(x));
      GEOSX_ASSERT_KERNEL(
         size1<=GEOSX_THREAD_SIZE(y),
         "The second dimension (%d) exceeds the number of x threads (%d).\n",
         size1, GEOSX_THREAD_SIZE(y));
      GEOSX_ASSERT_KERNEL(
         BatchSize==GEOSX_THREAD_SIZE(z),
         "The batchsize (%d) is not equal to the number of z threads (%d).\n",
         BatchSize, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamic2dThreadLayout(const Layout &rhs)
   : size0(rhs.template Size<0>()),
     size1(rhs.template Size<1>())
   {
      static_assert(
         2 == get_layout_rank<Layout>,
         "Can't copy-construct a layout of different rank.");
   }

   GEOSX_HOST_DEVICE inline
   int index(int idx0, int idx1) const
   {
      GEOSX_UNUSED_VAR(idx0, idx1);
      GEOSX_ASSERT_KERNEL(
         idx0==GEOSX_THREAD_ID(x),
         "The first index (%d) must be equal to the x thread index (%d)"
         " when using SizedDynamic2dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx0, GEOSX_THREAD_ID(x));
      GEOSX_ASSERT_KERNEL(
         idx1==GEOSX_THREAD_ID(y),
         "The second index (%d) must be equal to the y thread index (%d)"
         " when using SizedDynamic2dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx1, GEOSX_THREAD_ID(y));
      return 0;
   }

   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      static_assert(
         N>=0 && N<2,
         "Accessed size is higher than the rank of the Tensor.");
      return N==0? size0 : size1;
   }
};

// get_layout_rank
template <int BatchSize, int... Sizes>
struct get_layout_rank_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr int value = sizeof...(Sizes);
};

// is_dynamic_layout
template <int BatchSize, int... Sizes>
struct is_dynamic_layout_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr bool value = true;
};

// is_2d_threaded_layout
template <int BatchSize, int... Sizes>
struct is_2d_threaded_layout_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr bool value = true;
};

// is_serial_layout_dim
template <int BatchSize, int... Sizes>
struct is_serial_layout_dim_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>, 0>
{
   static constexpr bool value = false;
};

template <int BatchSize, int... Sizes>
struct is_serial_layout_dim_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>, 1>
{
   static constexpr bool value = false;
};

// is_threaded_layout_dim
template <int BatchSize, int... Sizes>
struct is_threaded_layout_dim_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>, 0>
{
   static constexpr bool value = true;
};

template <int BatchSize, int... Sizes>
struct is_threaded_layout_dim_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>, 1>
{
   static constexpr bool value = true;
};

// get_layout_size
template <int N, int BatchSize, int... Sizes>
struct get_layout_size_v<N, SizedDynamic2dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr int value = get_value<N,Sizes...>;
};

// get_layout_sizes
template <int BatchSize, int... Sizes>
struct get_layout_sizes_t<SizedDynamic2dThreadLayout<BatchSize,Sizes...>>
{
   using type = int_list<Sizes...>;
};

// get_layout_capacity
template <int BatchSize, int First, int Second, int... Rest>
struct get_layout_capacity_v<SizedDynamic2dThreadLayout<BatchSize,First,Second,Rest...>>
{
   static constexpr int value = get_layout_capacity<SizedDynamicLayout<Rest...>>;
};

template <int BatchSize, int First, int Second>
struct get_layout_capacity_v<SizedDynamic2dThreadLayout<BatchSize,First,Second>>
{
   static constexpr int value = 1;
};

template <int BatchSize, int First>
struct get_layout_capacity_v<SizedDynamic2dThreadLayout<BatchSize,First>>
{
   static constexpr int value = 1;
};

// get_layout_batch_size
template <int BatchSize, int... Sizes>
struct get_layout_batch_size_v<SizedDynamic2dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr int value = BatchSize;
};

// get_layout_result_type
template <int BatchSize, int... Sizes>
struct get_layout_result_type_t<SizedDynamic2dThreadLayout<BatchSize,Sizes...>>
{
   template <int... mySizes>
   using type = SizedDynamic2dThreadLayout<BatchSize,mySizes...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_DYNAMIC_2DTHREAD_LAYOUT_HPP_
