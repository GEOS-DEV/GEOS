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
 * @file dynamic_3dthread_layout.hpp
 */

#ifndef GEOSX_DYNAMIC_3DTHREAD_LAYOUT_HPP_
#define GEOSX_DYNAMIC_3DTHREAD_LAYOUT_HPP_

#include "tensor/utilities/error.hpp"
#include "dynamic_layout.hpp"
#include "layout_traits.hpp"

namespace geosx
{

namespace tensor
{

template <int BatchSize, int... Sizes>
class SizedDynamic3dThreadLayout;

template <int Rank, int BatchSize>
using Dynamic3dThreadLayout = instantiate<SizedDynamic3dThreadLayout,
                                          append<
                                             int_list<BatchSize>,
                                             int_repeat<Dynamic,Rank> > >;

template <int BatchSize, int FirstSize, int SecondSize, int ThirdSize, int... Sizes>
class SizedDynamic3dThreadLayout<BatchSize,FirstSize,SecondSize,ThirdSize,Sizes...>
{
private:
   const int size0;
   const int size1;
   const int size2;
   SizedDynamicLayout<Sizes...> layout;
public:
   template <typename... Ss> GEOSX_HOST_DEVICE inline
   SizedDynamic3dThreadLayout(int size0, int size1, int size2, Ss... sizes)
   : size0(size0), size1(size1), size2(size2), layout(sizes...)
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
         ThirdSize==Dynamic || ThirdSize==size2,
         "Compilation time (%d) and runtime sizes (%d) must be the same.\n",
         ThirdSize, size2);
      GEOSX_ASSERT_KERNEL(
         size0<=GEOSX_THREAD_SIZE(x),
         "The first dimension (%d) exceeds the number of x threads (%d).\n",
         size0, GEOSX_THREAD_SIZE(x));
      GEOSX_ASSERT_KERNEL(
         size1<=GEOSX_THREAD_SIZE(y),
         "The second dimension (%d) exceeds the number of y threads (%d).\n",
         size1, GEOSX_THREAD_SIZE(y));
      GEOSX_ASSERT_KERNEL(
         size2<=GEOSX_THREAD_SIZE(z),
         "The third dimension (%d) exceeds the number of z threads (%d).\n",
         size2, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamic3dThreadLayout(const Layout &rhs)
   : size0(rhs.template Size<0>()),
     size1(rhs.template Size<1>()),
     size2(rhs.template Size<2>()),
     layout( Get<0>( 0, Get<0>( 0, Get<0>( 0, rhs ) ) ) )
   {
      constexpr int Rank = sizeof...(Sizes) + 3;
      static_assert(
         Rank == get_layout_rank<Layout>,
         "Can't copy-construct a layout of different rank.");
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   int index(int idx0, int idx1, int idx2, Idx... idx) const
   {
      GEOSX_UNUSED_VAR(idx0, idx1, idx2);
      GEOSX_ASSERT_KERNEL(
         idx0==GEOSX_THREAD_ID(x),
         "The first index (%d) must be equal to the x thread index (%d)"
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx0, GEOSX_THREAD_ID(x));
      GEOSX_ASSERT_KERNEL(
         idx1==GEOSX_THREAD_ID(y),
         "The second index (%d) must be equal to the y thread index (%d)"
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx1, GEOSX_THREAD_ID(y));
      GEOSX_ASSERT_KERNEL(
         idx2==GEOSX_THREAD_ID(z),
         "The third index (%d) must be equal to the z thread index (%d)"
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx2, GEOSX_THREAD_ID(z));
      return layout.index(idx...);
   }

   // Can be constexpr if Tensor inherit from Layout
   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      constexpr int Rank = sizeof...(Sizes) + 3;
      static_assert(
         N>=0 && N<Rank,
         "Accessed size is higher than the rank of the Tensor.");
      return Dynamic3dThreadLayoutSize<N,Rank>::eval(size0,size1,size2,layout);
   }
};

template <int BatchSize, int FirstSize>
class SizedDynamic3dThreadLayout<BatchSize,FirstSize>
{
private:
   const int size0;
public:
   GEOSX_HOST_DEVICE inline
   SizedDynamic3dThreadLayout(int size0)
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
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamic3dThreadLayout(const Layout &rhs)
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
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx, GEOSX_THREAD_ID(x));
      return 0;
   }

   // Can be constexpr if Tensor inherit from Layout
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
class SizedDynamic3dThreadLayout<BatchSize,FirstSize,SecondSize>
{
private:
   const int size0;
   const int size1;
public:
   GEOSX_HOST_DEVICE inline
   SizedDynamic3dThreadLayout(int size0, int size1)
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
         "The second dimension (%d) exceeds the number of y threads (%d).\n",
         size1, GEOSX_THREAD_SIZE(y));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamic3dThreadLayout(const Layout &rhs)
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
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx0, GEOSX_THREAD_ID(x));
      GEOSX_ASSERT_KERNEL(
         idx1==GEOSX_THREAD_ID(y),
         "The second index (%d) must be equal to the y thread index (%d)"
         " when using SizedDynamic3dThreadLayout. Use shared memory"
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

template <int BatchSize, int FirstSize, int SecondSize, int ThirdSize>
class SizedDynamic3dThreadLayout<BatchSize,FirstSize,SecondSize,ThirdSize>
{
private:
   const int size0;
   const int size1;
   const int size2;
public:
   GEOSX_HOST_DEVICE inline
   SizedDynamic3dThreadLayout(int size0, int size1, int size2)
   : size0(size0), size1(size1), size2(size2)
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
         ThirdSize==Dynamic || ThirdSize==size2,
         "Compilation time (%d) and runtime sizes (%d) must be the same.\n",
         ThirdSize, size2);
      GEOSX_ASSERT_KERNEL(
         size0<=GEOSX_THREAD_SIZE(x),
         "The first dimension (%d) exceeds the number of x threads (%d).\n",
         size0, GEOSX_THREAD_SIZE(x));
      GEOSX_ASSERT_KERNEL(
         size1<=GEOSX_THREAD_SIZE(y),
         "The second dimension (%d) exceeds the number of y threads (%d).\n",
         size1, GEOSX_THREAD_SIZE(y));
      GEOSX_ASSERT_KERNEL(
         size2<=GEOSX_THREAD_SIZE(z),
         "The third dimension (%d) exceeds the number of z threads (%d).\n",
         size2, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   SizedDynamic3dThreadLayout(const Layout &rhs)
   : size0(rhs.template Size<0>()),
     size1(rhs.template Size<1>()),
     size2(rhs.template Size<2>())
   {
      static_assert(
         3 == get_layout_rank<Layout>,
         "Can't copy-construct a layout of different rank.");
   }

   GEOSX_HOST_DEVICE inline
   int index(int idx0, int idx1, int idx2) const
   {
      GEOSX_UNUSED_VAR(idx0, idx1, idx2);
      GEOSX_ASSERT_KERNEL(
         idx0==GEOSX_THREAD_ID(x),
         "The first index (%d) must be equal to the x thread index (%d)"
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx0, GEOSX_THREAD_ID(x));
      GEOSX_ASSERT_KERNEL(
         idx1==GEOSX_THREAD_ID(y),
         "The second index (%d) must be equal to the y thread index (%d)"
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx1, GEOSX_THREAD_ID(y));
      GEOSX_ASSERT_KERNEL(
         idx2==GEOSX_THREAD_ID(z),
         "The third index (%d) must be equal to the z thread index (%d)"
         " when using SizedDynamic3dThreadLayout. Use shared memory"
         " to access values stored in a different thread.\n",
         idx2, GEOSX_THREAD_ID(z));
      return 0;
   }

   // Can be constexpr if Tensor inherit from Layout
   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      static_assert(
         N>=0 && N<3,
         "Accessed size is higher than the rank of the Tensor.");
      return N==0? size0 : (N==1? size1 : size2);
   }
};

// get_layout_rank
template <int BatchSize, int... Sizes>
struct get_layout_rank_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr int value = sizeof...(Sizes);
};

// is_dynamic_layout
template <int BatchSize, int... Sizes>
struct is_dynamic_layout_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr bool value = true;
};

// is_3d_threaded_layout
template <int BatchSize, int... Sizes>
struct is_3d_threaded_layout_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr bool value = true;
};

// is_serial_layout_dim
template <int BatchSize, int... Sizes>
struct is_serial_layout_dim_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>, 0>
{
   static constexpr bool value = false;
};

template <int BatchSize, int... Sizes>
struct is_serial_layout_dim_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>, 1>
{
   static constexpr bool value = false;
};

template <int BatchSize, int... Sizes>
struct is_serial_layout_dim_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>, 2>
{
   static constexpr bool value = false;
};

// is_threaded_layout_dim
template <int BatchSize, int... Sizes>
struct is_threaded_layout_dim_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>, 0>
{
   static constexpr bool value = true;
};

template <int BatchSize, int... Sizes>
struct is_threaded_layout_dim_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>, 1>
{
   static constexpr bool value = true;
};

template <int BatchSize, int... Sizes>
struct is_threaded_layout_dim_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>, 2>
{
   static constexpr bool value = true;
};

// get_layout_size
template <int N, int BatchSize, int... Sizes>
struct get_layout_size_v<N, SizedDynamic3dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr int value = get_value<N,Sizes...>;
};

// get_layout_sizes
template <int BatchSize, int... Sizes>
struct get_layout_sizes_t<SizedDynamic3dThreadLayout<BatchSize,Sizes...>>
{
   using type = int_list<Sizes...>;
};

// get_layout_capacity
template <int BatchSize, int First, int Second, int Third, int... Rest>
struct get_layout_capacity_v<SizedDynamic3dThreadLayout<BatchSize,First,Second,Third,Rest...>>
{
   static constexpr int value = get_layout_capacity<SizedDynamicLayout<Rest...>>;
};

template <int BatchSize, int First, int Second, int Third>
struct get_layout_capacity_v<SizedDynamic3dThreadLayout<BatchSize,First,Second,Third>>
{
   static constexpr int value = 1;
};

template <int BatchSize, int First, int Second>
struct get_layout_capacity_v<SizedDynamic3dThreadLayout<BatchSize,First,Second>>
{
   static constexpr int value = 1;
};

template <int BatchSize, int First>
struct get_layout_capacity_v<SizedDynamic3dThreadLayout<BatchSize,First>>
{
   static constexpr int value = 1;
};

// get_layout_batch_size
template <int BatchSize, int... Sizes>
struct get_layout_batch_size_v<SizedDynamic3dThreadLayout<BatchSize,Sizes...>>
{
   static constexpr int value = BatchSize;
};

// get_layout_result_type
template <int BatchSize, int... Sizes>
struct get_layout_result_type_t<SizedDynamic3dThreadLayout<BatchSize,Sizes...>>
{
   template <int... mySizes>
   using type = SizedDynamic3dThreadLayout<BatchSize,mySizes...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_DYNAMIC_3DTHREAD_LAYOUT_HPP_
