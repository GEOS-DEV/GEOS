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
 * @file static_2dthread_layout.hpp
 */

#ifndef GEOSX_STATIC_2DTHREAD_LAYOUT_HPP_
#define GEOSX_STATIC_2DTHREAD_LAYOUT_HPP_

#include "tensor/utilities/error.hpp"
#include "static_layout.hpp"
#include "layout_traits.hpp"

namespace geosx
{

namespace tensor
{

/// Layout using a thread plane to distribute data
template <int... Dims>
class Static2dThreadLayout;

template <int DimX>
class Static2dThreadLayout<DimX>
{
public:
   GEOSX_HOST_DEVICE
   Static2dThreadLayout()
   {
    //   GEOSX_ASSERT_KERNEL(
    //      DimX<=GEOSX_THREAD_SIZE(x),
    //      "The first dimension (%d) exceeds the number of x threads (%d).\n",
    //      DimX, GEOSX_THREAD_SIZE(x));
    //   GEOSX_ASSERT_KERNEL(
    //      BatchSize==GEOSX_THREAD_SIZE(z),
    //      "The batchsize (%d) is not equal to the number of z threads (%d).\n",
    //      BatchSize, GEOSX_THREAD_SIZE(z));
   }

   GEOSX_HOST_DEVICE inline
   Static2dThreadLayout(int size0)
   {
      GEOSX_UNUSED_VAR(size0);
    //   GEOSX_ASSERT_KERNEL(
    //      size0==DimX,
    //      "The runtime first dimension (%d) is different to the compilation one (%d).\n",
    //      size0, DimX);
    //   GEOSX_ASSERT_KERNEL(
    //      DimX<=GEOSX_THREAD_SIZE(x),
    //      "The first dimension (%d) exceeds the number of x threads (%d).\n",
    //      DimX, GEOSX_THREAD_SIZE(x));
    //   GEOSX_ASSERT_KERNEL(
    //      BatchSize==GEOSX_THREAD_SIZE(z),
    //      "The batchsize (%d) is not equal to the number of z threads (%d).\n",
    //      BatchSize, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   Static2dThreadLayout(const Layout& rhs)
   {
      GEOSX_UNUSED_VAR(rhs);
      static_assert(
         1 == get_layout_rank<Layout>,
         "Can't copy-construct a layout of different rank.");
    //   GEOSX_ASSERT_KERNEL(
    //      rhs.template Size<0>() == DimX,
    //      "Layouts sizes don't match %d != %d.\n",
    //      DimX, rhs.template Size<0>());
   }

   GEOSX_HOST_DEVICE inline
   int index(int idx0) const
   {
      GEOSX_UNUSED_VAR(idx0);
    //   GEOSX_ASSERT_KERNEL(
    //      idx0==GEOSX_THREAD_ID(x),
    //      "The first index (%d) must be equal to the x thread index (%d)"
    //      " when using SizedStatic2dThreadLayout. Use shared memory"
    //      " to access values stored in a different thread.\n",
    //      idx0, GEOSX_THREAD_ID(x));
      return 0;
   }

   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      static_assert(
         N==0,
         "Accessed size is higher than the rank of the Tensor.");
      return DimX;
   }
};

template <int DimX, int DimY>
class Static2dThreadLayout<DimX, DimY>
{
public:
   GEOSX_HOST_DEVICE
   Static2dThreadLayout()
   {
    //   GEOSX_ASSERT_KERNEL(
    //      DimX<=GEOSX_THREAD_SIZE(x),
    //      "The first dimension (%d) exceeds the number of x threads (%d).\n",
    //      DimX, GEOSX_THREAD_SIZE(x));
    //   GEOSX_ASSERT_KERNEL(
    //      DimY<=GEOSX_THREAD_SIZE(y),
    //      "The second dimension (%d) exceeds the number of y threads (%d).\n",
    //      DimY, GEOSX_THREAD_SIZE(y));
    //   GEOSX_ASSERT_KERNEL(
    //      BatchSize==GEOSX_THREAD_SIZE(z),
    //      "The batchsize (%d) is not equal to the number of z threads (%d).\n",
    //      BatchSize, GEOSX_THREAD_SIZE(z));
   }

   GEOSX_HOST_DEVICE inline
   Static2dThreadLayout(int size0, int size1)
   {
      GEOSX_UNUSED_VAR(size0, size1);
    //   GEOSX_ASSERT_KERNEL(
    //      size0==DimX,
    //      "The runtime first dimension (%d) is different to the compilation one (%d).\n",
    //      size0, DimX);
    //   GEOSX_ASSERT_KERNEL(
    //      DimX<=GEOSX_THREAD_SIZE(x),
    //      "The first dimension (%d) exceeds the number of x threads (%d).\n",
    //      DimX, GEOSX_THREAD_SIZE(x));
    //   GEOSX_ASSERT_KERNEL(
    //      size1==DimY,
    //      "The runtime second dimension (%d) is different to the compilation one (%d).\n",
    //      size1, DimY);
    //   GEOSX_ASSERT_KERNEL(
    //      DimY<=GEOSX_THREAD_SIZE(y),
    //      "The second dimension (%d) exceeds the number of y threads (%d).\n",
    //      DimY, GEOSX_THREAD_SIZE(y));
    //   GEOSX_ASSERT_KERNEL(
    //      BatchSize==GEOSX_THREAD_SIZE(z),
    //      "The batchsize (%d) is not equal to the number of z threads (%d).\n",
    //      BatchSize, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   Static2dThreadLayout(const Layout& rhs)
   {
      GEOSX_UNUSED_VAR(rhs);
      static_assert(
         2 == get_layout_rank<Layout>,
         "Can't copy-construct a layout of different rank.");
    //   GEOSX_ASSERT_KERNEL(
    //      rhs.template Size<0>() == DimX,
    //      "Layouts sizes don't match %d != %d.\n",
    //      DimX, rhs.template Size<0>());
    //   GEOSX_ASSERT_KERNEL(
    //      rhs.template Size<1>() == DimY,
    //      "Layouts sizes don't match %d != %d.\n",
    //      DimY, rhs.template Size<1>());
   }

   GEOSX_HOST_DEVICE inline
   int index(int idx0, int idx1) const
   {
      GEOSX_UNUSED_VAR(idx0, idx1);
    //   GEOSX_ASSERT_KERNEL(
    //      idx0==GEOSX_THREAD_ID(x),
    //      "The first index (%d) must be equal to the x thread index (%d)"
    //      " when using SizedStatic2dThreadLayout. Use shared memory"
    //      " to access values stored in a different thread.\n",
    //      idx0, GEOSX_THREAD_ID(x));
    //   GEOSX_ASSERT_KERNEL(
    //      idx1==GEOSX_THREAD_ID(y),
    //      "The second index (%d) must be equal to the y thread index (%d)"
    //      " when using SizedStatic2dThreadLayout. Use shared memory"
    //      " to access values stored in a different thread.\n",
    //      idx1, GEOSX_THREAD_ID(y));
      return 0;
   }

   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      static_assert(
         N>=0 && N<2,
         "Accessed size is higher than the rank of the Tensor.");
      return get_value<N,DimX,DimY>;
   }
};

template <int DimX, int DimY, int... Dims>
class Static2dThreadLayout<DimX, DimY, Dims...>
{
private:
   StaticLayout<Dims...> layout;
public:
   GEOSX_HOST_DEVICE
   Static2dThreadLayout()
   {
    //   GEOSX_ASSERT_KERNEL(
    //      DimX<=GEOSX_THREAD_SIZE(x),
    //      "The first dimension (%d) exceeds the number of x threads (%d).\n",
    //      DimX, GEOSX_THREAD_SIZE(x));
    //   GEOSX_ASSERT_KERNEL(
    //      DimY<=GEOSX_THREAD_SIZE(y),
    //      "The second dimension (%d) exceeds the number of y threads (%d).\n",
    //      DimY, GEOSX_THREAD_SIZE(y));
    //   GEOSX_ASSERT_KERNEL(
    //      BatchSize==GEOSX_THREAD_SIZE(z),
    //      "The batchsize (%d) is not equal to the number of z threads (%d).\n",
    //      BatchSize, GEOSX_THREAD_SIZE(z));
   }

   template <typename... Sizes> GEOSX_HOST_DEVICE inline
   Static2dThreadLayout(int size0, int size1, Sizes... sizes)
   : layout(sizes...)
   {
      GEOSX_UNUSED_VAR(size0, size1);
    //   GEOSX_ASSERT_KERNEL(
    //      size0==DimX,
    //      "The runtime first dimension (%d) is different to the compilation one (%d).\n",
    //      size0, DimX);
    //   GEOSX_ASSERT_KERNEL(
    //      DimX<=GEOSX_THREAD_SIZE(x),
    //      "The first dimension (%d) exceeds the number of x threads (%d).\n",
    //      DimX, GEOSX_THREAD_SIZE(x));
    //   GEOSX_ASSERT_KERNEL(
    //      size1==DimY,
    //      "The runtime second dimension (%d) is different to the compilation one (%d).\n",
    //      size1, DimY);
    //   GEOSX_ASSERT_KERNEL(
    //      DimY<=GEOSX_THREAD_SIZE(y),
    //      "The second dimension (%d) exceeds the number of y threads (%d).\n",
    //      DimY, GEOSX_THREAD_SIZE(y));
    //   GEOSX_ASSERT_KERNEL(
    //      BatchSize==GEOSX_THREAD_SIZE(z),
    //      "The batchsize (%d) is not equal to the number of z threads (%d).\n",
    //      BatchSize, GEOSX_THREAD_SIZE(z));
   }

   template <typename Layout> GEOSX_HOST_DEVICE
   Static2dThreadLayout(const Layout& rhs)
   {
    //   GEOSX_UNUSED_VAR(rhs);
    //   static_assert(
    //      2 + sizeof...(Dims) == get_layout_rank<Layout>,
    //      "Can't copy-construct a layout of different rank.");
    //   GEOSX_ASSERT_KERNEL(
    //      rhs.template Size<0>() == DimX,
    //      "Layouts sizes don't match %d != %d.\n",
    //      DimX, rhs.template Size<0>());
    //   GEOSX_ASSERT_KERNEL(
    //      rhs.template Size<1>() == DimY,
    //      "Layouts sizes don't match %d != %d.\n",
    //      DimY, rhs.template Size<1>());
   }

   template <typename... Idx> GEOSX_HOST_DEVICE inline
   int index(int idx0, int idx1, Idx... idx) const
   {
      GEOSX_UNUSED_VAR(idx0, idx1);
    //   GEOSX_ASSERT_KERNEL(
    //      idx0==GEOSX_THREAD_ID(x),
    //      "The first index (%d) must be equal to the x thread index (%d)"
    //      " when using SizedStatic2dThreadLayout. Use shared memory"
    //      " to access values stored in a different thread.\n",
    //      idx0, GEOSX_THREAD_ID(x));
    //   GEOSX_ASSERT_KERNEL(
    //      idx1==GEOSX_THREAD_ID(y),
    //      "The second index (%d) must be equal to the y thread index (%d)"
    //      " when using SizedStatic2dThreadLayout. Use shared memory"
    //      " to access values stored in a different thread.\n",
    //      idx1, GEOSX_THREAD_ID(y));
      return layout.index(idx...);
   }

   template <int N> GEOSX_HOST_DEVICE inline
   constexpr int Size() const
   {
      static_assert(
         N>=0 && N<rank<DimX,DimY,Dims...>,
         "Accessed size is higher than the rank of the Tensor.");
      return get_value<N,DimX,DimY,Dims...>;
   }
};

// get_layout_rank
template <int... Dims>
struct get_layout_rank_v<Static2dThreadLayout<Dims...>>
{
   static constexpr int value = sizeof...(Dims);
};

// is_static_layout
template <int... Dims>
struct is_static_layout_v<Static2dThreadLayout<Dims...>>
{
   static constexpr bool value = true;
};

// is_2d_threaded_layout
template <int... Dims>
struct is_2d_threaded_layout_v<Static2dThreadLayout<Dims...>>
{
   static constexpr bool value = true;
};

// is_serial_layout_dim
template <int... Dims>
struct is_serial_layout_dim_v<Static2dThreadLayout<Dims...>, 0>
{
   static constexpr bool value = false;
};

template <int... Dims>
struct is_serial_layout_dim_v<Static2dThreadLayout<Dims...>, 1>
{
   static constexpr bool value = false;
};

// is_threaded_layout_dim
template <int... Dims>
struct is_threaded_layout_dim_v<Static2dThreadLayout<Dims...>, 0>
{
   static constexpr bool value = true;
};

template <int... Dims>
struct is_threaded_layout_dim_v<Static2dThreadLayout<Dims...>, 1>
{
   static constexpr bool value = true;
};

// get_layout_size
template <int N, int... Dims>
struct get_layout_size_v<N, Static2dThreadLayout<Dims...>>
{
   static constexpr int value = get_value<N, Dims...>;
};

// get_layout_sizes
template <int... Dims>
struct get_layout_sizes_t<Static2dThreadLayout<Dims...>>
{
   using type = int_list<Dims...>;
};

// get_layout_batch_size
template <int... Dims>
struct get_layout_batch_size_v<Static2dThreadLayout<Dims...>>
{
   static constexpr int value = 1;
};

// get_layout_capacity
template <int DimX>
struct get_layout_capacity_v<Static2dThreadLayout<DimX>>
{
   static constexpr int value = 1;
};

template <int DimX, int DimY>
struct get_layout_capacity_v<Static2dThreadLayout<DimX, DimY>>
{
   static constexpr int value = 1;
};

template <int DimX, int DimY, int... Dims>
struct get_layout_capacity_v<Static2dThreadLayout<DimX, DimY, Dims...>>
{
   static constexpr int value = prod(Dims...);
};

// get_layout_result_type
template <int... Sizes>
struct get_layout_result_type_t<Static2dThreadLayout<Sizes...>>
{
   template <int... Dims>
   using type = Static2dThreadLayout<Dims...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_STATIC_2DTHREAD_LAYOUT_HPP_
