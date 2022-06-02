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

/// A Class to compute the real index from the multi-indices of a DynamicLayout
template <int Dim, int N = 1>
struct DynamicLayoutIndex
{
   template <typename... Args> GEOSX_HOST_DEVICE inline
   static int eval(const int* sizes, int first, Args... args)
   {
      GEOSX_ASSERT_KERNEL(
         first<sizes[N-1],
         "Index (%d) greater than the dynamic size (%d).\n", first, sizes[N-1]);
      return first + sizes[N - 1] * DynamicLayoutIndex<Dim,N+1>::eval(sizes, args...);
   }
};

// Terminal case
template <int Dim>
struct DynamicLayoutIndex<Dim, Dim>
{
   GEOSX_HOST_DEVICE inline
   static int eval(const int* sizes, int first)
   {
      GEOSX_UNUSED_VAR(sizes);
      GEOSX_ASSERT_KERNEL(
         first<sizes[Dim-1],
         "Index (%d) greater than the dynamic size (%d).\n", first, sizes[Dim-1]);
      return first;
   }
};

/// A class to initialize the size of a DynamicLayout
template <int Dim, int N = 1>
struct InitDynamicLayout
{
   template <typename... Args> GEOSX_HOST_DEVICE inline
   static void result(int* sizes, int first, Args... args)
   {
      sizes[N - 1] = first;
      InitDynamicLayout<Dim,N+1>::result(sizes, args...);
   }

   template <typename Layout> GEOSX_HOST_DEVICE inline
   static void result(int* sizes, const Layout &rhs)
   {
      sizes[N - 1] = rhs.template Size<N-1>();
      InitDynamicLayout<Dim,N+1>::result(sizes, rhs);
   }
};

// Terminal case
template <int Dim>
struct InitDynamicLayout<Dim, Dim>
{
   template <typename... Args> GEOSX_HOST_DEVICE inline
   static void result(int* sizes, int first, Args... args)
   {
      GEOSX_UNUSED_VAR(args...);
      sizes[Dim - 1] = first;
   }

   template <typename Layout> GEOSX_HOST_DEVICE inline
   static void result(int* sizes, const Layout &rhs)
   {
      sizes[Dim - 1] = rhs.template Size<Dim-1>();
   }
};

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

template<int... Dims>
struct StaticELayoutIndex
{
   template <typename... Idx> GEOSX_HOST_DEVICE inline
   static int eval(Idx... args)
   {
      return StaticIndex<1,sizeof...(Dims)+1,Dims...,Dynamic>::eval(args...);
   }
};

// StaticELayoutSize
template <int StaticSize, int N, int... Sizes>
struct StaticELayoutSize
{
   GEOSX_HOST_DEVICE inline
   static int eval(int last_size)
   {
      GEOSX_UNUSED_VAR(last_size);
      return get_value<N,Sizes...>;
   }
};

template <int StaticSize, int... Sizes>
struct StaticELayoutSize<StaticSize, StaticSize, Sizes...>
{
   GEOSX_HOST_DEVICE inline
   static int eval(int last_size)
   {
      return last_size;
   }
};

template <int N, int Rank>
struct DynamicBlockLayoutSize
{
   template <typename Layout> GEOSX_HOST_DEVICE inline
   static int eval(int size0, int size1, const Layout &layout)
   {
      GEOSX_UNUSED_VAR(size0, size1);
      return layout.template Size<N-2>();
   }
};

template <int Rank>
struct DynamicBlockLayoutSize<0, Rank>
{
   template <typename Layout> GEOSX_HOST_DEVICE inline
   static int eval(int size0, int size1, const Layout &layout)
   {
      GEOSX_UNUSED_VAR(size1, layout);
      return size0;
   }
};

template <int Rank>
struct DynamicBlockLayoutSize<1, Rank>
{
   template <typename Layout> GEOSX_HOST_DEVICE inline
   static int eval(int size0, int size1, const Layout &layout)
   {
      GEOSX_UNUSED_VAR(size0, layout);
      return size1;
   }
};

template <int N, int Rank>
struct Dynamic3dThreadLayoutSize
{
   template <typename Layout> GEOSX_HOST_DEVICE inline
   static int eval(int size0, int size1, int size2,
                   const Layout &layout)
   {
      GEOSX_UNUSED_VAR(size0, size1, size2);
      return layout.template Size<N-3>();
   }
};

template <int Rank>
struct Dynamic3dThreadLayoutSize<0, Rank>
{
   template <typename Layout> GEOSX_HOST_DEVICE inline
   static int eval(int size0, int size1, int size2,
                   const Layout &layout)
   {
      GEOSX_UNUSED_VAR(size1, size2, layout);
      return size0;
   }
};

template <int Rank>
struct Dynamic3dThreadLayoutSize<1, Rank>
{
   template <typename Layout> GEOSX_HOST_DEVICE inline
   static int eval(int size0, int size1, int size2,
                   const Layout &layout)
   {
      GEOSX_UNUSED_VAR(size0, size2, layout);
      return size1;
   }
};

template <int Rank>
struct Dynamic3dThreadLayoutSize<2, Rank>
{
   template <typename Layout> GEOSX_HOST_DEVICE inline
   static int eval(int size0, int size1, int size2,
                   const Layout &layout)
   {
      GEOSX_UNUSED_VAR(size0, size1, layout);
      return size2;
   }
};

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_LAYOUT_IMPL_HPP_
