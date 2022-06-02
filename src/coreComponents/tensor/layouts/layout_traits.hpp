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
 * @file layout_traits.hpp
 */

#ifndef GEOSX_LAYOUT_TRAITS_HPP_
#define GEOSX_LAYOUT_TRAITS_HPP_

#include <type_traits>
#include "tensor/utilities/helper_constants.hpp"
#include "tensor/utilities/pow.hpp"

namespace geosx
{

namespace tensor
{

/// Forward declaration of the Tensor class
template <typename Container, typename Layout>
class TensorBase;

/////////////////
// Layout Traits

// get_layout_rank
template <typename Layout>
struct get_layout_rank_v;

template <typename Layout>
constexpr int get_layout_rank = get_layout_rank_v<Layout>::value;

template <typename Container, typename Layout>
struct get_layout_rank_v<TensorBase<Container,Layout>>
{
   static constexpr int value = get_layout_rank<Layout>;
};

// is_dynamic_layout
template <typename Layout>
struct is_dynamic_layout_v
{
   static constexpr bool value = false;
};

template <typename Layout>
constexpr bool is_dynamic_layout = is_dynamic_layout_v<Layout>::value;

// is_static_layout
template <typename Layout>
struct is_static_layout_v
{
   static constexpr bool value = false;
};

template <typename Layout>
constexpr bool is_static_layout = is_static_layout_v<Layout>::value;

// is_serial_layout
template <typename Layout>
struct is_serial_layout_v
{
   static constexpr bool value = false;
};

template <typename Layout>
constexpr bool is_serial_layout = is_serial_layout_v<Layout>::value;

// is_1d_threaded_layout
template <typename Layout>
struct is_1d_threaded_layout_v
{
   static constexpr bool value = false;
};

template <typename Layout>
constexpr bool is_1d_threaded_layout = is_1d_threaded_layout_v<Layout>::value;

// is_2d_threaded_layout
template <typename Layout>
struct is_2d_threaded_layout_v
{
   static constexpr bool value = false;
};

template <typename Layout>
constexpr bool is_2d_threaded_layout = is_2d_threaded_layout_v<Layout>::value;

// is_3d_threaded_layout
template <typename Layout>
struct is_3d_threaded_layout_v
{
   static constexpr bool value = false;
};

template <typename Layout>
constexpr bool is_3d_threaded_layout = is_3d_threaded_layout_v<Layout>::value;

// is_serial_layout_dim
template <typename Layout, int N>
struct is_serial_layout_dim_v
{
   static constexpr bool value = true;
};

template <typename Layout, int N>
constexpr bool is_serial_layout_dim = is_serial_layout_dim_v<Layout,N>::value;

// is_threaded_layout_dim
template <typename Layout, int N>
struct is_threaded_layout_dim_v
{
   static constexpr bool value = false;
};

template <typename Layout, int N>
constexpr bool is_threaded_layout_dim = is_threaded_layout_dim_v<Layout,N>::value;

// get_layout_size
template <int N, typename Layout>
struct get_layout_size_v;

template <int N, typename Layout>
constexpr int get_layout_size = get_layout_size_v<N,Layout>::value;

// get_layout_sizes
template <typename Layout>
struct get_layout_sizes_t;

template <typename Layout>
using get_layout_sizes = typename get_layout_sizes_t<Layout>::type;

// get_layout_batch_size
template <typename Layout>
struct get_layout_batch_size_v
{
   static constexpr int value = 1;
};

template <typename Layout>
constexpr int get_layout_batch_size = get_layout_batch_size_v<Layout>::value;

// get_layout_capacity
template <typename Layout>
struct get_layout_capacity_v;

template <typename Layout>
constexpr int get_layout_capacity = get_layout_capacity_v<Layout>::value;

// get_layout_result_type
template <typename Layout, int... Sizes>
struct get_layout_result_type_t;

template <typename Layout, int... Sizes>
using get_layout_result_type = typename get_layout_result_type_t<Layout>::template type<Sizes...>;

} // namespace tensor

} // namespace geosx

#endif // GEOSX_LAYOUT_TRAITS_HPP_
