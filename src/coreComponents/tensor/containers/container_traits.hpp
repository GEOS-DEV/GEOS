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
 * @file container_traits.hpp
 */

#ifndef GEOSX_CONTAINER_TRAITS_HPP_
#define GEOSX_CONTAINER_TRAITS_HPP_

namespace geosx
{

namespace tensor
{

////////////////////
// Container Traits

// get_container_type
template <typename Container>
struct get_container_type_t;

template <typename Container>
using get_container_type = typename get_container_type_t<Container>::type;

// get_container_sizes
template <typename Container>
struct get_container_sizes_t;

template <typename Container>
using get_container_sizes = typename get_container_sizes_t<Container>::type;

// get_unsized_container
template <typename Container>
struct get_unsized_container;

// is_pointer_container
template <typename Container>
struct is_pointer_container_v
{
   static constexpr bool value = false;
};

template <typename Tensor>
constexpr bool is_pointer_container = is_pointer_container_v<Tensor>::value;

} // namespace tensor

} // namespace geosx

#endif // GEOSX_CONTAINER_TRAITS_HPP_
