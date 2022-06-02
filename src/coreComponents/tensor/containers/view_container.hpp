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
 * @file read_container.hpp
 */

#ifndef GEOSX_VIEW_CONTAINER_HPP_
#define GEOSX_VIEW_CONTAINER_HPP_

#include "view_container.hpp"

namespace geosx
{

namespace tensor
{

/// A view Container
template <typename Container>
class ViewContainer
{
private:
   using T = get_container_type<Container>;
   Container &data;

public:
   GEOSX_HOST_DEVICE
   ViewContainer(Container &data): data(data) { }

   GEOSX_HOST_DEVICE
   const T& operator[](int i) const
   {
      return data[ i ];
   }

   GEOSX_HOST_DEVICE
   T& operator[](int i)
   {
      return data[ i ];
   }
};

/// A view Container
template <typename Container>
class ConstViewContainer
{
private:
   using T = get_container_type<Container>;
   const Container &data;

public:
   GEOSX_HOST_DEVICE
   ConstViewContainer(const Container &data): data(data) { }

   GEOSX_HOST_DEVICE
   const T& operator[](int i) const
   {
      return data[ i ];
   }
};

// get_container_type
template <typename Container>
struct get_container_type_t<ViewContainer<Container>>
{
   using type = get_container_type<Container>;
};

template <typename Container>
struct get_container_type_t<ConstViewContainer<Container>>
{
   using type = get_container_type<Container>;
};

// get_container_sizes
template <typename Container>
struct get_container_sizes_t<ViewContainer<Container>>
{
   using type = typename get_container_sizes_t<Container>::type;
};

template <typename Container>
struct get_container_sizes_t<ConstViewContainer<Container>>
{
   using type = typename get_container_sizes_t<Container>::type;
};

// get_unsized_container
template <typename Container>
struct get_unsized_container<ViewContainer<Container>>
{
   template <int... Sizes>
   using type = typename get_unsized_container<Container>::template type<Sizes...>;
};

template <typename Container>
struct get_unsized_container<ConstViewContainer<Container>>
{
   template <int... Sizes>
   using type = typename get_unsized_container<Container>::template type<Sizes...>;
};

// is_pointer_container
template <typename Container>
struct is_pointer_container_v<ViewContainer<Container>>
{
   static constexpr bool value = is_pointer_container<Container>;
};

template <typename Container>
struct is_pointer_container_v<ConstViewContainer<Container>>
{
   static constexpr bool value = is_pointer_container<Container>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_VIEW_CONTAINER_HPP_
