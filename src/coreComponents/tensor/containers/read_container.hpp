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

#ifndef GEOSX_READ_CONTAINER_HPP_
#define GEOSX_READ_CONTAINER_HPP_

#include "container_traits.hpp"

namespace geosx
{

namespace tensor
{

/// Non-owning const Container that can be moved between host and device.
template <typename T>
class ReadContainer
{
private:
   const T* data;

public:
   template <typename... Sizes> GEOSX_HOST_DEVICE
   ReadContainer(int size0, Sizes... sizes) : data(nullptr)
   {
      GEOSX_UNUSED_VAR(size0, sizes...);
      // static_assert(false,"Read Container are not supposed to be created like this");
   }

   GEOSX_HOST_DEVICE
   ReadContainer(const T* data) : data(data)
   { }

   GEOSX_HOST_DEVICE
   ReadContainer(const ReadContainer &rhs) : data(rhs.data)
   { }

   GEOSX_HOST_DEVICE
   const T& operator[](int i) const
   {
      return data[ i ];
   }
};

// get_container_type
template <typename T>
struct get_container_type_t<ReadContainer<T>>
{
   using type = T;
};

// is_pointer_container
template <typename T>
struct is_pointer_container_v<ReadContainer<T>>
{
   static constexpr bool value = true;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_READ_CONTAINER_HPP_
