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
 * @file stack_container.hpp
 */

#ifndef GEOSX_STACK_CONTAINER_HPP_
#define GEOSX_STACK_CONTAINER_HPP_

#include "container_traits.hpp"
#include "tensor/utilities/util.hpp"
#include "tensor/utilities/int_list.hpp"
#include "tensor/utilities/prod.hpp"

namespace geosx
{

namespace tensor
{

/// Owning Container statically sized and allocated on the stack.
template <typename T, int... Dims>
class StackContainer
{
private:
   T data[prod(Dims...)];

public:
   GEOSX_HOST_DEVICE
   StackContainer() { }

   template <typename... Sizes> GEOSX_HOST_DEVICE
   StackContainer(int size0, Sizes... sizes)
   {
      GEOSX_UNUSED_VAR( size0 );
      GEOSX_UNUSED_VAR( sizes... );
      // static_assert(
      //    sizeof...(Dims)==sizeof...(Sizes)+1,
      //    "Wrong number of dynamic sizes.");
      // TODO verify that Dims == sizes in Debug mode
   }

   GEOSX_HOST_DEVICE
   const T& operator[](int i) const
   {
      GEOSX_ASSERT_KERNEL(
         i<prod(Dims...),
         "Trying to access index %d in StackContainer of capacity %d).",
         i, prod(Dims...));
      return data[ i ];
   }

   GEOSX_HOST_DEVICE
   T& operator[](int i)
   {
      GEOSX_ASSERT_KERNEL(
         i<prod(Dims...),
         "Trying to access index %d in StackContainer of capacity %d).",
         i, prod(Dims...));
      return data[ i ];
   }
};

// get_container_type
template <typename T, int... Dims>
struct get_container_type_t<StackContainer<T,Dims...>>
{
   using type = T;
};

// get_container_sizes
template <typename T, int... Dims>
struct get_container_sizes_t<StackContainer<T, Dims...>>
{
   using type = int_list<Dims...>;
};

// get_unsized_container
template <typename T, int... Dims>
struct get_unsized_container<StackContainer<T, Dims...>>
{
   template <int... Sizes>
   using type = StackContainer<T, Sizes...>;
};

} // namespace tensor

} // namespace geosx

#endif // GEOSX_STACK_CONTAINER_HPP_
