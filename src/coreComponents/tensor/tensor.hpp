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
 * @file tensor.hpp
 */

#ifndef GEOSX_TENSOR_HPP_
#define GEOSX_TENSOR_HPP_

/// Contain compilation time functions used with tensors
#include "tensor_traits.hpp"

namespace geosx
{

namespace tensor
{

/** A generic tensor class using a linear memory container storing the values,
    and a layout mapping a rank N index to a linear index corresponding to the
    values indices in the container.
    @a Container is the type of data container, they can either be statically or
       dynamically allocated,
    @a Layout is a class that represents the data layout
       There is two main sub-categories of Layout, Static and Dynamic layouts.
       Dynamic Layout have the following signature:
       template <int Rank>,
       Static Layout have the following signature:
       template <int... Sizes>,
       where Sizes... is the list of the sizes of the dimensions of the Tensor.
   */
template <typename Container,
          typename Layout>
class TensorBase: public Container, public Layout
{
public:
   using T = get_tensor_value_type<TensorBase>;
   using container = Container;
   using layout = Layout;

   /// Default Constructor
   /** In order to use the Tensor default constructor, both the Container and
       the Layout need to have a default constructor.
   */
   GEOSX_HOST_DEVICE
   TensorBase() : Container(), Layout() { }

   /// Main Constructors
   /** Construct a Tensor by providing a pointer to its data, and the sizes of
       its Layout. The Container needs to have a constructor using a pointer for
       this constructor to be usable.
   */
   template <typename... Sizes> GEOSX_HOST_DEVICE
   TensorBase(T* ptr, Sizes... sizes)
   : Container(ptr), Layout(sizes...) { }

   /// Accessor
   /** This operator allows to access a value inside a Tensor by providing
       indices, the number of indices must be equal to the rank of the Tensor.
   */
   template <typename... Idx> GEOSX_HOST_DEVICE inline
   T& operator()(Idx... args)
   {
      static_assert(get_tensor_rank<TensorBase> == sizeof...(Idx),
                    "Wrong number of indices");
      return this->operator[]( this->index(args...) );
   }

   /// Const Accessor
   /** This operator allows to access a const value inside a const Tensor by
       providing indices, the number of indices must be equal to the rank of the
       Tensor.
   */
   template <typename... Idx> GEOSX_HOST_DEVICE inline
   const T& operator()(Idx... args) const
   {
      static_assert(get_tensor_rank<TensorBase> == sizeof...(Idx),
                    "Wrong number of indices");
      return this->operator[]( this->index(args...) );
   }
};

} // namespace tensor

} // namespace tensor

#endif // GEOSX_TENSOR
