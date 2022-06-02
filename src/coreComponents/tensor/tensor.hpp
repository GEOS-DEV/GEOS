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
/// Utility functions to abstract iterating over the tensor dimensions
#include "utilities/foreach.hpp"

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
   /** Construct a Tensor by providing the sizes of its different dimensions,
       both the Container and the Layout need to have similar constructors for
       this Tensor constructor to be usable.
   */
   template <typename... Sizes> GEOSX_HOST_DEVICE
   TensorBase(int size0, Sizes... sizes)
   : Container(size0,sizes...), Layout(size0,sizes...) { }

   /** Construct a Tensor by providing a pointer to its data, and the sizes of
       its Layout. The Container needs to have a constructor using a pointer for
       this constructor to be usable.
   */
   template <typename... Sizes> GEOSX_HOST_DEVICE
   TensorBase(T* ptr, Sizes... sizes)
   : Container(ptr), Layout(sizes...) { }

   /// Utility Constructors
   /** Construct a tensor based on a Layout, the Container needs to be default
       constructible. Note: A Tensor is a Layout (through inheritance).
   */
   GEOSX_HOST_DEVICE
   TensorBase(Layout index): Container(), Layout(index) { }

   /** Construct a Tensor by providing a Container object and a Layout object.
   */
   GEOSX_HOST_DEVICE
   TensorBase(Container data, Layout index): Container(data), Layout(index) { }

   /// Copy Constructors
   /** Copy a Tensor of the same type, the copy is deep or shallow depending on
       the Container.
   */
   GEOSX_HOST_DEVICE
   TensorBase(const TensorBase &rhs): Container(rhs), Layout(rhs) { }

   /** Deep copy of a Tensor of a different type. */
   template <typename OtherTensor,
             std::enable_if_t<
               is_tensor<OtherTensor>,
               bool> = true > GEOSX_HOST_DEVICE
   TensorBase(const OtherTensor &rhs): Container(), Layout(rhs)
   {
      ForallDims<TensorBase>::ApplyBinOp(*this, rhs, [&](auto... idx)
      {
         (*this)(idx...) = rhs(idx...);
      });
   }

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

   /// Initialization of a Tensor to a constant value.
   GEOSX_HOST_DEVICE inline
   TensorBase<Container,Layout>& operator=(const T &val)
   {
      ForallDims<TensorBase>::Apply(*this, [&](auto... idx)
      {
         (*this)(idx...) = val;
      });
      return *this;
   }

   /// operator=, compatible with other types of Tensors
   template <typename OtherTensor,
             std::enable_if_t<
               is_tensor<OtherTensor>,
               bool> = true > GEOSX_HOST_DEVICE inline
   TensorBase<Container,Layout>& operator=(const OtherTensor &rhs)
   {
      ForallDims<TensorBase>::ApplyBinOp(*this, rhs, [&](auto... idx)
      {
         (*this)(idx...) = rhs(idx...);
      });
      return *this;
   }

   /// operator+=, compatible with other types of Tensors
   template <typename OtherTensor,
             std::enable_if_t<
               is_tensor<OtherTensor>,
               bool> = true > GEOSX_HOST_DEVICE inline
   TensorBase<Container,Layout>& operator+=(const OtherTensor &rhs)
   {
      ForallDims<TensorBase>::ApplyBinOp(*this, rhs, [&](auto... idx)
      {
         (*this)(idx...) += rhs(idx...);
      });
      return *this;
   }

   /// operator-=, compatible with other types of Tensors
   template <typename OtherTensor,
             std::enable_if_t<
               is_tensor<OtherTensor>,
               bool> = true > GEOSX_HOST_DEVICE inline
   TensorBase<Container,Layout>& operator-=(const OtherTensor &rhs)
   {
      ForallDims<TensorBase>::ApplyBinOp(*this, rhs, [&](auto... idx)
      {
         (*this)(idx...) -= rhs(idx...);
      });
      return *this;
   }
};

} // namespace tensor

} // namespace tensor

#endif // GEOSX_TENSOR
