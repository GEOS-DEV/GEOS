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
 * @file get.hpp
 */

#ifndef GEOSX_TENSOR_GET
#define GEOSX_TENSOR_GET

#include "tensor/containers/container_traits.hpp"
#include "tensor/containers/view_container.hpp"
#include "tensor/layouts/layout_traits.hpp"
#include "tensor/layouts/restricted_layout.hpp"

namespace geosx
{

namespace tensor
{

/// Lazy accessor for the sub-Tensor extracted from idx in Nth dimension.
/** ex: `auto slice = Get<Dim>(idx, tensor);'
*/
template <int Dim, typename Container, typename Layout> GEOSX_HOST_DEVICE inline
auto Get(int idx, TensorBase<Container,Layout> &t)
{
   static_assert(Dim>=0 && Dim<get_layout_rank<Layout>,
      "Cannot access this dimension with Get");
   using C = ViewContainer<Container>;
   using L = RestrictedLayout<Dim,Layout>;
   using RestrictedTensor = TensorBase<C,L>;
   C data(t);
   L layout(idx,t);
   return RestrictedTensor(data,layout);
}

/// Lazy const accessor for the sub-Tensor extracted from idx in Nth dimension.
template <int Dim, typename Container, typename Layout> GEOSX_HOST_DEVICE inline
auto Get(int idx, const TensorBase<Container,Layout> &t)
{
   static_assert(Dim>=0 && Dim<get_layout_rank<Layout>,
      "Cannot access this dimension with Get");
   using C = ConstViewContainer<Container>;
   using L = RestrictedLayout<Dim,Layout>;
   using RestrictedTensor = TensorBase<C,L>;
   C data(t);
   L layout(idx,t);
   return RestrictedTensor(data,layout);
}

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_GET
