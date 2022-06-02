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
 * @file get_layout.hpp
 */

#ifndef GEOSX_TENSOR_UTIL_LAYOUT
#define GEOSX_TENSOR_UTIL_LAYOUT

#include "tensor/tensor.hpp"

namespace geosx
{

namespace tensor
{

template <typename Container, typename Layout> GEOSX_HOST_DEVICE
auto GetLayout(const TensorBase<Container,Layout>& t)
{
   return static_cast<const Layout&>(t);
}

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_UTIL_LAYOUT
