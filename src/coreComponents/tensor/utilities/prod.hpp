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
 * @file prod.hpp
 */

#ifndef GEOSX_TENSOR_UTIL_PROD
#define GEOSX_TENSOR_UTIL_PROD

namespace geosx
{

namespace tensor
{

/// Compute the product of a list of values
template <typename T> GEOSX_HOST_DEVICE
constexpr T prod(T first) {
   return first;
}

template <typename T, typename... D> GEOSX_HOST_DEVICE
constexpr T prod(T first, D... rest) {
   return first*prod(rest...);
}

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_UTIL_PROD
