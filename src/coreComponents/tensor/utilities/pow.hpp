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
 * @file pow.hpp
 */

#ifndef GEOSX_TENSOR_UTIL_POW
#define GEOSX_TENSOR_UTIL_POW

namespace geosx
{

namespace tensor
{

/// Compute x to the power n
template <typename T>
constexpr T pow(T x, unsigned int n)
{
   return n == 0 ? 1 : x * pow(x, n-1);
}

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_UTIL_POW
