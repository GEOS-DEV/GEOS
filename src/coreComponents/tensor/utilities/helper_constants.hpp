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
 * @file helper_constants.hpp
 */

#ifndef GEOSX_TENSOR_UTIL_CONSTANTS
#define GEOSX_TENSOR_UTIL_CONSTANTS

namespace geosx
{

namespace tensor
{

/// Compile time constant used to signify dynamic dimensions.
static constexpr int Dynamic = 0;

/// Compile time constant used to signify an error.
static constexpr int Error = -1;

/// The arbitrary maximum dynamic dimension size for stack allocated tensors.
static constexpr int DynamicMaxSize = 16;

#if (defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__))
#define GEOSX_DEVICE_COMPILE
#endif

/// Compile time constant indicating if the code being compiled is for device.
#ifdef GEOSX_DEVICE_COMPILE
static constexpr bool is_device = true;
#else
static constexpr bool is_device = false;
#endif

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_UTIL_CONSTANTS
