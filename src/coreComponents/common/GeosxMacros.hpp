/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GeosxMacros.hpp
 *
 * This file contains various macro definitions.
 */

#include "common/GeosxConfig.hpp"
#include "cxx-utilities/src/Macros.hpp"

#ifndef GEOSX_COMMON_GEOSXMACROS_HPP_
#define GEOSX_COMMON_GEOSXMACROS_HPP_

/// Specify default capture mode for GEOSX lambdas.
#define GEOSX_LAMBDA [=]

/**
 * @name Unused variable markers.
 *
 * These macros are used to explicitly mark a variable/argument as unused
 * and thus silence compiler warnings.
 */
///@{

/// Mark an unused argument and silence compiler warnings.
#define GEOSX_UNUSED_PARAM( X )

/// Mark an unused variable and silence compiler warnings.
#define GEOSX_UNUSED_VAR( X ) ( ( void ) X )

/// Mark a debug variable and silence compiler warnings.
#define GEOSX_DEBUG_VAR( X ) GEOSX_UNUSED_VAR( X )

#if defined(GEOSX_USE_OPENMP)
  #define PRAGMA_OMP( clause ) _Pragma(STRINGIZE(clause))
#else
  #define PRAGMA_OMP( clause )
#endif

/**
 * @name Host-device markers
 *
 * These macros are used to denote host/device/inline functions in a compiler-specific way.
 * They must be prepended to a function or lambda declaration/definition.
 * They will be defined differently when compiled by e.g. a CUDA compiler.
 */
///@{

#if defined(__CUDACC__)
  #define GEOSX_HOST __host__
  #define GEOSX_DEVICE __device__
  #define GEOSX_HOST_DEVICE __host__ __device__
  #define GEOSX_DEVICE_LAMBDA [=] __device__
  #define GEOSX_HOST_DEVICE_LAMBDA [=] __host__ __device__
  #define GEOSX_FORCE_INLINE __forceinline__
#else
/// Marks a host-only function.
  #define GEOSX_HOST
/// Marks a device-only function.
  #define GEOSX_DEVICE
/// Marks a host-device function.
  #define GEOSX_HOST_DEVICE
/// Marks a device-only lambda
  #define GEOSX_DEVICE_LAMBDA [=]
/// Marks a host-device lambda
  #define GEOSX_HOST_DEVICE_LAMBDA [=]
/// Marks a function or lambda for inlining
  #define GEOSX_FORCE_INLINE inline
#endif

///@}

/// preprocessor variable for the C99 restrict keyword for use with pointers
#define GEOSX_RESTRICT LVARRAY_RESTRICT

/// preprocessor variable for the C99 restrict keyword for use with the "this" pointer
#define GEOSX_RESTRICT_THIS LVARRAY_RESTRICT_THIS

#endif // GEOSX_COMMON_GEOSXMACROS_HPP_
