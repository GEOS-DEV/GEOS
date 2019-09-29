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
 * This file contains various macro definitions
 */

#include "common/GeosxConfig.hpp"

#ifndef COMMON_GEOSXMACROS_HPP_
#define COMMON_GEOSXMACROS_HPP_

#define GEOSX_LAMBDA [=]

// Use this to mark an unused argument and silence compiler warnings
#define GEOSX_UNUSED_ARG( X )

// Use this to mark an unused variable and silence compiler warnings.
#define GEOSX_UNUSED_VAR( X ) ( ( void ) X )

// Use this to mark a debug variable and silence compiler warnings.
#define GEOSX_DEBUG_VAR( X ) GEOSX_UNUSED_VAR( X )

// This will interpret A as a string
#define STRINGIZE_NX( A ) #A

// This will macro expand A and then interpret that as a string.
#define STRINGIZE( A ) STRINGIZE_NX( A )

#if defined(GEOSX_USE_OPENMP)
  #define PRAGMA_OMP( clause ) _Pragma(STRINGIZE(clause))
#else
  #define PRAGMA_OMP( clause )
#endif

#if defined(__CUDACC__)
  #define GEOSX_HOST __host__
  #define GEOSX_DEVICE __device__
  #define GEOSX_HOST_DEVICE __host__ __device__
  #define GEOSX_DEVICE_LAMBDA [=] __device__
  #define GEOSX_HOST_DEVICE_LAMBDA [=] __host__ __device__
  #define GEOSX_FORCE_INLINE __forceinline__
  #define GEOSX_CUDA_PARAM( X ) X
#else
  #define GEOSX_HOST
  #define GEOSX_DEVICE
  #define GEOSX_HOST_DEVICE
  #define GEOSX_DEVICE_LAMBDA [=]
  #define GEOSX_HOST_DEVICE_LAMBDA [=]
  #define GEOSX_FORCE_INLINE inline
  #define GEOSX_CUDA_PARAM( X )
#endif

#endif // COMMON_GEOSXMACROS_HPP_
