/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file GeosxMacros.hpp
 * This file contains various macro definitions
 */

#include "common/GeosxConfig.hpp"

#ifndef COMMON_GEOSXMACROS_HPP_
#define COMMON_GEOSXMACROS_HPP_

#define GEOSX_LAMBDA [=]

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
#else
  #define GEOSX_HOST
  #define GEOSX_DEVICE
  #define GEOSX_HOST_DEVICE
  #define GEOSX_DEVICE_LAMBDA [=]
  #define GEOSX_HOST_DEVICE_LAMBDA [=]
  #define GEOSX_FORCE_INLINE inline
#endif

#endif // COMMON_GEOSXMACROS_HPP_
