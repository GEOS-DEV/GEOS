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

#define STANDARD_ELEMENT_DNDX_LAYOUT 0
#define STANDARD_ELEMENT_DETJ_LAYOUT 0
#define STANDARD_ELEMENT_MEANSTRESS_LAYOUT 0
#define STANDARD_ELEMENT_DEVIATORSTRESS_LAYOUT 0
#define STANDARD_ELEMENT_TONODESRELATION_LAYOUT 0

#define STORE_NODE_DATA_LOCALLY 1

#define CALC_SHAPE_FUNCTION_DERIVATIVES 1

#define INLINE_STRESS_UPDATE 0

#define SSLE_USE_PATCH_KERNEL 1

#if SSLE_USE_PATCH_KERNEL
  #define SSLE_PATCH_KERNEL_VIZ_OUTPUT 0
  #define SSLE_PATCH_KERNEL_MAX_ELEMS 64
  #define SSLE_PATCH_KERNEL_MAX_NODES 128
  #define SSLE_PATCH_KERNEL_REORDER_NODES 0
#endif

#define GEOSX_LAMBDA [=]

#if defined(__CUDACC__)

#define FORCE_INLINE __forceinline__

#define GEOSX_HOST __host__
#define GEOSX_DEVICE __device__
#define GEOSX_HOST_DEVICE __host__ __device__

#define GEOSX_DEVICE_LAMBDA [=] __device__

#define GEOSX_FORCE_INLINE __forceinline__

#else

#define FORCE_INLINE inline

#define GEOSX_HOST
#define GEOSX_DEVICE
#define GEOSX_HOST_DEVICE

#define GEOSX_DEVICE_LAMBDA [=]

#define GEOSX_FORCE_INLINE
#endif

#endif // COMMON_GEOSXMACROS_HPP_
