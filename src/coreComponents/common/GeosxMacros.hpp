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

#if defined(USE_CUDA) && defined(__CUDACC__)
#define GEOSX_HOST __host__
#define GEOSX_DEVICE __device__
#define GEOSX_HOST_DEVICE __host__ __device__
#else
#define GEOSX_HOST
#define GEOSX_DEVICE
#define GEOSX_HOST_DEVICE
#endif

#endif
