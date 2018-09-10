/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

#ifndef SRC_COMPONENTS_CORE_SRC_COMMON_TIMING_MACROS_HPP_
#define SRC_COMPONENTS_CORE_SRC_COMMON_TIMING_MACROS_HPP_

#ifdef GEOSX_USE_CALIPER
#include <caliper/cali.h>

#define GEOS_MARK_FUNCTION CALI_CXX_MARK_FUNCTION

#define DO_STRINGIFY(arg) #arg
#define GEOS_CXX_MARK_LOOP_BEGIN(loop, loopName) CALI_CXX_MARK_LOOP_BEGIN(loop,DO_STRINGIFY(loopName))
#define GEOS_CXX_MARK_LOOP_END(loop) CALI_CXX_MARK_LOOP_END(loop)

#define GEOS_MARK_BEGIN(name) CALI_MARK_BEGIN(DO_STRINGIFY(name))
#define GEOS_MARK_END(name) CALI_MARK_END(DO_STRINGIFY(name))

#else

#define GEOS_MARK_FUNCTION
#define GEOS_CXX_MARK_LOOP_BEGIN(loop, loopName)
#define GEOS_CXX_MARK_LOOP_END(loop)
#define GEOS_MARK_BEGIN(name)
#define GEOS_MARK_END(name)
#endif



#endif
