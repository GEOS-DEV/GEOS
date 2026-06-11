/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file StackTraceData.hpp
 */

#ifndef GEOS_COMMON_LOGGER_STACKTRACEPARAMS_HPP
#define GEOS_COMMON_LOGGER_STACKTRACEPARAMS_HPP

#include "common/GeosxConfig.hpp"
#ifdef GEOS_USE_CPPTRACE
#include <cpptrace/cpptrace.hpp>
#include <cpptrace/from_current.hpp>
#endif

#include <string>

namespace geos
{

#ifdef GEOS_USE_CPPTRACE
using StackTraceType = cpptrace::stacktrace;
#else // fallback to LvArray stacktrace
using StackTraceType = std::string
#endif


/**
 * @struct StackTraceParams
 * @brief Dispatch struct to store the stacktrace type used by the current build.
 */
struct StackTraceParams
{
  StackTraceType trace;
};


#ifdef GEOS_USE_CPPTRACE
#define GEOS_STACKTRACE_PARAMETERS StackTraceParams{ cpptrace::from_current_exception() }
#else
#define GEOS_STACKTRACE_PARAMETERS StackTraceParams{ LvArray::system::stackTrace( true ) }
#endif


} /* namespace geos */

#endif /* GEOS_COMMON_LOGGER_STACKTRACEPARAMS_HPP */
