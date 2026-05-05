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
 * @file StackTrace.hpp
 */

#ifndef GEOS_COMMON_LOGGER_STACKTRACE_HPP
#define GEOS_COMMON_LOGGER_STACKTRACE_HPP

#include "common/GeosxConfig.hpp" // For the following guards
#ifdef GEOS_USE_CPPTRACE
#include <cpptrace/cpptrace.hpp>
#include <cpptrace/formatting.hpp>
#endif

#include <string>

namespace geos
{

/**
 * @brief Utility class to interact with stack traces.
 */
class StackTrace
{
public:

  /**
   * @brief Get a stack trace for the current thread.
   * @return The stack trace as a string.
   * @note Not signal-safe. Use signalSafeStackTrace() from inside a signal handler.
   */
  static std::string stackTrace();

  /**
   * @brief Get a stack trace from a context where signal-safety is required.
   * @return The stack trace as a string.
   */
  static std::string signalSafeStackTrace();

#ifdef GEOS_USE_CPPTRACE
  /**
   * @brief Access the configured cpptrace formatter.
   * @return The formatter instance.
   */
  static cpptrace::formatter & formatter();

  /**
   * @brief Format a cpptrace stacktrace using the configured formatter.
   * @param stacktrace The cpptrace stack trace to format.
   * @return Formatted stack trace string.
   */
  static std::string formatStackTrace( cpptrace::v1::stacktrace stacktrace );
#endif

};


} /* namespace geos */

#endif /* GEOS_COMMON_LOGGER_STACKTRACE_HPP */
