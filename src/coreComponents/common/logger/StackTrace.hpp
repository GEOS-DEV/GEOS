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

#include <string>
#include <vector>

namespace geos
{

struct StackTraceParams;

/**
 * @brief Utility class to interact with stack traces.
 */
class StackTrace
{
public:

  /**
   * @brief StackTrace constructor
   * @param params StackTraceParams class containing the type of stacktrace used,
   *               based on the build.
   */
  StackTrace( StackTraceParams const & params );

  /**
   * @brief Get a stack trace for the current thread.
   * @return The stack trace object containing the stack trace frames.
   */
  static StackTrace stackTrace();

  /**
   * @brief Stack trace frames accessor.
   * @return Container of frames from this stack trace.
   */
  std::vector< std::string > const & frames() const
  { return m_frames; }

  bool isValidStackTrace() const
  { return m_isValidStackTrace; }

private:

  std::vector< std::string > m_frames;

  bool m_isValidStackTrace = false;

};


} /* namespace geos */

#endif /* GEOS_COMMON_LOGGER_STACKTRACE_HPP */
