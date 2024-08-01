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
 * @file LogRouter.hpp
 */

#ifndef GEOS_COMMON_LOGROUTER_HPP
#define GEOS_COMMON_LOGROUTER_HPP

#include "common/LogOutput.hpp"
#include "common/LogErrorHistory.hpp"

namespace geos
{ // TODO document
namespace logging
{

class LogRouter {
public:

  static LogOutput s_main;


  LogRouter() = default;
  LogRouter( LogRouter && other ) = default;
  LogRouter( LogRouter const & other ) = delete;
  LogRouter & operator=( LogRouter const & other ) = delete;
  LogRouter & operator=( LogRouter && other ) = delete;


  void setHistoryMaxLevel( LogLevel maxLogLevel );

  void setTimeStepsToKeep( integer numberOfTimeSteps );

  // TODO : LoggingStrategy = logging outputs archetypes
  // void setOutputs( LoggingStrategy strategy );

  void addOutput( LoggingStrategy strategy );

private:

  LogErrorHistory m_errorHistory;

  LogLevel m_historyMaxLevel = LogLevel::Debug;

}

} /* namespace logging */
} /* namespace geos */

#endif /* GEOS_COMMON_LOGROUTER_HPP */