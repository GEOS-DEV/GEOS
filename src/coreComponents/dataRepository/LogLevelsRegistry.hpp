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
 * @file LogLevelsRegistry.hpp
 */

#ifndef GEOS_COMMON_LOGLEVELSREGISTRY_HPP
#define GEOS_COMMON_LOGLEVELSREGISTRY_HPP

#include "common/LogLevelsInfo.hpp"

namespace geos
{

class LogLevelsRegistry
{
public:
  
  /**
   * @brief Add a log description for a wrapper
   * @param level The minimum log level
   * @param description The description for the log level
   */
  void addEntry( integer level, std::string_view description );

  /**
   * @brief Construct the log level string description for a wrapper
   */
  string buildLogLevelDescription();

private:


  /// Map for building the log level string for each wrapper
  /// key : a logLevel condition, values : a set of description for a corresponding loglevel
  std::map< integer, std::vector< std::string > > m_logLevelsDescriptions;

};

/**
 * @brief Verify if a log level is active
 * @tparam LOG_LEVEL_INFO The structure containing log level information.
 * @param level Log level to be checked.
 * @return `true` if the log level is active, `false` otherwise.
 * @pre `LOG_LEVEL_INFO` must satisfy `logInfo::is_log_level_info`.
 *
 */
template< typename LOG_LEVEL_INFO >
std::enable_if_t< logInfo::is_log_level_info< LOG_LEVEL_INFO >, bool >
isLogLevelActive( int level )
{
  return level >= LOG_LEVEL_INFO::getMinLogLevel();
}
}

#endif
