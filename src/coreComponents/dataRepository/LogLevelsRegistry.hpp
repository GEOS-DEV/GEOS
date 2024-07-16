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

#include "common/LogLevels.hpp"

class LogLevelsRegistry
{
public:

  void addEntry( string level, string description );

  /**
   * @brief Construct the log level string description for a wrapper
   */
  void buildLogLevelDescription( std::map< std::string, std::vector< std::string > > const & logLevelsDescriptions );

private:


  /// Map for building the log level string for each wrapper
  /// key : a logLevel condition, values : a set of description for a corresponding loglevel
  std::map< std::string, std::vector< std::string > > m_logLevelsDescriptions;

}

template< typename LOG_LEVEL_INFO >
std::enable_if_t< is_log_level_info< LOG_LEVEL_INFO >, bool >
isLogLevelActive( int level )
{
  return level >= LOG_LEVEL_INFO::getMinLogLevel();
}
