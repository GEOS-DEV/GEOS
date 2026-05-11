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
 * @file RestartTiming.hpp
 */

#ifndef GEOS_DATAREPOSITORY_RESTARTTIMING_HPP_
#define GEOS_DATAREPOSITORY_RESTARTTIMING_HPP_

#include "common/DataTypes.hpp"
#include "common/Stopwatch.hpp"
#include "common/logger/Logger.hpp"

#include <cstdlib>
#include <cstring>

namespace geos
{
namespace dataRepository
{
namespace restartTiming
{

/**
 * @brief Return whether restart timing diagnostics are enabled.
 *
 * Set GEOS_RESTART_TIMING to any value other than 0, false, or off to enable
 * begin/end/progress messages during restart loading.
 */
inline bool enabled()
{
  static bool const value = []()
  {
    char const * const env = std::getenv( "GEOS_RESTART_TIMING" );
    return env != nullptr &&
           env[0] != '\0' &&
           std::strcmp( env, "0" ) != 0 &&
           std::strcmp( env, "false" ) != 0 &&
           std::strcmp( env, "FALSE" ) != 0 &&
           std::strcmp( env, "off" ) != 0 &&
           std::strcmp( env, "OFF" ) != 0;
  }();
  return value;
}

/**
 * @brief Scoped begin/end logger for restart timing diagnostics.
 */
class ScopedTimer
{
public:

  explicit ScopedTimer( string label ):
    m_label( std::move( label ) ),
    m_enabled( enabled() )
  {
    if( m_enabled )
    {
      GEOS_LOG_RANK( "RestartTiming BEGIN " << m_label );
    }
  }

  ~ScopedTimer()
  {
    if( m_enabled )
    {
      GEOS_LOG_RANK( "RestartTiming END   " << m_label << " elapsed=" << m_timer.elapsedTime() << " s" );
    }
  }

private:

  string m_label;
  Stopwatch m_timer;
  bool const m_enabled;
};

inline void logProgress( string const & label, real64 const elapsed )
{
  if( enabled() )
  {
    GEOS_LOG_RANK( "RestartTiming PROGRESS " << label << " elapsed=" << elapsed << " s" );
  }
}

} /* namespace restartTiming */
} /* namespace dataRepository */
} /* namespace geos */

#endif /* GEOS_DATAREPOSITORY_RESTARTTIMING_HPP_ */
