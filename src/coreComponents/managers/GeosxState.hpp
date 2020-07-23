/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GeosxState.hpp
 */

#ifndef GEOSX_MANAGERS_GEOSXSTATE_HPP_
#define GEOSX_MANAGERS_GEOSXSTATE_HPP_

#include "ProblemManager.hpp"

#include <chrono>

namespace geosx
{

std::string durationToString( std::chrono::system_clock::duration const duration );

class GeosxState
{
public:
  GeosxState();

  bool initializeDataRepository();

  bool run();

  dataRepository::Group * getGroupByPath( std::string const & path )
  { return m_problemManager.GetGroupByPath( path ); }

  std::chrono::system_clock::duration getInitTime() const
  { return m_initTime; }

  std::chrono::system_clock::duration getRunTime() const
  { return m_runTime; }

private:

  std::chrono::system_clock::time_point initialize();

  std::chrono::system_clock::time_point m_startTime;
  std::chrono::system_clock::duration m_initTime;
  std::chrono::system_clock::duration m_runTime;
  ProblemManager m_problemManager;
};

} // namespace geosx

#endif /* GEOSX_MANAGERS_GEOSXSTATE_HPP_ */
