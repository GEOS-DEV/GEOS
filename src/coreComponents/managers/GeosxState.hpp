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
#include <ostream>

namespace geosx
{

std::string durationToString( std::chrono::system_clock::duration const duration );

enum class State
{
  UNINITIALIZED = 0,
  INITIALIZED = 1,
  READY_TO_RUN = 2,
  COMPLETED = 3
};

std::ostream & operator<<( std::ostream & os, State const state );

class GeosxState
{
public:
  GeosxState();

  bool initializeDataRepository();

  void applyInitialConditions();

  void run();

  State getState() const
  { return m_state; }

  dataRepository::Group & getProblemManager()
  { return m_problemManager; }

  std::chrono::system_clock::duration getInitTime() const
  { return m_initTime; }

  std::chrono::system_clock::duration getRunTime() const
  { return m_runTime; }

private:

  std::chrono::system_clock::time_point initialize();

  std::chrono::system_clock::time_point m_startTime;
  std::chrono::system_clock::duration m_initTime;
  std::chrono::system_clock::duration m_runTime;
  State m_state;
  ProblemManager m_problemManager;
};

} // namespace geosx

#endif /* GEOSX_MANAGERS_GEOSXSTATE_HPP_ */
