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

// Source includes
#include "ProblemManager.hpp"

// TPL includes
#include <conduit.hpp>

// System includes
#include <functional>
#include <chrono>
#include <ostream>

// Forward declaration of conduit::Node
namespace conduit
{
  class Node;
}

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

  GeosxState( GeosxState const & ) = delete;
  GeosxState( GeosxState && ) = delete;
  GeosxState & operator=( GeosxState const & ) = delete;
  GeosxState & operator=( GeosxState && ) = delete;

  bool initializeDataRepository();

  void applyInitialConditions();

  void run();

  State getState() const
  { return m_state; }

  conduit::Node & getRootConduitNode()
  {
    GEOSX_ERROR_IF( m_rootNode == nullptr, "Not initialized." );
    return *m_rootNode;
  }

  ProblemManager & getProblemManager()
  {
    GEOSX_ERROR_IF( m_problemManager == nullptr, "Not initialized." );
    return *m_problemManager;
  }

  FieldSpecificationManager & getFieldSpecificationManager();

  FunctionManager & getFunctionManager();

  std::chrono::system_clock::duration getInitTime() const
  { return m_initTime; }

  std::chrono::system_clock::duration getRunTime() const
  { return m_runTime; }

private:
  std::unique_ptr< conduit::Node > m_rootNode;
  std::unique_ptr< ProblemManager > m_problemManager;
  State m_state;
  std::chrono::system_clock::duration m_initTime;
  std::chrono::system_clock::duration m_runTime;
};



GeosxState & getGlobalState();

void setGlobalStateAccessor( std::function< GeosxState * () > const & accessor );


} // namespace geosx

#endif /* GEOSX_MANAGERS_GEOSXSTATE_HPP_ */
