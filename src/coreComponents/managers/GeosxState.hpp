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
#include "common/DataTypes.hpp"

// System includes
#include <functional>
#include <chrono>
#include <memory>

// Forward declaration of conduit::Node.
namespace conduit
{
class Node;
}

#if defined( GEOSX_USE_CALIPER )
//Forward declaration of cali::ConfigManager.
namespace cali
{
class ConfigManager;
}
#endif


namespace geosx
{

// Forward declarations.
namespace dataRepository
{
class Group;
}

class ProblemManager;
class FieldSpecificationManager;
class FunctionManager;
class CommunicationTools;
struct CommandLineOptions;

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
  GeosxState( std::unique_ptr< CommandLineOptions > && commandLineOptions );

  /**
   * @brief Destructor.
   * @note This is the default constructor but it must be specified in the `cpp`, otherwise the use of
   *   `std::unique_ptr` with forward declared types won't work.
   */
  ~GeosxState();

  GeosxState( GeosxState const & ) = delete;
  GeosxState( GeosxState && ) = delete;
  GeosxState & operator=( GeosxState const & ) = delete;
  GeosxState & operator=( GeosxState && ) = delete;

  bool initializeDataRepository();

  void applyInitialConditions();

  void run();

  State getState() const
  { return m_state; }

  CommandLineOptions const & getCommandLineOptions()
  {
    GEOSX_ERROR_IF( m_commandLineOptions == nullptr, "Not initialized." );
    return *m_commandLineOptions;
  }

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

  /**
   * @brief Return the @c ProblemManager but as a @c dataRepository::Group.
   * @return The @c ProblemManager but as a @c dataRepository::Group.
   * @note This is useful if you only need at @c dataRepository::Group and don't want to
   *   include @c ProblemManager.hpp.
   */
  dataRepository::Group & getProblemManagerAsGroup();

  FieldSpecificationManager & getFieldSpecificationManager();

  FunctionManager & getFunctionManager();

  CommunicationTools & getCommunicationTools()
  {
    GEOSX_ERROR_IF( m_commTools == nullptr, "Not initialized." );
    return *m_commTools;
  }

  std::chrono::system_clock::duration getInitTime() const
  { return m_initTime; }

  std::chrono::system_clock::duration getRunTime() const
  { return m_runTime; }

private:
  State m_state;
  std::unique_ptr< CommandLineOptions > m_commandLineOptions;
  std::unique_ptr< conduit::Node > m_rootNode;
  std::unique_ptr< ProblemManager > m_problemManager;
  std::unique_ptr< CommunicationTools > m_commTools;

#if defined( GEOSX_USE_CALIPER )
  std::unique_ptr< cali::ConfigManager > m_caliperManager;
#endif

  std::chrono::system_clock::duration m_initTime;
  std::chrono::system_clock::duration m_runTime;
};

GeosxState & getGlobalState();

} // namespace geosx

#endif /* GEOSX_MANAGERS_GEOSXSTATE_HPP_ */
