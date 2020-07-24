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

#include "GeosxState.hpp"
#include "managers/initialization.hpp"

namespace geosx
{

/**
 * @class Timer
 * @brief Object that times the duration of its existence.
 */
class Timer
{
public:

  /**
   * @brief Constructor. The time the object is alive is added to @p duration.
   * @param duration A reference to the duration to add to.
   */
  Timer( std::chrono::system_clock::duration & duration ):
    m_start( std::chrono::system_clock::now() ),
    m_duration( duration )
  {}

  /// Destructor. Adds to the referenced duration.
  ~Timer()
  { m_duration += std::chrono::system_clock::now() - m_start; }

private:
  /// The time at which this object was constructed.
  std::chrono::system_clock::time_point const m_start;
  /// A reference to the duration to add to.
  std::chrono::system_clock::duration & m_duration;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string durationToString( std::chrono::system_clock::duration const duration )
{
  double const seconds = std::chrono::duration_cast< std::chrono::milliseconds >( duration ).count() / 1000.0;

  char buffer[ 32 ];
  std::snprintf( buffer, 32, "%20.3fs", seconds );

  return buffer;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::GeosxState():
  m_startTime( initialize() ),
  m_initTime(),
  m_runTime(),
  m_state( State::UNINITIALIZED ),
  m_problemManager( "Problem", nullptr )
{
  m_initTime = std::chrono::system_clock::now() - m_startTime;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GeosxState::initializeDataRepository()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOSX_ERROR_IF_NE( m_state, State::UNINITIALIZED );

  m_problemManager.ParseCommandLineInput();

  if( !m_problemManager.getSchemaFileName().empty() )
  {
    m_problemManager.GenerateDocumentation();
    m_state = State::INITIALIZED;
    return false;
  }

  m_problemManager.ParseInputFile();
  m_problemManager.ProblemSetup();

  m_state = State::INITIALIZED;

  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataRepository::Group * GeosxState::getGroupByPath( std::string const & path )
{ 
  GEOSX_ERROR_IF_EQ( m_state, State::UNINITIALIZED );
  return m_problemManager.GetGroupByPath( path );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::applyInitialConditions()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_initTime );
  
  GEOSX_ERROR_IF_NE( m_state, State::INITIALIZED );

  m_problemManager.ApplyInitialConditions();

  if ( getCommandLineOptions().beginFromRestart )
  { m_problemManager.ReadRestartOverwrite(); }

  m_state = State::READY_TO_RUN;
  MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::run()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_runTime );

  GEOSX_ERROR_IF_NE( m_state, State::READY_TO_RUN );

  if ( !m_problemManager.RunSimulation() )
  {
    m_state = State::COMPLETED;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::chrono::system_clock::time_point GeosxState::initialize()
{
  GEOSX_MARK_FUNCTION;
  std::chrono::system_clock::time_point const startTime = std::chrono::system_clock::now();

  std::string restartFileName;
  if( ProblemManager::ParseRestart( restartFileName ) )
  {
    GEOSX_LOG_RANK_0( "Loading restart file " << restartFileName );
    dataRepository::loadTree( restartFileName );
  }

  return startTime;
}

} // namespace geosx

