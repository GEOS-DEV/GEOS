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

// Source includes
#include "GeosxState.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/initialization.hpp"
#include "mpiCommunications/CommunicationTools.hpp"

// TPL includes
#include <conduit.hpp>

// System includes
#include <ostream>

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

std::ostream & operator<<( std::ostream & os, State const state )
{
  if ( state == State::UNINITIALIZED )
  { return os << "State::UNINITIALIZED"; }
  if ( state == State::INITIALIZED )
  { return os << "State::INITIALIZED"; }
  if ( state == State::READY_TO_RUN )
  { return os << "State::READY_TO_RUN"; }
  if ( state == State::COMPLETED )
  { return os << "State::COMPLETED"; }

  GEOSX_ERROR( "Unrecognized state. The integral value is: " << static_cast< int >( state ) );
  return os;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::GeosxState():
  m_state( State::UNINITIALIZED ),
  m_rootNode( new conduit::Node ),
  m_problemManager( nullptr ),
  m_commTools( std::make_unique< CommunicationTools >() ),
  m_initTime(),
  m_runTime()
{
  Timer timer( m_initTime );

  std::string restartFileName;
  if( ProblemManager::ParseRestart( restartFileName ) )
  {
    GEOSX_LOG_RANK_0( "Loading restart file " << restartFileName );
    dataRepository::loadTree( restartFileName, *m_rootNode );
  }

  m_problemManager = std::make_unique< ProblemManager >( "Problem", *m_rootNode );
  GEOSX_LOG_VAR( m_problemManager->getPath() );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::~GeosxState() = default;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GeosxState::initializeDataRepository()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOSX_ERROR_IF_NE( m_state, State::UNINITIALIZED );

  getProblemManager().ParseCommandLineInput();

  if( !getProblemManager().getSchemaFileName().empty() )
  {
    getProblemManager().GenerateDocumentation();
    m_state = State::INITIALIZED;
    return false;
  }

  getProblemManager().ParseInputFile();
  getProblemManager().ProblemSetup();

  m_state = State::INITIALIZED;

  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::applyInitialConditions()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_initTime );
  
  GEOSX_ERROR_IF_NE( m_state, State::INITIALIZED );

  getProblemManager().ApplyInitialConditions();

  if ( getCommandLineOptions().beginFromRestart )
  { getProblemManager().ReadRestartOverwrite(); }

  m_state = State::READY_TO_RUN;
  MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::run()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_runTime );

  GEOSX_ERROR_IF_NE( m_state, State::READY_TO_RUN );

  if ( !getProblemManager().RunSimulation() )
  {
    m_state = State::COMPLETED;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FieldSpecificationManager & GeosxState::getFieldSpecificationManager()
{ return getProblemManager().getFieldSpecificationManager(); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FunctionManager & GeosxState::getFunctionManager()
{ return getProblemManager().getFunctionManager(); }


static std::function< GeosxState * () > s_stateAccessor;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState & getGlobalState()
{
  GEOSX_ERROR_IF( s_stateAccessor == nullptr,
                  "A state accessor has not been set, set one with setGlobalStateAccessor." );

  GeosxState * const state = s_stateAccessor();
  GEOSX_ERROR_IF( state == nullptr,
                  "The state has not been created." );

  return *state;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setGlobalStateAccessor( std::function< GeosxState * () > const & accessor )
{ s_stateAccessor = accessor; }

} // namespace geosx

