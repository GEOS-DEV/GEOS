/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "GeosxState.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

// TPL includes
#include <conduit.hpp>

#if defined( GEOSX_USE_CALIPER )
  #include <caliper/cali-manager.h>
#endif

// System includes
#include <ostream>

namespace geosx
{

GeosxState * currentGlobalState = nullptr;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState & getGlobalState()
{
  GEOSX_ERROR_IF( currentGlobalState == nullptr,
                  "The state has not been created." );

  return *currentGlobalState;
}


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
string durationToString( std::chrono::system_clock::duration const duration )
{
  // If we want to print HH::MM::SS (maybe in addition to seconds-only):
  // return GEOSX_FMT( "{:%T}", duration );
  double const seconds = std::chrono::duration_cast< std::chrono::milliseconds >( duration ).count() / 1000.0;
  return GEOSX_FMT( "{:>20.3f}s", seconds );
}

std::ostream & operator<<( std::ostream & os, State const state )
{
  if( state == State::UNINITIALIZED )
  {
    return os << "State::UNINITIALIZED";
  }
  if( state == State::INITIALIZED )
  {
    return os << "State::INITIALIZED";
  }
  if( state == State::READY_TO_RUN )
  {
    return os << "State::READY_TO_RUN";
  }
  if( state == State::COMPLETED )
  {
    return os << "State::COMPLETED";
  }

  GEOSX_ERROR( "Unrecognized state. The integral value is: " << static_cast< int >( state ) );
  return os;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::GeosxState( std::unique_ptr< CommandLineOptions > && commandLineOptions ):
  m_state( State::UNINITIALIZED ),
  m_commandLineOptions( std::move( commandLineOptions ) ),
  m_rootNode( std::make_unique< conduit::Node >() ),
  m_problemManager( nullptr ),
  m_commTools( std::make_unique< CommunicationTools >() ),
#if defined( GEOSX_USE_CALIPER )
  m_caliperManager( std::make_unique< cali::ConfigManager >() ),
#endif
  m_initTime(),
  m_runTime()
{
  Timer timer( m_initTime );

#if defined( GEOSX_USE_CALIPER )
  setupCaliper( *m_caliperManager, getCommandLineOptions() );
#endif

  string restartFileName;
  if( ProblemManager::parseRestart( restartFileName, getCommandLineOptions() ) )
  {
    GEOSX_LOG_RANK_0( "Loading restart file " << restartFileName );
    dataRepository::loadTree( restartFileName, getRootConduitNode() );
  }

  m_problemManager = std::make_unique< ProblemManager >( getRootConduitNode() );

  GEOSX_ERROR_IF( currentGlobalState != nullptr, "Only one state can exist at a time." );
  currentGlobalState = this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::~GeosxState()
{
#if defined( GEOSX_USE_CALIPER )
  m_caliperManager->flush();
#endif

  GEOSX_ERROR_IF( currentGlobalState != this, "This shouldn't be possible." );
  currentGlobalState = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GeosxState::initializeDataRepository()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOSX_THROW_IF_NE( m_state, State::UNINITIALIZED, std::logic_error );

  getProblemManager().parseCommandLineInput();

  if( !getProblemManager().getSchemaFileName().empty() )
  {
    getProblemManager().generateDocumentation();
    m_state = State::INITIALIZED;
    return false;
  }

  getProblemManager().parseInputFile();
  getProblemManager().problemSetup();

  m_state = State::INITIALIZED;

  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::applyInitialConditions()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOSX_THROW_IF_NE( m_state, State::INITIALIZED, std::logic_error );

  getProblemManager().applyInitialConditions();

  if( getCommandLineOptions().beginFromRestart )
  {
    getProblemManager().readRestartOverwrite();
  }

  m_state = State::READY_TO_RUN;
  MpiWrapper::barrier( MPI_COMM_GEOSX );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::run()
{
  GEOSX_MARK_FUNCTION;
  Timer timer( m_runTime );

  GEOSX_THROW_IF_NE( m_state, State::READY_TO_RUN, std::logic_error );

  if( !getProblemManager().runSimulation() )
  {
    m_state = State::COMPLETED;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dataRepository::Group & GeosxState::getProblemManagerAsGroup()
{ return getProblemManager(); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FieldSpecificationManager & GeosxState::getFieldSpecificationManager()
{ return getProblemManager().getFieldSpecificationManager(); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FunctionManager & GeosxState::getFunctionManager()
{ return getProblemManager().getFunctionManager(); }

} // namespace geosx
