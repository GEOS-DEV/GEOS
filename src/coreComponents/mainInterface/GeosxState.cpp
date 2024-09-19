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
#include "dataRepository/Utilities.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/initialization.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "common/Timer.hpp"

// TPL includes
#include <conduit.hpp>

#if defined( GEOSX_USE_CALIPER )
  #include <caliper/cali-manager.h>
#endif

// System includes
#include <ostream>

namespace geos
{

GeosxState * currentGlobalState = nullptr;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState & getGlobalState()
{
  GEOS_ERROR_IF( currentGlobalState == nullptr,
                 "The state has not been created." );

  return *currentGlobalState;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
string durationToString( std::chrono::system_clock::duration const duration )
{
  // If we want to print HH::MM::SS (maybe in addition to seconds-only):
  // return GEOS_FMT( "{:%T}", duration );
  double const seconds = std::chrono::duration_cast< std::chrono::milliseconds >( duration ).count() / 1000.0;
  return GEOS_FMT( "{:>20.3f}s", seconds );
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

  GEOS_ERROR( "Unrecognized state. The integral value is: " << static_cast< int >( state ) );
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
    GEOS_LOG_RANK_0( "Loading restart file " << restartFileName );
    dataRepository::loadTree( restartFileName, getRootConduitNode() );
  }

  m_problemManager = std::make_unique< ProblemManager >( getRootConduitNode() );

  GEOS_ERROR_IF( currentGlobalState != nullptr, "Only one state can exist at a time." );
  currentGlobalState = this;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GeosxState::~GeosxState()
{
#if defined( GEOSX_USE_CALIPER )
  m_caliperManager->flush();
#endif

  GEOS_ERROR_IF( currentGlobalState != this, "This shouldn't be possible." );
  currentGlobalState = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GeosxState::initializeDataRepository()
{
  GEOS_MARK_FUNCTION;
  if ( ! initializeInputParsing( ) )
  {
    return false; // Early exit if initialization fails or is complete at the parsing stage.
  }

  // setup the problem after input parsing
  auto & problemManager = getProblemManager();
  problemManager.problemSetup();
  m_state = State::INITIALIZED;

  if( m_commandLineOptions->printMemoryUsage >= 0.0 )
  {
    dataRepository::printMemoryAllocation( problemManager, 0, m_commandLineOptions->printMemoryUsage );
  }

  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool GeosxState::initializeInputParsing()
{
  GEOS_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOS_THROW_IF_NE(m_state, State::UNINITIALIZED, std::logic_error);

  auto& problemManager = getProblemManager();

  problemManager.parseCommandLineInput(); // this currently merges multiple input files

  if (!problemManager.getSchemaFileName().empty())
  {
    problemManager.generateDocumentation();
    m_state = State::INITIALIZED;
    return false;
  }

  parseInputFiles();
  return false;

  if (m_commandLineOptions->preprocessOnly)
  {
    return false;
  }

  return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::parseInputFiles()
{
  using namespace dataRepository::inputProcessing;
  using document_type = inputParsing::input_document_type;
  auto & problemManager = getProblemManager();

  // Locate mergable Groups
  // GEOS_LOG_RANK_0( " -- generating skeleton" );
  // problemManager.generateDataStructureSkeleton( 0 );
  // GEOS_LOG_RANK_0( " -- skeleton generated" );
  // std::vector< dataRepository::Group const * > containerGroups;
  // problemManager.discoverGroupsRecursively( containerGroups, []( dataRepository::Group const & group ) { return group.numWrappers() == 0 && group.numSubGroups() > 0; } );
  // GEOS_LOG_RANK_0( " -- containers identified" );
  std::set< string > mergableNodes;
  // for( dataRepository::Group const * group : containerGroups )
  // {
  //   std::cout << group->getName() << std::endl;
  //   mergableNodes.insert( group->getCatalogName() );
  // }
  // containerGroups.clear( );
  // GEOS_LOG_RANK_0( " -- names stored" );
  // // problemManager.deregisterAllRecursive( );
  // GEOS_LOG_RANK_0( " -- skeleton destroyed" );

  std::set< string > allowedDeviations = { DomainPartition::groupKeysStruct::constitutiveManagerString(), MeshLevel::groupStructKeys::elemManagerString(), MeshLevel::groupStructKeys::particleManagerString() };

  // we get the input file name from the problem manager instead of the m_commandLineOptions because we have merged multiple -i options into a single file at this point
  string const & inputFileName = problemManager.getGroup< dataRepository::Group >( problemManager.groupKeys.commandLine ).getReference< string >( problemManager.viewKeys.inputFileName );
  document_type inputDocument;
  typename document_type::read_result_type readResult = inputDocument.loadFile( inputFileName, true );
  GEOS_THROW_IF( !readResult, GEOS_FMT( "Errors found while parsing input file\nDescription: {}\nOffset: {}",
                                       readResult.description(), readResult.offset ), InputError );

  // Extract the problem node and begin processing the user inputs
  typename document_type::node_type docRoot = inputDocument.getFirstChild();
  inputDocument.processIncludes( docRoot, 0 );
  if( m_commandLineOptions->preprocessOnly )
  {
    // by performing Declaration and TerseSyntax expansion we allow preprocessing to capture e.g. the application of Deprecations
    PreprocessOnly processor( inputDocument, mergableNodes, allowedDeviations );
    processor.execute( problemManager, docRoot );
    // Remove all the created Groups from preprocessing phases (only Declaration actually creates Groups)
    problemManager.deregisterAllRecursive( );
    if ( MpiWrapper::commRank() == 0 )
    {
      string const & outputDirectory = m_commandLineOptions->outputDirectory;
      inputDocument.saveFile( outputDirectory + "/" + inputFileName );
    }
  }
  else
  {
    AllProcessingPhases processor( inputDocument, mergableNodes, allowedDeviations );
    processor.execute( problemManager, docRoot );
    problemManager.processSchemaDeviations( inputDocument, mergableNodes );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GeosxState::applyInitialConditions()
{
  GEOS_MARK_FUNCTION;
  Timer timer( m_initTime );

  GEOS_THROW_IF_NE( m_state, State::INITIALIZED, std::logic_error );

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
  GEOS_MARK_FUNCTION;
  Timer timer( m_runTime );

  GEOS_THROW_IF_NE( m_state, State::READY_TO_RUN, std::logic_error );

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

} // namespace geos
