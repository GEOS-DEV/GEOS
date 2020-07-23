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

std::string durationToString( std::chrono::system_clock::duration const duration )
{
  double const seconds = std::chrono::duration_cast< std::chrono::milliseconds >( duration ).count() / 1000.0;

  char buffer[ 32 ];
  std::snprintf( buffer, 32, "%20.3fs", seconds );

  return buffer;
}

GeosxState::GeosxState():
  m_startTime( initialize() ),
  m_initTime(),
  m_runTime(),
  m_problemManager( "Problem", nullptr )
{
  m_initTime = std::chrono::system_clock::now() - m_startTime;
}

bool GeosxState::initializeDataRepository()
{
  GEOSX_MARK_FUNCTION;

  std::chrono::system_clock::time_point const begin = std::chrono::system_clock::now();

  m_problemManager.ParseCommandLineInput();

  if( !m_problemManager.getSchemaFileName().empty() )
  {
    m_problemManager.GenerateDocumentation();
    return false;
  }

  m_problemManager.ParseInputFile();
  m_problemManager.ProblemSetup();

  if ( getCommandLineOptions().beginFromRestart )
  { m_problemManager.ReadRestartOverwrite(); }

  MpiWrapper::Barrier( MPI_COMM_GEOSX );

  m_initTime += std::chrono::system_clock::now() - begin;

  return true;
}

bool GeosxState::run()
{
  GEOSX_MARK_FUNCTION;
  std::chrono::system_clock::time_point const begin = std::chrono::system_clock::now();
  bool const earlyReturn = m_problemManager.RunSimulation();
  m_runTime += std::chrono::system_clock::now() - begin;
  return earlyReturn;
}

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

