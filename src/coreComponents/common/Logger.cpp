/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Logger.cpp
 */

// Source includes
#include "Logger.hpp"
#include "Path.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geos
{

/**
 * @brief Insert an exception message in another one.
 * @param originalMsg original exception message (i.e. thrown from LVARRAY_THROW or GEOSX_THROW)
 * @param msgToInsert message to insert at the top of the originalMsg
 */
std::string InsertExMsg( std::string const & originalMsg, std::string const & msgToInsert )
{
  std::string newMsg( originalMsg );

  size_t insertPos = 0;
  // for readability purposes, we try to insert the message after the "***** Rank N: " or after "***** " instead of at the top.
  static auto constexpr rankLogStart =  "***** Rank ";
  static auto constexpr rankLogEnd =  ": ";
  static auto constexpr simpleLogStart =  "***** ";
  if( ( insertPos = newMsg.find( rankLogStart ) ) != std::string::npos )
  {
    insertPos = newMsg.find( rankLogEnd, insertPos + stringutilities::cstrlen( rankLogStart ) )
                + stringutilities::cstrlen( rankLogEnd );
  }
  else if( ( insertPos = newMsg.find_last_of( simpleLogStart ) ) != std::string::npos )
  {
    insertPos += stringutilities::cstrlen( simpleLogStart );
  }
  newMsg.insert( insertPos, msgToInsert );
  return newMsg;
}

InputError::InputError( std::exception const & subException, std::string const & msgToInsert ):
  std::runtime_error( InsertExMsg( subException.what(), msgToInsert ) )
{}


Logger logger;

Logger::Logger():
  rank( 0 ),
  ranksCount( 1 ),
  rankMsgPrefix( "" )
{
  reset();
}
void Logger::reset()
{
  globalLogLevel = 0;
  minLogLevel = 0;
  maxLogLevel = 5;

  setOutputToStd();
}
void Logger::initMpi( MPI_Comm mpiComm )
{
#ifdef GEOSX_USE_MPI
  comm = mpiComm;
  MPI_Comm_rank( mpiComm, &rank );
  MPI_Comm_size( mpiComm, &ranksCount );
  rankMsgPrefix = ranksCount > 0  ? GEOS_FMT( "Rank {}: ", rank ) : "";// TODO : choisir si l'on garde cette ternaire
#else
  GEOS_ERROR( "Trying to initialize MPI in serial build." );
#endif
}

void Logger::setOutputToFile( const std::string & rankOutputDir )
{
#ifdef GEOSX_USE_MPI
  MPI_Barrier( comm );
#endif
  std::string outputFolder;
  if( rankOutputDir != "" )
  {
    makeDirsForPath( rankOutputDir );
    outputFolder = rankOutputDir + '/';
  }

  fileOutStream = std::make_unique< std::ofstream >( GEOS_FMT( "{}rank_{}.out", outputFolder, rank ) );
  outStream = fileOutStream.get();
}
void Logger::setOutputToStd()
{
#ifdef GEOSX_USE_MPI
  if( ranksCount > 1 )
  {
    MPI_Barrier( comm );
  }
#endif

  fileOutStream = nullptr;
  outStream = &std::cout;
}

void Logger::setGlobalLogLevel( int value )
{ globalLogLevel = value; }
void Logger::setMinLogLevel( int value )
{ minLogLevel = value; }
void Logger::setMaxLogLevel( int value )
{ maxLogLevel = value; }


} // namespace geos
