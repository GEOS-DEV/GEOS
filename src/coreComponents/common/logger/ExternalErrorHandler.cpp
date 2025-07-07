/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ErrorHandling.cpp
 */

#include "ExternalErrorHandler.hpp"
#include "common/logger/Logger.hpp"

#include <fstream>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

namespace geos
{

OutputStreamDeviation::OutputStreamDeviation( int fileNo ):
  m_redirectedStream( fileNo ),
  m_originalStreamTarget( m_disabledPipe )
{
  //TODO : test, factorize, create PR
  // Create pipe with appropriate flags
  if( ::pipe( m_deviationPipe.fileDescriptorsArray ) == m_errorResult )
  {
    GEOS_WARNING( "Failed to create error pipe: " + std::string( ::strerror( errno ) ) );
  }

  // set O_CLOEXEC on both descriptors
  for( int fd : m_deviationPipe.fileDescriptorsArray )
  {
    int flags = ::fcntl( fd, F_GETFD );

    if( flags == m_errorResult || ::fcntl( fd, F_SETFD, flags | FD_CLOEXEC ) == m_errorResult )
    {
      ::close( m_deviationPipe.readEnd() );
      ::close( m_deviationPipe.writeEnd() );
      GEOS_WARNING( "Failed to set CLOEXEC: " + std::string( strerror( errno ) ) );
    }
  }

  // set read end to non blocking
  int flags = ::fcntl( m_deviationPipe.readEnd(), F_GETFL );
  if( flags == m_errorResult ||
      ::fcntl( m_deviationPipe.readEnd(), F_SETFL, flags | O_NONBLOCK ) == m_errorResult )
  {
    ::close( m_deviationPipe.readEnd() );
    ::close( m_deviationPipe.writeEnd() );
    GEOS_WARNING( "Failed to set non-blocking: " + std::string( strerror( errno ) ) );
  }

  // backup original descriptor
  m_originalStreamTarget = ::dup( m_redirectedStream );
  if( m_originalStreamTarget == m_errorResult )
  {
    ::close( m_deviationPipe.readEnd() );
    ::close( m_deviationPipe.writeEnd() );
    GEOS_WARNING( "Failed to duplicate original descriptor: " + std::string( strerror( errno ) ) );
  }

  // Redirect stderr to pipe
  if( ::dup2( m_deviationPipe.writeEnd(), m_redirectedStream ) == m_errorResult )
  {
    ::close( m_originalStreamTarget );
    ::close( m_deviationPipe.readEnd() );
    ::close( m_deviationPipe.writeEnd() );
    m_originalStreamTarget = m_disabledPipe;
    GEOS_WARNING( "Failed to redirect stream: " + std::string( strerror( errno ) ) );
  }

  // Close the write end of the parent pipe (we write through stderr now)
  ::close( m_deviationPipe.writeEnd() );
  m_deviationPipe.writeEnd() = m_disabledPipe;

  // various optimizations
  m_unprocessedData.reserve( 16384 );

  #ifdef __linux__
  ::fcntl( m_deviationPipe.readEnd(), F_SETPIPE_SZ, 1048576 );
  #endif

  #ifdef __APPLE__
  int bufsize = 65536;
  ::setsockopt( m_deviationPipe.readEnd(), SOL_SOCKET, SO_RCVBUF,
                &bufsize, sizeof(bufsize) );
  #endif
}

OutputStreamDeviation::~OutputStreamDeviation()
{
  if( m_originalStreamTarget != m_disabledPipe )
  {
    if( ::dup2( m_originalStreamTarget, m_redirectedStream ) == m_errorResult )
    {
      GEOS_WARNING( "Failed to restore pipe" );
    }

    ::close( m_originalStreamTarget );
    m_originalStreamTarget = m_disabledPipe;
  }

  if( m_deviationPipe.readEnd() != m_disabledPipe )
  {
    ::close( m_deviationPipe.readEnd() );
    m_deviationPipe.readEnd() = m_disabledPipe;
  }
}


void OutputStreamDeviation::flush( OutputStreamDeviation::LineHandlingFunctor const & lineFunctor )
{
  std::array< char, 8192 > readBuffer;
  ssize_t bytesRead;

  // read all pending data from the original stream & add it in the text buffer to process
  while( ( bytesRead = ::read( m_deviationPipe.readEnd(),
                               readBuffer.data(),
                               readBuffer.size() ) ) > 0 )
  {
    m_unprocessedData.append( readBuffer.data(), bytesRead );
  }

  {   // process each full lines
    size_t lineStart = 0;
    size_t lineEnd = 0;

    while((lineEnd = m_unprocessedData.find( '\n', lineStart )) != std::string::npos )
    {
      std::string_view line = std::string_view( m_unprocessedData.data() + lineStart,
                                                lineEnd - lineStart );
      lineFunctor( line );
      lineStart = lineEnd + 1;
    }

    // keep last line residual if it exists (incomplete line)
    if( lineStart < m_unprocessedData.size() )
    {
      m_unprocessedData.erase( 0, lineStart );
    }
    else
    {
      m_unprocessedData.clear();
    }
  }
}


void defaultErrorHandling( std::string_view errorMsg )
{
  std::cout << "External error: " << errorMsg << std::endl;
}

ExternalErrorHandler::ExternalErrorHandler():
  m_processErrorFunctor( defaultErrorHandling )
{}

ExternalErrorHandler::~ExternalErrorHandler()
{
  enableStderrPipe( false );
}

ExternalErrorHandler & ExternalErrorHandler::instance()
{
  static ExternalErrorHandler instance;
  return instance;
}

void ExternalErrorHandler::enableStderrPipe( bool enable )
{
  if( enable && !m_stderrDeviation )
  {
    m_stderrDeviation = std::make_unique< OutputStreamDeviation >( STDERR_FILENO );
  }
  else if( !enable && m_stderrDeviation )
  {
    m_stderrDeviation = nullptr;
  }
}

void ExternalErrorHandler::flush()
{
  if( m_stderrDeviation && m_processErrorFunctor )
  {
    m_stderrDeviation->flush( m_processErrorFunctor );
  }
}

} /* namespace geos */
