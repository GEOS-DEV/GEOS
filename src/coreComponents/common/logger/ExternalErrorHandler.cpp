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
 * @file ErrorHandler.cpp
 */

#include "ExternalErrorHandler.hpp"
#include "common/logger/Logger.hpp"

#include <fstream>
#include <stdio.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>
#include <fcntl.h>

namespace geos
{

void OutputStreamDeviation::Pipe::create()
{
  if( ::pipe( fileDescriptorsArray ) == m_errorResult )
  {
    GEOS_WARNING( "Failed to create error pipe: " + getLastPosixErrorString() );
  }
}

bool OutputStreamDeviation::Pipe::setDescriptorInheritanceMode( bool inherit )
{
  bool r = true;

  // set O_CLOEXEC flag on both descriptors
  for( PosixId pipeEndFD : fileDescriptorsArray )
  {
    // get current flags for pipe end
    PosixId flags = ::fcntl( pipeEndFD, F_GETFD );
    if( flags == m_errorResult )
    {
      r = false;
    }
    else
    {
      // enable O_CLOEXEC flag if inheritance is disabled
      flags = inherit ? ( flags & ~O_CLOEXEC ) : ( flags | FD_CLOEXEC );

      // apply new flags
      if( ::fcntl( pipeEndFD, F_SETFD, flags ) == m_errorResult )
        r = false;
    }
  }

  return r;
}

bool OutputStreamDeviation::Pipe::redirectWriteEnd( int targetFd )
{
  return ::dup2( writeEnd(), targetFd ) != m_errorResult;
}

void OutputStreamDeviation::Pipe::closePipe()
{
  closePipeEnd( readEnd() );
  closePipeEnd( writeEnd() );
}

OutputStreamDeviation::OutputStreamDeviation( PosixId fileNo ):
  m_redirectedStream( fileNo ),
  m_originalStreamTarget( m_disabledPipeEnd )
{
  // Create pipe with appropriate flags
  m_deviationPipe.create();
  m_deviationPipe.setDescriptorInheritanceMode( false );

  if( !setPipeEndBlockingMode( m_deviationPipe.readEnd(), true ) )
  {
    m_deviationPipe.closePipe();
    GEOS_WARNING( "Failed to set pipe end non-blocking: " + getLastPosixErrorString() );
  }

  // backup original descriptor for reset after destruction
  m_originalStreamTarget = duplicateDescriptor( m_redirectedStream );
  if( m_originalStreamTarget == m_disabledPipeEnd )
  {
    m_deviationPipe.closePipe();
    GEOS_WARNING( "Failed to duplicate original descriptor: " + getLastPosixErrorString() );
  }

  // Redirect stderr to pipe
  if( !m_deviationPipe.redirectWriteEnd( m_redirectedStream ) )
  {
    m_deviationPipe.closePipe();
    closePipeEnd( m_originalStreamTarget );
    GEOS_WARNING( "Failed to redirect stream: " + getLastPosixErrorString() );
  }

  // Close the write end of the parent pipe: no longer useful as the source "file" will be written externally
  // we will only use the readEnd of the pipe to get the content of it
  closePipeEnd( m_deviationPipe.writeEnd() );

  prepareStreamingBuffer();
}

bool OutputStreamDeviation::setPipeEndBlockingMode( PosixId pipeEnd, bool nonBlocking )
{
  // get current flags for pipe end
  PosixId flags = ::fcntl( pipeEnd, F_GETFL );
  if( flags == m_errorResult )
    return false;

  // set flags depending on nonBlocking value
  flags = nonBlocking ? ( flags | O_NONBLOCK ) : ( flags & ~O_NONBLOCK );

  // apply new flags
  return ::fcntl( pipeEnd, F_SETFL, flags ) != m_errorResult;
}

int OutputStreamDeviation::duplicateDescriptor( int pipeEnd )
{
  int fdDuplicata = ::dup( pipeEnd );
  // if ::dup() fails, it should returns -1, which is the same value as m_disabledPipeEnd
  return fdDuplicata;
}

void OutputStreamDeviation::prepareStreamingBuffer()
{
  // pre-allocate the message processing buffer
  m_unprocessedData.reserve( 16384 );

  { // resize the pipe stream buffer from 16ko to 1mo to be able to process the biggest error messages as one bit (platform specific)
    #ifdef __linux__
    ::fcntl( m_deviationPipe.readEnd(), F_SETPIPE_SZ, 1048576 );
    #endif

    #ifdef __APPLE__
    PosixId bufsize = 1048576;
    ::setsockopt( m_deviationPipe.readEnd(), SOL_SOCKET, SO_RCVBUF,
                  &bufsize, sizeof(bufsize) );
    #endif
  }
}

void OutputStreamDeviation::closePipeEnd( PosixId & pipeEnd )
{
  if( pipeEnd != m_disabledPipeEnd )
  {
    ::close( pipeEnd );
    pipeEnd = m_disabledPipeEnd;
  }
}

string OutputStreamDeviation::getLastPosixErrorString()
{
  return std::string( strerror( errno ) );
}

OutputStreamDeviation::~OutputStreamDeviation()
{
  if( m_originalStreamTarget != m_disabledPipeEnd )
  {
    if( ::dup2( m_originalStreamTarget, m_redirectedStream ) == m_errorResult )
    {
      GEOS_WARNING( "Failed to restore pipe" );
    }

    ::close( m_originalStreamTarget );
    m_originalStreamTarget = m_disabledPipeEnd;
  }

  if( m_deviationPipe.readEnd() != m_disabledPipeEnd )
  {
    ::close( m_deviationPipe.readEnd() );
    m_deviationPipe.readEnd() = m_disabledPipeEnd;
  }
}


void OutputStreamDeviation::flush( OutputStreamDeviation::LineHandlingFunctor const & lineFunctor,
                                   std::string_view detectionLocation )
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
      lineFunctor( line, detectionLocation );
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


ExternalErrorHandler::ExternalErrorHandler():
  m_processErrorFunctor( ExternalErrorHandler::defaultErrorHandling )
{}

ExternalErrorHandler::~ExternalErrorHandler()
{
  enableStderrPipeDeviation( false );
}

ExternalErrorHandler & ExternalErrorHandler::instance()
{
  static ExternalErrorHandler instance;
  return instance;
}

void ExternalErrorHandler::enableStderrPipeDeviation( bool enable )
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

void ExternalErrorHandler::flush( std::string_view detectionLocation )
{
  if( m_stderrDeviation && m_processErrorFunctor )
  {
    m_stderrDeviation->flush( m_processErrorFunctor, detectionLocation );
  }
}

void ExternalErrorHandler::defaultErrorHandling( std::string_view errorMsg,
                                                 std::string_view detectionLocation )
{
  std::cout << "External error, detected" << detectionLocation << ": " << errorMsg << std::endl;
}

bool ExternalErrorHandler::isNotAnErrorMsg( string_view errorMsg )
{
  return errorMsg.find( "INFO|" ) != string_view::npos || // VTK info messages
         errorMsg.find( "[HYPREDRV]" ) == 0; // hypredrive debug traces are informational
}

} /* namespace geos */
