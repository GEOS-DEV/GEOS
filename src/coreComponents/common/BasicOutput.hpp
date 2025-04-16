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
 * @file BasicOutput.hpp
 */

#ifndef GEOS_COMMON_BASICOUTPUT_HPP
#define GEOS_COMMON_BASICOUTPUT_HPP

#include <string_view>
#include "common/logger/Logger.hpp"

/**
 * @brief Verify a opened stream is open.
 *        Add appropriate messages to the error report when the operation fails.
 * @tparam StreamType Type of the stream, typically std::ofstream or std::ifstream.
 * @tparam ErrorReporter Type of the error reporter callable, signature: void( string_view msg )
 * @param stream The stream to verify.
 * @param errorReporter Callable that handles error reporting by collecting messages
 * @param isNewlyOpened Flag indicating if the stream was just opened before this call
 */
template< typename StreamType, typename ErrorReporter >
void validateStream( StreamType const & stream, bool isNewlyOpened, ErrorReporter && errorReporter )
{
  if( !stream.good() )
  {
    errorReporter( isNewlyOpened ?
                   "Output stream failed to open.\nPossible reasons: File doesn't exist in path / permissions / locking issue." :
                   "Output stream is in invalid state for writing.\nPossible reasons: Stream closed / buffer corruption / previous operation failed." );
  }
}

/**
 * @brief Verify a non-critical opened stream is open.
 *        Add appropriate warning in the log if the operation fails.
 * @tparam ErrorReporter Type of the error reporter callable, signature: void( string_view msg )
 * @param stream The stream to verify.
 * @param isNewlyOpened Flag indicating if the stream was just opened before this call
 */
template< typename StreamType >
void validateStream( StreamType const & stream, bool isNewlyOpened )
{
  validateStream( stream, isNewlyOpened,
                  [&]( std::string_view msg ) { GEOS_WARNING( msg ); } );
}

/**
 * @brief Verify a critical opened stream is open.
 *        Add appropriate errors in the log if the operation fails (terminate GEOS).
 * @tparam ErrorReporter Type of the error reporter callable, signature: void( string_view msg )
 * @param stream The stream to verify.
 * @param isNewlyOpened Flag indicating if the stream was just opened before this call
 */
template< typename StreamType >
void validateCriticalStream( StreamType const & stream, bool isNewlyOpened )
{
  validateStream( stream, isNewlyOpened,
                  [&]( std::string_view msg ) { GEOS_ERROR( msg ); } );
}

/**
 * @brief Helper function to write content to an output stream.
 *        Adds appropriate messages to the error report when the operation fails.
 * @tparam ErrorReporter Type of the error reporter callable, signature: void( string_view msg )
 * @param outputStream The stream to write the content to.
 * @param content The string view containing data to be written.
 * @param errorReporter Callable that handles error reporting by collecting messages
 * @param isNewlyOpened Flag indicating if the stream was just opened before this call
 * @note this method may be moved in a common/BasicOutput.xpp as a lot of output stream errors are not verified.
 */
template< typename ErrorReporter >
void toStream( std::ostream & outputStream,
               std::string_view content,
               bool isNewlyOpened,
               ErrorReporter && errorReporter )
{
  validateStream( outputStream, isNewlyOpened, errorReporter );

  const auto startPos = outputStream.tellp();
  errno = 0;

  outputStream << content;

  if( outputStream.bad() )
  {
    const auto bytesWritten = outputStream.tellp() - startPos;
    errorReporter( GEOS_FMT( "I/O error occurred while writing content, written {} / {} bytes.\n"
                             "Possible reasons: Insufficient disk space / read-only filesystem / disk disconnection",
                             bytesWritten, content.size() ) );
  }
  else if( !content.empty() && startPos >= 0 && outputStream.tellp() <= startPos )
  {
    errorReporter( "Export completed but no data was written\nPossible reasons: Disk quota exceeded / streaming logical error." );
  }

  if( errno != 0 )
    errorReporter( GEOS_FMT( "\n{}", std::strerror( errno ) ) );
}

/**
 * @brief Helper function to write content to an output stream.
 *        Adds appropriate messages to the log when the operation fails.
 * @param outputStream The stream to write the content to.
 * @param content The string view containing data to be written.
 * @param streamName Name of the stream to use in potencial streaming errors
 * @param critical Flag indicating if any writing error is critical
 * @param isNewlyOpened Flag indicating if the stream was just opened before this call
 * @note this method may be moved in a common/BasicOutput.xpp as a lot of output stream errors are not verified.
 */
template< typename ErrorReporter >
void toStream( std::ostream & outputStream, std::string_view content,
               std::string_view streamName, bool critical, bool isNewlyOpened )
{
  std::string msgs;
  toStream( outputStream, content, isNewlyOpened,
            [&]( std::string_view msg ) { msgs += msg; } );
  if( critical )
  {
    GEOS_ERROR( GEOS_FMT( "Error while writing to '{}':\n{}", streamName, msgs ) );
  }
  else
  {
    GEOS_WARNING( GEOS_FMT( "Error while writing to '{}':\n{}", streamName, msgs ) );
  }
}

#endif // GEOS_COMMON_BASICOUTPUT_HPP
