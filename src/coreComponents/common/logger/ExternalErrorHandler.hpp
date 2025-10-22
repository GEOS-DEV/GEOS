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
 * @file ErrorHandling.hpp
 */

#ifndef LOGGER_EXTERNALERRORHANDLER_HPP
#define LOGGER_EXTERNALERRORHANDLER_HPP

#include "ErrorHandling.hpp"

#include <functional>
#include <string_view>
#include <string>
#include <memory>

namespace geos
{

/**
 * @brief This class implements pipe redirection to allow to capture and process externally streamed messages
 */
class OutputStreamDeviation
{
public:

  /**
   * @brief A functor executed for each independant lines to process, taking the line as the 1st
   *        string_view parameter, and the detectionLocation as the 2nd.
   */
  using LineHandlingFunctor = std::function< void (std::string_view, std::string_view) >;

  /**
   * @brief Construct and enable a new pipe redirection.
   * @param fileNo The file descriptor number, as returned by fileno() or can be one of the
   *               following: "STDOUT_FILENO" (1), "STDERR_FILENO" (2).
   */
  OutputStreamDeviation( int fileNo );

  /**
   * @brief Destroy the OutputStreamDeviation object, restoring the original pipe state.
   */
  ~OutputStreamDeviation();

  /**
   * @brief Flush the buffer from the original output pipe in a string, allowing to log it where needed.
   * @param lineProcessingFunctor see LineHandlingFunctor.
   * @param detectionLocation A label to describe when the flush() operation is being made, thus
   *                          explaining to the user when the error has been detected.
   */
  void flush( LineHandlingFunctor const & lineProcessingFunctor, std::string_view detectionLocation );

private:
  struct Pipe
  {
    int fileDescriptorsArray[2];

    int & readEnd() { return fileDescriptorsArray[0]; }
    int & writeEnd() { return fileDescriptorsArray[1]; }
  };

  /// a special value to represant a disabled / not existing pipe.
  static constexpr int m_disabledPipe = -1;

  /// error values from POSIX functions.
  static constexpr int m_errorResult = -1;

  /// the original pipe to deviate
  int m_redirectedStream;

  /// Backup for restoring the original pipe at destruction.
  int m_originalStreamTarget;

  // the pipe that deviate the original pipe
  Pipe m_deviationPipe;

  /// a buffer to store the flush() results
  std::string m_unprocessedData;
};

/**
 * @brief Class to handle external error capture.
 *        This class role is to capture and process external error messages, using the geos logger for
*         better tracing, logging and handling of messages.
 */
class ExternalErrorHandler
{
public:

  /**
   * @brief A functor executed for each error mesage to process, taking the message as the 1st
   *        string_view parameter, and the detectionLocation as the 2nd.
   * @see defaultErrorHandling() for a default implementation.
   */
  using ErrorHandlingFunctor = OutputStreamDeviation::LineHandlingFunctor;

  /**
   * @brief Strinct singleton pattern has been choosen since we will only have single sources of external
   *        errors (stderr for now, we could extend that for HYPRE errors, or for more dependencies).
   * @return The unique global instance.
   */
  static ExternalErrorHandler & instance();

  /**
   * @brief Destructor, disable all error piping features.
   */
  ~ExternalErrorHandler();

  /**
   * @brief Set the function that process the external errors that have been captured. The processing
   *        typically consists in using the given error message, adding metadata, and logging the message.
   * @param errorHandlingFunctor see ErrorHandlingFunctor.
   * @note Implementation treat each independant lines as an single error.
   */
  void setErrorHandling( ErrorHandlingFunctor && errorHandlingFunctor )
  { m_processErrorFunctor = errorHandlingFunctor; }

  /**
   * @brief Enable capture of errors piped from the std::cerr stream.
   *        Helpful to capture GLIBC errors, or other errors from dependencies not managed by GEOS itself.
   * @param enable Enable the feature if true, disable it otherwise.
   * @note Disabled by default.
   */
  void enableStderrPipe( bool enable );

  /**
   * @brief Process all awaiting captured errors that were produced externally, then clear the error stream.
   * @param detectionLocation A label to describe when the flush() operation is being made, thus
   *                          explaining to the user when the error has been detected.
   * @see setErrorHandling() to set the error processing procedure.
   */
  void flush( std::string_view detectionLocation );

  /**
   * @brief Not designed for direct calls, error handling function in default use if never calling
   *        setErrorHandling().
   * @param errorMsg the error text message.
   * @param detectionLocation A label to describe to the user when the error has been detected.
   */
  static void defaultErrorHandling( std::string_view errorMsg, std::string_view detectionLocation );

private:
  std::unique_ptr< OutputStreamDeviation > m_stderrDeviation;

  ErrorHandlingFunctor m_processErrorFunctor;

  ExternalErrorHandler();
};

} /* namespace geos */

#endif
