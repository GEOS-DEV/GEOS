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
 * @brief This file provides the infrastructure to capture external errors.
 * @note Below is the architecture of the external error managment, in the scenario of a problematic
 *       infrastructure which deviates (thus breaks) stderr.
 * ```
 *  ________________________________________________________________________
 * |                            GEOS APPLICATION                            |
 * |------------------------------------------------------------------------|
 * |   _____________________         ____________________                   |
 * |  |  GEOS DEPENDANCIES  |       |   GEOS CORE        |                  |
 * |  |---------------------|       |--------------------|                  |
 * |  |   _________  _____  |       |   __________       |                  |
 * |  |  |         ||     | |       |  |          |      |                  |
 * |  |  | LvArray || std | |       |  | loggings |------------             |                     ________
 * |  |  |_________||_____| |       |  |__________|      |    |             |  log pipe          |        |
 * |  | _______   |    |    |       |   ______________   |    +--------------------------------->| stdout |
 * |  ||       |  |    |    |       |  |              |  |    |             |                    |________|
 * |  || Hypre |--+----+--------    |  | ErrorHandler |-------+             |
 * |  ||_______|            |  |    |  |______________|  |    |             |
 * |  |_____________________|  |    |____________________|    |             |
 * |                           |                              |             |
 * |                  external |                   ___________|___________  |
 * |                    errors |                  | OutputStreamDeviation | |      deviated       ________
 * |                    stream |                  |           |           | |     error pipe     |        |
 * |                           -------------------------------+   x-----------------x  +-------->| stderr |
 * |                                              |_______________________| |          |         |________|
 * |________________________________________________________________________|          |
 *                                                                               .............
 *                                                                               :   HPC X   :    potencial infratructure
 *                                                                               :   system  : <- which deviates the stderr
 *                                                                               : messaging :    for other reasons.
 *                                                                               :...........:
 * ```
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

  /// Posix identifier, can be a file handle, an error number... Must be consistent with posix functions.
  using PosixId = int;

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
  OutputStreamDeviation( PosixId fileNo );

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
    /// the file descriptors of the stream, read end first, then write end.
    PosixId fileDescriptorsArray[2];

    /**
     * @brief Initialize the file descriptors to open the pipe as a "virtual file".
     */
    void create();

    /**
     * @return the read end file descriptor of the pipe (serve to read data).
     */
    PosixId & readEnd()
    { return fileDescriptorsArray[0]; }

    /**
     * @return the write end file descriptor of the pipe (serve to write data).
     */
    PosixId & writeEnd()
    { return fileDescriptorsArray[1]; }

    /**
     * @brief Prevent or enable file descriptors of the pipe instance from being inherited by child processes.
     * @param inherit If true, enable file descriptors inheritance. If false, disable it.
     * @return True if successful, false otherwise. If the function fails, a warning is logged.
     * @details When disabled, this function sets the close-on-exec flag on a file descriptor,
     *          preventing it from being inherited by child processes when exec() is called.
     *          This helps prevent resource leaks and security issues in process management.
     */
    bool setDescriptorInheritanceMode( bool inherit );

    /**
     * @brief Redirect the write end of this pipe to a target stream, enabling capture of output
              written to that stream through the read end.
              The write end can still be used to write through the same pipe.
     * @param targetFd The file descriptor to redirect to (e.g., stdout, stderr).
     * @return True if successful, false otherwise.
     */
    bool redirectWriteEnd( int targetFd );

    /**
     * @brief Close the pipe ends.
     */
    void closePipe();
  };

  /// a special value to represant a disabled / not existing pipe.
  static constexpr PosixId m_disabledPipeEnd = -1;

  /// error values from POSIX functions.
  static constexpr PosixId m_errorResult = -1;

  /// the original pipe to deviate
  PosixId m_redirectedStream;

  /// Backup for restoring the original pipe at destruction.
  PosixId m_originalStreamTarget;

  /// the pipe that deviate the original pipe
  Pipe m_deviationPipe;

  /// a buffer to store the flush() results
  std::string m_unprocessedData;

  /**
   * @brief Set a file descriptor to non-blocking mode for asynchronous I/O operations.
   * @param pipeEnd The file descriptor to set the blocking mode.
   * @param nonBlocking If true, non-blocking mode is enabled. If false, blocking mode is used.
   * @return True if successful, false otherwise.
   * @details This function allows the read operation to return immediately even when no data is available,
   *          rather than blocking until data becomes available.
   */
  static bool setPipeEndBlockingMode( PosixId pipeEnd, bool nonBlocking );

  /**
   * @brief Duplicate a file descriptor, creating a new file descriptor that refers to the same underlying object.
   * @param pipeEnd The original file descriptor to duplicate.
   * @return The new duplicated file descriptor, or m_disabledPipeEnd on failure.
   * @details Used primarily for backup purposes to preserve the original stream state to restore that after cleanup.
   */
  static int duplicateDescriptor( int pipeEnd );

  /**
   * @brief Prepare the streaming buffer with a few optimisation by:
            - preallocating the intermediate message processing buffer,
            - growing the pipe stream to make it able to process the biggest error messages as one bit.
   */
  void prepareStreamingBuffer();

  /**
   * @brief Close the given pipe end
   * @param posixErrNo Reference to the POSIX pipe identifier. Set to m_disabledPipeEnd afterwards.
   */
  static void closePipeEnd( PosixId & pipeEnd );

  /**
   * @param posixErrNo POSIX error number
   * @return the error description of the given error number
   */
  static string getLastPosixErrorString();
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
  void enableStderrPipeDeviation( bool enable );

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

  /**
   * @param errorMsg The message generated
   * @return Return true if the message provided was generated by GEOS, false otherwise
   */
  bool static isNotAnErrorMsg( string_view errorMsg );

private:
  std::unique_ptr< OutputStreamDeviation > m_stderrDeviation;

  ErrorHandlingFunctor m_processErrorFunctor;

  ExternalErrorHandler();
};

} /* namespace geos */

#endif
