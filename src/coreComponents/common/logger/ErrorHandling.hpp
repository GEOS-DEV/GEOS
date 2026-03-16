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

#ifndef INITIALIZATION_ERROR_LOGGER_HPP
#define INITIALIZATION_ERROR_LOGGER_HPP

#include "common/DataTypes.hpp"
#include "common/format/LogPart.hpp"
#include "common/logger/LogHistory.hpp"
#include <mutex>

namespace geos
{

/**
 * @brief Logger for formatting and outputting diagnostics
 */
class ErrorLogger
{

public:

  /**
   * @return Global instance of the ErrorLogger class used for error/warning reporting.
   * @details This global instance is used across the codebase to log errors, warnings, and exceptions,
   *          and to write structured output of errors. It is used through the logging macros.
   * @note - local instances are possible for more specialized logging.
   *       - currently not available on GPU, use GEOS_WARNING/ERROR/ASSERT macros for this usecase.
   */
  GEOS_HOST static ErrorLogger & global();

  /**
   * @brief Create the YAML file or overwrite the contents if a YAML file of the same name already exists
   * And write its header when the command line option is enabled
   */
  void createFile();

  /**
   * @brief Enable the YAML file output, which is false by default
   * @param value A value of true enable the file writing
   */
  void enableFileOutput( bool value )
  { m_writeYaml = value; }

  /**
   * @return true if the YAML file output is enabled
   */
  bool isOutputFileEnabled() const
  { return m_writeYaml; }

  /**
   * @brief Set the name of the YAML file if specified by user
   * default is "errors.yaml"
   * @param filename the name of the YAML file
   */
  void setOutputFilename( std::string_view filename )
  { m_filename = filename; }

  /**
   * @return The file name of the output error file
   */
  std::string_view getOutputFilename()
  { return m_filename; }

  /**
   * @brief Convert a MsgType into a string
   * @param type the message type label
   * @return the string representation of the message type
   */
  static std::string toString( MsgType type );

  /**
   * @return Return the const general log stream
   */
  std::ostream const & getErrorStream() const
  { return m_stream; }

  /**
   * @return Return the const DiagnosticMsg
   */
  DiagnosticMsg const & getCurrentExceptionMsg() const
  { return m_getCurrentExceptionMsg;}

  /**
   * @brief Start building a new exception message
   * @param msgType Type of diagnostic (Warning, Error or Exception)
   * @param msgContent the message that can be completed
   * @param rank the rank(s) on which the diagnostic occured
   * @return Builder for the exception
   * @note One exception can exist at a time
   */
  DiagnosticMsgBuilder initCurrentExceptionMessage( MsgType msgType,
                                                    std::string_view msgContent,
                                                    integer rank );

  /**
   * @brief Modify/Continue building the current exception message
   * @return Builder for the exception
   */
  DiagnosticMsgBuilder modifyCurrentExceptionMessage()
  { return DiagnosticMsgBuilder::modify( m_getCurrentExceptionMsg ); }

  /**
   * @brief Write all the information retrieved about the current exception message into the instance
   * outputs (stream specified + optional yaml file)
   */
  void flushCurrentExceptionMessage();

  /**
   * @brief Write all the information retrieved about the diagnostic message into the instance
   * outputs (stream specified + optional yaml file)
   * @param errMsg a reference to the ErrorMsg to output, and will be re-initialized
   * @note Used for warnings and non-exception errors
   */
  void flushErrorMsg( DiagnosticMsg & errMsg );

  /**
   * @brief Format all information in ErrorMsg and write it to the specified output stream
   * @param errMsg The struct containing the error/warning object
   * @param os The output stream
   */
  static void formatMsgForLog( DiagnosticMsg const & errMsg, std::ostream & os );

  /**
   * @brief Write the ErrorMsg into the log stream output stream
   * @param errMsg The struct containing the error/warning object
   */
  void writeToLogStream( DiagnosticMsg & errMsg );

  /**
   * @brief Gets the current logger report data.
   * @return The current log part as a string.
   */
  LogHistory const & getLoggerReportData() const
  {return loggerMsgReportData;}

  /**
   * @brief Gets the current logger report data.
   * @return The current log part as a string.
   */
  LogHistory & getLoggerReportData()
  {return loggerMsgReportData;}

  /**
   * @brief Gets the current log part.
   * @return The current log part as a string.
   */
  string_view getCurrentLogPart() const
  {return m_currentLogPart;}

/**
 * @brief Sets the current log part.
 * @param currentLogPart The new log part to set.
 */
  void setCurrentLogPart( string_view currentLogPart )
  { m_currentLogPart = currentLogPart; }

private:

  /// The error constructed via exceptions
  DiagnosticMsg m_getCurrentExceptionMsg;
  /// The log history associated
  LogHistory loggerMsgReportData = {};

  /// Indicate whether the write to YAML command line option is enabled
  bool m_writeYaml = false;
  /// YAML file name
  std::string_view m_filename = "errors.yaml";
  /// The stream used for the log output. By default used std::cout
  std::ostream & m_stream = std::cout;
  /// The current log part being executed
  string m_currentLogPart;
  /// Avoid concurrent access between threads for log outputs
  std::mutex m_errorHandlerAsciiMutex;
  /// Avoid concurrent access between threads for yaml outputs
  std::mutex m_errorHandlerYamlMutex;

  /**
   * @brief Write all the information retrieved about the error/warning message into the YAML stream
   * @param errorMsg a reference to the diagnostic message
   */
  void writeToYamlStream( DiagnosticMsg & errMsg );

  /**
   * @brief Write the error message in the YAML file regarding indentation and line break
   * @param msg the message to write in the YAML
   * @param yamlFile The yaml file steam
   * @param indent The desired indentation
   */
  void streamMultilineYamlAttribute( std::string_view msg, std::ofstream & yamlFile,
                                     std::string_view indent );
};



} /* namespace geos */

#endif
