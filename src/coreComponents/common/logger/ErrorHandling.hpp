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
#include "common/format/Format.hpp"
#include "common/format/StringUtilities.hpp"
#include <mutex>

namespace geos
{

/**
 * @struct ErrorContext
 * Store contextual information about the error that occurred and assign it a priority
 * default is 0
 */
struct ErrorContext
{

  /**
   * @enum Attribute
   * Enumeration used to secure potential map keys
   */
  enum class Attribute
  {
    InputFile,
    InputLine,
    DataPath,
    DetectionLoc,
    Signal,
  };

  /// String containing the target object name followed by the the file and line declaring it.
  string m_formattedContext;

  /**
   * @brief The map contains contextual information about the error
   * It could be something like
   * "file" = "/path/to/file.xml"
   * "line" = "24"
   * or something like
   * "dataPath" = "/Functions/co2brine_philipsDensityTable
   * The key is a field of the Attribute enumeration and is converted to a string for writing in the YAML
   */
  map< Attribute, std::string > m_attributes;

  /**
   * @brief Priority level assigned to an error context.
   * @details Used to prioritize contexts (higher values = more relevant). Default is 0.
   */
  integer m_priority = 0;

  /**
   * @brief Construct to initialize ErrorContext
   * @param formattedContext String containing the target object name followed by the the file and line declaring it.
   * @param attributes Map containing contextual information about the error
   */
  ErrorContext( string formattedContext, map< Attribute, std::string > attributes ):
    m_formattedContext( formattedContext ),
    m_attributes( attributes ) {};

  /**
   * @brief Construct to initialize ErrorContext given a string containing the context and his priority
   * @param formattedContext String containing the target object name followed by the the file and line declaring it.
   * @param attributes Map containing contextual information about the error
   * @param priority Priority level assigned to an error context.
   */
  ErrorContext( string formattedContext, map< Attribute, std::string > attributes, integer priority ):
    m_formattedContext( formattedContext ),
    m_attributes( attributes ),
    m_priority( priority ) {};

  /**
   * @brief Set the priority value of the current error context information.
   * This way the different context information will appear in descending order during the error log output.
   * @param priority the new value to asign
   * @return the reference to the corresponding error
   */
  ErrorContext & setPriority( integer priority )
  { m_priority = priority; return *this; }

  /**
   * @brief Convert a value from the Attribute enumeration to a string
   * @param attribute the value of the enumeration to be converted
   * @return a string representation of the enumeration value
   */
  static std::string attributeToString( Attribute attribute );

};

/**
 * @enum MsgType
 * Enum listing the different types of possible diagnostics
 */
enum class MsgType
{
  Error,
  ExternalError,
  Warning,
  Exception,
  Undefined
};

/**
 * @brief Struct to construct the diagnostic message object
 */
struct DiagnosticMsg
{
  /// Type of diagnostic (Warning, Error or Exception)
  MsgType m_type = MsgType::Undefined;
  /// the message that can be completed
  std::string m_msg;
  /// the cause of the error (erroneous condition, failed assertion...) if identified (optional)
  std::string m_cause;
  /// the rank(s) on which the diagnostic occured
  std::set< int > m_ranksInfo;
  /// the source location file
  std::string m_file;
  /// the source location line (default is 0)
  integer m_line = 0;
  /// Additional information about the diagnostic in the input file
  stdVector< ErrorContext > m_contextsInfo;
  /// the stack trace
  stdVector< std::string > m_sourceCallStack;
  /// Indicates whether the stored call stack trace is valid and usable.
  bool m_isValidStackTrace = false;
};

/**
 * @brief Builder class for constructing DiagnosticMsg  objects
 */
class DiagnosticMsgBuilder
{
public:

/**
 * @brief Initialize a new DiagnosticMsg
 * @param msg The DiagnosticMsg being built
 * @param msgType Type of the diagnostic
 * @param msgContent The message of the diagnostic. It can be completed afterward
 * @param rank The rank on which the diagnostic occured
 * @return DiagnosticMsgBuilder
 */
  static DiagnosticMsgBuilder init( DiagnosticMsg & msg,
                                    MsgType msgType,
                                    std::string_view msgContent,
                                    integer rank );

  /**
   * @brief Modify an existing DiagnosticMsg
   * @param errorMsg The existing DiagnosticMsg
   * @return DiagnosticMsgBuilder
   */
  static DiagnosticMsgBuilder modify( DiagnosticMsg & errorMsg );

  /**
   * @brief Append exception text to the message
   * @param e The exception containing text to add
   * @param toEnd If true, append at end; otherwise prepend
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & addToMsg( std::exception const & e, bool toEnd  = false );

  /**
   * @brief Append text to the message
   * @param msg The text to add
   * @param toEnd If true, append at end; otherwise prepend
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & addToMsg( std::string_view msg, bool toEnd = false );

  /**
   * @brief Adds one or more context elements to the error
   * @tparam Args Variadic pack of compatible types (ErrorContext / DataContext)
   * @param args List of context data structures.
   * @return Reference to the current instance for method chaining.
   */
  template< typename ... Args >
  DiagnosticMsgBuilder & addContextInfo( Args && ... args )
  {
    ( this->addContextInfoImpl( ErrorContext( args ) ), ... );
    return *this;
  }

  /**
   * @brief Add where the detection occured
   * @param detectionLocation The context where the diagnostic happoned
   * @return The instance, for builder pattern.
   */
  DiagnosticMsgBuilder & addDetectionLocation( string_view detectionLocation );

  /**
   * @brief Add the signal to the DiagnosticMsg.
   *        - the signal can be one of the main error signals.
   *        - if the signal is SIGFPE, the nature of floating point error will be interpreted.
   * @param sig The signal, from ISO C99 or POSIX standard.
   * @param toEnd adds the message to the end if true, at the start otherwise.
   * @return The instance, for builder pattern.
   */
  DiagnosticMsgBuilder & addSignal( integer sig, bool toEnd = false );

  /**
   * @brief Set the source code location values (file and line where the error is detected)
   * @param msgFile Name of the source file location to add
   * @param msgLine Line of the source file location to add
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & setCodeLocation( std::string_view msgFile, integer msgLine );

  /**
   * @brief Set the type of the error, (amoung one of the MsgType)
   * @param msgType The type can be error, warning or exception
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & setType( MsgType msgType );

  /**
   * @brief Set the cause of the error
   * @param cause See documentation of m_cause.
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & setCause( std::string_view cause );

  /**
   * @brief Add a rank on which the error has been raised
   * @param rank The rank value
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & addRank( integer rank );

  /**
   * @brief Add the stack trace information about the error
   * @param stacktrace stack trace information to add
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & addCallStackInfo( std::string_view stacktrace );

  /**
   * @return Get the DiagnosticMsg
   */
  DiagnosticMsg & getDiagnosticMsg();

private:

  /**
   * @brief Private constructor - use init() or modify() instead
   * @param msg Reference to the DiagnosticMsg to build/modify
   */
  DiagnosticMsgBuilder( DiagnosticMsg & msg ):
    m_errorMsg( msg ){}

  /**
   * @brief Add contextual information about the error/warning
   * @param ctxInfo rvalue of the ErrorContext class
   */
  DiagnosticMsgBuilder & addContextInfoImpl( ErrorContext && ctxInfo );

  /// The diagnosticMsg being constructed
  DiagnosticMsg & m_errorMsg;
};

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
   * outputs (stream specified, std::cout by default + optional yaml file)
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

private:

  /// The error constructed via exceptions
  DiagnosticMsg m_getCurrentExceptionMsg;
  /// Indicate whether the write to YAML command line option is enabled
  bool m_writeYaml = false;
  /// YAML file name
  std::string_view m_filename = "errors.yaml";
  /// The stream used for the log output. By default used std::cout
  std::ostream & m_stream = std::cout;
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

class ErrorHandler
{
public:

  static ErrorHandler & instance();

  void setProgramAborter( std::function< void() > const & abortingFunctor )
  { m_abortingFunctor = abortingFunctor; }

  /**
   * @brief Post error-handling function that terminates the program.
   */
  void abortProgram()
  { m_abortingFunctor(); }

private:

  std::function< void() > m_abortingFunctor;

  /**
   * @brief Static class, no public constructor
   */
  ErrorHandler();

};

} /* namespace geos */

#endif
