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
 * @file DiagnosticMessage.hpp
 */

#include "common/format/EnumStrings.hpp"

namespace geos
{

/**
 * @enum MsgType
 * Enum listing the different types of possible errors
 */
enum class MsgType
{
  Error,
  ExternalError,
  Warning,
  Exception,
  Count // internal use, keep at last
};

/// Declare strings associated with output MsgType values.
ENUM_STRINGS( MsgType,
              "Error",
              "ExternalError",
              "Warning",
              "Exception",
              "Undefined" );


/**
 * @struct DiagnosticContext
 * Store contextual information about the error that occurred and assign it a priority
 * default is 0
 */
struct DiagnosticContext
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
   * @brief Construct to initialize DiagnosticContext
   * @param formattedContext String containing the target object name followed by the the file and line declaring it.
   * @param attributes Map containing contextual information about the error
   */
  DiagnosticContext( string formattedContext, map< Attribute, std::string > attributes ):
    m_formattedContext( formattedContext ),
    m_attributes( attributes ) {};

  /**
   * @brief Construct to initialize DiagnosticContext given a string containing the context and his priority
   * @param formattedContext String containing the target object name followed by the the file and line declaring it.
   * @param attributes Map containing contextual information about the error
   * @param priority Priority level assigned to an error context.
   */
  DiagnosticContext( string formattedContext, map< Attribute, std::string > attributes, integer priority ):
    m_formattedContext( formattedContext ),
    m_attributes( attributes ),
    m_priority( priority ) {};

  /**
   * @brief Set the priority value of the current error context information
   * This way the different context information will appear in descending order during the error log output
   * @param priority the new value to asign
   * @return the reference to the corresponding error
   */
  DiagnosticContext & setPriority( integer priority )
  { m_priority = priority; return *this; }

  /**
   * @brief Convert a value from the Attribute enumeration to a string
   * @param attribute the value of the enumeration to be converted
   * @return a string representation of the enumeration value
   */
  static std::string attributeToString( Attribute attribute );
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
  /// the path of the source location file
  std::string m_file;
  /// the source location line (default is 0)
  integer m_line = 0;
  /// The log part where the diagnostic occured
  string m_logPart;
  /// Additional information about the diagnostic in the input file
  std::vector< DiagnosticContext > m_contextsInfo;
  /// the stack trace
  std::vector< std::string > m_sourceCallStack;
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
   * @tparam Args Variadic pack of compatible types (DiagnosticContext / DataContext)
   * @param args List of context data structures.
   * @return Reference to the current instance for method chaining.
   */
  template< typename ... Args >
  DiagnosticMsgBuilder & addContextInfo( Args && ... args )
  {
    ( this->addContextInfoImpl( DiagnosticContext( args ) ), ... );
    return *this;
  }

  /**
   * @brief Add the detection location the DiagnosticMsg
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
   * @brief Set the type of the error
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
   * @param rank The value to add
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & addRank( integer rank );

  /**
   * @brief Add stack trace information about the error
   * @param stacktrace stack trace information to add
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & addCallStackInfo( std::string_view stacktrace );

  /**
   * @brief Set log part where the disgnostic occured
   * @param logPart The targetted log part
   * @return Reference to the current instance for method chaining.
   */
  DiagnosticMsgBuilder & setLogPart( std::string_view logPart );



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
   * @param ctxInfo rvalue of the DiagnosticContext class
   */
  DiagnosticMsgBuilder & addContextInfoImpl( DiagnosticContext && ctxInfo );

  /// The diagnosticMsg being constructed
  DiagnosticMsg & m_errorMsg;
};
}
