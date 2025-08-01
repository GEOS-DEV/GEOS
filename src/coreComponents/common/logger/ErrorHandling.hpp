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

namespace geos
{

/**
 * @class ErrorLogger
 * @brief Class to format and write different error/warning information that occured during the initialization
 */
class ErrorLogger
{

public:
  /**
   * @enum MsgType
   * Enum listing the different types of possible errors
   */
  enum class MsgType
  {
    Error,
    Warning,
    Exception,
    Undefined
  };

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
      DataPath
    };

    /// The map contains contextual information about the error
    /// It could be something like
    /// "file" = "/path/to/file.xml"
    /// "line" = "24"
    /// or something like
    /// "dataPath" = "/Functions/co2brine_philipsDensityTable
    /// The key is a field of the Attribute enumeration and is converted to a string for writing in the YAML
    map< Attribute, std::string > m_attributes;
    integer m_priority = 0;

    /**
     * @brief Set the priority value of the current error context information
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
   * @brief Struct to construct the error/warning object
   */
  struct ErrorMsg
  {
    /// the error type (Warning, Error or Exception)
    MsgType m_type = ErrorLogger::MsgType::Undefined;
    /// the erreur message that can be completed
    std::string m_msg;
    /// the source location file corresponding to the error in the code
    std::string m_file;
    /// the source location line corresponding to the error in the code (default is 0)
    integer m_line = 0;
    /// the rank(s) on which the error occured
    std::vector< int > m_ranksInfo;
    /// Additional information about the error in the input file
    std::vector< ErrorContext > m_contextsInfo;
    /// the stack trace
    std::vector< std::string > m_sourceCallStack;

    /**
     * @brief Construct a default Error Message without field specification
     */
    ErrorMsg() {};

    /**
     * @brief Construct a new Error Message from parameters
     * @param msgType the type of the message (error or warning)
     * @param msgContent the error/warning message content
     * @param msgFile the source file name where the error occcured
     * @param msgLine the line where the error occured
     */
    ErrorMsg( MsgType msgType, std::string_view msgContent, std::string_view msgFile, integer msgLine )
      : m_type( msgType ), m_msg( msgContent ), m_file( msgFile ), m_line( msgLine ) {}

    /**
     * @brief Add text to the current error msg
     * @param e the exception containing text to add
     * @param toEnd indicates whether to add the message at the beginning (true) or at the end (false)
     *              default is false
     * @return the reference to the current instance
     */
    ErrorMsg & addToMsg( std::exception const & e, bool toEnd = false );

    /**
     * @brief Add text to the current error msg
     * @param msg the text to add
     * @param toEnd indicates whether to add the message at the beginning (true) or at the end (false)
     *              default is false
     * @return the reference to the current instance
     */
    ErrorMsg & addToMsg( std::string_view msg, bool toEnd = false );

    /**
     * @brief Set the source code location values (file and line where the error is detected)
     * @param msgFile name of the source file location to add
     * @param msgLine line of the source file location to add
     * @return the reference to the current instance
     */
    ErrorMsg & setCodeLocation( std::string_view msgFile, integer msgLine );

    /**
     * @brief Set the type of the error
     * @param msgType the type can be error, warning or exception
     * @return the reference to the current instance
     */
    ErrorMsg & setType( MsgType msgType );

    /**
     * @brief Set the rank on which the error is raised
     * @param rank the value to asign
     * @return the reference to the current instance
     */
    ErrorMsg & setRank( int rank );

    /**
     * @brief Add stack trace information about the error
     * @param ossStackTrace stack trace information to add
     * @return the reference to the current instance
     */
    ErrorMsg & addCallStackInfo( std::string_view ossStackTrace );

    /**
     * @return true if the YAML file output is enabled
     */
    bool isValidStackTrace() const
    { return m_isValidStackTrace; }

    /**
     * @brief Adds one or more context elements to the error
     * @tparam Args variadic pack of argument types
     * @param args list of DataContexts
     */
    template< typename ... Args >
    void addContextInfo( Args && ... args );

private:
    /**
     * @brief Add contextual information about the error/warning
     * @param ctxInfo rvalue of the ErrorContext class
     */
    void addContextInfoImpl( ErrorContext && ctxInfo );

    bool m_isValidStackTrace = false;
  };

  /**
   * @return true if the YAML file output is enabled
   */
  bool isOutputFileEnabled() const
  { return m_writeYaml; }

  /**
   * @brief Enable the YAML file output, which is false by default
   * @param value A value of true enable the file writing
   */
  void enableFileOutput( bool value )
  { m_writeYaml = value; }

  /**
   * @brief Set the name of the YAML file if specified by user
   * default is "errors.yaml"
   * @param filename the name of the YAML file
   */
  void setOutputFilename( std::string_view filename )
  { m_filename = filename; }

  /**
   * @brief Gives acces to the error message that is currently being constructed,
   *        potencially at various application layers
   *        Use flushErrorMsg() when the message is fully constructed and you want it to be output
   * @return the reference to the current instance
   */
  ErrorMsg & currentErrorMsg()
  { return m_currentErrorMsg; }

  /**
   * @brief Create the YAML file or overwrite the contents if a YAML file of the same name already exists
   * And write its header when the command line option is enabled
   */
  void createFile();

  /**
   * @brief Convert a MsgType into a string
   * @param type the message type label
   * @return the string representation of the message type
   */
  static std::string toString( MsgType type );

  /**
   * @brief Write all the information retrieved about the error/warning message into the YAML file
   *        and reset the errorMsg instance to its initial state
   * @param errorMsg a constant reference to the error
   */
  void flushErrorMsg( ErrorMsg & errorMsg );

private:
  /// The error constructed via exceptions
  ErrorMsg m_currentErrorMsg;
  /// Indicate whether the write to YAML command line option is enabled
  bool m_writeYaml = false;
  /// YAML file name
  std::string_view m_filename = "errors.yaml";

  /**
   * @brief Write the error message in the YAML file regarding indentation and line break
   * @param msg the message to write in the YAML
   *            For the exception type, this message can be added as needed
   */
  void streamMultilineYamlAttribute( std::string_view msg, std::ofstream & yamlFile,
                                     std::string_view indent );
};

extern ErrorLogger g_errorLogger;

template< typename ... Args >
void ErrorLogger::ErrorMsg::addContextInfo( Args && ... args )
{
  ( this->addContextInfoImpl( ErrorContext( args ) ), ... );
}

} /* namespace geos */

#endif
