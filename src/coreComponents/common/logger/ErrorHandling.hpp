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

// Source includes
#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"

using namespace std;

namespace geos
{

/**
 * @class ErrorLogger
 * @brief Class to format and write the error/warning message that occured during the initialization
 */
class ErrorLogger
{

public:

  /**
   * @brief Construct a new Error Logger object
   *
   */
  ErrorLogger();

  /**
   * @brief Create the yaml file if the option is specified in the command line options
   * 
   */
  void createFile();

  /**
   * @enum MsgType
   * Enum listing the different types of possible errors
   */
  enum class MsgType
  {
    Error,
    Warning,
    Exception
  };

  /**
   * @brief Stores contextual information about the error that occurred and assigns it a priority (default is 0)
   * 
   */
  struct ContextInfo
  {
    std::map< std::string, std::string > m_ctxInfo;
    integer m_priority = 0;

    /**
     * @brief Set the priority of the current error context information
     * 
     * @param priority 
     * @return ContextInfo& 
     */
    ContextInfo & setPriority( integer priority )
    { m_priority = priority; return *this; }
  };

  /**
   * @brief Struct to define the error/warning message
   *
   */
  struct ErrorMsg
  {
    MsgType m_type;
    std::string m_msg;
    std::string m_file;
    integer m_line;
    std::vector< int > m_ranksInfo;
    std::vector< ContextInfo > m_contextsInfo;
    std::vector< std::string > m_sourceCallStack;

    int n = 0;

    /**
     * @brief Construct a new Error Msg object
     *
     */
    ErrorMsg() {};

    /**
     * @brief Construct a new Error Msg object
     *
     * @param msgType The type of the message (error or warning)
     * @param msgContent The error/warning message content
     * @param msgFile The file name where the error occcured
     * @param msgLine The line where the error occured
     */
    ErrorMsg( MsgType msgType, std::string msgContent, std::string msgFile, integer msgLine )
      : m_type( msgType ), m_msg( msgContent ), m_file( msgFile ), m_line( msgLine ) {}


    ErrorLogger * parent = nullptr;

    /**
     * @brief Fill the msg field of the structure with the error message
     *
     * @param e is the exception
     * @return ErrorMsg&
     */
    ErrorMsg & addToMsg( std::exception const & e );

    /**
     * @brief
     *
     * @param msg Add information about the error that occured to the msg field of the structure
     * @return ErrorMsg&
     */
    ErrorMsg & addToMsg( std::string const & msg );

    /**
     * @brief Set the Code Location object
     *
     * @param msgFile
     * @param msgLine
     * @return ErrorMsg&
     */
    ErrorMsg & setCodeLocation( string msgFile, integer msgLine );

    /**
     * @brief Set the Type object
     *
     * @param msgType
     * @return ErrorMsg&
     */
    ErrorMsg & setType( MsgType msgType );

    /**
     * @brief Adds one or more context elements to the error
     * 
     * @tparam Args 
     * @param args 
     */
    template< typename ... Args >
    void addContextInfo( Args && ... args );

    /**
     * @brief Set the rank on which the error is raised
     * 
     * @param rank 
     * @return ErrorMsg& 
     */
    ErrorMsg & setRank( int rank );

    /**
     * @brief Add stack trace information about the error/warning message to the ErrorMsg structure
     *
     * @param ossStackTrace stack trace information
     */
    ErrorLogger::ErrorMsg & addCallStackInfo( std::string const & ossStackTrace );

    private:
      /**
       * @brief Add contextual information about the error/warning message to the ErrorMsg structure
       *
       * @param info DataContext information  stored into a map
       */
      void addContextInfoImpl( ContextInfo && ctxInfo );
  };

  /**
   * @brief Return the error message information at the step where this getter is called
   * @return The current error msg
   */
  ErrorMsg & currentErrorMsg()
  { return m_currentErrorMsg; }

  /**
   * @brief Convert a MsgType into a string
   *
   * @param type
   * @return std::string
   */
  std::string toString( MsgType type );

  /**
   * @brief Write the error message in the yaml file regarding indentation and line break
   *
   * @param msg
   */
  void streamMultilineYamlAttribute( std::string_view msg, std::ofstream & yamlFile );

  /**
   * @brief Add the error/warning message into the yaml file
   *
   * @param errorMsg The error message informations formatted by the associated structure
   */
  void write( ErrorMsg const & errorMsg );

  /**
   * @brief Returns true whether the yaml file writing option is enabled by the user otherwise false 
   * 
   * @return true 
   * @return false 
   */
  bool writeFile() const 
  { return m_writeYaml; }

  /**
   * @brief Set the Write Value object
   * True whether the yaml file writing option is enabled by the user otherwise false
   * @param value 
   */
  void setWriteValue( bool value )
  { m_writeYaml = value; }

  /**
   * @brief Set the name of the yaml file if specified by user (default is "errors.yaml")
   * 
   * @param filename 
   */
  void setFilename( std::string filename )
  { m_filename = filename; }

private:
  // The error constructed via exceptions
  ErrorMsg m_currentErrorMsg;
  // Write in the yaml file
  bool m_writeYaml = false;
  // Yaml file name 
  std::string m_filename = "errors.yaml";
};

extern ErrorLogger errorLogger;

template< typename ... Args >
void ErrorLogger::ErrorMsg::addContextInfo( Args && ... args )
{
  ( this->addContextInfoImpl( ContextInfo( args ) ), ... );
}

} /* namespace geos */

#endif
