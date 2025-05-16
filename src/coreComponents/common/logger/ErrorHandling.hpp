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
   * @enum MsgType
   * Enum listing the different types of possible errors
   */
  enum class MsgType
  {
    Error,
    Warning
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
    std::vector< std::map< std::string, std::string > > m_contextsInfo;
    std::vector< std::string > m_sourceCallStack;

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
     * @param msgFile The file name where the error occured
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
     * @brief Add contextual information about the error/warning message to the ErrorMsg structure
     *
     * @param info DataContext information  stored into a map
     */
    void addContextInfo( std::map< std::string, std::string > && info );

    void addRankInfo( int rank );
    
    /**
     * @brief Add stack trace information about the error/warning message to the ErrorMsg structure
     *
     * @param ossStackTrace stack trace information
     */
    void addCallStackInfo( std::string const & ossStackTrace );
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
   * @brief Add the error/warning message into the yaml file
   *
   * @param errorMsg The error message informations formatted by the associated structure
   */
  void write( ErrorMsg const & errorMsg );

private:
  // The error constructed via exceptions
  ErrorMsg m_currentErrorMsg; 
};

extern ErrorLogger errorLogger;

} /* namespace geos */

#endif