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
   * @enum TypeMsg 
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
    MsgType type; 
    std::string msg; 
    std::string file; 
    integer line;
    std::vector< std::map< std::string, std::string > > contextsInfo;
    // std::vector< std::string > sourceCallStack;

    /**
     * @brief Construct a new Error Msg object
     * 
     * @param t The type of the message (error or warning)
     * @param m The error/warning message content
     * @param f The file name where the error occured
     * @param l The line where the error occured 
     */
    ErrorMsg( MsgType t, std::string m, std::string f, integer l ) : type( t ), msg( m ), file( f ), line( l ) {}
    
    /**
     * @brief Add contextual information about the error/warning message to the ErrorMsg structure
     * 
     * @param info 
     */
    void addContextInfo( std::map< std::string, std::string > && info );
  };

  std::string toString( MsgType type );

  /**
   * @brief Add the error/warning message into the yaml file
   * 
   * @param errorMsg The error message informations formatted by the associated structure
   */
  void write( ErrorMsg const & errorMsg );
};

extern ErrorLogger errorLogger;

} /* namespace geos */

# endif 