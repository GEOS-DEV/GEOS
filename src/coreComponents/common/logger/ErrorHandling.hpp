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
  enum class TypeMsg 
  {
      ERROR,
      WARNING
  };

  /**
   * @brief Struct to define the error/warning message
   * 
   */
  struct ErrorMsg 
  {
    TypeMsg type; 
    std::string msg; 
    std::string file; 
    integer line;
  };

  /**
   * @brief Structured the message based on the provided parameters
   * 
   * @param type The type of the message (error or warning)
   * @param msg The error/warning message content
   * @param file The file name where the error occured
   * @param line The line where the error occured 
   * @return ErrorMsg 
   */
  ErrorMsg errorMsgformatter( TypeMsg const & type, 
                              std::string const & msg,
                              std::string const & file, 
                              integer line )
  {
    return { type, msg, file, line };
  };

  std::string toString( TypeMsg type );

  /**
   * @brief Add the error/warning message into the yaml file
   * 
   * @param errorMsg The error message informations formatted by the associated structure
   */
  void errorMsgWritter( ErrorMsg const & errorMsg );
};

} /* namespace geos */

# endif 