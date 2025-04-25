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

#include "../DataTypes.hpp"

#include <yaml-cpp/yaml.h>
#include <fstream>

using namespace std;

class ErrorLogger 
{
  enum class TypeMsg 
  {
      ERROR,
      WARNING
  };

  struct ErrorMsg 
  {
    string msg; 
    TypeMsg type; 
    string file; 
    integer line;
  };

  // Fonction qui formatte 
  ErrorMsg errorMsgformatter( const string & type, 
                              const string & msg,
                              const string & file, 
                              integer line )
  {
    return { type, msg, file, line };
  };


  // Fonction qui écrit dans le yaml
  void errorMsgWritter( const ErrorMsg & errorMsg, const string filename )
  {
    YAML::Node newMsg;
    newMsg["message"] = errorMsg.msg;
    newMsg["type"] = errorMsg.type;

    YAML::Node location; 
    location["file"] = errorMsg.file; 
    location["line"] = errorMsg.line; 

    newMsg["location"] = location; 

    ofstream fout( filename, ios::app );
    fout << newMsg << "\n";
  };
};

# endif 