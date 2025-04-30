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
 * @file ErrorHandling.cpp
 */

// Source includes
#include "ErrorHandling.hpp"

// System includes
#include <fstream>
#include <filesystem>

namespace geos
{

  std::string ErrorLogger::toString( ErrorLogger::TypeMsg type )
  {
    switch ( type )
    {
      case ErrorLogger::TypeMsg::ERROR: return "Error";
      case ErrorLogger::TypeMsg::WARNING: return "Warning";
      default: return "Unknown";
    }
  }

  void ErrorLogger::errorMsgWritter( ErrorLogger::ErrorMsg const & errorMsg )
  {
    ErrorLogger logger;
    std::ofstream yamlFile( "errors.yaml", std::ios::app ); 
    if( yamlFile.is_open() )
    {
      std::ifstream yamlFileIn("errors.yaml");
      if( yamlFileIn.tellg() == 0 ) 
      {
        yamlFile << "errors: \n";
      }
      yamlFile << "     - message: " << errorMsg.msg << "\n";
      yamlFile << "       type: " << logger.toString( errorMsg.type ) << "\n";
      yamlFile << "       location: " << "\n";
      yamlFile << "             file: " << errorMsg.file << "\n";
      yamlFile << "             line: " << errorMsg.line << "\n\n";
      yamlFile.close();
      std::cout << "YAML file created successfully.\n";
    } 
    else 
    {
      std::cerr << "Unable to open file.\n";
    }
  }

} /* namespace geos */