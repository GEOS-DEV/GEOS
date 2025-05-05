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
#include <iostream>

namespace geos
{
  ErrorLogger errorLogger;

  std::string ErrorLogger::toString( ErrorLogger::MsgType type )
  {
    switch ( type )
    {
      case ErrorLogger::MsgType::Error: return "Error";
      case ErrorLogger::MsgType::Warning: return "Warning";
      default: return "Unknown";
    }
  }

  void ErrorLogger::write( ErrorLogger::ErrorMsg const & errorMsg )
  {
    std::string filename = "errors.yaml";
    std::ifstream checkYamlFile( filename );
    bool isEmpty = checkYamlFile.peek() == std::ifstream::traits_type::eof();
    checkYamlFile.close();

    std::ofstream yamlFile( filename, std::ios::app ); 
    if( yamlFile.is_open() )
    {
      if( isEmpty )
      {
        yamlFile << "errors: \n";
      }
      yamlFile << GEOS_FMT( "{:>2}- type: {}\n", " ", errorMsg.msg );
      yamlFile << GEOS_FMT( "{:>4}message: {}\n", " ", errorLogger.toString( errorMsg.type ) );
      yamlFile << GEOS_FMT( "{:>4}inputFileLocation:\n", " " );
      yamlFile << GEOS_FMT( "{:>6}- file: {}\n", " ", errorMsg.file );
      yamlFile << GEOS_FMT( "{:>8}line: {}\n\n", " ", errorMsg.line );
      yamlFile << GEOS_FMT( "{:>4}sourceLocation:\n", " " );
      yamlFile << GEOS_FMT( "{:>6}file: {}\n", " ", errorMsg.inputFileName );
      yamlFile << GEOS_FMT( "{:>6}line: {}\n\n", " ", errorMsg.inputFileLine );
      yamlFile.close();
      std::cout << "YAML file created successfully.\n";
    } 
    else 
    {
      std::cerr << "Unable to open file.\n";
    }
  }

} /* namespace geos */