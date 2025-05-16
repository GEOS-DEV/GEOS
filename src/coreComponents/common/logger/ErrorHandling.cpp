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
#include "common/logger/Logger.hpp"

// System includes
#include <fstream>
#include <iostream>
#include <utility>
#include <sstream>

namespace geos
{
static constexpr std::string_view m_filename = "errors.yaml";

ErrorLogger errorLogger{};

ErrorLogger::ErrorLogger()
{
  m_currentErrorMsg.parent = this;
  std::ofstream yamlFile( std::string( m_filename ), std::ios::out );
  if( yamlFile.is_open() )
  {
    yamlFile << "errors: \n";
    yamlFile.close();
  }
  else
  {
    GEOS_LOG( GEOS_FMT( "Unable to open error file for writing: {}", m_filename ) );
  }
}

void ErrorLogger::ErrorMsg::addContextInfo( std::map< std::string, std::string > && info )
{
  m_contextsInfo.emplace_back( std::move( info ) );
}

void ErrorLogger::ErrorMsg::addRankInfo( int rank )
{
  m_ranksInfo.push_back( rank );
}

void ErrorLogger::ErrorMsg::addCallStackInfo( std::string const & ossStackTrace )
{
  std::istringstream iss( ossStackTrace );
  std::string stackLine;
  std::size_t index; 

  while( std::getline( iss, stackLine) )
  {
    index = stackLine.find(':');
    m_sourceCallStack.push_back( stackLine.substr( index + 1 ) );
  }
}

std::string ErrorLogger::toString( ErrorLogger::MsgType type )
{
  switch( type )
  {
    case ErrorLogger::MsgType::Error: return "Error";
    case ErrorLogger::MsgType::Warning: return "Warning";
    default: return "Unknown";
  }
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addToMsg( std::exception const & e )
{
  parent->m_currentErrorMsg.m_msg = e.what(); 
  return parent->m_currentErrorMsg;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addToMsg( std::string const&  errorMsg )
{
  parent->m_currentErrorMsg.m_msg = GEOS_FMT( "{:>6}{}", " ", errorMsg ) + parent->m_currentErrorMsg.m_msg; // Inverser l'ordre FILO
  return parent->m_currentErrorMsg;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::setCodeLocation( string msgFile, integer msgLine )
{
  parent->m_currentErrorMsg.m_file = msgFile;
  parent->m_currentErrorMsg.m_line = msgLine;
  return parent->m_currentErrorMsg;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::setType( ErrorLogger::MsgType msgType )
{
  parent->m_currentErrorMsg.m_type = msgType;
  return parent->m_currentErrorMsg;
}
 
void ErrorLogger::write( ErrorLogger::ErrorMsg const & errorMsg ) //const
{
  std::ofstream yamlFile( std::string( m_filename ), std::ios::app );
  if( yamlFile.is_open() )
  {
    yamlFile << GEOS_FMT( "{:>2}- type: {}\n", " ", errorLogger.toString( errorMsg.m_type ) );
    yamlFile << GEOS_FMT( "{:>4}rank: ", " " );
    for( size_t i = 0; i < errorMsg.m_ranksInfo.size(); i++ )
    {
      yamlFile << errorMsg.m_ranksInfo[i];
    }
    yamlFile << "\n";
    yamlFile << GEOS_FMT( "{:>4}message: >-\n{} \n", " ", errorMsg.m_msg );
    if( !errorMsg.m_contextsInfo.empty() )
    {
      yamlFile << GEOS_FMT( "{:>4}contexts:\n", " " );

      for( size_t i = 0; i < errorMsg.m_contextsInfo.size(); i++ )
      {
        for( auto const & [key, value] : errorMsg.m_contextsInfo[i] )
        {
          if( key == "inputFileLine" )
          {
            yamlFile << GEOS_FMT( "{:>8}{}: {}\n", " ", key, value );
          }
          else 
          {
            yamlFile << GEOS_FMT( "{:>6}- {}: {}\n", " ", key, value );
          }
        }
      }
    }

    yamlFile << GEOS_FMT( "{:>4}sourceLocation:\n", " " );
    yamlFile << GEOS_FMT( "{:>6}file: {}\n", " ", errorMsg.m_file );
    yamlFile << GEOS_FMT( "{:>6}line: {}\n", " ", errorMsg.m_line );
    
    yamlFile << GEOS_FMT( "{:>4}sourceCallStack:\n", " " );
    
    for( size_t i = 0; i < errorMsg.m_sourceCallStack.size(); i++ )
    {
      if( i < 2 || i == errorMsg.m_sourceCallStack.size() - 1 ) continue; 
      yamlFile << GEOS_FMT( "{:>6}- {}: {}\n", " ", i-2, errorMsg.m_sourceCallStack[i] );
    }

    yamlFile.flush();
    GEOS_LOG( GEOS_FMT( "The error file {} was created successfully.", m_filename ) );
  }
  else
  {
    GEOS_LOG( GEOS_FMT( "Unable to open error file for writing: {}", m_filename ) );
  }
}

} /* namespace geos */
