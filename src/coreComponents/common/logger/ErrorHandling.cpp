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

#include "ErrorHandling.hpp"
#include "common/logger/Logger.hpp"

#include <fstream>
#include <string_view>

namespace geos
{
static constexpr std::string_view g_level1Start = "  - ";
static constexpr std::string_view g_level1Next =  "    ";
static constexpr std::string_view g_level2Start = "    - ";
static constexpr std::string_view g_level2Next =  "      ";
static constexpr std::string_view g_level3Start = "      - ";
static constexpr std::string_view g_level3Next =  "        ";

ErrorLogger g_errorLogger{};

ErrorLogger::ErrorLogger()
{
  m_currentErrorMsg.parent = this;
}

void ErrorLogger::createFile()
{
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

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addToMsg( std::exception const & e )
{
  parent->m_currentErrorMsg.m_msg = e.what();
  return parent->m_currentErrorMsg;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addToMsg( std::string errorMsg )
{
  parent->m_currentErrorMsg.m_msg = errorMsg + parent->m_currentErrorMsg.m_msg; 
  return parent->m_currentErrorMsg;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::setCodeLocation( std::string_view msgFile, integer msgLine )
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

void ErrorLogger::ErrorMsg::addContextInfoImpl( ErrorLogger::ContextInfo && ctxInfo )
{
  m_contextsInfo.emplace_back( std::move( ctxInfo ) );
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::setRank( int rank )
{
  m_ranksInfo.push_back( rank );
  return *this;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addCallStackInfo( std::string ossStackTrace )
{
  std::istringstream iss( ossStackTrace );
  std::string stackLine;
  std::size_t index;

  while( std::getline( iss, stackLine ) )
  {
    index = stackLine.find( ':' );
    m_sourceCallStack.push_back( stackLine.substr( index + 1 ) );
  }

  return *this;
}

std::string ErrorLogger::toString( ErrorLogger::MsgType type )
{
  switch( type )
  {
    case ErrorLogger::MsgType::Error: return "Error";
    case ErrorLogger::MsgType::Warning: return "Warning";
    case ErrorLogger::MsgType::Exception: return "Exception";
    default: return "Unknown";
  }
}

void ErrorLogger::streamMultilineYamlAttribute( std::string_view msg, std::ofstream& yamlFile )
{
  while( !msg.empty() )
  {
    const size_t index = msg.find( "\n" );
    std::string_view line = msg.substr( 0, index );
    yamlFile << g_level2Next << line << "\n";

    if( index == msg.npos )
      break;
    msg.remove_prefix( index + 1 );
  }
}

void ErrorLogger::write( ErrorLogger::ErrorMsg const & errorMsg ) //const
{
  std::ofstream yamlFile( std::string( m_filename ), std::ios::app );
  if( yamlFile.is_open() )
  {
    yamlFile << "\n" << g_level1Start << "type: " << g_errorLogger.toString( errorMsg.m_type ) << "\n";
    yamlFile << g_level1Next << "rank: ";
    for( size_t i = 0; i < errorMsg.m_ranksInfo.size(); i++ )
    {
      yamlFile << errorMsg.m_ranksInfo[i];
    }
    yamlFile << "\n";
    yamlFile << g_level1Next << "message: >-\n";
    streamMultilineYamlAttribute( errorMsg.m_msg, yamlFile );
    if( !errorMsg.m_contextsInfo.empty() )
    {
      yamlFile << g_level1Next << "contexts:\n";
      for( ContextInfo const & ctxInfo : errorMsg.m_contextsInfo )  
      {
        bool isFirst = true;
        for( auto const & [key, value] : ctxInfo.m_ctxInfo )
        {
          if( isFirst )
          {
              yamlFile << g_level3Start << key << ": " << value << "\n";
              isFirst = false;
          }
          else
          {
              yamlFile << g_level3Next << key << ": " << value << "\n";
          }
        }
        if( isFirst )
        {
          yamlFile << g_level3Start << "priority: " << ctxInfo.m_priority << "\n";
        }
        else 
        {
          yamlFile << g_level3Next << "priority: " <<ctxInfo.m_priority << "\n";
        }
      }
    }

    yamlFile << g_level1Next << "sourceLocation:\n";
    yamlFile << g_level2Next << "file: " << errorMsg.m_file << "\n";
    yamlFile << g_level2Next << "line: " << errorMsg.m_line << "\n";

    yamlFile << g_level1Next << "sourceCallStack:\n";

    for( size_t i = 0; i < errorMsg.m_sourceCallStack.size(); i++ )
    {
      if( i < 2 || i == errorMsg.m_sourceCallStack.size() - 1 )
        continue;
      yamlFile << g_level3Start << i-2 << errorMsg.m_sourceCallStack[i] << "\n"; 
    }

    yamlFile.flush();
    GEOS_LOG( GEOS_FMT( "The error file {} was appended.", m_filename ) );
  }
  else
  {
    GEOS_LOG( GEOS_FMT( "Unable to open error file for writing: {}", m_filename ) );
  }
}

} /* namespace geos */
