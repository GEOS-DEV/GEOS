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
#include "common/format/StringUtilities.hpp"

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
static constexpr std::string_view g_callStackMessage =
  "Callstack could not be retrieved. The format does not match the expected one.";

ErrorLogger g_errorLogger{};

void ErrorLogger::createFile()
{
  if( stringutilities::endsWith( m_filename, ".yaml" ) )
  {
    std::ofstream yamlFile( std::string( m_filename ), std::ios::out );
    if( yamlFile.is_open() )
    {
      yamlFile << "errors: \n\n";
      yamlFile.close();
    }
    else
    {
      GEOS_LOG_RANK( GEOS_FMT( "Unable to open error file for writing: {}", m_filename ) );
    }
  }
  else
  {
    enableFileOutput( false );
    GEOS_LOG_RANK( GEOS_FMT( "{} is a bad file name argument. The file must be in yaml format.", m_filename ) );
  }
}

std::string ErrorLogger::ErrorContext::attributeToString( ErrorLogger::ErrorContext::Attribute attribute )
{
  switch( attribute )
  {
    case ErrorLogger::ErrorContext::Attribute::InputFile: return "inputFile";
    case ErrorLogger::ErrorContext::Attribute::InputLine: return "inputLine";
    case ErrorLogger::ErrorContext::Attribute::DataPath: return "dataPath";
    default: return "Unknown";
  }
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addToMsg( std::exception const & e, bool toEnd )
{
  if( toEnd )
  {
    m_msg = m_msg + e.what();
  }
  else
  {
    m_msg = e.what() + m_msg;
  }
  return *this;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addToMsg( std::string_view errorMsg, bool toEnd )
{
  if( toEnd )
  {
    m_msg = m_msg + std::string( errorMsg );
  }
  else
  {
    m_msg = std::string( errorMsg ) + m_msg;
  }
  return *this;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::setCodeLocation( std::string_view msgFile, integer msgLine )
{
  m_file = msgFile;
  m_line = msgLine;
  return *this;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::setType( ErrorLogger::MsgType msgType )
{
  m_type = msgType;
  return *this;
}

void ErrorLogger::ErrorMsg::addContextInfoImpl( ErrorLogger::ErrorContext && ctxInfo )
{
  m_contextsInfo.emplace_back( std::move( ctxInfo ) );
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::setRank( int rank )
{
  m_ranksInfo.push_back( rank );
  return *this;
}

ErrorLogger::ErrorMsg & ErrorLogger::ErrorMsg::addCallStackInfo( std::string_view ossStackTrace )
{
  std::string str = std::string( ossStackTrace );
  std::istringstream iss( str );
  std::string stackLine;
  std::size_t index;

  std::regex pattern( R"(Frame \d+: \S+)" );

  while( std::getline( iss, stackLine ) )
  {
    if( std::regex_search( stackLine, pattern ))
    {
      m_isValidStackTrace = true;
      index = stackLine.find( ':' );
      m_sourceCallStack.push_back( stackLine.substr( index + 1 ) );
    }
  }

  if( !m_isValidStackTrace )
  {
    m_sourceCallStack.push_back( std::string( g_callStackMessage ) );
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

void ErrorLogger::streamMultilineYamlAttribute( std::string_view msg, std::ofstream & yamlFile,
                                                std::string_view indent )
{
  std::size_t i = 0;
  // Loop that runs through the string_view named msg
  while( i < msg.size() )
  {
    // Index of the next line break
    std::size_t index = msg.find( "\n", i );
    // If there is no line break, the entire string is taken
    if( index == std::string_view::npos )
    {
      index = msg.size();
    }
    // Writes the current line to the YAML file with the desired indentation
    std::string_view msgLine = msg.substr( i, index - i );
    yamlFile << indent << msgLine << "\n";
    // Move to the next line
    i = index + 1;
  }
}

void ErrorLogger::flushErrorMsg( ErrorLogger::ErrorMsg & errorMsg )
{
  std::ofstream yamlFile( std::string( m_filename ), std::ios::app );
  if( yamlFile.is_open() && g_errorLogger.isOutputFileEnabled() )
  {
    // General errors info (type, rank on which the error occured)
    yamlFile << g_level1Start << "type: " << g_errorLogger.toString( errorMsg.m_type ) << "\n";
    yamlFile << g_level1Next << "rank: ";
    for( auto const & info: errorMsg.m_ranksInfo )
    {
      yamlFile << info;
    }
    yamlFile << "\n";

    // Error message
    yamlFile << g_level1Next << "message: >-\n";
    streamMultilineYamlAttribute( errorMsg.m_msg, yamlFile, g_level2Next );

    // context information
    if( !errorMsg.m_contextsInfo.empty() )
    {
      // Sort contextual information by decreasing priority
      std::sort( errorMsg.m_contextsInfo.begin(), errorMsg.m_contextsInfo.end(),
                 []( const ErrorLogger::ErrorContext & a, const ErrorLogger::ErrorContext & b ) {
        return a.m_priority > b.m_priority;
      } );
      // Additional informations about the context of the error and priority information of each context
      yamlFile << g_level1Next << "contexts:\n";
      for( ErrorContext const & ctxInfo : errorMsg.m_contextsInfo )
      {
        yamlFile << g_level3Start << "priority: " << ctxInfo.m_priority << "\n";
        for( auto const & [key, value] : ctxInfo.m_attributes )
        {
          yamlFile << g_level3Next << ErrorContext::attributeToString( key ) << ": " << value << "\n";
        }
      }
    }

    // Location of the error in the code
    yamlFile << g_level1Next << "sourceLocation:\n";
    yamlFile << g_level2Next << "file: " << errorMsg.m_file << "\n";
    yamlFile << g_level2Next << "line: " << errorMsg.m_line << "\n";

    // Information about the stack trace
    yamlFile << g_level1Next << "sourceCallStack:\n";
    if( !errorMsg.isValidStackTrace() )
    {
      yamlFile << g_level3Start << "callStackMessage: " << errorMsg.m_sourceCallStack[0] << "\n";
    }
    else
    {
      for( size_t i = 0; i < errorMsg.m_sourceCallStack.size(); i++ )
      {
        yamlFile << g_level3Start << "frame" << i << ": " << errorMsg.m_sourceCallStack[i] << "\n";
      }
    }

    yamlFile << "\n";
    yamlFile.flush();
    errorMsg = ErrorMsg();
    GEOS_LOG_RANK( GEOS_FMT( "The error file {} was appended.", m_filename ) );
  }
  else
  {
    GEOS_LOG_RANK( GEOS_FMT( "Unable to open error file for writing.\n- Error file: {}\n- Error file enabled = {}.\n",
                             m_filename, g_errorLogger.isOutputFileEnabled() ) );
  }
}

} /* namespace geos */
