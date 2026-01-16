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
#include "common/DataTypes.hpp"
#include "common/logger/Logger.hpp"
#include "common/format/StringUtilities.hpp"

#include <fstream>
#include <regex>
#include <string_view>

// signal management
#include <csignal>
#include <cfenv>
#include <cstring>

namespace geos
{
static constexpr std::string_view g_level1Start = "  - ";
static constexpr std::string_view g_level1Next =  "    ";
// static constexpr std::string_view g_level2Start = "    - "; // unused for now
static constexpr std::string_view g_level2Next =  "      ";
static constexpr std::string_view g_level3Start = "      - ";
static constexpr std::string_view g_level3Next =  "        ";

ErrorLogger g_errorLogger{};

ErrorLogger & ErrorLogger::global()
{ return g_errorLogger; }

std::mutex ErrorLogger::m_errorHandlerMutex;

std::string ErrorContext::attributeToString( ErrorContext::Attribute attribute )
{
  switch( attribute )
  {
    case ErrorContext::Attribute::InputFile:    return "inputFile";
    case ErrorContext::Attribute::InputLine:    return "inputLine";
    case ErrorContext::Attribute::DataPath:     return "dataPath";
    case ErrorContext::Attribute::DetectionLoc: return "detectionLocation";
    case ErrorContext::Attribute::Signal:       return "signal";
    default:                                    return "unknown";
  }
}

DiagnosticMsgBuilder ErrorLogger::initCurrentExceptionMessage( MsgType msgType,
                                                               std::string_view msgContent,
                                                               integer rank )
{
  DiagnosticMsg diagnosticMsg;
  m_getCurrentExceptionMsg = DiagnosticMsgBuilder::init( diagnosticMsg,
                                                         msgType, msgContent,
                                                         rank )
                               .addCallStackInfo( LvArray::system::stackTrace( true ) )
                               .getDiagnosticMsg();
  return DiagnosticMsgBuilder::modify( m_getCurrentExceptionMsg );
}

DiagnosticMsgBuilder DiagnosticMsgBuilder::init( DiagnosticMsg & msg,
                                                 MsgType msgType,
                                                 std::string_view msgContent,
                                                 integer rank )
{
  return DiagnosticMsgBuilder( msg )
           .setType( msgType )
           .addToMsg( msgContent )
           .addRank( rank );
}

DiagnosticMsgBuilder DiagnosticMsgBuilder::modify( DiagnosticMsg & errorMsg )
{
  return DiagnosticMsgBuilder( errorMsg );
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::addContextInfoImpl( ErrorContext && ctxInfo )
{
  auto lowerBoundPos = std::lower_bound( m_errorMsg.m_contextsInfo.begin(), m_errorMsg.m_contextsInfo.end(),
                                         ctxInfo.m_priority,
                                         []( ErrorContext const & ctx, integer priority )
  { return ctx.m_priority >= priority; } );
  m_errorMsg.m_contextsInfo.insert( lowerBoundPos, std::move( ctxInfo ) );
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::addDetectionLocation( string_view detectionLocation )
{
  addContextInfo( ErrorContext{  string( detectionLocation ),
                                 { { ErrorContext::Attribute::DetectionLoc,
                                   string( detectionLocation ) } } } );
  return *this;
}


DiagnosticMsgBuilder & DiagnosticMsgBuilder::addToMsg( std::exception const & e, bool const toEnd )
{
  if( toEnd )
  {
    m_errorMsg.m_msg = m_errorMsg.m_msg + e.what();
  }
  else
  {
    m_errorMsg.m_msg = e.what() + m_errorMsg.m_msg;
  }
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::addToMsg( std::string_view errorMsg, bool toEnd )
{
  if( toEnd )
  {
    m_errorMsg.m_msg = m_errorMsg.m_msg + std::string( errorMsg );
  }
  else
  {
    m_errorMsg.m_msg = std::string( errorMsg ) + m_errorMsg.m_msg;
  }
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::addSignal( integer const sig, bool const toEnd )
{
  std::string errorMsg;
  if( sig == SIGFPE )
  {
    errorMsg = "Floating point error encountered: \n";

    if( std::fetestexcept( FE_DIVBYZERO ) )
      errorMsg += "- Division by zero operation.\n";

    if( std::fetestexcept( FE_INVALID ) )
      errorMsg += "- Domain error occurred in an earlier floating-point operation.\n";

    if( std::fetestexcept( FE_OVERFLOW ) )
      errorMsg += "- The result of the earlier floating-point operation was too large to be representable.\n";

    if( std::fetestexcept( FE_UNDERFLOW ) )
      errorMsg += "- The result of the earlier floating-point operation was subnormal with a loss of precision.\n";

    if( std::fetestexcept( FE_INEXACT ) )
      errorMsg += "- Inexact result.\n";
  }
  else
  {
    errorMsg =  GEOS_FMT( "Signal encountered (no. {}): {}",
                          sig, ::strsignal( sig ));
  }

  // standard messages
  addToMsg( errorMsg, toEnd );

  this->addContextInfo( ErrorContext{  "Signal (detected from Signal Handler)",
                                       { { ErrorContext::Attribute::Signal,
                                         std::to_string( sig )  },
                                         { ErrorContext::Attribute::DetectionLoc,
                                           string( "Signal handler" )  }
                                       } } );
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::setCodeLocation( std::string_view msgFile,
                                                              integer const msgLine )
{
  m_errorMsg.m_file = msgFile;
  m_errorMsg.m_line = msgLine;
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::setType( MsgType const msgType )
{
  m_errorMsg.m_type = msgType;
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::setCause( std::string_view cause )
{
  m_errorMsg.m_cause = cause;
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::addRank( integer const rank )
{
  m_errorMsg.m_ranksInfo.emplace( rank );
  return *this;
}

DiagnosticMsgBuilder & DiagnosticMsgBuilder::addCallStackInfo( std::string_view ossStackTrace )
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
      m_errorMsg.m_isValidStackTrace = true;
      index = stackLine.find( ':' );
      m_errorMsg.m_sourceCallStack.push_back( stackLine.substr( index + 1 ) );
    }
  }

  if( !m_errorMsg.m_isValidStackTrace )
  {
    m_errorMsg.m_sourceCallStack.push_back( str );
  }

  return *this;
}

DiagnosticMsg & DiagnosticMsgBuilder::getDiagnosticMsg()
{
  return m_errorMsg;
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

std::string ErrorLogger::toString( MsgType const type )
{
  switch( type )
  {
    case MsgType::Error: return "Error";
    case MsgType::Warning: return "Warning";
    case MsgType::Exception: return "Exception";
    case MsgType::ExternalError: return "ExternalError";
    default: return "Unknown";
  }
}

/**
 * @brief Retrieve all informations from the ErrorMsg and format and write into a stream.
 * @param errMsg Class containing all the error/warning information
 * @param os The output stream to write the content to.
 */
void ErrorLogger::writeToAscii( DiagnosticMsg const & errMsg, std::ostream & os )
{
  static constexpr string_view PREFIX = "***** ";
  std::lock_guard< std::mutex > guard( m_errorHandlerMutex );
  // --- HEADER ---
  os << PREFIX << ErrorLogger::toString( errMsg.m_type ) << "\n";
  if( !errMsg.m_file.empty())
  {
    os << PREFIX<< "LOCATION: " << errMsg.m_file;
    if( errMsg.m_line > 0 )
    {
      os << ":" << errMsg.m_line;
    }
    os << "\n";
  }
  if( !errMsg.m_cause.empty())
  {
    os << PREFIX << errMsg.m_cause << "\n";
  }
  os << PREFIX << "Rank " << stringutilities::join( errMsg.m_ranksInfo, ", " ) << "\n";
  // --- ERROR CONTEXT & MESSAGE ---
  std::vector< ErrorContext > const & contexts = errMsg.m_contextsInfo;
  if( contexts.empty() || contexts.front().m_formattedContext.empty())
  {
    os << PREFIX << "Message :\n";
  }
  else
  {
    os << PREFIX << "Message from " << contexts.front().m_formattedContext << ":\n";
  }
  if( !errMsg.m_msg.empty())
  {
    os << errMsg.m_msg << "\n";
  }

  if( contexts.size() > 1 )
  {
    os << PREFIX << "Additional contexts:\n";
    for( size_t i = 1; i < contexts.size(); ++i )
    {
      os << PREFIX << "- " << contexts[i].m_formattedContext << "\n";
    }
  }
  // --- STACKTRACE ---
  if( !errMsg.m_sourceCallStack.empty() )
  {
    os << "\n** StackTrace of "<< errMsg.m_sourceCallStack.size() << " frames **\n";
    for( size_t i = 0; i < errMsg.m_sourceCallStack.size(); i++ )
    {
      os << GEOS_FMT( "Frame {}: {}\n", i, errMsg.m_sourceCallStack[i] );
    }
    os << "=====\n";
    os.flush();
  }
}

void ErrorLogger::writeToYaml( DiagnosticMsg & errMsg )
{
  std::ofstream yamlFile( std::string( m_filename ), std::ios::app );
  std::lock_guard< std::mutex > guard( m_errorHandlerMutex );
  if( yamlFile.is_open() )
  {
    // General errors info (type, rank on which the error occured)
    yamlFile << g_level1Start << "type: " << ErrorLogger::toString( errMsg.m_type ) << "\n";
    yamlFile << g_level1Next << "rank: " << stringutilities::join( errMsg.m_ranksInfo, "," );
    yamlFile << "\n";

    // Error message
    yamlFile << g_level1Next << "message: >-\n";
    streamMultilineYamlAttribute( errMsg.m_msg, yamlFile, g_level2Next );
    // context information
    if( !errMsg.m_contextsInfo.empty() )
    {
      std::vector< ErrorContext > contextInfo = errMsg.m_contextsInfo;
      // Sort contextual information by decreasing priority
      std::sort( contextInfo.begin(), contextInfo.end(),
                 []( const ErrorContext & a, const ErrorContext & b ) {
        return a.m_priority > b.m_priority;
      } );
      // Additional informations about the context of the error and priority information of each context
      yamlFile << g_level1Next << "contexts:\n";
      for( ErrorContext const & ctxInfo : contextInfo )
      {
        yamlFile << g_level3Start << "priority: " << ctxInfo.m_priority << "\n";
        if( !ctxInfo.m_formattedContext.empty())
        {
          yamlFile << g_level3Next << "description: " << ctxInfo.m_formattedContext << "\n";
        }
        for( auto const & [key, value] : ctxInfo.m_attributes )
        {
          yamlFile << g_level3Next << ErrorContext::attributeToString( key ) << ": " << value << "\n";
        }
      }
    }

    // error cause
    if( !errMsg.m_cause.empty() )
    {
      yamlFile << g_level1Next << "cause: >-\n";
      streamMultilineYamlAttribute( errMsg.m_cause, yamlFile, g_level2Next );
    }

    // Location of the error in the code
    if( !errMsg.m_file.empty() )
    {
      yamlFile << g_level1Next << "sourceLocation:\n";
      yamlFile << g_level2Next << "file: " << errMsg.m_file << "\n";
      if( errMsg.m_line > 0 )
      {
        yamlFile << g_level2Next << "line: " << errMsg.m_line << "\n";
      }
    }
    // Information about the stack trace
    if( !errMsg.m_sourceCallStack.empty() )
    {
      yamlFile << g_level1Next << "sourceCallStack:\n";
      for( size_t i = 0; i < errMsg.m_sourceCallStack.size(); i++ )
      {
        yamlFile << ( errMsg.m_isValidStackTrace ?
                      GEOS_FMT( "{}frame{}: {}\n", g_level3Start, i, errMsg.m_sourceCallStack[i] ) :
                      GEOS_FMT( "{}{}\n", g_level3Start, errMsg.m_sourceCallStack[i] ) );
      }
    }

    yamlFile << "\n";
    yamlFile.flush();
    errMsg = DiagnosticMsg();
  }
  else
  {
    GEOS_LOG_RANK( GEOS_FMT( "Unable to open error file for writing.\n- Error file: {}\n", m_filename ) );
  }
}

void ErrorLogger::flushErrorMsg( DiagnosticMsg & errMsg )
{
  writeToAscii( errMsg, m_stream );
  if( isOutputFileEnabled() )
  {
    writeToYaml( errMsg );
  }
}

void ErrorLogger::flushCurrentExceptionMessage()
{
  writeToAscii( m_getCurrentExceptionMsg, m_stream );
  if( isOutputFileEnabled() )
  {
    writeToYaml( m_getCurrentExceptionMsg );
  }
}

} /* namespace geos */
