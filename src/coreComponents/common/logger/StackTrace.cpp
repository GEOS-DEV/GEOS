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
 * @file StackTrace.cpp
 */

#include "StackTrace.hpp"
#include "StackTraceParams.hpp"
#include "LvArray/src/system.hpp"

#ifdef GEOS_USE_CPPTRACE
#include <cpptrace/cpptrace.hpp>
#include <cpptrace/formatting.hpp>
#endif

#include <regex>
#include <sstream>

namespace geos
{


/**
 * @brief Collect the given stack trace into frames.
 * @param[in] trace The unified trace
 * @param[inout] frames Container for the collected frames
 * @param[inout] isValidStackTrace Result of whether the given trace is a valid trace
 */
void collectLvArrayTrace( std::string const & trace,
                          std::vector< std::string > & frames,
                          bool & isValidStackTrace )
{
  std::regex lvArrayPattern( R"(Frame \d+:\s*)" );

  std::istringstream iss( trace );
  std::string stackLine;

  while( std::getline( iss, stackLine ) )
  {
    std::smatch m;
    if( std::regex_search( stackLine, m, lvArrayPattern ) )
    {
      isValidStackTrace = true;
      frames.push_back( m.suffix().str() );
    }
  }

  if( !isValidStackTrace )
  {
    frames.push_back( trace );
  }
}

#ifdef GEOS_USE_CPPTRACE
/**
 * @brief Access the configured cpptrace formatter.
 * @return The formatter instance.
 */
static cpptrace::formatter const & formatter()
{
  static cpptrace::formatter const fmt = cpptrace::formatter{}
    .header( "" )
    .addresses( cpptrace::formatter::address_mode::none )
    .paths( cpptrace::formatter::path_mode::basename )
    .symbols( cpptrace::formatter::symbol_mode::pretty )
    .snippets( false )
    .columns( true )
    .filtered_frame_placeholders( false )
    .filter( []( cpptrace::stacktrace_frame const & frame )
  {
    static char const * const blacklist[] = {
      "__cxa_throw",
      "__cxa_rethrow",
      "_Unwind_RaiseException",
      "_Unwind_Resume",
      "cpptrace::detail::",
      "cpptrace::v1::detail::",
      "cpptrace::try_catch::",
      "cpptrace::v1::try_catch::",
    };
    for( char const * symbolToHide : blacklist )
    {
      if( frame.symbol.find( symbolToHide ) != std::string::npos )
      {
        return false;
      }
    }
    return true;
  } );

  return fmt;
}

/**
 * @brief Collect the given stack trace into frames.
 * @param[in] trace The unified trace
 * @param[inout] frames Container for the collected frames
 * @param[inout] isValidStackTrace Result of whether the given trace is a valid trace
 */
void collectCpptraceTrace( cpptrace::stacktrace const & trace,
                           std::vector< std::string > & frames,
                           bool & isValidStackTrace )
{
  std::string const formatted = formatter().format( trace );
  std::regex cpptracePattern( R"(^\s*#\d+\s+)" );

  std::istringstream iss( formatted );
  std::string stackLine;

  while( std::getline( iss, stackLine ) )
  {
    std::smatch m;
    if( std::regex_search( stackLine, m, cpptracePattern ) )
    {
      isValidStackTrace = true;
      frames.push_back( m.suffix().str() );
    }
  }

  if( !isValidStackTrace )
  {
    frames.push_back( formatted );
  }
}
#endif


StackTrace::StackTrace( StackTraceParams const & params )
{
#ifdef GEOS_USE_CPPTRACE
  collectCpptraceTrace( params.trace, m_frames, m_isValidStackTrace );
#else
  collectLvArrayTrace( params.trace, m_frames, m_isValidStackTrace );
#endif
}

StackTrace StackTrace::stackTrace()
{
#ifdef GEOS_USE_CPPTRACE
  return StackTrace( StackTraceParams{ cpptrace::generate_trace( /* skip = */ 1 ) } );
#else
  return StackTrace( LvArray::system::stackTrace( true ) );
#endif
}


} /* namespace geos */
