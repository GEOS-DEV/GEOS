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
#include "LvArray/src/system.hpp"

#include "common/GeosxConfig.hpp"
#ifdef GEOS_USE_CPPTRACE
#include <cpptrace/cpptrace.hpp>
#include <cpptrace/from_current.hpp>
#include <cpptrace/formatting.hpp>
#endif

namespace geos
{


std::string StackTrace::stackTrace()
{
#ifdef GEOS_USE_CPPTRACE
  return formatter().format( cpptrace::generate_trace( /* skip = */ 1 ) );
#else
  return LvArray::system::stackTrace( true );
#endif
}

std::string StackTrace::signalSafeStackTrace()
{
  return LvArray::system::stackTrace( true );
}

#ifdef GEOS_USE_CPPTRACE
cpptrace::formatter const & StackTrace::formatter()
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

std::string StackTrace::formatStackTrace( cpptrace::stacktrace const & stacktrace )
{
  return formatter().format( stacktrace );
}
#endif


} /* namespace geos */
