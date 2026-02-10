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
 * @file GeosExceptions.cpp
 */


#include "common/logger/GeosExceptions.hpp"

namespace geos
{

thread_local std::ostringstream Exception::m_formattingOSS;

void Exception::prepareWhat( DiagnosticMsg & msg ) noexcept
{
  m_formattingOSS.str( "" );
  m_formattingOSS.clear();

  ErrorLogger::formatMsgForLog( msg, m_formattingOSS );
  m_cachedWhat = m_formattingOSS.bad() ? "Exception formatting error!" : m_formattingOSS.str();
}

/**
 * @brief Insert an exception message in another one.
 * @param originalMsg original exception message (i.e. thrown from GEOS_THROW)
 * @param msgToInsert message to insert at the top of the originalMsg
 */
std::string insertExMsg( std::string const & originalMsg, std::string const & msgToInsert )
{
  std::string newMsg( originalMsg );
  size_t insertPos = 0;
  // for readability purposes, we try to insert the message after the "***** Rank N: " or after "***** " instead of at the top.
  static string_view constexpr rankLogStart =  "***** Rank ";
  static string_view constexpr rankLogEnd =  ": ";
  static string_view constexpr simpleLogStart =  "***** ";
  if( ( insertPos = newMsg.find( rankLogStart ) ) != std::string::npos )
  {
    insertPos = newMsg.find( rankLogEnd, insertPos + rankLogStart.size() )
                + rankLogEnd.size();
  }
  else if( ( insertPos = newMsg.rfind( simpleLogStart ) ) != std::string::npos )
  {
    insertPos += simpleLogStart.size();
  }
  else
  {
    insertPos = 0;
  }
  newMsg.insert( insertPos, msgToInsert );
  return newMsg;
}

InputError::InputError( std::exception const & subException, std::string const & msgToInsert ):
  geos::Exception( insertExMsg( subException.what(), msgToInsert ) )
{}

SimulationError::SimulationError( std::exception const & subException, std::string const & msgToInsert ):
  geos::Exception( insertExMsg( subException.what(), msgToInsert ) )
{}

}
