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

  ErrorLogger::writeToAscii( msg, m_formattingOSS );
  m_cachedWhat = m_formattingOSS.bad() ? "Exception formatting error!" : m_formattingOSS.str();
}
}
