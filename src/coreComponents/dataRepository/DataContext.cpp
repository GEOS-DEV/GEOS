/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DataContext.cpp
 */

#include "DataContext.hpp"

namespace geos
{
namespace dataRepository
{


DataContext::DataContext( string const & targetName ):
  m_targetName( targetName )
{}

std::ostream & operator<<( std::ostream & os, DataContext const & ctx )
{
  os << ctx.toString();
  return os;
}

DataContext::ToStringInfo::ToStringInfo( string const & targetName, string const & filePath, size_t line ):
  m_targetName( targetName ),
  m_filePath( filePath ),
  m_line( line )
{}

DataContext::ToStringInfo::ToStringInfo( string const & targetName ):
  m_targetName( targetName )
{}



} /* namespace dataRepository */
} /* namespace geos */
