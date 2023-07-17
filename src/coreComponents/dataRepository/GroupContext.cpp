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
 * @file GroupContext.cpp
 */

#include "GroupContext.hpp"

namespace geos
{
namespace dataRepository
{


GroupContext::GroupContext( Group & group, string const & objectName ):
  DataContext( objectName ),
  m_group( group )
{}
GroupContext::GroupContext( Group & group ):
  GroupContext( group, group.getName() )
{}

string GroupContext::toString() const
{
  // we add all parent names in a path string, showing only the first input file info encountered
  std::vector< string > parents;
  bool foundNearestLineInfo = false;
  for( Group const * parent = &m_group; parent->hasParent(); parent = &parent->getParent() )
  {
    ToStringInfo const info = parent->getDataContext().getToStringInfo();
    if( info.hasInputFileInfo() && !foundNearestLineInfo )
    {
      // avoiding spaces here is intended as we don't want any line return within the path.
      parents.push_back( GEOS_FMT( "{}({},l.{})", info.m_targetName, info.m_filePath, info.m_line ) );
      foundNearestLineInfo = true;
    }
    else
    {
      parents.push_back( string( info.m_targetName ) );
    }
  }
  return '/' + stringutilities::join( parents.rbegin(), parents.rend(), '/' );
}

DataContext::ToStringInfo GroupContext::getToStringInfo() const
{ return ToStringInfo( m_targetName ); }


WrapperContext::WrapperContext( WrapperBase & wrapper ):
  GroupContext( wrapper.getParent(), wrapper.getParent().getName() + '/' + wrapper.getName() ),
  m_typeName( wrapper.getName() )
{}

string WrapperContext::toString() const
{
  ToStringInfo const info = m_group.getDataContext().getToStringInfo();
  return info.hasInputFileInfo() ?
         GEOS_FMT( "{} ({}, l.{})", m_targetName, info.m_filePath, info.m_line ) :
         GEOS_FMT( "{}->{}", m_group.getDataContext().toString(), m_typeName );
}


} /* namespace dataRepository */
} /* namespace geos */
