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

#include "Group.hpp"
#include "WrapperBase.hpp"
#include "../mainInterface/ProblemManager.hpp"

namespace geos
{
namespace dataRepository
{


DataContext::DataContext( string const & objectName, bool const isDataFileContext ):
  m_objectName( objectName ),
  m_isDataFileContext( isDataFileContext )
{}

std::ostream & operator<<( std::ostream & os, DataContext const & sc )
{
  os << sc.toString();
  return os;
}


GroupContext::GroupContext( Group & group, string const & objectName ):
  DataContext( objectName, false ),
  m_group( group )
{}
GroupContext::GroupContext( Group & group ):
  GroupContext( group, group.getName() )
{}

string GroupContext::toString() const
{
  // it would be possible to insert the DataFileContext::toString() of parent objects when it exists, but is it relevant ?
  return m_group.getPath();
}


WrapperContext::WrapperContext( WrapperBase & wrapper ):
  GroupContext( wrapper.getParent(), wrapper.getName() )
{}

string WrapperContext::toString() const
{
  // if possible, we show the DataFileContext of the parent.
  if( m_group.getDataContext().isDataFileContext() )
  {
    return m_group.getDataContext().toString() + ", attribute " + m_objectName;
  }
  else
  {
    // it would be possible to insert the DataFileContext::toString() of parent objects when it exists, but is it relevant ?
    return m_group.getPath() + "/" + m_objectName;
  }
}


DataFileContext::DataFileContext( Group & group, xmlWrapper::xmlNodePos const & nodePos,
                                  string const & nodeTagName ):
  DataContext( group.getName(), true ),
  m_typeName( nodeTagName ),
  m_filePath( nodePos.filePath ),
  m_line( nodePos.line ),
  m_offsetInLine( nodePos.offsetInLine ),
  m_offset( nodePos.offset )
{}

DataFileContext::DataFileContext( WrapperBase & wrapper, xmlWrapper::xmlAttributePos const & attPos ):
  DataContext( wrapper.getParent().getName() + "/" + wrapper.getName(), true ),
  m_typeName( wrapper.getName() ),
  m_filePath( attPos.filePath ),
  m_line( attPos.line ),
  m_offsetInLine( attPos.offsetInLine ),
  m_offset( attPos.offset )
{}

string DataFileContext::toString() const
{
  std::ostringstream oss;
  oss << m_objectName;
  if( m_line != xmlWrapper::xmlDocument::npos )
  {
    oss << " (" << splitPath( m_filePath ).second << ", l." << m_line << ")";
  }
  else if( m_offset != xmlWrapper::xmlDocument::npos )
  {
    oss << " (" << splitPath( m_filePath ).second <<  ", offset " << m_offset << ")";
  }
  else
  {
    oss << " (Source file not found)";
  }
  return oss.str();
}


} /* namespace dataRepository */
} /* namespace geos */
