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
 * @file SourceContext.cpp
 */

//#include "codingUtilities/StringUtilities.hpp"
#include "Group.hpp"
#include "WrapperBase.hpp"
#include "../mainInterface/ProblemManager.hpp"

namespace geos
{
namespace dataRepository
{


SourceContext::SourceContext( string const & objectName, bool const isFileContext ):
  m_objectName( objectName ),
  m_isFileContext( isFileContext )
{}

std::ostream & operator<<( std::ostream & os, SourceContext const & sc )
{
  os << sc.toString();
  return os;
}


GroupContext::GroupContext( Group & group, string const & objectName ):
  SourceContext( objectName, false ),
  m_group( group )
{}
GroupContext::GroupContext( Group & group ):
  GroupContext( group, group.getName() )
{}

string GroupContext::toString() const
{
  // it would be possible to insert the FileContext::toString() of parent objects when it exists, but is it relevant ?
  return m_group.getPath();
}


WrapperContext::WrapperContext( WrapperBase & wrapper ):
  GroupContext( wrapper.getParent(), wrapper.getName() )
{}

string WrapperContext::toString() const
{
  // if possible, we show the FileContext of the parent.
  if( m_group.getSourceContext().isFileContext() )
  {
    return m_group.getSourceContext().toString() + ", attribute " + m_objectName;
  }
  else
  {
    // it would be possible to insert the FileContext::toString() of parent objects when it exists, but is it relevant ?
    return m_group.getPath() + "/" + m_objectName;
  }
}


FileContext::FileContext( Group & group, xmlWrapper::xmlNodePos const & nodePos,
                          string const & nodeTagName ):
  SourceContext( group.getName(), true ),
  m_typeName( nodeTagName ),
  m_filePath( nodePos.filePath ),
  m_line( nodePos.line ),
  m_offsetInLine( nodePos.offsetInLine ),
  m_offset( nodePos.offset )
{}

FileContext::FileContext( WrapperBase & wrapper, xmlWrapper::xmlAttributePos const & attPos ):
  SourceContext( wrapper.getParent().getName() + "/" + wrapper.getName(), true ),
  m_typeName( wrapper.getName() ),
  m_filePath( attPos.filePath ),
  m_line( attPos.line ),
  m_offsetInLine( attPos.offsetInLine ),
  m_offset( attPos.offset )
{}

string FileContext::toString() const
{
  std::ostringstream oss;
  oss << m_objectName << " (" << m_filePath;
  if( m_line != xmlWrapper::xmlDocument::npos )
  {
    oss << ": l." << m_line << ")";
  }
  else
  {
    // line hasn't been found, filename is probably wrong too, we just output the character offset.
    oss <<  ": offset " << m_offset << ")";
  }
  return oss.str();
}


} /* namespace dataRepository */
} /* namespace geos */
