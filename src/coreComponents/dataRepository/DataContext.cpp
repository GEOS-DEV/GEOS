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

std::ostream & operator<<( std::ostream & os, DataContext const & sc )
{
  os << sc.toString();
  return os;
}


/**
 * @return the node 'name' attribute if it exists, return the node tag name otherwise.
 * @param node the target node.
 */
string getNodeName( xmlWrapper::xmlNode const & node )
{
  xmlWrapper::xmlAttribute const nameAtt = node.getAttribute( "name" );
  if( !nameAtt.empty() )
  {
    return string( node.getAttribute( "name" ).value() );
  }
  else
  {
    return string( node.name() );
  }
}

DataFileContext::DataFileContext( xmlWrapper::xmlNode const & targetNode,
                                  xmlWrapper::xmlNodePos const & nodePos ):
  DataContext( getNodeName( targetNode ) ),
  m_typeName( targetNode.name() ),
  m_filePath( nodePos.filePath ),
  m_line( nodePos.line ),
  m_offsetInLine( nodePos.offsetInLine ),
  m_offset( nodePos.offset )
{}

DataFileContext::DataFileContext( xmlWrapper::xmlNode const & targetNode,
                                  xmlWrapper::xmlAttribute const & att,
                                  xmlWrapper::xmlAttributePos const & attPos ):
  DataContext( getNodeName( targetNode ) + '/' + att.name() ),
  m_typeName( att.name() ),
  m_filePath( attPos.filePath ),
  m_line( attPos.line ),
  m_offsetInLine( attPos.offsetInLine ),
  m_offset( attPos.offset )
{}

string DataFileContext::toString() const
{
  std::ostringstream oss;
  oss << m_targetName;
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

string DataFileContext::getTargetNameInPath( bool & foundNearestLine ) const
{
  std::ostringstream oss;
  oss << m_targetName;
  foundNearestLine = ( m_line != xmlWrapper::xmlDocument::npos );
  if( foundNearestLine )
  {
    oss << "(" << splitPath( m_filePath ).second << ",l." << m_line << ")";
  }
  return oss.str();
}

} /* namespace dataRepository */
} /* namespace geos */
