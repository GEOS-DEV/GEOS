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
  DataContext( objectName, false ),
  m_group( group )
{}
GroupContext::GroupContext( Group & group ):
  GroupContext( group, group.getName() )
{}

string GroupContext::toString() const
{
  string path;
  bool foundNearestLine = false;
  for( Group const * parentGroup = &m_group; parentGroup->hasParent(); parentGroup = &parentGroup->getParent() )
  {
    if( !foundNearestLine && parentGroup->getDataContext().isDataFileContext() )
    {
      DataFileContext const & parentContext =
        dynamic_cast< DataFileContext const & >( parentGroup->getDataContext() );
      if( parentContext.getLine() != xmlWrapper::xmlDocument::npos )
      {
        path.insert( 0, '/' + parentGroup->getName() + '(' +
                     splitPath( parentContext.getFilePath() ).second +
                     ",l." + std::to_string( parentContext.getLine() ) + ')' );
        foundNearestLine=true;
      }
      else
      {
        path.insert( 0, '/' + parentGroup->getName() );
      }
    }
    else
    {
      path.insert( 0, '/' + parentGroup->getName() );
    }
  }
  return path;
}


WrapperContext::WrapperContext( WrapperBase & wrapper ):
  GroupContext( wrapper.getParent(), wrapper.getName() )
{}

string WrapperContext::toString() const
{
  if( m_group.getDataContext().isDataFileContext() )
  {
    return m_group.getDataContext().toString() + ", attribute " + m_targetName;
  }
  else
  {
    return m_group.getDataContext().toString() + "/" + m_targetName;
  }
}


} /* namespace dataRepository */
} /* namespace geos */
