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
  string path;
  bool foundNearestLine = false;
  for( Group const * parentGroup = &m_group; parentGroup->hasParent(); parentGroup = &parentGroup->getParent() )
  {
    if( !foundNearestLine )
    {
      path.insert( 0, '/' + parentGroup->getDataContext().getTargetNameInPath( foundNearestLine ) );
    }
    else
    {
      path.insert( 0, '/' + parentGroup->getName() );
    }
  }
  return path;
}


WrapperContext::WrapperContext( WrapperBase & wrapper ):
  GroupContext( wrapper.getParent(), wrapper.getParent().getName() + '/' + wrapper.getName() ),
  m_typeName( wrapper.getName() )
{}

string WrapperContext::toString() const
{
  DataContext const & parentDC = m_group.getDataContext();
  return parentDC.toString() + parentDC.getWrapperSeparator() + m_typeName;
}


} /* namespace dataRepository */
} /* namespace geos */
