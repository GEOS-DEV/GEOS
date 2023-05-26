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
 * @file GeometricObjectManager.cpp
 */

#include "GeometricObjectManager.hpp"

#include "SimpleGeometricObjectBase.hpp"

namespace geos
{
GeometricObjectManager * GeometricObjectManager::m_instance = nullptr;

using namespace dataRepository;

GeometricObjectManager::GeometricObjectManager( string const & name,
                                                Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  GEOS_ERROR_IF( m_instance != nullptr, "Only one GeometricObjectManager can exist at a time." );
  m_instance = this;

}

GeometricObjectManager::~GeometricObjectManager()
{
  GEOS_ERROR_IF( m_instance != this, "m_instance != this should not be possible." );
  m_instance = nullptr;
}

GeometricObjectManager & GeometricObjectManager::getInstance()
{
  GEOS_ERROR_IF( m_instance == nullptr,
                 "GeometricObjectManager has not been constructed, or is already been destructed." );
  return *m_instance;
}

Group * GeometricObjectManager::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding Geometric Object: " << childKey << ", " << childName );
  std::unique_ptr< SimpleGeometricObjectBase > geometriObject = SimpleGeometricObjectBase::CatalogInterface::factory( childKey, childName, this );
  return &this->registerGroup< SimpleGeometricObjectBase >( childName, std::move( geometriObject ) );
}

/**
   * @brief This function is used to launch a unction over the target geometric objects with region type =
   * SimpleGeometricObjectBase.
   * @tparam LOOKUP_CONTAINER type of container of names or indices
   * @tparam LAMBDA type of the user-provided function
   * @param targetObjects target geometric objects names or indices
   * @param lambda kernel function
   */
template< typename OBJECTTYPE = SimpleGeometricObjectBase, typename ... OBJECTTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
void forGeometricObject( LOOKUP_CONTAINER const & targetObjects, LAMBDA && lambda )
{
  this->forSubGroups< OBJECTTYPE, OBJECTTYPES... >( targetObjects, std::forward< LAMBDA >( lambda ) );
}

void GeometricObjectManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from SimpleGeometricObjectBase here
  for( auto & catalogIter: SimpleGeometricObjectBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}



} /* namespace geos */
