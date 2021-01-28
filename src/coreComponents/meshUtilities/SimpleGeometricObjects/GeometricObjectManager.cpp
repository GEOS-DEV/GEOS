/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
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

namespace geosx
{

using namespace dataRepository;

GeometricObjectManager::GeometricObjectManager( string const & name,
                                                Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

GeometricObjectManager::~GeometricObjectManager()
{}

Group * GeometricObjectManager::createChild( string const & childKey, string const & childName )
{
  GEOSX_LOG_RANK_0( "Adding Geometric Object: " << childKey << ", " << childName );
  std::unique_ptr< SimpleGeometricObjectBase > geometriObject = SimpleGeometricObjectBase::CatalogInterface::factory( childKey, childName, this );
  return this->registerGroup< SimpleGeometricObjectBase >( childName, std::move( geometriObject ) );
}

void GeometricObjectManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from SimpleGeometricObjectBase here
  for( auto & catalogIter: SimpleGeometricObjectBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}



} /* namespace geosx */
