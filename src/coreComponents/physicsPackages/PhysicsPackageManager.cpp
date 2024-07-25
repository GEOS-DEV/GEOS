/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhysicsPackageManager.cpp
 */

#include "PhysicsPackageManager.hpp"

#include "PhysicsPackageBase.hpp"

namespace geos
{

using namespace dataRepository;

PhysicsPackageManager::PhysicsPackageManager( string const & name,
                                            Group * const parent ):
  Group( name, parent ),
  m_gravityVector( { 0.0, 0.0, -9.81 } )
{
  setInputFlags( InputFlags::REQUIRED );

  this->registerWrapper( viewKeyStruct::gravityVectorString(), &m_gravityVector ).
    setDefaultValue( m_gravityVector ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Gravity vector used in the physics solvers" );
}

PhysicsPackageManager::~PhysicsPackageManager()
{}


//START_SPHINX_INCLUDE_00
Group * PhysicsPackageManager::createChild( string const & childKey, string const & childName )
{
  Group * rval = nullptr;
  if( PhysicsPackageBase::CatalogInterface::hasKeyName( childKey ) )
  {
    GEOS_LOG_RANK_0( "Adding Solver of type " << childKey << ", named " << childName );
    rval = &registerGroup( childName,
                           PhysicsPackageBase::CatalogInterface::factory( childKey, childName, this ) );
  }
  return rval;
}


void PhysicsPackageManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from PhysicsPackageBase here
  for( auto & catalogIter: PhysicsPackageBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geos */
