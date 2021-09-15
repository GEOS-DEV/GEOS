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
 * @file PhysicsSolverManager.cpp
 */

#include "PhysicsSolverManager.hpp"

#include "SolverBase.hpp"

namespace geosx
{

using namespace dataRepository;

PhysicsSolverManager::PhysicsSolverManager( string const & name,
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

PhysicsSolverManager::~PhysicsSolverManager()
{}


//START_SPHINX_INCLUDE_00
Group * PhysicsSolverManager::createChild( string const & childKey, string const & childName )
{
  Group * rval = nullptr;
  if( SolverBase::CatalogInterface::hasKeyName( childKey ) )
  {
    GEOSX_LOG_RANK_0( "Adding Solver of type " << childKey << ", named " << childName );
    rval = &registerGroup( childName,
                           SolverBase::CatalogInterface::factory( childKey, childName, this ) );
  }
  return rval;
}


void PhysicsSolverManager::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from SolverBase here
  for( auto & catalogIter: SolverBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}


} /* namespace geosx */
