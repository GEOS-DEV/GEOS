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
 * @file MPMEventBase.cpp
 */

#include "MPMEventBase.hpp"

namespace geos
{

using namespace dataRepository;

// Events can optionally be triggered two ways
// - it queries whether it's dependencies have been met and then must run within 

MPMEventBase::MPMEventBase( string const & name,
                            Group * const parent ):
  Group( name, parent ),
  m_startTime( -1.0 ),
  m_endTime( DBL_MAX ),
  m_delay( 0.0 ),
  m_duration( 1e16 ),
  m_hasStarted( 0 ),
  m_isComplete( 0 )
{
  registerWrapper( viewKeyStruct::startTimeString(), &m_startTime ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_startTime ).
    setDescription( "Time at which event starts" );

  registerWrapper( viewKeyStruct::endTimeString(), &m_endTime ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_endTime ).
    setDescription( "Time at which event ends" );

  registerWrapper( viewKeyStruct::delayString(), &m_delay ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_delay ).
    setDescription( "Delay between dependencies being complete and starting the event" );

  registerWrapper( viewKeyStruct::durationString(), &m_duration ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_duration ).
    setDescription( "Time at which event ends" );

  registerWrapper( viewKeyStruct::dependenciesString(), &m_dependencies ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of the names for event dependencies" );

  registerWrapper( viewKeyStruct::hasStartedString(), &m_hasStarted ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Flag for whether event has started" );

  registerWrapper( viewKeyStruct::isCompleteString(), &m_isComplete ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Flag for whether event has been completed" );
}


MPMEventBase::~MPMEventBase()
{}


MPMEventBase::CatalogInterface::CatalogType & MPMEventBase::getCatalog()
{
  static MPMEventBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void MPMEventBase::postInputInitialization()
{
  if(m_dependencies.size() == 0)
  {
    GEOS_ERROR_IF( m_startTime < 0.0, getName() << " event start time must be less than or equal to 0.0 or was not specified");
    GEOS_ERROR_IF( m_startTime > m_endTime, getName() << " event start time must be less than end time!" );
  }
  else
  {
    GEOS_ERROR_IF( m_duration < 0.0, getName() << " event duration must be positive" );
  }
  
}

Group * MPMEventBase::createChild( string const & childKey, string const & childName )
{
  // GEOS_LOG_RANK_0( "Adding MPM Event: " << childKey << ", " << childName );
  // std::unique_ptr< MPMEventBase > event = MPMEventBase::CatalogInterface::factory( childKey, childName, this );
  // return &this->registerGroup< MPMEventBase >( childName, std::move( event ) );

  GEOS_LOG_RANK_0( GEOS_FMT( "{}: adding {} {}", getName(), childKey, childName ) );
  std::unique_ptr< MPMEventBase > event = MPMEventBase::CatalogInterface::factory( childKey, getDataContext(), childName, this );
  return &this->registerGroup< MPMEventBase >( childName, std::move( event ) );
}


void MPMEventBase::expandObjectCatalogs()
{
  // Only add children if the parent is of type EventManager
  // otherwise, this would fall into a loop
  if( strcmp( this->getParent().getName().c_str(), "MPMEvents" ) == 0 )
  {
    for( auto & catalogIter: MPMEventBase::getCatalog() )
    {
      createChild( catalogIter.first, catalogIter.first );
    }
  }
}


} /* namespace geos */
