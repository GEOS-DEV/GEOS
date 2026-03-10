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

MPMEventBase::MPMEventBase( string const & name,
                            Group * const parent ):
  Group( name, parent ),
  m_startTime( 0.0 ),
  m_endTime( 1e16 ), // Might overflow if set to DBL_MAX
  m_isComplete( 0 )
{
  registerWrapper( viewKeyStruct::startTimeString(), &m_startTime ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Time at which event starts" );

  registerWrapper( viewKeyStruct::endTimeString(), &m_endTime ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( DBL_MAX ).
    setDescription( "Time at which event ends" );

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
  GEOS_ERROR_IF( m_startTime > m_endTime, "Event start time must be less than end time!" );
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
