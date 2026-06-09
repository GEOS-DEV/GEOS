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
 * @file TemperatureProfileMPMEvent.cpp
 */

#include "TemperatureProfileMPMEvent.hpp"

namespace geos
{

using namespace dataRepository;

TemperatureProfileMPMEvent::TemperatureProfileMPMEvent( const string & name,
                                                        Group * const parent ):
  MPMEventBase( name, parent ),
  m_temperatureTable(),
  m_interpolationType( mpm::InterpolationOption::Linear )
{
  registerWrapper( viewKeyStruct::temperatureTableString(), &m_temperatureTable ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Event-local temperature table with rows { time, temperature } relative to the event start time" );

  registerWrapper( viewKeyStruct::interpolationTypeString(), &m_interpolationType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_interpolationType ).
    setRestartFlags( RestartFlags::WRITE_AND_READ ).
    setDescription( "Interpolation type of the event-local temperature table" );
}

TemperatureProfileMPMEvent::~TemperatureProfileMPMEvent()
{}

void TemperatureProfileMPMEvent::postInputInitialization()
{
  MPMEventBase::postInputInitialization();
  
  int numRows = m_temperatureTable.size( 0 );
  GEOS_ERROR_IF( numRows == 0, "Temperature table must have at least one entry" );
  for( int i = 0; i < numRows; ++i )
  {
    GEOS_ERROR_IF( m_temperatureTable[i].size() != 2, "Temperature table row " << i+1 << " must have 2 elements." );

    if( i == 0 )
    {
      GEOS_ERROR_IF( m_temperatureTable[0][0] > 0.0, "The first temperature-table time must be zero or earlier than the event start time." );
    }
    else
    {
      GEOS_ERROR_IF( ( m_temperatureTable[i][0] - m_temperatureTable[i-1][0] ) < 0.0, "Temperature table time entries must be monotonically increasing." );
    }
  }
}

REGISTER_CATALOG_ENTRY( MPMEventBase, TemperatureProfileMPMEvent, string const &, Group * const )

} /* namespace geos */