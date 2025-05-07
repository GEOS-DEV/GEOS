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
                                                    Group * const parent ) :
                                                    MPMEventBase(  name, parent ),
                                                    m_temperatureTable()
  {  
    registerWrapper( viewKeyStruct::temperatureTableString(), &m_temperatureTable ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Temperature profile specified as a table" );

    registerWrapper( viewKeyStruct::interpTypeString(), &m_interpType ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue( m_interpType ).
        setDescription( "Interpolation scheme: 0 (Linear), 1 (Cosine), 2 (Sigmoid)" ); // CC: double check this!
  }

  TemperatureProfileMPMEvent::~TemperatureProfileMPMEvent() 
  {}

  void TemperatureProfileMPMEvent::postProcessInput()
  {
    int numRows = m_temperatureTable.size( 0 );
    GEOS_ERROR_IF(numRows == 0, "Temperature profile event was included but temperature profile was specified.");
    for(int i = 0; i < numRows; ++i)
    {
      GEOS_ERROR_IF(m_temperatureTable[i].size() != 2, "Temperature profile must consists of a time column and temperature column" );
      GEOS_ERROR_IF(m_temperatureTable[i][0] < 0, "Temperature profile times must be positive.");
    }

    GEOS_LOG_RANK_0( "TemperatureProfileEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "interpType=" << m_interpType );
                     //TODO write temperature table to console
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, TemperatureProfileMPMEvent, string const &, Group * const )

} /* namespace geos */
